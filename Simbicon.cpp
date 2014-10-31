#include <Simbody.h>
#include "Humanoid.h"
#include "Simbicon.h"

#include <cassert>

using namespace SimTK;

// Normally Simbicon has 4 states per gait cycle (i.e. 2 per leg), 2 state is
// simplified and not as realistic.
#define TWO_STATE
// Torque angle to be flat relative to ground (not in original Simbicon paper)
//#define USE_GLOBAL_ANKLE
// Aim hip towards a global orientation (prevents meandering) (Not implemented
// for rigid contact because it depends on checking force above some maximum.)
//#define USE_GLOBAL_HIPROT
// Drop landing doesn't use controller so you can run with models other than
// the humanoid upon which the controller depends.
//#define DROP_LANDING

namespace { // file-scope symbols

// Controller gains, given by "strength" in N-m/radian. Then kp=strength
// and kd=2*sqrt(strength) for critical damping.
#define USE_ORIG_GAINS
#ifdef USE_ORIG_GAINS
const Real DefaultStrength              = 300;  
const Real NeckStrength                 = 100;  
const Real BackStrength                 = 300;  
const Real HipFlexionAdductionKp  = 1000;  // must be overdamped
const Real HipFlexionAdductionKd  = 100; 
const Real HipRotationStrength          = 300;  
const Real KneeStrength                 = 300;  
const Real ArmFlexionAdductionStrength  = 300;  
const Real ArmRotationStrength          = 300;  
const Real AnkleFlexionStrength         = 300;  
const Real AnkleInversionStrength       = 30;  //300
const Real ToeStrength                  = 30;   
#else
const Real DefaultStrength              = 300;  // was 300,30
const Real NeckStrength                 = 20;  // was 100,10
const Real BackStrength                 = 200;  // was 300,30
const Real HipFlexionAdductionStrength  = 500;  // was 1000,100
const Real HipRotationStrength          = 100;  // was 300,30
const Real KneeStrength                 = 300;  // was 300,30
const Real ArmFlexionAdductionStrength  = 10;  // was 300,30
const Real ArmRotationStrength          = 5;  // was 300,30
const Real AnkleFlexionStrength         = 100;  // was 300,30
const Real AnkleInversionStrength       = 50;   // was 300,30
const Real ToeStrength                  = 25;   // was 30,3
#endif

const Vec3 UnitX(1.0, 0.0, 0.0); 
const Vec3 UnitY(0.0, 1.0, 0.0); 
const Vec3 UnitZ(0.0, 0.0, 1.0);

// Convert muscle strength into critically damped control gains.
void calcGainsFromStrength(double strength, double& kp, double& kd) {
    kp = strength;
    kd = 2*std::sqrt(strength);
}

double clamp( double x, double maxTorque ) {
	if (x > maxTorque) 
		x = maxTorque;  
	
	if (x < -maxTorque) 
		x = -maxTorque;  

	return x; 
}
}

Simbicon::Simbicon(Humanoid& model) : _model(model) {
	_state = -1;
	_stateStartTime = -1.0;
	for (int i = 0; i < 2; i++) {  
		_lastSWTAngle[i] = -100.0; 
		_curSWTAngle[i] = -100.0; 
		_lastTrunkAngle[i] = -100.0; 
		_curTrunkAngle[i] = -100.0; 
		_lastRFootAngle[i] = -100.0; 
		_curRFootAngle[i] = -100.0; 
		_lastLFootAngle[i] = -100.0; 
		_curLFootAngle[i] = -100.0; 
	}
	_lastPelvisRotation = -100.0; 
	_curPelvisRotation = -100.0; 

}

void Simbicon::calcForce(const State&         s, 
                         Vector_<SpatialVec>& /*bodyForces*/, 
                         Vector_<Vec3>&       /*particleForces*/, 
                         Vector&              mobilityForces) const 
{
    Vector controls(NumActuators);
    computeControls(s,controls);

    for (int i=0; i<NumActuators; ++i)
        mobilityForces[_model.getUIndex(Actuator(i))] = controls[i];
}


// Project the pelvis z (right) and x (forward) directions onto the Ground
// (x-z) plane, and normalize. Projected z is the normal to the sagittal plane;
// projected x is the normal to the coronal plane.
void Simbicon::getSagCorNormals(const State& s, Vec3& sagN, Vec3& corN ) const {
	const MobilizedBody& pelvis = _model.pelvis;
	
    sagN = Vec3(pelvis.getBodyRotation(s).z());
    corN = Vec3(pelvis.getBodyRotation(s).x());
	sagN[YAxis] = 0; // project to y==0 
	sagN = sagN.normalize();  
	corN[YAxis] = 0; // project to y==0
	corN = corN.normalize();  
}

void Simbicon::getUpVectorInGround( const State& s, const MobilizedBody& b, 
	Vec3& up ) const { 
	up = Vec3(b.getBodyRotation(s).y());   
}

void Simbicon::getFrontVectorInGround( const State& s, const MobilizedBody& b, 
	Vec3& front ) const { 
	front = Vec3(b.getBodyRotation(s).x());
}

void Simbicon::fillInHipJointControls( const State& s, Vector& controls ) const {
	const SimbodyMatterSubsystem& matter = _model.matter;
	Vec3 sagN, corN; 	
	
    // Sagittal and coronal plane normals are right and front directions of
    // the pelvis, projected on the ground plane.
	getSagCorNormals(s, sagN, corN); 	
	
	int swh = hip_r_flexion;             // swing hip
	int sth = hip_l_flexion;             // stance hip
    MobilizedBody ankle = _model.foot_l; // stance ankle
	
	if (_state == 2 || _state == 3) { // right stance
		swh = hip_l_flexion; 
		sth = hip_r_flexion; 
        ankle = _model.foot_r;
	}
	int swhc = swh - 1; // coronal plane (hip adduction)
	int sthc = sth - 1; // stance hip adduction
	int sta = sth + 4;  // stance ankle dorsiflexion 

    const Transform& X_PF = ankle.getInboardFrame(s);
    Vec3 ankleLocInParent = X_PF.p();
    Vec3 ankleLoc = ankle.getParentMobilizedBody()
                         .findStationLocationInGround(s,ankleLocInParent);

	Vec3 com = matter.calcSystemMassCenterLocationInGround(s); 
	Vec3 d = com - ankleLoc; 
	Vec3 v_com = matter.calcSystemMassCenterVelocityInGround(s); 

	double d_sag = dot(corN, d);
	double v_sag = dot(corN, v_com);
	double d_cor = dot(sagN, d);
	double v_cor = dot(sagN, v_com);
	
	double thetad = 0.5; 
#ifndef TWO_STATE	
	if (_state == 1 || _state == 3) { 
		thetad = -0.1; 
	}
#endif
    double kp, kd; // position, derivative gains for hip flex/adduction
    #ifdef USE_ORIG_GAINS
    kp = HipFlexionAdductionKp; kd = HipFlexionAdductionKd;
    #else
    calcGainsFromStrength(HipFlexionAdductionStrength, kp, kd);
    #endif
	double cd = 0.2; // global tipping feedback
	double cv = 0.2;
	
	double trunkAngleVelEst[2] = {0, 0};
	double SWTAngleVelEst[2] = {0, 0};
	double RFootAngleVelEst[2] = {0, 0};
	double LFootAngleVelEst[2] = {0, 0};
	double PelvisRotationVelEst = 0;
	if (_lastTrunkAngle[0] > -100) { 
		// check there's a valid value for the _last*, 
		// otherwise just use 0 for vel
		for (int i = 0; i < 2; i++) {
			trunkAngleVelEst[i] = 
				(_curTrunkAngle[i] - _lastTrunkAngle[i])/STATE_UPD_STEPSIZE; 
			RFootAngleVelEst[i] = 
				(_curRFootAngle[i] - _lastRFootAngle[i])/STATE_UPD_STEPSIZE; 
			LFootAngleVelEst[i] = 
				(_curLFootAngle[i] - _lastLFootAngle[i])/STATE_UPD_STEPSIZE;
			if (_lastSWTAngle[i] > -100) { 
				SWTAngleVelEst[i] = 
					(_curSWTAngle[i] - _lastSWTAngle[i])/STATE_UPD_STEPSIZE; 
			}
		}
		PelvisRotationVelEst = 
			(_curPelvisRotation - _lastPelvisRotation)/STATE_UPD_STEPSIZE; 
	}
	 

	// sign change is needed for one of the stance legs in the coronal plane
	double sign = 1; 
	if (_state == 0 || _state == 1) // left stance	
		sign = -1;
	
	controls[swh] = clamp(  kp*(thetad + (cd*d_sag + cv*v_sag) - _curSWTAngle[0]) 
                          - kd*SWTAngleVelEst[0], kp);
	controls[swhc] = sign*(clamp(  kp*(0.0 + (cd*d_cor + cv*v_cor) - _curSWTAngle[1]) 
                                 - kd*SWTAngleVelEst[1], kp));
	
	// use stance hip to control the trunk
	controls[sth] =  -kp*(0. - _curTrunkAngle[0]) + kd*trunkAngleVelEst[0];   
	controls[sth] -= controls[swh];
	controls[sth] = clamp(controls[sth], kp);
		
	controls[sthc] = sign*(kp*(0. - _curTrunkAngle[1]) - kd*trunkAngleVelEst[1]);   
	controls[sthc] -= controls[swhc];
	controls[sthc] = clamp(controls[sthc], kp);
	
#ifdef USE_GLOBAL_ANKLE
    double kpaflex, kdaflex, kpainv, kdainv; // gains for ankle flex, inversion
    calcGainsFromStrength(AnkleFlexionStrength, kpaflex, kdaflex);
    calcGainsFromStrength(AnkleInversionStrength, kpainv, kdainv);
	controls[ankle_r_dorsiflexion] = 
        clamp(  kpaflex*(0. - _curRFootAngle[0]) 
              - 0.*kdaflex*RFootAngleVelEst[0], kpaflex);   
	controls[ankle_r_inversion] = 
        clamp( -kpainv*(0. - _curRFootAngle[1]) 
              + kdainv*RFootAngleVelEst[1], kpainv);   
	controls[ankle_l_dorsiflexion] = 
        clamp(  kpaflex*(0. - _curLFootAngle[0]) 
              - 0.*kdaflex*LFootAngleVelEst[0], kpaflex);   
	controls[ankle_l_inversion] = 
        clamp(  kpainv*(0. - _curLFootAngle[1]) 
              - kdainv*LFootAngleVelEst[1], kpainv);  
#endif

#ifdef USE_GLOBAL_HIPROT
    double kphrot, kdhrot; // gains for hip rotation
    calcGainsFromStrength(HipRotationStrength, kphrot, kdhrot);
	if (   (sth == hip_r_flexion && _curRFootContactForce > 100) 
        || (sth == hip_l_flexion  && _curLFootContactForce > 100)) 
    {
		controls[sth+1] = clamp(sign*( -kphrot*(0. - _curPelvisRotation) 
                                      + kdhrot*PelvisRotationVelEst ),
                                kphrot); 
	}
#endif

}

void Simbicon::computeControls(const State& s, Vector& controls) const
{
#ifdef DROP_LANDING
	for (int i = 0; i < controls.size(); i++) 
			controls[i] = 0.0; 
	return; 
#else
	int swh = hip_r_flexion; 
	int sth = hip_l_flexion; 
	
	if (_state == 2 || _state == 3) { // right stance
		swh = hip_l_flexion; 
		sth = hip_r_flexion; 
	}
	int swk = swh + 2; // swing knee
	int swa = swh + 4; // swing ankle
	int stk = sth + 2; // stance knee
	int sta = sth + 4; // stance ankle 
	for (int i = 0; i < NumActuators; i++) {
        double kp, kd;          // position gain, derivative gain
        calcGainsFromStrength(DefaultStrength, kp, kd);
		double thetad = 0.0;                    // desired angle
		
		if (i == neck_extension || i == neck_bending || i == neck_rotation) {
            calcGainsFromStrength(NeckStrength, kp, kd);
		}
		else if (i == back_tilt || i == back_list || i == back_rotation) {
            calcGainsFromStrength(BackStrength, kp, kd);
		}
		else if (   i == shoulder_r_flexion   || i == shoulder_l_flexion 
                 || i == shoulder_r_adduction || i == shoulder_l_adduction
                 || i == elbow_r_flexion      || i == elbow_l_flexion) {
            calcGainsFromStrength(ArmFlexionAdductionStrength, kp, kd);
		}
		else if (   i == shoulder_r_rotation || i == shoulder_l_rotation 
                 || i == elbow_r_rotation    || i == elbow_l_rotation) {
            calcGainsFromStrength(ArmRotationStrength, kp, kd);
		}
		else if (i == hip_r_rotation || i == hip_l_rotation) {
            calcGainsFromStrength(HipRotationStrength, kp, kd);
		}
		else if (i == knee_r_extension || i == knee_l_extension) {
            calcGainsFromStrength(KneeStrength, kp, kd);
		}
        else if (i == ankle_r_dorsiflexion || i == ankle_l_dorsiflexion) {
            calcGainsFromStrength(AnkleFlexionStrength, kp, kd);
        }
        else if (i == ankle_r_inversion || i == ankle_l_inversion) {
            calcGainsFromStrength(AnkleInversionStrength, kp, kd);
        }
		else if (i == mtp_r_dorsiflexion || i == mtp_l_dorsiflexion) { // toe
            calcGainsFromStrength(ToeStrength, kp, kd);
		}

		if (_state >= 0) {
			if (i == swk) 
				thetad = -1.1;  
			else if (i == swa) 
				thetad = 0.6;  
			else if (i == stk) 
				thetad = -0.05;

            #ifndef TWO_STATE	
			if (_state == 1 || _state == 3) { 
				if (i == swk) 
					thetad = -0.05;  
				else if (i == swa) 
					thetad = 0.15;  
				else if (i == stk) 
					thetad = -0.1; 
			}
            #endif
		}

        const Real qi = s.getQ()[_model.getQIndex(Actuator(i))];
        const Real ui = s.getU()[_model.getUIndex(Actuator(i))];
		controls[i] = clamp(kp*(thetad - qi) - kd*ui, kp); 
	}
	
	if (_state >= 0) 
		fillInHipJointControls(s, controls); 

	Vec3 com = _model.matter.calcSystemMassCenterLocationInGround(s); 
	if (com[1] < 0.7) {
		for (int i = 0; i < controls.size(); i++) 
			controls[i] = 0.0; 
	}
	return;
#endif
}
	
void Simbicon::
computeSecondaryStateVals(const State& s, Real lForce, Real rForce) {
	const MobilizedBody& pelvis = _model.pelvis;

	Vec3 upThigh;
	if (_state == 0 || _state == 1)  
		getUpVectorInGround(s, _model.thigh_r, upThigh); 
	else if (_state == 2 || _state == 3) 
		getUpVectorInGround(s, _model.thigh_l, upThigh); 
	
	Vec3 upPelvis; 
	getUpVectorInGround(s, pelvis, upPelvis); 
	
	Vec3 sagN, corN; 
	getSagCorNormals(s, sagN, corN); 	

	Vec3 XinSag = cross(UnitY, sagN); 
	Vec3 ZinCor = -cross(UnitY, corN); 

	// Store the current value, use these for velocity estimation
	for (int i = 0; i < 2; i++) {	
		_lastSWTAngle[i] = _curSWTAngle[i];
		_lastTrunkAngle[i] = _curTrunkAngle[i];
		_lastLFootAngle[i] = _curLFootAngle[i];
		_lastRFootAngle[i] = _curRFootAngle[i];
	}
	_lastPelvisRotation = _curPelvisRotation;

	// Update trunk and swing thigh global orientations in the saggital and 
	// coronal planes (by projecting the up vectors of the bodies in to the 
	// planes and calculate angles)
	Vec3 projUpThigh = upThigh - dot(upThigh, sagN)*sagN;
	_curSWTAngle[0] = 
		acos(dot(projUpThigh.normalize(), XinSag)) - Pi/2;

	Vec3 projUpPelvis = upPelvis - dot(upPelvis, sagN)*sagN;
	_curTrunkAngle[0] = 
		acos(dot(projUpPelvis.normalize(), XinSag)) - Pi/2;
	
	Vec3 projUpPelvisCor = upPelvis - dot(upPelvis, corN)*corN;
	_curTrunkAngle[1] = 
		acos(dot(projUpPelvisCor.normalize(), ZinCor)) - Pi/2;
	
	Vec3 projUpThighCor = upThigh - dot(upThigh, corN)*corN;
	_curSWTAngle[1] = 
		acos(dot(projUpThighCor.normalize(), ZinCor)) - Pi/2;
	

#ifdef USE_GLOBAL_ANKLE
	const MobilizedBody& rFoot = _model.foot_r;
	const MobilizedBody& lFoot = _model.foot_l;
	
	Vec3 upRFoot, upLFoot; 
	getUpVectorInGround(s, rFoot, upRFoot); 
	getUpVectorInGround(s, lFoot, upLFoot); 
	
	Vec3 projUpRFoot = upRFoot - dot(upRFoot, sagN)*sagN;
	_curRFootAngle[0] = 
		acos(dot(projUpRFoot.normalize(), XinSag)) - Pi/2;
	
	Vec3 projUpLFoot = upLFoot - dot(upLFoot, sagN)*sagN;
	_curLFootAngle[0] = 
		acos(dot(projUpLFoot.normalize(), XinSag)) - Pi/2;
	
	Vec3 projUpRFootCor = upRFoot - dot(upRFoot, corN)*corN;
	_curRFootAngle[1] = 
		acos(dot(projUpRFootCor.normalize(), ZinCor)) - Pi/2;
	
	Vec3 projUpLFootCor = upLFoot - dot(upLFoot, corN)*corN;
	_curLFootAngle[1] = 
		acos(dot(projUpLFootCor.normalize(), ZinCor)) - Pi/2;
	
#endif

#ifdef USE_GLOBAL_HIPROT
	Vec3 frontPelvis; 
	getFrontVectorInGround(s, pelvis, frontPelvis); 
	
	frontPelvis[1] = 0.0; // project to 2d
	_curPelvisRotation = acos(dot(frontPelvis.normalize(), Vec3(0.0, 0.0, 1.0))) - Pi/2;
#endif
 
}

SimbiconStateHandler::SimbiconStateHandler(Humanoid& model, Simbicon& simctrl,
    Real interval) 
        : PeriodicEventHandler(interval), _model(model), _simctrl(simctrl) {
	}

void SimbiconStateHandler::handleEvent(State& s, Real accuracy, bool& shouldTerminate) const
{
    shouldTerminate = false;
	Simbicon* simctrl = &_simctrl;

	_model.system.realize(s, Stage::Dynamics); 
    bool lContact, rContact; // total contact force magnitudes
    _model.findContactStatus(s, lContact, rContact);


#ifndef DROP_LANDING
	int curState = simctrl->getState();
	double duration = s.getTime() - simctrl->getStateStartTime();  

	if (curState < 0) {
		if (rContact) 
			simctrl->setState(2, s.getTime()); 
		else if (lContact)
			simctrl->setState(0, s.getTime()); 
	}
	else if (curState == 0) {
		if (duration > 0.3) {
			simctrl->setState(1, s.getTime()); 
//        	std::cout << "Trans 0 -> 1" << std::endl;;
		}
		else if (rContact && duration > 0.1) { 
			simctrl->setState(2, s.getTime()); 
  //      	std::cout << "Trans 0 -> 2" << std::endl;;
		}
	}
	else if (curState == 1) {
		if (rContact && duration > 0.1) { 
			simctrl->setState(2, s.getTime()); 
    //    	std::cout << "Trans 1 -> 2" << std::endl;;
		}
	}
	else if (curState == 2) {
		if (duration > 0.3) {
			simctrl->setState(3, s.getTime()); 
      //  	std::cout << "Trans 2 -> 3" << std::endl;;
		}
		else if (lContact && duration > 0.1) {
			simctrl->setState(0, s.getTime());
       // 	std::cout << "Trans 2 -> 0" << std::endl;;
		}
	}
	else if (curState == 3) {
		if (lContact && duration > 0.1) {
			simctrl->setState(0, s.getTime());
       // 	std::cout << "Trans 3 -> 0" << std::endl;;
		}
	}

	simctrl->computeSecondaryStateVals(s, 0, 0);
#endif
}

