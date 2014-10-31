#ifndef SIMBICON_H_
#define SIMBICON_H_

#include <Simbody.h>

class Humanoid;

#define STATE_UPD_STEPSIZE 0.005
//#define RIGID_CONTACT
 
class Simbicon : public SimTK::Force::Custom::Implementation
{
public:
	Simbicon(Humanoid& model); 
    void calcForce(const SimTK::State&                state, 
                   SimTK::Vector_<SimTK::SpatialVec>& bodyForces, 
                   SimTK::Vector_<SimTK::Vec3>&       particleForces, 
                   SimTK::Vector&                     mobilityForces) const 
                   OVERRIDE_11;

    SimTK::Real calcPotentialEnergy(const SimTK::State& state) const 
        OVERRIDE_11 {
        return 0;
    }

	int getState() const {
		return _state; 
	}
	void setState( int state, double time ) {
		_state = state; 
		_stateStartTime = time;
		if (state == 0 || state == 2) {
			// reset the swing thigh orientation when the swing leg changes
			for (int i = 0; i < 2; i++) {
				_lastSWTAngle[i] = -100.0;  
				_curSWTAngle[i] = -100.0;  
			}
		}
	}

	double getStateStartTime() const {
		return _stateStartTime; 
	}

	void computeSecondaryStateVals(const SimTK::State& s,
        SimTK::Real lForce, SimTK::Real rForce); 
	
private:
    void computeControls(const SimTK::State& s, SimTK::Vector& controls) const;

	void getSagCorNormals( const SimTK::State& s, 
		SimTK::Vec3& sagN, SimTK::Vec3& corN ) const; 
	void getUpVectorInGround( const SimTK::State& s, const SimTK::MobilizedBody& b, 
		SimTK::Vec3& up ) const; 
	void getFrontVectorInGround( const SimTK::State& s, const SimTK::MobilizedBody& b, 
		SimTK::Vec3& up ) const; 
	void fillInHipJointControls( const SimTK::State& s, 
		SimTK::Vector& controls ) const; 
	
    Humanoid& _model;
	int    _state;
	double _stateStartTime; 
	double _lastSWTAngle[2];  
	double _curSWTAngle[2];  
	double _lastTrunkAngle[2];  
	double _curTrunkAngle[2]; 
	double _lastRFootAngle[2];  
	double _curRFootAngle[2]; 
	double _lastLFootAngle[2];  
	double _curLFootAngle[2]; 
	double _lastPelvisRotation;  
	double _curPelvisRotation;  
	
};	

class SimbiconStateHandler : public SimTK::PeriodicEventHandler {
public:
    SimbiconStateHandler(Humanoid& m, Simbicon& simctrl, SimTK::Real interval);

    void handleEvent(SimTK::State& s, SimTK::Real accuracy, bool& shouldTerminate) const;
private:
    Humanoid& _model;
    Simbicon& _simctrl;
};

#endif // SIMBICON_H_


