#include "PredictiveOptimizationSystem.h"
#include "Simbicon.h"
#include "EventHandlers.h"
#include "Humanoid.h"

using namespace OpenSim;

SimTK::Visualizer* myVisualizer; 

int PredictiveOptimizationSystem::objectiveFunc(const SimTK::Vector& params,
		bool reporting, SimTK::Real& f) const {
	bool useViz = true; 
	SimTK::MultibodySystem& mbs = _model._system;

	const SimTK::SimbodyMatterSubsystem& matter = _model._matter;
	const SimTK::GeneralForceSubsystem&  forces = _model._forces;

	SimTK::State s = _model._system.realizeTopology();


	double CaptureVelocity = 0.001; 
	double CoefRest = 0.; 
	double mu_d = 2.0; 	
	double mu_s = 2.0; 	
	double mu_v = 0.05; 	
	SimTK::MyUnilateralConstraintSet unis(mbs, CaptureVelocity);
	
#ifdef RIGID_CONTACT
	const MarkerSet& markers = model->getMarkerSet(); 

	for (int i = 0; i < markers.getSize(); i++) {
		SimTK::MobilizedBody& mobBody = 
			matter.updMobilizedBody(markers[i].getBody().getIndex());
		const SimTK::Vec3& pointLoc = markers[i].getOffset();
		
		SimTK::MyPointContact* contact = 
			new SimTK::MyPointContact(mobBody, pointLoc, CoefRest);
		unis.addContactElement(contact);
		unis.addFrictionElement(
				new SimTK::MyPointContactFriction(*contact, mu_d, mu_s, mu_v,
					CaptureVelocity, // TODO: vtol?
					forces));
	}

	for (int i=0; i < unis.getNumContactElements(); ++i) {
		mbs.addEventHandler(new SimTK::ContactOn(mbs, unis,i, SimTK::Stage::Position));
		mbs.addEventHandler(new SimTK::ContactOn(mbs, unis,i, SimTK::Stage::Velocity));
		mbs.addEventHandler(new SimTK::ContactOn(mbs, unis,i, SimTK::Stage::Acceleration));
		mbs.addEventHandler(new SimTK::ContactOff(mbs, unis,i));

	}

	for (int i=0; i < unis.getNumFrictionElements(); ++i) {
		mbs.addEventHandler(new SimTK::StictionOn(mbs, unis, i));
		mbs.addEventHandler(new SimTK::StictionOff(mbs, unis, i));
	}
#endif
	 
//    SimTK::ExplicitEulerIntegrator integrator(mbs);
    //SimTK::SemiExplicitEuler2Integrator integrator(mbs);
//    SimTK::SemiExplicitEulerIntegrator integrator(mbs, 0.0001);
    //SimTK::RungeKuttaMersonIntegrator integrator(mbs);
    //SimTK::CPodesIntegrator integrator(mbs);
	//integrator.setOrderLimit(2); 
	
	SimTK::SemiExplicitEulerTimeStepper ts(mbs);

            
	if (useViz) {
		SimTK::Visualizer& viz = _model._viz;
		viz.setWindowTitle("Simbicon");
		viz.setBackgroundType(viz.GroundAndSky);
		viz.setGroundHeight(0.0);
		viz.setShowShadows(true);
		viz.setShowFrameRate(true);
		viz.setShowSimTime(true);
		viz.setShowFrameNumber(false);
//		model->updMatterSubsystem().setShowDefaultGeometry(true);
       	viz.addDecorationGenerator(new SimTK::ShowContact(unis));  
		SimTK::Visualizer::InputSilo& silo = *_model._userInput;
        UserInputHandler* userInput = new UserInputHandler(silo, 0.001);
        mbs.addEventHandler(userInput); 
//#ifndef RIGID_CONTACT
//		const ContactGeometrySet& contactSet = model->getContactGeometrySet();
//		double contactRad = 
//			static_cast<ContactSphere*>(&contactSet[1])->getRadius(); 
//		model->updVisualizer().getGeometryDecorationGenerator()->setDispMarkerRadius(contactRad);
//#endif
		myVisualizer = &viz; 
	}
//	integrator.setMaximumStepSize(0.005);
//	integrator.setMinimumStepSize(1e-10);
	integrator.setAccuracy(1e-2);
//	integrator.setConstraintTolerance(1e-3);
//	integrator.setAllowInterpolation(true);
	integrator.setReturnEveryInternalStep(true);
	integrator.setFinalTime(_finalTime);
        
	SimbiconStateHandler* simbiconHandler = new SimbiconStateHandler(_model, unis, STATE_UPD_STEPSIZE);
	_model._system.addEventHandler(simbiconHandler);

	s = _model._system.realizeTopology();
	
//	std::cout << s.updQ() << std::endl; 
//	std::cout << s.updU() << std::endl; 
	// TIMESTEPPER
	ts.setReportAllSignificantStates(true);
	ts.initialize(s);       // set integrator's initial state

	const double realStart = SimTK::realTime();    // start the timer.
	const double cpuStart = SimTK::cpuTime();    // start the timer.
	bool unknownException = false;

	int numEvents = 0; 
	SimTK::Integrator::SuccessfulStepStatus status;
	/*  (1) ReachedReportTime --------- stopped only to report; state might be interpolated.
		(2) ReachedEventTrigger ------- localized an event; this is the before state (interpolated).
		(3) ReachedScheduledEvent ----- reached the limit provided in stepTo() (scheduled event).
		(4) TimeHasAdvanced ----------- user requested control whenever an internal step is successful.
		(5) ReachedStepLimit ---------- took a lot of internal steps but didn't return control yet.
		(6) EndOfSimulation ----------- termination; don't call again.
		(7) StartOfContinuousInterval - the beginning of a continuous interval: either the start of
		the simulation, or t_high after an event handler has modified
		the state.
	 */
	try {
		while (!integrator.isSimulationOver()) {
			// Integrate until a significant step has been taken
			// (because both setReportAllSignificantStates and
			// integrator.setReturnEveryInternalStep is set to true).
			// A significant step includes all step statuses above.
			status = ts.stepTo(_finalTime);

			// Get the latest state and increment counter.
			s = ts.getState();

			// Do the following only when time has advanced in the simulation:
			//    1) Update the model cache;
			//    2) Append the state to a special time-advancing state storage.
			//
			// Note: advancing time is defined by the TimeHasAdvanced || 
			//       EndOfSimulation || StartOfContinuousInterval status
			//       and the number of states here should be the same as
			//       that reported by the integrator for the simulation.
			//       Note also that the STATE_POLLING_INTERVAL in Globals.h
			//       will affect the periodicity of the states reported
			//       due to StartOfContinuousInterval.
/*			if (status == SimTK::Integrator::TimeHasAdvanced 
					|| status == SimTK::Integrator::EndOfSimulation
					|| status == SimTK::Integrator::StartOfContinuousInterval) {
				std::cout << "Time advanced at t=" << s.getTime() << std::endl;
			}*/
			if (status == SimTK::Integrator::ReachedEventTrigger) 
				numEvents++; 
		}   
	}
	catch (const std::exception& e) {
		std::cout << "An exception occured during the simulation: \n\t" << e.what() << std::endl;
	}
	catch (...) {
		std::cout << "An exception occured during the simulation: \n\t" << std::endl;
	}

    double human_time = SimTK::realTime() - realStart;   
    double cpu_time = SimTK::cpuTime() - cpuStart;   
    const double sim_time = integrator.getTime();
    const double avrFunctEvalsPerstep = (double)integrator.getNumRealizations()/(double)integrator.getNumStepsTaken();

    std::cout << "Simulation (" << sim_time << " sec) completed in " << human_time << " sec." << std::endl;
    std::cout << "human_time/sim_time ratio = " << human_time/sim_time << "   (ratio < 1 indicates realtime)." << std::endl;

    std::cout << "Integrated using " << integrator.getMethodName() << " at accuracy " << integrator.getAccuracyInUse() << std::endl;
    std::cout << integrator.getNumStepsTaken() << " integration steps taken (" << integrator.getNumStepsTaken()/sim_time << " steps/sec)" << std::endl;

    std::cout << "Average step size: " << (1000*sim_time)/integrator.getNumStepsTaken() << " ms." << std::endl;
    std::cout << "Average step computation time: " << (1000*cpu_time)/integrator.getNumStepsTaken() << " ms." << std::endl;
    std::cout << "Average function evaluations (realizations) per step: " << avrFunctEvalsPerstep << std::endl;

    std::cout << "# STEPS/ATTEMPTS = "  << integrator.getNumStepsTaken() << " / " << integrator.getNumStepsAttempted() << std::endl;
    std::cout << "# ERR TEST FAILS = "  << integrator.getNumErrorTestFailures() << std::endl;
    std::cout << "# REALIZE/PROJECT = " << integrator.getNumRealizations() << " / " << integrator.getNumProjections() << std::endl;

	std::cout << "# Events triggered = " << numEvents << std::endl;
	f = 0; 
	return 0; 
}
