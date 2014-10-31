#include <Simbody.h>

#include "Humanoid.h"
#include "Simbicon.h"

using namespace SimTK;

namespace {
// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo, 
                     Real                   interval); 
    void handleEvent(State& state, Real accuracy,
                     bool& shouldTerminate) const OVERRIDE_11;
private:
    Visualizer::InputSilo& m_silo;
};

// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input. 
class OutputReporter : public PeriodicEventReporter {
public:
    OutputReporter(const Humanoid& model, 
                   Real            interval); 
    void handleEvent(const State& state) const OVERRIDE_11;
private:
    const Humanoid& _model;
};

// Write interesting integrator info to stdout at end of simulation.
void dumpIntegratorStats(double startCPU, double startTime,
                         const Integrator& integ);
}

//==============================================================================
// MAIN FUNCTION
//==============================================================================
int main(int argc, char **argv)
{
    try {

	double finalTime = 1000; 

    Humanoid model;
    model.system.addEventHandler
       (new UserInputHandler(*model.userInput, Real(0.1))); //check input every 100ms

    model.system.addEventReporter(new OutputReporter(model, .01));

    // Add the controller.
    Simbicon* simctrl = new Simbicon(model);
    Force::Custom simbicon(model.forces, simctrl); // takes ownership

    model.system.addEventHandler(new SimbiconStateHandler(model,*simctrl,
                                                          STATE_UPD_STEPSIZE));

    State s = model.system.realizeTopology();
    #ifdef USE_BALLS
        model.matter.setUseEulerAngles(s, true);
    #endif
    model.system.realizeModel(s);
    model.fillInActuatorMap(s);

    printf("Act: u\n");
    for (int i=0; i < NumActuators; ++i) {
        printf("%2d: %d\n", i, int(model.getUIndex(Actuator(i))));
    }

    //model.toes_r.lockAt(s, .2); model.toes_l.lockAt(s, .2);
    //model.foot_r.lockAt(s, Vec2(.2,0)); model.foot_l.lockAt(s, Vec2(.2,0));
    model.system.realize(s, Stage::Instance);

    printf("SIMBICON 3D:\n");
    printf("%d bodies, %d mobilities, -%d constraint equations -%d motions\n",
        model.matter.getNumBodies(), s.getNU(), s.getNMultipliers(), 
        model.matter.getKnownUDotIndex(s).size());

    model.trunk.setQToFitTranslation(s, Vec3(0,1.5,0));
    model.trunk.setUToFitLinearVelocity(s, Vec3(1,0,0));
    model.viz.report(s); 

    // Simulate.
    //CPodesIntegrator integ(model.system); integ.setOrderLimit(2); integ.setAccuracy(.01);
    
    //RungeKuttaMersonIntegrator integ(model.system); integ.setAccuracy(1e-3);
    //RungeKutta2Integrator integ(model.system); integ.setAccuracy(.1);
    //SemiExplicitEuler2Integrator integ(model.system); integ.setAccuracy(0.1);
    //SemiImplicitEulerIntegrator integ(model.system, .002);
    //integ.setConstraintTolerance(.001);
    //integ.setMaximumStepSize(.005);
    SemiExplicitEulerTimeStepper ts(model.system);
    //ts.setRestitutionModel(SemiExplicitEulerTimeStepper::Poisson);
    //ts.
    ts.initialize(s);
    model.viz.report(ts.getState());
    printf("Hit ENTER to simulate ... (ESC to quit)\n");
    model.userInput->waitForAnyUserInput(); model.userInput->clear();

    const double startCPU  = cpuTime(), startTime = realTime();
    do {
        model.viz.report(ts.getState());
        ts.stepTo(ts.getState().getTime() + 0.001);
    } while (ts.getTime() < 100);

    try {
        ts.stepTo(Infinity); // RUN
    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        model.system.realize(ts.getState());
        std::cout << "y=" << ts.getState().getY() << std::endl;
        std::cout << "ydot=" << ts.getState().getYDot() << std::endl;
        throw;
    }

    // dumpIntegratorStats(startCPU, startTime, integ);

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }	 
	return 0; 
}


//==============================================================================
//                           OUTPUT REPORTER
//==============================================================================
OutputReporter::OutputReporter(const Humanoid& model, 
                               Real            interval)
:   PeriodicEventReporter(interval), _model(model) {}

void OutputReporter::handleEvent(const State& state) const {
    //printf("OutputReporter @t=%g:\n", state.getTime());
    //Real fLeft, fRight;
    //_model.findContactForces(state, fLeft, fRight);
    //printf("  forces left=%g, right=%g\n", fLeft, fRight);
}

//==============================================================================
//                           USER INPUT HANDLER
//==============================================================================
UserInputHandler::UserInputHandler(Visualizer::InputSilo& silo, 
                                   Real                   interval) 
:   PeriodicEventHandler(interval), m_silo(silo) {}

void UserInputHandler::handleEvent(State& state, Real accuracy,
                                   bool& shouldTerminate) const  {
    while (m_silo.isAnyUserInput()) {
        unsigned key, modifiers;
        while (m_silo.takeKeyHit(key,modifiers))
            if (key == Visualizer::InputListener::KeyEsc) {
                shouldTerminate = true;
                m_silo.clear();
                return;
            }
    }  
}

//==============================================================================
//                        DUMP INTEGRATOR STATS
//==============================================================================
namespace {
void dumpIntegratorStats(double startCPU, double startTime, 
                         const Integrator& integ) {
    std::cout << "DONE: Simulated " << integ.getTime() << " seconds in " <<
        realTime()-startTime << " elapsed s, CPU="<< cpuTime()-startCPU << "s\n";
    #ifdef ANIMATE
    printf("***CAUTION: CPU time not accurate when animation is enabled.\n");
    #endif

    const int evals = integ.getNumRealizations();
    std::cout << "\nUsed "  << integ.getNumStepsTaken() << " steps, avg step=" 
        << (1000*integ.getTime())/integ.getNumStepsTaken() << "ms " 
        << (1000*integ.getTime())/evals << "ms/eval\n";

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n",  integ.getNumStepsTaken(), 
                                          integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n",     integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), 
                                          integ.getNumProjections());
    // Ratio of failed steps to successful ones is crude stiffness measure.
    printf("\nEstimated stiffness: %g\n",
       (double)integ.getNumErrorTestFailures()/integ.getNumStepsTaken());
}
}
