#ifndef EVENTHANDLERS_H_
#define EVENTHANDLERS_H_

#include <OpenSim/OpenSim.h>

//==============================================================================
// EVENT HANDLER: UserInputHandler
// Check for user input. If there has been some, process it.
//==============================================================================
class UserInputHandler : public SimTK::PeriodicEventHandler {
public:
    UserInputHandler(SimTK::Visualizer::InputSilo& silo, SimTK::Real interval); 

    void handleEvent(SimTK::State& s, SimTK::Real accuracy, bool& shouldTerminate) const;
private:
    SimTK::Visualizer::InputSilo& _silo;
};



#endif // EVENTHANDLERS_H_
