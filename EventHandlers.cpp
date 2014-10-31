#include "EventHandlers.h"
#include "Simbicon.h"

using namespace std;
using namespace SimTK;
using namespace OpenSim;

extern Visualizer* myVisualizer;  

UserInputHandler::UserInputHandler(
    Visualizer::InputSilo& silo,
    Real interval) 
        : PeriodicEventHandler(interval), _silo(silo) {
	}

void UserInputHandler::handleEvent(State& s, Real accuracy, bool& shouldTerminate) const
{
    shouldTerminate = false;
    while (_silo.isAnyUserInput()) {
        unsigned key, modifiers;
            
        while (_silo.takeKeyHit(key, modifiers)) {
            // End simulation: ESC
            if (key == Visualizer::InputListener::KeyEsc) {
                shouldTerminate = true;
                _silo.clear();
                return;
            }
            // Pause/unpause simulation: SPACE
            else if (key == ' ') {
				myVisualizer->report(s); 
                _silo.waitForAnyUserInput();
                return;
            }
        }
    }
}

