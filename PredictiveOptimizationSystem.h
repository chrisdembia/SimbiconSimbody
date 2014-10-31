#include <OpenSim/OpenSim.h>
class Humanoid;

namespace OpenSim {
class PredictiveOptimizationSystem : public SimTK::OptimizerSystem {
public: 
	PredictiveOptimizationSystem(int numParams, Humanoid& model, 
		double finalTime ) : 
		SimTK::OptimizerSystem(numParams), 
		_model(model), 
		_finalTime(finalTime)
		 {} 
	
	~PredictiveOptimizationSystem() {}
	
	int objectiveFunc(const SimTK::Vector& params,
			bool reporting, SimTK::Real& f) const OVERRIDE_11;  

private:
	Humanoid& _model;
	double    _finalTime;  
}; 
}
