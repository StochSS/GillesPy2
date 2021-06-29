#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <time.h>
#include <math.h>
#include "TauHybridSolver.h"
#include "template.h"
#include "hybrid_template.h"
#include "arg_parser.h"
using namespace Gillespy;

//Default values, replaced with command line args
unsigned int number_trajectories = 0;
unsigned int number_timesteps = 0;
int random_seed = 0;
double end_time = 100.0;
bool seed_time = true;
double increment = 0;
double tau_tol = 0.05;
class PropensityFunction : public IPropensityFunction{
public:

	double ODEEvaluate(int reaction_number, const std::vector <double> &S){
		return map_ode_propensity(reaction_number, S);
	}
    double TauEvaluate(unsigned int reaction_number, const std::vector<int> &S) {
		return map_propensity(reaction_number, S);
	}
    double evaluate(unsigned int reaction_number, unsigned int* S){return 1.0;}
};

double Gillespy::TauHybrid::HybridReaction::ode_propensity(
	ReactionId reaction_number,
	std::vector<double> &state)
{
	return map_ode_propensity(reaction_number, state);
}

double Gillespy::TauHybrid::HybridReaction::ssa_propensity(
	ReactionId reaction_number,
	std::vector<int> &state)
{
	return map_propensity(reaction_number, state);
}

int main(int argc, char* argv[]){
	ArgParser parser(argc, argv);
	random_seed = parser.seed;
	seed_time = (random_seed == -1);

	end_time = parser.end;
	number_trajectories = parser.trajectories;
	number_timesteps = parser.timesteps;
	tau_tol = parser.tau_tol;

	Model model(species_names, species_populations, reaction_names);
	add_reactions(model);

	if(seed_time){
		random_seed = time(NULL);
	}
	IPropensityFunction *propFun = new PropensityFunction();
	//Simulation INIT
	TauHybrid::HybridSimulation simulation(model);
	simulation.model = &model;
	simulation.end_time = end_time;
	simulation.random_seed = random_seed;
	simulation.number_timesteps = number_timesteps;
	simulation.number_trajectories = number_trajectories;
	simulation.propensity_function = propFun;
	init_simulation(&model, simulation);
	Gillespy::TauHybrid::map_species_modes(simulation.species_state);
	// Perform ODE  //
	TauHybrid::TauHybridCSolver(&simulation, tau_tol);
	simulation.output_results_buffer(std::cout);
	delete propFun;
	return 0;
}
