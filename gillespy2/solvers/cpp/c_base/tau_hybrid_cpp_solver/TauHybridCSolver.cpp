#include <iostream>
#include <csignal> //Included for timeout signal handling
#include <random>
#include "cvode.h" // prototypes for CVODE fcts., consts.
#include "nvector_serial.h"  // access to serial N_Vector
#include "sunlinsol_spgmr.h"  //access to SPGMR SUNLinearSolver
#include "cvode_spils.h" // access to CVSpils interface
#include "sundials_types.h"  // defs. of realtype, sunindextype
#include "sundials_math.h"  // contains the macros ABS, SUNSQR, EXP
#include "TauHybridCSolver.h"
#include "HybridModel.h"
#include "tau.h"
using namespace Gillespy;

// #define NV_Ith_S(v,i) (NV_DATA_S(v)[i]) // Access to individual components of data array, of N len vector

static int f(realtype t, N_Vector y, N_Vector y_dot, void *user_data); // forward declare function to be used to solve RHS of ODE

struct UserData {
  Gillespy::Simulation *my_sim;
  std::vector<double> reaction_states;
};
struct IntegratorOptions{
  // CVODE constants returned if bad output, or success output.
  // Constants: CV_SUCCESS,
  // CV_MEM_NULL: CVODE memory block not initialized through call to CVodeCreate
  // CV_NO_MALLOC: The allocation function CVodeInit not called
  // CV_ILL_Input: An input tolerance was negative
  int flag;
  // absolute tolerace of a system
  realtype abstol;
  // relative tolerance of system
  realtype reltol;
  // double max_step;
};
namespace Gillespy {
	bool interrupted = false;

	void signalHandler(int signum)
	{
		interrupted = true;
	}
	/**************************************
	 ***** POLICE TAPE : DO NOT CROSS *****
	 **************************************//*
	void init_species_mode(const Model &model, Simulation &simulation){
		int num_species = model.number_species;
		// int num_det_species = 0;
		for (int s = 0; s < num_species; s++){
			// if the user chooses discrete, initialise the partition flag to such.
			if (model.species[s].user_mode == DISCRETE){
				model.species[s].partition_mode = DISCRETE;
			}
			// otherwise, either the user chose continuous or dynamic (or null).
			// in any case, just initialise to continuous.
			else {
				model.species[s].partition_mode = CONTINUOUS;
				// num_det_species++;
			}
		}
		// simulation.number_det_species = num_det_species;
	}
	void partition_species(const Model &model, const std::vector<double> &propensity_values, std::vector<hybrid_state> curr_state, double tau_step, double current_time, std::map<int, bool> &det_species){
		// coefficient of variance- key:species id, value: cv
		std::map<int, double> cv;
		// means
		std::map<int, double> means;
		// standard deviation
		std::map<int, double> sd;
		//init means
		for (int i = 0; i < model.number_species; ++i){
			if (model.species[i].user_mode == DYNAMIC){
				if (model.species[i].partition_mode == CONTINUOUS){
					means.insert({i, curr_state[i].continuous});
				}else {
					means.insert({i, curr_state[i].discrete});
				}
			}
		}
		//init sd's
		for (int i = 0; i < model.number_species; ++i){
			if (model.species[i].user_mode == DYNAMIC){
				sd.insert({i, 0});
			}
		}
		// calculate means and standard deviations for dynamic-mode species involved in reactions
		for (int r = 0; r < model.number_reactions; ++r){
			for (int s = 0; s < model.number_species; ++s){
				// access list of species by accessing the correct element of the state-change vector (of the reaction)
				if (model.species[model.reactions[r].species_change[s]].user_mode == DYNAMIC ){
					// if less than 0, that means this is a reactant
					if (model.reactions[r].species_change[s] < 0){
						means[s] -= (tau_step * propensity_values[r] * model.reactions[r].species_change[s]);
						sd[s] += std::pow((tau_step * propensity_values[r] * model.reactions[r].species_change[s]),2);
					}
					// if greater than 0, that means this is a product
					if (model.reactions[r].species_change[s] > 0){
						means[s] += (tau_step * propensity_values[r] * model.reactions[r].species_change[s]);
						sd[s] += std::pow((tau_step * propensity_values[r] * model.reactions[r].species_change[s]),2);
					}
				}
			}
		}
		// calculate coefficient of variation using means and sd
		
		for (int s = 0; s < model.number_species; ++s){
			if (means.count(s) > 0){
				Species sref = model.species[s];
				if (sref.switch_min == 0) { // (default value means switch min not set, use switch tol)
					if (means[s] > 0){
						cv[s] = sd[s] / means[s];
					}else{
						cv[s] = 1;
					}
					if (cv[s] < sref.switch_tol){
						sref.partition_mode == CONTINUOUS;
					}else{
						sref.partition_mode == DISCRETE;
					}
				}else{
					if (means[s] > sref.switch_min){
						sref.partition_mode == CONTINUOUS;
					}else{
						sref.partition_mode == DISCRETE;
					}
				}
			}
		}
		return //TODO;
	}
		std::pair<std::map<std::string, int>, double> get_reactions(const Gillespy::Model *model, const std::vector<double> &propensity_values, double tau_step, double current_time, double save_time)
	{
		/*
     * Helper Function to get reactions fired from t to t+tau. Affects two values:
     * rxn_count - dict with key=Reaction channel value=number of times fired
     * curr_time - float representing current time
     *//*

		if (current_time + tau_step > save_time)
			tau_step = save_time - current_time;

		std::map<std::string, int> rxn_count; // map of how many times reaction is fired
		std::random_device rd;
		std::mt19937 generator(rd());
		std::pair<std::map<std::string, int>, double> values; // value pair to be returned, map of times {map of times reaction fired, current time}

		for (int i = 0; i < model->number_reactions; i++)
		{
			std::poisson_distribution<int> poisson(propensity_values[i] * tau_step);
			rxn_count[model->reactions[i].name] = poisson(generator);
		}
		current_time = current_time + tau_step;
		values.first = rxn_count;
		values.second = current_time;
		return values;
	}
	/******************************
	 ***** END OF POLICE TAPE *****
	 ******************************/

	

	void TauHybridCSolver(Gillespy::Simulation *simulation, const double tau_tol)
	{
		//timeouts not supported right now??
		signal(SIGINT, signalHandler);
		if (simulation) {
			int num_species = (simulation->model)->number_species;
			int num_reactions = (simulation->model)->number_reactions;
			int num_trajectories = simulation->number_trajectories;
			Model &model = *(simulation->model);
			std::unique_ptr<Species[]> &species = model.species;
			//TauArgs tau_args = initialize(*(simulation->model),tau_tol);
			double increment = simulation->timeline[1] - simulation->timeline[0];


			// Population/concentration state values for each species.
			// TODO: change back double -> hybrid_state, once we figure out how that works
			std::vector<double> current_state(num_species);
			//initialize propensity_values to 0 for each species
			std::vector<double> propensity_values(num_reactions);

			// Hybrid solver is highly dependent on random numbers.
			// In order to do this, a URN on the range [0,1) is generated.
			// log( uniform(rng) ) returns a real number on the range (-inf, 0).
			// TODO: either assign a seed or set seed to be configurable
			std::mt19937_64 rng;
			std::uniform_real_distribution<double> uniform(0, 1);

			// Represents the current "randomized state" for each reaction, used as a
			//   helper value to determine if/how many stochastic reactions fire.
			// This gets initialized to a random negative offset, and gets "less negative"
			//   during the integration step.
			// After each integration step, the reaction_state is used to count stochastic reactions.
			std::vector<double> reaction_state(num_reactions);
			for (int rxn_state = 0; rxn_state < num_reactions; ++rxn_state) {
				reaction_state[rxn_state] = uniform(rng);
			}

			//copy initial state for each trajectory
			for(int s = 0; s < num_species; s++){
				simulation->trajectoriesODE[0][0][s] = species[s].initial_population;
				current_state[s] = species[s].initial_population;
			}
			//Simulate for each trajectory
			for(int traj = 0; traj < num_trajectories; traj++){
				if (interrupted){
					break;
				}

				// Initialize the species population for the trajectory.
				for (int spec_i = 0; spec_i < num_species; ++spec_i) {
					current_state[spec_i] = species[spec_i].initial_population;
				}

				// Struct acts as container for simulation state.
				// This gets passed in to the integrator.
				UserData *data = new UserData {
					simulation,
					reaction_state
				};

				// Initialize integrator state.
				N_Vector y0 = N_VNew_Serial(num_species);
				for (int spec_i = 0; spec_i < num_species; ++spec_i) {
					NV_Ith_S(y0, spec_i) = species[spec_i].initial_population;
				}

				// Build the ODE memory object and initialize it.
				// Accepts initial integrator state y0, start time t0, and RHS f.
				void *cvode_mem = CVodeCreate(CV_BDF);
				realtype t0 = 0;
				int flag = 0;
				flag = CVodeInit(cvode_mem, f, t0, y0);
				flag = CVodeSStolerances(cvode_mem, GPY_HYBRID_RELTOL, GPY_HYBRID_ABSTOL);

				// Build the Linear Solver object and initialize it.
				SUNLinearSolver LS = SUNLinSol_SPGMR(y0, 0, 0);
				flag = CVodeSetUserData(cvode_mem, data);
				flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);

				// SIMULATION STEP LOOP
				double next_time;
				double tau_step = increment;
				int save_time = 0;
				while (simulation->current_time < simulation->end_time) {
					// Determine what the next time point is.
					// This will become current_time on the next iteration.
					// If a retry with a smaller tau_step is deemed necessary, this will change.
					next_time = simulation->current_time + tau_step;

					// Integration Step
					// For deterministic reactions, the concentrations are updated directly.
					// For stochastic reactions, integration updates the reaction_states vector.
					flag = CVode(cvode_mem, next_time, y0, &next_time, CV_NORMAL);

					// The newly-updated reaction_states vector may need to be reconciled now.
					// A positive reaction_state means reactions have potentially fired.
					// NOTE: it is possible for a population to swing negative, where a smaller Tau is needed.
					for (int rxn_i = 0; false && rxn_i < num_reactions; ++rxn_i) {
						// Temporary variable for the reaction's state.
						// Does not get updated unless the changes are deemed valid.
						double rxn_state = reaction_state[rxn_i];

						// Temporary array to store changes to dependent species, 0-initialized.
						unsigned int population_changes[num_species];
						for (int p_i = 0; p_i < num_species; ++p_i)
							population_changes[p_i] = 0;

						// Use the current reaction_state to count the number of firings.
						// If a negative population is detected, then the loop breaks prematurely.
						while (rxn_state > 0) {
							// "Fire" a reaction by recording changes in dependent species.
							// If a negative value is detected, break without saving changes.
							for (int spec_i = 0; spec_i < num_species; ++spec_i) {
								population_changes[spec_i] += model.reactions[rxn_i].species_change[spec_i];
								if (population_changes[spec_i] < 0) {
									break;
								}
							}

							// uniform(rng) is a random number on range (0,1), always fractional
							// This means that log(uniform(rng)) is always negative
							rxn_state += log(uniform(rng));
						}

						// Positive reaction state means a negative population was detected.
						// Only update state with the given population changes if valid.
						if (rxn_state <= 0) {
							for (int p_i = 0; p_i < num_species; ++p_i) {
								current_state[p_i] += population_changes[p_i];
							}
							reaction_state[rxn_i] = rxn_state;
						}
						else {
							// Invalid population state detected; try a smaller Tau step.
							next_time = simulation->current_time;
							tau_step *= 0.5;

							// TODO: Reset the integrator state to the previous time step.
						}
					}

					// Output the results for this time step.
					simulation->current_time = next_time;
					
					while (save_time <= next_time) {
						// Write each species, one at a time (from ODE solution)
						for (int spec_i = 0; spec_i < num_species; ++spec_i) {
							simulation->trajectoriesODE[traj][save_time][spec_i] = NV_Ith_S(y0, spec_i);
						}
						save_time += increment;
					}
					
				}

				// End of trajectory
				// Clean up integrator data structures
				N_VDestroy_Serial(y0);
				CVodeFree(&cvode_mem);
				SUNLinSolFree_SPGMR(LS);
				delete data;
			}
		}
	}
}

/**
 * Integrator function for ODE linear solver.
 * This gets passed directly to the Sundials ODE solver once initialized.
 */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	// Get y(t) vector and f(t, y) vector
	realtype *Y = N_VGetArrayPointer(y);
	realtype *dydt = N_VGetArrayPointer(ydot);
	realtype propensity;

	// Extract simulation data
	UserData *data = static_cast<UserData*>(user_data);
	Simulation *sim = data->my_sim;
	std::vector<double> reaction_states = data->reaction_states;
	std::vector<double> concentrations(sim->model->number_species);

	// Populate the current ODE state into the concentrations vector.
	// dy/dt results are initialized to zero, and become the change in propensity.
	unsigned int spec_i;
	for (spec_i = 0; spec_i < sim->model->number_species; ++spec_i) {
		concentrations[spec_i] = Y[spec_i];
		dydt[spec_i] = 0;
	}

	// Each species has a "spot" in the y and f(y,t) vector.
	// For each species, place the result of f(y,t) into dydt vector.
	unsigned int rxn_i;
	int species_change;
	for (rxn_i = 0; rxn_i < sim->model->number_reactions; ++rxn_i) {
		propensity = sim->propensity_function->ODEEvaluate(rxn_i, concentrations);

		for (spec_i = 0; spec_i < sim->model->number_species; ++spec_i) {
			// Use the evaluated propensity to update the concentration levels and reaction state.
			// Propensity is treated as positive if it's a product, negative if it's a reactant.
			species_change = sim->model->reactions[rxn_i].species_change[spec_i];
			if (species_change == 0)
				continue;

			// The product on the right evaluates to 1 if species_change is positive,
			//    and -1 if it's negative.
			// This is a branchless alternative to using an if-statement.
			// Saw a performance gain of ~20% by using this branchless method.
			dydt[spec_i] += propensity * (-1 + 2 * (species_change > 0));
		}
	}

	return 0;
};
