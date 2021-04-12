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
#include "model.h"
#include "tau.h"
using namespace Gillespy;

// #define NV_Ith_S(v,i) (NV_DATA_S(v)[i]) // Access to individual components of data array, of N len vector

static int f(realtype t, N_Vector y, N_Vector y_dot, void *user_data); // forward declare function to be used to solve RHS of ODE

struct UserData {
  Gillespy::Simulation *my_sim;
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
	void init_species_mode(const Model &model){
		int num_species = model.number_species;
		for (int s = 0; s < num_species; s++){
			// if the user chooses discrete, initialise the partition flag to such.
			if (model.species[s].user_mode == DISCRETE){
				model.species[s].partition_mode = DISCRETE;
			}
			// otherwise, either the user chose continuous or dynamic (or null).
			// in any case, just initialise to continuous.
			else {
				model.species[s].partition_mode = CONTINUOUS;
			}
		}
	}
	void partition_species(const Model &model, const std::vector<double> &propensity_values, double tau_step, double current_time)
	{
	}
		std::pair<std::map<std::string, int>, double> get_reactions(const Gillespy::Model *model, const std::vector<double> &propensity_values, double tau_step, double current_time, double save_time)
	{
		/*
     * Helper Function to get reactions fired from t to t+tau. Affects two values:
     * rxn_count - dict with key=Reaction channel value=number of times fired
     * curr_time - float representing current time
     */

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
			TauArgs tau_args = initialize(*(simulation->model),tau_tol);
			double increment = simulation->timeline[1] - simulation->timeline[0];

			//initialize current_state vector to 0 for each species
			std::vector<int> current_state(num_species);
			//initialize propensity_values to 0 for each species
			std::vector<double> propensity_values(num_reactions);

			//copy initial state for each trajectory
			for(int s = 0; s < num_species; s++){
				simulation->trajectories[0][0][s] = species[s].initial_population;
			}
			//Simulate for each trajectory
			//make new method here
			for(int traj = 0; traj < num_trajectories; traj++){
				if (interrupted){
					break;
				}

				for (int s = 0; s < num_species; s++) {
					current_state[s] = species[s].initial_population;
				}
				simulation->current_time = 0;
				//what is this?
				int entry_count = 0;
				//propensity sum is...
				double propensity_sum;
				//save time is...
				double save_time = 0;
				// steps rejected is...
				int steps_rejected = 0;
				//tau_step is...
				double tau_step;
				
				std::vector<int> prev_curr_state;

				// while (entry_count < simulation->number_timesteps)
			}
		}
	}
}

static int f(realtype t, N_Vector y, N_Vector y_dot, void *user_data) {
  	// N_VGetArrayPointer returns a pointer to the data in the N_Vector class.
  	realtype *ydata  = N_VGetArrayPointer(y); // pointer y vector
  	realtype *dydata = N_VGetArrayPointer(y_dot); // pointer ydot vec
  	UserData *sim_data;
  	sim_data = (UserData*) user_data; // void pointer magic

    std::vector <double> curr_state; // create vector of curr_state doubles, sent to propensity_function->ODEEvaluate()
  	int number_species = sim_data->my_sim->model->number_species; // for readability
  	int number_reacs = sim_data->my_sim->model->number_reactions; // for readability
   	std::vector <realtype> propensity; // Vector of propensities of type 'realtypes' (doubles used in SUNDIALS),

  	for (sunindextype i = 0; i < number_species; i++){
		dydata[i] = 0; // Initialize change in y to '0'
  		curr_state.push_back(ydata[i]); // add values found in our curr_state vector, 'ydata' to our curr_state <double> vector
  		// This vector is used for our propensity_function method "evaluate", defined in abstract in 'model.h'
  	}

  	for (sunindextype rxn = 0; rxn < number_reacs; rxn++){
		// Calculate propensity for each reaction, at current state
  		propensity.push_back((sim_data->my_sim)->propensity_function->ODEEvaluate((int)rxn, curr_state));

  	   	for (sunindextype spec = 0; spec < (sim_data->my_sim)->model->number_species; spec++){
			// if species is a product of this reaction, add propensity fxn
    			if ((sim_data->my_sim)->model->reactions[rxn].species_change[spec] > 0){
    				dydata[spec] += propensity[rxn];
    			}
    			// else if this species is reactant, subtract propensity fxn
    			else if ((sim_data->my_sim)->model->reactions[rxn].species_change[spec] < 0){
				dydata[spec] -= propensity[rxn];
    			}
  	  	  }
  	}
  return(0);
}

