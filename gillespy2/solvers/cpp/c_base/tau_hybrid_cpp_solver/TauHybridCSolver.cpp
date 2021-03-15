#include <iostream>
#include "cvode.h" // prototypes for CVODE fcts., consts.
#include "nvector_serial.h"  // access to serial N_Vector
#include "sunlinsol_spgmr.h"  //access to SPGMR SUNLinearSolver
#include "cvode_spils.h" // access to CVSpils interface
#include "sundials_types.h"  // defs. of realtype, sunindextype
#include "sundials_math.h"  // contains the macros ABS, SUNSQR, EXP
#include "TauHybridCSolver.h"
#include "model.h"
using namespace Gillespy;

#define NV_Ith_S(v,i) (NV_DATA_S(v)[i]) // Access to individual components of data array, of N len vector

static int f(realtype t, N_Vector y, N_Vector y_dot, void *user_data); // forward declare function to be used to solve RHS of ODE

struct UserData {
  Gillespy::Simulation *my_sim;
};

void TauHybridCSolver(Gillespy::Simulation* simulation, double increment){
	
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
