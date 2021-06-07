#include <iostream>
#include "cvode.h" // prototypes for CVODE fcts., consts.
#include "nvector_serial.h"  // access to serial N_Vector
#include "sunlinsol_spgmr.h"  //access to SPGMR SUNLinearSolver
#include "cvode_spils.h" // access to CVSpils interface
#include "sundials_types.h"  // defs. of realtype, sunindextype
#include "sundials_math.h"  // contains the macros ABS, SUNSQR, EXP
#include "ODECSolver.h"
#include "model.h"
using namespace Gillespy;

#define NV_Ith_S(v,i) (NV_DATA_S(v)[i]) // Access to individual components of data array, of N len vector

static int f(realtype t, N_Vector y, N_Vector y_dot, void *user_data); // forward declare function to be used to solve RHS of ODE

struct UserData {
  Gillespy::Simulation *my_sim;
};

void ODESolver(Gillespy::Simulation* simulation, double increment){
	int flag; // CVODE constants returned if bad output, or success output.
	// Constants: CV_SUCCESS,
	// CV_MEM_NULL: CVODE memory block not initialized through call to CVodeCreate
	// CV_NO_MALLOC: The allocation function CVodeInit not called
	// CV_ILL_Input: An input tolerance was negative


	// Allocate memory for data to be passed to RHS of ODE
    UserData *data = new UserData();
    data->my_sim = simulation;
    // my_sim points to a Gillespy::Simulation struct

	realtype abstol = 1e-5; // real tolerance of system
  	realtype reltol = 1e-5; // absolute tolerance of system

	// Initial conditions
	sunindextype N = (simulation -> model)->number_species; // length of problem, 'sunindextype' index's sundials N
	//  N_VECTOR is a custom vector type with various methods. API is located on Chapter 6 of CVode guide

	N_Vector y0; // Initialize initial condition vector as an N_Vector.
	y0 = N_VNew_Serial(N);

	for(unsigned int species_number = 0; species_number < ((simulation -> model) -> number_species); species_number++){
		NV_Ith_S(y0, species_number) = (simulation -> model) -> species[species_number].initial_population;
		simulation -> trajectoriesODE[0][0][species_number] = (simulation -> model) -> species[species_number].initial_population;
	} // Add species initial conditions to 'y0', our "current state vector"
	//Initialize CVODE solver object
	void* cvode_mem = NULL; // create cvode object ptr
	cvode_mem = CVodeCreate(CV_BDF); // CV_ADAMS for nonstiff, CV_BDF for stiff problems

	realtype t0 = 0; // time i.c, must use realtype
	flag = CVodeInit(cvode_mem, f, t0, y0); // Initalize ODE Solver with allocated memory, RHS function, t0, and y0.
	// Set tolerances defined in beginning
	flag = CVodeSStolerances(cvode_mem, reltol, abstol);

	// Create solver object
	SUNLinearSolver LS;
	// Choose linear solver module
	// SUNSPMR - Iterative Solver (compatible with serial, threadel, parallel, user supplied nvector)
	// SunLinearSolver_SUNSPGMR(N_Vector y, int pretype, intm axl)
	// N_Vector y = vector to be used in solver | int pretype = flag indicating desired precondition type. '0' = none
	// int maxl =  the number of Krylov basis vectors to use. Values <= 0 defaults to '5'
	LS = SUNLinSol_SPGMR(y0, 0, 0);

	// Attach linear solver module 
   	flag = CVodeSetUserData(cvode_mem, data);
    	// CVodeSetLinearSolver(cvode_mem, LS, J))
	// cvode_mem : pointer to CVODE memory block | LS : SUNLINSOL object to use for solving linear systems
	// J : SUNMATRIX object as template for the Jacobian, default NULL if not applicable
	flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);

	// For each point at which output is desired, call
	// ier = CVode(cvode_mem, tout, yout, &tret, itask). Here, itask specifies the return mode. The vector yout
	// (which can be the same as the vector y0 above) will contain y(t). More details at 4.5.7
	realtype tout; // the next time at which computed solution is desired
	realtype end_time = simulation->end_time;
	realtype step_length = increment;
	realtype tret = 0; // the time reached by the solver for output

	// realtype yout (in this program, called y0) : the computed solution vector
	// int itask : flag indicating the job of the solver for the NEXT user step
	// CV_NORMAL option causes the solver to take internal steps until it has reached or just passed the 'tout'
	// parameter. The solver interpolates in order to return an approximate value of y(tout).
	// Cvode() returns a vector, 'y0' (which is y(tout)), and corresponding  variable value 't' = tret (return time).
	// With CV_NORMAL, tret will be equal to tout, and y0 = y(tout)
	int curr_time = 0;
	for (tout = step_length; tout <= end_time; tout += step_length){
		flag = CVode(cvode_mem, tout, y0, &tret, CV_NORMAL);
		curr_time+=1;
		for (sunindextype species = 0; species < N; species++){
			simulation->trajectoriesODE[0][curr_time][(int)species] = NV_Ith_S(y0,species);
        	}
	}

	// Deallocate memory from solution vec
	N_VDestroy(y0);

	// Free solver mem
	CVodeFree(&cvode_mem);

	// Free linear solver/matrix mem
	SUNLinSolFree(LS);
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
