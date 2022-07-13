/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2022 GillesPy2 developers.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include "cvode.h"
#include "cvode_spils.h"
#include "sundials_math.h"
#include "sundials_types.h"
#include "nvector_serial.h"
#include "sunlinsol_spgmr.h"

#include "ODESolver.h"

/* Used to access individual components of a length N vector. */
#define NV_Ith_s(v, i) (NV_DATA_S(v)[i])

namespace Gillespy
{
    static volatile bool interrupted = false;

    GPY_INTERRUPT_HANDLER(signal_handler, {
        interrupted = true;
    })

    static int f(realtype t, N_Vector y, N_Vector y_dot, void *user_data);

    struct UserData
    {
        Simulation<double> *my_sim;
    };

    bool cmpf(float A, float B){
        return (fabs(A - B) < FLT_EPSILON);
    }

    void ODESolver(Simulation<double> *simulation, double increment, SolverConfiguration config)
    {
        GPY_INTERRUPT_INSTALL_HANDLER(signal_handler);

        // CVODE constants are returned on every success or failure.
        // CV_SUCCESS: Operation was successful.
        // CV_MEM_NULL: CVODE memory block was not initialized with CVodeCreate.
        // CV_NO_MALLOC: The allocation function CVodeInit was not called.
        // CV_ILL_INPUT: An input tolerance was negative.
        int flag;

        // Allocate memory for data passed into RHS.
        UserData *data = new UserData();
        data->my_sim = simulation;

        // Init the initial conditions.
        sunindextype N = (simulation->model)->number_species;
        N_Vector y0 = N_VNew_Serial(N);

        Model<double> *simulation_model = simulation->model;

        // Add species initial conditions to the current state vectory `y0`.
        for (unsigned int species_index = 0; species_index < simulation_model->number_species; species_index++)
        {
            double initial_population = simulation_model->species[species_index].initial_population;

            NV_Ith_S(y0, species_index) = initial_population;
            simulation->current_state[species_index] = initial_population;
        }
        simulation->output_buffer_range(std::cout);

        // Create and set CVODE object pointer. CV_ADAMS for nonstiff, CV_BDF for stiff.
        void *cvode_mem = CVodeCreate(CV_BDF);
        realtype t0 = 0;

        // Initialize the ODE solver and set tolerances.
        flag = CVodeInit(cvode_mem, f, t0, y0);
        flag = CVodeSStolerances(cvode_mem, config.rel_tol, config.abs_tol);
        switch (CVodeSetMaxStep(cvode_mem, config.max_step))
        {
            case CV_ILL_INPUT:
                std::cerr << "Bad step size: " << config.max_step << std::endl;
                break;
        }

        // Initialize and select the linear solver module.
        // SUNSPMR: Iterative Solver (compatible with serial, threaded, parallel, and user suppoed NVector).
        // SUNLinSol_SPGMR(N_Vector y, int pretype, intm axl)
        // - N_Vector y: Vector to be used in solver.
        // - int pretype: Flag indicating desired precondition type. `0` = None.
        // - int maxl: The number of Krylov basis vectors to use. Values <= 0 default to `5`.
        SUNLinearSolver linear_solver = SUNLinSol_SPGMR(y0, 0, 0);

        // Attach linear solver module.
        flag = CVodeSetUserData(cvode_mem, data);

        // CVodeSetLinearSolver(void *cvoid_mem, SUNLinearSolver LS, SUNMatrix A)
        // - void *cvode_mem: Pointer to CVODE memory block.
        // - SUNLinearSolver LS: SUNLINSOL object to use for solving linear systems.
        // - SUNMaxtrix A: Object to use as a template for Jacobian, defaults to NULL if not applicable.
        flag = CVodeSetLinearSolver(cvode_mem, linear_solver, NULL);

        // For each point at which output is desired, call `ier = CVode(cvode_mem, tout, yout, &tret, itask)`.
        // Here, itask specifies the return mode.
        // The vector `yout` (which can be the same as y0) will contain `y(t)`.

        // The next time at which the computed solution is desired.
        realtype tout;
        realtype end_time = simulation->end_time;
        realtype step_length = increment;

        for (tout = step_length; !interrupted && tout < end_time || cmpf(tout, end_time); tout += step_length)
        {
            // CV_NORMAL causes the solver to take internal steps until it has reached or just passed the `tout`
            // parameter. The solver interpolates in order to return an approximate value of `y(tout)`.
            // CVode() returns a vector `y0` (or `y(tout)`), and corresponding variable value `t` = `tret` (return time).
            // With CV_NORMAL `tret` is equal to `tout` and `y0` = `y(tout)`.
            flag = CVode(cvode_mem, tout, y0, &simulation->current_time, CV_NORMAL);

            for (sunindextype species = 0; species < N; species++)
            {
                simulation->current_state[(int)species] = NV_Ith_S(y0, species);
            }
            simulation->output_buffer_range(std::cout);
        }

        // Deallocate the solution vector.
        N_VDestroy(y0);

        // Deallocate solver memory.
        CVodeFree(&cvode_mem);

        // Deallocatae linear solver / matrix memory.
        SUNLinSolFree(linear_solver);
    }

    static int f(realtype t, N_Vector y, N_Vector y_dot, void *user_data)
    {
        UserData *sim_data = (UserData *)user_data;
        Simulation<double> *simulation = sim_data->my_sim;
        Model<double> *model = simulation->model;

        // N_VGetArrayPointer returns a pointer to the data property within N_Vector.
        realtype *ydata = N_VGetArrayPointer(y);
        realtype *dydata = N_VGetArrayPointer(y_dot);

        int number_species = model->number_species;
        int number_reactions = model->number_reactions;

        std::vector<double> current_state;
        std::vector<realtype> propensity;

        for (sunindextype species_index = 0; species_index < number_species; species_index++)
        {
            dydata[species_index] = 0;
            current_state.push_back(ydata[species_index]);
        }

        for (sunindextype reaction_index = 0; reaction_index < number_reactions; reaction_index++)
        {
            // Calculate propensity for each reaction at the current state.
            propensity.push_back(Reaction::propensity(reaction_index, ydata));

            for (sunindextype species_index = 0; species_index < number_species; species_index++)
            {
                // If the species is a product (positive) or reactant (negativ) of this reaction,
                // add the propensity function multiplied by the species_change value (e.g -2 for 
                // a bi-molecular reaction reactants). 
                if (model->reactions[reaction_index].species_change[species_index] != 0)
                {
                    dydata[species_index] += propensity[reaction_index] * model->reactions[reaction_index].species_change[species_index];
                }
            }
        }

        return 0;
    }
}
