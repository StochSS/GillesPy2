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

#pragma once

#include "HybridModel.h"
#include "cvode.h"
#include "sunlinsol_spgmr.h"
#include "sundials_types.h"
#include "nvector_serial.h"
#include <vector>
#include <random>

namespace Gillespy
{
    namespace TauHybrid
    {

    /* IntegratorStatus: represents the runtime state of the integrator.
     * OK indicates that no errors have occurred.
     */
    enum IntegrationStatus
    {
        // No errors have occurred.
        OK = 0,
        // Attempted to perform a SUNDIALS operation on a null CVODE object.
        NULL_POINTER,
        // A non-null object resulted in a memory error and must be initialized.
        BAD_MEMORY,
        // Could not perform integration, step size too small.
        BAD_STEP_SIZE
    };

    struct IntegratorData
    {
        HybridSimulation *simulation;
        std::vector<HybridSpecies> *species_state;
        std::vector<HybridReaction> *reaction_state;
        std::vector<Event> *events = nullptr;
        std::vector<std::function<double(double, const double*)>> active_triggers;
        // Container representing the rootfinder-enabled reactions.
        // Each integer at index i represents the reaction id corresponding to rootfinder element i.
        // In `rootfn`, this means that gout[i] is the "output" of reaction active_reaction_ids[i].
        // This is used to map the internal reaction number to the actual reaction id.
        std::vector<unsigned int> active_reaction_ids;
        std::vector<double> propensities;

        IntegratorData(HybridSimulation *simulation);
        IntegratorData(HybridSimulation *simulation, int num_species, int num_reactions);
        IntegratorData(IntegratorData &prev_data);
    };

    /* :IntegrationResults:
     * Organized data structure for accessing the integrator's output vector.
     * Contents are MUTABLE! Updating the values in any containing pointers
     *   will be permanently reflected in the integrator's vector.
     * 
     * All pointers in the structure point to different regions of the same vector.
     * N_Vector: [ --- concentrations --- | ---- rxn_offsets ---- ]
     */
    struct IntegrationResults
    {
        // concentrations: bounded by [0, num_species)
        realtype *concentrations;
        // reactions:      bounded by [num_species, num_species + num_reactions)
        realtype *reactions;
        int retcode;
    };

    struct URNGenerator
    {
    private:
        std::uniform_real_distribution<double> uniform;
        std::mt19937_64 rng;
        unsigned long long seed;
    public:
        double next();
        URNGenerator() = delete;
        explicit URNGenerator(unsigned long long seed);
    };

    class Integrator
    {
    private:
        void *cvode_mem;
        SUNLinearSolver solver;
        int num_species;
        int num_reactions;
        int *m_roots = nullptr;
        URNGenerator urn;
        Model<double> &model;
    public:
        // status: check for errors before using the results.
        IntegrationStatus status;
        N_Vector y;
        N_Vector y0;
        N_Vector y_save;
        realtype t;
        realtype t_save;

        /* save_state()
         * Creates a duplicate copy of the integrator's current solution vector.
         * Contents of the most recent duplicate will be restored when restore_state() is called.
         * 
         * Returns the time value of the integrator's saved state.
         */
        double save_state();

        /* restore_state()
         * Loads the most recent duplicated copy of the solution vector.
         * 
         * Returns the time value that the integrator was restored to.
         */
        double restore_state();

        /* refresh_state()
         * Loads any new changes to the solution vector without changing previous output.
         * Any new values assigned to the public N_Vector y will be recognized by the integrator.
         * The current time value remains the same. To change this, modify `t`.
         */
        void refresh_state();

        /* rereinitialize()
         * restore state to the values passed to the constructor.
         */
        void reinitialize();

        /// @brief Make events available to root-finder during integration.
        /// The root-finder itself is not activated until enable_root_finder() is called.
        ///
        /// @param events List of event objects to make available to the root-finder.
        /// The trigger functions of all given events are added as root-finder targets.
        void use_events(const std::vector<Event> &events);

        /// @brief Make events and reactions available to root-finder during integration.
        /// The root-finder itself is not activated until enable_root_finder() is called.
        ///
        /// @param events List of event objects to make available to the root-finder.
        /// @param reactions List of reaction objects to make available to the root-finder.
        void use_events(const std::vector<Event> &events, const std::vector<HybridReaction> &reactions);

        /// @brief Make reactions available to root-finder during integration.
        /// The root-finder itself is not activated until enable_root_finder() is called.
        ///
        /// @param reactions List of reaction objects to make available to the root-finder.
        void use_reactions(const std::vector<HybridReaction> &reactions);

        /// @brief Installs a CVODE root-finder onto the integrator.
        /// Any events or reactions provided by previous calls to use_events() or use_reactions()
        /// will cause the integrator to return early, which the integrate() method will indicate.
        bool enable_root_finder();

        /// @brief Removes the CVODE root-finder from the integrator.
        /// Early returns on root-finder events no longer happen,
        /// and the underlying SBML event data and reaction data are removed.
        bool disable_root_finder();

        /// @brief Configures CVODE to use the user-supplied configuration data.
        /// If all configurations were applied successfully, returns true. Otherwise, returns false.
        bool configure(SolverConfiguration config);

        void set_error_handler(CVErrHandlerFn error_handler);


		inline realtype *get_y_save_ptr()
		{
			return &N_VGetArrayPointer(y_save)[0];
		}
		inline realtype *get_y0_ptr()
		{
			return &N_VGetArrayPointer(y0)[0];
		}

		inline realtype *get_species_state()
		{
			return &N_VGetArrayPointer(y)[0];
		}
		inline realtype *get_reaction_state()
		{
			return &N_VGetArrayPointer(y)[num_species];
		}

        IntegrationResults integrate(double *t);
        IntegrationResults integrate(double *t, std::set<int> &event_roots, std::set<unsigned int> &reaction_roots);
        IntegratorData data;

        Integrator(HybridSimulation *simulation, Model<double> &model, URNGenerator urn, double reltol, double abstol);
        ~Integrator();
        void reset_model_vector();
    };

    N_Vector init_model_vector(Model<double> &model, URNGenerator urn);

    int rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data);
    int rootfn(realtype t, N_Vector y, realtype *gout, void *user_data);

    }
}
