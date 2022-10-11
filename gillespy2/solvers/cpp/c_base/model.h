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

#include "template.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) || defined(__WIN32__)
#include <windows.h>
#define GPY_PID_GET() ((int) GetCurrentProcessId())
#define GPY_INTERRUPT_HANDLER(handler_name, handler_code) \
static BOOL WINAPI handler_name(DWORD signum) {           \
    do handler_code while(0);                             \
    return TRUE;                                          \
}
#define GPY_INTERRUPT_INSTALL_HANDLER(handler) SetConsoleCtrlHandler(handler, TRUE)
#else
#include <unistd.h>
#include <csignal>
#define GPY_PID_GET() (getpid())
#define GPY_INTERRUPT_HANDLER(handler_name, handler_code) \
static void handler_name(int signum) {                    \
    do handler_code while(0);                             \
}
#define GPY_INTERRUPT_INSTALL_HANDLER(handler) signal(SIGINT, handler)
#endif

namespace Gillespy
{
    typedef unsigned int ReactionId;

    template <typename PType>
    struct Species {
        unsigned int id;
        PType initial_population;

        std::string name;

        // Needed by TauLeapingCSolver to hash into a set.
        bool operator < (const Species &other) const {
            return id < other.id;
        };
    };

    struct Reaction {
        unsigned int id;
        std::string name;

        // List of species that will change when this reaction fires.
        std::vector<unsigned int> affected_reactions;

        // List of reactions who's propensities will change when this reaction fires.
        std::unique_ptr<int[]> species_change;
        std::unique_ptr<int[]> reactants_change;
        std::unique_ptr<int[]> products_change;

        inline static double propensity(
                ReactionId reaction_id,
                double *state,
                double *parameters,
                const double *constants)
        {
            return map_propensity(reaction_id, state, parameters, constants);
        }

        inline static double propensity(
                ReactionId reaction_id,
                unsigned int *state,
                double *parameters,
                const double *constants)
        {
            return map_propensity(reaction_id, state, parameters, constants);
        }

        inline static double propensity(
                ReactionId reaction_id,
                int *state,
                double *parameters,
                const double *constants)
        {
            return map_propensity(reaction_id, state, parameters, constants);
        }

        inline static double propensity(ReactionId reaction_id, double *state)
        {
            return map_propensity(reaction_id, state, s_variables.get(), s_constants.get());
        }

        inline static double propensity(ReactionId reaction_id, int *state)
        {
            return map_propensity(reaction_id, state, s_variables.get(), s_constants.get());
        }

        inline static double propensity(ReactionId reaction_id, unsigned int *state)
        {
            return map_propensity(reaction_id, state, s_variables.get(), s_constants.get());
        }

        inline static void load_parameters()
        {
            s_variables = std::shared_ptr<double>(get_variables(&s_num_variables));
            s_constants = std::shared_ptr<const double>(get_constants(&s_num_constants));
        }

        static int s_num_variables;
        static int s_num_constants;
        static std::shared_ptr<double> s_variables;
        static std::shared_ptr<const double> s_constants;
    };

    template <typename PType>
    struct Model {
        unsigned int number_species;
        unsigned int number_reactions;

        std::unique_ptr<Species<PType>[]> species;
        std::unique_ptr<Reaction[]> reactions;

        void update_affected_reactions();

        Model(
            std::vector<std::string> species_names,
            std::vector<double> species_populations,
            std::vector<std::string> reaction_names
        );
    };

    template <typename PType>
    struct Simulation {
        int random_seed;

        unsigned int number_timesteps;
        unsigned int number_trajectories;
        // 0 is an invalid output interval and is instead used as a sentinel value.
        unsigned int output_interval = 0;

        double current_time;
        double end_time;
        double *timeline;

        PType *current_state;

        Model<PType> *model;

        template <class T> friend std::ostream &operator << (std::ostream &os, const Simulation<T> &simulation);

        // output_results_buffer: Writes the contents of the entire simulation trajectory.
        void output_results_buffer(std::ostream &os);

        /// \name output_buffer_range
        ///
        /// \brief Writes the contents of the simulation trajectory up to a certain index.
        /// The simulation maintains a "memory" of the last timestep it left off at.
        /// All timesteps between `last_timestep` (inclusive) and `next_timestep` (inclusive) are written.
        ///
        /// \param os Output stream to write to.
        /// \param next_timestep Which timestep index to stop writing from.
        void output_buffer_range(std::ostream &os, unsigned int next_timestep);

        /// \name output_buffer_range
        ///
        /// \brief Writes the contents of the next timestep of the simulation trajectory.
        /// The simulation maintains a "memory" of the last timestep it left off at.
        /// When no `next_timestep` is specified, it is assumed that only the next timestep is written.
        ///
        /// \param os Output stream to write to.
        void output_buffer_range(std::ostream &os);

        /// \name reset_output_buffer
        ///
        /// \brief Re-initializes the simulation's output buffer state to prepare for a new trajectory.
        /// When writing multiple trajectories, this should be called before each trajectory.
        ///
        /// \param trajectory_num Index pointing to the desired trajectory to output.
        void reset_output_buffer(unsigned int trajectory_num);

        /// \name output_buffer_final
        ///
        /// \brief Writes the final appending values of the buffer.
        /// Typically, this contains any necessary final data, like stop times or trajectory counts.
        ///
        /// \param os Output stream to write the final buffer contents to.
        void output_buffer_final(std::ostream &os);

        /// \name set_status
        ///
        /// \brief Sets the return code of the simulation program.
        /// \param status Opaque return code value.
        inline void set_status(uint8_t status) { m_status = status; }
        inline uint8_t get_status() { return m_status; }

        ~Simulation();

        enum SimulationStatus : uint8_t
        {
            OK = 0,
            PAUSED = 33,
        };

    private:
        unsigned int last_timestep;
        unsigned int trajectory_num;
        uint8_t m_status = 0;
    };

    /// \name SolverConfiguration
    /// \brief Container struct for ODE-specific configuration data.
    /// Used in ODE and Hybrid solvers to configure SUNDIALS solvers.
    struct SolverConfiguration
    {
        double rel_tol;
        double abs_tol;
        double max_step;
    };

    enum class LogLevel
    {
        INFO   = 0,
        WARN   = 1,
        ERR    = 2,
        CRIT   = 3,
        SILENT = 4,
    };

    class LogStream
    {
    public:
        virtual LogStream &operator<<(std::ostream &(*output)(std::ostream&)) { return *this; }
        virtual LogStream &operator<<(const char *output) { return *this; }
        virtual LogStream &operator<<(int output) { return *this; }
        virtual LogStream &operator<<(double output) { return *this; }
    };

    class StdErrLogStream : public LogStream
    {
    public:
        LogStream &operator<<(std::ostream &(*output)(std::ostream&)) override { std::cerr << output; return *this; }
        LogStream &operator<<(const char *output) override { std::cerr << output; return *this; }
        LogStream &operator<<(int output) override { std::cerr << output; return *this; }
        LogStream &operator<<(double output) override { std::cerr << output; return *this; }
    };

    class Logger
    {
        LogLevel m_log_level = LogLevel::CRIT;
        LogStream m_null = LogStream();
        StdErrLogStream m_stderr = StdErrLogStream();

    public:
        LogStream &info();
        LogStream &warn();
        LogStream &err();
        LogStream &crit();

        inline LogLevel get_log_level() { return m_log_level; }
        inline LogLevel set_log_level(LogLevel log_level)
        {
            LogLevel previous = m_log_level;
            m_log_level = log_level;
            return previous;
        }
    };

    template <typename TNum>
    void init_simulation(Model<TNum> *model, Simulation<TNum> &simulation);

    /* ================================= *
     * === DETERMINISTIC SIMULATIONS === *
     * ================================= */

    // Deterministic Species: species whose initial and runtime values are continuous.
    extern template struct Species<double>;
    // Deterministic Model: species state represented using concentration values.
    extern template struct Model<double>;
    // Deterministic Simulation: run using a continuous-valued model.
    extern template struct Simulation<double>;
    extern template void init_simulation<double>(Model<double> *model, Simulation<double> &simulation);

    /* ============================== *
     * === STOCHASTIC SIMULATIONS === *
     * ============================== */

    // Stochastic Species: species whose initial and runtime values are discrete.
    extern template struct Species<unsigned int>;
    extern template struct Species<int>;
    // Stochastic Model: species state represented using discrete population values.
    extern template struct Model<unsigned int>;
    extern template struct Model<int>;
    // Stochastic Simulation: run using a discrete-valued model.
    extern template struct Simulation<unsigned int>;
    extern template struct Simulation<int>;
    extern template void init_simulation<unsigned int>(Model<unsigned int> *model, Simulation<unsigned int> &simulation);
    extern template void init_simulation<int>(Model<int> *model, Simulation<int> &simulation);
}
