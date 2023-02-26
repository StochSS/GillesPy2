/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2023 GillesPy2 developers.
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

#include "model.h"

namespace Gillespy
{
    template <typename PType>
    class Simulation
    {
    public:
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

    template <typename TNum>
    void init_simulation(Model<TNum> *model, Simulation<TNum> &simulation);

    extern template struct Simulation<double>;
    extern template void init_simulation<double>(Model<double> *model, Simulation<double> &simulation);
    extern template struct Simulation<unsigned int>;
    extern template struct Simulation<int>;
    extern template void init_simulation<unsigned int>(Model<unsigned int> *model, Simulation<unsigned int> &simulation);
    extern template void init_simulation<int>(Model<int> *model, Simulation<int> &simulation);
}
