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

#include "simulation.h"

namespace Gillespy
{
    template <typename TNum>
    void init_timeline(Model<TNum> *model, Simulation<TNum> &simulation)
    {
        double timestep_size = simulation.end_time / (simulation.number_timesteps - 1);
        simulation.timeline = new double[simulation.number_timesteps];

        for (unsigned int i = 0; i < simulation.number_timesteps; ++i) {
            simulation.timeline[i] = timestep_size * i;
        }
    }

    template <typename TNum>
    void init_simulation(Model<TNum> *model, Simulation<TNum> &simulation)
    {
        init_timeline(model, simulation);

        simulation.current_state = new TNum[model->number_species];
        // Output interval must lie within the range (0, num_timesteps].
        // An output interval of 0 signifies to output entire trajectories.
        if (simulation.output_interval == 0 || simulation.output_interval > simulation.number_timesteps)
        {
            simulation.output_interval = simulation.number_timesteps;
        }
    }

    template <typename TNum>
    Simulation<TNum>::~Simulation()
    {
        delete[] timeline;
        delete[] current_state;
    }

    template <typename TNum>
    std::ostream &operator << (std::ostream &os, const Simulation<TNum> &simulation)
    {
        for (unsigned int timestep = 0; timestep < simulation.number_timesteps; timestep++) {
            os << simulation.timeline[timestep] << ' ';

            for (unsigned int trajectory = 0; trajectory < simulation.number_trajectories; trajectory++) {
                for (unsigned int species = 0; species < simulation.model->number_species; species++) {
                    os << simulation.trajectories[trajectory][timestep][species] << ' ';
                }
            }

            os << simulation.timeline[timestep] << ' ';
            os << '\n';
        }

        return os;
    }

    template <typename TNum>
    void Simulation<TNum>::output_results_buffer(std::ostream &os)
    {
        for (unsigned int trajectory = 0; trajectory < number_trajectories; trajectory++) {
            for (unsigned int timestep = 0; timestep < number_timesteps; timestep++) {
                os << timeline[timestep] << ',';

                for (unsigned int species = 0; species < model->number_species; species++) {
                    os << (double) current_state[species] << ',';
                }
            }
        }

        os << (int)current_time;
    }

    template <typename TNum>
    void Simulation<TNum>::output_buffer_range(std::ostream &os, unsigned int next_timestep)
    {
        next_timestep = std::min(number_timesteps - 1, next_timestep);
        // Each entry per timestep is a species population/concentration value.
        // If we have no species, then there's nothing to write!
        if (model->number_species == 0)
        {
            last_timestep = next_timestep;
            return;
        }

        unsigned int timestep;
        for (timestep = last_timestep; timestep <= next_timestep; ++timestep)
        {
            os << timeline[timestep];
            for (unsigned int species = 0; species < model->number_species; ++species)
            {
                os << ',' << current_state[species];
            }
            os << ',';

            if (timestep % output_interval == 0)
            {
                os.flush();
            }
        }

        last_timestep = timestep;
    }

    template <typename TNum>
    void Simulation<TNum>::reset_output_buffer(unsigned int trajectory_index)
    {
        last_timestep = 0;
        current_time = 0.0;
        trajectory_num = trajectory_index;

        for (unsigned int spec_i = 0; spec_i < model->number_species; ++spec_i)
        {
            current_state[spec_i] = model->species[spec_i].initial_population;
        }
    }

    template<typename PType>
    void Simulation<PType>::output_buffer_range(std::ostream &os)
    {
        output_buffer_range(os, last_timestep);
    }

    template<typename PType>
    void Simulation<PType>::output_buffer_final(std::ostream &os)
    {
        os << (int) current_time;
        os.flush();
    }

    template struct Simulation<double>;
    template struct Simulation<unsigned int>;
    template struct Simulation<int>;
    template void init_simulation<double>(Model<double> *model, Simulation<double> &simulation);
    template void init_simulation<unsigned int>(Model<unsigned int> *model, Simulation<unsigned int> &simulation);
    template void init_simulation<int>(Model<int> *model, Simulation<int> &simulation);
}
