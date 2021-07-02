/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2021 GillesPy2 developers.
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

#include "model.h"

namespace Gillespy {

	template <typename PType>
	Model<PType>::Model(
		std::vector<std::string> species_names,
		std::vector<double> species_populations,
		std::vector<std::string> reaction_names) :
		number_species(species_names.size()),
		number_reactions(reaction_names.size()) {

		species = std::make_unique<Species<PType>[]>(number_species);
		reactions = std::make_unique<Reaction[]>(number_reactions);

		for (unsigned int i = 0; i < number_species; i++) {
			species[i].id = i;
			species[i].initial_population = static_cast<PType>(species_populations[i]);
			species[i].name = species_names[i];
		}

		for (unsigned int reaction = 0; reaction < number_reactions; reaction++) {
			reactions[reaction].name = reaction_names[reaction];
			reactions[reaction].species_change = std::make_unique<int[]>(number_species);

			for (unsigned int species = 0; species < number_species; species++) {
				reactions[reaction].species_change[species] = 0;
			}

			reactions[reaction].affected_reactions = std::vector<unsigned int>();
		}
	}

	template <typename PType>
	void Model<PType>::update_affected_reactions() {
		// Clear affected_reactions for each reaction.
		for (unsigned int i = 0; i < number_reactions; i++) {
			reactions[i].affected_reactions.clear();
		}

		// Check all reactions for commong species changes -> affected reactions.
		for (unsigned int r1 = 0; r1 < number_reactions; r1++) {
			for (unsigned int r2 = 0; r2 < number_reactions; r2++) {
				for (unsigned int s = 0; s < number_species; s++) {
					if (reactions[r2].species_change[s] != 0) {
						reactions[r1].affected_reactions.push_back(r2);
					}
				}
			}
		}
	}

	template <typename TNum>
	void init_timeline(Model<TNum> *model, Simulation<TNum> &simulation) {
		double timestep_size = simulation.end_time / (simulation.number_timesteps - 1);
		simulation.timeline = new double[simulation.number_timesteps];

		for (unsigned int i = 0; i < simulation.number_timesteps; ++i) {
			simulation.timeline[i] = timestep_size * i;
		}
	}

	template <typename TNum>
	void init_simulation(Model<TNum> *model, Simulation<TNum> &simulation) {
		init_timeline(model, simulation);

		unsigned int trajectory_size = simulation.number_timesteps * (model->number_species);
		simulation.trajectories_1D = new TNum[simulation.number_trajectories * trajectory_size];
		simulation.trajectories = new TNum * *[simulation.number_trajectories];

		for (unsigned int trajectory = 0; trajectory < simulation.number_trajectories; trajectory++) {
			simulation.trajectories[trajectory] = new TNum * [simulation.number_timesteps];

			for (unsigned int timestep = 0; timestep < simulation.number_timesteps; timestep++) {
				simulation.trajectories[trajectory][timestep] =
					&(simulation.trajectories_1D[trajectory * trajectory_size + timestep * (model->number_species)]);
			}
		}
	}

	template <typename TNum>
	Simulation<TNum>::~Simulation() {
		delete[] timeline;
		delete[] trajectories_1D;

		for (unsigned int trajectory = 0; trajectory < number_trajectories; trajectory++) {
			delete[] trajectories[trajectory];
		}

		delete[] trajectories;
	}

	template <typename TNum>
	std::ostream &operator << (std::ostream &os, const Simulation<TNum> &simulation) {
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
	void Simulation<TNum>::output_results_buffer(std::ostream &os) {
		for (int trajectory = 0; trajectory < number_trajectories; trajectory++) {
			for (int timestep = 0; timestep < number_timesteps; timestep++) {
				os << timeline[timestep] << ',';

				for (int species = 0; species < model->number_species; species++) {
					os << trajectories[trajectory][timestep][species] << ',';
				}
			}
		}

		os << (int)current_time;
	}

	// DETERMINISTIC SIMULATIONS: explicit instantation of real-valued data structures.
	template struct Species<double>;
	template struct Model<double>;
	template struct Simulation<double>;

	// STOCHASTIC SIMULATIONS: explicit instantiation of discrete-valued data structures.
	template struct Species<unsigned int>;
	template struct Model<unsigned int>;
	template struct Simulation<unsigned int>;

	template void init_simulation<double>(Model<double> *model, Simulation<double> &simulation);
	template void init_simulation<unsigned int>(Model<unsigned int> *model, Simulation<unsigned int> &simulation);
}
