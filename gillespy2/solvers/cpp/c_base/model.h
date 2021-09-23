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

#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

namespace Gillespy {

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

	class IPropensityFunction {
	public:
		virtual double evaluate(unsigned int reaction_number, unsigned int *state) = 0;
		virtual double TauEvaluate(unsigned int reaction_number, const int *S) = 0;
		virtual double ODEEvaluate(int reaction_number, const std::vector<double> &S) = 0;

		virtual ~IPropensityFunction() {};
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

		PType *trajectories_1D;
		PType ***trajectories;
		PType *current_state;

		Model<PType> *model;

		IPropensityFunction *propensity_function;

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

		~Simulation();

	private:
		unsigned int last_timestep;
		unsigned int trajectory_num;
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
