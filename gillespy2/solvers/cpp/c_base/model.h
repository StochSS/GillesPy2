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

		virtual double TauEvaluate(unsigned int reaction_number, const std::vector<int> &S) = 0;
		virtual double ODEEvaluate(int reaction_number, const std::vector<double> &S) = 0;

		virtual ~IPropensityFunction() {};
	};

	template <typename PType>
	struct Simulation {
		int random_seed;

		unsigned int number_timesteps;
		unsigned int number_trajectories;

		double current_time;
		double end_time;
		double *timeline;

		PType *trajectories_1D;
		PType ***trajectories;

		Model<PType> *model;

		IPropensityFunction *propensity_function;

		template <class T> friend std::ostream &operator << (std::ostream &os, const Simulation<T> &simulation);

		void output_results_buffer(std::ostream &os);
		~Simulation();
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
	// Stochastic Model: species state represented using discrete population values.
	extern template struct Model<unsigned int>;
	// Stochastic Simulation: run using a discrete-valued model.
	extern template struct Simulation<unsigned int>;
	extern template void init_simulation<unsigned int>(Model<unsigned int> *model, Simulation<unsigned int> &simulation);
}
