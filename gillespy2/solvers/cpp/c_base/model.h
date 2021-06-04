#ifndef GILLESPY_MODEL
#define GILLESPY_MODEL
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

namespace Gillespy{

	//Represents info for a chemical reactant/product
	struct Species{
		unsigned int id; //useful for index id in arrays
		std :: string name;
		unsigned int initial_population;
		//Used for hashing into set, for TauLeapingCSolver
		bool operator < (const Species &other) const { return id < other.id; }
	};

	struct Reaction{
		unsigned int id; //useful for propensity function id associated
		std :: string name;
		std :: unique_ptr<int[]> species_change; //list of changes to species with this reaction firing
		std :: vector<unsigned int> affected_reactions; //list of which reactions have propensities that would change with this reaction firing
	};
  
	//Represents a model of reactions and species
	struct Model{
		void update_affected_reactions();
		unsigned int number_species;
		std :: unique_ptr<Species[]> species;
		unsigned int number_reactions;
		std :: unique_ptr<Reaction[]> reactions;
		Model(std :: vector<std :: string> species_names, std :: vector<unsigned int> species_populations, std :: vector<std :: string> reaction_names);
	};
  
	//Interface class to represent container for propensity functions
	class IPropensityFunction{
	public:
		virtual double evaluate(unsigned int reaction_number, unsigned int* state) = 0;
		virtual double TauEvaluate(unsigned int reaction_number, const std::vector<int> &S) = 0;
		virtual double ODEEvaluate(int reaction_number, const std::vector <double> &S) = 0;
		virtual ~IPropensityFunction() {};
	};

	template <typename PType>
	struct Simulation{
		Model* model;
		~Simulation();

		double* timeline;
		double end_time;
		double current_time;
		int random_seed;

		unsigned int number_timesteps;
		unsigned int number_trajectories;

		PType *trajectories_1D;
		PType ***trajectories;

		IPropensityFunction *propensity_function;
		template<class T> friend std :: ostream& operator<<(std :: ostream& os, const Simulation<T> &simulation);
		void output_results_buffer(std :: ostream& os);
	};

	/**
	 * Trajectory initializer function
	 * Populates the simulation object's data based on the given model
	 */
	template <typename TNum>
	void init_simulation(Model *model, Simulation<TNum> &simulation);

	// Explicitly instantiate each needed type of simulation
	extern template struct Simulation<double>;
	extern template struct Simulation<unsigned int>;
	extern template void init_simulation<double>(Model *model, Simulation<double> &simulation);
	extern template void init_simulation<unsigned int>(Model *model, Simulation<unsigned int> &simulation);

}
#endif
