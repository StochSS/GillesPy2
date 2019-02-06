#ifndef GILLESPY_MODEL
#define GILLESPY_MODEL
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

namespace Gillespy{

  //Represents info for a chemical reactant/product
  struct Species{
    unsigned int id; //useful for index id in arrays
    std :: string name;
    unsigned int initial_population;
  };
  
  struct Reaction{
    unsigned int id; //useful for propensity function id associated
    std :: string name;
    std :: unique_ptr<int[]> species_change; //list of changes to species with this reaction firing
    std :: vector<unsigned int> affected_reactions; //list of which reactions have propensities that would change with this reaction firing 
  };
  
  //Represents a model of reactions and species
  struct Model{
    unsigned int number_species;
    std :: unique_ptr<Species[]> species;
    unsigned int number_reactions;
    std :: unique_ptr<Reaction[]> reactions;
    Model(std :: vector<std :: string> species_names, std :: vector<unsigned int> species_populations, std :: vector<std :: string> reaction_names);
    void update_affected_reactions();
  };
  
  //Interface class to represent container for propensity functions
  class IPropensityFunction{
  public:
    virtual double evaluate(unsigned int reaction_number, unsigned int* state) = 0;
    virtual ~IPropensityFunction() {}; 
  };

  
  //Represents simulation return data
  struct Simulation{
    Model* model;
    double* timeline;
    double end_time;
    int random_seed;
    unsigned int number_timesteps;
    unsigned int number_trajectories;
    unsigned int* trajectories_1D;
    unsigned int*** trajectories;
    IPropensityFunction *propensity_function;
    Simulation(Model* model, unsigned int number_trajectories, unsigned int number_timesteps, double end_time, IPropensityFunction* propensity_function, int random_seed);
    ~Simulation();
    friend std :: ostream& operator<<(std :: ostream& os, const Simulation& simulation);
    void output_results_buffer(std :: ostream& os);
  };
}
#endif
