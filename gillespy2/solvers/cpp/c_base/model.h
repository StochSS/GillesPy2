#ifndef GILLESPY_MODEL
#define GILLESPY_MODEL
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <cstdint>

namespace Gillespy{

  //Represents info for a chemical reactant/product
  struct Species{
    uint32_t id; //useful for index id in arrays
    std :: string name;
    uint32_t initial_population;
  };
  
  struct Reaction{
    uint32_t id; //useful for propensity function id associated
    std :: string name;
    std :: unique_ptr<int[]> species_change; //list of changes to species with this reaction firing
    std :: vector<uint32_t> affected_reactions; //list of which reactions have propensities that would change with this reaction firing 
  };
  
  //Represents a model of reactions and species
  struct Model{
    uint32_t number_species;
    std :: unique_ptr<Species[]> species;
    uint32_t number_reactions;
    std :: unique_ptr<Reaction[]> reactions;
    Model(std :: vector<std :: string> species_names, std :: vector<uint32_t> species_populations, std :: vector<std :: string> reaction_names);
    void update_affected_reactions();
  };
  
  //Interface class to represent container for propensity functions
  class IPropensityFunction{
  public:
    virtual double evaluate(uint32_t reaction_number, uint32_t* state) = 0;
    virtual ~IPropensityFunction() {}; 
  };

  
  //Represents simulation return data
  struct Simulation{
    Model* model;
    double* timeline;
    double end_time;
    int random_seed;
    uint32_t number_timesteps;
    uint32_t number_trajectories;
    uint32_t* trajectories_1D;
    uint32_t*** trajectories;
    IPropensityFunction *propensity_function;
    Simulation(Model* model, uint32_t number_trajectories, uint32_t number_timesteps, double end_time, IPropensityFunction* propensity_function, int random_seed);
    ~Simulation();
    friend std :: ostream& operator<<(std :: ostream& os, const Simulation& simulation);
    void output_results_buffer(std :: ostream& os);
  };
}
#endif
