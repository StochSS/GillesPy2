#ifndef GILLESPY_MODEL
#define GILLESPY_MODEL
#include <string>
#include <vector>

namespace Gillespy{

  //Represents info for a chemical reactant/product
  struct Species{
    uint id; //useful for index id in arrays
    std :: string name;
    uint initial_population;
  };
  
  struct Reaction{
    uint id; //useful for propensity function id associated
    std :: string name;
    int* species_change; //list of changes to species with this reaction firing
    std :: vector<uint> affected_reactions; //list of which reactions have propensities that would change with this reaction firing 
  };
  
  //Represents a model of reactions and species
  struct Model{
    uint number_species;
    Species* species;
    uint number_reactions;
    Reaction* reactions;
  };
  
  //Build model with set number of species and reactions
  Model* build_model(std :: vector<std :: string> name_species, std :: vector<std :: string> name_reactions);
  
  //Free up space/clean up model
  void free_model(Model* model);


  //Interface class to represent container for propensity functions
  class IPropensityFunction{
  public:
    virtual double evaluate(uint reaction_number, uint* state) = 0;
    virtual ~IPropensityFunction() {}; 
  };

  
  //Represents simulation return data
  class Simulation{
  public:
    Model* model;
    double* timeline;
    double end_time;
    int random_seed;
    uint number_timesteps;
    uint number_trajectories;
    uint* trajectories_1D;
    uint*** trajectories;
    IPropensityFunction *propensity_function;
    Simulation(Model* model, uint number_trajectories, uint number_timesteps, double end_time, IPropensityFunction* propensity_function, int random_seed);
    ~Simulation();
  };  
}
#endif
