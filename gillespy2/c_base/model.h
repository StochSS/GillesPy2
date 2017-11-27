#ifndef GILLESPY_MODEL
#define GILLESPY_MODEL
#include <string>
#include <vector>

namespace Gillespy{

  //Represents info for a chemical reactant/product
  struct Species{
    int id; //useful for index id in arrays
    std :: string name;
    int initial_population;
  };
  
  struct Reaction{
    int id; //useful for propensity function id associated
    std :: string name;
    int* species_change; //list of changes to species with this reaction firing
    std :: vector<int> affected_reactions; //list of which reactions have propensities that would change with this reaction firing 
  };
  
  //Represents a model of reactions and species
  struct Model{
    int number_species;
    Species* species;
    int number_reactions;
    Reaction* reactions;
  };
  
  //Build model with set number of species and reactions
  Model* build_model(int number_species, std :: vector<std :: string> name_species, int number_reactions, std :: vector<std :: string> name_reactions);
  
  //Free up space/clean up model
  void free_model(Model* model);


  //Interface class to represent container for propensity functions
  class IPropensityFunction{
  public:
    virtual double evaluate(int reaction_number, int* state) = 0;
    virtual ~IPropensityFunction() {}; 
  };

  
  //Represents simulation return data
  class Simulation{
  public:
    Model* model;
    double* timeline;
    double end_time;
    int number_timesteps;
    int number_trajectories;
    int* trajectories_1D;
    int*** trajectories;
    IPropensityFunction *propensity_function;
    Simulation(Model* model, int number_timesteps, int number_trajectories, double end_time, IPropensityFunction* propensity_function);
    ~Simulation();
  };  
}
#endif
