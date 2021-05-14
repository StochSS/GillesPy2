#ifndef GILLESPY_MODEL
#define GILLESPY_MODEL
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

namespace Gillespy{
  #define CONTINUOUS 0
  #define DISCRETE 1
  #define DYNAMIC 2
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

  #define SSA 1
  #define ODE 2
  #define TAU 3
  #define HYBRID 4

  struct Simulation{
    Model* model;
    ~Simulation();

    // the type of simulation - SSA, ODE, TAU, or HYBRID
    int type;
    // array representing discrete time steps for the simulation
    double* timeline;
    
    double end_time;
    double current_time;
    int random_seed;
    unsigned int number_timesteps;
    unsigned int number_trajectories;

    unsigned int* trajectories_1D;
    // this is your results!!
    // first dimension: trajectory by number
    // second dimension: the associated timesteps for that trajectory 
    // third dimension: the species' population at the given timestep for the given trajectory (integer because discrete stochastic) 
    unsigned int*** trajectories;
    // this is your results!!
    // first dimension: trajectory by number
    // second dimension: the associated timesteps for that trajectory
    // third dimension: the species' population at the given timestep for the given trajectory (double because ODE)
    double* trajectories_1DODE;
    double*** trajectoriesODE;

    // this is essentially just an array that tracks whether or not a given species was computed continuously or discretely for a given timestep and trajectory (iteration of simulation)
    // CONTINUOUS 0
    // DISCRETE 1
    int*** trajectoriesHYBRID;

    IPropensityFunction *propensity_function;
    friend std :: ostream& operator<<(std :: ostream& os, const Simulation& simulation);
    void output_results_buffer(std :: ostream& os);
  };
  //Trajectory initializers for ODE and SSA solvers
  // void simulationODEINIT(Model* model, Simulation &simulation);
  // void simulationSSAINIT(Model* model, Simulation &simulation);
  void simulationINIT(Model* model, Simulation &simulation);

  
}
#endif
