#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <time.h>
#include <math.h>
#include "model.h"
#include "ODECSolver.h"
using namespace Gillespy;

//Default values, replaced with command line args
unsigned int number_trajectories = 1; // CHANGE BACK TO 0 WHEN NOT TESTING
unsigned int number_timesteps = 101; // ALSO CHANGE BACK TO 0,
int random_seed = 0;
double end_time = 100.0;
bool seed_time = true;
double increment = 1;

//Default constants
// __DEFINE_CONSTANTS__
const double V = 1.0;
std::string s_names[] = {
"Enzyme",
"Enzyme_Substrate_Complex",
"Product",
"Substrate",
};
unsigned int populations[] = {
120,
0,
0,
301
};
std::string r_names[] = {
"r1",
"r2",
"r3"
};
const double P0 = 0.0017;
const double P1 = 0.5;
const double P2 = 0.1;

class PropensityFunction : public IPropensityFunction{
public:

  double evaluate(int reaction_number, std::vector <double> S){
  switch(reaction_number){
            case 0:
                return P0*S[0]*S[3]/V;
            case 1:
                return P1*S[1];
            case 2:
                return P2*S[1];


    default: //Error
      return -1;
    }
  }

};

int main(int argc, char* argv[]){
  std :: vector<std :: string> species_names(s_names, s_names + sizeof(s_names)/sizeof(std :: string));
  std :: vector<unsigned int> species_populations(populations, populations + sizeof(populations)/sizeof(populations[0]));
  std :: vector<std :: string> reaction_names(r_names, r_names + sizeof(r_names)/sizeof(std :: string));

  Model model(species_names, species_populations, reaction_names);

  //Begin reaction species changes
//HERE WRITES MAN: __DEFINE_REACTIONS_
  model.reactions[0].species_change[0] = -1;

  model.reactions[0].species_change[1] = 1;

  model.reactions[0].species_change[3] = -1;

  model.reactions[1].species_change[0] = 1;

  model.reactions[1].species_change[1] = -1;

  model.reactions[1].species_change[3] = 1;

  model.reactions[2].species_change[0] = 1;

  model.reactions[2].species_change[1] = -1;

  model.reactions[2].species_change[2] = 1;

  //End reaction species changes
  model.update_affected_reactions();

  //Parse command line arguments
 std :: string arg;
 for(int i = 1; i < argc - 1; i++){
   arg = argv[i];
   if(argc > i+1 && arg.size() > 1 && arg[0] == '-'){
     std :: stringstream arg_stream(argv[i+1]);
     switch(arg[1]){
     case 's':
       arg_stream >> random_seed;
       seed_time = false;
       break;
     case 'e':
       arg_stream >> end_time;
       break;
      case 'i':
        arg_stream >> increment;
        break;
     case 't':
       if(arg[2] == 'r'){
         arg_stream >> number_trajectories;
       }else if(arg[2] == 'i'){
         arg_stream >> number_timesteps;
       }
       break;
     }
   }
 }

 if(seed_time){
   random_seed = time(NULL);
 }
  IPropensityFunction *propFun = new PropensityFunction();
  Simulation simulation(&model, number_trajectories, number_timesteps, end_time, propFun, random_seed, simulation.current_time, 1);

  ODESolver(&simulation,increment);
  //simulation.output_results_buffer(std :: cout);
  delete propFun;
  return 0;
}