#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <time.h>
#include <math.h>
#include "model.h"
#include "ssa.h"
#include "tau.h"
#include "tau_leaper.h"

using namespace Gillespy;

//Default values, replaced with command line args
unsigned int number_trajectories = 1; // CHANGE BACK TO 0 WHEN NOT TESTING GOOD LORD
unsigned int number_timesteps = 201; // ALSO CHANGE BACK TO 0, GOOD LORD
int random_seed = 0;
double end_time = 200.0; // PLEase< PLEAAELPEKAPKPASE CHANGE BACK TO 0
bool seed_time = true;

//Default constants
// __DEFINE_CONSTANTS__
const double V = 1.0;

std :: string s_names[] = {
"A",
"C",
"D1",
"D2",
"Da_prime",
"Dr_prime",
"Ma",
"Mr",
"R"
};

unsigned int populations[] = {
10,
10,
1,
1,
0,
0,
0,
0,
10
};

std :: string r_names[] = {
"S_A2",
"S_A3",
"S_Mr1",
"S_Mr2",
"a_A",
"a_Ma",
"a_Mr",
"a_R",
"s_A1",
"s_C",
"s_Da",
"s_Da_prime",
"s_Dr",
"s_Dr_prime",
"s_Ma1",
"s_Ma2",
"s_R1",
"s_r2"
};

const double P0 = 50.0;

const double P1 = 100.0;

const double P2 = 50.0;

const double P3 = 500.0;

const double P4 = 0.01;

const double P5 = 50.0;

const double P6 = 50.0;

const double P7 = 5.0;

const double P8 = 1.0;

const double P9 = 10.0;

const double P10 = 0.5;

const double P11 = 0.2;

const double P12 = 1.0;

const double P13 = 2.0;

const double P14 = 1.0;

class PropensityFunction : public IPropensityFunction{
public:
    
  double evaluate(unsigned int reaction_number, const std::vector<int> &S){
    switch(reaction_number){
                    case 0:
                        return P0*S[4];


                    case 1:
                        return P0*S[5];


                    case 2:
                        return P5*S[5];


                    case 3:
                        return P4*S[3];


                    case 4:
                        return P13*S[0];


                    case 5:
                        return P9*S[6];


                    case 6:
                        return P10*S[7];


                    case 7:
                        return P11*S[8];


                    case 8:
                        return P6*S[6];


                    case 9:
                        return P13*S[0]*S[8]/V;


                    case 10:
                        return P0*S[4];


                    case 11:
                        return P12*S[0]*S[2]/V;


                    case 12:
                        return P1*S[5];


                    case 13:
                        return P14*S[0]*S[3]/V;


                    case 14:
                        return P3*S[4];


                    case 15:
                        return P2*S[2];


                    case 16:
                        return P7*S[7];


                    case 17:
                        return P8*S[1];


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
  model.reactions[0].species_change[0] = 1;

  model.reactions[1].species_change[0] = 1;

  model.reactions[2].species_change[7] = 1;

  model.reactions[3].species_change[7] = 1;

  model.reactions[4].species_change[0] = -1;

  model.reactions[5].species_change[6] = -1;

  model.reactions[6].species_change[7] = -1;

  model.reactions[7].species_change[8] = -1;

  model.reactions[8].species_change[0] = 1;

  model.reactions[9].species_change[0] = -1;

  model.reactions[9].species_change[1] = 1;

  model.reactions[9].species_change[8] = -1;

  model.reactions[10].species_change[2] = 1;

  model.reactions[10].species_change[4] = -1;

  model.reactions[11].species_change[0] = -1;

  model.reactions[11].species_change[2] = -1;

  model.reactions[11].species_change[4] = 1;

  model.reactions[12].species_change[3] = 1;

  model.reactions[12].species_change[5] = -1;

  model.reactions[13].species_change[0] = -1;

  model.reactions[13].species_change[3] = -1;

  model.reactions[13].species_change[5] = 1;

  model.reactions[14].species_change[6] = 1;

  model.reactions[15].species_change[6] = 1;

  model.reactions[16].species_change[8] = 1;

  model.reactions[17].species_change[1] = -1;

  model.reactions[17].species_change[8] = 1;

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
  double tau_tol = 0.01;
  IPropensityFunction *propFun = new PropensityFunction();
  Simulation simulation(&model, number_trajectories, number_timesteps, end_time, propFun, random_seed, simulation.current_time);
  tau_leaper(&simulation, tau_tol);
  delete propFun;
  std::cout<<"SIM FINISHED!!"<<std::endl;
  return 0;
}
