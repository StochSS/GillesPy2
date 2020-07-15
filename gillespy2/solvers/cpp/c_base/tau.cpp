#include "tau.h"
#include <vector>
#include <iostream>
#include <memory>
#include <string>


int initialize(const Gillespy::Model &model, double tau_tol){

   TauArgs tau_args;

   std::cout << "Right before all the good stuff!" << std::endl;

   for (int i = 0; i<model.number_species; i++){
       //SEGFAULT HERE! tau_args->HOR[model.species[i].name] = 0;
       //std::cout << tau_args->HOR[model.species[i].name] << std::endl;

       std::cout << i << model.species[i].name << std::endl;
   }

   return 1;
}
