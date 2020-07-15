#ifndef TAU_H
#define TAU_H
#include "model.h"
#include <map>
#include <string>
#include <functional>
/*
    HOR                  Highest Order Reaction of species
    reactants            a list of all species in the model which act as reactants
    mu_i                 mu_i for each species
    sigma_i              sigma_i squared for each species
    g_i                  Relative species error allowance denominator
    epsilon_i            Relative error allowance of species
    critical_threshold   Reactant Population to be considered critical
 */
struct TauArgs{
    std::map<std::string,int> HOR;
    std::map<Gillespy::Species, int> reactants;
    std::map<std::string,int> mu_i;
    std::map<std::string,int> sigma_i;
    std::map<std::string, std::function<double(double)>> g_i;
    std::map<std::string,double> epsilon_i;
    int critical_threshold = 10;
};

int initialize(const Gillespy::Model &model, double tau_tol);
#endif // TAU_H
