#ifndef TAU_H
#define TAU_H
#include "model.h"
#include <map>
#include <string>
#include <functional>
#include <set>
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
    std::set<Gillespy::Species> reactants;
    //Below are g_i_lambdas, pop element when used
    std::map<std::string, std::function<double(double)>> g_i_lambdas;
    std::map<std::string,int> g_i;
    std::map<std::string,double> epsilon_i;
    int critical_threshold = 10;
};

TauArgs initialize(Gillespy::Model &model, const double tau_tol);
// SHOULD I ACTUALLY SEND a &TauARGS? DO TAU ARGS NEED TO BE MODIFIED? OR DO THEY STAY THE SAME AFTER INITIALIZING?
double select(Gillespy::Model &model, const TauArgs tau_args, const double tau_tol, const double current_time, const int save_time, double *propensity_values, int* current_state);
#endif // TAU_H
