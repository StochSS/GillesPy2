#ifndef TAU_H
#define TAU_H

#include "model.h"

#include <map>
#include <vector>
#include <set>
#include <functional>


namespace Gillespy {
    struct TauArgs
    {
        //Highest Order Reaction
        std::map<std::string, int> HOR;
        std::set<Gillespy::Species> reactants;
        //Below are g_i_lambdas, pop element when used
        std::map<std::string, std::function<double(double)>> g_i_lambdas;
        std::map<std::string, int> g_i;
        std::map<std::string, double> epsilon_i;
        std::map<int, std::vector<int>> reactions_reactants;
        std::map<int, std::vector<int>> products;
        int critical_threshold = 10;
    };
    TauArgs initialize(Gillespy::Model &model, double tau_tol);

    double select(Gillespy::Model &model, TauArgs &tau_args, const double &tau_tol, const double &current_time, const double &save_time, const std::vector<double> &propensity_values, const std::vector<int> &current_state);

    
}
#endif //TAU_H