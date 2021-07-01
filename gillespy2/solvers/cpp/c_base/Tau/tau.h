#ifndef TAU_H
#define TAU_H

#include "model.h"

#include <map>
#include <vector>
#include <set>
#include <functional>


namespace Gillespy {

    template<typename PType>
    struct TauArgs
    {
        //Highest Order Reaction
        std::map<std::string, int> HOR;
        std::set<Species<PType>> reactants;
        //Below are g_i_lambdas, pop element when used
        std::map<std::string, std::function<double(double)>> g_i_lambdas;
        std::map<std::string, int> g_i;
        std::map<std::string, double> epsilon_i;
        std::map<int, std::vector<int>> reactions_reactants;
        std::map<int, std::vector<int>> products;
        int critical_threshold = 10;
    };

    template<typename PType>
    TauArgs<PType> initialize(Model<PType> &model, double tau_tol);

    template<typename PType>
    double select(
        Model<PType> &model,
        TauArgs<PType> &tau_args,
        const double &tau_tol,
        const double &current_time,
        const double &save_time,
        const std::vector<double> &propensity_values,
        const std::vector<int> &current_state);

    // Continuous tau args is used by hybrid solver
    template struct TauArgs<double>;
    extern template TauArgs<double> initialize(Model<double> &model, double tau_tol);
    extern template double select<double>(
        Model<double> &model,
        TauArgs<double> &tau_args,
        const double &tau_tol,
        const double &current_time,
        const double &save_time,
        const std::vector<double> &propensity_values,
        const std::vector<int> &current_state);

    // Discrete tau args is used by tau leaping solver
    template struct TauArgs<unsigned int>;
    extern template TauArgs<unsigned int> initialize(Model<unsigned int> &model, double tau_tol);
    extern template double select<unsigned int>(
        Model<unsigned int> &model,
        TauArgs<unsigned int> &tau_args,
        const double &tau_tol,
        const double &current_time,
        const double &save_time,
        const std::vector<double> &propensity_values,
        const std::vector<int> &current_state);
}
#endif //TAU_H