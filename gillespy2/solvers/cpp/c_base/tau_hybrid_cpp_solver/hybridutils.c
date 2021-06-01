#include <iostream>
#include "TauHybridCSolver.h"
#include "HybridModel.h"
#include "tau.h"
#include "model.h"

// toggle_reactions()
namespace Gillespy {
    void create_diff_eqs()
    {
    }
    // Helper method to flag reactions that can be processed deterministically (continuous change)
    // without exceeding the user-supplied tolerance

    std::set<Reaction> flag_det_rxns(const Model &model, 
                                                         const std::map<int, bool> det_species,
                                                         std::map<int, bool> det_rxns, 
                                                         std::map<int, std::set<int>> dependent_species){
        int number_rxns = model.number_reactions;
        for (int rxn = 0; rxn < number_rxns; ++rxn){
            // start with the assumption that reaction is determinstic
            det_rxns[rxn] = true;
            // iterate through the dependent species of this reaction
            for (int s = 0; s < dependent_species[rxn].size(); ++s){
                // if any of the dependencies are set by the user as discrete OR
                // have been set as dynamic and has not been flagged as deterministic,
                // allow it to be modelled discretely
                if (model.species[s].user_mode == DISCRETE){
                    det_rxns[rxn] = false;
                    break;
                } 
                if (model.species[s].user_mode == DYNAMIC && det_species[s] == false) {
                    det_rxns[rxn] = false;
                    break;
                }
            }
        }
        //create a set of all deterministic reactions
        std::set<int> new_deterministic_reactions;
        for (int rxn = 0; rxn < det_rxns.size(); ++rxn) {
            if (det_rxns[rxn] == true){
                new_deterministic_reactions.insert(rxn);
            }
        }
        return new_deterministic_reactions;
    }
    void partition_species(const Model &model, const std::vector<double> &propensity_values, std::vector<hybrid_state> curr_state, double tau_step, double current_time, std::map<int, bool> &det_species)
    {
        // coefficient of variance- key:species id, value: cv
        std::map<int, double> cv;
        // means
        std::map<int, double> means;
        // standard deviation
        std::map<int, double> sd;
        //init means
        for (int i = 0; i < model.number_species; ++i)
        {
            if (model.species[i].user_mode == DYNAMIC)
            {
                if (model.species[i].partition_mode == CONTINUOUS)
                {
                    means.insert({i, curr_state[i].continuous});
                }
                else
                {
                    means.insert({i, curr_state[i].discrete});
                }
            }
        }
        //init sd's
        for (int i = 0; i < model.number_species; ++i)
        {
            if (model.species[i].user_mode == DYNAMIC)
            {
                sd.insert({i, 0});
            }
        }
        // calculate means and standard deviations for dynamic-mode species involved in reactions
        for (int r = 0; r < model.number_reactions; ++r)
        {
            for (int s = 0; s < model.number_species; ++s)
            {
                // access list of species by accessing the correct element of the state-change vector (of the reaction)
                if (model.species[model.reactions[r].species_change[s]].user_mode == DYNAMIC)
                {
                    // if less than 0, that means this is a reactant
                    if (model.reactions[r].species_change[s] < 0)
                    {
                        means[s] -= (tau_step * propensity_values[r] * model.reactions[r].species_change[s]);
                        sd[s] += std::pow((tau_step * propensity_values[r] * model.reactions[r].species_change[s]), 2);
                    }
                    // if greater than 0, that means this is a product
                    if (model.reactions[r].species_change[s] > 0)
                    {
                        means[s] += (tau_step * propensity_values[r] * model.reactions[r].species_change[s]);
                        sd[s] += std::pow((tau_step * propensity_values[r] * model.reactions[r].species_change[s]), 2);
                    }
                }
            }
        }
        // calculate coefficient of variation using means and sd
        for (int s = 0; s < model.number_species; ++s)
        {
            if (means.count(s) > 0)
            {
                Species sref = model.species[s];
                if (sref.switch_min == 0)
                { // (default value means switch min not set, use switch tol)
                    if (means[s] > 0)
                    {
                        cv[s] = sd[s] / means[s];
                    }
                    else
                    {
                        cv[s] = 1;
                    }
                    if (cv[s] < sref.switch_tol)
                    {
                        sref.partition_mode == CONTINUOUS;
                    }
                    else
                    {
                        sref.partition_mode == DISCRETE;
                    }
                }
                else
                {
                    if (means[s] > sref.switch_min)
                    {
                        sref.partition_mode == CONTINUOUS;
                    }
                    else
                    {
                        sref.partition_mode == DISCRETE;
                    }
                }
            }
        }
        return //TODO;
    }
}
