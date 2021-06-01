#include "statistics.h"
#include <random>

namespace Gillespy::TauHybrid::Statistics
{

    void init_species_mode(const Model &model, Simulation &simulation){
        int num_species = model.number_species;
        // int num_det_species = 0;
        for (int s = 0; s < num_species; s++){
            // if the user chooses discrete, initialise the partition flag to such.
            if (model.species[s].user_mode == DISCRETE){
                model.species[s].partition_mode = DISCRETE;
            }
            // otherwise, either the user chose continuous or dynamic (or null).
            // in any case, just initialise to continuous.
            else {
                model.species[s].partition_mode = CONTINUOUS;
                // num_det_species++;
            }
        }
        // simulation.number_det_species = num_det_species;
    }

	void partition_species(const Model &model, const std::vector<double> &propensity_values, std::vector<hybrid_state> curr_state, double tau_step, double current_time, std::map<int, bool> &det_species){
		// coefficient of variance- key:species id, value: cv
		std::map<int, double> cv;
		// means
		std::map<int, double> means;
		// standard deviation
		std::map<int, double> sd;
		//init means
		for (int i = 0; i < model.number_species; ++i){
			if (model.species[i].user_mode == DYNAMIC){
				if (model.species[i].partition_mode == CONTINUOUS){
					means.insert({i, curr_state[i].continuous});
				}else {
					means.insert({i, curr_state[i].discrete});
				}
			}
		}
		//init sd's
		for (int i = 0; i < model.number_species; ++i){
			if (model.species[i].user_mode == DYNAMIC){
				sd.insert({i, 0});
			}
		}
		// calculate means and standard deviations for dynamic-mode species involved in reactions
		for (int r = 0; r < model.number_reactions; ++r){
			for (int s = 0; s < model.number_species; ++s){
				// access list of species by accessing the correct element of the state-change vector (of the reaction)
				if (model.species[model.reactions[r].species_change[s]].user_mode == DYNAMIC ){
					// if less than 0, that means this is a reactant
					if (model.reactions[r].species_change[s] < 0){
						means[s] -= (tau_step * propensity_values[r] * model.reactions[r].species_change[s]);
						sd[s] += std::pow((tau_step * propensity_values[r] * model.reactions[r].species_change[s]),2);
					}
					// if greater than 0, that means this is a product
					if (model.reactions[r].species_change[s] > 0){
						means[s] += (tau_step * propensity_values[r] * model.reactions[r].species_change[s]);
						sd[s] += std::pow((tau_step * propensity_values[r] * model.reactions[r].species_change[s]),2);
					}
				}
			}
		}
		// calculate coefficient of variation using means and sd
		
		for (int s = 0; s < model.number_species; ++s){
			if (means.count(s) > 0){
				Species sref = model.species[s];
				if (sref.switch_min == 0) { // (default value means switch min not set, use switch tol)
					if (means[s] > 0){
						cv[s] = sd[s] / means[s];
					}else{
						cv[s] = 1;
					}
					if (cv[s] < sref.switch_tol){
						sref.partition_mode == CONTINUOUS;
					}else{
						sref.partition_mode == DISCRETE;
					}
				}else{
					if (means[s] > sref.switch_min){
						sref.partition_mode == CONTINUOUS;
					}else{
						sref.partition_mode == DISCRETE;
					}
				}
			}
		}

        return; // TODO
	}

    
	std::pair<std::map<std::string, int>, double> get_reactions(const Gillespy::Model *model, const std::vector<double> &propensity_values, double tau_step, double current_time, double save_time)
	{
	    /* Helper Function to get reactions fired from t to t+tau. Affects two values:
         * rxn_count - dict with key=Reaction channel value=number of times fired
         * curr_time - float representing current time
         */

		if (current_time + tau_step > save_time)
			tau_step = save_time - current_time;

		std::map<std::string, int> rxn_count; // map of how many times reaction is fired
		std::random_device rd;
		std::mt19937 generator(rd());
		std::pair<std::map<std::string, int>, double> values; // value pair to be returned, map of times {map of times reaction fired, current time}

		for (int i = 0; i < model->number_reactions; i++)
		{
			std::poisson_distribution<int> poisson(propensity_values[i] * tau_step);
			rxn_count[model->reactions[i].name] = poisson(generator);
		}
		current_time = current_time + tau_step;
		values.first = rxn_count;
		values.second = current_time;
		return values;
	}
}
