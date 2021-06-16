#include <set>
#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <random>
#include <csignal>
#include <iostream>
#include <algorithm>
#include <functional>

#include "TauLeapingSolver.h"

namespace Gillespy {
	volatile bool interrupted = false;

	void signalHandler(int signum) {
		interrupted = true;
	}

	struct TauArgs {
		int critical_threshold = 10;

		std::map<std::string, int> HOR;
		std::set<Gillespy::Species> reactants;

		std::map<std::string, std::function<double(double)>> g_i_lambdas;
		std::map<std::string, int> g_i;
		std::map<std::string, double> epsilon_i;
		std::map<int, std::vector<int>> reactions_reactants;
		std::map<int, std::vector<int>> products;
	};

	TauArgs initialize(Gillespy::Model &model, double tau_tol) {
		TauArgs tau_args;

		// Initialize highest order reactions to 0.
		for (int species_index = 0; species_index < model.number_species; species_index++) {
			tau_args.HOR[model.species[species_index].name] = 0;
		}

		for (int reaction_index = 0; reaction_index < model.number_reactions; reaction_index++) {
			int reaction_order = 0;

			for (int species_index = 0; species_index < model.number_species; species_index++) {
				if (model.reactions[reaction_index].species_change[species_index] > 0) {
					tau_args.products[reaction_index].push_back(species_index);
				}

				else if (model.reactions[reaction_index].species_change[species_index] < 0) {
					reaction_order++;

					tau_args.reactions_reactants[reaction_index].push_back(species_index);
					tau_args.reactants.insert(model.species[species_index]);
				}
			}

			if (tau_args.reactions_reactants[reaction_index].size() == 0) {
				continue;
			}

			// If this reaction's order is higher than the previous, initialize and set.
			for (auto const &reactant : tau_args.reactions_reactants[reaction_index]) {
				if (reaction_order <= tau_args.HOR[model.species[reactant].name]) {
					continue;
				}

				tau_args.HOR[model.species[reactant].name] = reaction_order;
				tau_args.g_i[model.species[reactant].name] = tau_args.HOR[model.species[reactant].name];

				int count = std::abs(model.reactions[reaction_index].species_change[reactant]);
				if (count == 2 && reaction_order == 2) {
					auto lambda = [](double x) {
						return (2 + (1 / (x - 1)));
					};
					tau_args.g_i_lambdas[model.species[reactant].name] = lambda;
				}

				else if (count == 2 && reaction_order == 3) {
					auto lambda = [](double x) {
						return ((3 / 2) * (2 + (1 / (x - 1))));
					};
					tau_args.g_i_lambdas[model.species[reactant].name] = lambda;
				}

				else if (count == 3) {
					auto lambda = [](double x) {
						return (3 + (1 / (x - 1)) + (2 / (x - 2)));
					};
					tau_args.g_i_lambdas[model.species[reactant].name] = lambda;
				}

				else {
					tau_args.g_i[model.species[reactant].name] = tau_args.HOR[model.species[reactant].name];
					tau_args.epsilon_i[model.species[reactant].name] = tau_tol / tau_args.g_i[model.species[reactant].name];
				}
			}
		}

		return tau_args;
	}

	double select(
		Gillespy::Model &model, 
		TauArgs &tau_args, 
		const double &tau_tol, 
		const double &current_time, 
		const double &save_time, 
		const std::vector<double> &propensity_values,
		const std::vector<int> &current_state) {

		bool critical = false;

		int v;

		// Tau time to step.
		double tau;

		// Smallest tau time for critical reactions.
		double critical_tau = 0;

		// Smallest tau time for non-critical reactions.
		double non_critical_tau = 0;

		std::map<std::string, double> critical_taus;
		std::map<std::string, double> mu_i;
		std::map<std::string, double> sigma_i;

		for (int species_index = 0; species_index < model.number_species; species_index++) {
			mu_i[model.species[species_index].name] = 0;
			sigma_i[model.species[species_index].name] = 0;
		}

		// Determine if there are any critical reactions and, if true, update mu_i and sigma_i.
		for (int reaction_index; reaction_index < model.number_reactions; reaction_index++) {
			for (auto const &reactant : tau_args.reactions_reactants[reaction_index]) {
				if (model.reactions[reaction_index].species_change[reactant] >= 0) {
					continue;
				}

				v = std::abs(model.reactions[reaction_index].species_change[reactant]);
				if ((double)current_state[reactant] / v < tau_args.critical_threshold && propensity_values[reaction_index] > 0) {
					critical = true;
				}

				int consumed = std::abs(model.reactions[reaction_index].species_change[reactant]);

				// Cao, Gillespie, Petzold 32a.
				mu_i[model.species[reactant].name] += consumed * propensity_values[reaction_index];
				sigma_i[model.species[reactant].name] += std::pow(consumed, 2) * propensity_values[reaction_index];
			}
		}

		// If a critical reaction is present, estimate tau for a single firing of each with a propensity > 0, taking the smallest tau.
		if (critical == true) {
			for (int reaction_index = 0; reaction_index < model.number_reactions; reaction_index++) {
				if (propensity_values[reaction_index] > 0) {
					critical_taus[model.reactions[reaction_index].name] = 1 / propensity_values[reaction_index];
				}
			}

			// Find the minimum of the critical taus.
			std::pair<std::string, double> min;
			min = *min_element(critical_taus.begin(), critical_taus.end(), [](const auto &lhs, const auto &rhs) {
				return lhs.second < rhs.second;
			});
			critical_tau = min.second;
		}

		if (tau_args.g_i_lambdas.size() > 0) {
			for (auto const &x : tau_args.g_i_lambdas) {
				tau_args.g_i[x.first] = tau_args.g_i_lambdas[x.first](tau_args.g_i[x.first]);
				tau_args.epsilon_i[x.first] = tau_tol / tau_args.g_i[x.first];

				tau_args.g_i_lambdas.erase(x.first);
			}
		}

		// Mapping of non-critical taus to be evaluated.
		std::map<std::string, double> tau_i;

		for (const auto &reactant : tau_args.reactants) {
			double calculated_max = tau_args.epsilon_i[reactant.name] * current_state[reactant.id];
			double max_pop_change_mean = std::max(calculated_max, 1.0);
			double max_pop_change_sd = std::pow(max_pop_change_mean, 2);

			if (mu_i[reactant.name] > 0) {
				tau_i[reactant.name] = std::min(std::abs(max_pop_change_mean / mu_i[reactant.name]), max_pop_change_sd / sigma_i[reactant.name]);
			}
		}

		if (tau_i.size() > 0) {
			// Find the minimum of tau_i.
			std::pair<std::string, double> min;
			min = *min_element(tau_i.begin(), tau_i.end(), [](const auto &lhs, const auto &rhs) {
				return lhs.second < rhs.second;
			});

			non_critical_tau = min.second;
		}

		// If all reactions are non-critical, use non-critical tau.
		if (critical == false) {
			tau = non_critical_tau;
		}

		// If all reactions are critical, use critical tau.
		else if (tau_i.size() == 0) {
			tau = critical_tau;
		}

		// If both critical and non-critical reactions exist, 
		// take the shortest tau between the two.
		else {
			tau = std::min(non_critical_tau, critical_tau);
		}

		if (tau > 0) {
			tau = std::max(tau, 1e-10);

			if (save_time - current_time > 0) {
				tau = std::min(tau, save_time - current_time);
			}
		}

		else {
			tau = save_time - current_time;
		}

		return tau;
	}

	std::pair<std::map<std::string, int>, double> get_reactions(
		const Gillespy::Model *model,
		const std::vector<double> &propensity_values,
		double tau_step,
		double current_time,
		double save_time) {

		/* 
		 * Helper function to get reactions fired from t to t + tau. Effects two values:
		 * - reaction_count: Map of items with key = reaction channel value and value = number of times fired.
		 * - current_time: A float value representing the current time.
		*/
		
		if (current_time + tau_step > save_time) {
			tau_step = save_time - current_time;
		}

		std::map<std::string, int> reaction_count;
		std::random_device rd;
		std::mt19937_64 generator(rd());

		// Value pairs to be returned.
		// The first is a map of times where the reaction fired, the second is the current time.
		std::pair<std::map<std::string, int>, double> values;

		for (int reaction_index = 0; reaction_index < model->number_reactions; reaction_index++) {
			std::poisson_distribution<int> poisson(propensity_values[reaction_index] * tau_step);
			reaction_count[model->reactions[reaction_index].name] = poisson(generator);
		}

		current_time = current_time + tau_step;
		values.first = reaction_count;
		values.second = current_time;

		return values;
	}

	void tau_leaper(Gillespy::Simulation<unsigned int> *simulation, const double tau_tol) {
		signal(SIGINT, signalHandler);

		if (!simulation) {
			return;
		}

		// Initialize the Tau args.
		TauArgs tau_args = initialize(*(simulation->model), tau_tol);
		double increment = simulation->timeline[1] - simulation->timeline[0];

		std::vector<int> current_state(simulation->model->number_species);
		std::vector<double> propensity_values(simulation->model->number_reactions);

		// Copy each trajectory's initial state.
		for (unsigned int species_index = 0; species_index < simulation->model->number_species; species_index++) {
			simulation->trajectories[0][0][species_index] = simulation->model->species[species_index].initial_population;
		}

		// Simulate each trajectory.
		for (unsigned int trajectory_index = 0; trajectory_index < simulation->number_trajectories; trajectory_index++) {
			if (interrupted) {
				break;
			}

			// Initialize the current state with species initial populations.
			for (int species_index = 0; species_index < simulation->model->number_species; species_index++) {
				current_state[species_index] = simulation->model->species[species_index].initial_population;
			}

			simulation->current_time = 0;

			int steps_rejected = 0;
			unsigned int entry_count = 0;
			
			// Will be assigned to later with tau::select().
			double tau_step;

			// Propensity sum initialization, will be added to later.
			double propensity_sum;
			double save_time = 0;

			std::vector<int> prev_curr_state;

			// Compute each save step.
			// Possible issue -- while less than end time? Could be incorrect.
			while (entry_count < simulation->number_timesteps) {
				if (interrupted) {
					break;
				}

				while (simulation->current_time < save_time) {
					if (interrupted) {
						break;
					}

					// Calculate propensities for each step.
					for (unsigned int reaction_index; reaction_index < simulation->model->number_reactions; reaction_index++) {
						propensity_values[reaction_index] = simulation->propensity_function->TauEvaluate(reaction_index, current_state);
					}

					tau_step = select(*(simulation->model), tau_args, tau_tol, simulation->current_time, save_time, propensity_values, current_state);

					int loop_count = 0;
					prev_curr_state = current_state;
					double prev_curr_time = simulation->current_time;

					while (true) {
						loop_count++;

						if (loop_count > 100) {
							throw std::runtime_error("Loop count exceeded 100, error");
						}

						std::map<std::string, int> reaction_count;
						std::pair<std::map<std::string, int>, double> values;

						values = get_reactions(simulation->model, propensity_values, tau_step, simulation->current_time, save_time);

						reaction_count = values.first;
						simulation->current_time = values.second;

						std::map<int, bool> species_modified;
						for (int reaction_index = 0; reaction_index < simulation->model->number_reactions; reaction_index++) {
							if (reaction_count[simulation->model->reactions[reaction_index].name] <= 0) {
								continue;
							}

							// Determine which species have been modified, and if valid, the resulting change in state.
							for (auto const &species : tau_args.reactions_reactants[reaction_index]) {
								species_modified[species] = true;

								// Add for both reactants and products because current_state is represented with negative number changes for reactants, positive for products.
								current_state[species] += 
									simulation->model->reactions[reaction_index].species_change[species] * reaction_count[simulation->model->reactions[reaction_index].name];
							}

							for (auto const &species : tau_args.products[reaction_index]) {
								species_modified[species] = true;

								// Ditto.
								current_state[species] += 
									simulation->model->reactions[reaction_index].species_change[species] * reaction_count[simulation->model->reactions[reaction_index].name];
							}
						}

						bool neg_state = false;
						for (auto const &x : species_modified) {
							if (current_state[x.first] < 0) {
								neg_state = true;
							}
						}

						if (neg_state == true) {
							current_state = prev_curr_state;
							simulation->current_time = prev_curr_time;

							tau_step = tau_step / 2;
						}

						else {
							break;
						}
					}
				}

				for (int species_index = 0; species_index < simulation->model->number_species; species_index++) {
					simulation->trajectories[trajectory_index][entry_count][species_index] = current_state[species_index];
				}

				save_time += increment;
				entry_count++;
			}
		}
	}
}
