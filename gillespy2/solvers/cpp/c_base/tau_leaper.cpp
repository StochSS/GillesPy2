#include "tau_leaper.h"
#include "tau.h"
#include <random>//Included for mt19937 random number generator, poisson distribution
#include <cmath>//Included for natural logarithm
#include <string.h>//Included for memcpy only
#include <csignal>//Included for timeout signal handling
#include <iostream> // For testing output

namespace Gillespy{
bool interrupted = false;

void signalHandler(int signum){
  interrupted = true;
}

std::pair<std::map<std::string,int>,double> get_reactions(const Gillespy::Model *model, const std::vector<double> &propensity_values, double tau_step, double current_time, double save_time){
    /*
     * Helper Function to get reactions fired from t to t+tau. Effects two values:
     *rxn_count - dict with key=Reaction channel value=number of times fired
     *curr_time - float representing current time
     */

    if (current_time + tau_step > save_time)
        tau_step = save_time - current_time;

    std::map<std::string, int> rxn_count; // map of how many times reaction is fired
    std::random_device rd;
    std::mt19937 generator(rd());
    std::pair<std::map<std::string,int>,double> values; // value pair to be returned, map of times {map of times reaction fired, current time}

    for (int i = 0; i < model->number_reactions; i++){
        std::poisson_distribution<int> poisson(propensity_values[i]*tau_step);
        rxn_count[model->reactions[i].name] = poisson(generator);
    }
    current_time = current_time+tau_step;
    values.first = rxn_count;
    values.second = current_time;
    return values;
}

void tau_leaper(Gillespy::Simulation* simulation, const double tau_tol){
    signal(SIGINT, signalHandler);
    if(simulation){
        //Initialize your tau args
        TauArgs tau_args = initialize(*(simulation->model),tau_tol);

        double increment = 1; // CHANGE THIS

        //Initialize current_state variables, propensity_values
        std::vector<int> current_state((simulation -> model) -> number_species);
        std::vector<double> propensity_values (simulation->model->number_reactions);


        //copy initial state for each trajectory
        for(unsigned int species_number = 0; species_number < ((simulation -> model) -> number_species); species_number++){
            simulation -> trajectories[0][0][species_number] = (simulation -> model) -> species[species_number].initial_population;
        }


        //Simulate for each trajectory
        for(unsigned int trajectory_number = 0; trajectory_number < simulation -> number_trajectories; trajectory_number++){
            if(interrupted){
                break;
            }

            for (int spec = 0; spec< simulation->model->number_species; spec++){
                current_state[spec] = (simulation->model->species[spec].initial_population);
            }



            //Initialize simulation variables
            double current_time = 0;
            unsigned int entry_count = 0;
            //Propensity sum initialization, to be added to later.
            double propensity_sum;
            //Start save time
            double save_time = 0;
            //Variable to keep track of rejected steps, debug
            int steps_rejected = 0;
            //Initialize tau_step, will be assigned using tau::select()
            double tau_step;
            std::vector <int> prev_curr_state;

            // Each save step
            while (entry_count < simulation->end_time+1){ // while less than end_time? Could be incorrect
                if (interrupted)
                    break;
                while(current_time < save_time){
                    if(interrupted) // If timeout, or keyboard interrupt
                        break;
                    //calculate propensities for each step
                    for(unsigned int reaction_number = 0; reaction_number < ((simulation -> model) -> number_reactions); reaction_number++){
                        propensity_values[reaction_number] = (simulation -> propensity_function -> evaluate(reaction_number, current_state));
                    }


                    tau_step = select(*(simulation->model), tau_args, tau_tol, current_time, save_time, propensity_values, current_state); // tau_step selection process

                    prev_curr_state = current_state;
                    double prev_curr_time = current_time;
                    int loop_cnt = 0;

                    while (true){

                        std::map<std::string, int> rxn_count;
                        std::pair<std::map<std::string,int>,double> values;
                        values = get_reactions(simulation -> model, propensity_values, tau_step, current_time, save_time);
                        rxn_count = values.first; // How many times reactions in model fired over a tau step
                        current_time  = values.second;

                        std::map<int,bool> species_modified;
                        for (int i = 0; i<simulation->model->number_reactions; i++){
                            if (rxn_count[simulation->model->reactions[i].name] > 0)
                                for (auto const &spec : tau_args.reactions_reactants[i]){
                                        species_modified[spec] = true;
                                        //+= for both reactants and products because current_state is represented with negative number changes for reactants, and positive for products.
                                        current_state[spec] += simulation->model->reactions[i].species_change[spec] * rxn_count[simulation->model->reactions[i].name];

                                    }
                                for (auto const &spec : tau_args.products[i]){
                                    species_modified[spec] = true;
                                    current_state[spec] += simulation->model->reactions[i].species_change[spec] * rxn_count[simulation->model->reactions[i].name];
                                }
                         }


                        bool neg_state = false;
                        for (auto const& x : species_modified)
                           {
                               if (current_state[x.first] < 0)
                                    neg_state = true;
                           }

                        if (neg_state == true){
                            current_state = prev_curr_state;
                            current_time = prev_curr_time;
                            tau_step /= 2;
                            }

                        else
                            break; // out of while true
                        }
                    }
                for (int i = 0; i<simulation->model->number_species; i++)
                    simulation->trajectories[trajectory_number][entry_count][i] = current_state[i];

                save_time += increment;
                entry_count += 1;
              }
                /*
                std::cout<<"TRAJECTORY: "<< trajectory_number<<std::endl;
                for (int entry = 0; entry<simulation->end_time+1; entry++){
                    std::cout<<"[";
                    for (int spec = 0; spec<simulation->model->number_species; spec++){
                        if (spec == 0)
                            std::cout<<simulation->trajectories[trajectory_number][entry][spec]<<"("<<simulation->model->species[spec].name<<")";
                        else
                            std::cout<<", "<<simulation->trajectories[trajectory_number][entry][spec]<<"("<<simulation->model->species[spec].name<<")";
                        }
                    std::cout<<"]";
                    std::cout<<std::endl;
                    }
                */


            }
        }
     }
  }


