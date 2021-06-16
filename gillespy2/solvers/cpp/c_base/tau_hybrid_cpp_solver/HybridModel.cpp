#include "HybridModel.h"

namespace Gillespy::TauHybrid
{

	void simulation_hybrid_init(HybridSimulation &simulation)
	{
		Model *model = simulation.model;
		init_timeline(simulation);
		simulation.type = HYBRID;

		unsigned int trajectory_size = simulation.number_timesteps * (model -> number_species);
		/* 1-dimensional arrays might be unnecessary; look into mapping 3D array directly */
		simulation.trajectories_hybrid1D = new hybrid_state[simulation.number_trajectories * trajectory_size];
		simulation.trajectories_hybrid = new hybrid_state**[simulation.number_trajectories];

		for(unsigned int i = 0; i < simulation.number_trajectories; i++){
			simulation.trajectories_hybrid[i] = new hybrid_state*[simulation.number_timesteps];
			for(unsigned int j = 0; j < simulation.number_timesteps; j++){
				simulation.trajectories_hybrid[i][j] = &(simulation.trajectories_hybrid1D[i * trajectory_size + j *  (model -> number_species)]);
			}
		}
	}

	void HybridSimulation::output_hybrid_results(std::ostream &os)
	{
		for (int i = 0 ; i < number_trajectories; i++){
			for (int j = 0; j < number_timesteps; j++){
				os << timeline[j] << ',';

				for (int k = 0; k < model->number_species; k++) {
					os << trajectories_hybrid[i][j][k].continuous << ',';
				}
			}

			os<<(int)current_time;
		}
	}

	HybridReaction::HybridReaction()
		: mode(SimulationState::CONTINUOUS),
		  base_reaction(nullptr)
	{
		// Empty constructor body
	}

	HybridSpecies::HybridSpecies()
		: user_mode(SimulationState::CONTINUOUS),
		  partition_mode(SimulationState::CONTINUOUS),
		  switch_tol(0.03),
		  switch_min(0)
	{
		// Empty constructor body
	}

	HybridSimulation::~HybridSimulation()
	{
		if (type == HYBRID) {
			for(unsigned int i = 0; i < number_trajectories; i++){
				delete trajectories_hybrid[i];
			}
			delete trajectories_hybrid;
		}
	}

}
