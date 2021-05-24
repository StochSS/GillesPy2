from gillespy2.solvers.cpp.c_decoder import BasicSimDecoder
from gillespy2.core import GillesPySolver

from .c_simulation import CSimulation, SimulationReturnCode
from ..utilities import cutils

class ODECSolver(GillesPySolver, CSimulation):
    type = "ODESimulation"

    def get_solver_settings(self):
        """
        :return: Tuple of strings, denoting all keyword argument for this solvers run() method.
        """
        return ('model', 't', 'number_of_trajectories', 'timeout', 'increment', 'seed', 'debug', 'profile')

    def run(self=None, model=None, t=20, number_of_trajectories=1, timeout=0,
            increment=0.05, seed=None, debug=False, profile=False, variables={}, resume=None, **kwargs):

        t = abs(t - int(resume["time"][-1]))
        number_timesteps = int(round(t / increment + 1))
        seed = int(seed)

        args = {
            "trajectories": number_of_trajectories,
            "timesteps": number_timesteps,
            "end": t,
            "increment": increment,
            "seed": seed
        }

        if self.variable:
            populations = cutils.update_species_init_values(model.listOfSpecies, self.species, variables, resume)
            parameter_values = cutils.change_param_values(model.listOfParameters, self.parameters, model.volume, variables)

            args.update({
                "initial_values": populations,
                "parameters": parameter_values
            })

        args = self._make_args(args)
        decoder = BasicSimDecoder()

        sim_status = self._run(self.type, args, decoder, timeout)

        if sim_status == SimulationReturnCode.DONE:
            return decoder.get_output(), 0