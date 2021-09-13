import gillespy2
from gillespy2.solvers.cpp.c_decoder import BasicSimDecoder
from gillespy2.solvers.utilities import solverutils as cutils
from gillespy2.core import GillesPySolver, gillespyError, Model

from .c_solver import CSolver, SimulationReturnCode

class TauHybridCSolver(GillesPySolver, CSolver):
    name = "TauHybridCSolver"
    target = "hybrid"

    def __init__(self, model: Model = None, output_directory: str = None, delete_directory: bool = True, resume=None, variable=False):
        options = None if model is None else TauHybridCSolver.__create_template_options(list(model.listOfSpecies.values()))
        super().__init__(model, output_directory, delete_directory, resume, variable, options)

    @classmethod
    def __create_template_options(cls, species: "list[gillespy2.Species]"):
        """
        Populate the given list of species modes into a set of template macro definitions.
        Generated options are specific to the Tau Hybrid solver,
          and get passed as custom definitons to the build engine.

        :param species: Ordered list of GillesPy2 species to generate options for.
        :return: Dictionary containing key-value pairs representing macro definitions.
        """
        species_mode_map = {
            "continuous": "CONTINUOUS_MODE",
            "discrete": "DISCRETE_MODE",
            "dynamic": "DYNAMIC_MODE",
        }

        species_mode_list = []
        for spec_id, spec in enumerate(species):
            # Continuous by default
            mode_keyword = species_mode_map.get(spec.mode, species_mode_map["dynamic"])
            species_mode_list.append(f"SPECIES_MODE({spec_id},{mode_keyword},{spec.switch_min})")

        return {
            f"GPY_HYBRID_SPECIES_MODES": " ".join(species_mode_list)
        }

    def get_solver_settings(self):
        """
        :return: Tuple of strings, denoting all keyword argument for this solvers run() method.
        """
        return ('model', 't', 'number_of_trajectories', 'timeout', 'increment', 'seed', 'debug', 'profile')

    def run(self=None, model: Model = None, t: int = 20, number_of_trajectories: int = 1, timeout: int = 0,
            increment: int = 0.05, seed: int = None, debug: bool = False, profile: bool = False, variables={}, 
            resume=None, tau_step: int = .03, tau_tol=0.03, **kwargs):

        if self is None or self.model is None:
            self = TauHybridCSolver(model, resume=resume)

        # Validate parameters prior to running the model.
        self._validate_type(variables, dict, "'variables' argument must be a dictionary.")
        self._validate_variables_in_set(variables, self.species + self.parameters)
        self._validate_resume(t, resume)
        self._validate_kwargs(**kwargs)
        self._validate_sbml_features({
            "Rate Rules": len(model.listOfRateRules),
            "Assignment Rules": len(model.listOfAssignmentRules),
            "Events": len(model.listOfEvents),
            "Function Definitions": len(model.listOfFunctionDefinitions)
        })

        if resume is not None:
            t = abs(t - int(resume["time"][-1]))

        number_timesteps = int(round(t / increment + 1))

        args = {
            "trajectories": number_of_trajectories,
            "timesteps": number_timesteps,
            "tau_tol": tau_tol,
            "end": t
        }

        if self.variable:
            populations = cutils.update_species_init_values(model.listOfSpecies, self.species, variables, resume)
            parameter_values = cutils.change_param_values(model.listOfParameters, self.parameters, model.volume, variables)

            args.update({
                "init_pop": populations,
                "parameters": parameter_values
            })

        seed = self._validate_seed(seed)
        if seed is not None:
            args.update({
                "seed": seed
            })


        args = self._make_args(args)
        decoder = BasicSimDecoder.create_default(number_of_trajectories, number_timesteps, len(self.model.listOfSpecies))

        sim_exec = self._build(model, self.target, self.variable, False)
        sim_status = self._run(sim_exec, args, decoder, timeout)

        if sim_status == SimulationReturnCode.FAILED:
            raise gillespyError.ExecutionError("Error encountered while running simulation C++ file:\n"
                f"Return code: {int(sim_status)}.\n")

        trajectories, time_stopped = decoder.get_output()

        simulation_data = self._format_output(trajectories)
        if sim_status == SimulationReturnCode.PAUSED:
            simulation_data = self._make_resume_data(time_stopped, simulation_data, t)
        if resume is not None:
            simulation_data = self._update_resume_data(resume, simulation_data, time_stopped)
        self.simulation_data = simulation_data

        return simulation_data, int(sim_status)
