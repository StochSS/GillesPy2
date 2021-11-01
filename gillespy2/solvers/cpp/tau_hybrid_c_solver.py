import gillespy2
from gillespy2.solvers.cpp.c_decoder import IterativeSimDecoder
from gillespy2.solvers.utilities import solverutils as cutils
from gillespy2.core import GillesPySolver, gillespyError, Model
from typing import Union
from gillespy2.core import Results

from .c_solver import CSolver, SimulationReturnCode
from gillespy2.solvers.cpp.build.template_gen import SanitizedModel

class TauHybridCSolver(GillesPySolver, CSolver):
    name = "TauHybridCSolver"
    target = "hybrid"

    @classmethod
    def __create_options(cls, model: "SanitizedModel") -> "SanitizedModel":
        """
        Populate the given list of species modes into a set of template macro definitions.
        Generated options are specific to the Tau Hybrid solver,
        and get passed as custom definitions to the build engine.

        :param model: Sanitized model containing sanitized species definitions.
        The GPY_HYBRID_SPECIES_MODES option will be set as an option for the model.
        :type model: SanitizedModel

        :returns: Pass-through of sanitized model object.
        :rtype: SanitizedModel
        """
        species_mode_map = {
            "continuous": "CONTINUOUS_MODE",
            "discrete": "DISCRETE_MODE",
            "dynamic": "DYNAMIC_MODE",
        }
        boundary_condition_types = [
            # When species.boundary_condition == False
            "STANDARD",
            # When species.boundary_condition == True
            "BOUNDARY",
        ]

        species_mode_list = []
        for spec_id, species in enumerate(model.species.values()):
            mode_keyword = species_mode_map.get(species.mode, species_mode_map["dynamic"])
            # Casting a bool to an int evaluates: False -> 0, and True -> 1
            # Explicit cast to bool for safety, in case boundary_condition is given weird values
            boundary_keyword = boundary_condition_types[int(bool(species.boundary_condition))]
            # Example: SPECIES_MODE(2, 10, CONTINUOUS_MODE, BOUNDARY)
            entry = f"SPECIES_MODE({spec_id},{species.switch_min},{mode_keyword},{boundary_keyword})"
            species_mode_list.append(entry)

        model.options["GPY_HYBRID_SPECIES_MODES"] = " ".join(species_mode_list)
        return model

    def _build(self, model: "Union[Model, SanitizedModel]", simulation_name: str, variable: bool, debug: bool = False,
               custom_definitions=None) -> str:
        sanitized_model = TauHybridCSolver.__create_options(SanitizedModel(model))
        for rate_rule in model.listOfRateRules.values():
            sanitized_model.use_rate_rule(rate_rule)
        return super()._build(sanitized_model, simulation_name, variable, debug)

    def get_solver_settings(self):
        """
        :return: Tuple of strings, denoting all keyword argument for this solvers run() method.
        """
        return ('model', 't', 'number_of_trajectories', 'timeout', 'increment', 'seed', 'debug', 'profile')

    def run(self=None, model: Model = None, t: int = 20, number_of_trajectories: int = 1, timeout: int = 0,
            increment: int = None, seed: int = None, debug: bool = False, profile: bool = False, variables={}, 
            resume=None, live_output: str = None, live_output_options: dict = {}, tau_step: int = .03, tau_tol=0.03, **kwargs):

        if self is None or self.model is None:
            self = TauHybridCSolver(model, resume=resume)

        increment = self.get_increment(model=model, increment=increment)

        # Validate parameters prior to running the model.
        self._validate_type(variables, dict, "'variables' argument must be a dictionary.")
        self._validate_variables_in_set(variables, self.species + self.parameters)
        self._validate_resume(t, resume)
        self._validate_kwargs(**kwargs)
        self._validate_sbml_features({
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
            "end": t,
            "interval": str(number_timesteps),
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

        if live_output is not None:
            live_output_options['type'] = live_output
            display_args = {
                "model": model, "number_of_trajectories": number_of_trajectories, "timeline": np.linspace(0, t, number_timesteps),
                "live_output_options": live_output_options, "resume": bool(resume)
            }
        else:
            display_args = None

        args = self._make_args(args)
        decoder = IterativeSimDecoder.create_default(number_of_trajectories, number_timesteps, len(self.model.listOfSpecies))

        sim_exec = self._build(model, self.target, self.variable, False)
        sim_status = self._run(sim_exec, args, decoder, timeout, display_args)

        if sim_status == SimulationReturnCode.FAILED:
            raise gillespyError.ExecutionError("Error encountered while running simulation C++ file:\n"
                f"Return code: {int(sim_status)}.\n")

        trajectories, time_stopped = decoder.get_output()

        simulation_data = self._format_output(trajectories)
        if sim_status == SimulationReturnCode.PAUSED:
            simulation_data = self._make_resume_data(time_stopped, simulation_data, t)
        if resume is not None:
            simulation_data = self._update_resume_data(resume, simulation_data, time_stopped)
        self.result = simulation_data
        self.rc = int(sim_status)

        return Results.build_from_solver_results(self, live_output_options)
