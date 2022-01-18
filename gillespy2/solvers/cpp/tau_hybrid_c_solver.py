import gillespy2
from gillespy2.solvers.cpp.c_decoder import IterativeSimDecoder
from gillespy2.solvers.utilities import solverutils as cutils
from gillespy2.core import GillesPySolver, Model
from gillespy2.core.gillespyError import *
from typing import Union
from gillespy2.core import Results

from .c_solver import CSolver, SimulationReturnCode
from gillespy2.solvers.cpp.build.template_gen import SanitizedModel

class TauHybridCSolver(GillesPySolver, CSolver):
    name = "TauHybridCSolver"
    target = "hybrid"

    @classmethod
    def __create_options(cls, sanitized_model: "SanitizedModel") -> "SanitizedModel":
        """
        Populate the given list of species modes into a set of template macro definitions.
        Generated options are specific to the Tau Hybrid solver,
        and get passed as custom definitions to the build engine.

        :param sanitized_model: Sanitized model containing sanitized species definitions.
        The GPY_HYBRID_SPECIES_MODES option will be set as an option for the model.
        :type sanitized_model: SanitizedModel

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
        trigger_mode_types = [
            # When event.use_values_from_trigger_time == False
            "USE_EVAL",
            # When event.use_values_from_trigger_time == True
            "USE_TRIGGER",
        ]
        persist_types = [
            # When event.trigger.persistent == False
            "IRREGULAR",
            # When event.trigger.persistent == True
            "PERSISTENT",
        ]
        initial_value_types = [
            # When event.trigger.initial_value == False
            "INIT_FALSE",
            # When event.trigger.initial_value == True
            "INIT_TRUE",
        ]

        species_mode_list = []
        for spec_id, species in enumerate(sanitized_model.species.values()):
            mode_keyword = species_mode_map.get(species.mode, species_mode_map["dynamic"])
            # Casting a bool to an int evaluates: False -> 0, and True -> 1
            # Explicit cast to bool for safety, in case boundary_condition is given weird values
            boundary_keyword = boundary_condition_types[int(bool(species.boundary_condition))]
            # Example: SPECIES_MODE(2, 10, CONTINUOUS_MODE, BOUNDARY)
            entry = f"SPECIES_MODE({spec_id},{species.switch_min},{mode_keyword},{boundary_keyword})"
            species_mode_list.append(entry)

        # EVENT(event_id, {targets}, trigger, delay, priority, use_trigger, use_persist)
        event_list = []
        # [SPECIES/VARIABLE]_ASSIGNMENT(assign_id, target_id, expr)
        event_assignment_list = []
        assign_id = 0
        for event_id, event in enumerate(sanitized_model.model.listOfEvents.values()):
            trigger = sanitized_model.expr.getexpr_cpp(event.trigger.expression)
            delay = sanitized_model.expr.getexpr_cpp(event.delay) \
                if event.delay is not None else "0"
            priority = sanitized_model.expr.getexpr_cpp(event.priority) \
                if event.priority is not None else "0"
            use_trigger = trigger_mode_types[int(bool(event.use_values_from_trigger_time))]
            use_persist = persist_types[int(bool(event.trigger.persistent))]
            initial_value = initial_value_types[int(bool(event.trigger.value or False))]

            assignments: "list[str]" = []
            for assign in event.assignments:
                variable = assign.variable
                expression = sanitized_model.expr.getexpr_cpp(assign.expression)

                if isinstance(variable, str):
                    if variable in sanitized_model.model.listOfSpecies:
                        variable = sanitized_model.model.listOfSpecies.get(variable)
                    elif variable in sanitized_model.model.listOfParameters:
                        variable = sanitized_model.model.listOfParameters.get(variable)
                    else:
                        raise ValueError(f"Invalid event assignment {assign}: received name {variable} "
                                         f"Must match the name of a valid Species or Parameter.")

                if isinstance(variable, gillespy2.Species):
                    assign_str = f"SPECIES_ASSIGNMENT(" \
                                 f"{assign_id},{sanitized_model.species_id.get(variable.name)},{expression})"
                elif isinstance(variable, gillespy2.Parameter):
                    assign_str = f"VARIABLE_ASSIGNMENT(" \
                                 f"{assign_id},{sanitized_model.parameter_id.get(variable.name)},{expression})"
                else:
                    raise ValueError(f"Invalid event assignment {assign}: received variable of type {type(variable)} "
                                     f"Must be of type str, Species, or Parameter")
                assignments.append(str(assign_id))
                event_assignment_list.append(assign_str)
                assign_id += 1
            # Check for "None"s
            for a in assignments:
                if a is None: raise Exception(f"assignment={a} is None in event={event}")
            if event_id is None: raise Exception(f"event_id is None in event={event}")
            if trigger is None: raise Exception(f"trigger is None in event={event}")
            if delay is None: raise Exception(f"delay is None in event={event}")
            if priority is None: raise Exception(f"priority is None in event={event}")
            if use_trigger is None: raise Exception(f"use_trigger is None in event={event}")
            if use_persist is None: raise Exception(f"use_persist is None in event={event}")
            if initial_value is None: raise Exception(f"initial_value is None in event={event}")

            assignments: "str" = " AND ".join(assignments)
            event_list.append(
                f"EVENT("
                f"{event_id},"
                f"{{{assignments}}},"
                f"{trigger},"
                f"{delay},"
                f"{priority},"
                f"{use_trigger},"
                f"{use_persist},"
                f"{initial_value}"
                f")"
            )

        sanitized_model.options["GPY_HYBRID_SPECIES_MODES"] = " ".join(species_mode_list)
        sanitized_model.options["GPY_HYBRID_EVENTS"] = " ".join(event_list)
        sanitized_model.options["GPY_HYBRID_NUM_EVENTS"] = str(len(event_list))
        sanitized_model.options["GPY_HYBRID_EVENT_ASSIGNMENTS"] = " ".join(event_assignment_list)
        sanitized_model.options["GPY_HYBRID_NUM_EVENT_ASSIGNMENTS"] = str(len(event_assignment_list))
        return sanitized_model

    def _build(self, model: "Union[Model, SanitizedModel]", simulation_name: str, variable: bool, debug: bool = False,
               custom_definitions=None) -> str:
        variable = variable or len(model.listOfEvents) > 0
        sanitized_model = TauHybridCSolver.__create_options(SanitizedModel(model, variable=variable))
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

        if self is None:
            self = TauHybridCSolver(model, resume=resume)
        if self.model is None:
            if model is None:
                raise SimulationError("A model is required to run the simulation.")
            self._set_model(model=model)
        if model is not None and model.get_json_hash() != self.model.get_json_hash():
            raise SimulationError("Model must equal TauHybridCSolver.model.")
        self.model.resolve_parameters()

        increment = self.get_increment(increment=increment)

        # Validate parameters prior to running the model.
        self._validate_type(variables, dict, "'variables' argument must be a dictionary.")
        self._validate_variables_in_set(variables, self.species + self.parameters)
        self._validate_resume(t, resume)
        self._validate_kwargs(**kwargs)
        self._validate_sbml_features({
            "Assignment Rules": len(self.model.listOfAssignmentRules),
            "Function Definitions": len(self.model.listOfFunctionDefinitions)
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
            populations = cutils.update_species_init_values(self.model.listOfSpecies, self.species, variables, resume)
            parameter_values = cutils.change_param_values(self.model.listOfParameters, self.parameters, self.model.volume, variables)

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
                "model": self.model, "number_of_trajectories": number_of_trajectories, "timeline": np.linspace(0, t, number_timesteps),
                "live_output_options": live_output_options, "resume": bool(resume)
            }
        else:
            display_args = None

        args = self._make_args(args)
        decoder = IterativeSimDecoder.create_default(number_of_trajectories, number_timesteps, len(self.model.listOfSpecies))

        sim_exec = self._build(self.model, self.target, self.variable, False)
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
