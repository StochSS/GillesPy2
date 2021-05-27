from collections import OrderedDict
from gillespy2.core import Species, Reaction, Parameter, Model

def write_template(path: str, model: Model, variable=False):
    """
    Write template definitions to a specified path. If the file does not exist, one is created.
    Definitions are written in `#define KEY VALUES` format.

    :param model: Model to get the macro definitions of.
    :type model: gillespy2.Model

    :param variable: Set to true to allow for non-constant parameter values.
    :type variable: bool
    """

    # Get a dictionary of model defines and transform into a list of strings in
    # `#define KEY VALUE` format.
    defines = get_model_defines(model, variable)
    template_lines = [(f"#define {key} {value}\n") for key, value in defines.items()]

    # Write generated lines to the template file.
    with open(path, "w") as template_file:
        template_file.writelines(template_lines)

def get_model_defines(model: Model, variable=False) -> "dict[str, str]":
    """
    Creates a dictionary of C++ macro definitions from the given model.
    The keys of the dictionary contain the name of the macro definition.
    The values of the dictionary are the values their corresponding macro should be defined to.

    :param model: Model to get the macro definitions of.
    :type model: gillespy2.Model

    :param variable: Set to true to allow for non-constant parameter values.
    :type variable: bool

    :return: Dictionary of fully-formatted macro definitions.
    """
    results = {}

    # Get definitions for variables
    parameter_mappings = model.sanitized_parameter_names()
    parameters = []
    parameter_names = []
    # System volume needs to be added manually.
    # It is assumed to ALWAYS be the first parameter.
    if model.volume:
        parameters.append(Parameter(name="V", expression=model.volume))
        parameter_names.append("V")
    for p_name, param in model.get_all_parameters().items():
        parameters.append(param)
        parameter_names.append(parameter_mappings[p_name])
    parameter_definitions = template_def_variables(parameters, parameter_names, variable)
    results.update(parameter_definitions)

    # Get definitions for species
    species = []
    species_names = []
    species_map = OrderedDict()
    for spec_id, spec_entry in enumerate(model.get_all_species().items()):
        spec_name, spec = spec_entry
        species.append(spec)
        species_names.append(spec_name)
        species_map[spec_name] = spec_id
    species_definitions = template_def_species(species, species_names)
    results.update(species_definitions)

    # Get definitions for reactions
    reaction_bases = model.get_all_reactions()
    reaction_names = []
    reactions = []
    for rxn_name, rxn in reaction_bases.items():
        reaction_names.append(rxn_name)
        reactions.append(rxn)
    reaction_definitions = template_def_reactions(reactions, reaction_names, species_map)
    results.update(reaction_definitions)
    
    # Get definitions for propensities
    stoch_propensities = []
    for rxn in reactions:
        stoch_propensity = rxn.sanitized_propensity_function(
            model.sanitized_species_names(),
            model.sanitized_parameter_names())
        stoch_propensities.append(stoch_propensity)
    stoch_propensity_definitions = template_def_propensities(stoch_propensities)
    results.update(stoch_propensity_definitions)

    return results

def template_def_variables(parameters: "list[Parameter]", sanitized_names: "list[str]", variable=False) -> "dict[str, str]":
    """
    Formats the relevant parameters to be passed to a C++ simulation template.
    Passed dictionaries/lists are assumed to be sanitized and sorted.
    
    :param parameters: Ordered list of runtime parameters.
    :type parameters: list[gillespy2.Parameter]

    :param sanitized_names: Ordered list of names for their corresponding runtime parameters.
    sanitized_names[i] should match parameters[i].
    :type sanitized_names: list[str]

    Returns the result as a list of tuples containing key-value pairs to be templated.
    """
    # Entries get defined as constant if variable is set to true.
    parameter_type = "VARIABLE" if variable else "CONSTANT"
    # Parameter entries, parsed and formatted
    parameter_set = []
    for param_id, parameter in enumerate(parameters):
        name = sanitized_names[param_id]

        parameter_set.append(f"{parameter_type}({name},{parameter.value})")

    return {
        "GPY_PARAMETER_VALUES": " ".join(parameter_set)
    }

def template_def_species(species: "list[Species]", sanitized_names: "list[str]") -> "dict[str, str]":
    """
    Passed dictionaries/lists are assumed to be sanitized and sorted.
    Formats the relevant species data to be passed to a C++ simulation template.

    :param species: Ordered list of species.
    :type species: gillespy2.Species

    :param sanitized_names: Ordered list of names corresponding to species.
    sanitized_names[i] should map to species[i].
    :type sanitized_names: list[str]

    :return: Dictionary of macro definitions for species and data related to species.
    """
    # Parse and format species initial populations
    populations = [str(specimen.initial_value) for specimen in species]
    num_species = str(len(populations))

    # Species names, parsed and formatted
    sanitized_names = [f"SPECIES_NAME({name})" for name in sanitized_names]
    populations = f"{{{','.join(populations)}}}"

    # Match each parameter with its macro definition name
    return {
        "GPY_INIT_POPULATIONS": populations,
        "GPY_NUM_SPECIES": num_species,
        "GPY_SPECIES_NAMES": " ".join(sanitized_names)
    }

def template_def_reactions(reactions: "list[Reaction]", sanitized_names: "list[str]", species_map: "OrderedDict[str, int]") -> "dict[str, str]":
    """
    Passed dictionaries/lists are assumed to be sanitized and sorted.
    Formats the relevant reactions and propensities to be passed to a C++ simulation template.

    :param reactions: Ordered list of reactions.
    The reaction's index in this list should correspond to its reaction id.
    For example, the reaction at reactions[3] has id 3.
    :type reactions: list[Reaction]

    :param sanitized_names: Ordered list of sanitized names for the reactions.
    The name's index in this list should match its corresponding reaction.
    sanitized_names[i] is matched to reactions[i].
    :type sanitized_names: list[str]

    :param species_map: Ordered dictionary mapping an unsanitized species name to its id.
    :type species_map: OrderedDict[str, int]

    :return: Dictionary of macro definitions for reactions.
    """
    num_reactions = str(len(reactions))
    reaction_set = []
    reaction_names = []

    for rxn_id, reaction in enumerate(reactions):
        name = sanitized_names[rxn_id]
        reaction_names.append(f"REACTION_NAME({name})")
        # Get stoichiometry matrix of reaction, turn into string, append to reacton_set
        species_change = [0] * len(species_map)

        for reactant, stoich in reaction.reactants.items():
            spec_id = species_map[reactant.name]
            species_change[spec_id] -= int(stoich)

        for product, stoich in reaction.products.items():
            spec_id = species_map[product.name]
            species_change[spec_id] += int(stoich)

        # Format the species changes as a stoichiometry set
        species_change = [str(dx) for dx in species_change]
        stoich = ",".join(species_change)
        reaction_set.append(f"{{{stoich}}}")

    reaction_set = f"{{{','.join(reaction_set)}}}"
    reaction_names = " ".join(reaction_names)

    return {
        "GPY_NUM_REACTIONS": num_reactions,
        "GPY_REACTIONS": reaction_set,
        "GPY_REACTION_NAMES": reaction_names
    }

def template_def_propensities(sanitized_propensities: "list[str]", ode=False) -> "dict[str, str]":
    """
    Formats the given list of pre-sorted, pre-sanitized propensities.

    :param sanitized_propensities: Ordered list of strings containing propensity functions.
    :type sanitized_propensities: list[str]

    :param ode: Boolean indicating whether the provided propensities are stochastic or deterministic.
    If set to True, then propensities will be assumed to be ODE propensities.
    :type ode: bool

    :return: Dictionary containing propensity macro definitions.
    """
    def_keyword = "ODE_PROPENSITIES" if ode else "PROPENSITIES"
    propensities = []

    for prop_id, prop in enumerate(sanitized_propensities):
        propensities.append(f"PROPENSITY({prop_id},{prop})")

    return {
        f"GPY_{def_keyword}": " ".join(propensities)
    }
