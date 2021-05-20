from typing import OrderedDict
from gillespy2.core import Species, Reaction, Parameter, Model

def template_def_variables(parameters: list[Parameter], sanitized_names = list[str], variable=False):
    """
    Formats the relevant parameters to be passed to a C++ simulation template.
    Passed dictionaries/lists are assumed to be sanitized and sorted,
      such that parameters[i] should match sanitized_names[i].

    Returns the result as a list of tuples containing key-value pairs to be templated.
    """
    # Entries get defined as constant if variable is set to true.
    parameter_type = "VARIABLE" if variable else "CONSTANT"
    # Parameter entries, parsed and formatted
    parameter_set = []
    for param_id, parameter in enumerate(parameters):
        name = sanitized_names[param_id]
        parameter_value = "{0}({1},{2})".format(parameter_type, name, parameter.value)
        parameter_set.append(parameter_value)

    return [
        ("GPY_PARAMETER_VALUES", " ".join(parameter_set))
    ]

def template_def_species(species: list[Species], sanitized_names: list[str]):
    """
    Passed dictionaries/lists are assumed to be sanitized and sorted.
    species[i] should map to sanitized_names[i].
    Formats the relevant species data to be passed to a C++ simulation template.

    Returns the result as a list of tuples containing key-value pairs to be templated.
    """
    # Parse and format species initial populations
    populations = [str(specimen.initial_value) for specimen in species]

    # Species names, parsed and formatted
    species_names = "{{{0}}}".format(",".join(sanitized_names))
    populations = "{{{0}}}".format(",".format(populations))
    num_species = str(len(populations))

    # Match each parameter with its macro definition name
    return [
        ("GPY_INIT_POPULATIONS", populations),
        ("GPY_NUM_SPECIES", num_species),
        ("GPY_SPECIES_NAMES", species_names),
    ]

def template_def_reactions(reactions: list[Reaction], sanitized_names: list[str], species_map: OrderedDict[str, Species]):
    """
    Passed dictionaries/lists are assumed to be sanitized and sorted.
    Formats the relevant reactions and propensities to be passed to a C++ simulation template.

    Returns the result as a list of tuples containing key-value pairs to be templated.
    """
    # TODO
    num_reactions = str(0)
    reaction_set = []
    reaction_names = []

    for rxn_id, reaction in enumerate(reactions):
        name = sanitized_names[rxn_id]
        reaction_names.append(name)
        # TODO: get stoichiometry matrix of reaction, turn into string, append to reacton_set
        species_change = ""

        for reactant in reaction.reactants:
            pass

        for product in reaction.products:
            pass

        reaction_set.append(species_change)

    reaction_set = "{{{0}}}".format(",".join(reaction_set))
    reaction_names = "{{{0}}}".format(",".join(reaction_names))

    return [
        ("GPY_NUM_REACTIONS", num_reactions),
        ("GPY_REACTIONS", reaction_set),
        ("GPY_REACTION_NAMES", reaction_names),
    ]
