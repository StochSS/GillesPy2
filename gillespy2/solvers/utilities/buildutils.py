from typing import OrderedDict
from gillespy2.core import Species, Reaction, Parameter, Model

def template_def_variables(parameters: list[Parameter], sanitized_names = list[str], variable=False) -> dict[str, str]:
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
    for param_id, parameter in parameters:
        name = sanitized_names[param_id]

        parameter_set.append(f"{parameter_type}({name},{parameter.value})")

    return {
        "GPY_PARAMETER_VALUES": " ".join(parameter_set)
    }

def template_def_species(species: list[Species], sanitized_names: list[str]) -> dict[str, str]:
    """
    Passed dictionaries/lists are assumed to be sanitized and sorted.
    species[i] should map to sanitized_names[i].
    Formats the relevant species data to be passed to a C++ simulation template.

    Returns the result as a list of tuples containing key-value pairs to be templated.
    """
    # Parse and format species initial populations
    populations = [str(specimen.initial_value) for specimen in species]

    # Species names, parsed and formatted
    species_names = f"{{{','.join(sanitized_names)}}}"
    populations = f"{{{','.join(populations)}}}"
    num_species = str(len(populations))

    # Match each parameter with its macro definition name
    return {
        "GPY_INIT_POPULATIONS": populations,
        "GPY_NUM_SPECIES": num_species,
        "GPY_SPECIES_NAMES": species_names
    }

def template_def_reactions(reactions: list[Reaction], sanitized_names: list[str], species_map: OrderedDict[str, int]) -> dict[str, str]:
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

    Returns the result as a list of tuples containing key-value pairs to be templated.
    """
    num_reactions = str(len(reactions))
    reaction_set = []
    reaction_names = []

    for rxn_id, reaction in enumerate(reactions):
        name = sanitized_names[rxn_id]
        reaction_names.append(name)
        # Get stoichiometry matrix of reaction, turn into string, append to reacton_set
        species_change = [0] * len(species_map)

        for reactant, stoich in reaction.reactants.items():
            spec_id = species_map[reactant.name]
            species_change[spec_id] -= stoich

        for product, stoich in reaction.products.items():
            spec_id = species_map[product.name]
            species_change[spec_id] += stoich

        # Format the species changes as a stoichiometry set
        species_change = [str(dx) for dx in species_change]
        stoich = ",".join(species_change)
        reaction_set.append(f"{{{stoich}}}")

    reaction_set = f"{{{','.join(reaction_set)}}}"
    reaction_names = f"{{{','.join(reaction_names)}}}"

    return {
        "GPY_NUM_REACTIONS": num_reactions,
        "GPY_REACTIONS": reaction_set,
        "GPY_REACTION_NAMES": reaction_names
    }
