# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2023 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import importlib

from gillespy2.core.gillespySolver import GillesPySolver


class DynamicCSolver(GillesPySolver):
    """
    Abstract class for an embedded C solver.

    :param lib_target: The member of the cgillespy library method to invoke
    :type lib_target: str
    """

    name = "DynamicCSolver"

    def __init__(self, lib_target: "str"):
        libcgillespy = importlib.import_module("libcgillespy")
        self.lib_target = getattr(libcgillespy, lib_target)

    def run(self, model, *args, **kwargs):
        deconstructed_species = {spec_name: float(spec.initial_value)
                                 for spec_name, spec in model.listOfSpecies.items()}
        deconstructed_reactions = {rx_name: {
            "reactants": [
                rx.reactants.get(spec, 0)
                for spec in model.listOfSpecies.values()
            ],
            "products": [
                rx.products.get(spec, 0)
                for spec in model.listOfSpecies.values()
            ],
            "rate_index": 0, # todo
        } for rx_name, rx in model.listOfReactions.items()}
        deconstructed_parameters = {param_key: float(eval(param.expression)) for param_key, param in model.listOfParameters.items()}
        result = self.lib_target(species=deconstructed_species,
                                 reactions=deconstructed_reactions,
                                 constants={},
                                 variables=deconstructed_parameters)
        # TODO: once return is implemented, convert to simulation trajectories object to conform to solver interface
        return result
