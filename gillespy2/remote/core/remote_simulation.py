'''
gillespy2.remote.core.remote_simulation
'''
# StochSS-Compute is a tool for running and caching GillesPy2 simulations remotely.
# Copyright (C) 2019-2023 GillesPy2 and StochSS developers.

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

from gillespy2.remote.client.endpoint import Endpoint
from gillespy2.remote.core.messages.simulation_run import SimulationRunRequest, SimulationRunResponse
from gillespy2.remote.core.messages.simulation_run_unique import SimulationRunUniqueRequest, SimulationRunUniqueResponse
from gillespy2.remote.core.messages.status import SimStatus
from gillespy2.remote.core.errors import RemoteSimulationError
from gillespy2.remote.core.remote_results import RemoteResults

class RemoteSimulation:
    '''
    An object representing a remote gillespy2 simulation. Requires a model and a host address.
    A solver type may be provided, but does not accept instantiated solvers.

    :param model: The model to simulate.
    :type model: gillespy2.Model

    :param server: A server to run the simulation. Optional if host is provided.
    :type server: gillespy2.remote.Server

    :param host: The address of a running instance of StochSS-Compute. Optional if server is provided.
    :type host: str

    :param port: The port to use when connecting to the host.
                 Only needed if default server port is changed. Defaults to 29681.
    :type port: int

    :param solver: The type of solver to use. Does not accept instantiated solvers.
    :type solver: gillespy2.GillesPySolver
    '''
    def __init__(self,
                 model,
                 server = None,
                 host: str = None,
                 port: int = 29681,
                 solver = None,
                 ):

        if server is not None and host is not None:
            raise RemoteSimulationError('Pass a ComputeServer/Cluster object or host but not both.')
        if server is None and host is None:
            raise RemoteSimulationError('Pass a ComputeServer/Cluster object or host.')
        if server is None and port is None:
            raise RemoteSimulationError('Pass a ComputeServer/Cluster object or port.')

        if server is None:
            from gillespy2.remote.client.compute_server import ComputeServer
            self.server = ComputeServer(host, port)
        else:
            self.server = server

        self.model = model

        if solver is not None:
            if hasattr(solver, 'is_instantiated'):
                raise RemoteSimulationError(
                    'RemoteSimulation does not accept an instantiated solver object. Pass a type.')
        self.solver = solver

    def is_cached(self, **params):
        '''
        Checks to see if a dummy simulation exists in the cache.

        :param params: Arguments for simulation.
        :type params: dict[str, Any]

        :returns: If the results are cached on the server.
        :rtype: bool
        '''
        if "solver" in params:
            if hasattr(params['solver'], 'is_instantiated'):
                raise RemoteSimulationError(
                    'RemoteSimulation does not accept an instantiated solver object. Pass a type.')
            params["solver"] = f"{params['solver'].__module__}.{params['solver'].__qualname__}"
        if self.solver is not None:
            params["solver"] = f"{self.solver.__module__}.{self.solver.__qualname__}"

        sim_request = SimulationRunRequest(model=self.model, **params)
        results_dummy = RemoteResults()
        results_dummy.id = sim_request.hash()
        results_dummy.server = self.server
        results_dummy.n_traj = params.get('number_of_trajectories', 1)
        return results_dummy.is_ready

    def run(self, ignore_cache=False, **params):
        # pylint:disable=line-too-long
        """
        Simulate the Model on the target ComputeServer, returning the results or a handle to a running simulation.
        See `here <https://stochss.github.io/GillesPy2/docs/build/html/classes/gillespy2.core.html#gillespy2.core.model.Model.run>`_.

        :param unique: When True, ignore cache completely and return always new results.
        :type unique: bool

        :param params: Arguments to pass directly to the Model#run call on the server.
        :type params: dict[str, Any]

        :returns: RemoteResults populated with Results if cached, otherwise and unpopulated RemoteResults
        :rtype: RemoteResults

        :raises RemoteSimulationError: In the case of SimStatus.ERROR
        """
        # pylint:enable=line-too-long

        if "solver" in params:
            if hasattr(params['solver'], 'is_instantiated'):
                raise RemoteSimulationError(
                    'RemoteSimulation does not accept an instantiated solver object. Pass a type.')
            params["solver"] = f"{params['solver'].__module__}.{params['solver'].__qualname__}"
        if self.solver is not None:
            params["solver"] = f"{self.solver.__module__}.{self.solver.__qualname__}"
        if ignore_cache is True:
            sim_request = SimulationRunUniqueRequest(self.model, **params)
            return self._run_unique(sim_request)
        if ignore_cache is False:
            sim_request = SimulationRunRequest(self.model, **params)
            return self._run(sim_request)

    def _run(self, request):
        '''
        :param request: Request to send to the server. Contains Model and related arguments.
        :type request: SimulationRunRequest
        '''
        response_raw = self.server.post(Endpoint.SIMULATION_GILLESPY2, sub="/run", request=request)
        if not response_raw.ok:
            raise Exception(response_raw.reason)

        sim_response = SimulationRunResponse.parse(response_raw.text)

        if sim_response.status == SimStatus.ERROR:
            raise RemoteSimulationError(sim_response.error_message)
        if sim_response.status == SimStatus.READY:
            remote_results =  RemoteResults(data=sim_response.results.data)
        else:
            remote_results =  RemoteResults()

        remote_results.id = sim_response.results_id
        remote_results.server = self.server
        remote_results.n_traj = request.kwargs.get('number_of_trajectories', 1)
        remote_results.task_id = sim_response.task_id

        return remote_results

    def _run_unique(self, request):
        '''
        Ignores the cache. Gives each simulation request a unique identifier.

        :param request: Request to send to the server. Contains Model and related arguments.
        :type request: SimulationRunUniqueRequest
        '''
        response_raw = self.server.post(Endpoint.SIMULATION_GILLESPY2, sub="/run/unique", request=request)

        if not response_raw.ok:
            raise Exception(response_raw.reason)
        sim_response = SimulationRunUniqueResponse.parse(response_raw.text)
        if not sim_response.status is SimStatus.RUNNING:
            raise Exception(sim_response.error_message)
        # non-conforming object creation ... possible refactor needed to solve, so left in.
        remote_results =  RemoteResults()
        remote_results.id = request.unique_key
        remote_results.task_id = request.unique_key
        remote_results.server = self.server
        remote_results.n_traj = request.kwargs.get('number_of_trajectories', 1)

        return remote_results
