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
from gillespy2.remote.core.messages.simulation_run_cache import SimulationRunCacheRequest, SimulationRunCacheResponse
from gillespy2.remote.core.messages.status import SimStatus, StatusRequest
from gillespy2.remote.core.errors import RemoteSimulationError
from gillespy2.remote.core.remote_results import RemoteResults

from gillespy2.remote.core.utils.log_config import init_logging
log = init_logging(__name__)

class RemoteSimulation:
    '''
    An object representing a remote gillespy2 simulation. Requires a model and a host address.
    A solver type may be provided, but does not accept instantiated solvers.

    :param model: The model to simulate.
    :type model: gillespy2.Model

    :param server: A server to run the simulation. Optional if host is provided.
    :type server: gillespy2.remote.Server

    :param host: The address of a running instance of gillespy2.remote. Optional if server is provided.
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

    def is_cached(self, namespace=None, **params):
        '''
        Checks to see if a dummy simulation exists in the cache.

        :param namespace: If provided, prepend to results path.
        :type namespace: str

        :param params: Arguments for simulation.
        :type params: dict[str, Any]

        :returns: If the results are cached on the server.
        :rtype: bool
        '''
        log.debug('is_cached()')
        n_traj_requested = params.get('number_of_trajectories', 1)
        n_traj_cached = self.get_n_traj_in_cache(namespace=namespace,
                                                 **params)
        if n_traj_requested <= n_traj_cached:
            return True
        else:
            return False

    def run(self, namespace=None, ignore_cache=False, parallelize=False, chunk_trajectories=False, **params):
        # pylint:disable=line-too-long
        """
        Simulate the Model on the target ComputeServer, returning the results or a handle to a running simulation.
        See `here <https://stochss.github.io/GillesPy2/docs/build/html/classes/gillespy2.core.html#gillespy2.core.model.Model.run>`_.

        :param namespace: If provided, prepend to results path.
        :type namespace: str

        :param ignore_cache: When True, ignore cache completely and return always new results.
        :type ignore_cache: bool

        :param params: Arguments to pass directly to the Model#run call on the server.
        :type params: dict[str, Any]

        :returns: RemoteResults populated with Results if cached, otherwise and unpopulated RemoteResults
        :rtype: RemoteResults

        :raises RemoteSimulationError: In the case of SimStatus.ERROR
        """
        # pylint:enable=line-too-long

        params = self._encode_solver_arg(**params)
        sim_request = SimulationRunCacheRequest(self.model,
                                                namespace=namespace,
                                                ignore_cache=ignore_cache,
                                                parallelize=parallelize,
                                                chunk_trajectories=chunk_trajectories,
                                                **params)
        log.debug('run()...')
        log.debug('sim_request.results_id')
        log.debug(sim_request.results_id)
        return self._run_cache(sim_request)

    def _run_cache(self, request):
        '''
        :param request: Request to send to the server. Contains Model and related arguments.
        :type request: SimulationRunRequest
        '''
        log.debug('_run_cache(request)...')
        log.debug('request.kwargs.get("solver", None)')
        log.debug(request.kwargs.get('solver', None))
        response_raw = self.server.post(Endpoint.SIMULATION_GILLESPY2, sub='/run/cache', request=request)

        if not response_raw.ok:
            raise Exception(response_raw.reason)

        sim_response = SimulationRunCacheResponse.parse(response_raw.text)

        if sim_response.status == SimStatus.ERROR:
            raise RemoteSimulationError(sim_response.error_message)
        if sim_response.status == SimStatus.READY:
            data = sim_response.results.data
        else:
            data = None
        return RemoteResults.factory(request.results_id,
                                     self.server,
                                     request.kwargs.get('number_of_trajectories', 1),
                                     request.id,
                                     request.namespace,
                                     data)

    def run_distributed(self, namespace=None, force_run=False, **params):
        log.debug('run_distributed()...')
        log.debug('namespace:')
        log.debug(namespace)
        log.debug('force_run:')
        log.debug(force_run)
        log.debug('params.get("solver", None)')
        log.debug(params.get('solver', None))
        params = self._encode_solver_arg(**params)
        log.debug('params.get("solver", None)')
        log.debug(params.get('solver', None))
        response_raw = self.server.get(Endpoint.DASK, sub='/number_of_workers')
        n_workers = int(response_raw.text)
        log.debug('n_workers:')
        log.debug(n_workers)
        n_traj_requested = params.get('number_of_trajectories', 1)
        results_collection = []
        if force_run is True:
            n_traj = n_traj_requested
        if force_run is False:
            n_traj_remote = self.get_n_traj_in_cache(namespace=namespace,
                                                    **params)
            log.debug('n_traj_remote:')
            log.debug(n_traj_remote)
            n_traj = n_traj_requested - n_traj_remote
                
        traj_per_worker = n_traj // n_workers
        log.debug('traj_per_worker:')
        log.debug(traj_per_worker)
        extra_traj = n_traj % n_workers
        log.debug('extra_traj:')
        log.debug(extra_traj)
        if traj_per_worker > 0:
            for _ in range(n_workers):
                params['number_of_trajectories'] = traj_per_worker
                sim_request = SimulationRunCacheRequest(self.model,
                                                        namespace=namespace,
                                                        force_run=True,
                                                        **params)
                results_collection.append(self._run_cache(sim_request))
        if extra_traj > 0:
            params['number_of_trajectories'] = extra_traj
            sim_request = SimulationRunCacheRequest(self.model,
                                                    namespace=namespace,
                                                    force_run=True,
                                                    **params)
            results_collection.append(self._run_cache(sim_request))
        
        return RemoteResults.factory(sim_request.results_id,
                                     self.server,
                                     n_traj_requested,
                                     None,
                                     sim_request.namespace,
                                     None)

    def get_n_traj_in_cache(self, namespace=None, **params) -> int:
        log.debug('get_n_traj_in_cache()...')
        log.debug('namespace:')
        log.debug(namespace)
        params = self._encode_solver_arg(**params)
        simulation_request = SimulationRunCacheRequest(self.model,
                                   namespace=namespace,
                                   **params)
        log.debug('simulation_request.results_id:')
        log.debug(simulation_request.results_id)
        n_traj = params.get('number_of_trajectories', 1)
        status_request = StatusRequest(results_id=simulation_request.results_id,
                                       n_traj=n_traj,
                                       namespace=namespace)
        response_raw = self.server.get(Endpoint.CACHE, sub='/number_of_trajectories', request=status_request)
        n_traj_in_cache = int(response_raw.text)
        log.debug('n_traj_in_cache:')
        log.debug(n_traj_in_cache)
        return n_traj_in_cache
    
    def _encode_solver_arg(self, **params):
        log.debug('_encode_solver_arg()...')
        value = params.get('solver', '') # Makes for a better fail if they pass something invalid.
        log.debug('params.get("solver", None)')
        log.debug(params.get('solver', None))
        if type(value) is not str:
            if hasattr(params['solver'], 'is_instantiated'):
                raise RemoteSimulationError(
                    'RemoteSimulation does not accept an instantiated solver object. Pass a type.')
            params['solver'] = f"{value.__module__}.{value.__qualname__}"
        if self.solver is not None and value == '':
            params["solver"] = f"{self.solver.__module__}.{self.solver.__qualname__}"
        log.debug('params.get("solver", None)')
        log.debug(params.get('solver', None))
        return params
