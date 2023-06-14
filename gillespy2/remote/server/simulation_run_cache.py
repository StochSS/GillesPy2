'''
gillespy2.remote.server.run
'''
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

import os
from secrets import token_hex

from tornado.web import RequestHandler
from tornado.ioloop import IOLoop
from distributed import Client
from gillespy2.remote.core.messages.status import SimStatus
from gillespy2.remote.core.messages.simulation_run_cache import SimulationRunCacheRequest, SimulationRunCacheResponse
from gillespy2.remote.server.cache import Cache

from gillespy2.remote.core.log_config import init_logging
log = init_logging(__name__)

class SimulationRunCacheHandler(RequestHandler):
    '''
    Endpoint for running Gillespy2 simulations as normal,
    except that trajectories are cached and reused
    if a particular model is run multiple times.
    '''

    scheduler_address = None
    cache_dir = None

    def initialize(self, scheduler_address, cache_dir):
        '''
        Sets the address to the Dask scheduler and the cache directory.

        :param scheduler_address: Scheduler address.
        :type scheduler_address: str

        :param cache_dir: Path to the cache.
        :type cache_dir: str
        '''

        self.scheduler_address = scheduler_address
        self.cache_dir = cache_dir

    async def post(self):
        '''
        Process simulation run request.
        '''
        sim_request = SimulationRunCacheRequest.parse(self.request.body)
        namespace = sim_request.namespace
        log.debug('%(namespace)s', locals())
        if namespace != '':
            namespaced_dir = os.path.join(namespace, self.cache_dir)
            self.cache_dir = namespaced_dir
            log.debug(namespaced_dir)
        if sim_request.ignore_cache is True:
            cache = Cache(self.cache_dir, sim_request.id)
        else:
            cache = Cache(self.cache_dir, sim_request.results_id)

        sim_hash = sim_request.results_id
        msg_0 = f'<{self.request.remote_ip}> | <{sim_hash}>'
        if not cache.exists():
            cache.create()
        empty = cache.is_empty()
        if not empty:
            # Check the number of trajectories in the request, default 1
            n_traj = sim_request.kwargs.get('number_of_trajectories', 1)
            # Compare that to the number of cached trajectories
            trajectories_needed =  cache.n_traj_needed(n_traj)
            if trajectories_needed > 0:
                sim_request.kwargs['number_of_trajectories'] = trajectories_needed
                msg = f'{msg_0} | Partial cache. Running {trajectories_needed} new trajectories.'
                log.info(msg)
                client = Client(self.scheduler_address)
                future = self._submit(sim_request, client)
                self._return_running(sim_hash, future.key)
                IOLoop.current().run_in_executor(None, self._cache, sim_hash, future, client)
            else:
                msg = f'{msg_0} | Returning cached results.'
                log.info(msg)
                results = cache.get_sample(n_traj)
                results_json = results.to_json()
                sim_response = SimulationRunCacheResponse(SimStatus.READY, results_id = sim_hash, results = results_json)
                self.write(sim_response.encode())
                self.finish()
        if empty:
            msg = f'{msg_0} | Results not cached. Running simulation.'
            log.info(msg)
            client = Client(self.scheduler_address)
            future = self._submit(sim_request, client)
            self._return_running(sim_hash, future.key)
            IOLoop.current().run_in_executor(None, self._cache, sim_hash, future, client)

    def _cache(self, results_id, future, client) -> None:
        '''
        :param results_id: Key to results.
        :type results_id: str

        :param future: Future that completes to gillespy2.Results.
        :type future: distributed.Future

        :param client: Handle to dask scheduler.
        :type client: distributed.Client

        '''
        results = future.result()
        client.close()
        cache = Cache(self.cache_dir, results_id)
        cache.save(results)

    def _submit(self, sim_request, client):
        '''
        :param sim_request: Incoming request.
        :type sim_request: SimulationRunCacheRequest

        :param client: Handle to dask scheduler.
        :type client: distributed.Client

        :returns: Future that completes to gillespy2.Results.
        :rtype: distributed.Future
        '''
        model = sim_request.model
        kwargs = sim_request.kwargs

        if "solver" in kwargs:
            from pydoc import locate
            kwargs["solver"] = locate(kwargs["solver"])

        return client.submit(model.run, **kwargs, key=sim_request.id)

    def _return_running(self, results_id, task_id):
        sim_response = SimulationRunCacheResponse(SimStatus.RUNNING, results_id=results_id, task_id=task_id)
        self.write(sim_response.encode())
        self.finish()
