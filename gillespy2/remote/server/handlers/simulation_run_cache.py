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
import random

from tornado.web import RequestHandler
from tornado.ioloop import IOLoop
from distributed import Client
from gillespy2.remote.core.messages.status import SimStatus
from gillespy2.remote.core.messages.simulation_run_cache import SimulationRunCacheRequest, SimulationRunCacheResponse
from gillespy2.remote.server.cache import Cache

from gillespy2.remote.core.utils.log_config import init_logging
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
        log.debug(sim_request.encode())
        log.debug('%(namespace)s', locals())

        results_id = sim_request.results_id
        cache = Cache(self.cache_dir, results_id, namespace=sim_request.namespace)
        msg_0 = f'<{self.request.remote_ip}> | <{results_id}>'
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
                self._run_cache(sim_request)
            else:
                msg = f'{msg_0} | Returning cached results.'
                log.info(msg)
                results = cache.get_sample(n_traj)
                results_json = results.to_json()
                sim_response = SimulationRunCacheResponse(SimStatus.READY, results_id = results_id, results = results_json)
                self.write(sim_response.encode())
                self.finish()
        if empty:
            msg = f'{msg_0} | Results not cached. Running simulation.'
            log.info(msg)
            self._run_cache(sim_request)

    def _run_cache(self, sim_request) -> None:
        '''
        :param results_id: Key to results.
        :type results_id: str

        :param future: Future that completes to gillespy2.Results.
        :type future: distributed.Future

        :param client: Handle to dask scheduler.
        :type client: distributed.Client

        '''
        results_id = sim_request.results_id
        client = Client(self.scheduler_address)
        log.debug('_run_cache():')
        log.debug(sim_request.parallelize)
        log.debug(sim_request.chunk_trajectories)
        if sim_request.parallelize is True:
            self._submit_parallel(sim_request, client)
        elif sim_request.chunk_trajectories is True:
            self._submit_chunks(sim_request, client)
        else:
            future = self._submit(sim_request, client)
            self._return_running(results_id, future.key)
            IOLoop.current().run_in_executor(None, self._cache, results_id, future, client)

    def _submit_parallel(self, sim_request, client: Client):
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
        results_id = sim_request.results_id
        if "solver" in kwargs:
            from pydoc import locate
            kwargs["solver"] = locate(kwargs["solver"])

        futures = []
        n_traj = kwargs['number_of_trajectories']
        del kwargs['number_of_trajectories']
        kwargs['seed'] = int(sim_request.id, 16)
        for i in range(n_traj):
            kwargs['seed'] += kwargs['seed']
            future = client.submit(model.run, **kwargs)
            futures.append(future)
        IOLoop.current().run_in_executor(None, self._cache_parallel, results_id, futures, client)
        self._return_running(results_id, results_id)

    def _submit_chunks(self, sim_request, client: Client):
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
        results_id = sim_request.results_id
        if "solver" in kwargs:
            from pydoc import locate
            kwargs["solver"] = locate(kwargs["solver"])
        n_traj = kwargs.get('number_of_trajectories',1)
        n_workers =  len(client.scheduler_info()['workers'])
        traj_per_worker = n_traj // n_workers
        extra_traj = n_traj % n_workers
        log.debug('_submit_chunks():')
        log.debug(traj_per_worker)
        log.debug(extra_traj)
        # return
        # traj_per_worker_list = []
        futures = []
        del kwargs['number_of_trajectories']
        kwargs['seed'] = int(sim_request.id, 16)
        for _ in range(n_workers):
            kwargs['seed'] += kwargs['seed']
            kwargs['number_of_trajectories'] = traj_per_worker
            future = client.submit(model.run, **kwargs)
            futures.append(future)
        if extra_traj > 0:
            kwargs['seed'] += kwargs['seed']
            kwargs['number_of_trajectories'] = extra_traj
            future = client.submit(model.run, **kwargs)
            futures.append(future)

        IOLoop.current().run_in_executor(None, self._cache_parallel, results_id, futures, client)
        self._return_running(results_id, results_id)
        


    def _cache_parallel(self, results_id, futures, client: Client) -> None:
        '''
        :param results_id: Key to results.
        :type results_id: str

        :param future: List of Futures that completes to gillespy2.Results.
        :type : List[distributed.Future]

        :param client: Handle to dask scheduler.
        :type client: distributed.Client

        '''
        results = client.gather(futures, asynchronous=False)
        results = sum(results)
        log.debug('_cache_parallel():')
        log.debug(results)
        client.close()
        # return
        cache = Cache(self.cache_dir, results_id)
        cache.save(results)
        
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

    def _submit(self, sim_request, client: Client):
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
