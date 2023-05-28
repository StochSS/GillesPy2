'''
stochss_compute.server.run
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

import random
from datetime import datetime
from secrets import token_hex

from tornado.web import RequestHandler
from tornado.ioloop import IOLoop
from distributed import Client, Future
from gillespy2.core import Results
from stochss_compute.core.messages.status import SimStatus
from stochss_compute.core.messages.simulation_run import SimulationRunRequest, SimulationRunResponse
from stochss_compute.server.cache import Cache


class RunHandler(RequestHandler):
    '''
    Endpoint for running Gillespy2 simulations.
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
        sim_request = SimulationRunRequest.parse(self.request.body)
        sim_hash = sim_request.hash()
        log_string = f'{datetime.now()} | <{self.request.remote_ip}> | Simulation Run Request | <{sim_hash}> | '
        cache = Cache(self.cache_dir, sim_hash)
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
                print(log_string +
                    f'Partial cache. Running {trajectories_needed} new trajectories.')
                client = Client(self.scheduler_address)
                future = self._submit(sim_request, sim_hash, client)
                self._return_running(sim_hash, future.key)
                IOLoop.current().run_in_executor(None, self._cache, sim_hash, future, client)
            else:
                print(log_string + 'Returning cached results.')
                results = cache.get()
                ret_traj = random.sample(results, n_traj)
                new_results = Results(ret_traj)
                new_results_json = new_results.to_json()
                sim_response = SimulationRunResponse(SimStatus.READY, results_id = sim_hash, results = new_results_json)
                self.write(sim_response.encode())
                self.finish()
        if empty:
            print(log_string + 'Results not cached. Running simulation.')
            client = Client(self.scheduler_address)
            future = self._submit(sim_request, sim_hash, client)
            self._return_running(sim_hash, future.key)
            IOLoop.current().run_in_executor(None, self._cache, sim_hash, future, client)

    def _cache(self, sim_hash, future: Future, client: Client):
        results = future.result()
        client.close()
        cache = Cache(self.cache_dir, sim_hash)
        cache.save(results)

    def _submit(self, sim_request, sim_hash, client: Client):
        model = sim_request.model
        kwargs = sim_request.kwargs
        n_traj = kwargs.get('number_of_trajectories', 1)
        if "solver" in kwargs:
            from pydoc import locate
            kwargs["solver"] = locate(kwargs["solver"])

        # keep client open for now! close?
        key = f'{sim_hash}:{n_traj}:{token_hex(8)}'
        future = client.submit(model.run, **kwargs, key=key)
        return future

    def _return_running(self, results_id, task_id):
        sim_response = SimulationRunResponse(SimStatus.RUNNING, results_id=results_id, task_id=task_id)
        self.write(sim_response.encode())
        self.finish()
