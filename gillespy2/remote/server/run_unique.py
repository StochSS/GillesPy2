'''
stochss_compute.server.run_unique
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

from datetime import datetime

from tornado.web import RequestHandler
from tornado.ioloop import IOLoop
from distributed import Client, Future
from gillespy2.remote.core.messages.status import SimStatus
from gillespy2.remote.core.exceptions import PRNGCollision
from gillespy2.remote.core.messages.simulation_run_unique import SimulationRunUniqueRequest, SimulationRunUniqueResponse
from gillespy2.remote.server.cache import Cache


class SimulationRunUniqueHandler(RequestHandler):
    '''
    Endpoint for running Gillespy2 simulations.
    '''

    def __init__(self, application, request, **kwargs):
        self.scheduler_address = None
        self.cache_dir = None
        self.unique_key = None
        super().__init__(application, request, **kwargs)

    def data_received(self, chunk: bytes):
        raise NotImplementedError()
    
    def initialize(self, scheduler_address, cache_dir):
        '''
        Sets the address to the Dask scheduler and the cache directory.
        Creates a new directory for one-off results files identifiable by token.

        :param scheduler_address: Scheduler address.
        :type scheduler_address: str

        :param cache_dir: Path to the cache.
        :type cache_dir: str
        '''
        self.scheduler_address = scheduler_address
        while cache_dir.endswith('/'):
            cache_dir = cache_dir[:-1]
        self.cache_dir = cache_dir + '/unique/'

    async def post(self):
        '''
        Process simulation run unique POST request.
        '''
        sim_request = SimulationRunUniqueRequest.parse(self.request.body)
        self.unique_key = sim_request.unique_key
        cache = Cache(self.cache_dir, self.unique_key)
        if cache.exists():
            self.set_status(404, reason='Try again with a different key, because that one is taken.')
            self.finish()
            raise PRNGCollision('Try again with a different key, because that one is taken.')
        cache.create()
        client = Client(self.scheduler_address)
        future = self._submit(sim_request, client)
        log_string = f'{datetime.now()} | <{self.request.remote_ip}> | Simulation Run Unique Request | <{self.unique_key}> | '
        print(log_string + 'Running simulation.')
        self._return_running()
        IOLoop.current().run_in_executor(None, self._cache, future, client)

    def _cache(self, future, client):
        '''
        Await results, close client, save to disk.

        :param future: Handle to the running simulation, to be awaited upon.
        :type future: distributed.Future

        :param client: Client to the Dask scheduler. Closing here for good measure, not sure if strictly necessary.
        :type client: distributed.Client
        '''
        results = future.result()
        client.close()
        cache = Cache(self.cache_dir, self.unique_key)
        cache.save(results)

    def _submit(self, sim_request, client):
        '''
        Submit request to dask scheduler.
        Uses pydoc.locate to convert str to solver class name object.

        :param sim_request: The user's request for a unique simulation.
        :type sim_request: SimulationRunUniqueRequest

        :param client: Client to the Dask scheduler.
        :type client: distributed.Client

        :returns: Handle to the running simulation and the results on the worker.
        :rtype: distributed.Future
        '''
        model = sim_request.model
        kwargs = sim_request.kwargs
        unique_key = sim_request.unique_key
        if "solver" in kwargs:
            # pylint:disable=import-outside-toplevel
            from pydoc import locate
            # pylint:enable=import-outside-toplevel
            kwargs["solver"] = locate(kwargs["solver"])

        future = client.submit(model.run, **kwargs, key=unique_key)
        return future

    def _return_running(self):
        '''
        Let the user know we submitted the simulation to the scheduler.
        '''
        sim_response = SimulationRunUniqueResponse(SimStatus.RUNNING)
        self.write(sim_response.encode())
        self.finish()
