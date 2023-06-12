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
from tornado.web import RequestHandler
from tornado.ioloop import IOLoop
from distributed import Client
from gillespy2.remote.core.messages.status import SimStatus
from gillespy2.remote.core.exceptions import PRNGCollision
from gillespy2.remote.core.messages.simulation_run import SimulationRunRequest, SimulationRunResponse
from gillespy2.remote.server.cache import Cache

from gillespy2.remote.core.log_config import init_logging
log = init_logging(__name__)

class SimulationRunHandler(RequestHandler):
    '''
    Endpoint for running GillesPy2 simulations.
    '''

    def __init__(self, application, request, **kwargs):
        self.scheduler_address = None
        self.cache_dir = None
        self.key = None
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
        self.cache_dir = cache_dir

    async def post(self):
        '''
        Process simulation run  POST request.
        '''
        sim_request = SimulationRunRequest.parse(self.request.body)
        if sim_request.namespace is not None:
            self.cache_dir = os.path.join(self.cache_dir, sim_request.namespace)
        cache = Cache(self.cache_dir, sim_request.id)
        if cache.exists():
            log.debug("This should not be happening.")
            self.finish()
            raise PRNGCollision
        try:
            cache.create()
            client = Client(self.scheduler_address)
            msg = f'<{self.request.remote_ip}> | <{sim_request.id}> | Running simulation.'
            log.info(msg)
            future = self._submit(sim_request, client)
            IOLoop.current().run_in_executor(None, self._cache, cache, future, client)
            self._return_running()
        except Exception as err:
            self._return_error(str(err))
            

    def _cache(self, cache, future, client):
        '''
        Await results, close client, save to disk.

        :param cache: Handle to the cache.
        :type cache: Cache

        :param future: Handle to the running simulation, to be awaited upon.
        :type future: distributed.Future

        :param client: Client to the Dask scheduler. Closing here for good measure, not sure if strictly necessary.
        :type client: distributed.Client
        '''
        results = future.result()
        client.close()
        cache.save(results)

    def _submit(self, sim_request, client):
        '''
        Submit request to dask scheduler.
        Uses pydoc.locate to convert str to solver class name object.

        :param sim_request: The user's request for a  simulation.
        :type sim_request: SimulationRunRequest

        :param client: Client to the Dask scheduler.
        :type client: distributed.Client

        :returns: Handle to the running simulation and the results on the worker.
        :rtype: distributed.Future
        '''
        model = sim_request.model
        kwargs = sim_request.kwargs
        key = sim_request.id
        kwargs['seed'] = int(sim_request.id, 16)
        if "solver" in kwargs:
            # pylint:disable=import-outside-toplevel
            from pydoc import locate
            # pylint:enable=import-outside-toplevel
            kwargs["solver"] = locate(kwargs["solver"])

        future = client.submit(model.run, **kwargs, key=key)
        return future

    def _return_running(self):
        '''
        Let the user know we submitted the simulation to the scheduler.
        '''
        sim_response = SimulationRunResponse(SimStatus.RUNNING)
        self.write(sim_response.encode())
        self.finish()

    def _return_error(self, error_message):
        '''
        Let the user know we submitted the simulation to the scheduler.
        '''
        sim_response = SimulationRunResponse(SimStatus.ERROR, error_message)
        self.write(sim_response.encode())
        self.finish()
