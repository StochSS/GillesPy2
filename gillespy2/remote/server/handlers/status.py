'''
gillespy2.remote.server.status
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
from distributed import Client
from tornado.web import RequestHandler
from gillespy2.remote.core.errors import RemoteSimulationError
from gillespy2.remote.core.messages.status import SimStatus, StatusRequest, StatusResponse
from gillespy2.remote.server.cache import Cache

from gillespy2.remote.core.utils.log_config import init_logging
log = init_logging(__name__)

class StatusHandler(RequestHandler):
    '''
    Endpoint for requesting the status of a simulation.
    '''
    def __init__(self, application, request, **kwargs):
        self.scheduler_address = None
        self.cache_dir = None
        super().__init__(application, request, **kwargs)

    def data_received(self, chunk: bytes):
        raise NotImplementedError()

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

    async def get(self):
        '''
        Process Status GET request.
        '''
        
        # status_request = StatusRequest.parse(self.request)
        log.debug(self.request.query_arguments)

        results_id = self.get_query_argument('results_id')
        n_traj = self.get_query_argument('n_traj', None)
        if n_traj is not None:
            n_traj = int(n_traj)
        task_id = self.get_query_argument('task_id', None)
        namespace = self.get_query_argument('namespace', None)

        if namespace is not None:
            self.cache_dir = os.path.join(self.cache_dir, namespace)
        # if results_id == task_id: # True iff call made using (ignore_cache=True)
        #     self.cache_dir = os.path.join(self.cache_dir, 'run/')
        
        cache = Cache(self.cache_dir, results_id)
        
        msg_0 = f'<{self.request.remote_ip}> | Results ID: <{results_id}> | Trajectories: {n_traj} | Task ID: {task_id}'
        log.info(msg_0)

        msg_1 = f'<{results_id}> | <{task_id}> | Status:'
        dne_msg = f'{msg_1} {SimStatus.DOES_NOT_EXIST.name}'
        ready_msg = f'{msg_1} {SimStatus.READY.name}'
        
        if cache.exists():
            log.debug('cache.exists(): True')
        
            if cache.is_empty():
        
                if task_id is not None:
        
                    state, err = await self._check_with_scheduler(task_id)

                    msg_2 = f'{msg_1} {SimStatus.RUNNING.name} | Task: {state} | Error: {err}'
                    log.info(msg_2)
        
                    if state == 'erred':
                        self._respond_error(err)
                    else:
                        self._respond_running(f'Scheduler Task State: {state}')
                
                else:

                    log.info(dne_msg)
                    self._respond_dne()
            
            else:
            
                if cache.is_ready(n_traj):
                    log.info(ready_msg)
                    self._respond_ready()

                else:

                    if task_id is not None:

                        state, err = await self._check_with_scheduler(task_id)

                        msg_2 = f'{msg_1} {SimStatus.RUNNING.name} | Task: {state} | Error: {err}'
                        log.info(msg_2)
                        
                        if state == 'erred':
                            self._respond_error(err)
                        else:
                            self._respond_running(f'Scheduler task state: {state}')
                    
                    else:

                        log.info(dne_msg)
                        self._respond_dne()

        else:
            log.debug('cache.exists(): False')
            log.info(dne_msg)
            self._respond_dne()

    def _respond_ready(self):
        status_response = StatusResponse(SimStatus.READY)
        self.write(status_response.encode())
        self.finish()

    def _respond_error(self, error_message):
        status_response = StatusResponse(SimStatus.ERROR, error_message)
        self.write(status_response.encode())
        self.finish()

    def _respond_dne(self):
        status_response = StatusResponse(SimStatus.DOES_NOT_EXIST, 'There is no record of that simulation.')
        self.write(status_response.encode())
        self.finish()

    def _respond_running(self, message):
        status_response = StatusResponse(SimStatus.RUNNING, message)
        self.write(status_response.encode())
        self.finish()

    async def _check_with_scheduler(self, task_id):
        '''
        Ask the scheduler for information about a task.
        '''
        client = Client(self.scheduler_address)

        # define function here so that it is pickle-able
        def scheduler_task_state(task_id, dask_scheduler=None):
            task = dask_scheduler.tasks.get(task_id)
            if task is None:
                return (None, None)
            if task.exception_text == "":
                return (task.state, None)
            return (task.state, task.exception_text)
        
        # Do not await. Reasons. It returns sync.
        _ = client.run_on_scheduler(scheduler_task_state, task_id)
        client.close()
        return _
