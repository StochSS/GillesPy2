'''
gillespy2.remote.server.status
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

from distributed import Client
from tornado.web import RequestHandler
from gillespy2.remote.core.errors import RemoteSimulationError
from gillespy2.remote.core.messages.status import SimStatus, StatusResponse
from gillespy2.remote.server.cache import Cache

from gillespy2.remote.core.log_config import init_logging
log = init_logging(__name__)

class StatusHandler(RequestHandler):
    '''
    Endpoint for requesting the status of a simulation.
    '''
    def __init__(self, application, request, **kwargs):
        self.scheduler_address = None
        self.cache_dir = None
        self.task_id = None
        self.results_id = None
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

    async def get(self, results_id, n_traj, task_id):
        '''
        Process GET request.

        :param results_id: Hash of the simulation. Required.
        :type results_id: str

        :param n_traj: Number of trajectories in the request. Default 1.
        :type n_traj: str

        :param task_id: ID of the running simulation. Required.
        :type task_id: str
        '''
        if '' in (results_id, n_traj):
            self.set_status(404, reason=f'Malformed request: {self.request.uri}')
            self.finish()
            raise RemoteSimulationError(f'Malformed request: {self.request.uri}')
        self.results_id = results_id
        self.task_id = task_id
        n_traj = int(n_traj)
        unique = results_id == task_id
        log.debug('unique: %(unique)s', locals())
        cache = Cache(self.cache_dir, results_id, unique=unique)
        log_string = f'<{self.request.remote_ip}> | Results ID: <{results_id}> | Trajectories: {n_traj} | Task ID: {task_id}'
        log.info(log_string)

        msg = f'<{results_id}> | <{task_id}> | Status: '

        exists = cache.exists()
        log.debug('exists: %(exists)s', locals())
        if exists:
            empty = cache.is_empty()
            if empty:
                if self.task_id not in ('', None):
                    state, err = await self._check_with_scheduler()
                    log.info(msg + SimStatus.RUNNING.name + f' | Task: {state} | Error: {err}')
                    if state == 'erred':
                        self._respond_error(err)
                    else:
                        self._respond_running(f'Scheduler task state: {state}')
                else:
                    log.info(msg+SimStatus.DOES_NOT_EXIST.name)
                    self._respond_dne()
            else:
                ready = cache.is_ready(n_traj)
                if ready:
                    log.info(msg+SimStatus.READY.name)
                    self._respond_ready()
                else:
                    if self.task_id not in ('', None):
                        state, err = await self._check_with_scheduler()
                        log.info(msg+SimStatus.RUNNING.name+f' | Task: {state} | error: {err}')
                        if state == 'erred':
                            self._respond_error(err)
                        else:
                            self._respond_running(f'Scheduler task state: {state}')
                    else:
                        log.info(msg+SimStatus.DOES_NOT_EXIST.name)
                        self._respond_dne()
        else:
            log.info(msg+SimStatus.DOES_NOT_EXIST.name)
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

    async def _check_with_scheduler(self):
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
        ret = client.run_on_scheduler(scheduler_task_state, self.task_id)
        client.close()
        return ret
