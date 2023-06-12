'''
gillespy2.remote.core.remote_results
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

from time import sleep
from gillespy2 import Results
from gillespy2.remote.client.endpoint import Endpoint
from gillespy2.remote.core.errors import RemoteSimulationError
from gillespy2.remote.core.messages.results import ResultsResponse
from gillespy2.remote.core.messages.status import StatusRequest, StatusResponse, SimStatus

from gillespy2.remote.core.log_config import init_logging
log = init_logging(__name__)

class RemoteResults(Results):
    '''
    Wrapper for a gillespy2.Results object that exists on a remote server and which is then downloaded locally.
    A Results object is: A List of Trajectory objects created by a gillespy2 solver, extends the UserList object.

    These four fields must be initialized manually: id, server, n_traj, task_id.

    :param data: A list of trajectory objects.
    :type data: UserList

    :param id: ID of the cached Results object.
    :type id: str

    :param server: The remote instance of StochSS-Compute where the Results are cached.
    :type server: gillespy2.remote.ComputeServer

    :param task_id: Handle for the running simulation.
    :type task_id: str
    '''
    # These four fields are initialized after object creation.
    id = None
    server = None
    n_traj = None
    task_id = None
    namespace = None

    # pylint:disable=super-init-not-called
    def __init__(self, data = None):
        self._data = data
    # pylint:enable=super-init-not-called

    @property
    def data(self):
        """
        The trajectory data.

        :returns: self._data
        :rtype: UserList
        """
        if None in (self.id, self.server, self.n_traj):
            raise Exception('RemoteResults must have a self.id, self.server and self.n_traj.')

        if self._data is None:
            self._resolve()
        return self._data

    @property
    def sim_status(self):
        '''
        Fetch the simulation status.

        :returns: Simulation status enum as a string.
        :rtype: str
        '''
        if self._data is not None:
            return self._status().status.name

    def get_gillespy2_results(self):
        """
        Get the GillesPy2 results object from the remote results.

        :returns: The generated GillesPy2 results object.
        :rtype: gillespy.Results
        """
        return Results(self.data)


    @property
    def is_ready(self):
        """
        True if results exist in cache on the server.

        :returns: status == SimStatus.READY
        :rtype: bool
        """
        if self._data is not None:
            return self._status().status == SimStatus.READY
        return True

    def _status(self):
        '''
        It is undefined/illegal behavior to call this function if self._data is not None.        
        '''
        if self._data is not None:
            raise Exception('TODO Name this exception class. Cant call status on a finished simulation.')
        # Request the status of a submitted simulation.
        status_request = StatusRequest(self.id, self.namespace)
        response_raw = self.server.post(Endpoint.SIMULATION_GILLESPY2,
                                       f"/status")
        if not response_raw.ok:
            raise RemoteSimulationError(response_raw.reason)

        status_response = StatusResponse.parse(response_raw.text)
        return status_response

    def _resolve(self):
        '''
        It is undefined behavior to call this function if self._data is not None.        
        '''
        status_response = self._status()
        status = status_response.status

        if status == SimStatus.RUNNING:
            log.info('Simulation is running. Downloading results when complete......')
            while True:
                sleep(5)
                status_response = self._status()
                status = status_response.status
                if status != SimStatus.RUNNING:
                    break

        if status in (SimStatus.DOES_NOT_EXIST, SimStatus.ERROR):
            raise RemoteSimulationError(status_response.message)

        if status == SimStatus.READY:
            log.info('Results ready. Fetching.......')
            if self.id == self.task_id:
                response_raw = self.server.get(Endpoint.SIMULATION_GILLESPY2, f"/{self.id}/results")
            else:
                response_raw = self.server.get(Endpoint.SIMULATION_GILLESPY2, f"/{self.id}/{self.n_traj}/results")
            if not response_raw.ok:
                raise RemoteSimulationError(response_raw.reason)

            response = ResultsResponse.parse(response_raw.text)

            self._data = response.results.data # Fully initialized
