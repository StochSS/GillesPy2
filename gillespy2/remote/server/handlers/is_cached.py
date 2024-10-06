'''
gillespy2.remote.server.is_cached
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
from gillespy2.remote.core.errors import RemoteSimulationError
from gillespy2.remote.core.messages.status import SimStatus, StatusResponse
from gillespy2.remote.server.cache import Cache

class IsCachedHandler(RequestHandler):
    '''
    Endpoint that will determine if a particular simulation exists in the cache.
    '''
    def initialize(self, cache_dir):
        '''
        Set cache path.

        :param cache_dir: Path to the cache.
        :type cache_dir: str
        '''
        self.cache_dir = cache_dir

    async def get(self, results_id = None, n_traj = None):
        '''
        Process GET request.

        :param results_id: Hash of the simulation.
        :type results_id: str

        :param n_traj: Number of trajectories to check for.
        :type n_traj: str
        '''
        if None in (results_id, n_traj):
            raise RemoteSimulationError('Malformed request')
        n_traj = int(n_traj)
        cache = Cache(self.cache_dir, results_id)
        print(f'\
{datetime.now()} |\
 Source: <{self.request.remote_ip}> |\
 Cache Check |\
 <{results_id}> |\
 Trajectories: {n_traj} ')
        msg = f'{datetime.now()} | <{results_id}> | Trajectories: {n_traj} | Status: '
        exists = cache.exists()
        if exists:
            empty = cache.is_empty()
            if empty:
                print(msg+SimStatus.DOES_NOT_EXIST.name)
                self._respond_dne('That simulation is not currently cached.')
            else:
                ready = cache.is_ready(n_traj)
                if ready:
                    print(msg+SimStatus.READY.name)
                    self._respond_ready()
                else:
                    print(msg+SimStatus.DOES_NOT_EXIST.name)
                    self._respond_dne(f'Not enough trajectories in cache. \
                                      Requested: {n_traj}, Available: {cache.n_traj_in_cache()}')
        else:
            print(msg+SimStatus.DOES_NOT_EXIST.name)
            self._respond_dne('There is no record of that simulation')

    def _respond_ready(self):
        status_response = StatusResponse(SimStatus.READY)
        self.write(status_response.encode())
        self.finish()

    def _respond_dne(self, msg):
        status_response = StatusResponse(SimStatus.DOES_NOT_EXIST, msg)
        self.write(status_response.encode())
        self.finish()
