'''
gillespy2.remote.server.results
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

from tornado.web import RequestHandler
from gillespy2.remote.core.errors import RemoteSimulationError
from gillespy2.remote.core.messages.results import ResultsResponse
from gillespy2.remote.server.cache import Cache

from gillespy2.remote.core.log_config import init_logging
log = init_logging(__name__)

class ResultsHandler(RequestHandler):
    '''
    Endpoint for Results objects.
    '''
    def __init__(self, application, request, **kwargs):
        self.cache_dir = None
        super().__init__(application, request, **kwargs)

    def data_received(self, chunk: bytes):
        raise NotImplementedError()

    def initialize(self, cache_dir):
        '''
        Set the cache directory.

        :param cache_dir: Path to the cache.
        :type cache_dir: str
        '''
        while cache_dir.endswith('/'):
            cache_dir = cache_dir[:-1]
        self.cache_dir = cache_dir + '/run/'

    async def get(self, results_id = None):
        '''
        Process GET request.

        :param results_id: Hash of the simulation.
        :type results_id: str
        
        :param n_traj: Number of trajectories in the request.
        :type n_traj: str
        '''
        if results_id in ('', '/'):
            self.set_status(404, reason=f'Malformed request: {self.request.uri}')
            self.finish()
            raise RemoteSimulationError(f'Malformed request | <{self.request.remote_ip}>')
        msg = f' <{self.request.remote_ip}> | Results Request | <{results_id}>'
        log.info(msg)
        cache = Cache(self.cache_dir, results_id)
        if cache.is_ready():
            results = cache.read()
            results_response = ResultsResponse(results)
            self.write(results_response.encode())
        else:
            # This should not happen!
            self.set_status(404, f'Results "{results_id}" not found.')
        self.finish()
