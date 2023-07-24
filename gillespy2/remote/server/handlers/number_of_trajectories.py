'''
gillespy2.remote.server.number_of_trajectories
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

from tornado.web import RequestHandler

from gillespy2.remote.core.utils.log_config import init_logging
from gillespy2.remote.server.cache import Cache
log = init_logging(__name__)

class NumberOfTrajectoriesHandler(RequestHandler):
    '''
    Responds with the number of trajectories in the cache.
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
        self.cache_dir = cache_dir

    def get(self):
        '''
        Process GET request.
        
        :returns: Number of trajectories in the cache.
        :rtype: str
        '''
        msg = f'Request from {self.request.remote_ip}'
        log.info(msg)
        results_id = self.get_query_argument('results_id')
        namespace = self.get_query_argument('namespace', None)
        msg = f'results_id: <{results_id}> | namespace: <{namespace}>'
        log.info(msg)
        cache = Cache(self.cache_dir, results_id, namespace=namespace)
        n_traj_in_cache = cache.n_traj_in_cache()
        self.write(str(n_traj_in_cache))
        self.finish()
