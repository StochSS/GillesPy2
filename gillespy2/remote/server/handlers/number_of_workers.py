'''
gillespy2.remote.server.sourceip
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
from distributed import Client

from gillespy2.remote.core.utils.log_config import init_logging
log = init_logging(__name__)

class NumberOfWorkersHandler(RequestHandler):
    '''
    Responds with the number of workers attached to the scheduler.
    '''
    scheduler_address = None
    def initialize(self, scheduler_address):
        '''
        Sets the address to the Dask scheduler and the cache directory.

        :param scheduler_address: Scheduler address.
        :type scheduler_address: str

        :param cache_dir: Path to the cache.
        :type cache_dir: str
        '''

        self.scheduler_address = scheduler_address

    def get(self):
        '''
        Process POST request.
        
        :returns: request.remote_ip
        :rtype: str
        '''
        msg = f' <{self.request.remote_ip}>'
        log.info(msg)
        client = Client(self.scheduler_address)
        n_workers =  len(client.scheduler_info().get('workers', []))
        self.write(str(n_workers))
        self.finish()
        client.close()
