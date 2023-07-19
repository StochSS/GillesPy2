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

import os
from tornado.web import RequestHandler
from gillespy2.remote.core.messages.source_ip import SourceIpRequest, SourceIpResponse

class SourceIpHandler(RequestHandler):
    '''
    Responds with the IP address associated with the request.
    Used only by cloud api.
    '''

    def post(self):
        '''
        Process POST request.
        
        :returns: request.remote_ip
        :rtype: str
        '''
        source_ip = self.request.remote_ip
        print(f'[SourceIp Request] | Source: <{source_ip}>')
        source_ip_request = SourceIpRequest.parse(self.request.body)
        # could possibly also check just to see if request is valid?
        if source_ip_request.cloud_key == os.environ.get('CLOUD_LOCK'):
            source_ip_response = SourceIpResponse(source_ip=source_ip)
            self.write(source_ip_response.encode())
        else:
            self.set_status(403, 'Access denied.')
        self.finish()
