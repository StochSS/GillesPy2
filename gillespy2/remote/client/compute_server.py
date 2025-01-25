'''
ComputeServer(Server)
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

from gillespy2.remote.client.server import Server

class ComputeServer(Server):
    '''
    Simple object representing a remote instance of StochSS-Compute.

    :param host: Address of the remote server.
    :type host: str

    :param port: Port on which to connect. Defaults to 29681.
    :type port: int
    '''
    # pylint: disable=super-init-not-called
    def __init__(self, host: str, port: int = 29681):
        host = host.replace('http://','')
        host = host.split(':')[0]
        self._address = f"http://{host}:{port}"
    # pylint: enable=super-init-not-called

    @property
    def address(self) -> str:
        """
        The server's IP address and port.

        :returns: "http://{host}:{port}"
        """
        return self._address
