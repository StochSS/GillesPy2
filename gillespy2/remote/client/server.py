'''
Server(ABC)
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

from time import sleep
from abc import ABC, abstractmethod
import requests
from gillespy2.remote.client.endpoint import Endpoint
from gillespy2.remote.core.messages.base import Request

from gillespy2.remote.core.utils.log_config import init_logging
log = init_logging(__name__)

class Server(ABC):
    '''
    Abstract Server class with hard coded endpoints.

    :raises TypeError: Server cannot be instantiated directly. Must be ComputeServer or Cluster.
    '''

    _endpoints = {
        Endpoint.SIMULATION_GILLESPY2: "/api/v3/simulation/gillespy2",
        Endpoint.CLOUD: "/api/v3/cloud",
        Endpoint.CACHE: "/api/v3/cache",
        Endpoint.DASK: '/api/v3/dask',
    }

    def __init__(self) -> None:
        raise TypeError('Server cannot be instantiated directly. Must be ComputeServer or Cluster.')

    @property
    @abstractmethod
    def address(self):
        '''
        NotImplemented
        '''
        return NotImplemented

    def get(self, endpoint: Endpoint, sub: str, request: Request = None):
        '''
        Send a GET request to endpoint.

        :param endpoint: The API endpoint.
        :type endpoint: Endpoint

        :param sub: Final part of url string.
        :type sub: str

        :param request: An object that inherits from Request.
        :type request: Request

        :returns: The HTTP response.
        :rtype: requests.Response
        '''
        if request is not None:
            log.debug(request.encode())
        url = f"{self.address}{self._endpoints[endpoint]}{sub}"
        log.debug(url)
        n_try = 1
        sec = 3
        while n_try <= 3:
            try:
                if request is not None:
                    return requests.get(url, timeout=30, params=request.encode())
                return requests.get( url, timeout=30)

            except ConnectionError:
                print(f"Connection refused by server. Retrying in {sec} seconds....")
                sleep(sec)
                n_try += 1
                sec *= n_try
            except Exception as err:
                print(f"Unknown error: {err}. Retrying in {sec} seconds....")
                sleep(sec)
                n_try += 1
                sec *= n_try

    def post(self, endpoint: Endpoint, sub: str, request: Request = None):
        '''
        Send a POST request to endpoint.

        :param endpoint: The API endpoint.
        :type endpoint: Endpoint

        :param sub: Final part of url string.
        :type sub: str

        :param request: An object that inherits from Request.
        :type request: Request

        :returns: The HTTP response.
        :rtype: requests.Response
        '''
        log.debug(request.encode())

        if self.address is NotImplemented:
            raise NotImplementedError

        url = f"{self.address}{self._endpoints[endpoint]}{sub}"
        log.debug(url)
        n_try = 1
        sec = 3
        while n_try <= 3:
            try:
                if request is None:
                    print(f"[POST] {url}")
                    return requests.post(url, timeout=30)
                print(f"[{type(request).__name__}] {url}")
                return requests.post(url, json=request.encode(), timeout=30)

            except ConnectionError:
                print(f"Connection refused by server. Retrying in {sec} seconds....")
                sleep(sec)
                n_try += 1
                sec *= n_try
            except Exception as err:
                print(f"Unknown error: {err}. Retrying in {sec} seconds....")
                sleep(sec)
                n_try += 1
                sec *= n_try
