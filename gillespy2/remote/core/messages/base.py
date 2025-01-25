'''
gillespy2.remote.core.messages
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

from abc import ABC, abstractmethod

class Request(ABC):
    '''
    Base class.
    '''
    @abstractmethod
    def encode(self):
        '''
        Encode self for http.
        '''
    @staticmethod
    @abstractmethod
    def parse(raw_request):
        '''
        Parse http for python.
        '''

class Response(ABC):
    '''
    Base class.
    '''
    @abstractmethod
    def encode(self):
        '''
        Encode self for http.
        '''
    @staticmethod
    @abstractmethod
    def parse(raw_response):
        '''
        Parse http for python.
        '''
