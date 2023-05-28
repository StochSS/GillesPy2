'''
stochss_compute.cloud.exceptions
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

class ResourceException(Exception):
    '''
    Misconfigured or out-of-date resources detected in the cluster setup.
    '''
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
        print('Missing or misconfigured resources.')

class EC2ImportException(Exception):
    '''
    Some extra dependencies are required for EC2.
    '''
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
        print('StochSS-Compute on EC2 requires boto3 and paramiko to be installed.\nTry: pip install boto3 paramiko')

class EC2Exception(Exception):
    '''
    General exception class.
    '''
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
