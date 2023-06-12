'''
gillespy2.remote.server.cache
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
from json.decoder import JSONDecodeError
from datetime import datetime
import random
from filelock import SoftFileLock
from gillespy2 import Results

class Cache:
    '''
    Cache

    :param cache_dir: The root cache directory.
    :type cache_dir: str

    :param results_id: Simulation hash.
    :type results_id: str
    '''
    def __init__(self, cache_dir, results_id):
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
        self.results_path = os.path.join(cache_dir, f'{results_id}.results')

    def create(self) -> None:
        '''
        Create the results file if it does not exist.
        '''
        try:
            with open(self.results_path, 'x', encoding='utf-8') as file:
                file.close()
        except FileExistsError:
            pass

    def exists(self) -> bool:
        '''
        Check if the results file exists.

        :returns: os.path.exists(self.results_path)
        :rtype: bool
        '''
        return os.path.exists(self.results_path)

    def is_empty(self) -> bool:
        '''
        Check if the results are empty.

        :returns: filesize == 0 or self.exists()
        :rtype: bool
        '''
        lock = SoftFileLock(f'{self.results_path}.lock')
        with lock:
            if self.exists():
                filesize = os.path.getsize(self.results_path)
                return filesize == 0
            return True

    def is_ready(self, n_traj_wanted=0) -> bool:
        '''
        Check if the results are ready to be retrieved from the cache.

        :param n_traj_wanted: The number of requested trajectories.
        :type: int

        :returns: n_traj_wanted <= len(<Results in cache>)
        :rtype: bool
        '''
        results = self.get()
        if results is None or n_traj_wanted > len(results):
            return False
        return True

    def n_traj_needed(self, n_traj_wanted) -> int:
        '''
        Calculate the difference between the number of trajectories the user has requested
         and the number of trajectories currently in the cache.

        :param n_traj_wanted: The number of requested trajectories.
        :type: int

        :returns: A number greater than or equal to zero.
        :rtype: int
        '''
        if self.is_empty():
            return n_traj_wanted
        results = self.get()
        if results is None:
            return n_traj_wanted
        diff = n_traj_wanted - len(results)
        if diff > 0:
            return diff
        return 0

    def n_traj_in_cache(self) -> int:
        '''
        Check the number of trajectories in the cache.

        :returns: `len()` of the gillespy2.Results
        :rtype: int
        '''
        if self.is_empty():
            return 0
        results = self.get()
        if results is not None:
            return len(results)
        return 0

    def get(self):
        '''
        Retrieve a gillespy2.Results object from the cache or None if error.

        :returns: Results.from_json(results_json)
        :rtype: gillespy2.Results or None
        '''
        try:
            results_json = self.read()
            return Results.from_json(results_json)
        except JSONDecodeError:
            return None

    def get_sample(self, n_traj) -> Results or None:
        '''
        Retrieve a gillespy2.Results by sampling from the cache or None if error.

        :returns: Results.from_json(results_json)
        :rtype: gillespy2.Results or None
        '''
        results = self.get()
        if results is None or len(results) <= n_traj:
            return results
        ret_traj = random.sample(results, n_traj)
        return Results(ret_traj)

    def read(self) -> str:
        '''
        Retrieve a gillespy2.Results object as a JSON-formatted string.

        :returns: The output of reading the file.
        :rtype: str
        '''
        lock = SoftFileLock(f'{self.results_path}.lock')
        with lock:
            with open(self.results_path,'r', encoding='utf-8') as file:
                return file.read()

    def save(self, results: Results) -> None:
        '''
        Save a newly processed gillespy2.Results object to the cache.

        :param results: The new Results.
        :type: gillespy2.Results
        '''
        msg = f'{datetime.now()} | Cache | <{self.results_path}> | '
        lock = SoftFileLock(f'{self.results_path}.lock')
        with lock:
            with open(self.results_path, 'r+', encoding='utf-8') as file:
                try:
                    old_results = Results.from_json(file.read())
                    combined_results = results + old_results
                    print(msg+'Add')
                    file.seek(0)
                    file.write(combined_results.to_json())
                except JSONDecodeError:
                    print(msg+'New')
                    file.seek(0)
                    file.write(results.to_json())
