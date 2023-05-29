'''
gillespy2.remote.server.api
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
import asyncio
import subprocess
from logging import DEBUG, INFO
from tornado.web import Application
from gillespy2.remote.server.simulation_run import SimulationRunHandler
from gillespy2.remote.server.simulation_run_cache import SimulationRunCacheHandler
from gillespy2.remote.server.sourceip import SourceIpHandler
from gillespy2.remote.server.status import StatusHandler
from gillespy2.remote.server.results import ResultsHandler
from gillespy2.remote.core.log_config import init_logging, set_global_log_level
log = init_logging(__name__)

def _make_app(dask_host, dask_scheduler_port, cache):
    scheduler_address = f'{dask_host}:{dask_scheduler_port}'
    return Application([
        (r"/api/v2/simulation/gillespy2/run", SimulationRunHandler,
            {'scheduler_address': scheduler_address, 'cache_dir': cache}),
        (r"/api/v2/simulation/gillespy2/(?P<results_id>[0-9a-zA-Z][0-9a-zA-Z]*?)/results",
            ResultsHandler, {'cache_dir': cache}),
        (r"/api/v2/simulation/gillespy2/run/cache", SimulationRunCacheHandler,
            {'scheduler_address': scheduler_address, 'cache_dir': cache}),
        (r"/api/v2/simulation/gillespy2/(?P<results_id>[0-9a-zA-Z][0-9a-zA-Z]*?)/(?P<n_traj>[1-9][0-9]*?)/(?P<task_id>[0-9a-zA-Z][0-9a-zA-Z]*?)/status",
            StatusHandler, {'scheduler_address': scheduler_address, 'cache_dir': cache}),
        (r"/api/v2/simulation/gillespy2/(?P<results_id>.*?)/(?P<n_traj>[1-9]\d*?)/results",
            ResultsHandler, {'cache_dir': cache}),
        # (r"/api/v2/cache/gillespy2/(?P<results_id>.*?)/(?P<n_traj>[1-9]\d*?)/is_cached",
        #     IsCachedHandler, {'cache_dir': cache}),
        (r"/api/v2/cloud/sourceip", SourceIpHandler),
    ])

async def start_api(
        port = 29681,
        cache_trajectories = True,
        cache_path = 'cache/',
        dask_host = 'localhost',
        dask_scheduler_port = 8786,
        rm = False,
        logging_level = DEBUG,
        ):
    """
    Start the REST API with the following arguments.

    :param port: The port to listen on.
    :type port: int

    :param cache_trajectories: If True, default behavior is to cache trajectories. If False, trajectory cacheing is turned off by default. Can be overridden on client side.
    :type cache_trajectories: bool

    :param cache_path: The cache directory path.
    :type cache_path: str

    :param dask_host: The address of the dask cluster.
    :type dask_host: str

    :param dask_scheduler_port: The port of the dask cluster.
    :type dask_scheduler_port: int

    :param rm: Delete the cache when exiting this program.
    :type rm: bool

    :param logging_level: Set log level for gillespy2.remote.
    :type debug: logging._Level
    """

    set_global_log_level(logging_level)
    # TODO clean up lock files here

    cache_path = os.path.abspath(cache_path)
    app = _make_app(dask_host, dask_scheduler_port, cache_path)
    app.listen(port)
    msg='''
=========================================================================
  StochSS-Compute listening on port: %(port)d                                 
  Cache directory: %(cache_path)s                                        
  Connecting to Dask scheduler at: %(dask_host)s:%(dask_scheduler_port)d 
=========================================================================
'''
    log.info(msg, locals())

    try:
        await asyncio.Event().wait()
    except asyncio.exceptions.CancelledError as error:
        log.error(error)
    finally:
        if rm and os.path.exists(cache_path):
            log.info('Removing cache...')
            subprocess.Popen(['rm', '-r', cache_path])
            log.info('Cache Removed OK')
            