'''
gillespy2.remote.launch
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

import sys
import asyncio
from argparse import ArgumentParser, Namespace
from distributed import LocalCluster
from gillespy2.remote.server.api import start_api

from logging import INFO, getLevelName
from gillespy2.remote.core.log_config import init_logging
log = init_logging(__name__)

def _add_shared_args(parser):
    '''
    :type parser: ArgumentParser 
    
    :rtype: ArgumentParser
    '''

    server = parser.add_argument_group('Server')
    server.add_argument("-p", "--port", default=29681, type=int, required=False,
        help="The port to use for the server. Defaults to 29681.")
    server.add_argument("-l", "--logging-level", default=INFO, required=False,
        help='Set the logging level threshold. Str or int. Defaults to INFO (20).')

    cache = parser.add_argument_group('Cache')
    cache.add_argument('-c', '--cache_path', default='cache/', required=False,
        help='Path to use for the cache.')
    cache.add_argument('--rm', default=False, action='store_true', required=False,
        help='Whether to delete the cache upon exit. Default False.')
    
    return parser


def launch_server():
    '''
    Start the REST API. Alias to script "gillespy2-remote".

    `gillespy2-remote --help`
    OR
    `python -m gillespy2.remote.launch --help`
    '''
    def _parse_args() -> Namespace:
        desc = '''
GillesPy2 Remote allows you to run simulations remotely on your own Dask cluster.
To launch both simultaneously, use `gillespy2-remote-cluster` instead.
Trajectories are automatically cached to support multiple users running the same model.
'''
        parser = ArgumentParser(description=desc, add_help=True, conflict_handler='resolve')

        parser = _add_shared_args(parser)

        dask = parser.add_argument_group('Dask')
        dask.add_argument("-H", "--dask-host", default='localhost', required=False,
                            help="The host to use for the dask scheduler. Defaults to localhost.")
        dask.add_argument("-P", "--dask-scheduler-port", default=8786, type=int, required=False,
                            help="The port to use for the dask scheduler. Defaults to 8786.")
        return parser.parse_args()

    args = _parse_args()
    args.logging_level = getLevelName(args.logging_level)

    asyncio.run(start_api(**args.__dict__))


def launch_with_cluster():
    '''
    Start up a Dask cluster along with gillespy2.remote REST API. Alias to script "gillespy2-remote-cluster".

    `gillespy2-remote-cluster --help`
    OR
    `python -m gillespy2.remote.launch cluster --help`
    '''

    def _parse_args() -> Namespace:
        desc = '''
        Startup script for a GillesPy2 Remote and pre-configured Dask Distributed cluster.
        Command-line options allow you to override automatic cluster configuration.
        Your trajectories are automatically cached to support multiple users running the same model.'''
        parser = ArgumentParser(description=desc, add_help=True, conflict_handler='resolve')

        parser = _add_shared_args(parser)

        dask = parser.add_argument_group('Dask')
        dask.add_argument("-H", "--dask-host", default=None, required=False,
            help="The host to use for the dask scheduler. Defaults to localhost.")
        dask.add_argument("-P", "--dask-scheduler-port", default=0, type=int, required=False,
            help="The port to use for the dask scheduler. 0 for a random port. \
            Defaults to a random port.")
        dask.add_argument('-W', '--dask-n-workers', default=None, type=int, required=False,
            help='Configure the number of workers. Defaults to one per core.')
        dask.add_argument('-T', '--dask-threads-per-worker', default=None, required=False, type=int,
            help='Configure the threads per worker. \
            Default will let Dask decide based on your CPU.')
        dask.add_argument('--dask-processes', default=None, required=False, type=bool,
            help='Whether to use processes (True) or threads (False). \
            Defaults to True, unless worker_class=Worker, in which case it defaults to False.')
        dask.add_argument('-D', '--dask-dashboard-address', default=':8787', required=False,
            help='Address on which to listen for the Bokeh diagnostics server \
            like ‘localhost:8787’ or ‘0.0.0.0:8787’. Defaults to ‘:8787’. \
            Set to None to disable the dashboard. Use ‘:0’ for a random port.')
        dask.add_argument('-N', '--dask-name', default=None, required=False,
            help='A name to use when printing out the cluster, defaults to the type name.')        
        
        return  parser.parse_args()


    args = _parse_args()
    args.logging_level = getLevelName(args.logging_level)

    dask_args = {}
    for (arg, value) in vars(args).items():
        if arg.startswith('dask_'):
            dask_args[arg[5:]] = value
    log.info('Launching Dask Cluster.')
    cluster = LocalCluster(**dask_args)
    tokens = cluster.scheduler_address.split(':')
    dask_host = tokens[1][2:]
    dask_port = int(tokens[2])
    msg = f'Scheduler Address: <{cluster.scheduler_address}>'
    log.info(msg)
    for i, worker in cluster.workers.items():
        msg = f'Worker {i}: {worker}'
        log.info(msg)

    msg = f'Dashboard Link: <{cluster.dashboard_link}>\n'
    log.info(msg)

    try:
        asyncio.run(start_api(port=args.port, cache_path=args.cache_path,
            dask_host=dask_host, dask_scheduler_port=dask_port, rm=args.rm, logging_level=args.logging_level))
    except asyncio.exceptions.CancelledError:
        pass
    finally:
        log.info('Shutting down cluster.')
        asyncio.run(cluster.close())
        log.info('Cluster terminated.')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] == 'cluster':
            del sys.argv[1]
            launch_with_cluster()
    launch_server()
