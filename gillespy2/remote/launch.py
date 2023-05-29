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

def launch_server():
    '''
    Start the REST API. Alias to script "gillespy2-remote".

    `gillespy2-remote --help`
    OR
    `python -m gillespy2.remote.launch --help`
    '''
    def _parse_args() -> Namespace:
        desc = '''
            StochSS-Compute is a server and cache that anonymizes StochSS simulation data.
        '''
        parser = ArgumentParser(description=desc, add_help=True, conflict_handler='resolve')

        server = parser.add_argument_group('Server')
        server.add_argument("-p", "--port", default=29681, type=int, required=False,
                            help="The port to use for the server. Defaults to 29681.")

        cache = parser.add_argument_group('Cache')
        cache.add_argument('-c', '--cache', default='cache/', required=False,
            help='Path to use for the cache. Default ./cache')
        cache.add_argument('--rm', '--rm-cache', default=False, required=False,
            help='Whether to delete the cache upon exit. Default False.')

        dask = parser.add_argument_group('Dask')
        dask.add_argument("-H", "--dask-host", default='localhost', required=False,
                            help="The host to use for the dask scheduler. Defaults to localhost.")
        dask.add_argument("-P", "--dask-scheduler-port", default=8786, type=int, required=False,
                            help="The port to use for the dask scheduler. Defaults to 8786.")
        return parser.parse_args()

    args = _parse_args()
    asyncio.run(start_api(**args.__dict__))


def launch_with_cluster():
    '''
    Start up a Dask Cluster and StochSS-Compute REST API. Alias to script "stochss-compute-cluster".

    `gillespy2-remote-cluster --help`
    OR
    `python -m gillespy2.remote.launch cluster --help`
    '''

    def _parse_args() -> Namespace:
        usage = '''
            gillespy2-remote-cluster -p PORT
        '''
        desc = '''
            Startup script for a StochSS-Compute cluster.
            StochSS-Compute is a server and cache that anonymizes StochSS simulation data.
            Uses Dask, a Python parallel computing library.   
        '''
        parser = ArgumentParser(description=desc, add_help=True, usage=usage,
            conflict_handler='resolve')

        server = parser.add_argument_group('Server')
        server.add_argument("-p", "--port", default=29681, type=int, required=False,
            help="The port to use for the server. Defaults to 29681.")

        cache = parser.add_argument_group('Cache')
        cache.add_argument('-c', '--cache', default='cache/', required=False,
            help='Path to use for the cache.')
        cache.add_argument('--rm', default=False, action='store_true', required=False,
            help='Whether to delete the cache upon exit. Default False.')

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
            help='A name to use when printing out the cluster, defaults to type name.')
        args =  parser.parse_args()
        return args


    args = _parse_args()

    dask_args = {}
    for (arg, value) in vars(args).items():
        if arg.startswith('dask_'):
            dask_args[arg[5:]] = value
    print('Launching Dask Cluster...')
    cluster = LocalCluster(**dask_args)
    tokens = cluster.scheduler_address.split(':')
    dask_host = tokens[1][2:]
    dask_port = int(tokens[2])
    print(f'Scheduler Address: <{cluster.scheduler_address}>')
    for i, worker in cluster.workers.items():
        print(f'Worker {i}: {worker}')

    print(f'Dashboard Link: <{cluster.dashboard_link}>\n')

    try:
        asyncio.run(start_api(port=args.port, cache=args.cache,
            dask_host=dask_host, dask_scheduler_port=dask_port, rm=args.rm))
    except asyncio.exceptions.CancelledError:
        pass
    finally:
        print('Shutting down cluster...', end='')
        asyncio.run(cluster.close())
        print('OK')

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] == 'cluster':
            del sys.argv[1]
            launch_with_cluster()
    launch_server()
