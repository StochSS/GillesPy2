'''
gillespy2.remote.core.messages.status
'''
from enum import Enum
from tornado.escape import json_decode
from gillespy2.remote.core.exceptions import MessageParseException
from gillespy2.remote.core.messages.base import Request, Response

class SimStatus(Enum):
    '''
    Status describing a remote simulation.
    '''
    RUNNING = 'The simulation is currently running.'
    READY = 'Simulation is done and results exist in the cache.'
    ERROR = 'The Simulation has encountered an error.'
    DOES_NOT_EXIST = 'There is no evidence of this simulation either running or on disk.'

    @staticmethod
    def from_str(name):
        '''
        Convert str to Enum.
        '''
        if name == 'RUNNING':
            return SimStatus.RUNNING
        if name == 'READY':
            return SimStatus.READY
        if name == 'ERROR':
            return SimStatus.ERROR
        if name == 'DOES_NOT_EXIST':
            return SimStatus.DOES_NOT_EXIST
        # pylint: disable=no-member
        raise ValueError(f'Not a valid status.\n{SimStatus._member_names_}')
        # pylint: enable=no-member

class StatusRequest(Request):
    '''
    A request for simulation status.

    :param results_id: Hash of the SimulationRunCacheRequest or key from SimulationRunRequest
    :type results_id: str

    :param n_traj: Number of requested trajectories. Defaults to 1.
    :type n_traj: int | None

    :param task_id: Handle to a currently running simulation.
    :type task_id: str | None

    :param namespace: Optional namespace to prepend to results directory.
    :type namespace: str | None
    '''
    def __init__(self, results_id, n_traj=None, task_id=None, namespace=None):
        self.results_id = results_id
        self.namespace = namespace
        self.n_traj = n_traj
        self.task_id = task_id
        
    def encode(self):
        '''
        :returns: self.__dict__
        :rtype: dict
        '''
        return self.__dict__

    @staticmethod
    def parse(raw_request):
        '''
        Parse HTTP request.

        :param raw_request: The request.
        :type raw_request: dict[str, str]

        :returns: The decoded object.
        :rtype: StatusRequest
        '''
        request_dict = json_decode(raw_request)
        results_id = request_dict.get('results_id', None)
        namespace = request_dict.get('namespace', None)
        n_traj = request_dict.get('n_traj', None)
        task_id = request_dict.get('task_id', None)
        return StatusRequest(results_id, namespace, n_traj, task_id)

class StatusResponse(Response):
    '''
    A response from the server about simulation status.

    :param status: Status of the simulation
    :type status: SimStatus
    
    :param message: Possible error message or otherwise
    :type message: str | None
    '''
    def __init__(self, status, message=None):
        self.status = status
        self.message = message

    def encode(self):
        '''
        Encodes self.
        :returns: self as dict
        :rtype: dict
        '''
        return {'status': self.status.name,
                'message': self.message}

    @staticmethod
    def parse(raw_response):
        '''
        Parse HTTP response.

        :param raw_response: The response.
        :type raw_response: dict[str, str]

        :returns: The decoded object.
        :rtype: StatusResponse
        '''
        
        response_dict = json_decode(raw_response)
        
        try:
            status = SimStatus.from_str(response_dict.get('status'))
        except ValueError as err:
            raise MessageParseException from err
        
        message = response_dict.get('message', None)
        return StatusResponse(status, message=message)