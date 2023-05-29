'''
gillespy2.remote.core.messages.status
'''
from enum import Enum
from tornado.escape import json_decode
from gillespy2.remote.core.messages.base import Request, Response

class SimStatus(Enum):
    '''
    Status describing a remote simulation.
    '''
    PENDING = 'The simulation is pending.'
    RUNNING = 'The simulation is still running.'
    READY = 'Simulation is done and results exist in the cache.'
    ERROR = 'The Simulation has encountered an error.'
    DOES_NOT_EXIST = 'There is no evidence of this simulation either running or on disk.'

    @staticmethod
    def from_str(name):
        '''
        Convert str to Enum.
        '''
        if name == 'PENDING':
            return SimStatus.PENDING
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

    :param results_id: Hash of the SimulationRunRequest
    :type results_id: str
    '''
    def __init__(self, results_id):
        self.results_id = results_id
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
        return StatusRequest(request_dict['results_id'])

class StatusResponse(Response):
    '''
    A response from the server about simulation status.

    :param status: Status of the simulation
    :type status: SimStatus
    
    :param message: Possible error message or otherwise
    :type message: str
    '''
    def __init__(self, status, message = None):
        self.status = status
        self.message = message

    def encode(self):
        '''
        Encodes self.
        :returns: self as dict
        :rtype: dict[str, str]
        '''
        return {'status': self.status.name,
                'message': self.message or ''}

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
        status = SimStatus.from_str(response_dict['status'])
        message = response_dict['message']
        if not message:
            return StatusResponse(status)
        else:
            return StatusResponse(status, message)