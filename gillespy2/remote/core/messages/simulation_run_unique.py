'''
gillespy2.remote.core.messages.simulation_run_unique
'''
from secrets import token_hex
from tornado.escape import json_decode
from gillespy2 import Model
from gillespy2.remote.core.messages.base import Request, Response
from gillespy2.remote.core.messages.status import SimStatus

class SimulationRunUniqueRequest(Request):
    '''
    A one-off simulation request identifiable by a unique key.

    :param model: A model to run.
    :type model: gillespy2.Model

    :param kwargs: kwargs for the model.run() call.
    :type kwargs: dict[str, Any]
    '''
    def __init__(self, model, **kwargs):
        self.model = model
        self.kwargs = kwargs
        self.unique_key = token_hex(7)

    def encode(self):
        '''
        JSON-encode model and then encode self to dict.
        '''
        return {'model': self.model.to_json(),
                'kwargs': self.kwargs,
                'unique_key': self.unique_key,
                }

    @staticmethod
    def parse(raw_request):
        '''
        Parse raw HTTP request. Done server-side.

        :param raw_request: The request.
        :type raw_request: dict[str, str]

        :returns: The decoded object.
        :rtype: SimulationRunUniqueRequest
        '''
        request_dict = json_decode(raw_request)
        model = Model.from_json(request_dict['model'])
        kwargs = request_dict['kwargs']
        _ = SimulationRunUniqueRequest(model, **kwargs)
        _.unique_key = request_dict['unique_key'] # apply correct token (from raw request) after object construction.
        return _

class SimulationRunUniqueResponse(Response):
    '''
    A response from the server regarding a SimulationRunUniqueRequest.

    :param status: The status of the simulation.
    :type status: SimStatus

    :param error_message: Possible error message.
    :type error_message: str | None

    '''
    def __init__(self, status, error_message = None):
        self.status = status
        self.error_message = error_message

    def encode(self):
        '''
        Encode self to dict.
        '''
        return {'status': self.status.name,
                'error_message': self.error_message or '',
                }

    @staticmethod
    def parse(raw_response):
        '''
        Parse HTTP response.

        :param raw_response: The response.
        :type raw_response: dict[str, str]

        :returns: The decoded object.
        :rtype: SimulationRunResponse
        '''
        response_dict = json_decode(raw_response)
        status = SimStatus.from_str(response_dict['status'])
        error_message = response_dict['error_message']
        return SimulationRunUniqueResponse(status, error_message)
