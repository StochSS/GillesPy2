'''
gillespy2.remote.core.messages.simulation_run
'''
from secrets import token_hex
from tornado.escape import json_decode
from gillespy2 import Model
from gillespy2.remote.core.messages.base import Request, Response
from gillespy2.remote.core.messages.status import SimStatus

class SimulationRunRequest(Request):
    '''
    A simulation request identifiable by a key.

    :param model: A model to run.
    :type model: gillespy2.Model

    :param namespace: Optional namespace for the results.
    :type namespace: str | None

    :param kwargs: kwargs for the model.run() call.
    :type kwargs: dict[str, Any]
    '''
    def __init__(self, model, namespace=None, **kwargs):
        self.model = model
        self.kwargs = kwargs
        self.key = token_hex(16)
        self.namespace = namespace

    def encode(self):
        '''
        JSON-encode model and then encode self to dict.
        '''
        return {'model': self.model.to_json(),
                'kwargs': self.kwargs,
                'key': self.key,
                'namespace': self.namespace or ''
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
        namespace = request_dict['namespace']
        _ = SimulationRunRequest(model, namespace=namespace,**kwargs)
        _.key = request_dict['key'] # apply correct token (from raw request) after object construction.
        return _

class SimulationRunResponse(Response):
    '''
    A response from the server regarding a SimulationRunRequest.

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
        return SimulationRunResponse(status, error_message)
