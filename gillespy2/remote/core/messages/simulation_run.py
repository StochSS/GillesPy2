'''
gillespy2.remote.core.messages.simulation_run
'''
import copy
from hashlib import md5
from json import JSONDecodeError
from secrets import token_hex
from tornado.escape import json_decode, json_encode
from gillespy2 import Model
from gillespy2.remote.core.exceptions import MessageParseException
from gillespy2.remote.core.messages.base import Request, Response
from gillespy2.remote.core.messages.status import SimStatus

class SimulationRunRequest(Request):
    '''
    A simulation request message object.

    :param model: A model to run.
    :type model: gillespy2.Model

    :param namespace: Optional namespace for the results.
    :type namespace: str | None

    :param kwargs: kwargs for the model.run() call.
    :type kwargs: dict[str, Any]
    '''
    def __init__(self, model: Model, namespace=None, ignore_cache=False, **kwargs):
        self.model = model
        self.kwargs = kwargs
        self.id = token_hex(32)
        self.results_id = self.hash()
        self.namespace = namespace
        self.ignore_cache = ignore_cache

    def encode(self):
        '''
        JSON-encode model and then encode self to dict.
        '''

        return {
                'model': self.model.to_json(),
                'kwargs': self.kwargs,
                'id': self.id,
                'results_id': self.results_id,
                'namespace': self.namespace,
                'ignore_cache': self.ignore_cache
                }

    @staticmethod
    def parse(raw_request):
        '''
        Parse raw HTTP request. Done server-side.

        :param raw_request: The request.
        :type raw_request: dict[str, str]

        :returns: The decoded object.
        :rtype: SimulationRunRequest
        '''
        try:
            request_dict = json_decode(raw_request)
        except JSONDecodeError as err:
            raise MessageParseException from err

        model = Model.from_json(request_dict.get('model', None))
        kwargs = request_dict.get('kwargs', {})
        id = request_dict.get('id', None) # apply correct token (from raw request) after object construction.
        results_id = request_dict.get('results_id', None) # apply correct token (from raw request) after object construction.
        namespace = request_dict.get('namespace', None)
        ignore_cache = request_dict.get('ignore_cache', None)
        if None in (model, id, results_id, ignore_cache):
            raise MessageParseException        
        _ = SimulationRunRequest(model, namespace=namespace, ignore_cache=ignore_cache, **kwargs)
        _.id = id # apply correct token (from raw request) after object construction.
        return _
    
    def hash(self):
        '''
        Generate a unique hash of this simulation request.
        Does not include number_of_trajectories in this calculation.

        :returns: md5 hex digest.
        :rtype: str
        '''
        anon_model_string = self.model.to_anon().to_json(encode_private=False)
        popped_kwargs = {kw:self.kwargs[kw] for kw in self.kwargs if kw!='number_of_trajectories'}
        # Explanation of line above:
        # Take 'self.kwargs' (a dict), and add all entries to a new dictionary,
        # EXCEPT the 'number_of_trajectories' key/value pair.
        kwargs_string = json_encode(popped_kwargs)
        request_string =  f'{anon_model_string}{kwargs_string}'
        _hash = md5(str.encode(request_string)).hexdigest()
        return _hash

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
                'error_message': self.error_message}

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
