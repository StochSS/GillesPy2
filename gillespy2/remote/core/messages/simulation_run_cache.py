'''
gillespy2.remote.core.messages.simulation_run_cache
'''
from hashlib import md5
from json import JSONDecodeError
from secrets import token_hex
from tornado.escape import json_decode, json_encode
from gillespy2 import Results
from gillespy2.core.model import Model
from gillespy2.remote.core.exceptions import MessageParseException
from gillespy2.remote.core.messages.base import Request, Response
from gillespy2.remote.core.messages.status import SimStatus

from gillespy2.remote.core.utils.log_config import init_logging
log = init_logging(__name__)

class SimulationRunCacheRequest(Request):
    '''
    A simulation request message object.

    :param model: A model to run.
    :type model: gillespy2.Model

    :param namespace: Optional namespace for the results.
    :type namespace: str | None

    :param kwargs: kwargs for the model.run() call.
    :type kwargs: dict[str, Any]
    '''
    def __init__(self,
                 model,
                 namespace=None,
                 force_run=False,
                 ignore_cache=False,
                 parallelize=False,
                 chunk_trajectories=False,
                 **kwargs):
        self.model = model
        self.kwargs = kwargs
        self.id = token_hex(16)
        if ignore_cache is True:
            self.results_id = self.id
        else:
            self.results_id = self.hash()
        self.namespace = namespace
        self.force_run = force_run
        self.ignore_cache = ignore_cache
        self.parallelize = parallelize
        self.chunk_trajectories = chunk_trajectories

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
                'force_run': self.force_run,
                'ignore_cache': self.ignore_cache,
                'parallelize': self.parallelize,
                'chunk_trajectories': self.chunk_trajectories,
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
        force_run = request_dict.get('force_run', None)
        ignore_cache = request_dict.get('ignore_cache', None)
        parallelize = request_dict.get('parallelize', None)
        chunk_trajectories = request_dict.get('chunk_trajectories', None)
        if None in (model, id, results_id, force_run, ignore_cache, parallelize, chunk_trajectories):
            raise MessageParseException        
        _ = SimulationRunCacheRequest(model,
                                      namespace=namespace,
                                      force_run=force_run,
                                      ignore_cache=ignore_cache,
                                      parallelize=parallelize,
                                      chunk_trajectories=chunk_trajectories,
                                      **kwargs)
        _.id = id # apply correct token (from raw request) after object construction.
        _.results_id = results_id # apply correct token (from raw request) after object construction.
        return _
    
    def hash(self):
        '''
        Generate a unique hash of this simulation request.
        Does not include number_of_trajectories in this calculation.

        :returns: md5 hex digest.
        :rtype: str
        '''
        log.debug('SimulationRunCacheRequest.hash()...')
        anon_model_string = self.model.to_anon().to_json(encode_private=False)
        popped_kwargs = {kw:self.kwargs[kw] for kw in self.kwargs if kw!='number_of_trajectories'}
        # Explanation of line above:
        # Take 'self.kwargs' (a dict), and add all entries to a new dictionary,
        # EXCEPT the 'number_of_trajectories' key/value pair.
        kwargs_string = json_encode(popped_kwargs)
        log.debug('kwargs_string:')
        log.debug(kwargs_string)
        request_string =  f'{anon_model_string}{kwargs_string}'
        _hash = md5(str.encode(request_string)).hexdigest()
        log.debug('_hash:')
        log.debug(_hash)
        return _hash


class SimulationRunCacheResponse(Response):
    '''
    A response from the server regarding a SimulationRunCacheRequest.

    :param status: The status of the simulation.
    :type status: SimStatus

    :param error_message: Possible error message.
    :type error_message: str | None

    :param results_id: Hash of the simulation request. Identifies the results.
    :type results_id: str | None

    :param task_id: Task ID. Handle to a running simulation.
    :type task_id: str | None

    :param results: JSON-Encoded gillespy2.Results
    :type results: str | None
    '''
    def __init__(self, status, error_message = None, results_id = None, results = None, task_id = None):
        self.status = status
        self.error_message = error_message
        self.results_id = results_id
        self.results = results
        self.task_id = task_id

    def encode(self):
        '''
        Encode self to dict.
        '''
        return {'status': self.status.name,
                'error_message': self.error_message,
                'results_id': self.results_id,
                'results': self.results,
                'task_id': self.task_id,}

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
        results_id = response_dict.get('results_id', None)
        error_message = response_dict.get('error_message', None)
        task_id = response_dict.get('task_id', None)
        results = response_dict.get('results', None)
        if results is not None:
            results = Results.from_json(results)
        return SimulationRunCacheResponse(status, error_message, results_id, results, task_id)
