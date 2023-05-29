'''
gillespy2.remote.core.messages.simulation_run_cache
'''
from hashlib import md5
from tornado.escape import json_decode, json_encode
from gillespy2 import Results
from gillespy2.remote.core.messages.base import Response
from gillespy2.remote.core.messages.simulation_run import SimulationRunRequest
from gillespy2.remote.core.messages.status import SimStatus

class SimulationRunCacheRequest(SimulationRunRequest):
    '''
    A simulation request.

    :param model: A model to run.
    :type model: gillespy2.Model

    :param kwargs: kwargs for the model.run() call.
    :type kwargs: dict[str, Any]
    '''
    def __init__(self, model, namespace=None, **kwargs):
        return super().__init__(model, namespace=namespace, **kwargs)

    def encode(self):
        '''
        JSON-encode model and then encode self to dict.
        '''
        return super().encode()

    @staticmethod
    def parse(raw_request):
        '''
        Parse HTTP request.

        :param raw_request: The request.
        :type raw_request: dict[str, str]

        :returns: The decoded object.
        :rtype: SimulationRunRequest
        '''
        # return SimulationRunCacheRequest()
        _ = SimulationRunRequest.parse(raw_request)
        return SimulationRunCacheRequest(_.model, namespace=_.namespace, **_.kwargs)

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
                'error_message': self.error_message or '',
                'results_id': self.results_id or '',
                'results': self.results or '',
                'task_id': self.task_id or '',}

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
        results_id = response_dict['results_id']
        error_message = response_dict['error_message']
        task_id = response_dict['task_id']
        if response_dict['results'] != '':
            results = Results.from_json(response_dict['results'])
        else:
            results = None
        return SimulationRunCacheResponse(status, error_message, results_id, results, task_id)
