'''
gillespy2.remote.core.messages.results
'''
from tornado.escape import json_decode
from gillespy2 import Results

from gillespy2.remote.core.messages.base import Request, Response

class ResultsRequest(Request):
    '''
    Request results from the server.

    :param results_id: Hash of the SimulationRunCacheRequest
    :type results_id: str

    :param namespace: Optional namespace to prepend to results directory.
    :type namespace: str

    '''
    def __init__(self, results_id, namespace=None, n_traj=None):
        self.results_id = results_id
        self.namespace = namespace
        self.n_traj = n_traj

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
        :rtype: ResultsRequest
        '''
        request_dict = json_decode(raw_request)
        results_id = request_dict.get('results_id', None)
        namespace = request_dict.get('namespace', None)
        n_traj = request_dict.get('n_traj', None)
        return ResultsRequest(results_id, namespace, n_traj)

class ResultsResponse(Response):
    '''
    A response from the server about the Results.

    :param results: The requested Results from the cache. (JSON)
    :type results: str
    
    '''
    def __init__(self, results = None):
        self.results = results

    def encode(self):
        '''
        :returns: self.__dict__
        :rtype: dict
        '''
        return {'results': self.results or ''}

    @staticmethod
    def parse(raw_response):
        '''
        Parse HTTP response.

        :param raw_response: The response.
        :type raw_response: dict[str, str]

        :returns: The decoded object.
        :rtype: ResultsResponse
        '''
        response_dict = json_decode(raw_response)
        if response_dict['results'] != '':
            results = Results.from_json(response_dict['results'])
        else:
            results = None
        return ResultsResponse(results)