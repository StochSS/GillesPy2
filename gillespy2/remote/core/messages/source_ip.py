'''
gillespy2.remote.core.messages.source_ip
'''
from tornado.escape import json_decode
from gillespy2.remote.core.messages.base import Request, Response

class SourceIpRequest(Request):
    '''
    Restrict server access.

    :param cloud_key: Random key generated locally during launch.
    :type cloud_key: str
    '''

    def __init__(self, cloud_key):
        self.cloud_key = cloud_key

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
        :rtype: SourceIpRequest
        '''
        request_dict = json_decode(raw_request)
        return SourceIpRequest(request_dict['cloud_key'])

class SourceIpResponse(Response):
    '''
    Response from server containing IP address of the source.

    :param source_ip: IP address of the client.
    :type source_ip: str
    '''
    def __init__(self, source_ip):
        self.source_ip = source_ip

    def encode(self):
        '''
        :returns: self.__dict__
        :rtype: dict
        '''
        return self.__dict__

    @staticmethod
    def parse(raw_response):
        '''
        Parses a http response and returns a python object.

        :param raw_response: A raw http SourceIpResponse from the server.
        :type raw_response: str

        :returns: The decoded object.
        :rtype: SourceIpResponse
        '''
        response_dict = json_decode(raw_response)
        return SourceIpResponse(response_dict['source_ip'])
