'''
stochss_compute.core.log_config

Global Logging Configuration
'''

from logging import getLogger

def init_logging(name):
    '''
    Call after import to initialize logs in a module file.
    To follow convention, use predefined __name__.

    Like so:

    from stochss_compute.core.log_config import init_logs
    log = init_logs(__name__)

    :param name: Name for the logger. Use the dot-separated module path string.
    :type name: str

    :returns: A module specific logger.
    :rtype: logging.Logger
    '''
    logger = getLogger(name)
    return logger

def set_global_log_level(level):
    '''
    Sets the root logger log level.

    :param level: NOTSET:0, DEBUG:10, INFO:20, WARNING:30, ERROR:40, CRITICAL:50, etc.
    :type level: int | logging._Level
    '''
    getLogger(__name__.split('.', maxsplit=1)[0]).setLevel(level)
    