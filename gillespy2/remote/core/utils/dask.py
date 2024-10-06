import os
message = '''distributed:
  version: 2
  logging:
    distributed: error
    distributed.client: error
    distributed.worker: error
    distributed.nanny: error
    distributed.core: error
    '''
def silence_dask_logs():
    config_path_string = '~/.config/dask/distributed.yaml'
    config_path = os.path.expanduser(config_path_string)
    if os.path.exists(config_path) is True:
        print(f'To silence dask logs, edit the file {config_path} to contain:')
        print()
        print(message)
    else:
        config_dir_string = '~/.config/dask/'
        config_dir = os.path.expanduser(config_dir_string)
        if os.path.exists(config_dir) is False:
            os.makedirs(config_dir)
        with open(config_path, 'x', encoding='utf-8') as file:
            file.write(message)
            file.close()
