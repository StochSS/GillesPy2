'''
stochss_compute.cloud.ec2
'''
# StochSS-Compute is a tool for running and caching GillesPy2 simulations remotely.
# Copyright (C) 2019-2023 GillesPy2 and StochSS developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import logging
from time import sleep
from secrets import token_hex
from stochss_compute.client.server import Server
from stochss_compute.cloud.ec2_config import EC2LocalConfig, EC2RemoteConfig
from stochss_compute.core.messages.source_ip import SourceIpRequest, SourceIpResponse
from stochss_compute.cloud.exceptions import EC2ImportException, ResourceException, EC2Exception
from stochss_compute.client.endpoint import Endpoint
try:
    import boto3
    from botocore.config import Config
    from botocore.session import get_session
    from botocore.exceptions import ClientError
    from paramiko import SSHClient, AutoAddPolicy
except ImportError as err:
    raise EC2ImportException from err


def _ec2_logger():
    log = logging.getLogger("EC2Cluster")
    log.setLevel(logging.INFO)
    log.propagate = False

    if not log.handlers:
        _formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        _handler = logging.StreamHandler()
        _handler.setFormatter(_formatter)
        log.addHandler(_handler)

    return log


class EC2Cluster(Server):
    """
    Attempts to load a StochSS-Compute cluster. Otherwise just initializes a new cluster.

    :param local_config: Optional. Allows configuration of local cluster resources.
    :type local_config: EC2LocalConfig

    :param remote_config: Optional. Allows configuration of remote cluster resource identifiers.
    :type remote_config: EC2RemoteConfig

    :raises EC2Exception: possible boto3 ClientError from AWS calls. See `here <https://boto3.amazonaws.com/v1/documentation/api/latest/guide/error-handling.html#aws-service-exceptions>`_.
    """
    log = _ec2_logger()

    _init = False
    _client = None
    _resources = None
    _restricted: bool = False
    _subnets = {
        'public': None,
        'private': None
    }
    _default_security_group = None
    _server_security_group = None
    _vpc = None
    _server = None
    _ami = None

    _local_config = EC2LocalConfig()
    _remote_config = EC2RemoteConfig()

    def __init__(self, local_config=None, remote_config=None) -> None:
        if local_config is not None:
            self._local_config = local_config
        if remote_config is not None:
            self._remote_config = remote_config

        if self._remote_config.region is not None:
            config = Config(region_name=self._remote_config.region)
            region = self._remote_config.region
            # Overrides any underlying configurationz
            self._client = boto3.client('ec2', config=config)
            self._resources = boto3.resource('ec2', config=config)
        else:
            region = get_session().get_config_variable('region')
            self._client = boto3.client('ec2')
            self._resources = boto3.resource('ec2')

        if self._remote_config.ami is not None:
            self._ami = self._remote_config.ami
        else:
            try:
                self._ami = self._remote_config._AMIS[region]
            except KeyError as err2:
                self._set_status('region error')
                raise EC2Exception(f'Unsupported region. Currently Supported: \
                                   {list(self._remote_config._AMIS.keys())}. \
                                   Try providing an AMI identifier.') from err2

        try:
            self._load_cluster()
        except ClientError as c_e:
            self._set_status(c_e.response['Error']['Code'])
            raise EC2Exception(c_e.response['Error']['Message']) from c_e
        except ResourceException:
            self.clean_up()

    @property
    def address(self) -> str:
        """
        The server's IP address and port.

        :returns: "http://{ip}:{port}"
        :rtype: str

        :raises EC2Exception: Do not call before launching a cluster.
        """
        if self._server is None:
            raise EC2Exception('No server found. First launch a cluster.')
        if self._server.public_ip_address is None:
            self._server.reload()
        if self._server.public_ip_address is None:
            raise EC2Exception('No public address found.')

        return f'http://{self._server.public_ip_address}:{self._remote_config.api_port}'

    @property
    def status(self) -> str:
        '''
        Return the EC2 instance status.

        :returns: A status set locally, or, if connected, a status fetched from the instance.
        :rtype: str
        '''
        if self._server is None:
            return self._status
        else:
            return self._server.state['Name']

    def _set_status(self, status):
        self._status = status
        if self._local_config.status_file is not None:
            with open(self._local_config.status_file, 'w', encoding='utf-8') as file:
                file.write(status)

    def launch_single_node_instance(self, instance_type):
        """
        Launches a single node StochSS-Compute instance. Make sure to check instance_type pricing before launching.

        :param instance_type: Example: 't3.nano' See full list `here <https://aws.amazon.com/ec2/instance-types/>`_.
        :type instance_type: str
        
        :raises EC2Exception: possible boto3 ClientError from AWS calls. See `here <https://boto3.amazonaws.com/v1/documentation/api/latest/guide/error-handling.html#aws-service-exceptions>`_.
         """
        if self._init is True:
            raise EC2Exception('You cannot launch more than one \
                               StochSS-Compute cluster instance \
                               per EC2Cluster object.')

        self._set_status('launching')
        try:
            self._launch_network()
            self._create_root_key()
            self._launch_head_node(instance_type=instance_type)
        except ClientError as c_e:
            self._set_status(c_e.response['Error']['Code'])
            raise EC2Exception(c_e.response['Error']['Message']) from c_e
        self._set_status(self._server.state['Name'])

    def clean_up(self):
        """
        Terminates and removes all cluster resources.

        :raises EC2Exception: possible boto3 ClientError from AWS calls. See `here <https://boto3.amazonaws.com/v1/documentation/api/latest/guide/error-handling.html#aws-service-exceptions>`_.
        """
        self._set_status('terminating')
        self._init = False

        vpc_search_filter = [
            {
                'Name': 'tag:Name',
                'Values': [
                    self._remote_config.vpc_name
                ]
            }
        ]
        try:
            vpc_response = self._client.describe_vpcs(
                Filters=vpc_search_filter)
            for vpc_dict in vpc_response['Vpcs']:
                vpc_id = vpc_dict['VpcId']
                vpc = self._resources.Vpc(vpc_id)
                for instance in vpc.instances.all():
                    instance.terminate()
                    self.log.info(
                        'Terminating "%s". This might take a minute.......', instance.id)
                    instance.wait_until_terminated()
                    self._server = None
                    self.log.info('Instance "%s" terminated.', instance.id)
                for s_g in vpc.security_groups.all():
                    if s_g.group_name == self._remote_config.security_group_name:
                        self.log.info('Deleting "%s".......', s_g.id)
                        s_g.delete()
                        self._server_security_group = None
                        self.log.info('Security group "%s" deleted.', s_g.id)
                    elif s_g.group_name == 'default':
                        self._default_security_group = None
                for subnet in vpc.subnets.all():
                    self.log.info('Deleting %s.......', subnet.id)
                    subnet.delete()
                    self._subnets['public'] = None
                    self.log.info('Subnet %s deleted.', subnet.id)
                for igw in vpc.internet_gateways.all():
                    self.log.info('Detaching %s.......', igw.id)
                    igw.detach_from_vpc(VpcId=vpc.vpc_id)
                    self.log.info('Gateway %s detached.', igw.id)
                    self.log.info('Deleting %s.......', igw.id)
                    igw.delete()
                    self.log.info('Gateway %s deleted.', igw.id)
                self.log.info('Deleting %s.......', vpc.id)
                vpc.delete()
                self._vpc = None
                self.log.info('VPC %s deleted.', vpc.id)
            try:
                self._client.describe_key_pairs(
                    KeyNames=[self._remote_config.key_name])
                key_pair = self._resources.KeyPair(
                    self._remote_config.key_name)
                self.log.info(
                    'Deleting "%s".', self._remote_config.key_name)
                self.log.info(
                    'Key "%s" deleted.', self._remote_config.key_name)
                key_pair.delete()
            except:
                pass
        except ClientError as c_e:
            self._set_status(c_e.response['Error']['Code'])
            raise EC2Exception(c_e.response['Error']['Message']) from c_e
        self._delete_root_key()
        self._set_status('terminated')

    def _launch_network(self):
        """
        Launches required network resources.
        """
        self.log.info("Launching Network.......")
        self._create_sssc_vpc()
        self._create_sssc_subnet(public=True)
        self._create_sssc_subnet(public=False)
        self._create_sssc_security_group()
        self._vpc.reload()

    def _create_root_key(self):
        """
        Creates a key pair for SSH login and instance launch.
        """

        response = self._client.create_key_pair(
            KeyName=self._remote_config.key_name,
            KeyType=self._local_config.key_type,
            KeyFormat=self._local_config.key_format)

        waiter = self._client.get_waiter('key_pair_exists')
        waiter.wait(KeyNames=[self._remote_config.key_name])
        os.makedirs(self._local_config.key_dir, exist_ok=True)
        with open(self._local_config.key_path, 'x', encoding='utf-8') as key:
            key.write(response['KeyMaterial'])
            os.chmod(self._local_config.key_path, 0o400)

    def _delete_root_key(self) -> None:
        """
        Deletes key from local filesystem if it exists.
        """
        if os.path.exists(self._local_config.key_path):
            self.log.info(
                'Deleting "%s".', self._local_config.key_path)
            os.remove(self._local_config.key_path)
            self.log.info('"%s" deleted.', self._local_config.key_path)

    def _create_sssc_vpc(self):
        """
        Creates a vpc.
        """
        vpc_cidr_block = '172.31.0.0/16'
        vpc_tag = [
            {
                'ResourceType': 'vpc',
                'Tags': [
                    {
                        'Key': 'Name',
                        'Value': self._remote_config.vpc_name
                    }
                ]
            }
        ]

        vpc_response = self._client.create_vpc(
            CidrBlock=vpc_cidr_block, TagSpecifications=vpc_tag)
        vpc_id = vpc_response['Vpc']['VpcId']
        vpc_waiter_exist = self._client.get_waiter('vpc_exists')
        vpc_waiter_exist.wait(VpcIds=[vpc_id])
        vpc_waiter_avail = self._client.get_waiter('vpc_available')
        vpc_waiter_avail.wait(VpcIds=[vpc_id])
        self._vpc = self._resources.Vpc(vpc_id)
        self._default_security_group = list(
            sg for sg in self._vpc.security_groups.all())[0]

        self._client.modify_vpc_attribute(
            VpcId=vpc_id, EnableDnsSupport={'Value': True})
        self._client.modify_vpc_attribute(
            VpcId=vpc_id, EnableDnsHostnames={'Value': True})

        igw_response = self._client.create_internet_gateway()
        igw_id = igw_response['InternetGateway']['InternetGatewayId']
        igw_waiter = self._client.get_waiter('internet_gateway_exists')
        igw_waiter.wait(InternetGatewayIds=[igw_id])

        self._vpc.attach_internet_gateway(InternetGatewayId=igw_id)
        for rtb in self._vpc.route_tables.all():
            if rtb.associations_attribute[0]['Main'] is True:
                rtb_id = rtb.route_table_id
        self._client.create_route(
            RouteTableId=rtb_id, GatewayId=igw_id, DestinationCidrBlock='0.0.0.0/0')

        self._vpc.reload()

    def _create_sssc_subnet(self, public: bool):
        """
        Creates a public or private subnet.
        """
        if public is True:
            label = 'public'
            subnet_cidr_block = '172.31.0.0/20'
        else:
            label = 'private'
            subnet_cidr_block = '172.31.16.0/20'

        subnet_tag = [
            {
                'ResourceType': 'subnet',
                'Tags': [
                    {
                        'Key': 'Name',
                        'Value': f'{self._remote_config.subnet_name}-{label}'
                    }
                ]
            }
        ]
        self._subnets[label] = self._vpc.create_subnet(
            CidrBlock=subnet_cidr_block, TagSpecifications=subnet_tag)
        waiter = self._client.get_waiter('subnet_available')
        waiter.wait(SubnetIds=[self._subnets[label].id])
        self._client.modify_subnet_attribute(
            SubnetId=self._subnets[label].id, MapPublicIpOnLaunch={'Value': True})
        self._subnets[label].reload()

    def _create_sssc_security_group(self):
        """
        Creates a security group for SSH and StochSS-Compute API access.
        """
        description = 'Default Security Group for StochSS-Compute.'
        self._server_security_group = self._vpc.create_security_group(
            Description=description, GroupName=self._remote_config.security_group_name)
        sshargs = {
            'CidrIp': '0.0.0.0/0',
            'FromPort': 22,
            'ToPort': 22,
            'IpProtocol': 'tcp',
        }
        self._server_security_group.authorize_ingress(**sshargs)
        sgargs = {
            'CidrIp': '0.0.0.0/0',
            'FromPort': self._remote_config.api_port,
            'ToPort': self._remote_config.api_port,
            'IpProtocol': 'tcp',
            'TagSpecifications': [
                {
                    'ResourceType': 'security-group-rule',
                    'Tags': [
                        {
                            'Key': 'Name',
                            'Value': 'api-server'
                        },
                    ]
                },
            ]
        }
        self._server_security_group.authorize_ingress(**sgargs)
        self._server_security_group.reload()

    def _restrict_ingress(self, ip_address: str = ''):
        """
        Modifies the security group API ingress rule to
        only allow access on the specified port from the given ip address.
        """
        rule_filter = [
            {
                'Name': 'group-id',
                'Values': [
                    self._server_security_group.id,
                ]
            },
            {
                'Name': 'tag:Name',
                'Values': [
                    'api-server',
                ]
            },
        ]
        sgr_response = self._client.describe_security_group_rules(
            Filters=rule_filter)
        sgr_id = sgr_response['SecurityGroupRules'][0]['SecurityGroupRuleId']
        new_sg_rules = [
            {
                'SecurityGroupRuleId': sgr_id,
                'SecurityGroupRule': {
                    'IpProtocol': 'tcp',
                    'FromPort': self._remote_config.api_port,
                    'ToPort': self._remote_config.api_port,
                    'CidrIpv4': f'{ip_address}/32',
                    'Description': 'Restricts cluster access.'
                }
            },
        ]
        self._client.modify_security_group_rules(
            GroupId=self._server_security_group.id, SecurityGroupRules=new_sg_rules)
        self._server_security_group.reload()

    def _launch_head_node(self, instance_type):
        """
        Launches a StochSS-Compute server instance.
        """
        cloud_key = token_hex(32)

        launch_commands = f'''#!/bin/bash
sudo yum update -y
sudo yum -y install docker
sudo usermod -a -G docker ec2-user
sudo service docker start
sudo chmod 666 /var/run/docker.sock 
docker run --network host --rm -t -e CLOUD_LOCK={cloud_key} --name sssc stochss/stochss-compute:cloud stochss-compute-cluster -p {self._remote_config.api_port} > /home/ec2-user/sssc-out 2> /home/ec2-user/sssc-err &
'''
        kwargs = {
            'ImageId': self._ami,
            'InstanceType': instance_type,
            'KeyName': self._remote_config.key_name,
            'MinCount': 1,
            'MaxCount': 1,
            'SubnetId': self._subnets['public'].id,
            'SecurityGroupIds': [self._default_security_group.id, self._server_security_group.id],
            'TagSpecifications': [
                {
                    'ResourceType': 'instance',
                    'Tags': [
                        {
                            'Key': 'Name',
                            'Value': self._remote_config.server_name
                        },
                    ]
                },
            ],
            'UserData': launch_commands,
        }

        self.log.info(
            'Launching StochSS-Compute server instance. This might take a minute.......')
        try:
            response = self._client.run_instances(**kwargs)
        except ClientError as c_e:
            raise EC2Exception from c_e

        instance_id = response['Instances'][0]['InstanceId']
        # try catch
        self._server = self._resources.Instance(instance_id)
        self._server.wait_until_exists()
        self._server.wait_until_running()

        self.log.info('Instance "%s" is running.', instance_id)

        self._poll_launch_progress(['sssc'])

        self.log.info('Restricting server access to only your ip.')
        source_ip = self._get_source_ip(cloud_key)

        self._restrict_ingress(source_ip)
        self._init = True
        self.log.info('StochSS-Compute ready to go!')

    def _poll_launch_progress(self, container_names, mock=False):
        """
        Polls the instance to see if the Docker container is running.

        :param container_names: A list of Docker container names to check against.
        :type container_names: List[str]
        """
        if mock is True:
            from test.unit_tests.mock_ssh import MockSSH
            ssh = MockSSH()
        else:
            ssh = SSHClient()
        ssh.set_missing_host_key_policy(AutoAddPolicy())
        sshtries = 0
        while True:
            try:
                ssh.connect(self._server.public_ip_address, username='ec2-user',
                            key_filename=self._local_config.key_path, look_for_keys=False)
                break
            except Exception as err2:
                if sshtries >= 5:
                    raise err2
                self._server.reload()
                sleep(5)
                sshtries += 1
                continue
        for container in container_names:
            sshtries = 0
            while True:
                sleep(60)
                _, stdout, stderr = ssh.exec_command(
                    "docker container inspect -f '{{.State.Running}}' " + f'{container}')
                rc = stdout.channel.recv_exit_status()
                out = stdout.readlines()
                err2 = stderr.readlines()
                if rc == -1:
                    ssh.close()
                    raise EC2Exception(
                        "Something went wrong connecting to the server. No exit status provided by the server.")
                # Wait for yum update, docker install, container download
                if rc == 1 or rc == 127:
                    self.log.info('Waiting on Docker daemon.')
                    sshtries += 1
                    if sshtries >= 5:
                        ssh.close()
                        raise EC2Exception(
                            f"Something went wrong with Docker. Max retry attempts exceeded.\nError:\n{''.join(err2)}")
                if rc == 0:
                    if 'true\n' in out:
                        sleep(10)
                        self.log.info('Container "%s" is running.', container)
                        break
        ssh.close()

    def _get_source_ip(self, cloud_key):
        """
        Ping the server to find the IP address associated with the request.

        :param cloud_key: A secret key which must match the random key used to launch the instance.
        :type cloud_key: str
        """
        source_ip_request = SourceIpRequest(cloud_key=cloud_key)
        response_raw = self.post(
            Endpoint.CLOUD, sub='/sourceip', request=source_ip_request)
        if not response_raw.ok:
            raise EC2Exception(response_raw.reason)
        response = SourceIpResponse.parse(response_raw.text)
        return response.source_ip

    def _load_cluster(self):
        '''
        Reload cluster resources. Returns False if no vpc named sssc-vpc.
        '''

        vpc_search_filter = [
            {
                'Name': 'tag:Name',
                'Values': [
                    self._remote_config.vpc_name
                ]
            }
        ]
        vpc_response = self._client.describe_vpcs(Filters=vpc_search_filter)

        if len(vpc_response['Vpcs']) == 0:
            if os.path.exists(self._local_config.key_path):
                self._set_status('key error')
                raise ResourceException
            else:
                try:
                    keypair = self._client.describe_key_pairs(
                        KeyNames=[self._remote_config.key_name])
                    if keypair is not None:
                        self._set_status('key error')
                        raise ResourceException
                except:
                    pass
            return False
        if len(vpc_response['Vpcs']) == 2:
            self.log.error('More than one VPC named "%s".',
                           self._remote_config.vpc_name)
            self._set_status('VPC error')
            raise ResourceException
        vpc_id = vpc_response['Vpcs'][0]['VpcId']
        self._vpc = self._resources.Vpc(vpc_id)
        vpc = self._vpc
        errors = False
        for instance in vpc.instances.all():
            for tag in instance.tags:
                if tag['Key'] == 'Name' and tag['Value'] == self._remote_config.server_name:
                    self._server = instance
        if self._server is None:
            self.log.warn('No instances named "%s".',
                          self._remote_config.server_name)
            self._set_status('server error')
            errors = True
        for s_g in vpc.security_groups.all():
            if s_g.group_name == 'default':
                self._default_security_group = s_g
            if s_g.group_name == self._remote_config.security_group_name:
                for rule in s_g.ip_permissions:
                    if rule['FromPort'] == self._remote_config.api_port \
                            and rule['ToPort'] == self._remote_config.api_port \
                            and rule['IpRanges'][0]['CidrIp'] == '0.0.0.0/0':
                        self.log.warn('Security Group rule error.')
                        self._set_status('security group error')
                        errors = True
                self._server_security_group = s_g
        if self._server_security_group is None:
            self.log.warn('No security group named "%s".',
                          self._remote_config.security_group_name)
            self._set_status('security group error')
            errors = True
        for subnet in vpc.subnets.all():
            for tag in subnet.tags:
                if tag['Key'] == 'Name' and tag['Value'] == f'{self._remote_config.subnet_name}-public':
                    self._subnets['public'] = subnet
                if tag['Key'] == 'Name' and tag['Value'] == f'{self._remote_config.subnet_name}-private':
                    self._subnets['private'] = subnet
        if None in self._subnets.values():
            self.log.warn('Missing or misconfigured subnet.')
            self._set_status('subnet error')
            errors = True
        if errors is True:
            raise ResourceException
        else:
            self._init = True
            self.log.info('Cluster loaded.')
            self._set_status(self._server.state['Name'])
            return True
