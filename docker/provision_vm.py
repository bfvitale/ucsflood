#!/opt/python-3.11.5/bin/python3

from collections import defaultdict
from collections.abc import Iterable
import json
import re
import requests
import sys
from typing import Any
import warnings

from google.api_core.extended_operation import ExtendedOperation
import google.auth
import google.auth.transport.requests
from google.cloud import compute_v1

PROJECT = 'ucs1-00'
GCP_REGION = 'us-central1'
ZONE = 'us-central1-a'

def wait_for_extended_operation(
    operation: ExtendedOperation, verbose_name: str = "operation",
    timeout: int = 300
) -> Any:
    result = operation.result(timeout=timeout)

    if operation.error_code:
        print(
            f'Error during {verbose_name}: [Code: {operation.error_code}]: {operation.error_message}',
            file=sys.stderr,
            flush=True,
        )
        print(f"Operation ID: {operation.name}", file=sys.stderr, flush=True)
        raise operation.exception() or RuntimeError(operation.error_message)

    if operation.warnings:
        print(f"Warnings during {verbose_name}:\n", file=sys.stderr, flush=True)
        for warning in operation.warnings:
            print(f" - {warning.code}: {warning.message}", file=sys.stderr, flush=True)

    return result

def list_all_instances(project_id: str,
                       )-> dict[str, Iterable[compute_v1.Instance]]:
    instance_client = compute_v1.InstancesClient()
    request = compute_v1.AggregatedListInstancesRequest()
    request.project = project_id
    request.max_results = 100

    agg_list = instance_client.aggregated_list(request=request)

    all_instances = defaultdict(list)
    for zone, response in agg_list:
        if response.instances:
            all_instances[zone].extend(response.instances)
    return all_instances

def new_instance_name():
    instances_by_zone = list_all_instances(PROJECT)
    names = [instance.name
             for zone, instances in instances_by_zone.items()
             for instance in instances]
    nums = []
    PREFIX='worker'
    for name in names:
        pattern = f'{PREFIX}(\d+)'
        m = re.match(pattern, name)
        if not m:
            continue
        num = int(m.group(1))
        nums.append(num)
    if len(nums) == 0:
        next_num = 0
    else:
        for next_num in range(0, max(nums) + 2):
            if next_num not in nums:
                break
    assert next_num not in nums
    assert next_num < 100
    return f'{PREFIX}{next_num:02d}'

def build_create_instance_json(vm_name):
    disk_image_name = 'projects/debian-cloud/global/images/family/debian-12'
    service_account_email = f'vmsa01@{PROJECT}.iam.gserviceaccount.com'
    
    d = {
        'machineType': f'zones/{ZONE}/machineTypes/c3d-standard-16',
        'name': vm_name,
        'disks': [
            {
                'autoDelete': True,
                'boot': True,
                'deviceName': vm_name,
                'initializeParams': {
                    'diskSizeGb': '10',
                    'diskType': f'zones/{ZONE}/diskTypes/pd-balanced',
                    'sourceImage': disk_image_name,
                },
                'mode': 'READ_WRITE',
                'type': 'PERSISTENT'
            },
            {
                'autoDelete': True,
                'boot': False,
                'deviceName': 'scratch',
                'initializeParams': {
                    'diskSizeGb': '3200',
                    'diskType': f'zones/{ZONE}/diskTypes/pd-balanced',
                },
                'mode': 'READ_WRITE',
                'type': 'PERSISTENT'
            }
        ],
        'networkInterfaces': [
            {
                'accessConfigs': [
                    {
                        'name': 'External NAT',
                        'networkTier': "Standard"
                    }
                ],
                'nicType': 'GVNIC',
                'stackType': 'IPV4_ONLY'
            }
        ],
        'scheduling': {
            'automaticRestart': False,
            'onHostMaintenance': 'MIGRATE',
            'provisioningModel': 'STANDARD',
            'maxRunDuration': {
                'seconds': 86400
            }
        },
        'serviceAccounts': [
            {
                'email': service_account_email,
                'scopes': [
                    'https://www.googleapis.com/auth/cloud-platform',
                    'https://www.googleapis.com/auth/devstorage.read_write',
                    'https://www.googleapis.com/auth/logging.write',
                    'https://www.googleapis.com/auth/monitoring.write',
                    'https://www.googleapis.com/auth/servicecontrol',
                    'https://www.googleapis.com/auth/service.management.readonly',
                    'https://www.googleapis.com/auth/trace.append'
                ]
            }
        ],
        'metadata': {
            'items': [
                {
                    'key': 'startup-script-url',
                    'value': 'gs://ucsadmin/vm_startup_script'
                }
            ]
        },
        'params': {
            'resourceManagerTags': {
                'tagKeys/281483994258228': 'tagValues/281483844204535'
            }
        }
    }
    return json.loads(json.dumps(d))

def provision_vm_rest():
    urlpath = f'/compute/v1/projects/{PROJECT}/zones/{ZONE}/instances'
    url = 'https://compute.googleapis.com' + urlpath

    vm_name = new_instance_name()
    message = json.dumps(build_create_instance_json(vm_name))
    
    creds, project = google.auth.default()
    auth_req = google.auth.transport.requests.Request()
    creds.refresh(auth_req)
    
    response = requests.post(
        url,
        headers={'Authorization': f'Bearer {creds.token}'},
        data=message
    )
    print(response.status_code)
    print(response.json())
    assert response.status_code == 200
    print(vm_name)

def main(argv):
    provision_vm_rest()
    
if __name__ == '__main__':
    sys.exit(main(sys.argv))
