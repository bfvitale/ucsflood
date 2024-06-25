#!/opt/python-3.11.5/bin/python3

# This was very hastilty assembled. See file 'TODO' for ideas for improvement.

import argparse
import dataclasses
import json
import os
import re
import subprocess
import sys
import time

from google.api_core.extended_operation import ExtendedOperation
import google.auth
import google.auth.transport.requests
from google.cloud import compute_v1
from google.cloud import storage as gcs

GCP_PROJECT = 'ucs1-00'
GCS_BUCKET = 'ucsflood'  # do not include leading 'gs://'
OUTPUT = f'gs://{GCS_BUCKET}/output'

def auth():
    creds, project = google.auth.default()
    auth_req = google.auth.transport.requests.Request()
    creds.refresh(auth_req)  # do we need to do this periodically?

@dataclasses.dataclass
class VM:
    name: str
    ip_address: str
    status: str = None
    scenario: str = None

def ssh_command(ip_address, remote_cmd):
    c = ['ssh',
         '-o', 'StrictHostKeyChecking=no',
         '-o', 'KexAlgorithms=ecdh-sha2-nistp521',
         '-p', '2222',
         ip_address ]
    c.append(remote_cmd)
    return c
    
def inventory_workers():
    instance_client = compute_v1.InstancesClient()
    request = compute_v1.AggregatedListInstancesRequest()
    request.project = GCP_PROJECT
    request.max_results = 200
    vms = []
    
    for zone, response in instance_client.aggregated_list(request=request):
        if not response.instances: continue
        for i in response.instances:
            vm = VM(
                i.name,
                status=i.status,
                ip_address=i.network_interfaces[0].access_configs[0].nat_i_p)
            vms.append(vm)
    
    for vm in vms:
        if vm.status != 'RUNNING':
            continue
        scenario = None
        try:
            out = subprocess.run(
                ssh_command(vm.ip_address, 'cat ucs/status'),
                capture_output=True)
        except subprocess.CalledProcessError as e:
            print(e)
            vm.status = 'NOT_RESPONDING'
        else:
            if (m := re.search(r'No such file or directory', str(out.stderr))):
                vm.status = 'READY'
            elif (m := re.search(r'Connection refused',
                                 str(out.stderr))):
                vm.status = 'CONN_REFUSED'
            else:
                if len(out.stdout) == 0:
                    print(out.stderr)
                record = json.loads(out.stdout)
                host = record['host']
                vm_status = record['status']
                scenario = record['scenario']
                region = record['region']
            
                assert host == f'{vm.name}-container'
            
                vm.scenario = f'{scenario}-{region}'
                if vm_status == 'done':
                    vm.status = 'READY'
                    vm.scenario = None

    def extract_num_from_name(name):
        if (m := re.search(r'\d+$', name)):
            return int(m.group(0))
        return name

    vms.sort(key=lambda v: extract_num_from_name(v.name))
    return vms

def inventory_output():
    client = gcs.Client()
    prefix = 'output'
    fields = 'items(name),nextPageToken'
    blobs = client.list_blobs(GCS_BUCKET, prefix=prefix, fields=fields)
    r = []
    for blob in blobs:
        m = re.search(r'is_flooded_sea-(.+).tif$', blob.name)
        if not m: continue
        scenario_region = m.group(1)
        r.append(scenario_region)
    return r

def read_worklist():
    """Return an ORDERED list of (scenario,region) we want."""
    
    want = []
    for line in open('worklist'):
        line = line.strip()
        fields = line.split()
        assert len(fields) == 2
        scenario, region = fields
        want.append(f'{scenario}-{region}')
    return want

def start(scenario_region, vm, dry_run):
    fields = scenario_region.split('-')
    assert len(fields) == 4
    year, severity, frequency, region = fields
    remote_cmd = (f'cd /home/ben/ucs/code; ./run_one_scenario.py '
                  f'--year={year} --severity={severity} '
                  f'--frequency={frequency} --region={region} --daemonize')
    cmd = ssh_command(vm.ip_address, remote_cmd)
    if dry_run:
        action = 'recommend starting'
    else:
        action = 'starting'
    date = time.strftime('%F %T')
    print(f'{date} {action} {scenario_region} on {vm.name}')
    if not dry_run:
        subprocess.run(cmd, check=True)                       

def main(argv):
    parser = argparse.ArgumentParser(
            description=(
                'run jobs from worklist (if not already done) on any free VMs'
            )
    )
    parser.add_argument('--dry_run', action='store_true')
    parser.add_argument('--report_vm', action='store_true')
    parser.add_argument('--report_jobs', action='store_true')
    args = parser.parse_args()

    auth()
    vms = inventory_workers()
    have = set(inventory_output())
    want = read_worklist()
    in_progress = set([vm.scenario for vm in vms if vm.status == 'RUNNING'])
    remain = [x for x in want if x not in have and x not in in_progress]
    ready_vms = [vm for vm in vms if vm.status == 'READY']

    print(f'{len(ready_vms)} VMs free {len(in_progress)} running')
    # 'have' may contain scenarios not on our list.
    print(f'{len(have)} of {len(want)} scenarios done {len(remain)} remain')

    if args.report_vm:
        for vm in vms:
            msg = f'vm {vm.name} {vm.status}'
            if vm.scenario:
                msg += f' {vm.scenario}'
            print(msg)
            
    if args.report_jobs:
        for scenario_region in have:
            print(f'scenario_region {scenario_region} done')

    if remain and ready_vms:
        start(remain[0], ready_vms[0], args.dry_run)
    
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
