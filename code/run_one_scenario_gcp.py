!/opt/python-3.11.5/bin/python3

import argparse
import dataclasses
import fcntl
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time

import google.cloud.storage as gcs

import rasterio

import regionlib
import tile_lib
from ucs_constants import *

GCS_BUCKET = 'ucsflood'  # do not include leading 'gs://'
STATUS_PATH = pathlib.Path('/home/ben/ucs/status')

@dataclasses.dataclass(kw_only=True)
class Options(regionlib.BaseOptions):
    daemonize: bool = False

    @classmethod
    def from_args(cls, args):
        options = super(Options, cls).from_args(args)
        options.daemonize = args.daemonize
        return options

def parse_options():
    parser = argparse.ArgumentParser(
        description=(
            'run end-to-end'
        )
    )
    regionlib.BaseOptions.setup_argparser(parser)
    parser.add_argument('--daemonize', action='store_true', help='daemonize')
    args = parser.parse_args()
    return Options.from_args(args)

def ensure_os_path(d: pathlib.Path):
    """Adds 'd' to os.environ['PATH'] if not already present."""
    path = os.environ[PATH_KEY].split(':')
    if str(d) not in path:
        os.environ[PATH_KEY] += f':{d}'

class Runner:
    def __init__(self, exparms, status_writer):
        self._exparms = exparms
        self._status_writer = status_writer

    def scenario_args(self):
        return [
            self.year_arg(),
            self.severity_arg(),
            self.frequency_arg()
        ]

    def year_arg(self):
        return f'--year={self._exparms.year}'

    def severity_arg(self):
        return f'--severity={self._exparms.severity}'

    def frequency_arg(self):
        return f'--frequency={self._exparms.frequency}'

    def run(self, command, args=[]):
        command_args = [command] + args
        print(f'running {" ".join(command_args)}')
        self._status_writer.write(command_args[0])
        subprocess.run(command_args, check=True)

def vrt_to_tif(runner, vrt_path, tif_path):
    ds = rasterio.open(vrt_path)
    match ds.dtypes[0]:
      case 'float32' | 'float64':
        predictor = GDAL_PREDICTOR_FLOAT
      case 'uint8' | 'int8' | 'byte' :
        predictor = GDAL_PREDICTOR_DELTA
      case _:
        predictor = GDAL_PREDICTOR_NONE
    
    runner.run(
        'gdal_translate', [
            '-co', 'COMPRESS=ZSTD',
            '-co', 'ZSTD_LEVEL=3',
            '-co', f'{PREDICTOR_KEY}={predictor}',
            '-co', 'TILED=YES',
            '-co', 'BIGTIFF=YES',
            '-co', 'NUM_THREADS=8',
            '-colorinterp', 'gray',
            '-stats',
            # Place stats in-band, rather than in 'sidecar' XML file.
            '--config', 'GDAL_PAM_ENABLED', 'NO',  
            '--config', 'GDAL_CACHEMAX', '512',
            str(vrt_path), str(tif_path)
        ]
    )
    print(f'wrote {tif_path}')

def zip_flood_depth_sea(runner, regions, exparms):
    for region in regions:
        zip_path = regionlib.flood_depth_sea_zip_path(region, exparms)
        vrt_path = regionlib.flood_depth_sea_vrt_path(region, exparms)
        tile_paths = tile_lib.sub_rasters_from_vrt(vrt_path)
        index_path = regionlib.tile_index_path(region, exparms)

        old_wd = os.getcwd()
        wd = vrt_path.parent
        os.chdir(wd)

        # TODO: handle this case by copying files to a staging area.
        assert index_path.parent == wd
        assert all([p.parent.parent == wd for p in tile_paths])

        contents = [vrt_path, index_path]
        contents.extend(sorted(tile_paths))
        rel_contents = [str(p.relative_to(wd)) for p in contents]

        tmp_dir = tempfile.TemporaryDirectory(dir=zip_path.parent)
        tmp_zip_path = pathlib.Path(tmp_dir.name) / 'tmp.zip'

        # Some of our tiles are quite compressible, despite being ZSTD TIFFs,
        # because they are all NODATA. Zip makes short work even of already-
        # compressed TIFFs, so we don't use the '-n .tif' option.
        zip_args = [str(tmp_zip_path)]
        zip_args.extend(rel_contents)
        runner.run('zip', zip_args)
        
        tmp_zip_path.rename(zip_path)
        tmp_dir.cleanup()
        os.chdir(old_wd)
        print(f'wrote {zip_path}')
        
# TODO: move this into extract_sea_clumps.py.
def vrts_to_tifs(runner, regions, exparms):
    for region in regions:
        depth_vrt_path = regionlib.flood_depth_sea_vrt_path(region, exparms)
        depth_tif_path = regionlib.flood_depth_sea_tif_path(region, exparms)
               
        flooded_vrt_path = regionlib.is_flooded_sea_vrt_path(region, exparms)
        flooded_tif_path = regionlib.is_flooded_sea_tif_path(region, exparms)

        # TODO: run these in parallel.
        vrt_to_tif(runner, depth_vrt_path, depth_tif_path)
        vrt_to_tif(runner, flooded_vrt_path, flooded_tif_path)

def upload_output(runner, regions, exparms):
    scenario = regionlib.exparms_suffix(exparms)
    dst_dir = f'gs://{GCS_BUCKET}/output/{scenario}'
    
    for region in regions:
        srcs = [
            regionlib.is_flooded_sea_tif_path(region, exparms),
            regionlib.flood_depth_sea_tif_path(region, exparms)
        ]
        ovrs = [src.with_suffix('.tif.ovr') for src in srcs]
        srcs.extend(ovrs)

        for src in srcs:
            dst = os.path.join(dst_dir, src.name)
            # TODO: use google.cloud.storage module.
            runner.run('gcloud',
                       ['storage', 'cp', str(src), dst])
            src.unlink()

def delete_scratch(runner, exparms):
    path = regionlib.slr_scenario_path(exparms)
    shutil.rmtree(path)    

def add_overview_one(runner, tif_path):
    assert tif_path.suffix == '.tif'
    overview_path = tif_path.with_suffix('.tif.ovr')
    overview_path.unlink(True)
    runner.run(
        'gdaladdo', [
            '--config', 'GDAL_CACHEMAX', '512',
            '--config', 'GDAL_NUM_THREADS', 'ALL_CPUS',
            '--config', 'COMPRESS_OVERVIEW', 'ZSTD',
            '--config', 'BIGTIFF_OVERVIEW', 'YES',
            '-r', 'bilinear',
            '-ro',
            str(tif_path)
        ])
    assert overview_path.exists()
    
def add_overviews(runner, regions, exparms):
    for region in regions:
        for tif_path in [regionlib.is_flooded_sea_tif_path(region, exparms),
                         regionlib.flood_depth_sea_tif_path(region, exparms)]:
            add_overview_one(runner, tif_path)
    
class StatusWriter:
    def __init__(self, region, exparms):
        self._client = gcs.Client()
        self._bucket = self._client.get_bucket(GCS_BUCKET)
        scenario_path = regionlib.output_scenario_path(exparms)
        self._scenario = regionlib.exparms_suffix(exparms)
        self._region = region
        self._path = f'output/{scenario_path.name}/status.{region}'
        self._blob = self._bucket.blob(self._path)
        self._hostname = socket.gethostname().split('.')[0]

    def write(self, status: str):
        record = {
            'host': self._hostname,
            'time': time.strftime('%F-%T'),
            'scenario': self._scenario,
            'region': str(self._region),
            'pid': os.getpid(),
            'status': status
        }
        msg = '%s\n' % json.dumps(record)

        # Write to local file.
        with open(STATUS_PATH, 'w') as fd:
            fd.write(msg)

        # Write to Google Cloud Storage
        try:
            self._blob.upload_from_string(msg)
        except e:
            print(e)
            print(f'while updating status {msg}')
        # GCS does not like us writing the same bucket too often.
        time.sleep(1)
                
def run_one_scenario(region, exparms):
    status_writer = StatusWriter(region, exparms)
    exclusive()
    
    runner = Runner(exparms, status_writer)
    scenario_args = runner.scenario_args()
    scenario_region_args = scenario_args[:]
    scenario_region_args.append(f'--region={region}')

    status_writer.write('start')
    runner.run('join_tides_transect.py')
    runner.run('gen_transects_z.py', scenario_args)
    runner.run('interpolate_slr.py', scenario_args)
    runner.run('add_mhhw_slr.py', scenario_region_args)
    runner.run('subtract.py', scenario_region_args)
    runner.run('clump.py', scenario_region_args)
    runner.run('extract_sea_clumps.py', scenario_region_args)
    runner.run('make_tile_index.py', scenario_region_args)
    vrts_to_tifs(runner,  [region], exparms)
    add_overviews(runner, [region], exparms)
    upload_output(runner, [region], exparms)
    delete_scratch(runner, exparms)
    status_writer.write('done')

def exclusive():
    exists = STATUS_PATH.exists()
    fd = open(STATUS_PATH, 'a+')
    if not exists:
        record = { 'pid': os.getpid() }
        fd.write(json.dumps(record))
        fd.flush()
        fd.seek(0)
        
    try:
        fcntl.lockf(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        msg = fd.read()
        status = json.loads(msg)
        pid = status['pid']
        os.stderr.write(f'pid {pid} has lock')

def daemonize():
    if os.fork():
        os._exit(0)

    os.setsid()

    if os.fork():
        os._exit(0)
        
    # Redirect stdin stdout stderr.
    ofd = os.open('/var/tmp/out', os.O_RDWR | os.O_APPEND | os.O_CREAT)
    os.dup2(ofd, 1)
    os.dup2(ofd, 2)
    os.close(ofd)
    
    ifd = os.open(os.devnull, os.O_RDONLY)
    os.dup2(ifd, 0)
    os.close(ifd)

def main(argv):
    options = parse_options()
    assert len(options.years) == 1
    assert len(options.severities) == 1
    assert len(options.frequencies) == 1
    assert len(options.regions) == 1

    ensure_os_path('/opt/python-3.11.5/bin')  # for gdal_translate
    ensure_os_path(CODE_DIR)

    if options.daemonize:
        daemonize()

    exparms = regionlib.ExParms(
        options.years[0], options.severities[0], options.frequencies[0])
    return run_one_scenario(options.regions[0], exparms)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
