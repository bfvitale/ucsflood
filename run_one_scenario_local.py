#!/opt/python-3.11.5/bin/python3

import argparse
import json
import os
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

CLOUD_BUCKET = 'ucsflood'  # do not include leading 'gs://'

def parse_options():
    parser = argparse.ArgumentParser(
        description=(
            'run end-to-end'
        )
    )
    regionlib.BaseOptions.setup_argparser(parser)
    args = parser.parse_args()
    return regionlib.BaseOptions.from_args(args)

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

def upload_output(runner, exparms):
    src = regionlib.output_scenario_path(exparms)
    
    dst = f'gs://{CLOUD_BUCKET}/output/'
    print(f'copying {src} to {dst}')
    runner.run('gcloud', ['storage', 'cp', '--recursive', str(src), str(dst)])

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

class DummyStatusWriter:
    def __init(self, exparms):
        pass

    def write(self, status: str):
        print(status)
        pass
    
class StatusWriter:
    def __init__(self, exparms):
        self._client = gcs.Client()
        self._bucket = self._client.get_bucket(CLOUD_BUCKET)
        scenario_path = regionlib.output_scenario_path(exparms)
        self._path = f'output/{scenario_path.name}/status'
        self._blob = self._bucket.blob(self._path)
        self._hostname = socket.gethostname().split('.')[0]
        
    def write(self, status: str):
        record = {
            'host': self._hostname,
            'time': time.strftime('%F-%T'),
            'status': status
        }
        msg = '%s\n' % json.dumps(record)
        try:
            self._blob.upload_from_string(msg)
        except e:
            print(e)
            print(f'while updating status {msg}')
        # GCS does not like us writing the same bucket too often.
        time.sleep(1)
                
def run_one_scenario(regions, exparms):
    # status_writer = StatusWriter(exparms)
    status_writer = DummyStatusWriter()
    
    runner = Runner(exparms, status_writer)
    scenario_args = runner.scenario_args()
    scenario_region_args = scenario_args[:]
    if sorted(regions) != sorted(regionlib.regions_to_process()):
        assert len(regions) == 1
        scenario_region_args.append(f'--region={regions[0]}')

    #status_writer.write('start')
    runner.run('join_tides_transect.py')
    runner.run('gen_transects_z.py', scenario_args)
    runner.run('interpolate_slr.py', scenario_args)
    runner.run('add_mhhw_slr.py', scenario_region_args)
    runner.run('subtract.py', scenario_region_args)
    runner.run('clump.py', scenario_region_args)
    runner.run('extract_sea_clumps.py', scenario_region_args)
    runner.run('make_tile_index.py', scenario_region_args)
    # zip_flood_depth_sea(runner, regions, exparms)
    vrts_to_tifs(runner, regions, exparms)
    add_overviews(runner, regions, exparms)
    # upload_output(runner, exparms)
    #status_writer.write('done')
    
def main(argv):
    options = parse_options()
    assert len(options.years) == 1
    assert len(options.severities) == 1
    assert len(options.frequencies) == 1

    ensure_os_path('/opt/python-3.11.5/bin')  # for gdal_translate
    ensure_os_path(CODE_DIR)
    
    exparms = regionlib.ExParms(
        options.years[0], options.severities[0], options.frequencies[0])
    return run_one_scenario(options.regions, exparms)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
