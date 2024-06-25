#!/opt/python-3.11.5/bin/python3

import os
import pathlib
import re
import sys
import zipfile

from osgeo import gdal

import regionlib
from ucs_constants import *

def check_driver():
    ver = list(map(int, (gdal.__version__).split('.')))
    assert len(ver) > 1
    major, minor = ver[0:2]
    assert major >= 3 and minor >= 7
    drivers = [gdal.GetDriver(i).ShortName
               for i in range(gdal.GetDriverCount())]
    assert 'OpenFileGDB' in drivers
    
def extract_tiff_from_gdb(in_dir, out_dir):
    dataset = gdal.Open(str(in_dir))

    options = gdal.TranslateOptions(format = 'GTiff', creationOptions = [
        'TILED=YES', 'NUM_THREADS=4', 'COMPRESS=DEFLATE', 'PREDICTOR=2'])
    subs = dataset.GetMetadata_List("SUBDATASETS")
    sub_datasets = {}
    for sub in subs:
        m = re.match(
            r'SUBDATASET_(?P<num>\d+)_(?P<fname>NAME|DESC)=(?P<fval>.*)$',
            sub)
        assert m
        if m.group('fname') == 'DESC': continue
        fval = m.group('fval')
        fields = fval.split(':')
        assert len(fields) == 3
        driver, filename, sub_name = fields
        sub_datasets[sub_name] = fval
    dataset = None  # release GDAL resources

    for sub_name, sub_path in sub_datasets.items():
        print(sub_name)
        ofn = str(out_dir / (sub_name + '.tif'))
        dataset = gdal.Open(sub_path)
        r = gdal.Translate(ofn, dataset, options=options)
        assert r

def extract_gdb_from_zip(src_zip, dst_dir):
    zf = zipfile.ZipFile(src_zip)
    zf.extractall(MHHW_DIR)

def main():
    check_driver()
    # TODO: re-cast this using Makefile or similar, so things are cached and
    # happen lazily, only when needed.

    MHHW_DIR.mkdir(parents=True, exist_ok=True)
    extract_gdb_from_zip(
        MHHW_INPUT_DIR / 'NOAA_OCM_MHHWInterpolatedSurface.zip',
        MHHW_DIR)

    MHHW_TIFF_DIR.mkdir(parents=True, exist_ok=True)
    extract_tiff_from_gdb(MHHW_NOAA_GDB, MHHW_TIFF_DIR)

    for fn in ['aev_tidal_surface_hi_coast_mhhw_msl.tif',
               'aev_tidal_surface_gu_coast_mhhw_guvd.tif']:
        print(fn)
        target = MHHW_INPUT_DIR / fn
        link = MHHW_TIFF_DIR / fn
        link.unlink(missing_ok=True)
        link.symlink_to(os.path.relpath(target, start=MHHW_TIFF_DIR))
    
if __name__ == '__main__':
    gdal.UseExceptions()
    main()
