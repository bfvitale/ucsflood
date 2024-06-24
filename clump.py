#!/opt/python-3.11.5/bin/python3

import argparse
import atexit
import os
import shutil
import sys
import tempfile
import time

import grass_session
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
import numpy as np
from osgeo import gdal, gdalconst
import rasterio

import regionlib
from ucs_constants import *

def clump(region, exparms: regionlib.ExParms):
    input_path = regionlib.is_flooded_vrt_path(region, exparms)
    ds = rasterio.open(input_path)
    print(f'clumping {input_path.stem}')
    
    session = grass_session.Session()
    grass_root = pathlib.Path('/home/ben/ucs/scratch/grass')

    location_path = grass_root / region
    shutil.rmtree(location_path, ignore_errors=True)
    
    session.open(gisdb=str(grass_root), location=region,
                 create_opts=ds.crs)
    grass_raster_in = input_path.stem
    start_time = time.time()
    r.external(input=str(input_path), output=grass_raster_in,
               overwrite=True, flags='r')
    print(f'r.external took {time.time() - start_time} secs')
    g.region(raster=grass_raster_in)
    output_path = regionlib.is_flooded_clumped_path(region, exparms)
    r.external_out(directory=str(output_path.parent), format='GTiff',
                   options=('BIGTIFF=YES,COMPRESS=ZSTD,ZSTD_LEVEL=3,'
                            'TILED=YES,PROFILE=GeoTIFF'))
    r.clump(input=grass_raster_in, output=output_path.name, flags='d',
            overwrite=True)

    # Fix no data value: if GRASS wrote MAX_INT32 + 1, which is out of range,
    # change to MIN_INT32.

    min_int32 = np.iinfo(np.int32).min
    ds = gdal.Open(str(output_path), gdalconst.GA_Update)
    band = ds.GetRasterBand(1)
    ds_nodata = band.GetNoDataValue()
    if isinstance(ds_nodata, float) and ds_nodata == -float(min_int32):
        band.SetNoDataValue(min_int32)
    del ds  # close the file (in gdal fashion)

    print(f'wrote {output_path}')
    return 0

def parse_options():
    parser = argparse.ArgumentParser(
        description=(
            'label connected components of is_flooded'
        )
    )
    regionlib.BaseOptions.setup_argparser(parser)
    args = parser.parse_args()
    return regionlib.BaseOptions.from_args(args)

def main(argv):
    options = parse_options()
    assert len(options.years) == 1
    assert len(options.severities) == 1
    assert len(options.frequencies) == 1

    gdal.UseExceptions()
    exparms = regionlib.ExParms(
        options.years[0], options.severities[0], options.frequencies[0])
    for region in options.regions:
        clump(region, exparms)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
