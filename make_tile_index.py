#!/opt/python-3.11.5/bin/python3

import argparse
import math
import os
import re
import sys
import time

import fiona
import numpy as np
from osgeo import gdal
import pathlib
import rasterio
import shapely

import regionlib
import tile_lib
from ucs_constants import *

def make_tile_index(region, exparms):
    index_path = regionlib.tile_index_path(region, exparms)
    
    vrt_path = regionlib.flood_depth_sea_vrt_path(region, exparms)
    tile_paths = tile_lib.sub_rasters_from_vrt(vrt_path)
    
    ds0 = rasterio.open(tile_paths[0])
    bbox0 = shapely.geometry.box(*ds0.bounds)
    schema = {GEOMETRY_KEY: bbox0.geom_type, PROPERTIES_KEY: {'tile': 'str'} }
    shape_file = fiona.open(index_path, 'w', crs=ds0.crs, schema=schema)
    ds0.close()

    sorted_paths = [None] * len(tile_paths)
    for path in tile_paths:
        m = re.match(r'.*(\d{3})-of-(\d{3}).tif$', str(path))
        assert m, path
        tile_num = int(m.group(1))
        ntiles = int(m.group(2))
        assert ntiles == len(tile_paths)
        sorted_paths[tile_num] = path

    for tile_num, path in enumerate(sorted_paths):
        ds = rasterio.open(path)
        bbox = shapely.geometry.box(*ds.bounds)
        feature = regionlib.shapely_to_feature(bbox, {'tile': tile_num})
        shape_file.write(feature)
    shape_file.close()
    print(f'wrote {index_path}')

def parse_options():
    parser = argparse.ArgumentParser(
        description=(
            'create an index of tiles (sub rasters) in flood_depth_sea VRT'
        )
    )
    regionlib.BaseOptions.setup_argparser(parser)
    args = parser.parse_args()
    return regionlib.BaseOptions.from_args(args)
    
def main(argv):
    gdal.UseExceptions()

    options = parse_options()
    assert len(options.years) == 1
    assert len(options.severities) == 1
    assert len(options.frequencies) == 1    
    exparms = regionlib.ExParms(
        options.years[0], options.severities[0], options.frequencies[0])

    GIGABYTE = 2 ** 30
    gdal_cachemax = 2 * GIGABYTE
    with rasterio.Env(GDAL_CACHEMAX=gdal_cachemax):
        for region in options.regions:
            make_tile_index(region, exparms)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
