#!/opt/python-3.11.5/bin/python3

import argparse
import concurrent.futures
import logging
import math
import multiprocessing
import os
import sys
import tempfile
import time
from typing import Iterable, Optional

import fiona
import numpy as np
from osgeo import gdal
import pathlib
import rasterio
import shapely
import tqdm
import xml.etree.ElementTree as xet

import regionlib
import tile_lib
from ucs_constants import *

PROG_NAME = 'subtract'

def parse_options():
    parser = argparse.ArgumentParser(
        description=(
            'compute sea_surface by summing MHHW and SLR data'
        )
    )
    regionlib.BaseOptions.setup_argparser(parser)
    args = parser.parse_args()
    return regionlib.BaseOptions.from_args(args)


no_data_block_cache = tile_lib.NoDataBlockCache()

def subtract_one_tile(dem_path, surf_path, tile_num, ntiles, dem_window,
                      depth_path, is_flooded_path):
    logging.info(f'subtracting tile {tile_num}')
    dem_ds = rasterio.open(dem_path)
    surf_ds = rasterio.open(surf_path)
    block_size = 128
    base_profile = dem_ds.profile | {
        'bigtiff': 'YES',
        'compress': COMPRESS_ALG,
        'driver': 'GTiff',
        TILED_KEY: True,
        'num_threads': 8,
        'ZSTD_LEVEL': ZSTD_LEVEL_SMALL,
        'width': dem_window.width,
        'height': dem_window.height,
        'transform': dem_ds.window_transform(dem_window),
        BLOCKXSIZE_KEY: block_size,
        BLOCKYSIZE_KEY: block_size
    }
    depth_profile = base_profile | {
        PREDICTOR_KEY: GDAL_PREDICTOR_FLOAT
    }
    is_flooded_profile = base_profile | {
        DTYPE_KEY: INT8_KEY,
        # NO_DATA is same as 'not flooded'. For purposes of later ('clump')
        # stage, this is good, because we don't want to distinguish unflooded
        # pixels from unknown pixels.
        NODATA_KEY: 0,
    }

    depth_tmp = regionlib.TempFileWrapper(depth_path)
    is_flooded_tmp = regionlib.TempFileWrapper(is_flooded_path)
    
    depth_ds = rasterio.open(depth_tmp.temp_path, 'w', **depth_profile)
    is_flooded_ds = rasterio.open(
        is_flooded_tmp.temp_path, 'w', nbits=1, **is_flooded_profile)
    logging.info(f'writing {depth_path.stem}')
    logging.info(f'writing {is_flooded_path.stem}')
    BAND = 1
    blocks = tile_lib.blocks_in_tile(dem_window, block_size)
    for block in blocks:
        dst_window = rasterio.windows.Window(
            block.col_off - dem_window.col_off,
            block.row_off - dem_window.row_off,
            block.width, block.height)
        bounds = dem_ds.window_bounds(block)
        surf_window = surf_ds.window(*bounds)
        shape = (block.height, block.width)
        surf_data = surf_ds.read(
            indexes=BAND, window=surf_window, masked=True,
            out_shape=shape,
            resampling=rasterio.enums.Resampling.bilinear)

        # This is complex, but saves work when one or more of the input blocks
        # is all NO_DATA. Reduces run time by 50%.
        if surf_data.mask.all():
            depth_data = no_data_block_cache.get(shape, np.float32(fNO_DATA))
            is_flooded_data = no_data_block_cache.get(shape, np.uint8(0))
        else:
            dem_data = dem_ds.read(indexes=BAND, window=block, masked=True)
            if dem_data.mask.all():
                depth_data = no_data_block_cache.get(
                    shape, np.float32(fNO_DATA))
                is_flooded_data = no_data_block_cache.get(shape, np.uint8(0))
            else:
                depth_data = np.ma.where(
                    dem_data <= surf_data, surf_data - dem_data,
                    np.float32(fNO_DATA))
                is_flooded_data = (dem_data <= surf_data).astype(np.uint8)
        
        depth_ds.write(depth_data, indexes=BAND, window=dst_window)
        is_flooded_ds.write(is_flooded_data, indexes=BAND, window=dst_window)

    depth_ds.close()
    depth_tmp.finish()
    is_flooded_ds.close()
    is_flooded_tmp.finish()
    return True

def subtract_dem_from_slr(region, exparms: regionlib.ExParms):
    depth_vrt_path = regionlib.flood_depth_vrt_path(region, exparms)
    is_flooded_vrt_path = regionlib.is_flooded_vrt_path(region, exparms)
    depth_vrt_path.unlink(True)
    is_flooded_vrt_path.unlink(True)

    dem_path = regionlib.regional_dem_vrt_path(region)
    surf_path = regionlib.sea_surface_path(region, exparms)
    dem_ds = rasterio.open(dem_path)
    surf_ds = rasterio.open(surf_path)
    assert not rasterio.coords.disjoint_bounds(dem_ds.bounds, surf_ds.bounds)

    logging.info(f'reading {dem_path}')
    logging.info(f'reading {surf_path}')
    tiles = tile_lib.nonempty_tiles_in_dem_vrt(region, TILE_SIZE)
    depth_paths = []
    is_flooded_paths = []

    executor = concurrent.futures.ProcessPoolExecutor(4)
    futures = []

    for tile_num, window in enumerate(tiles):
        depth_path = regionlib.flood_depth_tile_path(
            region, exparms, tile_num, len(tiles))
        is_flooded_path = regionlib.is_flooded_tile_path(
            region, exparms, tile_num, len(tiles))
        depth_paths.append(depth_path)
        is_flooded_paths.append(is_flooded_path)

        futures.append(executor.submit(
            subtract_one_tile,
            dem_path, surf_path, tile_num, len(tiles), window, depth_path,
            is_flooded_path))

    for tile_num, future in enumerate(
            concurrent.futures.as_completed(futures)):
        logging.info(f'completed tile {tile_num}')
        assert future.result()

    for depth_path in depth_paths:
        assert depth_path.exists()
            
    executor.shutdown()

    gdal.BuildVRT(str(depth_vrt_path), list(map(str, depth_paths)))
    logging.info(f'wrote {depth_vrt_path}')
    gdal.BuildVRT(str(is_flooded_vrt_path), list(map(str, is_flooded_paths)))
    logging.info(f'wrote {is_flooded_vrt_path}')

def main(argv):
    gdal.UseExceptions()
    tile_lib.setup_log(PROG_NAME)
    options = parse_options()

    assert len(options.severities) == 1
    assert len(options.years) == 1
    assert len(options.frequencies) == 1
    year = options.years[0]
    severity = options.severities[0]
    frequency = options.frequencies[0]
    
    for region in options.regions:
        GIGABYTE = 2 ** 30
        gdal_cachemax = 2 * GIGABYTE
        with rasterio.Env(GDAL_CACHEMAX=gdal_cachemax):
            exparms = regionlib.ExParms(year, severity, frequency)
            subtract_dem_from_slr(region, exparms)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
