#!/opt/python-3.11.5/bin/python3

PROG_NAME = 'extract_sea_clumps'

"""
Outputs two rasters:

 flood_depth_max_clumps:
   Covers entire DEM area for region, but valid data only for "interesting"
   connected flooded areas. "Interesting" areas consist of the *largest* area
   from:

     - The overall region (typically, this is the area connected to the sea)

     - Each of the areas defined as polygons in shapefile 'interesting.gpkg'

 is_flooded_max_clumps:
   Just like above 'flood_depth_max_clumps', but only 1-bit per pixel,
   indicating there is flooding at that pixel.
"""

import argparse
import collections
import concurrent.futures
import dataclasses
import logging
import os
import sys
import tempfile

import fiona
from osgeo import gdal
from osgeo import gdalconst
import numpy as np
import numpy.ma
import pyproj
import rasterio
import shapely
import tqdm

import regionlib
import tile_lib
from ucs_constants import *

def histogram_of_clumped_one_tile(clumped_path, tile: rasterio.windows.Window):
    block_size = 128
    BAND = 1
 
    ds = rasterio.open(clumped_path)
    no_data = ds.nodata
    assert no_data
    
    bins_tile = np.array([0], dtype=np.int64)
    blocks = tile_lib.blocks_in_tile(tile, block_size)
    
    for block in blocks:
        data = ds.read(indexes=BAND, window=block)
        bins_block = np.bincount(data[data != no_data])
        if bins_block.size == 0:
            continue
        size = max(bins_block.size, bins_tile.size)
        bins_block.resize(size)
        bins_tile.resize(size)
        bins_tile += bins_block
        
    ds.close()
    return bins_tile

# TODO: Move this step to a separate stage, because this is slow and (a wrapper
# def) writes a file caching the output.
def histogram_of_clumped(region, exparms: regionlib.ExParms) -> dict[int, int]:
    """Returns mapping from clump_id to clump size."""
    
    clumped_path = regionlib.is_flooded_clumped_path(region, exparms)
    logging.info(f'reading {clumped_path}')
    
    # Ensure we can use DEM tiling on 'clumped' raster.
    dem_path = regionlib.regional_dem_vrt_path(region)
    dem_ds = rasterio.open(dem_path)
    clumped_ds = rasterio.open(clumped_path)
    assert dem_ds.shape == clumped_ds.shape
    dem_ds.close()
    clumped_ds.close()
    
    bins_raster = np.array([0], np.int64)

    executor = concurrent.futures.ProcessPoolExecutor(4)
    futures = []

    tiles = tile_lib.nonempty_tiles_in_dem_vrt(region, TILE_SIZE)
    for tile in tiles:
        futures.append(executor.submit(
            histogram_of_clumped_one_tile, clumped_path, tile))
        
    for tile_num, future in enumerate(
            concurrent.futures.as_completed(futures)):
        bins_tile = future.result()

        if bins_tile.size == 0:
            continue
        size = max(bins_tile.size, bins_raster.size)
        bins_tile_copy = bins_tile.copy()
        bins_tile_copy.resize(size)
        bins_raster.resize(size)
        bins_raster += bins_tile_copy
        logging.info(f'completed histogram tile {tile_num} of {len(tiles)}')

    return bins_raster

def histogram_of_clumped_cached(region, exparms: regionlib.ExParms):
    hist_path = regionlib.clump_size_histogram_path(region, exparms)
    clumped_path = regionlib.is_flooded_clumped_path(region, exparms)
    
    if (hist_path.exists()
        and hist_path.stat().st_mtime > clumped_path.stat().st_mtime):
        # Read from cache file.
        fd = open(hist_path)
        size = int(fd.readline())
        histogram = np.zeros(size, dtype=np.int64)
        for line in fd:
            fields = line.strip().split()
            assert len(fields) == 2
            clump_id, clump_size = map(int, fields)
            histogram[clump_id] = clump_size
        logging.info(f'using cached histogram {hist_path}')
        return histogram
    
    histogram = histogram_of_clumped(region, exparms)

    # Write cache file for next time.
    tmp = tempfile.NamedTemporaryFile('w', prefix=str(hist_path.parent),
                                      encoding='utf-8', delete=False)
    tmp.write(f'{histogram.size}\n')
    for clump_id in range(histogram.size):
        s = f'{clump_id} {histogram[clump_id]}'
        tmp.write(f'{s}\n')
    tmp.close()
    os.rename(tmp.name, hist_path)

    logging.info(f'wrote {hist_path}')
    return histogram

def compute_stats(path):
    ds = gdal.Open(str(path), gdal.GA_Update)
    try:
        ds.GetRasterBand(1).ComputeStatistics(False)
    except RuntimeError as e:
        if 'no valid pixels' not in e.args[0]:
            raise
    ds = None

no_data_block_cache = tile_lib.NoDataBlockCache()

def write_sea_clumps_one_tile(
        region, exparms, clumps, tile_num, ntiles, window):
    # Inputs.
    depth_in_ds = rasterio.open(regionlib.flood_depth_vrt_path(region, exparms))
    clump_ds = rasterio.open(
        regionlib.is_flooded_clumped_path(region, exparms))
    
    # Outputs.
    depth_sea_tmp = regionlib.TempFileWrapper(
        regionlib.flood_depth_sea_tile_path(region, exparms, tile_num, ntiles))
    flooded_sea_tmp = regionlib.TempFileWrapper(
        regionlib.is_flooded_sea_tile_path(region, exparms, tile_num, ntiles))
    
    block_size = 128
    profile = depth_in_ds.profile | {
        'bigtiff': 'YES',
        'compress': COMPRESS_ALG,
        ZSTD_LEVEL_KEY: ZSTD_LEVEL_SMALL,
        'driver': 'GTiff',
        TILED_KEY: True,
        'num_threads': 8,
        'width': window.width,
        'height': window.height,
        'transform': depth_in_ds.window_transform(window),
        BLOCKXSIZE_KEY: block_size,
        BLOCKYSIZE_KEY: block_size
    }
    depth_sea_profile = profile | {
        PREDICTOR_KEY: GDAL_PREDICTOR_FLOAT,
    }
    flooded_sea_profile = profile | {
        DTYPE_KEY: INT8_KEY,
        NODATA_KEY: 0,
        'photometric': 'MINISWHITE'
    }

    depth_sea_ds = rasterio.open(
        depth_sea_tmp.temp_path, 'w', **depth_sea_profile)
    flooded_sea_ds = rasterio.open(
        flooded_sea_tmp.temp_path, 'w', nbits=1, **flooded_sea_profile)  

    # Copy data, retaining only pixels with clump_id in 'clumps'.
    BAND = 1
    blocks = tile_lib.blocks_in_tile(window, block_size)

    for block in blocks:
        dst_window = rasterio.windows.Window(
            block.col_off - window.col_off, block.row_off - window.row_off,
            block.width, block.height)
        # Read input we need in all cases.
        clump_data = clump_ds.read(indexes=BAND, window=block)

        # Compute outputs, reading if needed.
        flooded_sea_data = None
        depth_sea_data = None
        if np.any(clump_data != clump_ds.nodata):
            flooded_sea_data = np.isin(clump_data, clumps).astype(np.uint8)
            if np.any(flooded_sea_data):
                depth_in_data = depth_in_ds.read(indexes=BAND, window=block)
                depth_sea_data = np.where(
                    flooded_sea_data, depth_in_data, depth_sea_ds.nodata)

        # If we did not compute output (because some input was NO_DATA),
        # fill output with NO_DATA.
        shape = (block.height, block.width)
        if flooded_sea_data is None:
            flooded_sea_data = no_data_block_cache.get(
                shape, np.uint8(flooded_sea_ds.nodata))
        if depth_sea_data is None:
            depth_sea_data = no_data_block_cache.get(
                shape, np.float32(fNO_DATA))

        # Write outputs.
        flooded_sea_ds.write(flooded_sea_data, indexes=BAND, window=dst_window)
        depth_sea_ds.write(depth_sea_data, indexes=BAND, window=dst_window)

    flooded_sea_ds.close()
    depth_sea_ds.close()
    
    compute_stats(depth_sea_tmp.temp_path)
    for tmp in [flooded_sea_tmp, depth_sea_tmp]:
        logging.info(f'wrote {tmp.final_path.stem}')
        tmp.finish()

    return flooded_sea_tmp.final_path, depth_sea_tmp.final_path

def write_sea_clumps(region, exparms, clumps):
    tiles = tile_lib.nonempty_tiles_in_dem_vrt(region, TILE_SIZE)
    ntiles = len(tiles)
    flooded_paths = []
    depth_paths = []

    flooded_vrt_path = regionlib.is_flooded_sea_vrt_path(region, exparms)
    flooded_vrt_path.unlink(True)
    depth_vrt_path = regionlib.flood_depth_sea_vrt_path(region, exparms)
    depth_vrt_path.unlink(True)
    
    executor = concurrent.futures.ProcessPoolExecutor(4)
    futures = []

    for tile_num, window in enumerate(tiles):
        futures.append(executor.submit(
            write_sea_clumps_one_tile,
            region, exparms, clumps, tile_num, len(tiles), window))

    with tqdm.tqdm(total=len(futures), position=0) as progress:
        progress.set_description('ALL')
        for future in concurrent.futures.as_completed(futures):
            progress.update()
            flooded_path, depth_path = future.result()
            flooded_paths.append(flooded_path)
            depth_paths.append(depth_path)

    executor.shutdown()

    for vrt_path, raster_paths in [(flooded_vrt_path, flooded_paths),
                                   (depth_vrt_path, depth_paths)]:
        gdal.BuildVRT(str(vrt_path), list(map(str, raster_paths)))
        #logging.info(f'computing stats {vrt_path.stem}')
        #compute_stats(vrt_path)
        #logging.info(f'wrote {vrt_path}')

def isolated_areas_for_region(features: fiona.Collection, region
                              ) -> list[shapely.Polygon]:
    """Return 'isolated_areas' for given region."""
    polys = []
    for feature in features:
        shape = shapely.geometry.shape(feature.geometry)
        if region == regionlib.region_from_point(shape.representative_point()):
            polys.append(shape)
    return polys
    
def identify_sea_clumps(region: regionlib.Region,
                                exparms: regionlib.ExParms):
    clumped_path = regionlib.is_flooded_clumped_path(region, exparms)
    clumped_ds = rasterio.open(clumped_path)

    # Compute separate histogram for each isolated area.
    
    clumps = []
    isolated = fiona.open(ISOLATED_SHAPEFILE)
    crs_region = regionlib.choose_utm(region)
    projector = pyproj.Transformer.from_crs(
        isolated.crs, crs_region, always_xy=True, only_best=True).transform

    for poly_unprojected in isolated_areas_for_region(isolated, region):
        poly = shapely.ops.transform(projector, poly_unprojected)
        window = rasterio.windows.round_window_to_full_blocks(
            clumped_ds.window(*poly.bounds),
            clumped_ds.block_shapes)
        histogram = histogram_of_clumped_one_tile(clumped_path, window)
        clumps.append(np.argmax(histogram))

    # Compute histogram for entire region.
    histogram_region = histogram_of_clumped_cached(region, exparms)
    
    clumps.append(np.argmax(histogram_region))
    return clumps

def extract_sea_clumps(region, exparms):
    clumps = identify_sea_clumps(region, exparms)
    write_sea_clumps(region, exparms, np.array(clumps))
    
def parse_options():
    parser = argparse.ArgumentParser(
        description=(
            'extract flood_depth for non-trivial connected clumps'
        )
    )
    regionlib.BaseOptions.setup_argparser(parser)
    args = parser.parse_args()
    return regionlib.BaseOptions.from_args(args)

def main(argv):
    tile_lib.setup_log(PROG_NAME)
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
            extract_sea_clumps(region, exparms)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
