#!/opt/python-3.11.5/bin/python3

import argparse
import atexit
import glob
import json
import os
import sys
import tempfile
import time
from typing import Optional

import fiona
import numpy as np
from osgeo import gdal
from osgeo import gdalconst
import pathlib
import pyproj
import rasterio
import rasterio.crs
import rasterio.io
import rasterio.mask
import rasterio.warp
import shapely
import shapely.ops

import regionlib
from ucs_constants import *

BAND = 1
NUM_THREADS = 8  # compression threads in gdal writer. Consider ALL_CPUS or %.

# TODO: factor to a core def, and separate wrappers for memfile or filename.
def set_no_data(src_path: pathlib.Path, dst_path: pathlib.Path) -> None:
    """
      This is overkill. Of 89 DEM files, only these need repair:

      Metadata nodata is -3.402823e+38, but pixels are -9999:
        CA_EKA2_GCS_5m_NAVD88m 
        CA_EKA1_GCS_5m_NAVD88m
    
      Metadata nodata unset, but pixels are -9999:
        OR_PQR2_GCS_5m_NAVD88m
        WA_SEW1_GCS_5m_NAVD88m

      Not yet checked:
        FL_Pan_East
    """
    
    sys.stderr.write(f'setting no_data {src_path.stem} ... ')
    sys.stderr.flush()
    start_time = time.monotonic()
    
    src = rasterio.open(src_path)
    if not src.profile['tiled']:
        sys.stderr.write(f'WARNING: {src_path.stem} not tiled\n')
    dst_profile = src.profile.copy()  # copy appears to be unnecessary
    dst_profile.update({
        'bigtiff': 'YES',
        'compress': COMPRESS_ALG,
        'driver': 'GTiff',
        'nodata': NO_DATA_CANONICAL,
        'tiled': True,
        'NUM_THREADS': NUM_THREADS,
        'PREDICTOR': GDAL_PREDICTOR_FLOAT,
        # speed more important than size for intermediate file
        'ZSTD_LEVEL': ZSTD_LEVEL_FAST
    })

    fNO_DATA = float(NO_DATA_CANONICAL)
    with rasterio.open(dst_path, 'w', **dst_profile) as dst:
        for ji, window in src.block_windows(BAND):
            data = src.read(BAND, window=window)
            # 49 DEMs also have -99 value pixels, which are water.
            np.place(data, np.isin(data, [-9999, -9999.9, -32767]), fNO_DATA)
            dst.write(data, indexes=BAND, window=window)
    elapsed = time.monotonic() - start_time
    sys.stderr.write(f'{elapsed:.1f} s\n')
    return None

def reproject(src_path: pathlib.Path, dst_path: pathlib.Path):
    sys.stderr.write(f're-projecting {src_path.stem} to {dst_path} ... ')
    sys.stderr.flush()
    dst_tmp = regionlib.TempFileWrapper(dst_path)

    start_time = time.monotonic()
    src = rasterio.open(src_path)

    region = regionlib.region_from_tiff_ds(src)
    dst_epsg = regionlib.choose_utm(region)
    dst_crs = rasterio.crs.CRS.from_string(dst_epsg)

    if src.crs == dst_crs:
        # TODO: link or copy?
        sys.stderr.write('WARNING: raster does not require reprojection\n')

    bounds = src.bounds
    RESOLUTION = 3
    dst_transform, dst_width, dst_height = (
        rasterio.warp.calculate_default_transform(
            src_crs=src.crs, dst_crs=dst_crs, width=src.width,
            height=src.height, left=bounds.left, bottom=bounds.bottom,
            right=bounds.right, top=bounds.top, resolution=RESOLUTION))
    
    dst_profile = src.meta.copy()
    dst_profile.update({
        'bigtiff': True,
        'compress': COMPRESS_ALG,
        'crs': dst_crs,
        'driver': 'GTiff',
        'height': dst_height,
        'width': dst_width,
        'tiled': True,
        'NUM_THREADS': NUM_THREADS,
        'PREDICTOR': GDAL_PREDICTOR_FLOAT,
        'transform': dst_transform,
        'warp_mem_limit': WARP_MEM_LIMIT,
        'ZSTD_LEVEL': ZSTD_LEVEL_FAST,
    })
    with rasterio.open(dst_tmp.temp_path, 'w', **dst_profile) as dst_ds:
        rasterio.warp.reproject(
            source=rasterio.band(src, BAND),
            destination=rasterio.band(dst_ds, BAND),
            src_crs=src.crs,
            dst_crs=dst_crs,
            src_transform=src.transform,
            num_threads=NUM_THREADS,
            # 'bilinear' is 11% slower vs default 'nearest', but better quality
            resampling=rasterio.enums.Resampling.nearest,
            warp_mem_limit=WARP_MEM_LIMIT
        )

    dst_ds.close()
    dst_tmp.finish()
    elapsed = time.monotonic() - start_time
    sys.stderr.write(f'{elapsed:.1f} s\n')

MhhwExtentType = dict[regionlib.Region, shapely.MultiPolygon]

def get_mhhw_extents() -> MhhwExtentType:
    '''Returns dict mapping Region to MultiPolygon in the regional CRS, having
       bounds of regional MHHW. Runtime ~7 secs.'''

    # TODO: Cache this in a file? rasterio.features.shapes is expensive.
    shapes_by_region = {}
    tiffs = list(MHHW_TIFF_DIR.glob('*.tif'))
    assert tiffs

    for tiff in tiffs:
        ds = rasterio.open(tiff)
        region = regionlib.region_from_tiff_ds(ds)
        regional_crs = regionlib.choose_utm(region)
        transformer = pyproj.Transformer.from_crs(ds.crs, regional_crs,
                                                  always_xy=True).transform
        mask = ds.read_masks(1)
        shapes = [
            shapely.ops.transform(transformer, shapely.geometry.shape(poly))
            for poly, val in rasterio.features.shapes(
                    mask, mask=mask, transform=ds.transform)
            if round(val) != 0
        ]
        
        if region not in shapes_by_region:
            shapes_by_region[region] = []
        shapes_by_region[region].extend(shapes)

    extent_by_region = {}
    for region, shapes in shapes_by_region.items():
        shape = shapely.union_all(shapes)
        if isinstance(shape, shapely.Polygon):
            shape = shapely.MultiPolygon([shape])
        extent_by_region[region] = shape
    return extent_by_region

def fill_multi_polygon(multi) -> shapely.MultiPolygon:
    if isinstance(multi, shapely.Polygon):
        multi = shapely.MultiPolygon([multi])
    filled = shapely.MultiPolygon(
        [shapely.Polygon(p.exterior) for p in multi.geoms])
    assert (multi - filled).area == 0.
    return filled

def crop_dem_to_mhhw_extent(mhhw_extents, src_path, cropped_path) -> None:
    sys.stderr.write(f'cropping {src_path.stem} ... ')
    sys.stderr.flush()

    cropped_tmp = regionlib.TempFileWrapper(cropped_path)

    start_time = time.monotonic()
    src = rasterio.open(src_path)
    region = regionlib.region_from_tiff_ds(src)
    assert region in mhhw_extents

    # Figure out cutline.
    mhhw_extent = mhhw_extents[region]    
    dem_bbox = shapely.geometry.box(*src.bounds)
    assert mhhw_extent.intersects(dem_bbox)
    extent = fill_multi_polygon(mhhw_extent).intersection(dem_bbox)
    
    cutline = tempfile.NamedTemporaryFile(dir='/tmp', suffix='.gpkg')
    schema = {GEOMETRY_KEY: extent.geom_type, PROPERTIES_KEY: {} }
    fh = fiona.open(cutline, 'w', crs=src.crs, schema=schema, driver='GPKG',
                    overwrite=True)
    fh.write(regionlib.shapely_to_feature(extent))
    fh.close()

    # 'extent' already in regional CRS.
    # use gdal.Warp for cropToCutline. No re-projection.
    warp_options = gdal.WarpOptions(
        creationOptions=[
            'BIGTIFF=YES',
            f'COMPRESS={COMPRESS_ALG}',
            f'NUM_THREADS={NUM_THREADS}',
            f'PREDICTOR={GDAL_PREDICTOR_FLOAT}',
            'TILED=YES',
            f'ZSTD_LEVEL={ZSTD_LEVEL_SMALL}'
        ],
        dstNodata=NO_DATA_CANONICAL,
        cutlineDSName=cutline.name,
        cropToCutline=True,
        format='GTiff',
        multithread=True,
        warpMemoryLimit=WARP_MEM_LIMIT,
        warpOptions=[]  # NUM_THREADS > 1 slows cropping of region PRVI
    )

    gdal.Warp(str(cropped_tmp.temp_path), str(src_path), options=warp_options)
    cropped_tmp.finish()
    elapsed = time.monotonic() - start_time
    sys.stderr.write(f'{elapsed:.1f} s\n')

class ProgressReporter:
    def __init__(self):
        self.start_time = time.time()
        self.last_time = None
        self.calls = 0

    def report(self, complete, message, data):
        self.calls += 1
        now = time.time()
        if (int(100 * complete) % 10 > 0 and self.last_time
            and now - self.last_time < 10):
            return
        self.last_time = now
        elapsed = now - self.start_time
        print(f'{elapsed:.1f} secs {round(complete * 100)}%')

def merge_dem_to_region(vrt_path: pathlib.Path, srcs: list[pathlib.Path]):
    """Merge all DEMs in each region into a single raster for that region.
    
       If region is set, process only that region. Otherwise, process all
       regions. Assumes inputs are from output of prepare_dem_for_merge.
    """

    # Write the VRT.
    vrt_options = gdal.BuildVRTOptions(resolution='highest', strict=True)
    sys.stderr.write(f'writing {vrt_path}\n')
    gdal.BuildVRT(str(vrt_path), sorted(map(str, srcs)), options=vrt_options)

def prepare_dem_for_merge(region: regionlib.Region,
                          mhhw_extents: MhhwExtentType, dem_path: pathlib.Path):
    """Returns path to prepared file."""
    # TODO: check uncompressed size; use memfiles for rasters that fit in RAM.
    # 10 - 15% faster.
    
    nodata_path = (DEM_TMP_DIR / dem_path.stem).with_suffix('.nodata.tif')
    warped_path = (DEM_TMP_DIR / dem_path.stem).with_suffix('.warped.tif')
    cropped_dir = DEM_REGIONAL_DIR / str(region) 
    cropped_path = (cropped_dir / dem_path.stem).with_suffix('.cropped.tif')
    if not cropped_dir.exists():
        cropped_dir.mkdir()
        
    set_no_data(dem_path, nodata_path)
    reproject(nodata_path, warped_path)
    
    # Last step is to crop. This is output to a tmp file and atomically moved
    # into place. TODO: move atomic logic to a make(1)-like wrapper.

    # TODO: crop_dem_to_mhhw_extent is implemented with gdal.Warp.
    # 'reproject' uses rasterio.warp.reproject, which internally uses
    # gdal.Warp. gdal.Warp can reproject and crop in one pass. This crop step
    # could be folded into def 'reproject' if the gdal.Warp there was exposed,
    # saving 1-2 hours of combined runtime from EAST and WEST regions.

    # TODO: should we skip this crop, and instead align all three rasters
    # (SLR, MHHW, and DEM) after we've created them?
    # The 5 most-cropped DEMs are: AL 40%, FL_Keys 28%, TX_South2 13%, NH 9%,
    # ME_East 7%. Perhaps crop only the top 5 or 10.

    crop_dem_to_mhhw_extent(mhhw_extents, warped_path, cropped_path)
    #nodata_path.unlink()
    #warped_path.unlink()
    return cropped_path

def inventory_dem() -> dict[regionlib.Region, list[pathlib.Path]]:
    dem_paths = sorted(list(DEM_UNZIP_DIR.glob('*.tif')))
    assert dem_paths
    assert len(dem_paths) == 89
    dems_by_region = {}
    for path in dem_paths:
        if path.stem.startswith('CNMI_GCS'):
            # Guam MHHW data does not extend to CNMI (Mariana Islands).
            sys.stderr.write(f'skipping DEM {path.stem}\n')
            continue
        with rasterio.open(path) as ds:
            region = regionlib.region_from_tiff_ds(ds)
            dems_by_region.setdefault(region, []).append(path)
    return dems_by_region

def process_region(mhhw_extents, region: regionlib.Region,
                   dems: list[pathlib.Path]):
    """
      'dems' should contain only the DEMs for 'region'.
    """
    GIGABYTE = 2 ** 30
    gdal_cachemax = 1 * GIGABYTE
    with rasterio.Env(GDAL_CACHEMAX=gdal_cachemax,
                      CPL_LOG_ERRORS=True):
        vrt_path = regionlib.regional_dem_vrt_path(region)
        vrt_path.unlink(True)
        prepared_dems = []
        for dem in dems:
            p = prepare_dem_for_merge(region, mhhw_extents, dem)
            prepared_dems.append(p)

        merge_dem_to_region(vrt_path, prepared_dems)

def parse_options():
    parser = argparse.ArgumentParser(description=(
        'prepare and combine all DEMs into one VRT and TIFF per region'))
    parser.add_argument('--region', type=str,
                        choices=regionlib.Region.all_string(),
                        help='region to process. By default, all regions')
    options = parser.parse_args()
    # Convert options.region from string to enum.
    if options.region:
        options.region = regionlib.Region(options.region)
    return options
        
def main(argv):
    gdal.UseExceptions()

    options = parse_options()
    DEM_REGIONAL_DIR.mkdir(parents=True, exist_ok=True)
    mhhw_extents = get_mhhw_extents()
    dems_by_region = inventory_dem()

    # If cmd-line args specify regions, only process those.
    if options.region:
        regions = {options.region}
    else:
        regions = set(dems_by_region.keys())

    for region in regions.intersection(regionlib.regions_to_skip()):
        sys.stderr.write(f'skipping {region.upper()}\n')
    regions -= regionlib.regions_to_skip()
    if not regions:
        sys.stderr.write(f'no regions to process')

    for region in regions:
        process_region(mhhw_extents, region, dems_by_region[region])

if __name__ == '__main__':
    sys.exit(main(sys.argv))
