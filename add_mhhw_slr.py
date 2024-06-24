#!/opt/python-3.11.5/bin/python3

# 105 minutes for all (region, year, severity, frequencies).

import argparse
import concurrent.futures
import dataclasses
import os
import sys
import tempfile
from typing import List, Optional

import fiona
import numpy as np
import rasterio
import rasterio.features
import rasterio.mask
import rasterio.merge
import shapely
import shapely.geometry
import tqdm

import regionlib
from ucs_constants import *

@dataclasses.dataclass(kw_only=True)
class Options(regionlib.BaseOptions):
    pass

    @classmethod
    def from_args(cls, args):
        options = super(Options, cls).from_args(args)
        return options

def parse_options():
    parser = argparse.ArgumentParser(
        description=(
            'compute sea_surface by summing MHHW and SLR data'
        )
    )
    regionlib.BaseOptions.setup_argparser(parser)
    args = parser.parse_args()
    return Options.from_args(args)

def slr_mostly_contains_mhhw(slr_ds, mhhw_shapes):
    slr_mask = slr_ds.read_masks(1)
    slr_shapes = [
        shapely.geometry.shape(poly)
        for poly, val in rasterio.features.shapes(
                slr_mask, mask=slr_mask, transform=slr_ds.transform)
        if round(val) != 0
    ]
    slr_union = shapely.union_all(slr_shapes)
    mhhw_union = shapely.union_all(mhhw_shapes)

    # TODO: use .contains instead of .difference. Currently does not work due
    # to edge artifacts from rasterio.features.shapes, and version of shapely
    # we are using does not yet have the 'set_precision' method.
    diff = shapely.difference(mhhw_union, slr_union)
    ratio = diff.area / mhhw_union.area
    if diff == 0 or ratio < 0.0005:
        return True

    print(f'diff.area {diff.area:.2f} ratio {ratio:.1g}')
    
    # Save diff to disk for post-mortem inspection.
    schema = {GEOMETRY_KEY: diff.geom_type, PROPERTIES_KEY: {} }
    diff_fd = fiona.open('diff.gpkg', 'w', crs=slr_ds.crs, schema=schema)
    diff_fd.write(regionlib.shapely_to_feature(diff))
    diff_fd.close()

    # Also save SLR shape.
    schema = { GEOMETRY_KEY: slr_union.geom_type, PROPERTIES_KEY: {} }
    slr_fd = fiona.open('slr-shape.gpkg', 'w', crs=slr_ds.crs, schema=schema)
    slr_fd.write(regionlib.shapely_to_feature(slr_union))
    slr_fd.close()
    
    return False

def add_mhhw_slr_one_raster(region, exparms) -> None:
    """
      Add regional rasters SLR and MHHW, placing output in regional
      'sea_surface'.

      SLR is Interpolated Tide Gauge Surface (include predicted Sea Level Rise).
      MHHW is Inundation Mapping Tidal Surface (Mean Higher High Water).

      The sum is masked to include only the extent of the MHHW raster, to
      ensure we have data only at the intersection of the two rasters.
    """
    
    slr_path = regionlib.regional_slr_path(region, exparms)
    mhhw_path = regionlib.regional_mhhw_path(region)
    sum_path = regionlib.sea_surface_path(region, exparms)
    
    slr_ds = rasterio.open(slr_path)
    mhhw_ds = rasterio.open(mhhw_path)
    CRS_KEY = 'crs'
    mhhw_crs = mhhw_ds.meta[CRS_KEY]
    
    assert mhhw_crs == slr_ds.meta[CRS_KEY]
    assert mhhw_ds.nodata == slr_ds.nodata, (
        f'mhhw no_data {mhhw_ds.nodata} != slr no_data {slr_ds.nodata}')
    res = mhhw_ds.res[0]
    assert res == mhhw_ds.res[1]
    assert (np.isclose(100, res) or np.isclose(50, res)
            or (region == regionlib.Region.GU
                and np.isclose(246, res, .01))
            ), (f'region {region} mhhw res {mhhw_ds.res}')
    assert mhhw_ds.count == 1

    # This reads entire mask into RAM. TODO: try gdal_polygonize, on a VRT which
    # extracts a 'has_data' (i.e. mask) band.
    #mhhw_mask = mhhw_ds.read_masks(1)
    #mhhw_shapes = [
    #    shapely.geometry.shape(poly)
    #    for poly, val in rasterio.features.shapes(
    #            mhhw_mask, mask=mhhw_mask, transform=mhhw_ds.transform)
    #    if round(val) != 0
    #]
    #assert slr_mostly_contains_mhhw(slr_ds, mhhw_shapes)

    assert mhhw_ds.profile[TILED_KEY]
    sum_profile = mhhw_ds.profile | {
        'compress': COMPRESS_ALG,
        'ZSTD_LEVEL': ZSTD_LEVEL_FAST,
        'PREDICTOR': GDAL_PREDICTOR_FLOAT,
        'num_threads': 8,
    }
    # numpy MaskedArray operator+ warns on overflow of masked elements.
    old_error_settings = np.seterr(over='ignore')
    
    sum_ds = rasterio.open(sum_path, 'w', **sum_profile)
    blocks = [block for ij, block in mhhw_ds.block_windows()]
    for block in blocks:
        mhhw_data = mhhw_ds.read(indexes=1, window=block, masked=True)
        slr_window = slr_ds.window(*(mhhw_ds.window_bounds(block)))
        shape = (block.height, block.width)
        slr_data = slr_ds.read(
            indexes=1, window=slr_window, out_shape=shape, masked=True,
            resampling=rasterio.enums.Resampling.bilinear)
        sum_data = mhhw_data + slr_data
        # Use in-band no_data pixels, rather than a separate .msk file.
        sum_data = np.where(sum_data.mask, sum_ds.nodata, sum_data.data)
        sum_ds.write(sum_data, window=block, indexes=1)
    sum_ds.close()
    
    np.seterr(**old_error_settings)
    
    return sum_path

def add_mhhw_slr(options):
    executor = concurrent.futures.ProcessPoolExecutor(4, max_tasks_per_child=1)
    futures = []
    
    for year in options.years:
        for severity in options.severities:
            for frequency in options.frequencies:
                exparms = regionlib.ExParms(year, severity, frequency)
                for region in options.regions:
                    futures.append(executor.submit(
                        add_mhhw_slr_one_raster, region, exparms))

    # progress = tqdm.tqdm(total=len(futures))
    for i, future in enumerate(concurrent.futures.as_completed(futures)):
        # progress.update()
        sum_path = future.result()
        print(f'{i} / {len(futures)} wrote {sum_path.stem}', flush=True)
    
def main() -> int:
    options = parse_options()
    add_mhhw_slr(options)
    return 0

if __name__ == '__main__':
    sys.exit(main())
