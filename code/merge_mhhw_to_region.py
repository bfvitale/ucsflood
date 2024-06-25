#!/opt/python-3.11.5/bin/python3

"""Merge all rasters for each region into one raster per region."""

# Re-projects all rasters to a single common CRS (coordinate reference system)
# for the region, which is chosen (in regionlib.py) as the best UTM projection
# for the center of the region.

# Note, not quite the same as original permanent_inundatation code. The TIFFs
# overlap slightly. For those areas, that code took the MAXIMUM Z value. This
# code uses gdal.Warp, which chooses the Z value from the input file listed
# last.

# TODO: consider using rasterio.warp instead of gdal.Warp.
#       do rasters have any other issues? e.g. is no_data set? -9999 pixels?

import glob
import os.path
import pathlib
import sys

import fiona
from osgeo import gdal, ogr, osr
import pyproj
import rasterio
import rasterio.features
import rasterio.mask
import shapely
import shapely.geometry

import regionlib
from ucs_constants import *

def build_mhhw_index() -> dict[regionlib.Region, pathlib.Path]:
    mhhw_by_region = {}
    tiffs = list(MHHW_TIFF_DIR.glob('*.tif'))
    assert tiffs
    for path in tiffs:
        path = path
        assert path not in mhhw_by_region  # check for duplicate names
        region = regionlib.region_from_tiff_fn(path)
        mhhw_by_region.setdefault(region, []).append(path)

    with open(MHHW_DIR / 'mhhw_by_region', 'w') as index:
        for region, paths in mhhw_by_region.items():
            for path in paths:
                index.write(f'{path} {region}\n')

    return mhhw_by_region

def compute_slr_hulls() -> dict[regionlib.Region, shapely.Polygon]:
    """Compute convex hull of each region's SLR surface."""
    transects_by_region: dict[regionlib.Region, list[shapely.shape]] = {}
    transects = fiona.open(TRANSECT_SEGMENT_FILENAME)

    for feature in transects:
        shape = shapely.geometry.shape(feature.geometry)
        region = feature.properties[regionlib.REGION_KEY].upper()
        transects_by_region.setdefault(region, []).append(shape)

    hulls_by_region = {}
    for region, shapes in transects_by_region.items():
        union = shapely.union_all(shapes)
        hull = shapely.convex_hull(union)
        utm = regionlib.choose_utm(region)
        projector = pyproj.Transformer.from_crs(
            transects.crs, utm, always_xy=True, only_best=True).transform
        projected = shapely.ops.transform(projector, hull)
        hulls_by_region[region] = projected

    return hulls_by_region

"""A more flexible approach, robust to later pipeline-order changes,
would be to have the add/subtract raster steps tolerate any inputs, and
align them (e.g. crop, project, match resolution) just before the
arithmetic step. That helps decouple the SLR, MHHW, and DEM
preparation. Cropping here is a performance optimization. Not a big
deal for low-ish resolution MHHW and very low resolution SLRs, but
significant for DEMs. It also permits us to measure the area and warn if
masking discards too much data."""

def area_of_raster_hull(filename):
    ds = rasterio.open(filename)
    mask = ds.read_masks(1)
    shapes = [
        shapely.geometry.shape(poly)
        for poly, val in rasterio.features.shapes(
                mask, mask=mask, transform=ds.transform)
        if round(val) != 0
    ]
    union = shapely.union_all(shapes)
    hull = union.convex_hull
    area = hull.area
    return area

def crop_mhhw_to_slr(merged_fn, slr_polygon, cropped_fn):
    merged_ds = rasterio.open(merged_fn)

    # Crop the MHHW by an additional buffer of 500m, to ensure it fits strictly
    # within the SLR boundary, even with edge artifacts when the SLR is
    # rasterized with 500m resolution.
    # TODO: revisit this, and or extend the transects.
    slr_polygon = slr_polygon.buffer(-500)
    
    # TODO: Revisit: this reads entire raster into memory. gdal.Warp can crop
    # and mosaic in one call, and handle rasters larger than RAM.
    image, transform = rasterio.mask.mask(merged_ds, [slr_polygon], crop=True)
    dst_profile = merged_ds.profile | {
        'height': image.shape[1],
        'width': image.shape[2],
        'transform': transform        
    }
    with rasterio.open(cropped_fn, 'w', **dst_profile) as cropped_ds:
        cropped_ds.write(image)

    # Measure how much area we cropped.
    merged_area = area_of_raster_hull(merged_fn)
    cropped_area = area_of_raster_hull(cropped_fn)
    lost = 100. * (1 - cropped_area / merged_area)
    print(f'masking away {lost:.1f} %')

def main(argv):
    gdal.UseExceptions()

    hulls_by_region = compute_slr_hulls()
    mhhw_by_region = build_mhhw_index()
    for path in [MHHW_UNCROPPED_DIR, MHHW_REGIONAL_DIR]:
        path.mkdir(parents=True, exist_ok=True)
    
    for region, paths in mhhw_by_region.items():
        # TODO: We save the uncropped file, but perhaps we could delete it.
        merged_fn = regionlib.regional_mhhw_uncropped_path(region)
        
        cropped_fn = regionlib.regional_mhhw_cropped_path(region)
        crs_name = regionlib.choose_utm(region)
        
        warp_options = gdal.WarpOptions(
            dstSRS=crs_name,
            warpMemoryLimit=2000,
            creationOptions=[
                'BIGTIFF=YES',
                f'COMPRESS={COMPRESS_ALG}',
                'NUM_THREADS=ALL_CPUS',
                'PREDICTOR=3',
                'TILED=YES',
                'ZSTD_LEVEL=3'
            ],
            dstNodata=NO_DATA_CANONICAL,
            multithread=True
            # warpOptions=[NUM_THREADS=x] does not seem to help.
        )
        print(merged_fn.stem, crs_name)
        for path in paths:
            sr = osr.SpatialReference()
            sr.ImportFromWkt(gdal.Open(str(path)).GetProjection())
            assert sr.GetAuthorityName(None) == 'EPSG'
            print('\t', path.stem, f'epsg:{sr.GetAuthorityCode(None)}')
            
        gdal.Warp(str(merged_fn), list(map(str, paths)), options=warp_options)
        # TODO: consider warp_options 'cutlineDSName' and 'cropToCutline'.
        crop_mhhw_to_slr(merged_fn, hulls_by_region[region], cropped_fn)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
