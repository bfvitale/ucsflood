#!/opt/python-3.11.5/bin/python3

# Outputs a Sea-Level-Rise raster (TIF), from a vector file of points (GPKG).
# The output is a low-resolution (e.g. 500 meter) raster giving estimated Sea
# Level Rise in meters for a given experimental scenario (year, severity).

# The output is computed from the output using Natural Neighbor interpolation.

# The input file similarly provides the estimated Sea Level Rise for a given
# (year, severity), but only at each Tide Gauge in the region.

# The input file is typically the output of the gen_transects_z stage.
# The output file is typically fed to the add_mhhw_slr stage to yield a sea
# surface raster.

# Sample invocation:
# interpolate_slr.py \
#   --pixel_size=1000  # faster. use 500 (default) for production
#   --in_file=transect_points_z-2040-high.gpkg \
#   --out_prefix=/var/tmp/sea-level-rise-2040-high

# Outputs 240 raster files in 201 seconds, one raster for each
# (region, year, severity) combination.

import argparse
import concurrent.futures
import dataclasses
import os
import pathlib
import sys
import tempfile
import time
from typing import List, Optional

import fiona
import cgalinterp
import numpy as np
import numpy.typing as npt
import pyproj
import rasterio
import rasterio.features
import shapely
import shapely.ops
from typing import List

import regionlib
import tile_lib
from ucs_constants import *

PROG_NAME = 'interpolate_slr'

@dataclasses.dataclass(kw_only=True)
class Options(regionlib.BaseOptions):
    """
    Input and output filenames are derived by regionlib to at well-known
    locations from, as a function of 'year', 'severity', 'frequency', and
    'region'. 

    A separate output file is generated for each element in the cartesian
    product of 'regions', 'years', and 'severities', 'frequencies'.
    """
    pixel_size: Optional[int] = None

    @classmethod
    def from_args(cls, args):
        options = super(Options, cls).from_args(args)
        options.pixel_size = args.pixel_size

        return options

def setup_options() -> Options:
    parser = argparse.ArgumentParser(
        description=(
            'interpolate an SLR (Sea Level Rise) surface from transect points'
        )
    )
    
    regionlib.BaseOptions.setup_argparser(parser)
    parser.add_argument('--pixel_size', type=int, default=500,
                        help='size of each raster pixel, in meters')
    
    args = parser.parse_args()
    return Options.from_args(args)

def choose_utm(points: List[shapely.geometry.Point]) -> pyproj.CRS:

    """Given a list of lat-long points, return CRS for best UTM zone.
       Returns as follows:

        west EPSG:32610 (UTM Zone 10N)
        hi EPSG:32604 (4N)
        east EPSG:32617 (17N; compromise including New England and Gulf Coast)
        pr EPSG:32619 (19N)
        vi EPSG:32620 (20N)
        as EPSG:32702 (2N)
        gu EPSG:32655 (55S)
"""
    multi = shapely.geometry.MultiPoint(points)
    box = shapely.geometry.box(*multi.bounds)
    c = shapely.centroid(box)

    region = regionlib.region_from_point(c)
    return regionlib.choose_utm(region)

# NaturalNeighbor outputs NO_DATA for points outside convex hull, so this
# function is UNUSED.
def trim_region(utm_points, grid: npt.NDArray, transform: rasterio.Affine):
    """Returns copy of 'grid' masked by the convex hull of 'utm_points'.
    Points outside the hull are set to NO_DATA.
    'transform' could be output from rasterio.transform.from_origin()."""
    multi = shapely.geometry.MultiPoint(utm_points)
    hull = multi.convex_hull.exterior
    polygon = shapely.geometry.Polygon(hull)
    mask = rasterio.features.geometry_mask([polygon], grid.shape, transform)
    masked = np.where(mask, np.float32(fNO_DATA), grid)
    return masked

def interpolate_region(out_filename, region, utm_points, features, crs,
                       pixel_size):
    n = len(utm_points)
    assert n == len(features)
    xin = np.zeros(n, dtype=np.float32)
    yin = np.zeros(n, dtype=np.float32)
    zin = np.zeros(n, dtype=np.float32)
    for i, point in enumerate(utm_points):
        assert len(point.coords) == 1
        coord = point.coords[0]
        xin[i] = coord[0]
        yin[i] = coord[1]
        zin[i] = features[i].properties[TIDE_HEIGHT_KEY]

    # https://coast.noaa.gov/data/digitalcoast/pdf/slr-inundation-methods.pdf
    # recommends 500 - 3000 meters.
    step = pixel_size

    grid = cgalinterp.NaturalNeighbour(xin, yin, zin, step, fNO_DATA)
    assert grid.dtype == np.float32
   
    transform = rasterio.transform.from_origin(np.min(xin), np.max(yin), step,
                                               step)

    dataset = rasterio.open(
        out_filename, 'w', driver='GTiff', compress='ZSTD', tiled=True,
        predictor=GDAL_PREDICTOR_FLOAT, num_threads=4,
        zstd_level=ZSTD_LEVEL_FAST, count=1, height=grid.shape[0],
        width=grid.shape[1],
        dtype=str(grid.dtype), crs=crs, transform=transform, nodata=fNO_DATA)
    dataset.write(grid, indexes=1)
    dataset.close()
    
    return grid.shape

class Processor:
    def __init__(self, pixel_size: int):
        self._pixel_size = pixel_size

        self.in_crs = None
        self.geo_shapes_by_region = {}  # geographic coords i.e. long/lat
        self.features_by_region = {}
        self.utm_shapes_by_region = {}  # coords projected to regional UTM
        self.crs_by_region = {}

    def read_input(self, in_filename):
        """Read shapes from input file (file contains all regions).
           We leave shapes in their original coordinate reference, which we
           save in self.in_crs.
        """
        
        in_fd = fiona.open(in_filename)
        self.in_crs = in_fd.crs
        
        # TODO(vitale): instead of "parallel arrays", consider id(shape) trick.
        for feature in in_fd:
            shape = shapely.geometry.shape(feature.geometry)
            region = regionlib.region_from_point(shape)
            self.geo_shapes_by_region.setdefault(region, []).append(shape)
            self.features_by_region.setdefault(region, []).append(feature)

        assert self.geo_shapes_by_region
        num_regions = len(self.geo_shapes_by_region)
        num_features = sum(
            [len(v) for v in self.geo_shapes_by_region.values()])
        assert num_features == 1865
        assert num_regions == 6

    def project_input(self):
        """Choose best projection for each region, and project all points."""
       
        for region, shapes in self.geo_shapes_by_region.items():
            utm = choose_utm(shapes)
            self.crs_by_region[region] = utm
            projector = pyproj.Transformer.from_crs(
                self.in_crs, utm, always_xy=True, only_best=True).transform
        
            for shape in shapes:
                projected = shapely.ops.transform(projector, shape)
                self.utm_shapes_by_region.setdefault(region, []).append(
                    projected)

    def process_region(self, region, exparms):
        out_path = regionlib.regional_slr_path(region, exparms)
        start = time.time()

        with tempfile.NamedTemporaryFile(
                dir=out_path.parent, suffix='.tif') as tmp:
            
            shape = interpolate_region(
                tmp, region, self.utm_shapes_by_region[region],
                self.features_by_region[region], self.crs_by_region[region],
                self._pixel_size)

            os.replace(tmp.name, out_path)
            tmp._closer.delete = False

        end = time.time()
        print('%s %s %.1f' % (out_path.name, shape, round(end - start, 1)))

def interpolate_slr_one_year(options, year):
    for severity in options.severities:
        for frequency in options.frequencies:
            exparms = regionlib.ExParms(year, severity, frequency)
            in_path = regionlib.slr_transect_path(exparms) 
            processor = Processor(options.pixel_size)
            processor.read_input(in_path)
            processor.project_input()

            for region in options.regions:
                processor.process_region(region, exparms)

    return True

def interpolate_slr(options):
    executor = concurrent.futures.ProcessPoolExecutor(4)
    futures = []

    for year in options.years:
        futures.append(executor.submit(
            interpolate_slr_one_year, options, year))

    for future in concurrent.futures.as_completed(futures):
        # must call future.result to check for exceptions
        assert future.result()
    
    return 0
    
def main():
    options = setup_options()
    tile_lib.setup_log(PROG_NAME)
    return interpolate_slr(options)

if __name__ == '__main__':
    sys.exit(main())
