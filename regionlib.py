#!/opt/python-3.11.5/bin/python3

# This file contains code shared by many pipeline stage scripts.
# TODO: rename from 'regionlib' to something less specific, and/or refactor
# into separate libs with more-descriptive names.

import argparse
import enum
import dataclasses
import os
import pathlib
import re
import sys
import tempfile
from typing import Self

import fiona
from osgeo import gdal
from osgeo import osr
import pyproj
import rasterio
import shapely
import shapely.geometry
import shapely.ops

from ucs_constants import *

# TODO: consider using plain string constants instead of Enum.

class Region(enum.StrEnum):
    AS = 'AS'
    EAST = 'EAST'
    GU = 'GU'
    HI = 'HI'
    PRVI = 'PRVI'
    WEST = 'WEST'

    @classmethod
    def all(cls):
        return list(cls._member_map_.values())
    
    @classmethod
    def all_string(cls):
        return list(map(lambda r: r.upper(), cls.all()))
    
REGION_NAMES = {
    Region.AS: 'American Samoa',
    Region.EAST: 'East',
    Region.GU: 'Guam', 
    Region.HI: 'Hawaii',
    Region.PRVI: 'Puerto Rico and U.S. Virgin Islands',
    Region.WEST: 'West'
}

class Severity(enum.StrEnum):
    LOW = enum.auto()
    INT_LOW = enum.auto() #'INT_LOW'
    INT = enum.auto() # 'INT'
    INT_HIGH = enum.auto() #'INT_HIGH'
    HIGH = enum.auto() #'HIGH'

    @classmethod
    def all(cls):
        return list(cls._member_map_.values())

    @classmethod
    def all_string(cls):
        return list(map(lambda r: r.lower(), cls.all()))

def choose_utm(region: str | Region) -> str:
    """Returns best NAD83 UTM zone for region."""
    
    WGS84_EPSGs = {
        Region.AS:   32702,  # 2N
        Region.EAST: 32617,  # 17N; compromise from New England to Gulf Coast
        Region.GU:   32655,  # 55S
        Region.HI:   32604,  # 4N
        Region.PRVI: 32619,
        Region.WEST: 32610   # UTM Zone 10N
    }
    
    NAD83_EPSGs = {
        # TODO: instead of 2195 for AS, should we use EPSG 6636 (NAD83(PA11)) ?
        Region.AS:   2195,   # NAD83(HARN) UTM 2N http://epsg.io/2195
        # East is a compromise, including New England and Gulf Coast.
        Region.EAST: 26917,  # NAD83 UTM 17N http://epsg.io/26917
        Region.GU:   8693,   # NAD83(MA11) UTM 55S http://epsg.io/8693
        Region.HI:   26904,  # NAD83 UTM 4N http://epsg.io/26904
        Region.PRVI: 26920,  # NAD83 UTM 20N http://epsg.io/26920
        Region.WEST: 26910  # NAD83 UTM 10N http://epsg.io/26910
    }

    if isinstance(region, str):
        region = Region(region)
    return f'epsg:{NAD83_EPSGs[region]}'


@dataclasses.dataclass
class RegionBound:
    region: str
    west: float   # longitude (WGS84)
    north: float  # latitude 
    east: float   # longitude
    south: float  # latitude 

    def contains(self, point: shapely.Point) -> bool:
        return ((point.y >= self.north and point.y <= self.south or
                 point.y >= self.south and point.y <= self.north) and
                (point.x >= self.west and point.x <= self.east))

# Crude WGS84 boundaries for partitioning our analysis.
REGION_BOUNDS = [
    RegionBound(Region.AS, -171, -14, -169, -14.6),
    RegionBound(Region.EAST, -99, 47, -66, 24),
    #RegionBound(Region.GU, 144.4, 15.4, 146, 12.9), # includes Mariana Islands
    # We have DEM for Mariana Islands, but not MHHW data.
    RegionBound(Region.GU, 144.4, 14.0, 146, 12.9), # excludes Mariana Islands
    RegionBound(Region.HI, -161, 18.5, -154.5, 22.8),
    RegionBound(Region.PRVI, -68.5, 17, -64.4, 19.5),
    RegionBound(Region.WEST, -125.3, 49.1, -115, 31.8)
]

def region_from_point(point: shapely.Point) -> Region:
    for rb in REGION_BOUNDS:
        if rb.contains(point):
            return rb.region
    raise RuntimeError(f'region not found for point {point}')

def wgs84_epsg():
    return f'EPSG:{WGS84_EPSG}'

def canonicalize_region(region: str) -> Region:
    region = region.upper()
    if region in ['PR', 'VI']:
        return Region.PRVI
    return Region(region)

def region_from_tiff_ds(ds: rasterio.io.DatasetReader) -> Region:
    box = shapely.geometry.box(*(ds.bounds))
    transform_to_wgs84 = pyproj.Transformer.from_crs(
        ds.crs, wgs84_epsg(), always_xy=True).transform
    box = shapely.ops.transform(transform_to_wgs84, box)
    return region_from_point(box.centroid)    

def region_from_tiff_fn(fn: str | pathlib.Path) -> Region:
    ds = rasterio.open(str(fn))
    return region_from_tiff_ds(ds)

@dataclasses.dataclass(frozen=True)
class ExParms:
    """Experimental Parameters. Models one 'scenario'. """
    year: int
    severity: Severity
    frequency: int  # how often does it flood annually? twice, four times, 26..

def available_flood_frequencies() -> list[int]:
    freqs = []
    for fn in TIDES_DIR.glob(f'ucs_level_*days.nc'):
        assert (m := re.match(r'ucs_level_(\d+)days.nc$', fn.name))
        freqs.append(int(m.group(1)))
    assert freqs
    return sorted(freqs)

def tide_data_path(frequency: int) -> pathlib.Path:
    assert 1 < frequency < 100
    return TIDES_DIR / f'ucs_level_{frequency:02d}days.nc'

@dataclasses.dataclass(kw_only=True)
class BaseOptions:
    regions: list[Region] = dataclasses.field(default_factory=list)
    years: list[int] = dataclasses.field(default_factory=list)
    severities: list[Severity] = dataclasses.field(default_factory=list)
    frequencies: list[int] = dataclasses.field(default_factory=list)

    @classmethod
    def setup_argparser(cls, parser: argparse.ArgumentParser):
        parser.add_argument('--region', type=str,
                            choices=Region._member_names_,
                            help='region to process. By default, all regions')
        parser.add_argument(
            '--year', type=int,
            help='year [2000 - 2100) to process. By default, each decade')
        parser.add_argument('--severity', type=str,
                            choices=list(map(lambda x: x.lower(),
                                             Severity._member_names_)),
                            help='severity to process. By default, all.')
        parser.add_argument('--frequency', type=int,
                            choices=available_flood_frequencies(),
                            help='flood frequency to process. By default, all.')

    @classmethod
    def from_args(cls, args) -> Self:
        """
          args is the return value of argparse.ArgumentParser.parse_args().
        """
        options = cls()
        if args.region:
            # Convert string to Enum Region.
            region = Region(args.region)
            if region in regions_to_skip():
                sys.stderr.write(f'region {args.region} not supported\n')
                sys.exit(1)
            options.regions = [region]
        else:
            # Default is *all* regions.
            options.regions = list(regions_to_process())

        if args.year:
            if args.year not in range(2020, 2100 + 1):
                sys.stderr.write('--year must be in [2020, 2100]\n')
                sys.exit(1)
            options.years = [args.year]
        else:
            # Default is decades between 2300 and 2100 inclusive.
            options.years = list(range(2030, 2100 + 1, 10))

        if args.severity:
            try:
                severity = Severity(args.severity)
            except ValueError:
                severities = Severity._member_names_
                sys.stderr.write(f'--severity must be one of {severities}\n')
                sys.exit(1)
            options.severities = [severity]
        else:
            options.severities = list(Severity)

        available_frequencies = available_flood_frequencies()
        if args.frequency:
            if args.frequency not in available_frequencies:
                sys.stderr.write('--frequency must be in {available}\n')
                sys.exit(1)
            options.frequencies = [args.frequency]
        else:
            # Default is all frequencies for which we have data.
            options.frequencies = available_frequencies

        return options

def regional_mhhw_uncropped_path(region: Region) -> pathlib.Path:
    return MHHW_UNCROPPED_DIR / f'mhhw-{region}-uncropped.tif'

def regional_mhhw_cropped_path(region: Region) -> pathlib.Path:
    """Cropped to align with SLR surface."""
    return MHHW_REGIONAL_DIR / f'mhhw-{region}.tif'

def regional_mhhw_path(region: Region) -> pathlib.Path:
    return regional_mhhw_cropped_path(region)

def slr_scenario_path(exparms):
    """Return a per-experimental-scenario subdirectory of SLR_DIR."""
    assert 0 < exparms.frequency < 100
    d = SLR_DIR / exparms_suffix(exparms)
    if not d.exists():
        d.mkdir(parents=True, exist_ok=True)
    return d

def output_scenario_path(exparms):
    """Return a per-experimental-scenario subdirectory of OUTPUT_DIR."""
    assert 0 < exparms.frequency < 100
    p = OUTPUT_DIR / exparms_suffix(exparms)
    if not p.exists():
        p.mkdir(parents=True, exist_ok=True)
    return p

def exparms_suffix(exparms: ExParms):
    assert 0 < exparms.frequency < 100
    return f'{exparms.year}-{exparms.severity}-{exparms.frequency:02d}'

def exparms_suffix_region(region: Region, exparms: ExParms):
    suffix = exparms_suffix(exparms)
    return f'{suffix}-{region}'

def regional_slr_path(region: Region, exparms: ExParms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'slr-{suffix}.tif'

def regional_dem_vrt_path(region: Region) -> pathlib.Path:
    return DEM_REGIONAL_DIR / f'dem-{region}.vrt'

def format_tile_name(tile_num: int, ntiles: int) -> str:
    assert ntiles < 1000
    assert tile_num < ntiles
    return f'{tile_num:03d}-of-{ntiles:03d}'

def regional_dem_tile_path(region: Region, tile_num: int, num_tiles: int
                           ) -> pathlib.Path:
    assert len(str(num_tiles)) < 1000
    tile = format_tile_name(tile_num, num_tiles)
    name = f'dem-{region}-tile-{tile}.tif'
    return DEM_REGIONAL_DIR / name

def sea_surface_path(region, exparms: ExParms) -> pathlib.Path:
    ex = exparms
    assert isinstance(region, Region)
    assert isinstance(ex.year, int)
    assert isinstance(ex.severity, Severity)
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'sea_surface-{suffix}.tif'

def flood_depth_path(region: Region, exparms: ExParms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'flood_depth-{suffix}.tif'

def flood_depth_tile_path(region: Region, exparms: ExParms, tile_num: int,
                          ntiles: int) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    tile = format_tile_name(tile_num, ntiles)
    return slr_scenario_path(exparms) / f'flood_depth-{suffix}-{tile}.tif'

def flood_depth_vrt_path(region: Region, exparms: ExParms) -> pathlib.Path:
    assert isinstance(region, Region)
    assert isinstance(exparms.year, int)
    assert isinstance(exparms.severity, Severity)
    assert isinstance(exparms.frequency, int)
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'flood_depth-{suffix}.vrt'

def is_flooded_path(region: Region, exparms: ExParms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'is_flooded-{suffix}.tif'

def is_flooded_tile_path(region: Region, exparms: ExParms,
                         tile_num: int, ntiles: int) -> pathlib.Path:
    tile = format_tile_name(tile_num, ntiles)
    suffix = exparms_suffix_region(region, exparms)
    return (slr_scenario_path(exparms) / f'is_flooded-{suffix}-{tile}.tif')

def is_flooded_vrt_path(region, exparms: ExParms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'is_flooded-tiled-{suffix}.vrt'

def is_flooded_clumped_path(region: Region, exparms: ExParms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'is_flooded_clumped-{suffix}.tif'

def slr_transect_path(exparms: ExParms) -> pathlib.Path:
    assert isinstance(exparms.severity, Severity)
    assert exparms.year >= 2000 and exparms.year <= 2100
    suffix = exparms_suffix(exparms)
    return slr_scenario_path(exparms) / f'transects_points_z-{suffix}.gpkg'

def clump_size_histogram_path(region: Region, exparms: ExParms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'clump_size_histogram-{suffix}.txt'

def flood_depth_sea_vrt_path(region, exparms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'flood_depth_sea-{suffix}.vrt'

def tile_index_path(region, exparms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'tile_index-{suffix}.gpkg'
    
def flood_depth_sea_tif_path(region, exparms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return output_scenario_path(exparms) / f'flood_depth_sea-{suffix}.tif'

def flood_depth_sea_tile_path(region: Region, exparms: ExParms,
                              tile_num: int, num_tiles: int) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    tile = format_tile_name(tile_num, num_tiles)
    d = slr_scenario_path(exparms) / f'{region}'
    if not d.exists():
        d.mkdir(parents=True, exist_ok=True)
    return d / f'flood_depth_sea-{suffix}-{tile}.tif'

def flood_depth_sea_zip_path(region: Region, exparms: ExParms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return (output_scenario_path(exparms)
            / f'flood_depth_sea-tiles-{suffix}.zip')

def is_flooded_sea_vrt_path(region, exparms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return slr_scenario_path(exparms) / f'is_flooded_sea-{suffix}.vrt'

def is_flooded_sea_tile_path(region: Region, exparms: ExParms,
                              tile_num: int, num_tiles: int) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    tile = format_tile_name(tile_num, num_tiles)
    return (slr_scenario_path(exparms)
            / f'is_flooded_sea-{suffix}-{tile}.tif')

def is_flooded_sea_tif_path(region, exparms) -> pathlib.Path:
    suffix = exparms_suffix_region(region, exparms)
    return output_scenario_path(exparms) / f'is_flooded_sea-{suffix}.tif'

def decades():
    return list(range(2030, 2100 + 1, 10))

def regions_to_skip():
    return set([Region.AS, Region.HI])

def regions_to_process():
    return sorted(set(Region) - regions_to_skip())

class TempFileWrapper:
    def __init__(self, path: pathlib.Path, suffix=None):
        self._final_path = path
        path.unlink(True)
        temp_file = tempfile.NamedTemporaryFile(
            prefix=path.stem, dir=path.parent, delete=False, suffix=suffix)
        self._temp_path = pathlib.Path(temp_file.name)
        temp_file.close()

    @property
    def temp_path(self) -> pathlib.Path:
        return self._temp_path

    @property
    def final_path(self) -> pathlib.Path:
        return self._final_path
    
    def finish(self):
        self._temp_path.rename(self._final_path)

def shapely_to_feature(shape, properties={}):
    return fiona.Feature(
        fiona.Geometry.from_dict(
            shapely.geometry.mapping(shape)),
        properties=properties)

if __name__ == '__main__':
    # TODO: write real tests.
    #for fn in pathlib.Path(DEM_UNZIP_DIR).glob('*.tif'):
    #    print(fn.stem, region_from_tiff_fn(fn))
    # print(REGION_NAMES)
    # print([r for r in Region])
    # print(regional_mhhw_path(Region.EAST))
    # print(list(Severity))
    # print(regions_to_process())
    # print(regional_dem_tile_path(Region.EAST, 2, 5))
    # print(is_flooded_clumped_path(Region.EAST, ExParms(2030, Severity.INT)))
    pass
