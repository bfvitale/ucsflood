import pathlib

# TODO: consider changing from constants to methods / defs.
# We may want to pass these as arguments from build system (e.g. Makefile) to
# each 'stage' program.

WGS84_EPSG = 4326

UCS_DIR = pathlib.Path('/home/ben/ucs')
INPUT_DIR = UCS_DIR / 'input'
OUTPUT_DIR = UCS_DIR / 'output'
CODE_DIR = UCS_DIR / 'code'

DEM_INPUT_DIR = INPUT_DIR / 'dem/zips'
MHHW_INPUT_DIR = INPUT_DIR / 'mhhw'
TIDES_DIR = INPUT_DIR / 'tides'

TRANSECT_INPUT_DIR = INPUT_DIR / 'transect'
TRANSECT_SHAPEFILE = TRANSECT_INPUT_DIR / 'transect.gpkg'

# Path of file of polygons covering areas with DEM coverage, which are isolated
# from largest connected area in the regional DEM.
ISOLATED_SHAPEFILE = INPUT_DIR / 'isolated_areas.gpkg'

SCRATCH_DIR = UCS_DIR / 'scratch'

TRANSECT_SCRATCH_DIR = SCRATCH_DIR / 'transect'
TRANSECT_SEGMENT_FILENAME = (
    TRANSECT_SCRATCH_DIR / 'transect-segments-with-stations.gpkg')
TRANSECT_POINT_FILENAME = (
    TRANSECT_SCRATCH_DIR / 'transect-points-with-stations.gpkg')

MHHW_NOAA_GDB = SCRATCH_DIR / 'mhhw/NOAA_OCM_MHHW.gdb/'
MHHW_DIR = SCRATCH_DIR / 'mhhw'
MHHW_TIFF_DIR = MHHW_DIR / 'tiff'
MHHW_REGIONAL_DIR = MHHW_DIR / 'regional'
MHHW_UNCROPPED_DIR = MHHW_DIR / 'regional-uncropped'

DEM_DIR = SCRATCH_DIR / 'dem'
DEM_UNZIP_DIR = DEM_DIR / 'unzipped'
DEM_TMP_DIR = DEM_DIR / 'tmp'
DEM_PREPARED_DIR = DEM_DIR / 'prepared'
DEM_REGIONAL_DIR = DEM_DIR / 'regional'
DEM_BOUNDS_FN = DEM_DIR / 'dem-bounds.gpkg'

SLR_DIR = SCRATCH_DIR / 'slr'

# GDAL datatype for this constant is 'string', not 'float'!
NO_DATA_CANONICAL = '-3.40282306073709653e+38'
fNO_DATA = float(NO_DATA_CANONICAL)

# Dictionary keys; help detect typos.
# Strings used as keys in shapes' properties dictionaries.
TIDE_HEIGHT_KEY = 'tide_height'
TIDE_STATION_NAME_KEY = 'tide_station_name'
PROPERTIES_KEY = 'properties'
GEOMETRY_KEY = 'geometry'
REGION_KEY = 'region'

# Other keys used in dictionaries.
# TODO: are these available in some library, e.g. GDAL?
BLOCKXSIZE_KEY = 'blockxsize'
BLOCKYSIZE_KEY = 'blockysize'
TILED_KEY = 'tiled'
PATH_KEY = 'PATH'
DTYPE_KEY = 'dtype'
INT8_KEY = 'int8'
NODATA_KEY = 'nodata'
PREDICTOR_KEY = 'PREDICTOR'
ZSTD_LEVEL_KEY = 'ZSTD_LEVEL'

# ZSTD faster and smaller vs DEFLATE, LZMA, LZW
COMPRESS_ALG = 'ZSTD'
GDAL_PREDICTOR_NONE = 1
GDAL_PREDICTOR_DELTA = 2
GDAL_PREDICTOR_FLOAT = 3
ZSTD_LEVEL_FAST = 1  # used for temporary files
ZSTD_LEVEL_SMALL = 3
WARP_MEM_LIMIT = 2048  # MegaBytes. TODO: consider using 'percentage' setting.

TILE_SIZE = 32768
