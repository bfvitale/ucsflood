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

def inventory_dem() -> dict[regionlib.Region, list[pathlib.Path]]:
    dem_paths = sorted(list(DEM_UNZIP_DIR.glob('*.tif')))
    assert dem_paths
    for path in dem_paths:
        ds = rasterio.open(str(path))
        print(path.stem, ds.crs.units_factor, ds.res, ds.crs.is_geographic)
        #print(ds.res)
        #if not ds.crs.is_geographic:
        #    print(path.stem)
        
inventory_dem()

# 4.499999999999997e-05


# 2.9668008751160784

2.6949458523585537e-05
