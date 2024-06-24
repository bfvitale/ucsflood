#!/opt/python-3.11.5/bin/python3

import argparse
import atexit
import dataclasses
import glob
import itertools
import json
import math
import os
import sys
import tempfile
import time
import tqdm
from typing import Optional
import xml.etree.ElementTree as xet

import affine
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
import rasterio.sample
import rasterio.warp
import shapely
import shapely.ops

import regionlib
from ucs_constants import *

orig_vrt = '/home/ben/ucs/scratch/slr/flood_depth-2030-int-EAST.tif'
tiled_vrt = '/home/ben/ucs/scratch/slr/flood_depth-2030-int-EAST.vrt'

@dataclasses.dataclass
class SubRaster:
    fn: str
    win: rasterio.windows.Window
    bounds: shapely.Polygon
    transform: affine.Affine
    vrt_transform: affine.Affine
    ds: rasterio.DatasetReader

def sub_rasters_from_vrt(vrt_path):
    vrt_ds = rasterio.open(vrt_path)
    
    tree = xet.parse(vrt_path)
    sub_rasters = {}
    for elem in tree.getroot().iter('ComplexSource'):
        fn = elem.find('SourceFilename').text
        if not pathlib.Path(fn).is_absolute():
            fn = pathlib.Path(vrt_path).parent.joinpath(fn)

        rect = elem.find('DstRect')
        win = rasterio.windows.Window(math.ceil(float(rect.get('xOff'))),
                                      float(rect.get('yOff')),
                                      float(rect.get('xSize')),
                                      float(rect.get('ySize')))
        sub_ds = rasterio.open(fn)
        sub = SubRaster(fn, win, shapely.box(*sub_ds.bounds), sub_ds.transform,
                        vrt_ds.transform, sub_ds)
        sub_rasters[fn] = sub
    return sub_rasters
    
def find_sub_raster(xy, sub_rasters):
    p = shapely.Point(*xy)
    r = []
    for sub in sub_rasters.values():
        if sub.bounds.contains(p):
            r.append(sub)
    assert len(r) == 1
    return r[0]

#orig_sub_rasters = sub_rasters_from_vrt(orig_vrt)
tiled_sub_rasters = sub_rasters_from_vrt(tiled_vrt)

orig_ds = rasterio.open(orig_vrt)
tiled_ds = rasterio.open(tiled_vrt)

# EAST mismatch 53815125 /74521520 block_ij 5748, 2906

window_list1 = list(orig_ds.block_windows(1))
# window_list1 = window_list1[7214984:]

# block_ij, w = window_list1[i]
# for i in range(len(window_list1)):
#for block_ij, w in itertools.islice(orig_ds.block_windows(1), 8066675, None):

progress_bar = tqdm.tqdm(total=len(window_list1))
for block_ij, w in window_list1:
    progress_bar.update()
    orig_data = orig_ds.read(1, window=w)
    tiled_data = tiled_ds.read(1, window=w)
    if not np.array_equal(orig_data, tiled_data):
        print(f'mismatch block_ij {block_ij}')
        iter = np.nditer(orig_data, flags=['multi_index'])
        for orig_pixel in iter:
            tiled_pixel = tiled_data[iter.multi_index]
            row = w.row_off + iter.multi_index[0]
            col = w.col_off + iter.multi_index[1]
            #w2 = rasterio.windows.Window(col, row, 1, 1)
            #orig_pixel_check = orig_ds.read(1, window=w2)[0][0]
            #tiled_pixel_check = tiled_ds.read(1, window=w2)[0][0]
            #assert (orig_pixel == orig_pixel_check
            #        and tiled_pixel == tiled_pixel_check)
            orig_xy = orig_ds.xy(row, col)
            tiled_xy = tiled_ds.xy(row, col)
            assert orig_xy == tiled_xy
            tile = find_sub_raster(orig_xy, tiled_sub_rasters)
            if not tile:
                print('cannot find tile')
                sys.exit(0)
            #if '-004-' in os.path.basename(tile.fn):
            #    continue
            tile_row, tile_col = rasterio.transform.rowcol(tile.transform,
                                                           *orig_xy)
            tile_xy = rasterio.transform.xy(tile.transform,
                                            tile_row, tile_col)
            assert np.allclose(tile_xy, orig_xy)
            tile_win = rasterio.windows.Window(tile_col, tile_row, 1, 1)
            tile_pixel = tile.ds.read(1, window=tile_win)[0][0]
            # if orig_pixel != tiled_pixel:
                # sub = find_sub_raster(orig_xy, orig_sub_rasters)
                #sub_row, sub_col = rasterio.transform.rowcol(
                #    sub.transform, *orig_xy)
                #sub_win = rasterio.windows.Window(sub_col, sub_row, 1, 1)
                #assert sub_row < sub.ds.shape[0]
                #assert sub_col < sub.ds.shape[1]
                #sub_pixel = sub.ds.read(1, window=sub_win)[0][0]
                #print(os.path.basename(tile.fn), iter.multi_index, (row, col),
                #      orig_xy, orig_pixel, tiled_pixel, tile_pixel, sub_pixel)
                #for dy in (-1, 0, 1):
                #    for dx in (-1, 0, 1):
                #        pixel = tile.ds.read(
                #            1, window=rasterio.windows.Window(
                #                tile_win.col_off + dx, tile_win.row_off + dy,
                #                1, 1))[0][0]
                #        print(dx, dy, pixel)
                #sys.exit(1)
#for ji, window in window_list1:
#    data = src.read(BAND, window=window)
