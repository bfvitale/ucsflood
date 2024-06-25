#!/opt/python-3.11.5/bin/python3

import logging
import os
import pathlib
import sys
from typing import Iterable

import numpy as np
import pathlib
import rasterio
import shapely
import xml.etree.ElementTree as xet

import regionlib
from ucs_constants import *

def sub_rasters_from_vrt(vrt_path: pathlib.Path) -> list[pathlib.Path]:
    # TODO: instead use gdal.Open().GetFileList(). But, that includes the VRT
    # itself.
    tree = xet.parse(vrt_path)
    paths = []
    for elem in tree.getroot().iter('SourceFilename'):
        path = pathlib.Path(elem.text)
        if elem.get('relativeToVRT'):
            path = vrt_path.resolve(strict=True).parent / path
        paths.append(path)
    return paths
        
def get_vrt_sub_extent(vrt_path: pathlib.Path) -> shapely.MultiPolygon:
    """Returns a MultiPolygon which is the union of the bounding boxes of each
    sub-raster in the input VRT.

    Pre-condition: The VRT is simply a mosaic with no other processing."""

    boxes = []
    for raster_path in sub_rasters_from_vrt(vrt_path):
        ds = rasterio.open(raster_path)
        boxes.append(shapely.geometry.box(*ds.bounds))
    r = shapely.union_all(boxes)
    if r.geom_type == 'Polygon':
        r = shapely.MultiPolygon([r])
    return r

def tiles_in_window(src_window: rasterio.windows.Window, tile_size: int
                    ) -> Iterable[rasterio.windows.Window]:
    """A generator yielding a set of tiles, each a square Window of size
    (tile_size, tile_size), which cover the src_window."""
    windows = []
    for row in range(0, src_window.height, tile_size):
        tile_height = tile_size
        if row + tile_height > src_window.height:
            tile_height = src_window.height - row
        for col in range(0, src_window.width, tile_size):
            tile_width = tile_size
            if col + tile_width > src_window.width:
                tile_width = src_window.width - col
            window = rasterio.windows.Window(
                src_window.col_off + col, src_window.row_off + row,
                tile_width, tile_height)
            windows.append(window)
    return windows

def blocks_in_tile(tile: rasterio.windows.Window, block_size: int
                   ) -> Iterable[rasterio.windows.Window]:
    """An alias for def 'tiles_in_window', for clarity of exposition."""
    return tiles_in_window(tile, block_size)

def nonempty_tiles_in_dem_vrt(region: regionlib.Region, tile_size: int
                              ) -> Iterable[rasterio.windows.Window]:
    """Returns List of Tiles covering VRT for region, excluding any that are
    empty, i.e. entirely NO_DATA, i.e. aren't backed by any sub-raster in VRT.

    Note, a returned tile might still be empty, because it intersects with an
    empty part of a sub-raster.
    """
    vrt_path = regionlib.regional_dem_vrt_path(region)
    vrt_ds = rasterio.open(vrt_path)
    vrt_window = rasterio.windows.Window(0, 0, vrt_ds.width, vrt_ds.height)
    all_tiles = tiles_in_window(vrt_window, tile_size)

    # TODO: Why do this in coordinate space? Can we stay in pixel space by
    # consulting the DstRect each SourceFilename sub raster in the VRT XML?
    # But, violating encapsulation of VRT XML is not ideal.
    boxes = [
        shapely.box(*rasterio.windows.bounds(tile, vrt_ds.transform))
        for tile in all_tiles
    ]
    
    DEBUG = True
    if DEBUG:
        # Check that tiles exactly cover VRT.
        assert shapely.equals(shapely.union_all(boxes),
                              shapely.box(*vrt_ds.bounds))

    sub_extent = get_vrt_sub_extent(vrt_path)
    
    if DEBUG:
        # Ensure area is mostly contiguous. WEST has gap between coast and
        # Channel Islands. PRVI has gap betwen PR and VI.
        if region in [regionlib.Region.WEST, regionlib.Region.PRVI]:
            assert len(sub_extent.geoms) == 2
        else:
            assert len(sub_extent.geoms) == 1

    non_empty_tiles = [tile
                       for tile, box in zip(all_tiles, boxes)
                       if box.intersects(sub_extent)]
    logging.info(
        f'{len(non_empty_tiles)} non-empty of total {len(all_tiles)} tiles')
    return non_empty_tiles

class NoDataBlockCache:
    def __init__(self):
        self._blocks = {}
        
    def get(self, shape: tuple[int, int], val: np.float32 | np.uint8):
        key = (val, shape)
        if key not in self._blocks:
            self._blocks[key] = np.full(shape, val)
        return self._blocks[key]

def setup_log(prog_name, is_child=False):
    handler_stderr = logging.StreamHandler(sys.stderr)
    if is_child:
        file_pid = os.getppid()
    else:
        file_pid = os.getpid()
    filename = pathlib.Path(f'/var/tmp/log/{prog_name}.{file_pid}')
    log_dir = filename.parent
    if not log_dir.exists():
        log_dir.mkdir()
    log_file = open(filename, 'a')
    handler_stderr = logging.StreamHandler(sys.stderr)
    handler_file = logging.StreamHandler(log_file)

    formatter = logging.Formatter('%(asctime)s %(process)d %(message)s',
                                  datefmt='%F %T')
    handler_stderr.setFormatter(formatter)
    handler_file.setFormatter(formatter)
    
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.addHandler(handler_stderr)
    root.addHandler(handler_file)
    if not is_child:
        sys.stderr.write(f'logging to {filename}\n')
        logging.info(f'starting')
