#!/opt/python-3.11.5/bin/python3

"""
  Joins transects to their nearest tide gauge station.
  That is, for each transect, find the nearest tide gauge station.
  Transects were hand-drawn by humans.

  The output is a .gpkg file (similar to a "shape" or "vector" file).

  See def main(), below, for the schema. Roughly, each record (shape)
  in the file is a 'LineString' geometry (one transect), annotated
  with 'tide_station_name'.

  GIS community refers to this operation as a "Spatial Join".

"""

import math
from math import acos, asin, atan2, cos, degrees, pi, radians, sin, sqrt
import sys
from typing import Dict, List, Sequence, Tuple
import unittest

import fiona
import numpy as np
from numpy.testing import assert_allclose
import shapely
import shapely.ops
from sklearn.neighbors import BallTree

import regionlib
import tides
from ucs_constants import *

WGS84 = fiona.crs.from_epsg(4326)

class Point:
    R = 6371000  # mean earth radius in meters
    
    def __init__(self, point: fiona.model.Feature):
        assert point.geometry['type'] == 'Point'
        self._fiona_point = point
        
    @staticmethod
    def from_lon_lat(lon, lat, properties={}):
        return Point(
            regionlib.shapely_to_feature(shapely.Point(lon, lat), properties))

    @property
    def coords(self):
        """Returns (x, y) (or (lon, lat)) tuple."""
        return self._fiona_point.geometry.coordinates

    @property
    def lat(self):
        return self.coords[1]

    @property
    def lon(self):
        return self.coords[0]

    def Sdistance(self, other: 'Point'):
        """Returns distance in meters from this Point to other Point, using a
           Spherical Earth.
        """
        lat0 = radians(self.lat)
        lat1 = radians(other.lat)
        dlon = radians(other.lon - self.lon)

        # Haversine distance.
        dlat = lat1 - lat0
        ang = sin(dlat / 2) ** 2 + cos(lat0) * cos(lat1) * sin(dlon / 2) ** 2
        drad = 2 * atan2(sqrt(ang), sqrt(1 - ang))
        return drad * self.R
        
        # Spherical law of cosines.
        #ang = sin(lat0) * sin(lat1) + cos(lat0) * cos(lat1) * cos(dlon)
        #return acos(ang) * self.R

    def Sbearing(self, other: 'Point') -> float:
        """Bearing on Spherical Earth."""
        lat0 = radians(self.lat)
        lat1 = radians(other.lat)
        dlon = radians(other.lon - self.lon)
        dy = sin(dlon) * cos(lat1)
        dx = cos(lat0) * sin(lat1) - sin(lat0) * cos(lat1) * cos(dlon)
        angle = degrees(atan2(dy, dx))
        return (angle + 360) % 360
 
    def __str__(self):
        return f'Point({self.lon:.6f}, {self.lat:.6f})'

def intermediate_point(p0: Point, p1: Point, fraction: float) -> Point:
    """Returns a point that is 'fraction' along the great-circle path from p0
    to p1."""
    lon0, lat0, lon1, lat1 = map(radians, (p0.lon, p0.lat, p1.lon, p1.lat))
    delta = p0.Sdistance(p1) / Point.R
    a = sin((1 - fraction) * delta) / sin(delta)
    b = sin(fraction * delta) / sin(delta)    
    x = a * cos(lat0) * cos(lon0) + b * cos(lat1) * cos(lon1)
    y = a * cos(lat0) * sin(lon0) + b * cos(lat1) * sin(lon1)
    z = a * sin(lat0) + b * sin(lat1)

    lat = atan2(z, sqrt(x ** 2 + y ** 2))
    lon = atan2(y, x)
    return Point.from_lon_lat(degrees(lon), degrees(lat))

def points_from_line(line: shapely.geometry.LineString,
                     num_points: int) -> List[Point]:
    """We don't use shapely.segmentize(), because we are working on sphere
    using unprojected coordinates. (We could project, segmentize,
    unproject.)
    """
    assert len(line.coords) == 2  # only 1 line in LineString
    assert num_points > 1
    lc0, lc1 = line.coords
    p0 = Point.from_lon_lat(*lc0)
    p1 = Point.from_lon_lat(*lc1)
    interval = 1 / (num_points - 1)
    return [intermediate_point(p0, p1, i * interval) for i in range(num_points)]
    
class TestPointsFromLine(unittest.TestCase):
    def test_intermediate_point(self):
        cambridge = Point.from_lon_lat(0.119, 52.205)
        paris = Point.from_lon_lat(2.351, 48.857)
        ip1 = intermediate_point(cambridge, paris, 0.25)
        self.assertAlmostEqual(ip1.lon, 0.7073, 4)
        self.assertAlmostEqual(ip1.lat, 51.3721, 4)

        lax = Point.from_lon_lat(118.4, 33.95)
        jfk = Point.from_lon_lat(73.783, 40.633)
        ip2 = intermediate_point(lax, jfk, 100 / 2144)
        self.assertAlmostEqual(ip2.lon, 116.552, 3)
        self.assertAlmostEqual(ip2.lat, 34.617, 3)

        ip3 = intermediate_point(lax, jfk, 0.4)
        self.assertAlmostEqual(ip3.lon, 101.626, 3)
        self.assertAlmostEqual(ip3.lat, 38.669, 3)
        
    def test_points_from_line(self):
        lax = Point.from_lon_lat(118.4, 33.95)
        jfk = Point.from_lon_lat(73.783, 40.633)
        line = shapely.geometry.LineString([lax.coords, jfk.coords])
        points1 = [p.coords for p in points_from_line(line, 3)]
        assert_allclose(points1, [lax.coords, (97.137, 39.455), jfk.coords],
                        atol=1e-3)
    
def cross_arc_distance(e0: Point, e1: Point, p: Point) -> float:
    """Returns distance in meters between p and great-circle-segment arc
    (e0, e1).
    
    N.B. If cross track point is not on the minor arc between
    endpoints, instead return distance from p to closest endpoint.

    See http://stackoverflow.com/questions/32771458 

    """
    R = Point.R
    dist0p_meters = e0.Sdistance(p)
    dist0p_rad = dist0p_meters / R  # radians
    theta0p = np.deg2rad(e0.Sbearing(p))
    theta01 = np.deg2rad(e0.Sbearing(e1))

    dtheta = abs(theta0p - theta01)
    if pi / 2 < dtheta < 3 * pi / 2:
        # Obtuse angle. Closest point to p is e0.
        return dist0p_meters

    # Cross-track ("XT") distance.
    dist_XT_rad = asin(sin(dist0p_rad) * sin(theta0p - theta01))
    dist_XT_meters = dist_XT_rad * R

    # Compute along-track ("AT") distance from e0 to p.
    dist_AT_meters = acos(cos(dist0p_rad) / cos(dist_XT_rad)) * R

    # Is cross-track intersection point ("3") beyond the arc segment?
    if dist_AT_meters > e0.Sdistance(e1):
        # Yes, point 3 is beyond arc segment. Closest point to p is e1.
        return e1.Sdistance(p)

    # Point 3 is _on_ arc segment, cross-arc distance cross-track distance.
    return abs(dist_XT_meters)

class TestCrossTrack(unittest.TestCase):
    def test_base(self):
        bangor = Point.from_lon_lat(-68.828, 44.807)
        montreal = Point.from_lon_lat(-73.741, 45.471)
        halifax = Point.from_lon_lat(-63.510, 44.880)
        ottawa = Point.from_lon_lat(-75.667, 45.323)
        lawrencetown = Point.from_lon_lat(-63.352, 44.647)
        
        self.assertAlmostEqual(cross_arc_distance(montreal, halifax, bangor),
                               54962.862,  3)
        self.assertAlmostEqual(cross_arc_distance(halifax, montreal, bangor),
                               54962.862,  3)

        # Test cases when cross-track point lies outside shorter arc.

        # Ottawa is west of the (montreal - halifax) arc.
        self.assertAlmostEqual(cross_arc_distance(montreal, halifax, ottawa),
                               ottawa.Sdistance(montreal), 3)
        self.assertAlmostEqual(cross_arc_distance(halifax, montreal, ottawa),
                               ottawa.Sdistance(montreal), 3)

        # Lawrencetown is east of the (montreal - halifax) arc.
        self.assertAlmostEqual(
            cross_arc_distance(halifax, montreal, lawrencetown),
            lawrencetown.Sdistance(halifax), 3)
        self.assertAlmostEqual(
            cross_arc_distance(montreal, halifax, lawrencetown),
            lawrencetown.Sdistance(halifax), 3)

        # Some tests from stackoverflow.
        def PLL(lon, lat): return Point.from_lon_lat(lon, lat)
        
        self.assertAlmostEqual(
            cross_arc_distance(PLL(-55.5, -10.1), PLL(-45.1, -15.2),
                               PLL(-62.5, -10.5)),
            767094.768, 3)
        self.assertAlmostEqual(
            cross_arc_distance(PLL(-45.1, -15.2), PLL(-55.5, -10.1),
                               PLL(-62.5, -10.5)),
            767094.768, 3)

        self.assertAlmostEqual(
            cross_arc_distance(PLL(60.5, 40.5), PLL(80.5, 50.5), PLL(69, 51)),
            479609.299, 3)
        self.assertAlmostEqual(
            cross_arc_distance(PLL(80.5, 50.5), PLL(60.5, 40.5), PLL(69, 51)),
            479609.299, 3)
        
        self.assertAlmostEqual(
            cross_arc_distance(PLL(35.61, 21.72), PLL(40.7, 23.65),
                               PLL(42, 25)),
            199706.839, 3)
        self.assertAlmostEqual(
            cross_arc_distance(PLL(40.7, 23.65), PLL(35.61, 21.72),
                               PLL(42, 25)),
            199706.839, 3)

        self.assertAlmostEqual(cross_arc_distance(PLL(30.98883,-29.7762),
                                                  PLL(30.98877, -29.77483),
                                                  PLL(30.988693, -29.775664)),
                               10.949, 3)
        self.assertAlmostEqual(cross_arc_distance(PLL(30.98877, -29.77483),
                                                  PLL(30.98883,-29.7762),
                                                  PLL(30.988693, -29.775664)),
                               10.949, 3)

# unittest.main(); sys.exit(0)

def cross_arc_shapely(line: shapely.geometry.LineString,
                      point: Point):
    """Cross-arc distance from Point to one-line shapely LineString."""
    assert len(line.coords) == 2  # only 1 line in LineString
    lc0, lc1 = line.coords
    e0 = Point.from_lon_lat(*lc0)
    e1 = Point.from_lon_lat(*lc1)
    return cross_arc_distance(e0, e1, point)

DISTANCE_TO_TIDE_STATION_KEY = 'distance_to_tide_station'
POINT_NUM_KEY = 'point_num'
TIDE_STATION_NAME_KEY = 'tide_station_name'
TRANSECT_ID_KEY = 'transect_id'

def rad_rev(p: Tuple[float,float]) -> Tuple[float,float]:
    """We use (lon, lat) degrees; scikit.haversine uses (lat,lon) radians."""
    return (radians(p[1]), radians(p[0]))

def process_transect_segment(ball: BallTree, stations: List[tides.TideStation],
                             segment_file, feature: fiona.Feature):
    """Annotates transect with nearest TideStation, and writes to segment_file.
       N.B. 'ball' contains Tide Station coords in (lat,lon) radians."""
    
    transect = shapely.geometry.shape(feature.geometry)
   
    # Initial guess: get nearest 10 TideStations to _one endpoint_ of transect.
    indices_nearest = ball.query([rad_rev(transect.coords[0])], k=10)[1][0]
    assert(len(indices_nearest) == 10)

    # Find exact distance from each station to _entire_ transect (not just an
    # endpoint), and choose the closest station.
    closest_dist: float = math.inf
    closest: tides.TideStation = None

    for index in indices_nearest:
        station = stations[index]
        station_point = Point.from_lon_lat(*station.coords())
        distance = cross_arc_shapely(transect, station_point)
        if distance < closest_dist:
            closest_dist = distance
            closest = station

    props = dict(feature.properties.items())
    props[TIDE_STATION_NAME_KEY] = closest.name
    shape = shapely.geometry.shape(feature.geometry)
    region = regionlib.region_from_point(shape.representative_point())
    assert region == regionlib.canonicalize_region(props[regionlib.REGION_KEY])
    props[regionlib.REGION_KEY] = region
    
    segment_file.write(
        fiona.Feature(geometry=feature.geometry, properties=props))

def process_transect_points(ball: BallTree, stations: List[tides.TideStation],
                            point_file, feature: fiona.Feature,
                            transect_id: int) -> int:
    """Factors transect into points, annotates each Point with nearest
       TideStation. Writes each Point to output. Returns number of points."""

    transect = shapely.geometry.shape(feature.geometry)
    num_points = int(feature.properties['num_points'])
    points = points_from_line(transect, num_points)

    # Find closest tide station to each transect point.
    records = []
    for point_num, point in enumerate(points):
        distances, indices = ball.query([rad_rev(point.coords)], k=1)
        # TODO: Consider recomputing distance using elliptical earth.
        distance = distances[0][0] * Point.R  # convert radians to meters
        station = stations[indices[0][0]]

        props = dict(feature.properties.items())
        props.update({
            POINT_NUM_KEY: point_num,
            DISTANCE_TO_TIDE_STATION_KEY: distance,
            TIDE_STATION_NAME_KEY: station.name,
            TRANSECT_ID_KEY: transect_id
        })
        shape = shapely.geometry.Point(point.coords)
        records.append(regionlib.shapely_to_feature(shape,props))
        
    point_file.writerecords(records)
    return len(records)

def main():
    # Read Tide data.
    # We hard-code flood_frequency = 26 here, because we aren't reading the
    # actual tide heights, at this pipeline stage. We only care about the
    # coordinates and ID of each tide-station. However, the tides.get_tide_data
    # function requires a frequency parameter. TODO: read tide data without Z.
    # The actual frequency is used to read the tide data in gen_transects_z.py
    # when it calls TideData.get_level(station, severity, year)
    FLOOD_FREQUENCY = 26
    tide_data = tides.get_tide_data(FLOOD_FREQUENCY)

    # Open transect input file.
    infile = fiona.open(TRANSECT_SHAPEFILE)

    # Build BallTree of Tide Stations.
    stations = list(tide_data.stations.values())
    ball_arg = [rad_rev(s.coords()) for s in stations]
    ball = BallTree(ball_arg, metric='haversine')
    
    # Open two output 'vector' files.
    TRANSECT_SCRATCH_DIR.mkdir(parents=True, exist_ok=True)

    # The first output file models transects as line segments.
    schema = infile.schema
    assert schema[GEOMETRY_KEY] in ('LineString', 'MultiLineString')
    schema[PROPERTIES_KEY][TIDE_STATION_NAME_KEY] = 'str:80'
    segment_file = fiona.open(TRANSECT_SEGMENT_FILENAME, "w", schema=schema,
                              crs=infile.crs)

    # The other output file models each transect as a series of points.
    schema[GEOMETRY_KEY] = 'Point'
    props = schema[PROPERTIES_KEY]
    props[POINT_NUM_KEY] = 'int'
    props[DISTANCE_TO_TIDE_STATION_KEY] = 'float'
    props[TRANSECT_ID_KEY] = 'int'

    point_file = fiona.open(TRANSECT_POINT_FILENAME, "w", schema=schema,
                            crs=infile.crs)

    # Join tides with transects (annotate each transect line and transect point
    # with the ID of its nearest Tide Station).
    # Write output file records as we go.
    total_points = 0
    for num_transects, feature in enumerate(infile):
        transect_id = len(segment_file)
        process_transect_segment(ball, stations, segment_file, feature)
        num_pts = process_transect_points(ball, stations, point_file, feature,
                                          transect_id)
        assert num_pts
        total_points += num_pts
    assert num_transects
    print(f'wrote {total_points} points {TRANSECT_POINT_FILENAME.name}')
    print(f'wrote {num_transects} transects {TRANSECT_SEGMENT_FILENAME.name}')
    return 0

if __name__ == '__main__':
    sys.exit(main())
