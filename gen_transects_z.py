#!/opt/python-3.11.5/bin/python3

"""
  Add tide heights to transects-with-stations file, given year and Sea Level
  Rise severity. Output to file transects-points-z-{year}-{severity}.

  TODO: fold this simple fast step into interpolate_slr.py.
"""

# Sample invocation
# gen_transects_z.py --year=2040 --severity=high --frequency=26

import argparse
import logging
import sys

import fiona
import numpy
import shapely

import regionlib
import tides
from ucs_constants import *

def setup_opts():
    parser = argparse.ArgumentParser(description=(
        'Add tide heights to transects-with-stations file, given year and '
        'Sea Level Rise severity.'
    ))
    parser.add_argument(
        '--year', type=int,
        help='Sea Level projection year [2030..2100]. Default: each decade'
    )
    parser.add_argument(
        '--severity', type=str,
        choices=[x.lower() for x in regionlib.Severity._member_names_],
        help='Sea Level Rise severity to process. By default, all.'
    )
    parser.add_argument(
        '--frequency', type=int, help='Number of floods per year',
        choices=regionlib.available_flood_frequencies()
    )

    args = parser.parse_args()
    if args.year:
        if args.year not in range(2020, 2100 + 1):
            sys.stderr.write('if set, --year must be in [2020..2100]\n')
            sys.exit(1)

    if (args.year or args.severity or args.frequency) and not (
            args.year and args.severity and args.frequency):
        sys.stderr.write('if any of --year, --severity, or --frequency are set,'
                         ' all must be set\n')
        sys.exit(1)
        
    return args

def process_one(tide_data: tides.TideData, exparms: regionlib.ExParms):
    transect_path = TRANSECT_POINT_FILENAME
    out_path = regionlib.slr_transect_path(exparms)
    
    in_fd = fiona.open(transect_path)
    out_schema = in_fd.schema.copy()
    out_schema[PROPERTIES_KEY][TIDE_HEIGHT_KEY] = 'float'
    out_fd = fiona.open(out_path, "w", schema=out_schema, crs=in_fd.crs)

    out_features = []
    for feature in in_fd:
        point = shapely.geometry.shape(feature.geometry)
        station = feature.properties[TIDE_STATION_NAME_KEY]
        z = tide_data.get_level(station, exparms.severity, exparms.year)
        assert numpy.ndim(z) == 0
        z = float(z)

        # TideData units are int centimeters above MHHW; convert to float
        # meters.
        z /= 100.
        
        # We'd like to store the Z value as the z-coord of a 3D Point. Probably
        # in meters of elevation relative to some reference. But, it's unclear
        # how, when X and Y are in lon / lat units, rather than meters.
        # zpoint = shapely.geometry.Point(point.x, point.y, z).

        # Store Z as separate property; ArcGIS Pro seems to prefer this.
        out_properties = dict(feature.properties)
        out_properties[TIDE_HEIGHT_KEY] = z

        out_features.append(regionlib.shapely_to_feature(point, out_properties))

    out_fd.writerecords(out_features)
    print(f'wrote {out_path.stem}')
    return 0

def main():
    opts = setup_opts()
    if opts.year:
        years = [opts.year]
    else:
        years = range(2030, 2100 + 1, 10)

    if opts.severity:
        severities = [regionlib.Severity(opts.severity)]
    else:
        severities = list(regionlib.Severity)

    if opts.frequency:
        frequencies = [opts.frequency]
    else:
        frequencies = regionlib.available_flood_frequencies()

    for frequency in frequencies:
        tide_data = tides.get_tide_data(frequency)
        for year in years:
            for severity in severities:
                exparms = regionlib.ExParms(year, severity, frequency)
                process_one(tide_data, exparms)

if __name__ == '__main__':
    sys.exit(main())
