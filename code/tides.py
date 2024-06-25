#!/opt/python-3.11.5/bin/python3

"""Read the NetCDF file with Tide Gauge data and predictions; parse
into data structures. Print some of them.

"""

import dataclasses
import os
import re
import sys

import fiona
import netCDF4
import numpy as np
import shapely

import american_states as america
import regionlib

from ucs_constants import *

@dataclasses.dataclass
class TideStation:
    num: int
    name: str
    noaa_id: str
    lat: float
    lon: float

    def __post_init__(self):
        # TODO(vitale): reconsider this validation.
        if self.name in ('Apra Harbor, Guam',
                         'Wake Island',
                         'Kwajalein, RMI',
                         'Sand Island (Midway)'):
            self.state = 'territory'
            return
        m = re.match(r'(.+), (?P<state>[A-Z]{2})', self.name)
        if not m:
            raise ValueError(f'tide station name "{self.name}" must end in '
                             '", XX" where XX is a 2-letter U.S. state')
        state = m.group('state')
        if state not in america.state_dict:
            raise ValueError(f'tide station "{self.name}" unknown state '
                             '"{state}"')
        self.state = state

    def coords(self):
        return (self.lon, self.lat)
    
def normalize_longitude(x):
    if x > 180:
        x -= 360
    return x
    
class TideData:
    def __init__(self, netcdf_filename: str, level_var_name: str):
        self.netcdf_filename = netcdf_filename
        self.dataset = ds = netCDF4.Dataset(netcdf_filename)
        # .stations is keyed by TideStation.name.
        self.stations: dict[str, TideStation] = self.read_stations(ds)
        self.severities: dict[str, int] = self.read_severities(ds)
        # .years maps year (e.g. 2025) to record number (e.g. 3).
        self.years: dict[int, int] = self.read_years(ds)
        self.level = ds.variables[level_var_name]
        
    def read_stations(self, ds: netCDF4.Dataset):
        """Reads NetCDF dataset, builds map from name to Tide Station."""
        v_station = ds.variables['station']
        v_lat = ds.variables['latitude']
        v_lon = ds.variables['longitude']
        v_name = ds.variables['station_name']
        v_noaa_id = ds.variables['station_id']        

        stations = {}
        
        for i, station_num in enumerate(v_station):
            assert(i == station_num)
            lat = np.float64(v_lat[i])
            lon = normalize_longitude(np.float64(v_lon[i]))
            noaa_id = v_noaa_id[i]
            name = v_name[i]     
            stations[name] = TideStation(i, name, noaa_id, lat, lon)
            
        return stations

    def read_severities(self, ds: netCDF4.Dataset):
        """Builds dictionary mapping severity name to small integer."""
        severities = {}
        # We obnoxiously rename the NetCDF variable 'scenario' to 'severity',
        # as in "Severity of Sea Level Rise projection". We later use
        # 'scenario' to refer to the (year, severity, flood_frequency) tuple.
        for i, severity in enumerate(ds.variables['scenario']):
            severities[severity] = i

        # Our enum regionlib.Severity exactly models the NetCDF variable.
        assert list(severities.keys()) == regionlib.Severity.all_string()
        return severities

    def read_years(self, ds: netCDF4.Dataset):
        """Reads NetCDF dataset, builds map from year to position in data."""
        years = {}
        for i, year in enumerate(ds.variables['year']):
            years[year[()]] = i  # "[()]" gets scalar from numpy 0-dim array
        return years

    def get_level(self, station_name, severity_name, year) -> int:
        """Returns:
             highest level exceeded at least freq times per year on average,
             for each year, in centimeters above MHHW."""

        station_num = self.stations[station_name].num
        severity_num = self.severities[severity_name]
        year_num = self.years[year]
        return self.level[station_num, severity_num, year_num]

def write_stations(tides: TideData):
    WGS84 = fiona.crs.CRS.from_epsg(4326)
    schema = {
        GEOMETRY_KEY: 'Point',
        PROPERTIES_KEY: {
            'name': 'str',
            'noaa_id': 'str',
            'state': 'str'
        }
    }
    f = fiona.open('/tmp/stations.gpkg', "w", crs=WGS84, schema=schema)
    for station in tides.stations.values():
        point = shapely.geometry.Point(station.lon, station.lat)
        feature = {
            GEOMETRY_KEY: shapely.geometry.mapping(point),
            PROPERTIES_KEY: {
                'name': station.name,
                'noaa_id': station.noaa_id,
                'state': station.state
            }
        }
        f.write(feature)

def write_severities(tides: TideData):
    for severity in tides.severities:
        print(severity)

def print_all_tides(tides: TideData):
    for severity in tides.severities:
        for year in tides.years:
            for station in tides.stations.values():
                v = tides.get_level(station.name, severity, year)
                print(f'{year} {severity:5s} {station.name:20s} {v}')

def print_select_tides(tides: TideData, station_names: list[str]):
    for severity in tides.severities:
        for year in tides.years:
            for name in station_names:
                v = tides.get_level(name, severity, year)
                print(f'{year} {severity:5s} {name:20s} {v}')

def get_tide_data(frequency: int) -> TideData:
    fn = regionlib.tide_data_path(frequency)
    level_variable_name = f'level_{frequency:02d}'
    return TideData(fn, level_variable_name)

# TODO: write some tests! This main() program exercises things a bit.
def main():
    tides = get_tide_data(4)
    # print_all_tides(tides)
    # print_select_tides(tides, ['San Juan, PR', 'Charlotte Amalie, VI'])
    # write_stations(tides)
    write_severities(tides)
    
if __name__ == "__main__":
    main()
