"""
Helpers for getting data programatically.

"""

import datetime
import requests
import numpy as np


def get_geonet_quakes(
    min_latitude=-49.0, max_latitude=-40.0,
    min_longitude=164.0, max_longitude=182.0,
    min_magnitude=0.0, max_magnitude=9.0,
    min_depth=0.0, max_depth=500.0,
    start_time=datetime.datetime(1960, 1, 1),
    end_time=datetime.datetime(2020, 1, 1),
):
    """
    Get a dataframe of the earthquakes in the GeoNet catalogue.
    
    Parameters
    ----------
    min_latitude
        Minimum latitude in degrees for search
    max_latitude
        Maximum latitude in degrees for search
    min_longitude
        Minimum longitude in degrees for search
    max_longitude
        Maximum longitude in degrees for search
    min_depth
        Minimum depth in km for search
    max_depth
        Maximum depth in km for search
    min_magnitude
        Minimum magnitude for search
    max_magnitude
        Maximum magnitude for search
    start_time
        Start date and time for search
    end_time
        End date and time for search
        
    Returns
    -------
    dictionary of columns
    """
    type_mapper = {
        'publicid': str,
         'eventtype': str,
         'longitude': float,
         'latitude': float,
         'magnitude': float,
         'depth': float,
         'magnitudetype': str,
         'depthtype': str,
         'evaluationmethod': str,
         'evaluationstatus': str,
         'evaluationmode': str,
         'earthmodel': str,
         'usedphasecount': int,
         'usedstationcount': int,
         'magnitudestationcount': int,
         'minimumdistance': float,
         'azimuthalgap': float,
         'originerror': float,
         'magnitudeuncertainty': float}
    
    # Convert start_time and end_time to strings
    start_time = start_time.strftime("%Y-%m-%dT%H:%M:%S")
    end_time = end_time.strftime("%Y-%m-%dT%H:%M:%S")
    # Use the more efficient f-string formatting
    query_string = (
        "https://quakesearch.geonet.org.nz/csv?bbox="
        f"{min_longitude},{min_latitude},{max_longitude},"
        f"{max_latitude}&minmag={min_magnitude}"
        f"&maxmag={max_magnitude}&mindepth={min_depth}"
        f"&maxdepth={max_depth}&startdate={start_time}"
        f"&enddate={end_time}")
    print(f"Using query: {query_string}")
    response = requests.get(query_string)
    lines = response.content.decode().splitlines()
    columns = lines[0].split(',')
    result = dict()
    for column_no, column in enumerate(columns):
        column = column.strip()
        if column.endswith("time"):
            vals = np.array([datetime.datetime.strptime(
                row.split(',')[column_no],
                "%Y-%m-%dT%H:%M:%S.%fZ"
            ) for row in lines[1:]])
        else:
            vals = np.array([type_mapper[column](row.split(',')[column_no]) 
                             for row in lines[1:]])
        result.update({column: vals})
    # Ideally this should be ordered by origin time
    order = np.argsort(result["origintime"])
    for key, val in result.items():
        result[key] = val[order]
    return result


def get_gnss_for_station(
    station: str, 
    fits_url: str = "http://fits.geonet.org.nz/observation",
    starttime: datetime.datetime = None,
    endtime: datetime.datetime = None,
) -> dict:
    """
    Get GNSS data from GeoNet for the station
    
    Parameters
    ----------
    station
        The name of the station you want to get data for
    fits_url
        URL of the FITS data service you want to query.
    starttime
        Earliest timestamp to return
    endtime
        Latest timestamp to return
        
    Returns
    -------
    Dictionary with keys:
        time 
            list of timestamps of observations
        north
            list of offsets in mm in the north direction
        east
            list of offsets in mm in the east direction
        up          
            list of vertical offsets in mm
        north_error
            list of errors in mm for north
        east_error
            list of errors in mm for east
        up_error
            list of erros in mm for up
    
    """
    # Initialise an empty dictionary that we will append to
    out = dict(time=[],
               north=[],
               east=[],
               up=[],
               north_error=[],
               east_error=[],
               up_error=[])
    for channel in {"north", "east", "up"}:
        parameters = {"typeID": channel[0], "siteID": station}
        response = requests.get(fits_url, params=parameters)
        assert response.status_code == 200, "Bad request"
        payload = response.content.decode("utf-8").split("\n")
        # payload is a csv with header
        # This is a list-comprehension, a type of fast, one-line for loop
        payload = [p.split(',') for p in payload]
        # Check that this is what we expect
        assert payload[0][0] == 'date-time', "Unkown format"
        assert len(payload[0]) == 3, "Unknown format"
        times, displacements, errors = zip(*[
            (datetime.datetime.strptime(p[0], '%Y-%m-%dT%H:%M:%S.%fZ'),
             float(p[1]), float(p[2])) for p in payload[1:-1]])
        if len(out["time"]) == 0:
            out.update({"time": times})
        else:
            assert out["time"] == times, "Different time sampling for different components."
        out.update({channel: displacements, f"{channel}_error": errors})
    starttime = starttime or min(out["time"])
    endtime = endtime or max(out["time"])
    mask = [starttime <= t <=endtime for t in out["time"]]
    for key, value in out.items():
        out.update({key: [v for v, m in zip(value, mask) if m]})
    return out