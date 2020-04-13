"""
Code to extract a subset of a catalog adjacent to an arbitrary cross-section.

:author: Calum Chamberlain
:date: 23/03/2018
"""

import math
from datetime import datetime as dt

EARTHRADIUS = 6371  # Global definition of earth radius in km.


class Location():
    """
    Location in x, y, z in some co-ordinate system.
    """
    def __init__(self, x, y, z, origin, strike, dip, 
                 time=None, magnitude=None):
        self.x = x
        self.y = y
        self.z = z
        self.origin = origin
        self.strike = strike
        self.dip = dip
        if time and not isinstance(time, dt):
            raise IOError("time is not a datetime object.")
        self.time = time
        self.magnitude = magnitude

    def __str__(self):
        return ("Location: x: {0} y: {1} z: {2}\n\tOrigin: "
                "{3} Strike: {4} Dip: {5}".format(
                    self.x, self.y, self.z, self.origin, self.strike, 
                    self.dip))

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        for key in self.__dict__.keys():
            if not self.__dict__[key] == other.__dict__[key]:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other) 

    def to_geographic(self):
        """
        Convert to Geographic reference frame.

        :returns: Geographic
        """
        # Rotate back by dip
        s = math.radians(self.strike)
        d = math.radians(90 - self.dip)
        x1 = (self.x * math.cos(-d)) + (self.z * math.sin(-d))
        z = (-self.x * math.sin(-d)) + (self.z * math.cos(-d))
        # Rotate back by strike into North, East, Down reference
        x = (x1 * math.cos(s)) + (self.y * math.sin(s))
        y = (-x1 * math.sin(s)) + (self.y * math.cos(s)) 
        # Convert to geographic co-ordinates
        latitude = y / EARTHRADIUS # Degrees north of origin in radians
        latitude += math.radians(self.origin.latitude)
        latitude = math.degrees(latitude)
        mean_lat = math.radians((latitude + self.origin.latitude) / 2) 
        # is this calculation better done in radians?
        longitude = x / (math.cos(mean_lat) * EARTHRADIUS)
        longitude += math.radians(self.origin.longitude)
        longitude = math.degrees(longitude)
        depth = z + self.origin.depth
        geog = Geographic(latitude=round(latitude, 6),
                          longitude=round(longitude, 6),
                          depth=round(depth, 6), time=self.time,
                          magnitude=self.magnitude)
        return geog


class Geographic():
    """
    Geographic position in lat, long and depth as deg, deg, km (+ve down)
    """
    def __init__(self, latitude, longitude, depth, time=None, magnitude=None):
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        if time and not isinstance(time, dt):
            raise IOError("time is not a datetime object.")
        self.time = time
        self.magnitude = magnitude

    def __str__(self):
        return "Geographic: Lat {0}, Long {1}, Depth (km) {2}".format(
            self.latitude, self.longitude, self.depth)
    
    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        for key in self.__dict__.keys():
            if not self.__dict__[key] == other.__dict__[key]:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other) 

    def to_xyz(self, origin, strike, dip):
        """
        Convert geographic location to arbitrary co-ordinate system.

        :type origin: Location
        :param origin: Origin for new co-ordinate system.
        :type strike: float
        :param strike: Degrees clockwise from north to rotate system.
        :type dip: float
        :param dip: Degrees down from horizontal to rotate system.

        :returns: Location
        """
        # Calculate x, y, z from North, East, Down co-ordinate system
        # x is km East, y is km North, z is km up
        z = self.depth - origin.depth
        mean_lat = math.radians((self.latitude + origin.latitude) / 2)
        x = math.radians(self.longitude - origin.longitude) # Degrees east
        x *= math.cos(mean_lat) * EARTHRADIUS
        y = math.radians(self.latitude - origin.latitude) # Degrees north
        y *= EARTHRADIUS

        s = math.radians(strike)
        d = math.radians(90 - dip)
        # Rotate through strike (clockwise from North)
        x1 = (x * math.cos(-s)) + (y * math.sin(-s))
        y1 = (-x * math.sin(-s)) + (y * math.cos(-s)) 
        """
        x1 is horizontal distance perpendicular to strike,
        y1 is horizontal distance along strike - this needs no further rotation.
        """
        # Rotate z and y1 through dip.
        x2 = (x1 * math.cos(d)) + (z * math.sin(d))
        z1 = (-x1 * math.sin(d)) + (z * math.cos(d))
        return Location(round(x2, 6), round(y1, 6), round(z1, 6),
                        origin, strike, dip, time=self.time, 
                        magnitude=self.magnitude)


if __name__ == '__main__':
    print("Nothing to see here...")