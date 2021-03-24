""" 
Simple earthquake location inversion tools.

Author: Calum Chamberlain
"""

import numpy as np
import math
from datetime import datetime as dt
import logging

LOGGER = logging.getLogger("Geiger")

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


def update_model(p_times, p_locations, s_times, s_locations, model, vp, vs,
                 norm: int = 2):
    """
    Recompute model to minimise residuals.

    :type p_times: np.ndarray
    :param p_times: Observed P arrival times
    :type p_locations: list of dict with keys "x", "y", "z"
    :param p_locations: The locations for the p_times observations
    :type s_times: np.ndarray
    :param s_times: Observed S arrival times
    :type s_locations: list of dict with keys "x", "y", "z"
    :param s_locations: The locations for the s_times observations
    :type model: dict keyed by "x", "y", "z", "time"
    :param model: The current model to be updated
    :type vp: float
    :param vp: P-wave velocity
    :type vs: float
    :param vs: S-wave velocity

    :returns: model
    """
    data = np.concatenate([p_times, s_times])
    residuals = np.zeros_like(data)
    condition = np.zeros((len(residuals), 4))  # G matrix in Stein and Wysession

    for i, p_loc in enumerate(p_locations):
        residuals[i], distance = _residual(model, p_loc, p_times[i], vp, True)
        condition[i][0] = 1
        condition[i][1] = (model["x"] - p_loc["x"]) / (vp * distance)
        condition[i][2] = (model["y"] - p_loc["y"]) / (vp * distance)
        condition[i][3] = (model["z"] - p_loc["z"]) / (vp * distance)

    for i, s_loc in enumerate(s_locations):
        j = i + len(p_times)
        residuals[j], distance = _residual(model, s_loc, s_times[i], vs, True)
        condition[j][0] = 1
        condition[j][1] = (model["x"] - s_loc["x"]) / (vs * distance)
        condition[j][2] = (model["y"] - s_loc["y"]) / (vs * distance)
        condition[j][3] = (model["z"] - s_loc["z"]) / (vs * distance)

    model_change = np.dot(
        np.linalg.inv(np.dot(np.transpose(condition), condition)),
        np.dot(np.transpose(condition), residuals))
    for i, key in enumerate(["time", "x", "y", "z"]):
        model[key] += model_change[i]
    LOGGER.info(
        "Model updated to: "
        "(time {0:.2f}, x {1:.2f}, y {2:.2f} z {3:.2f})".format(
        model["time"], model["x"], model["y"], model["z"]))
    return model, np.linalg.norm(model_change, norm)


def _residual(model, obs_loc, obs_time, velocity, return_distance=False):
    distance = (
            (model["x"] - obs_loc["x"]) ** 2 + 
            (model["y"] - obs_loc["y"]) ** 2 +
            (model["z"] - obs_loc["z"]) ** 2) ** 0.5
    time = model["time"] + distance / velocity
    residual = obs_time - time
    LOGGER.debug(
        "Distance to station: {0:6.2f}\tTime: {1:4.2f}\tVelocity: {2:.2f}\tResidual: {3:.2f}".format(
            distance, time, velocity, residual))
    if return_distance:
        return residual, distance
    return residual


def geiger_locate_lat_lon(p_times, p_locations, s_times, s_locations, vp, vs,
                          max_it=10, convergence=0.1, starting_depth=5.0,
                          plot=False, l2norm: bool = True):
    """
    Locate a seismic source using a simple linearised inversion.

    :type p_times: np.ndarray
    :param p_times: Observed P arrival times
    :type p_locations: list of dict with keys "lat", "lon", "z"
    :param p_locations: The locations for the p_times observations, z should be in km
    :type s_times: np.ndarray
    :param s_times: Observed S arrival times
    :type s_locations: list of dict with keys "lat", "lon", "z"
    :param s_locations: The locations for the s_times observations z should be in km
    :type vp: float
    :param vp: P-wave velocity
    :type vs: float
    :param vs: S-wave velocity
    :type max_it: int
    :param max_int: Maximum iterations for location
    :type convergence: float
    :param convergence: Model change to define convergence
    :type starting_depth: float
    :param starting_depth: Depth for starting location
    :type plot: bool
    :param plot: Whether or not to plot the location and residual history

    :returns dict
    """
    p_xyz = [Geographic(latitude=p["lat"], longitude=p["lon"], depth=p["z"])
             for p in p_locations]
    s_xyz = [Geographic(latitude=p["lat"], longitude=p["lon"], depth=p["z"])
             for p in s_locations]
    origin = Geographic(latitude=p_xyz[0].latitude,
                        longitude=p_xyz[0].longitude,
                        depth=0.0)

    p_xyz = [p.to_xyz(origin, 0, 90) for p in p_xyz]
    s_xyz = [p.to_xyz(origin, 0, 90) for p in s_xyz]
    p_xyz = [{"x": p.x, "y": p.y, "z": -p.z} for p in p_xyz]
    s_xyz = [{"x": p.x, "y": p.y, "z": -p.z} for p in s_xyz]

    # Convert from UTCDateTime to float
    _time = min(np.concatenate([p_times, s_times]))
    _p_times = np.array([_t - _time for _t in p_times])
    _s_times = np.array([_t - _time for _t in s_times])

    model = geiger_locate(
        p_times=_p_times, p_locations=p_xyz, s_times=_s_times,
        s_locations=s_xyz, vp=vp, vs=vs, max_it=max_it, convergence=convergence,
        starting_depth=starting_depth, plot=plot, l2norm=l2norm)
    
    # convert location back to lat lon
    model_time = _time + model["time"]
    model = Location(x=model["x"], y=model["y"], z=model["z"], origin=origin,
                     strike=0, dip=90)
    model = model.to_geographic()
    return {"lat": model.latitude, "lon": model.longitude, "z": model.depth,
            "time": model_time}


def geiger_locate(p_times, p_locations, s_times, s_locations, vp, vs,
                  max_it=10, convergence=0.1, starting_depth=5.0, plot=False,
                  l2norm: bool = True):
    """
    Locate a seismic source using a simple linearised inversion.

    :type p_times: np.ndarray
    :param p_times: Observed P arrival times
    :type p_locations: list of dict with keys "x", "y", "z"
    :param p_locations: The locations for the p_times observations
    :type s_times: np.ndarray
    :param s_times: Observed S arrival times
    :type s_locations: list of dict with keys "x", "y", "z"
    :param s_locations: The locations for the s_times observations
    :type vp: float
    :param vp: P-wave velocity
    :type vs: float
    :param vs: S-wave velocity
    :type max_it: int
    :param max_int: Maximum iterations for location
    :type convergence: float
    :param convergence: Model change to define convergence
    :type starting_depth: float
    :param starting_depth: Depth for starting location
    :type plot: bool
    :param plot: Whether or not to plot the location and residual history

    :returns dict
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from copy import deepcopy

    if l2norm:
        norm = 2
    else:
        norm = 1

    def norm_residual(model):
        r = []
        for i, p_loc in enumerate(p_locations):
            r.append(_residual(
                model=model, obs_loc=p_loc, obs_time=p_times[i], velocity=vp))
        for i, s_loc in enumerate(s_locations):
            r.append(_residual(
                model=model, obs_loc=s_loc, obs_time=s_times[i], velocity=vs))
        return np.linalg.norm(r, norm)

    model = deepcopy(p_locations[np.where(p_times == min(p_times))[0][0]])
    model.update({"z": starting_depth, "time": min(p_times)})
    LOGGER.debug(
        "Starting model at "
        "(time {0:.2f}, x {1:.2f}, y {2:.2f} z {3:.2f})".format(
        model["time"], model["x"], model["y"], model["z"]))
    
    models = [deepcopy(model)]
    residuals = [norm_residual(model)]
    LOGGER.info(
        "Starting model has residual of {0:.2f}".format(residuals[-1]))

    for i in range(max_it):
        model, model_change = update_model(
            p_times=p_times, p_locations=p_locations, s_times=s_times,
            s_locations=s_locations, model=model, vp=vp, vs=vs, norm=norm)
        if model["z"] < 0:
            print("AIRQUAKE! - flipping sign")
            model.update({"z": -1 * model["z"]})
        models.append(deepcopy(model))
        residuals.append(norm_residual(model))
        LOGGER.info(
            "Updated model has residual of {0:.2f}".format(residuals[-1]))
        if model_change < convergence:
            LOGGER.info("Model has reached convergence threshold")
            break

    if plot:
        fig = plt.figure(figsize=(5, 10))
        scatter_ax = fig.add_subplot(2, 1, 1, projection='3d')
        x = [m["x"] for m in models]
        y = [m["y"] for m in models]
        z = [m["z"] for m in models]
        t = [m["time"] for m in models]
        iteration = np.arange(len(models))
        scatter_ax.plot(x, y, z, "k", linewidth=0.5)
        cmplasma = plt.get_cmap("jet")
        scatter_ax.scatter(x, y, z, c=iteration, cmap=cmplasma, label="locations")
        scatter_ax.set_xlabel("x")
        scatter_ax.set_ylabel("y")
        scatter_ax.set_zlabel("z")

        # Plot stations
        sta_x = [s["x"] for s in p_locations + s_locations]
        sta_y = [s["y"] for s in p_locations + s_locations]
        sta_z = [s["z"] for s in p_locations + s_locations]
        scatter_ax.scatter(sta_x, sta_y, sta_z, color="r", marker="v",
                           label="stations")
        scatter_ax.legend()

        scatter_ax.invert_zaxis()

        residual_ax = fig.add_subplot(2, 1, 2)
        residual_ax.plot(residuals, "k")
        residual_ax.set_xlabel("Iteration")
        residual_ax.set_ylabel("Residual (s)")
        plt.show()
    return model