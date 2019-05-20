"""
Simple linearised inversion of seismic phase arrival times to get source
location in a homogeneous half-space. Used for teaching, not production ready.

:author: Calum Chamberlain
:date: 13/5/2019
"""

import numpy as np
import logging

LOGGER = logging.getLogger("Geiger")


def seisan_hash(event, inventory, velocities, vpvs, search_angle=2, 
                min_incorrect=1, ratio_error=0.2, angle_prob=10,
                multiple_threshold=0.1, clean=True):
    """
    Use SEISAN's implementation of HASH to compute the focal mechanism of an
    event.

    :note: 
        Requires command line version of driver to be installed. This is
        provided alongside this code as gphs445_hash.for and should be copied
        to your PRO directory as hash_seisan.for - this can be built by running
        `make all` in your PRO directory.
    """
    import subprocess
    from obspy.core.event import (
        QuantityError, NodalPlane, NodalPlanes, FocalMechanism,
        ResourceIdentifier)

    event_back = seisan_hyp(
        event=event, inventory=inventory, velocities=velocities, vpvs=vpvs,
        clean=False)
    subprocess.call(
        ["hash_seisan", str(search_angle), str(min_incorrect),
         str(ratio_error), str(angle_prob), str(multiple_threshold)])
    focal_mechanisms = _read_hash()
    for mechanism in focal_mechanisms:
        nodal_p = NodalPlane(
            strike=mechanism["strike"], dip=mechanism["dip"],
            rake=mechanism["rake"])
        fm = FocalMechanism(nodal_planes=NodalPlanes(nodal_plane_1=nodal_p))
        fm.method_id = ResourceIdentifier(
            "smi:nc.anss.org/focalMehcanism/HASH")
        event_back.focal_mechanisms.append(fm)
    if clean:
        _cleanup()
    return event_back


def _read_hash():
    import os
    if not os.path.isfile("hash_seisan.out"):
        Logger.error("hash_seisan.out not found")
        return []
    with open("hash_seisan.out", "r") as f:
        lines = f.read().splitlines()
    focal_mechanisms, focal_mechanism = ([], None)
    for line in lines:
        if line.startswith("Strike,dip,rake"):
            if focal_mechanism:
                focal_mechanisms.append(focal_mechanism)
            focal_mechanism = {
                "strike": float(line.split()[-3]),
                "dip": float(line.split()[-2]),
                "rake": float(line.split()[-1])}
        elif line.startswith("Fault+aux plane uncertainty"):
            focal_mechanism.update({
                "uncertainty 1": float(line.split()[-2]),
                "uncertainty_2": float(line.split()[-1])})
    if focal_mechanism:
        focal_mechanisms.append(focal_mechanism)
    return focal_mechanisms


def seisan_hyp(event, inventory, velocities, vpvs, clean=True):
    """
    Use SEISAN's Hypocentre program to locate an event.
    """
    import warnings
    import subprocess
    from obspy.core.event import Origin
    from nordic_io.core import write_select, read_nordic
    
    # Write STATION0.HYP file
    _write_station0(inventory, velocities, vpvs)
    
    subprocess.call(['remodl'])
    subprocess.call(['setbrn'])
    
    event_unlocated = event.copy()
    try:
        old_origin = event.preferred_origin() or event.origins[0]
        origin = Origin(time=old_origin.time)
    except IndexError:
        origin = Origin(time=sorted(event.picks, key=lambda p: p.time)[0].time)
    event_unlocated.origins = [origin]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        event_unlocated.write(format="NORDIC", filename="to_be_located")
    subprocess.call(['hyp', "to_be_located"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        event_back = read_nordic("hyp.out")
    event_out = event.copy()
    # We lose some info in the round-trip to nordic
    event_out.origins[0] = event_back[0].origins[0]
    event_out.picks = event_back[0].picks
    if clean:
        _cleanup()
    return event_out

def _cleanup():
    import os

    # Clean up
    files_to_remove = [
        "hyp.out", "to_be_located", "remodl.tbl", "remodl1.lis", "remodl2.lis",
        "print.out", "gmap.cur.kml", "hypmag.out", "hypsum.out", "remodl.hed",
        "IASP91_linux.HED", "IASP91_linux.TBL", "setbrn1.lis", "setbrn2.lis",
        "setbrn3.lis", "STATION0.HYP", "focmec.dat", "focmec.inp", "fort.17",
        "fps.out", "hash_seisan.out", "pspolar.inp"]
    for f in files_to_remove:
        if os.path.isfile(f):
            os.remove(f)


def _stationtoseisan(station):
    """
    Convert obspy inventory to string formatted for Seisan STATION0.HYP file.

    :type station: obspy.core.inventory.station.Station
    :param station: Inventory containing a single station.

    :returns: str

    .. note:: Only works to the low-precision level at the moment (see seisan \
        manual for explanation).
    """

    if station.latitude < 0:
        lat_str = 'S'
    else:
        lat_str = 'N'
    if station.longitude < 0:  # Stored in =/- 180, not 0-360
        lon_str = 'W'
    else:
        lon_str = 'E'
    if len(station.code) > 4:
        sta_str = station.code[0:4]
    else:
        sta_str = station.code.ljust(4)
    if len(station.channels) > 0:
        depth = station.channels[0].depth
    else:
        msg = 'No depth found in station.channels, have you set the level ' +\
              'of stationXML download to channel if using obspy.get_stations?'
        raise IOError(msg)
    elev = str(int(round(station.elevation - depth))).rjust(4)
    # lat and long are written in STATION0.HYP in deg,decimal mins
    lat = abs(station.latitude)
    lat_degree = int(lat)
    lat_decimal_minute = (lat - lat_degree) * 60
    lon = abs(station.longitude)
    lon_degree = int(lon)
    lon_decimal_minute = (lon - lon_degree) * 60
    lat = ''.join([str(int(abs(lat_degree))),
                   '{0:.2f}'.format(lat_decimal_minute).rjust(5)])
    lon = ''.join([str(int(abs(lon_degree))),
                   '{0:.2f}'.format(lon_decimal_minute).rjust(5)])
    station_str = ''.join(['  ', sta_str, lat, lat_str, lon, lon_str, elev])
    return station_str


def _write_station0(inventory, velocities, vpvs):
    out = (
        "RESET TEST(02)=500.0\nRESET TEST(07)=-3.0\nRESET TEST(08)=2.6\n"
        "RESET TEST(09)=0.001\nRESET TEST(11)=99.0\nRESET TEST(13)=5.0\n"
        "RESET TEST(34)=1.5\nRESET TEST(35)=2.5\nRESET TEST(36)=0.0\n"
        "RESET TEST(41)=20000.0\nRESET TEST(43)=5.0\nRESET TEST(51)=3.6\n"
        "RESET TEST(50)=1.0\nRESET TEST(56)= 1.0\nRESET TEST(58)= 99990.0\n"
        "RESET TEST(40)=0.0\nRESET TEST(60)=0.0\nRESET TEST(71)=1.0\n"
        "RESET TEST(75)=1.0\nRESET TEST(76)=0.910\nRESET TEST(77)=0.00087\n"
        "RESET TEST(78)=-1.67\nRESET TEST(79)=1.0\nRESET TEST(80)=3.0\n"
        "RESET TEST(81)=1.0\nRESET TEST(82)=1.0\nRESET TEST(83)=1.0\n"
        "RESET TEST(88)=1.0\nRESET TEST(85)=0.1\nRESET TEST(91)=0.1\n")
    for network in inventory:
        for station in network:
            out += "\n" + _stationtoseisan(station)
    out += "\n\n"
    # Add velocity model
    for layer in velocities:
        if "moho" in layer.keys() and layer["moho"]:
            out += "{0:7.3f}   {1:7.3f}    N     \n".format(
                layer["velocity"], layer["top"])
        else:
            out += "{0:7.3f}   {1:7.3f}\n".format(
                layer["velocity"], layer["top"])
    out += "\n15.0 1100.2200. {0:.2f} \nTES\n".format(vpvs)
    with open("STATION0.HYP", "w") as f:
        f.write(out)
    return


def update_model(p_times, p_locations, s_times, s_locations, model, vp, vs):
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
    return model, np.linalg.norm(model_change, 2)


def _residual(model, obs_loc, obs_time, velocity, return_distance=False):
    distance = (
            (model["x"] - obs_loc["x"]) ** 2 + 
            (model["y"] - obs_loc["y"]) ** 2 +
            (model["z"] - obs_loc["z"]) ** 2) ** 0.5
    time = model["time"] + distance / velocity
    residual = obs_time - time
    LOGGER.debug(
        "Distance to station: {0:.2f}\tResidual: {1:.2f}".format(
            distance, residual))
    if return_distance:
        return residual, distance
    return residual


def geiger_locate_lat_lon(p_times, p_locations, s_times, s_locations, vp, vs,
                          max_it=10, convergence=0.1, starting_depth=5.0,
                          plot=False):
    """
    Locate a seismic source using a simple linearised inversion.

    :type p_times: np.ndarray
    :param p_times: Observed P arrival times
    :type p_locations: list of dict with keys "lat", "lon", "z"
    :param p_locations: The locations for the p_times observations
    :type s_times: np.ndarray
    :param s_times: Observed S arrival times
    :type s_locations: list of dict with keys "lat", "lon", "z"
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
    from coordinates import Geographic, Location

    p_xyz = [Geographic(latitude=p["lat"], longitude=p["lon"], depth=p["z"])
             for p in p_locations]
    s_xyz = [Geographic(latitude=p["lat"], longitude=p["lon"], depth=p["z"])
             for p in s_locations]
    origin = p_xyz[0]

    p_xyz = [p.to_xyz(origin, 0, 90) for p in p_xyz]
    s_xyz = [p.to_xyz(origin, 0, 90) for p in s_xyz]
    p_xyz = [{"x": p.x, "y": p.y, "z": p.z} for p in p_xyz]
    s_xyz = [{"x": p.x, "y": p.y, "z": p.z} for p in s_xyz]

    # Convert from UTCDateTime to float
    _time = min(np.concatenate([p_times, s_times]))
    _p_times = np.array([_t - _time for _t in p_times])
    _s_times = np.array([_t - _time for _t in s_times])

    model = geiger_locate(
        p_times=_p_times, p_locations=p_xyz, s_times=_s_times,
        s_locations=s_xyz, vp=vp, vs=vs, max_it=max_it, convergence=convergence,
        starting_depth=starting_depth, plot=plot)
    
    # convert location back to lat lon
    model_time = _time + model["time"]
    model = Location(x=model["x"], y=model["y"], z=model["z"], origin=origin,
                     strike=0, dip=90)
    model = model.to_geographic()
    return {"lat": model.latitude, "lon": model.longitude, "z": model.depth,
            "time": model_time}


def geiger_locate(p_times, p_locations, s_times, s_locations, vp, vs,
                  max_it=10, convergence=0.1, starting_depth=5.0, plot=False):
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

    def norm_residual(model):
        r = []
        for i, p_loc in enumerate(p_locations):
            r.append(_residual(
                model=model, obs_loc=p_loc, obs_time=p_times[i], velocity=vp))
        for i, s_loc in enumerate(s_locations):
            r.append(_residual(
                model=model, obs_loc=s_loc, obs_time=s_times[i], velocity=vs))
        return np.linalg.norm(r, 2)

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
            s_locations=s_locations, model=model, vp=vp, vs=vs)
        models.append(deepcopy(model))
        residuals.append(norm_residual(model))
        LOGGER.info(
            "Updated model has residual of {0:.2f}".format(residuals[-1]))
        if model_change < convergence:
            LOGGER.info("Model has reached convergence threshold")
            break

    if plot:
        fig = plt.figure()
        scatter_ax = fig.add_subplot(2, 1, 1, projection='3d')
        x = [m["x"] for m in models]
        y = [m["y"] for m in models]
        z = [m["z"] for m in models]
        t = [m["time"] for m in models]
        scatter_ax.plot(x, y, z, "k", linewidth=0.5)
        cmplasma = plt.get_cmap("jet")
        scatter_ax.scatter(x, y, z, c=t, cmap=cmplasma, label="locations")
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


def simple_test():
    p_times = np.array([3, 2.1354, 1.1662, 1.3416])
    p_locations = [
        {"x": 0, "y": 0, "z": 0},
        {"x": 2, "y": 15, "z": 0},
        {"x": 10, "y": 7, "z": 0},
        {"x": 14, "y": 12, "z": 0}]
    s_times = np.array([5, 3.559, 1.9437, 2.2361])
    s_locations = p_locations
    vp = 5
    vs = 3
    model = geiger_locate(
        p_times=p_times, p_locations=p_locations, s_times=s_times,
        s_locations=s_locations, vp=vp, vs=vs, max_it=10, 
        convergence=0.1, starting_depth=50.0, plot=True)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s", 
        level=logging.INFO)
    simple_test()