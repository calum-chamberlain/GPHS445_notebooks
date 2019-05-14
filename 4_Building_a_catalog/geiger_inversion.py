"""
Simple linearised inversion of seismic phase arrival times to get source
location in a homogeneous half-space. Used for teaching, not production ready.

:author: Calum Chamberlain
:date: 13/5/2019
"""

import numpy as np
import logging

LOGGER = logging.getLogger("Geiger")


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