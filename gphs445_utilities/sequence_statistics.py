"""
Calum's attempt at calculating earthquake catalog stats - adapted from Zmap

"""

from typing import List, Tuple
from collections import Counter
import logging

import numpy as np

from obspy import UTCDateTime
from obspy.core.event import Catalog, Event


Logger = logging.getLogger(__name__)

# Set some globals for comparison - there is probably a better way to do this.
MIN_MAG = -2 ** 64
MIN_TIME, MAX_TIME = (UTCDateTime(0), UTCDateTime(ns=2 ** 64))


def event_time(event: Event) -> UTCDateTime:
    """ Get a time for an event. """
    try:
        timed_obj = event.preferred_origin() or event.origins[0]
    except IndexError:
        try:
            timed_obj = sorted(event.picks, key=lambda p: p.time)[0]
        except IndexError:
            print("Neither origin nor pick found")
            return UTCDateTime(0)
    return timed_obj.time


def event_magnitude(event: Event) -> float:
    """ Get a magnitude for the event. """
    try:
        magnitude = event.preferred_magnitude() or event.magnitudes[0]
        magnitude = magnitude.mag
    except IndexError:
        magnitude = None
    return magnitude


def _same_sign(a: float, b: float) -> bool:
    """ Check if the sign of two numbers is the same. """
    return abs(a + b) == abs(a) + abs(b)


def _different_sign(a: float, b: float) -> bool:
    """ Reverse of _same_sign. """
    return not _same_sign(a, b)


def _estimate_c_error(
    duration: float, 
    event_times: np.ndarray, 
    k: float, 
    p: float, 
    c: float
) -> float:
    """ Estimate error in p-value """
    qsum = k * ((1 / (duration + c) ** p) - (1 / c ** p))
    psum = (1.0 / (event_times + c)).sum()
    return qsum + p * psum


def _estimate_p_error(
    duration: float, 
    event_times: np.ndarray, 
    k: float, 
    p: float, 
    c: float
) -> float:
    """ Estimate error in c-value """
    qp = 1 - p
    sum_log = np.sum(np.log(event_times + c))
    qsum_log = k / qp ** 2
    qsum_log = qsum_log * (
        (duration + c) ** qp * (
            1 - qp * np.log(duration + c)
        ) - c ** qp * (1 - qp * np.log(c)))
    return qsum_log + sum_log


def _calculate_standard_deviations(
    k: float, 
    c: float, 
    p: float, 
    duration: float,
    elapsed: float = 0.0
) -> Tuple[float, float, float]:
    """
    """
    from math import log
    matrix = np.empty((3, 3))

    f = ((duration + c) ** (-p + 1)) / (-p + 1)
    h = ((elapsed + c) ** (-p + 1)) / (-p + 1)
    matrix[0][0] = (1 / k) * (f - h)

    f = (duration + c) ** -p
    h = (elapsed + c) ** -p
    matrix[0][1] = f - h
    
    f = (-(duration + c) ** (-p + 1)) * (
        ((log(duration + c)) / (-p + 1)) - (1 / ((-p + 1) ** 2)))
    h = (-(elapsed + c) ** (-p + 1)) * (
        ((log(elapsed + c)) / (-p + 1)) - (1 / ((-p + 1) ** 2)))
    matrix[0][2] = f - h
    
    matrix[1][0] = matrix[0][1]
    
    f = ((duration + c) ** (-p - 1)) / (p + 1)
    h = ((elapsed + c) ** (-p - 1)) / (p + 1)
    matrix[1][1] = (-k) * (p ** 2) * (f - h)
    
    f = ((duration + c) ** (-p)) * (
        ((log(duration + c)) / (-p)) - (1 / (p ** 2)))
    h = ((elapsed + c) ** (-p)) * (((log(elapsed + c)) / (-p)) - (1 / (p ** 2)))
    matrix[1][2] = (k * p) * (f - h)
    
    matrix[2][0] = matrix[0][2]
    matrix[2][1] = matrix[1][2]
    
    f10 = ((duration + c) ** (-p + 1)) * ((log(duration + c)) ** 2) / (-p + 1)
    f11 = (2 * ((duration + c) ** (-p + 1))) / ((-p + 1) ** 2)
    f12 = (log(duration + c)) - (1 / (-p + 1))
    f9 = f10 - (f11 * f12)
    
    h10 = ((elapsed + c) ** (-p + 1)) * ((log(elapsed + c)) ** 2) / (-p + 1)
    h11 = (2 * ((elapsed + c) ** (-p + 1))) / ((-p + 1) ** 2)
    h12 = (log(elapsed + c)) - (1 / (-p + 1))
    h9 = h10 - (h11 * h12)
    matrix[2][2] = (k) * (f9 - h9)

    matrix = np.linalg.inv(matrix)
    k_std = matrix[0][0] ** .5
    c_std = matrix[1][1] ** .5
    p_std = matrix[2][2] ** .5
    return k_std, p_std, c_std


    def _estimate_p_error(
        duration: float, 
        event_times: np.ndarray, 
        k: float, 
        p: float, 
        c: float
    ) -> float:
        """ Estimate error in p-value """
        qsum = k * ((1 / (duration + c) ** p) - (1 / c ** p))
        psum = (1.0 / event_times + c).sum()
        return qsum + p * psum

    def _estimate_c_error(
        duration: float, 
        event_times: np.ndarray, 
        k: float, 
        p: float, 
        c: float
    ) -> float:
        """ Estimate error in c-value """
        qp = 1 - p
        sum_log = np.sum(np.log(event_times + c))
        qsum_log = k / qp ** 2
        qsum_log = qsum_log * (
            (duration + c) ** qp * (
                1 - qp * np.log(duration + c)
            ) - c ** qp * (1 - qp * np.log(c)))
        return qsumln + sumln


    def _estimate_p_error(
        duration: float, 
        event_times: np.ndarray, 
        k: float, 
        p: float, 
        c: float
    ) -> float:
        """ Estimate error in p-value """
        qsum = k * ((1 / (duration + c) ** p) - (1 / c ** p))
        psum = (1.0 / event_times + c).sum()
        return qsum + p * psum

    def _estimate_c_error(
        duration: float, 
        event_times: np.ndarray, 
        k: float, 
        p: float, 
        c: float
    ) -> float:
        """ Estimate error in c-value """
        qp = 1 - p
        sum_log = np.sum(np.log(event_times + c))
        qsum_log = k / qp ** 2
        qsum_log = qsum_log * (
            (duration + c) ** qp * (
                1 - qp * np.log(duration + c)
            ) - c ** qp * (1 - qp * np.log(c)))
        return qsumln + sumln


class Aftershocks(object):
    # Parameters for Omori solving
    _p_initial, _c_initial = (1.1, 0.1)
    _minimum_p_error, _minimum_c_error = (1e-3, 1e-3)
    _pstep_initial, _cstep_initial = (0.05, 0.05)
    _minimum_pstep, _minimum_cstep = (1e-4, 1e-4)
    _step_reduction_factor = .9
    _maxloops = 500

    def __init__(self, catalog: Catalog, mainshock: Event):
        self.catalog = sorted(catalog, key=lambda e: event_time(e))
        self.mainshock = mainshock

    def __repr__(self):
        return ("Aftershocks(catalog=<Catalog of {0} events>, "
                "mainshock=<M {1:.2f}>)".format(
                    len(self.catalog), self.mainshock_magnitude))

    @property
    def mainshock_time(self):
        return event_time(self.mainshock)

    @property
    def mainshock_magnitude(self):
        return event_magnitude(self.mainshock)

    @property
    def mag_time(self):
        for event in self.catalog:
            magnitude = event_magnitude(event)
            ev_time = event_time(event)
            if ev_time == UTCDateTime(0):
                ev_time = None
            yield (magnitude, ev_time)

    def plot_rate(
        self,
        endtime: UTCDateTime = None, 
        min_magnitude: float = None,
        k: float = None, 
        c: float = None, 
        p: float = None):
        """
        Plot the cumulative number of events in the catalogue.

        Parameters
        ----------
        endtime
            The endtime for the catalogue window to plot for.
        min_magnitude
            The smallest magnitude to plot
        k
            Omori parameter K
        c
            Omori parameter c
        p
            Omori parameter p
        
        Returns
        -------
        A Figure of the plot.
        """
        import matplotlib.pyplot as plt

        min_magnitude = min_magnitude or MIN_MAG
        endtime = endtime or MAX_TIME

        fig, ax = plt.subplots(1)
        eq_times = [(eq_time - self.mainshock_time) / 86400.0
                    for mag, eq_time in self.mag_time 
                    if mag and mag >= min_magnitude and eq_time <= endtime]
        ax.scatter(eq_times, [i for i in range(len(eq_times))], 
                   color='k', marker="o")
        ax.set_xlabel("Time since mainshock (days)")
        ax.set_ylabel("Cumulative number of events")

        # Plot the fitted Omori decay
        fitted_times = np.linspace(0, max(eq_times), 200)
        fitted_numbers = k * (
            c ** (1 - p) - (fitted_times + c) ** (1 -p)) / (p - 1)
        ax.plot(fitted_times, fitted_numbers, "r-", 
                label="Omori k={0:.2f}, c={1:.2f}, p={2:.2f}".format(k, c, p))
        ax.legend()
        return fig

    def plot_magnitude_frequency(
        self, 
        starttime: UTCDateTime = None, 
        endtime: UTCDateTime = None, 
        mc: float = None, 
        a: float = None, 
        b: float = None, 
        bin_width: float = 0.1
    ):
        """
        Plot the magnitude-frequency distribution for the catalogue.

        Parameters
        ----------
        starttime
            The starttime for the catalogue window to plot for.
        endtime
            The endtime for the catalogue window to plot for.
        mc
            Magnitude of completeness of the catalogue
        a
            Gutenberg-Richter a-value in the cumulative form
        b
            Gutenberg-Richter b-value
        bin_width
            Width of magnitude bin for histogram plot.
        
        Returns
        -------
        A Figure of the plot.
        """
        import matplotlib.pyplot as plt

        starttime = starttime or MIN_TIME
        endtime = endtime or MAX_TIME

        fig, ax = plt.subplots(1)
        mags = self._get_magnitudes(starttime, endtime)
        mags.sort(reverse=True)
        max_mag = int(mags[0] + 1)
        min_mag = int(mags[-1])
        lower_bins = np.arange(max_mag, min_mag, -bin_width)
        cumulative = []
        for lower_bound in lower_bins:
            n_events = 0
            for mag in mags:
                if mag >= lower_bound:
                    n_events += 1
                else:
                    break
            cumulative.append(n_events)
        bin_points = [b - (bin_width * 0.5) for b in lower_bins]
        ax.scatter(bin_points, cumulative, color='k', marker="x")
        ax.set_yscale('log')
        ax.set_ylim((0.9, max(cumulative) + 1))
        ax.set_ylabel("Cumulative events")
        ax.set_xlabel("Magnitude")
        # Plot the fit
        if a and b:
            predicted_n = 10 ** (a - b * lower_bins)
            ax.plot(lower_bins, predicted_n, 'k', 
                    label="b={0:.2f}, a={1:.2f}".format(b, a))
        # Plot the completeness
        if mc:
            ax.vlines([mc], 0, 1, transform=ax.get_xaxis_transform(),
                      colors='r', linestyles="dashed",
                      label="$M_C$ = {0:.02f}".format(mc))
            ax.legend()
        return fig

    def _get_magnitudes(
        self,
        starttime: UTCDateTime = None,
        endtime: UTCDateTime = None
    ) -> List:
        starttime = starttime or MIN_TIME
        endtime = endtime or MAX_TIME
        return [mag for mag, ev_time in self.mag_time 
                if starttime <= ev_time <= endtime and mag is not None]


    def goodness_of_fit(
        self,
        max_mag: float = None,
        starttime: UTCDateTime = None,
        endtime: UTCDateTime = None,
        completeness_step: float = 0.1,
        completeness_range: Tuple[float, float] = None,
    ) -> float:
        """
        Calculate magnitude of completeness using goodness-of-fit method.

        Parameters
        ----------
        starttime
            Earliest event-time to include in completeness calculation. 
            If None, then all events up to `endtime` will be used.
        endtime
            Last event-time to include in completeness calculation. If None 
            then all events after `starttime` will be used.
        completeness_step
            Magnitude interval for computing completeness across.

        Returns
        -------
        Magnitude of completeness.
        """
        m_c, r = None, 9999
        magnitudes = self._get_magnitudes(starttime, endtime)
        completeness_range = completeness_range or (
            min(magnitudes), max(magnitudes))
        cut = completeness_range[0]
        while cut <= completeness_range[1]:
            _, _r, _ = self.b_value(
                m_c=cut, max_mag=max_mag, starttime=starttime, endtime=endtime)
            if _r is not None and _r < r:
                m_c, r = cut, _r
            cut += completeness_step
        return m_c


    def max_curv(
        self,
        max_mag: float = None,
        starttime: UTCDateTime = None,
        endtime: UTCDateTime = None, 
    ) -> float:
        """
        Calculate magnitude of completeness using maximum curvature method.

        Parameters
        ----------
        starttime
            Earliest event-time to include in completeness calculation. 
            If None, then all events up to `endtime` will be used.
        endtime
            Last event-time to include in completeness calculation. If None 
            then all events after `starttime` will be used.

        Returns
        -------
        Magnitude of completeness.

        Warning
        -------
            This method often under-estimates completeness.  It is common to 
            add a constant to this.
        Note
        ----
        Code adapted from eqcorrscan.utils.mag_calc.calc_max_curv
        """
        max_mag = max_mag or 100
        magnitudes = [m for m in self._get_magnitudes(starttime, endtime) 
                      if m < max_mag]
        counts = Counter(magnitudes)
        df = np.zeros(len(counts))
        mag_steps = np.zeros_like(df)
        grad = np.zeros(len(counts) - 1)
        grad_points = grad.copy()
        for i, magnitude in enumerate(sorted(counts.keys(), reverse=True)):
            mag_steps[i] = magnitude
            if i > 0:
                df[i] = counts[magnitude] + df[i - 1]
            else:
                df[i] = counts[magnitude]
        for i, val in enumerate(df):
            if i > 0:
                grad[i - 1] = (val - df[i - 1]) / (mag_steps[i] - mag_steps[i - 1])
                grad_points[i - 1] = mag_steps[i] - ((mag_steps[i] -
                                                    mag_steps[i - 1]) / 2.0)
        # Need to find the second order derivative
        curvature = np.zeros(len(grad) - 1)
        curvature_points = curvature.copy()
        for i, _grad in enumerate(grad):
            if i > 0:
                curvature[i - 1] = (_grad - grad[i - 1]) / (
                    grad_points[i] - grad_points[i - 1])
                curvature_points[i - 1] = grad_points[i] - (
                    (grad_points[i] - grad_points[i - 1]) / 2.0)
        return curvature_points[np.argmax(abs(curvature))]

    def b_value(
        self, 
        m_c: float, 
        max_mag: float = None,
        starttime: UTCDateTime = None,
        endtime: UTCDateTime = None, 
    ) -> Tuple[float, float, float]:
        """ 
        Calculate b-value 
        
        Parameters
        ----------
        m_c
            Magnitude of completeness.
        max_mag
            Maximum magnitude to fit up to
        starttime
            Earliest event-time to include in completeness calculation. 
            If None, then all events up to `endtime` will be used.
        endtime
            Last event-time to include in completeness calculation. If None 
            then all events after `starttime` will be used.

        Returns
        -------
        b-value, residual, a-value
        """
        magnitudes = self._get_magnitudes(starttime, endtime)
        counts = Counter(magnitudes)
        cdf = np.zeros(len(counts))
        mag_steps = np.zeros(len(counts))
        for i, magnitude in enumerate(sorted(counts.keys(), reverse=True)):
            mag_steps[i] = magnitude
            if i > 0:
                cdf[i] = cdf[i - 1] + counts[magnitude]
            else:
                cdf[i] = counts[magnitude]

        if not max_mag:
            max_mag = max(magnitudes)
        if m_c >= max_mag or m_c >= max(magnitudes):
            Logger.warning("Can't compute, m_c is above max_mag")
            return None, None, None
        complete_mags = []
        complete_freq = []
        for i, mag in enumerate(mag_steps):
            if m_c <= mag <= max_mag:
                complete_mags.append(mag)
                complete_freq.append(np.log10(cdf[i]))
        if len(complete_mags) < 4:
            Logger.warning('Not computing, fewer than 4 events')
            return None, None, None
        fit = np.polyfit(complete_mags, complete_freq, 1, full=True)
        # Calculate the residuals according to the Wiemer & Wys 2000 definition
        predicted_freqs = [fit[0][1] - abs(fit[0][0] * M)
                        for M in complete_mags]
        r = 100 - ((np.sum([abs(complete_freq[i] - predicted_freqs[i])
                        for i in range(len(complete_freq))]) * 100) /
                np.sum(complete_freq))
        return abs(fit[0][0]), r, fit[0][1]

    def omori(
        self, 
        endtime: UTCDateTime = None, 
        plot: bool = False, 
        min_magnitude: int = None,
    ) -> Tuple[float, float, float, float, float, float]:
        """
        Calculate the Omori decay parameters for a section of catalogue.

        Follows the modified Omori law:

        n(t) = K / (t + c) ** p

        Where n(t) is the rate of earthquakes at time t.  K, c and p are 
        constants.

        Parameters
        ----------
        endtime
            The endtime for the catalogue window to calculate for.
        plot
            Whether to plot the fit or not
        min_magnitude
            The smallest magnitude to use for calculating the fit.

        Returns
        -------
        Paramters K, c and p, and their associated stanard deviations
        K = productivity
        c = time adjustment
        p = falloff rate
        
        Example usage
        -------------
        > K, K_std, c, c_std, p, p_std = atershocks.omori()
        """
        from math import copysign

        endtime = endtime or MAX_TIME
        min_magnitude = min_magnitude or MIN_MAG

        event_times = [(ev_time - self.mainshock_time) / 86400.0
                       for mag, ev_time in self.mag_time
                       if mag and mag >= min_magnitude and ev_time <= endtime]
        event_times = np.array(event_times)
        n = len(event_times)
        duration = event_times.max()

        # Code adapted from Zmap (internals were working in zmap, but couldn't call it.)
        def estimate_k(p, c):
            return (1 - p) * n / ((duration + c) ** (1 - p) - c ** (1 - p))

        p, c = (self._p_initial, self._c_initial)
        pstep, cstep = (self._pstep_initial, self._cstep_initial)
        k = estimate_k(p, c)  # Apparently K doesn't change?

        # Iteratively solve inversion
        previous_p_error, previous_c_error = (None, None)
        i = 0
        while i < self._maxloops:
            if p == 1.0:
                Logger.warning("P == 1.0, adjusting")
            
            c_error = _estimate_c_error(duration, event_times, k, p, c)
            p_error = _estimate_p_error(duration, event_times, k, p, c)

            # Break conditions
            if abs(c_error) < self._minimum_c_error:
                break
            elif abs(p_error) < self._minimum_p_error:
                break
            elif cstep <= self._minimum_cstep:
                break
            elif pstep <= self._minimum_pstep:
                break
            
            # Adjust P
            if previous_p_error and _different_sign(previous_p_error, p_error) and pstep >= self._minimum_pstep:
                pstep *= self._step_reduction_factor
            p = p - copysign(pstep, p_error)
            previous_p_error = p_error
            # Adjust C
            if previous_c_error and _different_sign(previous_c_error, c_error) and cstep >= self._minimum_cstep:
                cstep *= self._step_reduction_factor
            c = c - copysign(cstep, c_error)
            if c <= 0:
                c = cstep
            previous_c_error = c_error
            i += 1
        else:
            Logger.error("Did not converge")
        
        # Calculate standard deviations in k, p, c
        k_std, p_std, c_std = _calculate_standard_deviations(
            k=k, p=p, c=c, duration=duration)

        if plot:
            fig = self.plot_rate(
                endtime=endtime, min_magnitude=min_magnitude, k=k, c=c, p=p)
            fig.show()
        return k, k_std, c, c_std, p, p_std


def make_plots_for_catalog(catalog: Catalog, cutoff: float = None):
    mainshock = catalog[0]
    mainshock_mag = (
        mainshock.preferred_magnitude() or mainshock.magnitudes[0]).mag
    for ev in catalog:
        try:
            mag = (ev.preferred_magnitude() or ev.magnitudes[0]).mag
        except IndexError:
            continue
        if mag > mainshock_mag:
            mainshock = ev
            mainshock_mag = mag

    aft = Aftershocks(catalog, mainshock)
    magnitudes = [event_magnitude(ev) for ev in catalog if event_magnitude(ev)]
    magnitudes.sort()
    max_mag = magnitudes[-2]  # Do not try to fit all the way up to the mainshock
    mc = aft.max_curv()
    # mc = aft.goodness_of_fit(
    #     max_mag=max_mag, completeness_range=(0., 3.))
    b, r, a = aft.b_value(m_c=mc, max_mag=max_mag)
    cutoff = cutoff or mc
    k, k_std, c, c_std, p, p_std = aft.omori(min_magnitude=cutoff)
    Logger.info(
        "Stats: mc {mc}, b {b}, residual {r}, a {a}, k {k}, c {c}, "
        "p {p}".format(mc=mc, b=b, r=r, a=a, k=k, c=c, p=p))
    mf_fig = aft.plot_magnitude_frequency(mc=mc, a=a, b=b)
    omori_fig = aft.plot_rate(k=k, c=c, p=p, min_magnitude=cutoff)
    return mf_fig, omori_fig


if __name__ == "__main__":
    import argparse
    import glob
    from obspy import read_events
    import matplotlib.pyplot as plt

    plt.ioff()
    logging.basicConfig(
        level="INFO",
        format="%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")

    parser = argparse.ArgumentParser(
        description="Make plots of Omori and GR for a given sequence")
    parser.add_argument(
        "-s", "--sequence", type=str, help="Sequence to make plots for.",
        required=True)
    parser.add_argument(
        "-c", "--cutoff-magnitude", type=float, 
        help="Cutoff magnitude for Omori plots", default=2.5)
    
    args = vars(parser.parse_args())

    catalog_files = glob.glob(
        "../Sequences/{0}/repicked_catalog_end_*.xml".format(args["sequence"]))
    catalog_files.sort()
    geonet_file = "../Sequences/{0}/{0}_geonet.xml".format(args["sequence"])

    # Make GeoNet plots
    geonet_cat = read_events(geonet_file)
    geonet_cat.events.sort(key=lambda e: e.origins[0].time)
    mainshock = geonet_cat[0]
    for ev in geonet_cat:
        mainshock_mag = mainshock.preferred_magnitude() or mainshock.magnitudes[0]
        try:
            mag = ev.preferred_magnitude() or ev.magnitudes[0]
        except IndexError:
            continue
        if mag.mag > mainshock_mag.mag:
            mainshock = ev
    geonet_aft = Aftershocks(geonet_cat, mainshock=mainshock)
    min_time, max_time = (
        geonet_cat[0].origins[0].time, geonet_cat[-1].origins[0].time)
    mf_fig, omori_fig = make_plots_for_catalog(
        geonet_cat, cutoff=args["cutoff_magnitude"])
    omori_fig.suptitle("GeoNet $M>{0}$ {1} {2}-{3}".format(
        args["cutoff_magnitude"], args["sequence"], 
        min_time.strftime("%Y/%m/%dT%H:%M:%S"), 
        max_time.strftime("%Y/%m/%dT%H:%M:%S")))
    mf_fig.suptitle("GeoNet {0} {1}-{2}".format(
        args["sequence"], min_time.strftime("%Y/%m/%dT%H:%M:%S"), 
        max_time.strftime("%Y/%m/%dT%H:%M:%S")))
    omori_fig.savefig("../Sequences/{0}/GeoNet_omori.png".format(
        args["sequence"]))
    mf_fig.savefig("../Sequences/{0}/GeoNet_mf.png".format(args["sequence"]))

    # Make Matched plots
    for catalog_file in catalog_files:
        Logger.info("Working on file {0}".format(catalog_file))
        file_end = catalog_file.split(
            "repicked_catalog_end_")[-1].split(".xml")[0]
        cat = read_events(catalog_file)
        cat.events.sort(key=lambda e: e.origins[0].time)
        min_time, max_time = cat[0].origins[0].time, cat[-1].origins[0].time
        mainshock = cat[0]
        for ev in cat:
            mainshock_mag = mainshock.preferred_magnitude() or mainshock.magnitudes[0]
            try:
                mag = ev.preferred_magnitude() or ev.magnitudes[0]
            except IndexError:
                continue
            if mag.mag > mainshock_mag.mag:
                mainshock = ev
        aft = Aftershocks(cat, mainshock=mainshock)
        try:
            mf_fig, omori_fig = make_plots_for_catalog(
                cat, cutoff=args["cutoff_magnitude"])
            mf_fig.suptitle("Matched {0} {1}-{2}".format(
                args["sequence"], min_time.strftime("%Y/%m/%dT%H:%M:%S"), 
                max_time.strftime("%Y/%m/%dT%H:%M:%S")))
            mf_fig.savefig("../Sequences/{0}/MF_{1}.png".format(
                args["sequence"], file_end))
            omori_fig.suptitle("Matched $M>{0}$ {1} {2}-{3}".format(
                args["cutoff_magnitude"], args["sequence"], 
                min_time.strftime("%Y/%m/%dT%H:%M:%S"), 
                max_time.strftime("%Y/%m/%dT%H:%M:%S")))
            omori_fig.savefig("../Sequences/{0}/Omori_{1}.png".format(
                args["sequence"], file_end))
        except Exception as e:
            print(e)
        plt.close("all")
