"""
Functions to help with earthquake statistics.
"""

import numpy as np
import matplotlib.pyplot as plt

from collections import Counter


def calc_b_value(magnitudes, completeness, max_mag=None, plotvar=True):
    """
    Calculate the b-value for a range of completeness magnitudes.
    Take from EQcorrscan.

    Calculates a power-law fit to given magnitudes for each completeness
    magnitude.  Plots the b-values and residuals for the fitted catalogue
    against the completeness values. Computes fits using numpy.polyfit,
    which uses a least-squares technique.

    :type magnitudes: list
    :param magnitudes: Magnitudes to compute the b-value for.
    :type completeness: list
    :param completeness: list of completeness values to compute b-values for.
    :type max_mag: float
    :param max_mag: Maximum magnitude to attempt to fit in magnitudes.
    :type plotvar: bool
    :param plotvar: Turn plotting on or off.

    :rtype: list
    :return:
        List of tuples of (completeness, b-value, residual, number of
        magnitudes used)

    .. Note::
        High "residuals" indicate better fit. Residuals are calculated
        according to the Wiemer & Wyss 2000, Minimum Magnitude of Completeness
        in Earthquake Catalogs: Examples from Alaska, the Western United
        States, and Japan, BSSA.

    .. rubric:: Example

    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.mag_calc import calc_b_value
    >>> client = Client('IRIS')
    >>> t1 = UTCDateTime('2012-03-26T00:00:00')
    >>> t2 = t1 + (3 * 86400)
    >>> catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3)
    >>> magnitudes = [event.magnitudes[0].mag for event in catalog]
    >>> b_values = calc_b_value(magnitudes, completeness=np.arange(3, 7, 0.2),
    ...                         plotvar=False)
    >>> round(b_values[4][1], 1)
    1.0
    >>> # We can set a maximum magnitude:
    >>> b_values = calc_b_value(magnitudes, completeness=np.arange(3, 7, 0.2),
    ...                         plotvar=False, max_mag=5)
    >>> round(b_values[4][1], 1)
    1.0
    """
    b_values = []
    # Calculate the cdf for all magnitudes
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
    for m_c in completeness:
        if m_c >= max_mag or m_c >= max(magnitudes):
            print('Not computing completeness at %s, above max_mag' %
                  str(m_c))
            break
        complete_mags = []
        complete_freq = []
        for i, mag in enumerate(mag_steps):
            if mag >= m_c <= max_mag:
                complete_mags.append(mag)
                complete_freq.append(np.log10(cdf[i]))
        if len(complete_mags) < 4:
            print('Not computing completeness above ' + str(m_c) +
                  ', fewer than 4 events')
            break
        fit = np.polyfit(complete_mags, complete_freq, 1, full=True)
        # Calculate the residuals according to the Wiemer & Wys 2000 definition
        predicted_freqs = [fit[0][1] - abs(fit[0][0] * M)
                           for M in complete_mags]
        r = 100 - ((np.sum([abs(complete_freq[i] - predicted_freqs[i])
                           for i in range(len(complete_freq))]) * 100) /
                   np.sum(complete_freq))
        b_values.append((m_c, abs(fit[0][0]), r, len(complete_mags)))
    if plotvar:
        fig, ax1 = plt.subplots()
        b_vals = ax1.scatter(list(zip(*b_values))[0], list(zip(*b_values))[1],
                             c='k')
        resid = ax1.scatter(list(zip(*b_values))[0],
                            [100 - b for b in list(zip(*b_values))[2]], c='r')
        ax1.set_ylabel('b-value and residual')
        plt.xlabel('Completeness magnitude')
        ax2 = ax1.twinx()
        ax2.set_ylabel('Number of events used in fit')
        n_ev = ax2.scatter(list(zip(*b_values))[0], list(zip(*b_values))[3],
                           c='g')
        fig.legend(handles=(b_vals, resid, n_ev),
                   labels=('b-values', 'residuals', 'number of events'),
                   loc='lower right')
        ax1.set_title('Possible completeness values')
        plt.show()
    return b_values


def calc_max_curv(magnitudes, bin_size=0.5, plotvar=False):
    """
    Calculate the magnitude of completeness using the maximum curvature method.

    :type magnitudes: list or numpy array
    :param magnitudes:
        List of magnitudes from which to compute the maximum curvature which
        will give an estimate of the magnitude of completeness given the
        assumption of a power-law scaling.
    :type bin_size: float
    :param bin_size:
        Width of magnitude bins used to compute the non-cumulative distribution
    :type plotvar: bool
    :param plotvar: Turn plotting on and off

    :rtype: float
    :return: Magnitude at maximum curvature

    .. Note:: Should be used as a guide, often under-estimates Mc.

    .. rubric:: Example

    >>> import numpy as np
    >>> mags = np.arange(3, 6, .1)
    >>> N = 10 ** (5 - 1 * mags)
    >>> magnitudes = [0, 2, 3, 2.5, 2.2, 1.0]  # Some below completeness
    >>> for mag, n in zip(mags, N):
    ...     magnitudes.extend([mag for _ in range(int(n))])
    >>> calc_max_curv(magnitudes, plotvar=False)
    3.0
    """
    min_bin, max_bin = int(min(magnitudes)), int(max(magnitudes) + 1)
    bins = np.arange(min_bin, max_bin + bin_size, bin_size)
    df, bins = np.histogram(magnitudes, bins)
    grad = (df[1:] - df[0:-1]) / bin_size
    # Need to find the second order derivative
    curvature = (grad[1:] - grad[0:-1]) / bin_size
    max_curv = bins[np.argmax(np.abs(curvature))] + bin_size
    if plotvar:
        fig, ax = plt.subplots()
        ax.scatter(bins[:-1] + bin_size / 2, df, color="k",
                   label="Magnitudes")
        ax.axvline(x=max_curv, color="red", label="Maximum curvature")
        ax1 = ax.twinx()
        ax1.plot(bins[:-1] + bin_size / 2, np.cumsum(df[::-1])[::-1],
                 color="k", label="Cumulative distribution")
        ax1.scatter(bins[1:-1], grad, color="r", label="Gradient")
        ax2 = ax.twinx()
        ax2.scatter(bins[1:-2] + bin_size, curvature, color="blue",
                    label="Curvature")
        # Code borrowed from https://matplotlib.org/3.1.1/gallery/ticks_and_
        # spines/multiple_yaxis_with_spines.html#sphx-glr-gallery-ticks-and-
        # spines-multiple-yaxis-with-spines-py
        ax2.spines["right"].set_position(("axes", 1.2))
        ax2.set_frame_on(True)
        ax2.patch.set_visible(False)
        for sp in ax2.spines.values():
            sp.set_visible(False)
        ax2.spines["right"].set_visible(True)

        ax.set_ylabel("N earthquakes in bin")
        ax.set_xlabel("Magnitude")
        ax1.set_ylabel("Cumulative events and gradient")
        ax2.set_ylabel("Curvature")
        fig.legend()
        fig.show()
    return float(max_curv)


def _finalise_figure(fig, **kwargs):  # pragma: no cover
    """
    Internal function to wrap up a figure.
    {plotting_kwargs}
    """
    import matplotlib.pyplot as plt

    title = kwargs.get("title")
    show = kwargs.get("show", True)
    save = kwargs.get("save", False)
    savefile = kwargs.get("savefile", "EQcorrscan_figure.png")
    return_fig = kwargs.get("return_figure", False)
    size = kwargs.get("size", (10.5, 7.5))
    fig.set_size_inches(size)
    if title:
        fig.suptitle(title)
    if save:
        fig.savefig(savefile, bbox_inches="tight")
        print("Saved figure to {0}".format(savefile))
    if show:
        plt.show(block=True)
    if return_fig:
        return fig
    fig.clf()
    plt.close(fig)
    return None


def freq_mag(magnitudes, completeness, max_mag, binsize=0.2, **kwargs):
    """
    Plot a frequency-magnitude histogram and cumulative density plot.

    Currently this will compute a b-value, for a given completeness.
    B-value is computed by linear fitting to section of curve between
    completeness and max_mag.

    :type magnitudes: list
    :param magnitudes: list of float of magnitudes
    :type completeness: float
    :param completeness: Level to compute the b-value above
    :type max_mag: float
    :param max_mag: Maximum magnitude to try and fit a b-value to
    :type binsize: float
    :param binsize: Width of histogram bins, defaults to 0.2
    {plotting_kwargs}

    :returns: :class:`matplotlib.figure.Figure`

    .. Note::
        See :func:`eqcorrscan.utils.mag_calc.calc_b_value` for a least-squares
        method of estimating completeness and b-value. For estimating maximum
        curvature see :func:`eqcorrscan.utils.mag_calc.calc_max_curv`.

    .. rubric:: Example

    >>> from obspy.clients.fdsn import Client
    >>> from obspy import UTCDateTime
    >>> from eqcorrscan.utils.plotting import freq_mag
    >>> client = Client('IRIS')
    >>> t1 = UTCDateTime('2012-03-26T00:00:00')
    >>> t2 = t1 + (3 * 86400)
    >>> catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3)
    >>> magnitudes = [event.preferred_magnitude().mag for event in catalog]
    >>> freq_mag(magnitudes, completeness=4, max_mag=7) # doctest: +SKIP

    .. plot::

        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        from eqcorrscan.utils.plotting import freq_mag
        client = Client('IRIS')
        t1 = UTCDateTime('2012-03-26T00:00:00')
        t2 = t1 + (3 * 86400)
        catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=3)
        magnitudes = [event.preferred_magnitude().mag for event in catalog]
        freq_mag(magnitudes, completeness=4, max_mag=7)
    """
    import matplotlib.pyplot as plt
    # Ensure magnitudes are sorted
    magnitudes.sort()
    # Check that there are no nans or infs
    if np.isnan(magnitudes).any():
        print('Found nan values, removing them')
        magnitudes = [mag for mag in magnitudes if not np.isnan(mag)]
    if np.isinf(magnitudes).any():
        print('Found inf values, removing them')
        magnitudes = [mag for mag in magnitudes if not np.isinf(mag)]
    fig, ax1 = plt.subplots()
    # Set up the bins, the bin-size could be a variables
    bins = np.arange(int(min(magnitudes) - 1), int(max(magnitudes) + 1),
                     binsize)
    n, bins, patches = ax1.hist(magnitudes, bins, facecolor='Black',
                                alpha=0.5, label='Magnitudes')
    ax1.set_ylabel('Frequency')
    ax1.set_ylim([0, max(n) + 0.5 * max(n)])
    plt.xlabel('Magnitude')
    # Now make the cumulative density function
    counts = Counter(magnitudes)
    cdf = np.zeros(len(counts))
    mag_steps = np.zeros(len(counts))
    for i, magnitude in enumerate(sorted(counts.keys(), reverse=True)):
        mag_steps[i] = magnitude
        if i > 0:
            cdf[i] = cdf[i - 1] + counts[magnitude]
        else:
            cdf[i] = counts[magnitude]
    ax2 = ax1.twinx()
    # ax2.scatter(magnitudes, np.log10(cdf), c='k', marker='+', s=20, lw=2,
    ax2.scatter(mag_steps, np.log10(cdf), c='k', marker='+', s=20, lw=2,
                label='Magnitude cumulative density')
    # Now we want to calculate the b-value and plot the fit
    x = []
    y = []
    for i, magnitude in enumerate(mag_steps):
        if completeness <= magnitude <= max_mag:
            x.append(magnitude)
            y.append(cdf[i])
    fit = np.polyfit(x, np.log10(y), 1)
    fit_fn = np.poly1d(fit)
    ax2.plot(magnitudes, fit_fn(magnitudes), '--k',
             label='GR trend, b-value = ' + str(abs(fit[0]))[0:4] +
             '\n $M_C$ = ' + str(completeness))
    ax2.set_ylabel('$Log_{10}$ of cumulative density')
    plt.xlim([min(magnitudes) - 0.1, max(magnitudes) + 0.2])
    plt.ylim([min(np.log10(cdf)) - 0.5, max(np.log10(cdf)) + 1.0])
    plt.legend(loc=2)
    fig = _finalise_figure(fig=fig, **kwargs)  # pragma: no cover
    return fig