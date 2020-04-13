""" Functions for event detection in seismograms. """

import numpy as np


def moving_rms(a, n=3):
    """
    Compute the moving root-mean-square
    of a in windows of length n
    
    :type a: numpy.ndarray
    :type n: int
    """
    ret = np.cumsum(a ** 2, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return np.sqrt(ret[n - 1:] / n)


def sta_lta(tr, sta_len, lta_len):
    """
    Compute the STA/LTA ratio for a seismic trace
    
    :type tr: obspy.core.Trace
    :param tr: The trace to compute the STA/LTA for
    :type sta_len: float
    :param sta_len: Length of STA window in seconds
    :type lta_len: float
    :param lta_len: Length of LTA window in seconds
    
    :returns: numpy.ndarray of the detection statistic
    """
    n = tr.stats.npts
    lta_len_samples = int(round(lta_len * tr.stats.sampling_rate))
    sta_len_samples = int(round(sta_len * tr.stats.sampling_rate))
    if sta_len_samples * tr.stats.delta != sta_len:
        print("STA: {0} closest sample length is {1}".format(
            sta_len, sta_len_samples))
    if lta_len_samples * tr.stats.delta != lta_len:
        print("LTA: {0} closest sample length is {1}".format(
            lta_len, lta_len_samples))
    # lta starts one window length in to the trace
    lta = np.zeros(n)
    lta[lta_len_samples - 1:] =  moving_rms(
        a=tr.data, n=lta_len_samples)
    sta = np.zeros(n)
    sta[sta_len_samples - 1:] = moving_rms(
        a=tr.data, n=sta_len_samples)
    detector = np.zeros(n)
    detector[lta_len_samples - 1:] = sta[lta_len_samples - 1:] / lta[lta_len_samples - 1:]
    return detector


def plot_detector(tr, detector):
    """
    Plot the trace with the detection statistic overlaid
    
    :type tr: obspy.core.Trace
    :type detector: numpy.ndarray
    """
    from obspy import Stream
    st = Stream()
    st += tr
    detector_tr = tr.copy()
    detector_tr.data = detector
    detector_tr.stats.network = "STALTA"
    st += detector_tr
    fig = st.plot(equal_scale=False)
    return fig