"""
Plot picks onto waveforms for GeoNet events.


"""

import numpy as np

from obspy import Stream
from obspy.core.event import Event
from matplotlib.pyplot import Figure


def plot_polarities(event):
    import warnings
    import matplotlib.pyplot as plt
    import mplstereonet
    from obspy.imaging.beachball import aux_plane

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='stereonet')
    
    up_toa, up_baz, down_toa, down_baz = ([], [], [], [])

    for pick in event.picks:
        try:
            polarity = pick.polarity
        except AttributeError:
            continue
        if polarity == "undecidable":
            continue
        origin = event.preferred_origin() or event.origins[-1]
        # We want the take-off angle and azimuth
        toa, baz = (None, None)
        for arrival in origin.arrivals:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                _pick = arrival.pick_id.get_referred_object()
            if _pick is not None and _pick == pick:
                toa = arrival.takeoff_angle
                baz = arrival.azimuth
        if not (toa and baz):
            continue
        if 0. <= toa < 90.:
            toa = 90. - toa  # complement for downward angles
        elif 90. <= toa <= 180.:
            toa = 270. - toa  # project upward angles
        baz -= 90 # Calculate strike azi from direct (dip-pointing) azi
        if polarity == "positive":
            up_toa.append(toa)
            up_baz.append(baz)
        elif polarity == "negative":
            down_toa.append(toa)
            down_baz.append(baz)
    ax.rake(up_baz, up_toa, 90, "ro", label="Compressional")
    ax.rake(down_baz, down_toa, 90, "bo", label="Dilatational")
    for fm in event.focal_mechanisms:
        if fm.nodal_planes:
            if fm.nodal_planes.nodal_plane_1:
                _np = fm.nodal_planes.nodal_plane_1
                ax.plane(_np.strike, _np.dip, "k")
            if fm.nodal_planes.nodal_plane_2:
                _np = fm.nodal_planes.nodal_plane_2
                ax.plane(_np.strike, _np.dip, "k")
            elif fm.nodal_planes.nodal_plane_1:
                # Calculate the aux plane
                _np = fm.nodal_planes.nodal_plane_1
                _str, _dip, _rake = aux_plane(_np.strike, _np.dip, _np.rake)
                ax.plane(_str, _dip, "k")
    ax.legend()
    return fig


def get_geonet_waveforms(
    event: Event,
    delay: float = 0.,
    length: float = 120.,
    all_components: bool = False,
) -> Stream:
    """
    Get picked GeoNet waveforms for an event.

    Parameters
    ----------
    event:
        The event to get waveforms for
    delay:
        Delay in seconds relative to the origin-time of the event to download 
        waveforms from
    length:
        Total length in seconds for each channel to download
    all_components:
        Whether to download all components (True), or just those picked (False)

    Returns
    -------
    Stream downloaded.
    """
    from obspy.clients.fdsn import Client

    client = Client("GEONET")
    
    seed_ids = {pick.waveform_id.get_seed_string() for pick in event.picks}
    if all_components:
        seed_ids = {sid[0:-1] + "?" for sid in seed_ids}
    try:
        origin_time = (event.preferred_origin() or event.origins[0]).time
    except IndexError:
        try:
            origin_time = sorted(event.picks, key=lambda p: p.time)[0].time
        except IndexError:
            raise NotImplementedError("No origin time or pick times found")    
    st = Stream()
    for seed_id in seed_ids:
        print(f"Downloading for {seed_id}")
        try:
            st += client.get_waveforms(
                *seed_id.split('.'), starttime=origin_time + delay, 
                endtime=origin_time + length + delay)
        except Exception as e:
            print(f"Some error downloading {seed_id}: \n{e}")
    st.trim(origin_time + delay, origin_time + delay + length)
    st.merge().sort(["starttime"])
    return st


def plot_picked(event: Event, st: Stream) -> Figure:
    """
    Plot picks onto waveforms for an event.
    
    Parameters
    ----------
    event:
        The event with picks to plot
    st:
        The stream to plot onto

    Returns
    -------
    Figure of the plot
    """
    import matplotlib.pyplot as plt

    pick_color = {'P': 'r', 'S': 'b'}
    
    fig, axes = plt.subplots(nrows=len(st), figsize=(15, 8), sharex=True)
    for tr, ax in zip(st, axes):
        x = np.arange(0, tr.stats.npts)
        x = x * tr.stats.delta
        x = [(tr.stats.starttime + _).datetime for _ in x]
        ax.plot(x, tr.data, 'k')
        ax.set_ylabel(tr.id, rotation=0)
        ax.set_yticks([])
        # Find matching picks and plot them
        for pick in event.picks:
            if pick.waveform_id.get_seed_string() == tr.id:
                if x[0] < pick.time.datetime < x[-1]:
                    ax.axvline(x=pick.time.datetime, 
                               color=pick_color[pick.phase_hint])
    fig.subplots_adjust(wspace=0, hspace=0)
    return fig