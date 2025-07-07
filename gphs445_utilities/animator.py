"""
Animator for obspy catalogs

:author: Calum J Chamberlain
:date: 30/5/2019
:licence: LGPL v3
"""

from obspy.imaging.cm import obspy_sequential
from obspy import Catalog, UTCDateTime
import numpy as np
import datetime
import warnings

from matplotlib.dates import AutoDateFormatter, AutoDateLocator, date2num
from matplotlib.cm import get_cmap

import cartopy.crs as ccrs
try:
    from progressbar import ProgressBar
    HAS_PROGRESS = True
except:
    HAS_PROGRESS = False


def _get_plot_coords(catalog):
    lats, lons, mags, colors = ([], [], [], [])
    if len(catalog) == 0:
        return [None], [None], [], []
    for event in catalog:
        if not event.origins:
            msg = ("Event '%s' does not have an origin and will not be "
                   "plotted." % str(event.resource_id))
            warnings.warn(msg)
            continue
        origin = event.preferred_origin() or event.origins[0]
        lats.append(origin.latitude)
        lons.append(origin.longitude)
        try:
            magnitude = event.preferred_magnitude() or event.magnitudes[0]
            mag = magnitude.mag
        except (IndexError, AttributeError):
            mag = 0.1
        mags.append(mag)
        colors.append((origin.get('depth') or np.nan) / 1e3)
    return lats, lons, mags, colors


def _get_event_time(event):
    """
    Get the time of an event - origin preferred, then pick.

    :type event: `obspy.core.event.Event`
    :return: UTCDateTime
    """
    try:
        origin = event.preferred_origin() or event.origins[0]
        t = origin.time
    except (IndexError, AttributeError):
        warnings.warn(
            "No origin time for event, defaulting to pick time")
        try:
            t = min([pick.time for pick in event.picks])
        except ValueError:
            raise AttributeError("No time for event")
    return t


def _blank_map(lons, lats, color, projection="global",
               resolution='110m', continent_fill_color='0.8',
               water_fill_color='1.0', colormap=None, colorbar=None,
               title=None, colorbar_ticklabel_format=None,
               color_label="Depth (km)", logarithmic_color=False,
               proj_kwargs=None, figsize=(20, 20), pad=0.05, norm=None,
               fig=None, map_ax=None, cm_ax=None):
    """
    Plot a map for the region appropriate for the catalog, but do not plot the
    events themselves.

    Used to set-up a figure before animation.

    :type lons: list/tuple of floats
    :param lons: Longitudes of the data points.
    :type lats: list/tuple of floats
    :param lats: Latitudes of the data points.
    :type color: list/tuple of floats (or objects that can be
        converted to floats, like e.g.
        :class:`~obspy.core.utcdatetime.UTCDateTime`)
    :param color: Color information of the individual data points to be
        used in the specified color map (e.g. origin depths,
        origin times).
    :type projection: str, optional
    :param projection: The map projection.
        Currently supported are:

            * ``"global"`` (Will plot the whole world using
              :class:`~cartopy.crs.Mollweide`.)
            * ``"ortho"`` (Will center around the mean lat/long using
              :class:`~cartopy.crs.Orthographic`.)
            * ``"local"`` (Will plot around local events using
              :class:`~cartopy.crs.AlbersEqualArea`.)
            * Any other Cartopy :class:`~cartopy.crs.Projection`. An instance
              of this class will be created using the supplied ``proj_kwargs``.

        Defaults to "global"
    :type resolution: str, optional
    :param resolution: Resolution of the boundary database to use. Will be
        passed directly to the Cartopy module. Possible values are:

            * ``"110m"``
            * ``"50m"``
            * ``"10m"``

        Defaults to ``"110m"``. If you specify another resolutoin, 
        GSHHG will be used
    :type continent_fill_color: Valid matplotlib color, optional
    :param continent_fill_color:  Color of the continents. Defaults to
        ``"0.9"`` which is a light gray.
    :type water_fill_color: Valid matplotlib color, optional
    :param water_fill_color: Color of all water bodies.
        Defaults to ``"white"``.
    :type colormap: str, any matplotlib colormap, optional
    :param colormap: The colormap for color-coding the events as provided
        in `color` kwarg.
        The event with the smallest `color` property will have the
        color of one end of the colormap and the event with the highest
        `color` property the color of the other end with all other events
        in between.    if colorbar is not None:
        show_colorbar = colorbar
    else:
        if len(lons) > 1 and hasattr(color, "__len__") and \
                not isinstance(color, (str, 
                )):
            show_colorbar = True
        else:
            show_colorbar = False
        Defaults to None which will    if colorbar is not None:
        show_colorbar = colorbar
    else:
        if len(lons) > 1 and hasattr(color, "__len__") and \
                not isinstance(color, (str, native_str)):
            show_colorbar = True
        else:
            show_colorbar = False use the default matplotlib colormap.
    :type colorbar: bool, optional    if colorbar is not None:
        show_colorbar = colorbar
    else:
        if len(lons) > 1 and hasattr(color, "__len__") and \
                not isinstance(color, (str, native_str)):
            show_colorbar = True
        else:
            show_colorbar = False
    :param colorbar: When left `Non    if colorbar is not None:
        show_colorbar = colorbar
    else:
        if len(lons) > 1 and hasattr(color, "__len__") and \
                not isinstance(color, (str, native_str)):
            show_colorbar = True
        else:
            show_colorbar = Falsee`, a colorbar is plotted if more than one
        object is plotted. Using `True`/`False` the colorbar can be forced
        on/off.
    :type title: str
    :param title: Title above plot.
    :type colorbar_ticklabel_format: str or function or
        subclass of :class:`matplotlib.ticker.Formatter`
    :param colorbar_ticklabel_format: Format string or Formatter used to format
        colorbar tick labels.
    :type proj_kwargs: dict
    :param proj_kwargs: Keyword arguments to pass to the Cartopy
        :class:`~cartopy.ccrs.Projection`. In this dictionary, you may specify
        ``central_longitude='auto'`` or ``central_latitude='auto'`` to have
        this function calculate the latitude or longitude as it would for other
        projections. Some arguments may be ignored if you choose one of the
        built-in ``projection`` choices.
    :return: Figure
    """
    from obspy.imaging.maps import (
        _CARTOPY_FEATURES, _CARTOPY_RESOLUTIONS, mean_longitude)
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt
    from matplotlib.cm import get_cmap
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize, LogNorm
    from matplotlib.ticker import (
        FormatStrFormatter, Formatter, FuncFormatter, MaxNLocator)

    if isinstance(color[0], (datetime.datetime, UTCDateTime)):
        datetimeplot = True
        color = [date2num(getattr(t, 'datetime', t)) for t in color]
    else:
        datetimeplot = False

    if fig is None:
        fig = plt.figure(figsize=figsize)

    # The colorbar should only be plotted if more then one event is
    # present.
    if colorbar is not None:
        show_colorbar = colorbar
    else:
        if len(lons) > 1 and hasattr(color, "__len__") and \
                not isinstance(color, str):
            show_colorbar = True
        else:
            show_colorbar = False

    if projection == "local":
        ax_x0, ax_width = 0.10, 0.80
    elif projection == "global":
        ax_x0, ax_width = 0.01, 0.98
    else:
        ax_x0, ax_width = 0.05, 0.90

    proj_kwargs = proj_kwargs or {}
    if projection == 'global':
        proj_kwargs['central_longitude'] = np.mean(lons)
        proj = ccrs.Mollweide(**proj_kwargs)
    elif projection == 'ortho':
        proj_kwargs['central_latitude'] = np.mean(lats)
        proj_kwargs['central_longitude'] = mean_longitude(lons)
        proj = ccrs.Orthographic(**proj_kwargs)
    elif projection == 'local':
        if min(lons) < -150 and max(lons) > 150:
            max_lons = max(np.array(lons) % 360)
            min_lons = min(np.array(lons) % 360)
        else:
            max_lons = max(lons)
            min_lons = min(lons)
        lat_0 = max(lats) / 2. + min(lats) / 2.
        lon_0 = max_lons / 2. + min_lons / 2.
        if lon_0 > 180:
            lon_0 -= 360
        deg2m_lat = 2 * np.pi * 6371 * 1000 / 360
        deg2m_lon = deg2m_lat * np.cos(lat_0 / 180 * np.pi)
        if len(lats) > 1:
            height = (max(lats) - min(lats)) * deg2m_lat
            width = (max_lons - min_lons) * deg2m_lon
            margin = pad * (width + height)
            height += margin
            width += margin
        else:
            height = 2.0 * deg2m_lat
            width = 5.0 * deg2m_lon
        # Do intelligent aspect calculation for local projection
        # adjust to figure dimensions
        w, h = fig.get_size_inches()
        aspect = w / h
        if show_colorbar:
            aspect *= 1.2
        if width / height < aspect:
            width = height * aspect
        else:
            height = width / aspect

        proj_kwargs['central_latitude'] = lat_0
        proj_kwargs['central_longitude'] = lon_0
        proj_kwargs['standard_parallels'] = [lat_0, lat_0]
        proj = ccrs.AlbersEqualArea(**proj_kwargs)

    # User-supplied projection.
    elif isinstance(projection, type):
        if 'central_longitude' in proj_kwargs:
            if proj_kwargs['central_longitude'] == 'auto':
                proj_kwargs['central_longitude'] = mean_longitude(lons)
        if 'central_latitude' in proj_kwargs:
            if proj_kwargs['central_latitude'] == 'auto':
                proj_kwargs['central_latitude'] = np.mean(lats)
        if 'pole_longitude' in proj_kwargs:
            if proj_kwargs['pole_longitude'] == 'auto':
                proj_kwargs['pole_longitude'] = np.mean(lons)
        if 'pole_latitude' in proj_kwargs:
            if proj_kwargs['pole_latitude'] == 'auto':
                proj_kwargs['pole_latitude'] = np.mean(lats)

        proj = projection(**proj_kwargs)

    else:
        msg = "Projection '%s' not supported." % projection
        raise ValueError(msg)

    # TODO: Add option to show cumulative plot
    if show_colorbar:
        if map_ax is None:
            map_ax = fig.add_axes([ax_x0, 0.13, ax_width, 0.77], 
                                  projection=proj)
        elif not hasattr(map_ax, "projection") or map_ax.projection != proj:
            geom = map_ax.get_geometry()
            # Remove and rebuild
            map_ax.remove()
            map_ax = fig.add_subplot(*geom, projection=proj)
        if cm_ax is None:
            cm_ax = fig.add_axes([ax_x0, 0.05, ax_width, 0.05])
        plt.sca(map_ax)
    else:
        ax_y0, ax_height = 0.05, 0.85
        if projection == "local":
            ax_y0 += 0.05
            ax_height -= 0.05
        if map_ax is None:
            map_ax = fig.add_axes([ax_x0, ax_y0, ax_width, ax_height],
                                  projection=proj)
        elif not hasattr(map_ax, "projection") or map_ax.projection != proj:
            geom = map_ax.get_geometry()
            # Remove and rebuild
            map_ax.remove()
            map_ax = fig.add_subplot(*geom, projection=proj)

    if projection == 'local':
        x0, y0 = proj.transform_point(lon_0, lat_0, proj.as_geodetic())
        map_ax.set_xlim(x0 - width / 2, x0 + width / 2)
        map_ax.set_ylim(y0 - height / 2, y0 + height / 2)
    else:
        map_ax.set_global()

    # Pick features at specified resolution.
    if resolution in ("10m", "50m", "110m"):
        # Use NaturalEarthFeature interface
        resolution = _CARTOPY_RESOLUTIONS[resolution]
        try:
            borders, land, ocean = _CARTOPY_FEATURES[resolution]
        except KeyError:
            borders = cfeature.NaturalEarthFeature(cfeature.BORDERS.category,
                                                cfeature.BORDERS.name,
                                                resolution,
                                                edgecolor='none',
                                                facecolor='none')
            land = cfeature.NaturalEarthFeature(cfeature.LAND.category,
                                                cfeature.LAND.name, resolution,
                                                edgecolor='face', facecolor='none')
            ocean = cfeature.NaturalEarthFeature(cfeature.OCEAN.category,
                                                cfeature.OCEAN.name, resolution,
                                                edgecolor='face',
                                                facecolor='none')
            _CARTOPY_FEATURES[resolution] = (borders, land, ocean)

        # Draw coast lines, country boundaries, fill continents.
        map_ax.set_facecolor(water_fill_color)
        map_ax.add_feature(ocean, facecolor=water_fill_color)
        map_ax.add_feature(land, facecolor=continent_fill_color)
        map_ax.add_feature(borders, edgecolor='0.75')
        map_ax.coastlines(resolution=resolution, color='0.4')
    else:
        # Use the GSHHG interface
        coast = cfeature.GSHHSFeature(
            scale=resolution, levels=[1], facecolor=continent_fill_color, 
            edgecolor="0.4")
        map_ax.set_facecolor(water_fill_color)
        map_ax.add_feature(coast)
        

    # Draw grid lines - TODO: draw_labels=True doesn't work yet.
    if projection == 'local':
        map_ax.gridlines()
    else:
        # Draw lat/lon grid lines every 30 degrees.
        map_ax.gridlines(xlocs=range(-180, 181, 30), ylocs=range(-90, 91, 30))

    # scatter = map_ax.scatter(lons, lats, marker=marker, s=size, c=color,
    #                          zorder=10, cmap=colormap,
    #                          transform=ccrs.Geodetic())
    if norm is None:
        if logarithmic_color:
            norm = LogNorm(vmin=min(color), vmax=max(color))
        else:
            norm = Normalize(vmin=min(color), vmax=max(color))

    if title:
        plt.suptitle(title)

    # Only show the colorbar for more than one event.
    if show_colorbar:
        if colorbar_ticklabel_format is not None:
            if isinstance(colorbar_ticklabel_format, (str, native_str)):
                formatter = FormatStrFormatter(colorbar_ticklabel_format)
            elif hasattr(colorbar_ticklabel_format, '__call__'):
                formatter = FuncFormatter(colorbar_ticklabel_format)
            elif isinstance(colorbar_ticklabel_format, Formatter):
                formatter = colorbar_ticklabel_format
            locator = MaxNLocator(5)
        else:
            if datetimeplot:
                locator = AutoDateLocator()
                formatter = AutoDateFormatter(locator)
                # Compat with old matplotlib versions.
                if hasattr(formatter, "scaled"):
                    formatter.scaled[1 / (24. * 60.)] = '%H:%M:%S'
            else:
                locator = None
                formatter = None
        cmap = get_cmap(name=colormap)
        cb = ColorbarBase(
            cm_ax, norm=norm, cmap=cmap, orientation='horizontal',
            ticks=locator, format=formatter)
        cb.set_label(color_label)
        # Compat with old matplotlib versions.
        if hasattr(cb, "update_ticks"):
            cb.update_ticks()
    if show_colorbar:
        return fig, map_ax, cm_ax, cb
    return fig, map_ax


class AnimatedCatalog(Catalog):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def time_sort(self, *args, **kwargs):
        """
        Sort the catalog by time - origin-time preferred, but if not found,
        pick time will be used
        """
        _times = []
        for event in self.events:
            t = _get_event_time(event=event)
            _times.append((t, event))
        _times.sort(key=lambda _t: _t[0], *args, **kwargs)
        _, self.events = zip(*_times)
        return self

    def animate(self, projection='global', resolution='l',
                continent_fill_color='0.9', water_fill_color='1.0',
                colormap=None, show=True, title=None, time_step=86400,
                decay=10, interval=10, figsize=(15, 15), fig=None, **kwargs):
        """
        Animate the catalog in time.

        :type time_step: float
        :param time_step: Frame size in seconds.
        :type decay: int
        :param decay:
            How long to keep old events on the plot in frames - will fade out
            (decrease alpha) linearly through this period.
        :type interval: int
        :param interval: Time between frames (ms)

        :return:
        """
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation
        from matplotlib.colors import Normalize
        from collections import deque

        """################## Set up catalog chunks #######################"""
        self.time_sort()
        sub_catalogs, chunk_ends = ([], [])
        catalog_start = _get_event_time(event=self.events[0])
        catalog_end = _get_event_time(event=self.events[-1])
        chunk_start = catalog_start
        while chunk_start < catalog_end:
            chunk_end = chunk_start + time_step
            sub_catalog = Catalog([
                e for e in self.events
                if chunk_start <= _get_event_time(e) < chunk_end])
            sub_catalogs.append(sub_catalog)
            chunk_ends.append(chunk_end)
            chunk_start += time_step
        catalog_deck = deque([Catalog() for _ in range(decay)], maxlen=decay)
        alphas = np.linspace(0., 1., decay)

        frames = int(round(
            ((catalog_end - catalog_start) / time_step) + decay / 2, 0))

        """#################### Set up the empty figure ####################"""
        # lat/lon coordinates, magnitudes, dates
        lats, lons, mags, colors = _get_plot_coords(self)

        # Create the colormap for date based plotting.
        if colormap is None:
            colormap = obspy_sequential
        elif isinstance(colormap, str):
            colormap = get_cmap(colormap)
        norm = Normalize(vmin=min(colors), vmax=max(colors))

        if title is None:
            title = "Animated Catalog"

        min_size = 0.1
        max_size = 20
        min_size_ = min(mags) - 1
        max_size_ = max(mags) + 1

        if fig is None:
            fig, map_ax, _, _ = _blank_map(
                lons, lats, colors, projection=projection,
                resolution=resolution,
                continent_fill_color=continent_fill_color,
                water_fill_color=water_fill_color, colormap=colormap,
                title=title, color_label="Depth (km)", figsize=figsize)
        else:
            map_ax = fig.gca()

        """ ############# Set up the initial scatters ##################### """
        scatters = []
        for _, alpha in zip(catalog_deck, alphas):
            scatters.append(map_ax.scatter(
                [], [], marker="o", s=[], c=[], zorder=10,
                cmap=colormap, transform=ccrs.PlateCarree(), alpha=alpha))
        frame_time = catalog_start - interval
        timestamp = map_ax.text(
            0.05, 0.05, frame_time.strftime("%Y/%m/%d %H:%M:%S.%d"),
            horizontalalignment="left", verticalalignment="bottom",
            transform=map_ax.transAxes)

        if HAS_PROGRESS:
            bar = ProgressBar(max_value=frames)
        """ ####################### Animation function #####################"""
        def update(frame):
            if len(sub_catalogs) > 0:
                catalog_deck.append(sub_catalogs.pop(0))
            else:
                catalog_deck.append(Catalog())

            for i, cat in enumerate(catalog_deck):
                lats, lons, mags, colors = _get_plot_coords(cat)
                frac = [(0.2 + (_i - min_size_)) / (max_size_ - min_size_)
                        for _i in mags]
                size_plot = [(_i * (max_size - min_size)) ** 2 for _i in frac]
                # lons, lats needs to be transformed
                scatters[i].set_offsets(list(zip(lons, lats)))
                # Colors needs to be sequence of rgba tuples
                scatters[i].set_color(colormap(norm(colors), alphas[i]))
                scatters[i].set_sizes(size_plot)
            frame_time = catalog_start + (frame * time_step)
            timestamp.set_text(frame_time.strftime("%Y/%m/%d %H:%M:%S.%d"))
            if HAS_PROGRESS:
                bar.update(frame)
            else:
                print(f"\r{frame}")
            artists = [timestamp, *scatters]
            return artists

        anim = FuncAnimation(
            fig, update, frames=frames, interval=interval, repeat=False,
            blit=True)
        if HAS_PROGRESS:
            bar.finish()

        if show:
            plt.show()
        return anim


if __name__ == '__main__':
    from obspy.clients.fdsn import Client
    from matplotlib import animation

    client = Client("GEONET")
    cat = client.get_events(
        starttime=UTCDateTime(2019, 1, 1), endtime=UTCDateTime(2019, 1, 10),
        maxdepth=120, maxlatitude=-32., minlatitude=-49, minlongitude=164.2,
        maxlongitude=179.)
    print("Downloaded catalog of {0} events".format(len(cat)))
    ani_cat = AnimatedCatalog(cat)
    step = 1800
    fade_out = 3 * 86400
    decay = int(round(fade_out / step))
    fig = ani_cat.animate(
        projection="local", title="GeoNet Catalog", show=False, time_step=step,
        decay=decay, resolution="h")

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=45, metadata=dict(artist="Me"), bitrate=1800)

    fig.save("Test_animation.mp4", writer=writer)


