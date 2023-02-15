"""
Functions to help interact with focal mechanisms.

"""

import matplotlib.pyplot as plt
import mplstereonet
from mplstereonet import stereonet_math

from obspy.core.event import Event

from typing import List
from collections import namedtuple

# Hack to cope with numpy change
import numpy as np

np.float = float

Polarity = namedtuple("Polarity", ("azimuth", "toa", "polarity", "station"))

UPS = ("U", "up", "positive", "compressional")
DOWNS = ("D", "down", "negative", "dilatational")

POLARITY_PROPS = {
    "up": {"label": "Compressional", "color": "red", "marker": "^"},
    "down": {"label": "Dilatational", "marker": "v", 
             "markeredgecolor": 'black', "markerfacecolor": "white"},
    "unknown": {"label": "Unknown", "color": "black", "marker": "x"}}


class NodalPlane():
    def __init__(self, strike, dip, rake):
        self.strike = strike % 360
        # TODO: I could take care of this being out of bounds.
        assert 0 <= dip <= 90, "Dip must be between 0 and 90"
        self.dip = dip
        rake = rake % 360
        if rake > 180:
            rake = rake - 360
        self.rake = rake

    def __repr__(self):
        return (
            f"NodalPlane(strike={self.strike}, dip={self.dip}, "
            f"rake={self.rake})")
    # TODO: Implement finding aux plane and P,T,B axes.

    def plot(self, ax = None, show: bool = False, color: str = "k", 
             label: str = None, markeredgecolor: str = None):
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="stereonet")
        rake_label, pole_label = None, None
        if label:
            rake_label, pole_label = f"{label} rake", f"{label} pole"
        artists = []
        artists.extend(ax.plane(
            self.strike, self.dip, color=color, label=label))
        # Note rake is measured anticlockwise from horizontal for seismologists,
        # but clockwise for geologists (as in mplstereonet)
        rake = -1 * self.rake
        # if self.strike < 180:
        #     rake = -1 * rake
        artists.extend(ax.rake(
            self.strike, self.dip, rake, color=color, label=rake_label,
            markeredgecolor=markeredgecolor))
        artists.extend(ax.pole(
            self.strike, self.dip, color=color, marker="s", label=pole_label,
            markeredgecolor=markeredgecolor))
        if show:
            plt.show()
        return artists


class FocalMechanism():
    def __init__(
        self, 
        nodal_plane_1: NodalPlane = None,
        nodal_plane_2: NodalPlane = None,
        polarities: List[Polarity] = None
    ):
        self.nodal_plane_1 = nodal_plane_1 or NodalPlane(45, 45, 45)
        self.nodal_plane_2 = nodal_plane_2 or NodalPlane(35, 35, 35)
        self.polarities = polarities or []

    @classmethod
    def from_event(cls, event: Event):
        try:
            origin = event.preferred_origin() or event.origins[-1]
        except IndexError:
            raise NotImplementedError("Event needs an origin")
        polarities = []

        for pick in event.picks:
            if pick.polarity and pick.phase_hint.startswith("P"):
                # Get the arrival
                pick_seed_id = pick.waveform_id.get_seed_string()
                print(f"Found polarity of {pick.polarity} for {pick_seed_id}")
                for arr in origin.arrivals:
                    arr_pick = arr.pick_id.get_referred_object()
                    if arr_pick and arr_pick.waveform_id.get_seed_string() == pick_seed_id:
                        if arr.phase == "P":
                            if arr.takeoff_angle < 0:
                                toa = abs(arr.takeoff_angle)
                                az = arr.azimuth % 360
                            else:
                                toa = arr.takeoff_angle
                                az = (arr.azimuth + 180) % 360
                            polarity = Polarity(az, toa,
                                                pick.polarity,
                                                station=pick_seed_id)
                            polarities.append(polarity)
                            break
                else:
                    print(f"No arrival found for polarity pick on {pick_seed_id}")
        return cls(polarities=polarities)

    def plot(self, show: bool = True, label: bool = True) -> plt.Figure:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="stereonet")

        self.plot_polarities(ax=ax, show=False, label=label)
        for nodal_plane in (self.nodal_plane_1, self.nodal_plane_2):
            nodal_plane.plot(ax=ax, show=False)
        ax.legend()
        ax.grid()
        if show:
            plt.show(block=True)
        return fig

    def plot_polarities(self, ax = None, show: bool = False, label: bool = True) -> plt.Axes:
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="stereonet")
        for polarity in self.polarities:
            toa, baz = polarity.toa, polarity.azimuth
            assert 0 <= toa <= 180, f"Take off angle ({toa}) does not lie between 0 and 180"
            if 0. <= toa < 90.:
                toa = 90. - toa  # complement for downward angles
            elif 90. <= toa <= 180.:
                toa = 270. - toa  # project upward angles
            baz -= 90 # Calculate strike azi from direct (dip-pointing) azi
            if polarity.polarity in UPS:
                plot_kwargs = POLARITY_PROPS["up"]
            elif polarity.polarity in DOWNS:
                plot_kwargs = POLARITY_PROPS["down"]
            else:
                plot_kwargs = POLARITY_PROPS["unknown"]
            # Ensure boundedness
            baz = baz % 360
            ax.rake(baz, toa, 90, **plot_kwargs)
            if label:
                lon, lat = mplstereonet.rake(baz, toa, 90)
                ax.text(lon[0], lat[0], polarity.station)
            plot_kwargs.pop("label", None)  # Only label once
        if show:
            plt.show(block=True)
        return ax

    def tweak_np2(self, increment: float = 0.5, show: bool = True):
        from matplotlib.widgets import Slider
        np2_color = "turquoise"
        np2_outline = "orange"

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="stereonet")
        self.plot_polarities(ax=ax, show=False)
        self.nodal_plane_1.plot(ax=ax, show=False, label="NP1")
        self._artists = self.nodal_plane_2.plot(
            ax=ax, show=False, color=np2_color, label="NP2", 
            markeredgecolor=np2_outline)

        fig.subplots_adjust(bottom=0.25)
        fig.legend()
        ax.grid()
        # Add sliders
        x_start = 0.2
        # Strike
        axstrike = fig.add_axes([x_start, 0.2, 0.65, 0.03])
        strike_slider = Slider(
            axstrike, "NP-2 Strike", 0.0, 360.0, 
            valinit=self.nodal_plane_2.strike, valstep=increment,
            facecolor=np2_color)
        # Dip
        axdip = fig.add_axes([x_start, 0.15, 0.65, 0.03])
        dip_slider = Slider(
            axdip, "NP-2 Dip", 0.0, 90.0, 
            valinit=self.nodal_plane_2.dip, valstep=increment,
            facecolor=np2_color)
        # Rake
        axrake = fig.add_axes([x_start, 0.1, 0.65, 0.03])
        rake_slider = Slider(
            axrake, "NP-2 Rake", -180.0, 180.0, 
            valinit=self.nodal_plane_2.rake, valstep=increment,
            facecolor=np2_color)

        def update(val):
            self.nodal_plane_2.strike = strike_slider.val
            self.nodal_plane_2.dip = dip_slider.val
            self.nodal_plane_2.rake = rake_slider.val
            for artist in self._artists:
                artist.remove()
            self._artists = self.nodal_plane_2.plot(
                ax=ax, color=np2_color, markeredgecolor=np2_outline)
            fig.canvas.draw_idle()

        strike_slider.on_changed(update)
        dip_slider.on_changed(update)
        rake_slider.on_changed(update)

        if show:
            plt.show(block=True)
        return fig

    def find_planes(self, increment: float = 0.5, show: bool = True):
        from matplotlib.widgets import Slider
        np2_color = "turquoise"
        np2_outline = "orange"

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="stereonet")
        self.plot_polarities(ax=ax, show=False)
        self._np1_artists = self.nodal_plane_1.plot(
            ax=ax, show=False, label="NP1")
        self._np2_artists = self.nodal_plane_2.plot(
            ax=ax, show=False, color=np2_color, label="NP2", 
            markeredgecolor=np2_outline)

        fig.subplots_adjust(bottom=0.35)
        fig.legend()
        ax.grid()
        # Add sliders
        x_start = 0.2
        # -------------------- NP1 --------------------

        # Strike
        np1axstrike = fig.add_axes([x_start, 0.26, 0.65, 0.03])
        np1strike_slider = Slider(
            np1axstrike, "NP-1 Strike", 0.0, 360.0, 
            valinit=self.nodal_plane_1.strike, valstep=increment,
            facecolor="k")
        # Dip
        np1axdip = fig.add_axes([x_start, 0.21, 0.65, 0.03])
        np1dip_slider = Slider(
            np1axdip, "NP-1 Dip", 0.0, 90.0, 
            valinit=self.nodal_plane_1.dip, valstep=increment,
            facecolor="k")
        # Rake
        np1axrake = fig.add_axes([x_start, 0.16, 0.65, 0.03])
        np1rake_slider = Slider(
            np1axrake, "NP-1 Rake", -180.0, 180.0, 
            valinit=self.nodal_plane_1.rake, valstep=increment,
            facecolor="k")

        def update_np1(val):
            self.nodal_plane_1.strike = np1strike_slider.val
            self.nodal_plane_1.dip = np1dip_slider.val
            self.nodal_plane_1.rake = np1rake_slider.val
            for artist in self._np1_artists:
                artist.remove()
            self._np1_artists = self.nodal_plane_1.plot(ax=ax)
            fig.canvas.draw_idle()

        np1strike_slider.on_changed(update_np1)
        np1dip_slider.on_changed(update_np1)
        np1rake_slider.on_changed(update_np1)
        # -------------------- NP2 --------------------

        # Strike
        np2axstrike = fig.add_axes([x_start, 0.11, 0.65, 0.03])
        np2strike_slider = Slider(
            np2axstrike, "NP-2 Strike", 0.0, 360.0, 
            valinit=self.nodal_plane_2.strike, valstep=increment,
            facecolor=np2_color)
        # Dip
        np2axdip = fig.add_axes([x_start, 0.06, 0.65, 0.03])
        np2dip_slider = Slider(
            np2axdip, "NP-2 Dip", 0.0, 90.0, 
            valinit=self.nodal_plane_2.dip, valstep=increment,
            facecolor=np2_color)
        # Rake
        np2axrake = fig.add_axes([x_start, 0.01, 0.65, 0.03])
        np2rake_slider = Slider(
            np2axrake, "NP-2 Rake", -180.0, 180.0, 
            valinit=self.nodal_plane_2.rake, valstep=increment,
            facecolor=np2_color)

        def update_np2(val):
            self.nodal_plane_2.strike = np2strike_slider.val
            self.nodal_plane_2.dip = np2dip_slider.val
            self.nodal_plane_2.rake = np2rake_slider.val
            for artist in self._np2_artists:
                artist.remove()
            self._np2_artists = self.nodal_plane_2.plot(
                ax=ax, color=np2_color, markeredgecolor=np2_outline)
            fig.canvas.draw_idle()

        np2strike_slider.on_changed(update_np2)
        np2dip_slider.on_changed(update_np2)
        np2rake_slider.on_changed(update_np2)

        if show:
            plt.show(block=True)
        return fig

"""
TODO make sliders for strike, dip and rake of conjugate plane so for jupyter.
The idea is to have:
- another version where they input the azimuth, toa and polarity into a table,
  then it plots those data and provides two sets of s/d/r sliders and they have to work out the nodal planes.
"""

if __name__ == "__main__":
    np1 = NodalPlane(0, 30, 120)
    np2 = NodalPlane(146, 64, 73)
    polarities = [
        Polarity(90, 45, "positive"),
        Polarity(50, 60, "up"),
        Polarity(40, 75, "down"),
        Polarity(200, 45, "positive"),
        Polarity(250, 55, "down"),
    ]

    fm = FocalMechanism(
        nodal_plane_1=np1, nodal_plane_2=np2, polarities=polarities)
    fm.find_planes()

        