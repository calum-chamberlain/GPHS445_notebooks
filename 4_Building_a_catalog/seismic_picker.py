"""
Simple seismic phase picking applet for the GPHS445 course at VUW

:author: Calum Chamberlain
:date: 7 May 2019

License: LGPL v.3
"""

from matplotlib.dates import num2date, date2num
from obspy import UTCDateTime
import numpy as np


def seismic_picker(st, event_in=None):
    """
    A simple GUI for picking seismic phase arrivals.

    :type st: obspy.core.Stream
    :param st: The stream to pick within
    :type event_in: obspy.core.event.Event
    :param event_in: Optional event to edit.

    :returns: obspy.core.event
    """
    import getpass
    from obspy import UTCDateTime
    from obspy.core.event import (
        Event, Pick, WaveformStreamID, CreationInfo, Amplitude)

    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib.lines import Line2D
    from matplotlib import style

    style.use("ggplot")
    st = st.merge()
    st.sort(
        ["network", "station", "location", "channel", "starttime"])
    plot_start = min([tr.stats.starttime for tr in st]).datetime
    plot_end = max([tr.stats.endtime for tr in st]).datetime
    n_traces = len(st)
    fig, axes = plt.subplots(
        figsize=(8, 8), dpi=100, sharex=True, nrows=n_traces, ncols=1,
        squeeze=False)
    p_picks, s_picks, amplitude_picks, duration_picks = ({}, {}, {}, {})
    for row, tr in zip(axes, st):
        # Get linked pick
        p_pick_time, s_pick_time, amplitude, duration, polarity = (
            None, None, None, None, None)
        if event_in is not None:
            linked_picks = [
                p for p in event_in.picks
                if p.waveform_id.get_seed_string() == tr.id]
            # Get linked P phases
            linked_p_picks = [p for p in linked_picks if p.phase_hint == "P"]
            if len(linked_p_picks) > 1:
                print(
                    "Found {0} P picks for {1}, using the first at {2}".format(
                        len(linked_p_picks), tr.id, linked_p_picks[0].time))
            if len(linked_p_picks) > 0:
                p_pick_time = linked_p_picks[0].time
                try:
                    polarity = linked_p_picks[0].polarity
                    if polarity == "positive":
                        polarity = "up"
                    elif polarity == "negative":
                        polarity = "down"
                    else:
                        polarity = None
                except AttributeError:
                    polarity = None
            # Get linked S phases
            linked_s_picks = [p for p in linked_picks if p.phase_hint == "S"]
            if len(linked_s_picks) > 1:
                print(
                    "Found {0} S picks for {1}, using the first at {2}".format(
                        len(linked_s_picks), tr.id, linked_s_picks[0].time))
            if len(linked_s_picks) > 0:
                s_pick_time = linked_s_picks[0].time
            # Get linked amplitudes
            linked_amplitudes = [
                a for a in event_in.amplitudes
                if a.waveform_id.get_seed_string() == tr.id]
            linked_amplitude_points = [
                a for a in linked_amplitudes if a.type != 'END']
            if len(linked_amplitude_points) > 1:
                print("Found multiple amplitude picks, using the first")
            if len(linked_amplitude_points) > 0:
                amplitude = {
                    "amplitude": linked_amplitude_points[0].generic_amplitude,
                    "time": linked_amplitude_points[0].
                    pick_id.get_referred_object().time}
            # Get durations
            linked_durations = [
                a for a in linked_amplitudes if a.type == 'END']
            if len(linked_durations) > 1:
                print("Found multiple durations, using the first")
            if len(linked_durations) > 0:
                duration = {
                    "duration": linked_durations[0].generic_amplitude,
                    "time": linked_durations[0].pick_id.
                    get_referred_object().time}
        x = np.arange(0, tr.stats.npts)
        x = x * tr.stats.delta
        x = [(tr.stats.starttime + _).datetime for _ in x]
        for col in row:
            col.plot(x, tr.data, color="k", linewidth=0.5)
            # col.set_title(tr.id, rotation=0)
            col.text(0.02, 0.95, tr.id, transform=col.transAxes,
                     fontdict=dict(fontsize="small", ha="left", va="top"),
                     bbox=dict(boxstyle="round", fc="w", alpha=0.8))
            col.format_xdata = mdates.DateFormatter('%Y-%m-%dT%H:%M:%S.%f')
            col.set_xlabel("Time")
            col.set_xlim([plot_start, plot_end])
            if p_pick_time is not None:
                p_pick_line = col.add_line(
                    Line2D(xdata=[p_pick_time.datetime, p_pick_time.datetime],
                           ydata=list(col.get_ylim()), color='r'))
            else:
                p_pick_line = col.add_line(
                    Line2D(xdata=[], ydata=[], color="r"))
            if s_pick_time is not None:
                s_pick_line = col.add_line(
                    Line2D(xdata=[s_pick_time.datetime, s_pick_time.datetime],
                           ydata=list(col.get_ylim()), color="b"))
            else:
                s_pick_line = col.add_line(
                    Line2D(xdata=[], ydata=[], color="b"))
            if amplitude is not None:
                amplitude_pick = col.add_line(
                    Line2D(xdata=[amplitude["time"].datetime],
                           ydata=[amplitude["amplitude"]],
                           marker="o", markerfacecolor="r"))
            else:
                amplitude_pick = col.add_line(
                    Line2D(xdata=[], ydata=[], marker="+",
                           markerfacecolor="r", alpha=0.5))
            if duration is not None:
                duration_pick = col.add_line(
                    Line2D(xdata=[duration["time"].datetime,
                                  duration["time"].datetime],
                           ydata=list(col.get_ylim()),
                           color="k", linestyle="--"))
            else:
                duration_pick = col.add_line(
                    Line2D(xdata=[], ydata=[], color="k", linestyle="--"))

            p_picks.update(
                {tr.id: Picker(p_pick_line, button=1, polarity=polarity,
                               allow_polarity=True)})
            s_picks.update(
                {tr.id: Picker(s_pick_line, button=3)})
            amplitude_picks.update(
                {tr.id: Picker(amplitude_pick, button='a')})
            duration_picks.update(
                {tr.id: Picker(duration_pick, button='e')})
    print(
        "Make your picks using:\n"
        "\tleft mouse button: P\n\tright mouse button: S\n\t'a': amplitude"
        " at mouse location\n\t'e': end duration at mouse location - needs a "
        "P-pick to calculate duration\n\tP-wave polarity can be picked using "
        "the up and down arrows while hovering over a P-pick.\n\n"
        "Picks can be deleted by hovering over them and pressing"
        " the middle mouse button")
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.canvas.draw()
    plt.show()
    event_out = Event()
    for trace_id, picker in p_picks.items():
        if picker.time is not None:
            if picker.polarity == "up":
                polarity = "positive"
            elif picker.polarity == "down":
                polarity = "negative"
            else:
                polarity = "undecidable"
            event_out.picks.append(Pick(
                phase_hint="P", time=picker.time,
                waveform_id=WaveformStreamID(seed_string=trace_id),
                evaluation_mode="manual", polarity=polarity,
                creation_info=CreationInfo(author=getpass.getuser())))
    for trace_id, picker in s_picks.items():
        if picker.time is not None:
            event_out.picks.append(Pick(
                phase_hint="S", time=picker.time,
                waveform_id=WaveformStreamID(seed_string=trace_id),
                evaluation_mode="manual",
                creation_info=CreationInfo(author=getpass.getuser())))
    for trace_id, picker in amplitude_picks.items():
        if picker.time is not None:
            amp_pick = Pick(
                time=picker.time, phase_hint="IAML",
                waveform_id=WaveformStreamID(seed_string=trace_id),
                evaluation_mode="manual",
                creation_info=CreationInfo(author=getpass.getuser()))
            event_out.picks.append(amp_pick)
            event_out.amplitudes.append(Amplitude(
                generic_amplitude=picker.amplitude, type="A",
                pick_id=amp_pick.resource_id,
                waveform_id=WaveformStreamID(seed_string=trace_id),
                evaluation_mode="manual",
                creation_info=CreationInfo(author=getpass.getuser())))
    for trace_id, picker in duration_picks.items():
        if picker.time is None:
            continue
        # Get linked P pick
        duration_start = [
            p for p in event_out.picks 
            if p.phase_hint == "P" and
            p.waveform_id.get_seed_string() == trace_id]
        if len(duration_start) == 0:
            print("No matching P for duration on {0}".format(trace_id))
            continue
        duration_start = sorted([p.time for p in duration_start])[0]
        print("Duration: {0:.4f}s".format(picker.time - duration_start))
        duration_pick = Pick(
            time=picker.time, phase_hint="END",
            waveform_id=WaveformStreamID(seed_string=trace_id),
            evaluation_mode="manual",
            creation_info=CreationInfo(author=getpass.getuser()))
        event_out.picks.append(duration_pick)
        event_out.amplitudes.append(Amplitude(
            generic_amplitude=picker.time - duration_start, type="END",
            pick_id=duration_pick.resource_id,
            waveform_id=WaveformStreamID(seed_string=trace_id),
            evaluation_mode="manual",
            creation_info=CreationInfo(author=getpass.getuser())))
    return event_out


class Picker:
    def __init__(self, line, button=1, delete_threshold=0.2,
                 allow_polarity=False, polarity=None):
        if button == 2:
            raise IOError("Middle mouse button reserved for pick deletion")
        self.line = line
        self.button = button
        self.delete_threshold = delete_threshold
        self.allow_polarity = allow_polarity
        self.polarity = polarity
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        if len(self.ys) == 1:
            self.amplitude = self.ys[0]
        else:
            self.amplitude = None
        if len(self.xs) > 0:
            self.time = UTCDateTime(self.xs[0])
        else:
            self.time = None

        self.press = False
        self.move = False
        self.c1 = line.figure.canvas.mpl_connect(
            'button_press_event', self.onpress)
        self.c2 = line.figure.canvas.mpl_connect(
            'button_release_event', self.onrelease)
        self.c3 = line.figure.canvas.mpl_connect(
            'motion_notify_event', self.onmove)
        self.amp = line.figure.canvas.mpl_connect(
            'key_press_event', self.on_key)

    def on_key(self, event):
        """Register a key press"""
        if event.inaxes != self.line.axes:
            return
        if event.key == self.button:
            if self.button == "a":
                self.xs = [num2date(event.xdata)]
                self.ys = [event.ydata]
                self.time = UTCDateTime(num2date(event.xdata))
                self.amplitude = event.ydata
                print("Amplitude pick made at {0} with amplitude "
                      "{1:.2g}".format(self.time, self.amplitude))
            elif self.button == "e":
                self.xs = [num2date(event.xdata), num2date(event.xdata)]
                self.ys = list(self.line.axes.get_ylim())
                self.time = UTCDateTime(num2date(event.xdata))
                print("Duration end pick made at {0}".format(self.time))
            else:
                print("I only know what to do with a and e")
        # elif event.key == 2:
        #     # Delete the pick
        #     if self.time is not None:
        #         diff = abs((self.xs[0] - num2date(event.xdata)).total_seconds())
        #         if diff < self.delete_threshold:
        #             self.xs = []
        #             self.ys = []
        #             print("Deleted pick at time {0}".format(self.time))
        #             self.time = None
        #             self.line.set_data(self.xs, self.ys)
        #             self.line.axes.draw_artist(self.line)
        #             self.line.figure.canvas.draw()
        #             return
        elif event.key in ["up", "down"] and self.allow_polarity:
            if self.time is not None:
                diff = abs((self.xs[0] - num2date(event.xdata)).total_seconds())
                if diff < self.delete_threshold:
                    print("Polarity {0} recorded".format(event.key))
                    self.polarity = event.key
        else:
            return
        self.line.set_data(self.xs, self.ys)
        self.line.axes.draw_artist(self.line)
        self.line.figure.canvas.blit(self.line.axes.bbox)

    def onclick(self, event):
        """Make the pick"""
        if event.inaxes != self.line.axes:
            return
        if event.button == self.button:
            self.xs = [num2date(event.xdata), num2date(event.xdata)]
            self.ys = list(self.line.axes.get_ylim())
            self.time = UTCDateTime(num2date(event.xdata))
            print("Pick made at {0}".format(self.time))
        elif event.button == 2:
            # Delete the pick
            if self.time is not None:
                diff = abs((self.xs[0] - num2date(event.xdata)).total_seconds())
                if diff < self.delete_threshold:
                    self.xs = []
                    self.ys = []
                    print("Deleted pick at time {0}".format(self.time))
                    self.time = None
                    self.line.set_data(self.xs, self.ys)
                    self.line.axes.draw_artist(self.line)
                    self.line.figure.canvas.draw()
                    return
        else:
            return
        self.line.set_data(self.xs, self.ys)
        self.line.axes.draw_artist(self.line)
        self.line.figure.canvas.blit(self.line.axes.bbox)

    def onpress(self, event):
        self.press = True

    def onmove(self, event):
        if self.press:
            self.move = True

    def onrelease(self, event):
        if self.press and not self.move:
            self.onclick(event)
        self.press = False
        self.move = False

    def __repr__(self):
        return "Picker object set to time {0} with amplitude {1}".format(
            self.time, self.amplitude)


if __name__ == "__main__":
    import argparse
    import glob
    from obspy import read, read_events

    parser = argparse.ArgumentParser(description="Seismic Picker")
    parser.add_argument(
        "-f", "--file", type=str, required=False,
        help="The file you want to read in and pick. Supports wildcards")
    parser.add_argument(
        "-i", "--infile", type=str, required=False,
        help="A previous event file to plot the picks from")
    parser.add_argument(
        "-o", "--outfile", type=str, required=False,
        default="picked_event.xml", help="File to save picks to")
    parser.add_argument(
        "-l", "--lowcut", type=float, required=False,
        help="Lowcut filter to apply before picking")
    parser.add_argument(
        "-t", "--highcut", type=float, required=False,
        help="Highcut filter to apply before picking")

    args = vars(parser.parse_args())
    if args['file'] is None:
        st = read()
    else:
        st = read(args['file'])
    if args["infile"] is not None:
        cat = read_events(args["infile"])
        if len(cat) > 1:
            print("{0} events in catalog, using the first".format(len(cat)))
        event_in = cat[0]
    else:
        event_in = None
    
    if args['highcut'] is not None and args['lowcut'] is not None:
        st = st.filter(
            "bandpass", freqmin=args['lowcut'], freqmax=args['highcut'],
            corners=4)
    elif args['highcut'] is not None:
        st = st.filter("lowpass", freq=args["highcut"], corners=4)
    elif args['lowcut'] is not None:
        st = st.filter("highpass", freq=args["lowcut"], corners=4)

    event_out = seismic_picker(st, event_in=event_in)
    if len(event_out.picks) > 0:
        print("Saving picks to {0}".format(args['outfile']))
        event_out.write(args["outfile"], format="QUAKEML")
    else:
        print("No picks made, no outfile")
