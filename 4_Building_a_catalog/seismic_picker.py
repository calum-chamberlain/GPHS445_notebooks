"""
Simple seismic phase picking applet for the GPHS445 course at VUW

:author: Calum Chamberlain
:date: 7 May 2019

License: LGPL v.3
"""

import numpy as np

from bokeh.io import show
from bokeh.layouts import column, row
from bokeh.plotting import figure
from bokeh.models import CustomJS, Div
from bokeh import events
from bokeh.models.formatters import DatetimeTickFormatter


def display_pick(div, attributes=[], style = 'float:left;clear:left;font_size=10pt'):
    "Build a suitable CustomJS to display the current event in the div model."
    return CustomJS(args=dict(div=div), code="""
        var attrs = %s; var args = [];
        for (var i = 0; i<attrs.length; i++) {
            args.push(attrs[i] + '=' + Number(cb_obj[attrs[i]]).toFixed(2));
        }
        var line = "<span style=%r><b>Pick: </b>(" + args.join(", ") + ")</span>\\r";
        var text = div.text.concat(line);
        var lines = text.split("\\r")
        if (lines.length > 35)
            lines.shift();
        div.text = lines.join("\\r");
    """ % (attributes, style))

# Needs datetime and plot location of last pick - TapTool?

def seismic_picker(st, plot_width=1000, plot_height=200):
    point_attributes = ['x'] 
    plots = []
    for i, tr in enumerate(st):
        x = np.arange(0, tr.stats.npts)
        x = x * tr.stats.delta
        x = [(tr.stats.starttime + _).datetime for _ in x]
        if i == 0:
            p = figure(
                plot_width=plot_width, plot_height=plot_height, title=tr.id,
                tools="pan, wheel_zoom")
            p1 = p
        else:
            p = figure(
                plot_width=plot_width, plot_height=plot_height, title=tr.id,
                tools="pan, wheel_zoom", x_range=p1.x_range)
        p.line(x, tr.data, color='black')
        div = Div(width=400, height=p.plot_height)
        datetick_formatter = DatetimeTickFormatter(
            days=["%m/%d"], months=["%m/%d"],
            hours=["%m/%d %H:%M:%S"], minutes=["%m/%d %H:%M:%S"], 
            seconds=["%m/%d %H:%M:%S"], hourmin=["%m/%d %H:%M:%S"],
            minsec=["%m/%d %H:%M:%S"])
        p.xaxis.formatter = datetick_formatter
        p.js_on_event(
            events.Tap, display_pick(div, attributes=point_attributes))
        plots.append(row(p, div))
    show(column(plots))