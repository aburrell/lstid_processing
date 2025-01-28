#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
# ----------------------------------------------------------------------------
"""Functions that support diagnostic plotting."""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from importlib import resources


def close_figures(figs):
    """Close a list of figures, freeing up memory

    Parameters
    ----------
    figs : list-like
        List of figure handles

    """
    for fig in figs:
        plt.close(fig)

    return


def nrl_colors(num):
    """Determine color codes for NRL in a standard order

    Parameters
    ----------
    num : int
        Index number, will cycle to select an index between 0 and 10

    Returns
    -------
    color : str
        HEX string or named color (for unofficial colors)

    """
    # There are 11 defined colors, adjust the index to the appropriate range
    cnum = num

    while(cnum < 0):
        cnum += 10

    while(cnum > 10):
        cnum -= 10

    # Find the HEX code for the adjusted index
    if cnum == 0:
        color = '#162956'  # Navy
    elif cnum == 1:
        color = '#6DA9CF'  # Light Blue
    elif cnum == 2:
        color = '#2E5590'  # Blue
    elif cnum == 3:
        color = '#4C4F59'  # Gray
    elif cnum == 4:
        color = '#FAB208'  # Yellow
    elif cnum == 5:
        color = 'lightsteelblue'  # Light Steel Blue
    elif cnum == 6:
        color = 'orange'  # Orange
    elif cnum == 7:
        color = 'lightgray'  # Light Gray
    elif cnum == 8:
        color = 'k'  # Black
    elif cnum == 9:
        color = 'gold'  # Gold
    elif cnum == 10:
        color = 'lightskyblue'  # Light sky blue

    return color


def get_marker(num, small=False, line=None):
    """Select a matplotlib standard plotting marker in a specified order.

    Parameters
    ----------
    num : int
        Index number, starting at zero
    small : bool
        Use the small markers (default=False)
    line : bool or NoneType
        Replace marker with line (True), add marker to line (False), or
        do not include any line (None) (default=None)

    Returns
    -------
    marker : str
        String denoting a marker to be used in matplotlib formatting

    """

    big_markers = ["o", "s", "^", "v", "*", "X", "+", "d", "p", "<", ">", "8",
                   "D", "P", "x"]
    small_markers = [".", ",", "|", "_"]
    line_markers = ["-", "--", "-.", ":"]

    if line:
        mark = line_markers
    elif small:
        if line is None:
            mark = small_markers
        else:
            mark = ["{:s}{:s}".format(mm, line_markers[-1])
                    for mm in small_markers]
    else:
        if line is None:
            mark = big_markers
        else:
            mark = ["{:s}{:s}".format(mm, line_markers[0])
                    for mm in big_markers]

    i = num % len(mark)

    return mark[i]


def add_colorbar(figure_handle, contour_handle, zmin, zmax, zinc=6, name=None,
                 units=None, orient='vertical', scale='linear', width=1,
                 loc=[0.9, 0.1, 0.03, 0.8], overflow='neither', ticks=None,
                 ax=None, pad=.05):
    """Create a colorbar and add it to a figure.

    Parameters
    ----------
    figure_handle : matplotlib.figure.Figure
        handle to the figure
    contour_handle : matplotlib.collections.PathCollection
        handle to the plotted contour collection
    zmin : float
        minimum z value
    zmax : float
        maximum z value
    zinc : int
        z tick incriment (default=6)
    name : str or NoneType
        z variable name (default=None)
    units : str or NoneType
        z variable units (default=None)
    orient : str
        bar orientation, see mpl.colorbar (default='vertical')
    scale : str
        Tick scaling type, one of 'default', 'linear', 'datetime',
        'scientific', 'log', or 'exponential' (default='linear')
    width : float
        fraction of width (0-1), (default=1)
    loc : list-like
        location of colorbar in figure coordinates, see also mpl.colorbar
        (default=[0.95, 0.1, 0.03, 0.8])
    overflow : str
        include triangles showing overflow colors? Options include 'neither',
        'both', more? (default="neither")
    ticks : list-like or NoneType
        List of custom tick values, only used if value is not None
        (default=None)
    ax : matplotlib.axes._subplots.AxesSubplot or NoneType
        Axis to which the colorbar will be attached or None to create
        a unique axis for the colorbar (default=None)
    pad : float
        If `ax` is not None, fraction of the original Axes between the colorbar
        and new image Axes (default=0.05)

    Returns
    -------
    ax2 : matplotlib.axes._subplots.AxesSubplot
        Handle to the colorbar axis
    cb : matplotlib.colorbar.Colorbar
        Handle to the colorbar

    """
    # Set the colorbar options
    cb_kwargs = {'orientation': orient, 'extend': overflow}

    # Set the z range and output the colorbar
    if scale.lower() in ["linear", 'datetime', 'scientific']:
        cb_kwargs['ticks'] = np.linspace(zmin, zmax, zinc, endpoint=True)
    elif scale.lower() != 'default':
        cb_kwargs['ticks'] = np.logspace(np.log10(zmin), np.log10(zmax), zinc,
                                         endpoint=True)

    cax = None
    if ax is None:
        try:
            cax = figure_handle.add_axes(loc)
        except Exception:
            ax = figure_handle.axes[-1]

    # Change the z scale, if necessary
    if scale.lower() in ["log", "exponential"]:
        cb_kwargs['format'] = mpl.ticker.FormatStrFormatter('%7.2E')
    elif scale.lower() in ['scientific']:
        cb_kwargs['format'] = mpl.ticker.ScalarFormatter(useMathText=True)
        cb_kwargs['format'].set_powerlimits((0, 0))

    # Set up the colorbar
    if cax is None:
        cb = figure_handle.colorbar(contour_handle, ax=ax, pad=pad,
                                    **cb_kwargs)
    else:
        cb = figure_handle.colorbar(contour_handle, cax=cax, **cb_kwargs)

    # See if custom tick labels are required
    custom_ticks = False

    # See if the data is custom, temporal, or if the limits are multiples of pi
    if ticks is not None:
        cb_kwargs['ticks'] = ticks
        custom_ticks = True
    elif scale.lower() == 'datetime':
        # Requires max and min in datetime epoch days
        days_inc = (zmax - zmin) / zinc
        if days_inc > 365.0:
            ff = "%j/%Y"
            if units is None:
                units = "DDD/YYYY"
        elif days_inc > 1.0:
            ff = "%j"
            if units is None:
                units = "DDD"
        elif days_inc > 0.0007:
            ff = "%H:%M"
            if units is None:
                units = "HH:MM"
        else:
            ff = "%H:%M:%S"
            if units is None:
                units = "HH:MM:SS"

        cb_kwargs['ticks'] = [
            "{:s}".format(mpl.dates.num2date(wval).strftime(ff))
            for wval in cb_kwargs['ticks']]
        custom_ticks = True
    elif zmin % np.pi == 0 and zmax % np.pi == 0:
        wfac = cb_kwargs['ticks'] / np.pi
        cb_kwargs['ticks'] = list(cb_kwargs['ticks'])
        for i, wval in enumerate(wfac):
            if wval == 0.0:
                cb_kwargs['ticks'][i] = "{:.0f}".format(wval)
            elif wval == 1.0:
                cb_kwargs['ticks'][i] = r"$\pi$"
            elif wval == -1.0:
                cb_kwargs['ticks'][i] = r"-$\pi$"
            elif wval == int(wval):
                cb_kwargs['ticks'][i] = r"{:.0f}$\pi$".format(wval)
            else:
                cb_kwargs['ticks'][i] = r"{:.2f}$\pi$".format(wval)
        custom_ticks = True

    if custom_ticks:
        cb_kwargs['ticks'] = np.asarray(cb_kwargs['ticks'])
        if orient == "vertical":
            cb.ax.set_yticklabels(cb_kwargs['ticks'])
        else:
            cb.ax.set_xticklabels(cb_kwargs['ticks'])

    # Set the label and update the ticks
    if name is not None:
        if units is not None:
            cb.set_label(r'{:s} (${:s}$)'.format(name, units))
        else:
            cb.set_label(r'{:s}'.format(name))

    # Return the handle for the colorbar (which is treated as a subplot)
    # to allow for additional adjustments
    return cax, cb


def add_magnetic_equator(ax, meq_file=None, line_style='-', line_width=2,
                         line_color='k'):
    """Add a line showing the location of the magnetic equator.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        Subplot handle
    meq_file : str or NoneType
        Name of file containing location of magnetic equator with longitude in
        the first column and latitude in the second, or None to use local file
        (default=None)
    line_style : str
        Line style to use (default='-')
    line_width : int
        Line width to plot (default=2)
    line_color : str
        Line color to plot, may use named colors (default='k')

    Raises
    ------
    IOError
        If unable to find or open provided filename

    """
    # Test the magnetic equator filename
    if meq_file is None:
        meq_file = os.path.join(resources.files(__package__), 'map_files',
                                'igrf10.inc_equator')

    if not os.path.isfile(meq_file):
        raise IOError('unable to find magnetic equator file: {:}'.format(
            meq_file))

    # Load the magnetic equator data
    meq_data = np.loadtxt(meq_file, unpack=True)

    # Ensure it spans -180 to 360
    lon = list(meq_data[0])
    lat = list(meq_data[1])
    lat.extend(list(meq_data[1]))
    if meq_data[0].min() < 0.0:
        lon.extend(list(meq_data[0] + 360.0))
    else:
        lon.extend(list(meq_data[0] - 360.0))

    ilon = np.argsort(lon)
    lon = np.array(lon)[ilon]
    lat = np.array(lat)[ilon]

    # Plot the magnetic equator data
    ax.plot(lon, lat, line_style, color=line_color, lw=line_width)

    return
