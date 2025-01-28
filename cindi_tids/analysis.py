#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
# ----------------------------------------------------------------------------
"""TID analysis and plotting routines for CINDI data."""

import cmocean
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

import pysat
import pysatNASA

from cindi_tids import logger
from cindi_tids import plot_rout
from cindi_tids.smoothing import filter_rout


def init_cindi_tid_data(min_period=150, max_period=300):
    """Initialize a pysat object for loading CINDI TID data.

    Parameters
    ----------
    min_period : int or float
        Minimum period in seconds, -1 for high-pass filter (default=150)
    max_period : int or float
        Maximum period in seconds, -1 for low-pass filter (default=300)

    Returns
    -------
    cindi : pysat.Instrument
        pysat Instrument object with additional data parameters needed for
        TID analysis

    See Also
    --------
    Analysis.Smooth_Data.filter_rout.pysat_rel_data_butter

    Notes
    -----
    Period bandpass was determined using the orbital period and the wavelength
    of LSTIDs

    """
    # Initialize the CINDI instrument with no cleaning (temporary)
    # and orbital information
    orbit_info = {'kind': 'longitude', 'index': 'glon'}

    try:
        cindi = pysat.Instrument("cnofs", "ivm", clean_level='none',
                                 orbit_info=orbit_info)
    except KeyError:
        cindi = pysat.Instrument(inst_module=pysatNASA.instruments.cnofs_ivm,
                                 clean_level='none', orbit_info=orbit_info)

    # Set the relative data calculations
    relative = {'ionDensity': True, 'ionVelparallel': False,
                'ionVelmeridional': False}
    for dkey in relative.keys():
        rkey = '{:s}_rel_butter_Tmin{:.0f}s_Tmax{:.0f}s'.format(
            dkey, min_period, max_period)
        cindi.custom_attach(filter_rout.pysat_rel_data_butter,
                            args=[dkey],
                            kwargs={'samp_period': 1, 'min_period': min_period,
                                    'max_period': max_period, 'out_key': rkey,
                                    'nan_val': 0.0, 'relative': relative[dkey]})

    return cindi


def rel_data_orbit_plot(cindi, rel_key, vel_keys=None, vel_rel_keys=None,
                        vel_labels=None, fig=None, figname=None, ni_lim=None,
                        rel_lim=None, vel_max=None, vrel_lim=None):
    """Plot a standard figure for CINDI relative ion density.

    Parameters
    ----------
    cindi : pysat.Instrument
        C/NOFS CINDI IVM data, with a single orbit loaded
    rel_key : str
        Data key pointing to the relative ion density
    vel_keys : list or NoneType
        List of data keys pointing to additonal data, meridional and
        field-aligned drifts are encouraged (default=None)
    vel_rel_keys : list or NoneType
        List of data keys pointing to additonal relative quantities, meridional
        and field-aligned drifts are encouraged (default=None)
    vel_labels : list or NoneType
        List of label strings for additonal relative quantities, of the same
        order as the vel_rel_keys and vel_keys input (default=None)
    fig : matplotlib.figure.Figure or NoneType
        Figure handle, or None to initialize locally (default=None)
    figname : str or NoneType
        Figure name for saving output or None to not save to file
        (default=None)
    ni_lim : tuple or NoneType
        y-axis limits for the ion density data or None to use defaults
        (default=None)
    rel_lim : tuple or NoneType
        y-axis limits for the relative ion density data or None to use
        defaults (default=None)
    vel_max : float
        y-axis limits (symmetric) for the velocity data or None to use defaults
        (default=None)
    vrel_lim : tuple or NoneType
        y-axis limits for the addtinal relative data or None to use defaults
        (default=None)

    Returns
    -------
    fig : matplotlib.figure.Figure or None
        Figure handle or None, if figure was closed

    """
    # Set the standard data keys and data sets
    ni_key = 'ionDensity'
    lon_key = 'glon'
    lat_key = 'glat'
    alt_key = 'altitude'
    slt_key = 'slt'
    ctime = mpl.dates.date2num(cindi.index)

    # Get the data limits if not provided
    if ni_lim is None:
        ni_lim = (np.nanmin(cindi[ni_key]), np.nanmax(cindi[ni_key]))

    if rel_lim is None:
        rel_lim = (np.nanmin(cindi[rel_key]), np.nanmax(cindi[rel_key]))

    if vel_keys is not None and vel_max is None:
        vel_max = max([np.nanmax(abs(cindi[vkey])) for vkey in vel_keys])

    if vel_rel_keys is not None and vrel_lim is None:
        vrel_lim = max([np.nanmax(abs(cindi[vkey])) for vkey in vel_rel_keys])

    # Initialize the figure, if one is not supplied
    if fig is None:
        fig = plt.figure(figsize=(8.0, 7.5))

    # Add the desired axes in two groups
    grid_orbit = fig.add_gridspec(2, 1, top=.9, bottom=.55)
    grid_ni = fig.add_gridspec(2, 1, top=.45, bottom=.1)

    ax_alt = fig.add_subplot(grid_orbit[0, 0])
    ax_lat = fig.add_subplot(grid_orbit[1, 0])
    ax_ni = fig.add_subplot(grid_ni[0, 0])
    ax_rel = fig.add_subplot(grid_ni[1, 0])

    # Add the labels to the axes
    ax_alt.set_ylabel('Altitude (km)')
    ax_lat.set_ylabel(r'Latitude ($^\circ$)')
    ax_lat.set_xlabel(r'Longitude ($^\circ$)')
    ax_ni.set_ylabel('N$_i$ (cm$^{-3}$)')
    ax_rel.set_ylabel('Filtered N$_i$ / N$_i$')
    ax_rel.set_xlabel('Universal Time (HH:MM)')

    if vel_keys is not None:
        ax_vel = ax_ni.twinx()
        ax_vel.set_ylabel('v$_i$ (m s$^{-1}$)')

    if vel_rel_keys is not None:
        ax_vrel = ax_rel.twinx()
        ax_vrel.set_ylabel('Filtered v$_i$')

    # Plot the data
    con_t = ax_alt.scatter(cindi[lon_key], cindi[alt_key], c=ctime,
                           marker='o', s=1, cmap=cmocean.cm.thermal,
                           vmin=ctime[0], vmax=ctime[-1])
    plot_rout.add_magnetic_equator(ax_lat)
    ax_lat.scatter(cindi[lon_key], cindi[lat_key], c=cindi[slt_key],
                   marker='o', s=1, cmap=cmocean.cm.phase, vmin=0, vmax=24,
                   zorder=10)
    con_slt = ax_ni.scatter(cindi.index, cindi[ni_key], c=cindi[slt_key],
                            marker='o', s=1, cmap=cmocean.cm.phase, vmin=0,
                            vmax=24)
    ax_rel.plot(cindi.index, np.zeros(shape=cindi.index.shape), 'k:', lw=2)
    ax_rel.scatter(cindi.index, cindi[rel_key], c=cindi[slt_key], marker='o',
                   s=1, cmap=cmocean.cm.phase, vmin=0, vmax=24, label="N$_i$")

    ncol = 1
    if vel_keys is not None:
        ax_vel.plot(cindi.index, np.zeros(shape=cindi.index.shape), 'k:', lw=2)
        ax_vel.set_ylim(-vel_max, vel_max)

    if vel_rel_keys is not None and vel_keys is not None:
        if vel_labels is None:
            vel_labels = vel_keys
        ncol += len(vel_labels)

        for i, vkey in enumerate(vel_rel_keys):
            ax_vel.plot(cindi.index, cindi[vel_keys[i]],
                        plot_rout.get_marker(i), ms=1,
                        color=plot_rout.nrl_colors(i))
            ax_vrel.plot(cindi.index, cindi[vkey], plot_rout.get_marker(i),
                         ms=1, color=plot_rout.nrl_colors(i),
                         label=vel_labels[i])
        ax_vrel.legend(loc=1, fontsize="medium", ncol=ncol, scatterpoints=3,
                       bbox_to_anchor=[1.15, -.15], markerscale=5)
        ax_vrel.set_ylim(vrel_lim)

    # Format the axes
    ax_alt.set_xlim(0, 360)
    ax_alt.set_ylim(350, 860)
    ax_alt.xaxis.set_major_locator(mpl.ticker.MultipleLocator(60))
    ax_alt.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(''))
    ax_lat.set_xlim(0, 360)
    ax_lat.set_ylim(-15, 15)
    ax_lat.xaxis.set_major_locator(mpl.ticker.MultipleLocator(60))
    ax_ni.set_yscale('log')
    ax_ni.xaxis.set_major_formatter(mpl.dates.DateFormatter(""))
    ax_ni.set_xlim(cindi.index[0], cindi.index[-1])
    ax_ni.set_ylim(ni_lim)
    ax_rel.xaxis.set_major_formatter(mpl.dates.DateFormatter("%H:%M"))
    ax_rel.set_xlim(cindi.index[0], cindi.index[-1])
    ax_rel.set_ylim(rel_lim)

    # Format the figure
    fig.subplots_adjust(hspace=0.0, right=0.85)
    fig.suptitle(
        '{:} {:} on {:}\n{:}'.format(
            cindi.platform.upper(), cindi.name.upper(),
            cindi.index[0].strftime('%d %b %Y'),
            cindi.meta[rel_key, cindi.meta.labels.desc].split('. ')[-1]),
        fontsize='medium')

    # Add colorbars
    plot_rout.add_colorbar(fig, con_t, ctime[0], ctime[-1], zinc=5, name="UT",
                           scale='datetime', loc=[0.88, 0.75, 0.01, 0.15])
    plot_rout.add_colorbar(fig, con_slt, 0, 24, zinc=7, name="SLT", units='h',
                           scale='linear', loc=[0.88, 0.55, 0.01, 0.15])

    # Save the figure, if output name supplied
    if figname is not None:
        if os.path.isfile(figname):
            logger.warning(''.join(['file [', figname, '] already exists, ',
                                    'returning figure handle instead of ',
                                    'writing over existing file']))
        else:
            fig.savefig(figname)
            plot_rout.close_figures([fig])
            fig = None

    return fig


def plot_cindi_tid_orbits(plot_dir, stime, etime=None, min_period=150,
                          max_period=300):
    """Create CINDI LSTID Ni variation orbit plots.

    Parameters
    ----------
    plot_dir : str or NoneType
        Plot directory or None to return all created figure handles
    stime : dt.datetime
        Start time
    etime : dt.datetime or NoneType
        End time or None if only 1 day is required (default=None)
    min_period : int or float
        Minimum period in seconds, -1 for high-pass filter (default=150)
    max_period : int or float
        Maximum period in seconds, -1 for low-pass filter (default=300)

    Returns
    -------
    figs : list
        List of figure handles (empty if all saved to file)

    Raises
    ------
    IOError
        If plot directory doesn't exist

    """
    # Test the plot directory
    if plot_dir is not None and not os.path.isdir(plot_dir):
        logger.warn('specified plot directory does not exist: {:}'.format(
            plot_dir))
        plot_dir = None

    # Initialize the output
    figs = list()

    # Initialize the CINDI object and data limits
    cindi = init_cindi_tid_data(min_period=min_period, max_period=max_period)
    rkey = cindi.custom_kwargs[0]['out_key']
    rel_vel_keys = [cindi.custom_kwargs[i]['out_key'] for i in [1, 2]]
    vel_keys = [rvkey.split('_')[0] for rvkey in rel_vel_keys]
    vel_labels = ["v$_{:s}$".format("{||}" if vkey.find("parallel") > 0
                                    else "{mer}") for vkey in vel_keys]
    ni_lim = [1000, 4.0e6]
    rel_lim = [-0.3, 0.3]
    vrel_lim = [-30, 30]

    # Set the file name root
    if plot_dir is not None:
        figroot = "_".join([cindi.platform, cindi.name, rkey])

    # Load the data and test the ion density limits
    cindi.load(date=stime, end_date=etime)

    if np.nanmin(cindi['ionDensity']) < ni_lim[0]:
        ni_lim[0] = np.nanmin(cindi['ionDensity'])

    if np.nanmax(cindi['ionDensity']) > ni_lim[1]:
        ni_lim[1] = np.nanmax(cindi['ionDensity'])

    # Cycle by orbit
    if etime is None:
        etime = stime
    else:
        etime -= dt.timedelta(days=1)

    cindi.bounds = (stime, etime)
    cindi.orbits[0]
    for iorb in range(cindi.orbits.num):
        # Cycle to the next orbit, if needed
        if iorb > 0:
            cindi.orbits.next()

        # Set the output filename
        if plot_dir is None:
            figname = None
        else:
            figname = os.path.join(plot_dir, "{:s}_orbit{:03d}.png".format(
                figroot, cindi.orbits.current))

        # Create and save the figure for this orbit
        fig = rel_data_orbit_plot(cindi, rkey, vel_keys=vel_keys,
                                  vel_rel_keys=rel_vel_keys,
                                  vel_labels=vel_labels, figname=figname,
                                  ni_lim=ni_lim, rel_lim=rel_lim, vel_max=200,
                                  vrel_lim=vrel_lim)

        if fig is not None:
            figs.append(fig)

    return figs


def identify_tid(cindi, dens_pert_thresh=0.2, dens_quiet_thresh=0.1,
                 vel_pert_thresh=10.0, vel_sec=60, join_sec=300,
                 dens_var='ionDensity_rel_butter_Tmin150s_Tmax300s',
                 vmer_var='ionVelmeridional_rel_butter_Tmin150s_Tmax300s',
                 vpar_var='ionVelparallel_rel_butter_Tmin150s_Tmax300s'):
    """Identify periods of TID activity in CINDI IVM data.

    Parameters
    ----------
    cindi : pysat.Instrument
        CINDI IVM instrument object with perturbation data
    dens_pert_thresh : float
        Absolute value for plasma density threshold, above which wave activity
        is present (default=0.2)
    dens_quiet_thresh : float
        Absolute value for plasma density threshold, below which wave activity
        is absent (default=0.1)
    vel_pert_thresh : float
        Absolute value for plasma velocity threshold, above which wave activity
        is present (default=10.0)
    vel_sec : int
        Maximum number of seconds allowed between a perturbed density
        observation and an inquiet velocity observation (default=60)
    join_sec : int
        Maximum number of seconds between triggers to identify a cohesive event
        period (default=300)
    dens_var : str
        Perturbed plasma density variable name
        (default='ionDensity_rel_butter_Tmin150s_Tmax300s')
    vmer_var : str
        Perturbed meridional plasma drift variable name
        (default='ionVelmeridional_rel_butter_Tmin150s_Tmax300s')
    vpar_var : str
        Perturbed field-aligned plasma drift variable name
        (default='ionVelparallel_rel_butter_Tmin150s_Tmax300s')

    Returns
    -------
    event_start : list
        List of TID event start times
    event_end : list
        List of TID event end times, same length as `event_start`

    """

    # Test to ensure the necessary variable names are present
    if np.any([var not in cindi.variables for var in [dens_var, vmer_var,
                                                      vpar_var]]):
        raise KeyError('cindi Instrument missing required data variable')

    # Initalize the output data
    event_start = []
    event_end = []

    # Get the perturbed and active (not quiet) indexes
    ipert_dens = np.where(abs(cindi[dens_var]) > dens_pert_thresh)[0]
    iactive = np.where((abs(cindi[dens_var]) >= dens_quiet_thresh)
                       | (abs(cindi[vmer_var]) > vel_pert_thresh)
                       | (abs(cindi[vpar_var]) > vel_pert_thresh))[0]
    ipert_vel = np.where((abs(cindi[vmer_var]) > vel_pert_thresh)
                         | (abs(cindi[vpar_var]) > vel_pert_thresh))[0]

    # Cycle through each perturbed time, checking to see if the velocity is
    # sufficently active
    for ptime in cindi.index[ipert_dens]:
        if len(event_start) == 0 or ptime >= event_end[-1]:
            # Only evaluate times not currently encompassed by an event
            del_vel_sec = abs(cindi.index[ipert_vel]
                              - ptime).min().total_seconds()

            if del_vel_sec <= vel_sec:
                # This is not a plasma bubble, since velocities and densities
                # both have waves. Find the start of the event.
                del_dens_sec = np.array([abs(ctime - ptime).total_seconds()
                                         for ctime in cindi.index[iactive]])
                ievent = np.where(del_dens_sec < join_sec)[0]
                estart = cindi.index[iactive[ievent]].min()
                eend = cindi.index[iactive[ievent]].max()

                # If this is a new event boundary, it may need to be expanded
                expand_start = False
                expand_end = False

                # Update the event time boundaries
                if(len(event_end) > 0
                   and estart < event_end[-1] + dt.timedelta(seconds=join_sec)):
                    # This event starts within the join period of the previous
                    # event.  Combine by expanding the event time period.
                    if eend > event_end[-1]:
                        event_end[-1] = eend
                        expand_end = True
                    if estart < event_start[-1]:
                        event_start[-1] = estart
                        expand_start = True
                else:
                    # This is a new event, save the start and end times
                    event_start.append(estart)
                    event_end.append(eend)
                    expand_start = True
                    expand_end = True

                while expand_start:
                    # Test to see if there are more active times in range of
                    # the start of the event
                    del_dens_sec = np.array([
                        abs(ctime - event_start[-1]).total_seconds()
                        for ctime in cindi.index[iactive]])
                    ievent = np.where(del_dens_sec < join_sec)[0]

                    if len(ievent) > 0:
                        estart = cindi.index[iactive[ievent]].min()
                        if estart < event_start[-1]:
                            event_start[-1] = estart
                        else:
                            expand_start = False
                    else:
                        expand_start = False

                while expand_end:
                    # Test to see if there are more active times in range of
                    # the end of the event
                    del_dens_sec = np.array([
                        abs(ctime - event_end[-1]).total_seconds()
                        for ctime in cindi.index[iactive]])
                    ievent = np.where(del_dens_sec < join_sec)[0]

                    if len(ievent) > 0:
                        eend = cindi.index[iactive[ievent]].max()
                        if eend > event_end[-1]:
                            event_end[-1] = eend
                        else:
                            expand_end = False
                    else:
                        expand_end = False

    # Return the start and end times of the events
    return event_start, event_end
