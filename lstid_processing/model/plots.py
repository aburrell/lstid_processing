#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
#
# DISTRIBUTION STATEMENT A: Approved for public release. Distribution is
# unlimited.
# ----------------------------------------------------------------------------
"""Create some useful plots for the SAMI3 concatonated data."""

import cmocean
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pds

import lstid_processing.plot_rout as pr
from lstid_processing.model import analysis


def plot_f2_peak_char_and_diff(sami, f2_inds, nt=None, dat_keys=None,
                               start=dt.datetime(2014, 3, 25, 22),
                               stop=dt.datetime(2014, 3, 26, 6)):
    """Plot the F2 peak characteristics and hemispheric differences.

    Parameters
    ----------
    sami : xr.Dataset
        Dataset with concatonated SAMI3 data
    f2_inds : dict
        Dict with north and south indices for the F2 peaks along a field line
    nt : int or NoneType
        Event time index or None to not include a vertical line marking this
        point (default=None)
    dat_keys : list-like
        List of data keys to plot, if None uses altitude and field-aligned
        neutral wind (default=None)
    start : dt.datetime
        Starting time for plot (default=dt.datetime(2014, 3, 25, 22))
    stop : dt.datetime
        Ending time for plot (default=dt.datetime(2014, 3, 26, 6))

    Returns
    -------
    fig : plt.Figure
        Figure handle

    """
    # Set the defaults if needed
    if dat_keys is None:
        dat_keys = ['zalt', 'vnq']

    # Initialize the figure
    fig = plt.figure(figsize=(6.51, 7.68))  # Set in size for 2 data variables
    axes = [fig.add_subplot(len(dat_keys), 1, 1 + i)
            for i in range(len(dat_keys))]

    # Cycle through each data variable, plotting the North, South, and N-S
    stimes = sami['datetime'].values
    for i, ax in enumerate(axes):
        ax.plot(stimes, sami[dat_keys[i]].values[f2_inds['north']], '-',
                color=pr.nrl_colors(0), lw=2, label='North')
        ax.plot(stimes, sami[dat_keys[i]].values[['south']], '--',
                color=pr.nrl_colors(1), lw=2, label='South')
        ax.plot(stimes, sami[dat_keys[i]].values[f2_inds['north']]
                - sami[dat_keys[i]].values[f2_inds['south']], '-.',
                color=pr.nrl_colors(4), lw=2, label='North-South')

        # Format the axes
        ax.set_xlim(start, stop)
        ax.set_ylabel('{:s} ({:s})'.format(dat_keys[i],
                                           sami[dat_keys[i]].units))
        ax.grid()

        if i > len(axes) - 1:
            ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(''))
        else:
            ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%j\n%H:%M'))
            ax.set_xlabel('Universal Time')

        if nt is not None:
            ylim = ax.get_ylim()
            ax.plot([stimes[nt], stimes[nt]], ylim, 'k-', lw=2)

        if i == 0:
            ax.legend(fontsize='medium', ncols=3, loc=1,
                      bbox_to_anchor=(.94, 1.23))

    # Format the figure
    fig.subplots_adjust(top=.84, hspace=.05)
    fig.suptitle(
        'Field Line F$_2$ Peak Characteristics\n{:} - {:} UT'.format(
            start.strftime('%d %b %Y %H:%M'), stop.strftime('%d %b %Y %H:%M')),
        fontsize='medium')

    return fig


def plot_dom_acceleration(sami, nt, nlind, dom_acc=None, nfinds=None,
                          nz_sat=None, nz_north=None, nz_south=None):
    """Plot the dominate acceleration terms for a single time and meridian.

    Parameters
    ----------
    sami : xr.Dataset
        Dataset with concatonated SAMI3 data
    nt : int
        Time index
    nlind : int
        'nl' meridian index
    dom_acc : array-like or NoneType
        Array of values specifying the dominant acceleration term or None
        to calculate (default=None)
    nfinds : list-like or None
        List of 'nf' indices for field lines to be plotted (default=None)
    nz_sat : list-like or None
        List of 'nz' indices corresponding to the `nfinds` indices for
        satellite locations (default=None)
    nz_north : list-like or None
        List of 'nz' indices corresponding to the `nfinds` indices for northern
        F2 peak locations (default=None)
    nz_south : list-like or None
        List of 'nz' indices corresponding to the `nfinds` indices for southern
        F2 peak locations (default=None)

    Returns
    -------
    fig : plt.Figure
        Figure handle

    """
    # Get the acceleration terms, if needed
    if dom_acc is None:
        dom_acc = analysis.get_dominant_acceleration(sami)

    # Initialize the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Set the contour levels
    levels = np.linspace(0, 3, 4)

    # Plot the acceleration terms
    con = ax.contourf(sami['glat'][nt, nlind], sami['zalt'][nt, nlind],
                      dom_acc[nt], cmap=cmocean.cm.deep, levels=levels)

    # Format the plot
    ax.set_xlim(-15, 30)
    ax.set_ylim(100, 1000)
    ax.set_ylabel('Altitude (km)')
    ax.set_xlabel(r'Geodetic Latitude ($^\circ$)')

    pr.add_colorbar(fig, con, 0, 3, zinc=4, name='Strongest O$^+$ Variations',
                    loc=[0.88, 0.11, 0.02, 0.77])

    # Format the figure
    fig.subplots_adjust(right=.87)
    fig.suptitle('{:}'.format(sami['datetime'].values[nt].strftime(
        '%d %b %Y %H:%M UT')), fontsize='medium')

    # Add lines, if desired
    if nfinds is not None:
        for i, nfind in enumerate(nfinds):
            ax.plot(sami['glat'][nt, nlind, nfind],
                    sami['zalt'][nt, nlind, nfind], 'w--', lw=2)

            if nz_sat is not None:
                ax.plot(sami['glat'][nt, nlind, nfind, nz_sat[i]],
                        sami['zalt'][nt, nlind, nfind, nz_sat[i]], 'wo', ms=10,
                        mfc='none')

            if nz_north is not None:
                ax.plot(sami['glat'][nt, nlind, nfind, nz_north[i]],
                        sami['zalt'][nt, nlind, nfind, nz_north[i]], 'wX',
                        ms=10)

            if nz_south is not None:
                ax.plot(sami['glat'][nt, nlind, nfind, nz_south[i]],
                        sami['zalt'][nt, nlind, nfind, nz_south[i]], 'wX',
                        ms=10)

    return fig


def plot_dens_var(sami, nlind, nfind, nzind, sat_key, nt=None, nt_color='k',
                  title='', start=dt.datetime(2014, 3, 25, 22),
                  stop=dt.datetime(2014, 3, 26, 6)):
    """Plot the electron, O+, and H+ density variations at a desired location.

    Parameters
    ----------
    sami : xr.Dataset
        Dataset with concatonated SAMI3 data
    nlind : int
        'nl' index
    nfind : int
        'nf' index
    nzind : int
        'nf' index
    nt : int or None
        Time index to plot a vertical line marking time or None (default=None)
    nt_color : str
        Color for vertical line (default='k')
    sat_key : str
        Single character string used to specify the satellite/meridian for
        the data variations (e.g., 'c' or 'd')
    title ; str
         Figure title string (default='')
    start : dt.datetime
        Starting time for plot (default=dt.datetime(2014, 3, 25, 22))
    stop : dt.datetime
        Ending time for plot (default=dt.datetime(2014, 3, 26, 6))

    Returns
    -------
    fig : plt.Figure
        Figure handle

    """
    # Set the data keys
    dkeys = ['rel_dene_{:s}'.format(sat_key), 'rel_deni1_{:s}'.format(sat_key),
             'rel_deni2_{:s}'.format(sat_key)]
    labels = [r'$\Delta N_e/N_e$', r'$\Delta N_{H^+}/N_{H^+}$',
              r'$\Delta N_{O^+}/N_{O^+}$']
    percent = [100, 100 * sami['deni1'].values[nt, nlind, nlind, nzind]
               / sami['dene'].values[nt, nlind, nlind, nzind],
               100 * sami['deni2'].values[nt, nlind, nlind, nzind]
               / sami['dene'].values[nt, nlind, nlind, nzind]]

    # Initialize figure
    fig = plt.figure(figsize=(8.34, 6.75))
    axes = [fig.add_subplot(len(dkeys), 1, 1 + i) for i in range(len(dkeys))]

    # Cycle through each data type
    for i, ax in enumerate(axes):
        if nt is not None:
            ax.plot([sami['datetime'].values[nt], sami['datetime'].values[nt]],
                    [-3000, 3000], '-', color=nt_color, lw=2)

        # Plot the data
        ax.plot(sami['datetime'].values, sami[dkeys[i]][:, nzind], '--',
                color=pr.nrl_colors(0), lw=2)

        # Format the axis
        ax.set_ylim(-.05, .05)
        ax.set_xlim(start, stop)
        ax.grid()
        ax.set_ylabel(labels[i])

        if i < len(axes) - 1:
            ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(''))
        else:
            ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%j\n%H:%M'))
            ax.set_xlabel('Universal Time')

        if i == 0:
            ax.text(16154.85, 0.06, 'Percent Ion or Electron Composition')

        ax.text(16155.12, 0.3, '{:.1f}%'.format(percent[i]))

    # Set the title
    fig.suptitle(
        r'SAMI3 for {:d} at {:.1f}$^\circ$N, {:.1f}$^\circ$E, {:.1f} km'.format(
            start.year, sami['glat'].values[nt, nlind, nfind, nzind],
            sami['glon'].values[nt, nlind, nfind, nzind],
            sami['zalt'].values[nt, nlind, nfind, nzind]), fontsize='medium')

    return fig


def plot_field_lines_w_linear_var(sami, nt, nlind, nfinds, nz_sat, nz_north,
                                  nz_south, colors, sat_labels):
    """Plot the 6-panel conjunction figure.

    Parameters
    ----------
    sami : xr.Dataset
        SAMI3 standard or concatonated data set
    nt : int
        Time index
    nlind : int
        'nl' index
    nfinds : list-like
        List of two 'nl' indexes
    nz_sat : list-like
        List of 'nz' indices corresponding to the `nfinds` indices for
        satellite locations
    nz_north : list-like
        List of 'nz' indices corresponding to the `nfinds` indices for northern
        F2 peak locations
    nz_south : list-like or None
        List of 'nz' indices corresponding to the `nfinds` indices for southern
        F2 peak locations
    colors : list-like
        This that contains two colors for each satellite as a list, for example
        [['orange', 'tan'], ['k', 'grey']]
    sat_labels : list-like
        List of satellite labels

    Returns
    -------
    fig : plt.Figure
        Figure handle

    """
    # Make plot with many variables
    fig = plt.figure(figsize=(14.5, 4.8))
    axes = [fig.add_subplot(2, 3, 1 + i) for i in range(6)]
    twins = [None for axes in range(6)]

    # Cycle through each satellite
    for i, nfind in enumerate(nfinds):
        # Create the field-line plot with satellite and peak density locations.
        # Start by plotting the field line
        axes[0].plot(sami['glat'].values[nlind, nfind, :],
                     sami['zalt'].values[nlind, nfind, :], '-',
                     color=colors[i][0], label=sat_labels[i], lw=2)

        # Plot the F2 peaks on top of the field line
        axes[0].plot(sami['glat'].values[nlind, nfind, nz_north[i]],
                     sami['zalt'].values[nlind, nfind, nz_north[i]], 'X',
                     color=colors[i][1], ms=10)
        axes[0].plot(sami['glat'].values[nlind, nfind, nz_south[i]],
                     sami['zalt'].values[nlind, nfind, nz_south[i]], 'X',
                     color=colors[i][1], ms=10)

        # Plot the satellite location on top of the field line
        axes[0].plot(sami['glat'].values[nlind, nfind, nz_sat[i]],
                     sami['zalt'].values[nlind, nfind, nz_sat[i]], 'o',
                     color='w', ms=10, markeredgecolor=colors[i][0],
                     markeredgewidth=2, zorder=1)

        # Create the electron density plot
        axes[3].plot(sami['glat'].values[nlind, nfind, :],
                     sami['dene'].values[nt, nlind, nfind, :], '-',
                     color=colors[i][0], label=sat_labels[i], lw=2)

        # Plot the F2 peaks on top of the electron density
        axes[3].plot(sami['glat'].values[nlind, nfind, nz_north[i]],
                     sami['dene'].values[nt, nlind, nfind, nz_north[i]], 'X',
                     color=colors[i][1], ms=10)
        axes[3].plot(sami['glat'].values[nlind, nfind, nz_south[i]],
                     sami['dene'].values[nt, nlind, nfind, nz_south[i]], 'X',
                     color=colors[i][1], ms=10)

        # Plot the satellite location on top of the electron density
        axes[3].plot(sami['glat'].values[nlind, nfind, nz_sat[i]],
                     sami['dene'].values[nt, nlind, nfind, nz_sat[i]], 'o',
                     color='w', ms=10, markeredgecolor=colors[i][0],
                     markeredgewidth=2, zorder=1)

        # Create the DMSP time variation plot
        iax = 1 + i * 3
        axes[iax].plot(
            np.full(shape=(2,), fill_value=sami['hrut'].values[nt]),
            [-300, 300], 'k-', lw=2)
        axes[iax].plot(sami['hrut'].values,
                       sami['rel_vnq_d'].values[:, nz_sat[i]] / 100.0, '--',
                       color='blue', label='u$_{||}$', lw=2)
        axes[iax].plot(sami['hrut'].values,
                       sami['rel_vnp_d'].values[:, nz_sat[i]] / 100.0, ':',
                       color='cornflowerblue', label='u$_{exb_{mer}}$', lw=2)
        axes[iax].plot(sami['hrut'].values,
                       sami['rel_vsi2_d'].values[:, nz_sat[i]] / 100.0, '--',
                       color='green', label='v$_{||_{O^+}}$', lw=2)
        axes[iax].plot(sami['hrut'].values,
                       sami['rel_u1p_d'].values[:, nz_sat[i]] / 100.0, ':',
                       color='yellowgreen', label='v$_{exb_{mer}}$', lw=2)
        axes[iax].plot([sami['hrut'].values[nt], sami['hrut'].values[nt]],
                       [-200, 200], '-', color=colors[i][0], lw=2)

        # Set the common formatting
        axes[iax].set_ylabel('{:s}\n{:s}'.format(sat_labels[i],
                                                 'filtered vel (m s$^{-1}$)'))
        axes[iax].set_xlim(0, 4)
        axes[iax].yaxis.grid()

        # Plot along the twin
        if twins[iax] is None:
            twins[iax] = axes[iax].twinx()
            twins[iax].yaxis.label.set_color('darkviolet')
            twins[iax].tick_params(axis='y', colors='darkviolet')

        twins[iax].plot(sami['hrut'].values,
                        sami['rel_dene_d'].values[:, nz_sat[i]], '-',
                        color='darkviolet', lw=2, label='Electron Density')
        twins[iax].plot(sami['hrut'].values,
                        sami['rel_denn2_d'].values[:, nz_sat[i]], '-.',
                        color='darkviolet', lw=2, label='O Density')

        # Create the lat variation plot
        iax = 2 + i * 3
        axes[iax].plot(
            np.full(shape=2,
                    fill_value=sami['glat'].values[nlind, nfind, nz_north[i]]),
            [-500, 500], '-', color=colors[i][1], lw=2)
        axes[iax].plot(
            np.full(shape=2,
                    fill_value=sami['glat'].values[nlind, nfind, nz_south[i]]),
            [-500, 500], '-', color=colors[i][1], lw=2)
        axes[iax].plot(
            np.full(shape=2,
                    fill_value=sami['glat'].values[nlind, nfind, nz_sat[i]]),
            [-500, 500], '-', color=colors[i][0], lw=2)
        axes[iax].plot(sami['glat'].values[nlind, nfind, :],
                       sami['rel_vnq_d'].values[nt, :] / 100.0, '--',
                       color='blue', label='u$_{||}$', lw=2)
        axes[iax].plot(sami['glat'].values[nlind, nfind, :],
                       sami['rel_vnp_d'].values[nt, :] / 100.0, ':',
                       color='cornflowerblue', label='u$_{exb_{mer}}$', lw=2)
        axes[iax].plot(sami['glat'].values[nlind, nfind, :],
                       sami['rel_vsi2_d'].values[nt, :] / 100.0, '--',
                       color='green', label='v$_{||_{O^+}}$', lw=2)
        axes[iax].plot(sami['glat'].values[nlind, nfind, :],
                       sami['rel_u1p_d'].values[nt, :] / 100.0, ':',
                       color='yellowgreen', label='v$_{exb_{mer}}$', lw=2)

        # Intialize the twin
        if twins[iax] is None:
            twins[iax] = axes[iax].twinx()
            twins[iax].set_ylabel(r'$\Delta$ N / N')
            twins[iax].yaxis.label.set_color('darkviolet')
            twins[iax].tick_params(axis='y', colors='darkviolet')

        # Plot along the latitude twin
        twins[iax].plot(sami['glat'].values[nlind, nfind, :],
                        sami['rel_dene_d'].values[nt, :], '-',
                        color='darkviolet', lw=2, label='Electron Density')
        twins[iax].plot(sami['glat'].values[nlind, nfind, :],
                        sami['rel_denn2_d'].values[nt, :], '-.',
                        color='darkviolet', lw=2, label='O Density')

        # Set the common formatting
        axes[iax].set_xlim(
            sami['glat'].values[nlind, nfind, nz_south[i]] - 1.0,
            sami['glat'].values[nlind, nfind, nz_north[i]] + 1.0)
        axes[iax].yaxis.grid()

    # Format the field-line axis
    axes[0].set_ylim(0, 1000)
    axes[0].set_ylabel('Altitude (km)')
    axes[0].xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(''))

    # Format the electron density axis
    axes[3].ticklabel_format(style='sci', scilimits=(-3, 3), axis='y')
    axes[3].yaxis.major.formatter._useMathText = True
    axes[3].set_xlabel(r'Geographic Latitude ($^\circ$)')
    axes[3].set_ylabel('N$_e$ (cm$^{-3}$)')

    # Format the time plot axes
    axes[1].set_ylim(-40, 40)
    axes[1].xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(''))
    axes[4].set_ylim(-15, 15)
    axes[4].set_xlabel('Universal Time (h)')
    twins[1].set_ylim(-.2, .2)
    twins[4].set_ylim(-.03, .03)

    # Format the latitude plot axes
    axes[2].set_ylim(-20, 20)
    axes[5].set_ylim(-15, 15)
    twins[2].set_ylim(-.1, .1)
    twins[5].set_ylim(-.15, .15)
    axes[2].xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(''))
    axes[5].set_xlabel(r'Geographic Latitude ($^\circ$)')

    # Add the legend and adjust the spacing
    fig.subplots_adjust(left=.07, right=.93, wspace=.35)
    axes[0].legend(loc=1, ncol=2, bbox_to_anchor=[1.05, 1.36])
    axes[2].legend(ncol=4, loc=1, bbox_to_anchor=[0.0, 1.36])
    twins[2].legend(ncol=2, loc=1, bbox_to_anchor=[1.3, 1.36])

    return fig


def get_plot_tid_peaks(sami, nt_start, nt_stop, nlind, nfind, nzinds,
                       dat_keys, dat_labels=None, dat_scale=None,
                       peak_height=None, add_lines=True, add_line_labels=False,
                       min_lat=None, min_sec=None, min_lat_break=None,
                       min_lin_fit=None, max_lat=None, max_sec=None,
                       max_lat_break=None, max_lin_fit=None):
    """Calculate and plot the TID peaks for a given time and altitude range.

    Parameters
    ----------
    sami : xr.Dataset
        SAMI3 concatonated data set
    nt_start: int
        Starting time index
    nt_stop : int
        Ending time index
    nlind : int
        'nl' index
    nfind : int
         'nf' index
    nzinds : list-like
        List of 'nz' indices corresponding to, e.g., the topside ionosphere
    dat_keys : list-like or str
        List of data keys to plot or a satellite string to use defaults
    dat_labels : list-like or NoneType
        List of data labels, will be overwritten if `dat_keys` is a satellite
        string (default=None)
    dat_scale : list-like or NoneType
        List of data scales for keys, if a satellite string is used for
        `dat_keys` this will be reset (default=None)
    peak_height : list-like or NoneType
        List of peak height minima for keys, if a satellite string is used for
        `dat_keys` this will be reset (default=None)
    add_lines : bool
        Add linear fits to plot (default=True)
    add_line_labels : bool
        Add the linear fits to the legend (default=False)
    min_lat : dict or NoneType
        Latitudes corresponding to the minima indices with keys corresponding
        to data variables or None to calculate (default=None)
    min_sec : dict or NoneType
        Seconds from the starting time corresponding to the minima indices with
        keys corresponding to data variables or None to calculate (default=None)
    min_lat_break : dict or NoneType
        Indices of the latitude breaks with keys corresponding to data
        variables for the minima or None to calculate (default=None)
    min_lin_fit : dict or NoneType
        Output from scipy.stats.linregress for each valid fit period or None to
        calculate (default=None)
    max_lat : dict or NoneType
        Latitudes corresponding to the maxima indices with keys corresponding
        to data variables or None to calculate (default=None)
    max_sec : dict or NoneType
        Seconds from the starting time corresponding to the maxima indices with
        keys corresponding to data variables or None to calculate (default=None)
    max_lat_break : dict or NoneType
        Indices of the latitude breaks with keys corresponding to data
        variables for the maxim aor None to calculate (default=None)
    max_lin_fit : dict or NoneType
        Output from scipy.stats.linregress for each valid fit period or None to
        calculate (default=None)

    Returns
    -------
    min_lat : dict
        Latitudes corresponding to the minima indices with keys corresponding
        to data variables
    min_sec : dict
        Seconds from the starting time corresponding to the minima indices with
        keys corresponding to data variables
    min_lat_break : dict
        Indices of the latitude breaks with keys corresponding to data
        variables for the minima
    min_lin_fit : dict
        Output from scipy.stats.linregress for each valid fit period
    max_lat : dict
        Latitudes corresponding to the maxima indices with keys corresponding
        to data variables
    max_sec : dict
        Seconds from the starting time corresponding to the maxima indices with
        keys corresponding to data variables
    max_lat_break : dict
        Indices of the latitude breaks with keys corresponding to data
        variables for the maxima
    max_lin_fit : dict
        Output from scipy.stats.linregress for each valid fit period
    fig : plt.Figure
        Figure handle

    """
    # Get a subset in time of the SAMI3 dataset
    sel_sami = sami.sel(num_times=sami['datetime'].num_times[nt_start:nt_stop])
    start = pds.to_datetime(sami['datetime'][nt_start].values).to_pydatetime()

    # Set the data keys, if needed
    if dat_keys in ['c', 'd']:
        if dat_keys == 'd':
            peak_heights = [5.0, 0.025, 5.0]
        else:
            peak_heights = [None, 0.025, None]

        dat_keys = ["_".join([dkey, dat_keys]) for dkey in [
            'rel_vnq', 'rel_dene', 'rel_vsi2']]
        dat_scale = [100.0, 1.0, 100.0]
        dat_labels = ['$\Delta$ u$_{||}$\nGeo Lat ($^\circ$N)',
                      '$\Delta N_e/N_e$\nGeo Lat ($^\circ$N)',
                      '$\Delta$ v$_{{||}_{O^+}}$\nGeo Lat ($^\circ$N)']

    if None in [dat_scale, dat_labels, peak_heights]:
        raise ValueError('must provide dat_scale, dat_labels, and peak_heights')

    # Calculate any data that is needed
    if None in [min_lat, min_sec, min_lat_break]:
        # Get the peak indices
        min_ind, max_ind = analysis.get_topside_peaks(
            sel_sami, nzinds, dat_keys, dat_scale, peak_heights)

        # Get the data along the minima
        min_lat, min_sec, min_lat_break = analysis.find_linear_breaks(
            sel_sami, nlind, nfind, nzinds, min_ind)
    else:
        max_ind = None

    if None in [max_lat, max_sec, max_lat_break]:
        if max_ind is None:
            # Get the peak indices
            min_ind, max_ind = analysis.get_topside_peaks(
                sel_sami, nzinds, dat_keys, dat_scale, peak_heights)

        # Get the data along the minima
        max_lat, max_sec, max_lat_break = analysis.find_linear_breaks(
            sel_sami, nlind, nfind, nzinds, max_ind)

    if min_lin_fit is None:
        min_lin_fit = analysis.fit_lines_to_peaks(min_lat, min_sec,
                                                  min_lat_break)

    if max_lin_fit is None:
        max_lin_fit = analysis.fit_lines_to_peaks(max_lat, max_sec,
                                                  max_lat_break)

    # Create the figure
    fig = plt.figure(figsize=(7.98, 7.3))
    ax = {dkey: fig.add_subplot(len(dat_keys), 1, 1 + i)
          for i, dkey in enumerate(dat_keys)}

    # Cycle through each of the axes
    ylim = [None, None]
    for i, dkey in enumerate(dat_keys):
        # Plot the minima
        ax[dkey].plot([start + dt.timedelta(seconds=tval)
                       for tval in min_sec[dkey]], min_lat[dkey], 'o',
                      color=pr.nrl_colors(1), label='Minima')

        # Plot the maxima
        ax[dkey].plot([start + dt.timedelta(seconds=tval)
                       for tval in max_sec[dkey]], max_lat[dkey], 's',
                      color=pr.nrl_colors(0), label='Maxima')

        # Add the line fits, if desired
        if add_lines:
            # Cycle through each minima fit
            for fit in min_lin_fit[dkey]:
                xcalc = np.linspace(fit[1], fit[2], 100)
                xplot = [start + dt.timedelta(seconds=x) for x in xcalc]
                yplot = fit[0].slope * xcalc + fit[0].intercept
                if add_line_labels:
                    label = "{:.3f}x {:s} {:.3f}; r={:.2f}".format(
                        fit[0].slope, "+" if fit[0].intercept >= 0 else "-",
                        fit[0].intercept, fit[0].rvalue)
                else:
                    label = None
                ax[dkey].plot(xplot, yplot, "--", color=pr.nrl_colors(3),
                              label=label, lw=2)

            # Cycle through each maxima fit
            for fit in max_lin_fit[dkey]:
                xcalc = np.linspace(fit[1], fit[2], 100)
                xplot = [start + dt.timedelta(seconds=x) for x in xcalc]
                yplot = fit[0].slope * xcalc + fit[0].intercept
                if add_line_labels:
                    label = "{:.3f}x {:s} {:.3f}; r={:.2f}".format(
                        fit[0].slope, "+" if fit[0].intercept >= 0 else "-",
                        fit[0].intercept, fit[0].rvalue)
                else:
                    label = None
                ax[dkey].plot(xplot, yplot, "-.", color=pr.nrl_colors(4),
                              label=label, lw=2)

        # Format the axis
        ax[dkey].set_ylabel(dat_labels[i])
        ax[dkey].yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
        ax[dkey].grid()
        ax[dkey].set_xlim(sel_sami['datetime'].values[0],
                          sel_sami['datetime'].values[-1])
        spec = ax[dkey].get_subplotspec()

        if spec.is_last_row():
            ax[dkey].xaxis.set_major_formatter(
                mpl.dates.DateFormatter('%j\n%H:%M'))
            ax[dkey].set_xlabel('Universal Time')
        else:
            ax[dkey].xaxis.set_major_formatter(
                mpl.ticker.FormatStrFormatter(''))

            if spec.is_first_row():
                ax[dkey].legend(fontsize='small', numpoints=1, loc=1,
                                bbox_to_anchor=(.75, 1.25), ncols=2)

        # Get the y-axis limits
        ax_ylim = ax[dkey].get_ylim()
        if ylim[0] is None or ylim[0] < ax_ylim[0]:
            ylim[0] = ax_ylim[0]

        if ylim[1] is None or ylim[1] > ax_ylim[1]:
            ylim[1] = ax_ylim[1]

    for aa in ax.values():
        aa.set_ylim(ylim)

    return(min_lat, min_sec, min_lat_break, min_lin_fit, max_lat, max_sec,
           max_lat_break, max_lin_fit, fig)
