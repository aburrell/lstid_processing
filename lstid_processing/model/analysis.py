#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
# ----------------------------------------------------------------------------
"""Analyse concatonated SAMI3 files."""

import datetime as dt
import numpy as np
import pandas as pds
from scipy import signal
from scipy import stats


def get_default_indices(time_array=None, lat_array=None):
    """Get the default indices for the C/NOFS and DMSP F16 field lines.

    Parameters
    ----------
    time_array : array-like or NoneType
        Array of times, from which the conjugate time index may be selected.
    lat_array : array-like or NoneType
        Array of latitudes with dimensions of time x nl x nf x nz.

    Returns
    -------
    nlind : int
        SAMI3 'nl' index for the C/NOFS and DMSP F16 location
    nfindc : int
        SAMI3 'nf' index for the C/NOFS location
    nfindc : int
        SAMI3 'nf' index for the DMSP F16 location
    nt : int or NoneType
        None if `time_array` is None, time index for the satelite conjuction
        otherwise
    nzindc : int
        None if `lat_array` or `time_array` is None; otherwise the altitude
        index for the C/NOFS location
    nzindd : int
        None if `lat_array` is None or `time_array` is None; otherwise altitude
        index for the DMSP F16 location

    Notes
    -----
    C/NOFS is at: 11.3 N, 33.2 E, 479.4 km
    DMSP F16 is at: 11.3 N, 33.2 E, 845.8 km
    Conjunction at: 26 March 2014 at 02:49:59 UT

    """
    # CINDI field line: 11.3 N, 33.2 E, 479.4 km
    # DMSP field line: 11.3 N, 33.2 E, 845.8 km
    nlind = 26
    nfindc = 44
    nfindd = 58

    if time_array is None:
        nt = None
        nzindc = None
        nzindd = None
    else:
        nt = get_time_index(time_array, dt.datetime(2014, 3, 26, 2, 49, 59))

        if lat_array is None:
            nzindc = None
            nzindd = None
        else:
            nzindc = abs(lat_array[nt, nlind, nfindc] - 11.3).argmin()
            nzindd = abs(lat_array[nt, nlind, nfindd] - 11.3).argmin()

    return nlind, nfindc, nfindd, nt, nzindc, nzindd


def get_time_index(time_array, dtime):
    """Retrieve the index corresponding to the closest time value.

    Parameters
    ----------
    time_array : array-like
        Array of times, from which the desired time index may be selected.
    dtime : dt.datetime
        Desired time as a datetime object

    Returns
    -------
    nt : int
        Time index

    """
    nt = int(abs(pds.to_datetime(time_array).to_pydatetime() - dtime).argmin())
    return nt


def get_f2_peaks(nlind, nfind, sami):
    """Get the F2 peak locations for a specified field line.

    Parameters
    ----------
    nlind : int
        SAMI3 `nl` index for the field line
    nfind : int
        SAMI3 `nf` index for the field line
    sami : xr.Dataset
        Dataset with SAMI3 values

    Returns
    -------
    f2_inds : dict
       Dict with keys 'south' and 'north' that contain arrays of indicies
       corresponding to the altitude indices for the F2 peak across all times.
       Note that these will only be valid if the field line reaches into the
       topside ionosphere.

    """
    # Get the halfway point along the field line
    phalf = int(sami['glat'][0, nlind, nfind].shape[0] / 2)

    # Get the index of the density peak along each half of the field line
    f2_inds = {
        'south': sami['dene'].values[:, nlind, nfind, :phalf].argmax(axis=1),
        'north': phalf + sami['dene'].values[:, nlind, nfind,
                                             phalf:].argmax(axis=1)}

    return f2_inds


def get_dominant_acceleration(sami):
    """Get the dominant acceleration variation term indices.

    Parameters
    ----------
    sami : xr.Dataset
        SAMI3 concatonated data with acceleration variation terms

    Returns
    -------
    dom_acc : array-like
        Array of numbers indicating which acceleration term has the greatest
        variability; has values of 1, 2, and 3

    Notes
    -----
    1 is the plasma pressure gradient, 2 is the ion-neutral collisions, and
    3 is the ion-ion collisions

    """
    # Get the indexes where the collisions are the highest
    ineu = np.where((sami['rel_u2_max'].values >= sami['rel_u1_max'].values)
                    & (sami['rel_u2_max'].values >= sami['rel_u3_max'].values))

    iion = np.where((sami['rel_u3_max'].values >= sami['rel_u1_max'].values)
                    & (sami['rel_u3_max'].values >= sami['rel_u2_max'].values))

    # Initalize the output to specify that plasma pressure is highest
    dom_acc = np.ones(shape=sami['rel_u3_max'].values.shape)

    # Reassign indices where collisions are greater
    dom_acc[ineu] = 2
    dom_acc[iion] = 3

    return dom_acc


def get_topside_peaks(sami, nzinds, dat_keys, dat_scale, peak_height):
    """Get the variation peaks across a range of altitudes.

    Parameters
    ----------
    sami : xr.Dataset
        Dataset with SAMI3 values for the desired time period
    nzinds : list-like
        Lists of altitude indices at which peaks will be calculated
    dat_keys : list-like
        List of data keys to evaluate, expects relative variation keys like
        ['rel_dene_d']
    dat_scale : list-like
        List of scaling parameters to divide from the data
    peak_height : list-like
        List of heights that peaks must meet in order to qualify in the
        detection, can be None to include all peaks regardless of significance

    Returns
    -------
    min_ind : dict
        Dict with keys corresponding to data keys and values contianing a list
        of indices that locate significant peak minima
    max_ind : dict
        Dict with keys corresponding to data keys and values contianing a list
        of indices that locate significant peak maxima

    """
    # Initialize the output
    min_ind = {dkey: list() for dkey in dat_keys}
    max_ind = {dkey: list() for dkey in dat_keys}

    # Cycle through each data type
    for i, dkey in enumerate(dat_keys):
        # Cycle through all altitudes
        for nzind in nzinds:
            # Find the maxima across the desired time period
            peak_ind, _ = signal.find_peaks(
                sami[dkey][:, nzind].values / dat_scale[i],
                height=peak_height[i])
            max_ind[dkey].append(peak_ind)

            # Find the minima across the desired time period
            peak_ind, _ = signal.find_peaks(
                -1.0 * sami[dkey][:, nzind].values / dat_scale[i],
                height=peak_height[i])
            min_ind[dkey].append(peak_ind)

    return min_ind, max_ind


def find_linear_breaks(sami, nlind, nfind, nzinds, peak_ind, lat_break=6.0,
                       sec_break=900, max_lat=20.0):
    """Identify the breaks in peaks that indicate a change in wavefront.

    Parameters
    ----------
    sami : xr.Dataset
        Dataset with SAMI3 values for the desired time period
    nlind : int
        'nl' index for the relative data included in the minima and maxima
        data keys
    nfind : int
        'nf' index for the relative data included in the minima and maxima
        data keys
    nzinds : list-like
        Lists of altitude indices at which peaks were calculated
    peak_ind : dict
        Dict with keys corresponding to data keys and values contianing a list
        of indices that locate significant peak minima or maxima
    lat_break : float
        Change in latitude that indicates a break (default=6.0)
    sec_break : int
        Change in time (seconds) that indicates a break (default=900)
    max_lat : float
        Maximum latitude for data, starting point for breaks (default=20.0)

    Returns
    -------
    peak_lat : dict
        Latitudes corresponding to the peak indices with keys corresponding to
        data variables
    peak_sec : dict
        Seconds from the starting time corresponding to the peak indices with
        keys corresponding to data variables
    lat_break_inds : dict
        Indices of the latitude breaks with keys corresponding to data
        variables

    Notes
    -----
    Only identifies North-to-South trends automatically

    """
    # Initialize the output
    peak_lat = dict()
    peak_sec = dict()
    lat_break_inds = dict()

    # Cycle through all the data keys
    start = pds.to_datetime(sami['datetime'].values[0]).to_pydatetime()
    for dkey in peak_ind.keys():
        # Initalize the output lists
        peak_lat[dkey] = list()
        peak_sec[dkey] = list()
        lat_break_inds[dkey] = list()

        # Cycle through the altitudes to gather the peak locations
        for i, nzind in enumerate(nzinds):
            # Get the latitudes and times for the peak locations
            if len(peak_ind[dkey][i]) > 0:
                lat = sami['glat'].values[
                    peak_ind[dkey][i], nlind, nfind, nzind]
                sec = [tval.total_seconds() for tval in (
                    pds.to_datetime(sami['datetime'].values[
                        peak_ind[dkey][i]]).to_pydatetime() - start)]

                # Save the output
                peak_lat[dkey].extend(list(lat))
                peak_sec[dkey].extend(sec)

        # Find breaks in the locations. Start by sorting data by location
        # relative to the maximum in the North and time
        lat = max_lat - np.array(peak_lat[dkey])
        isort = np.lexsort((lat, peak_sec[dkey]))

        peak_lat[dkey] = np.array(peak_lat[dkey])[isort]
        peak_sec[dkey] = np.array(peak_sec[dkey])[isort]

        # Calculate the latitude and time differences of the sorted data
        lat_diff = abs(peak_lat[dkey][1:] - peak_lat[dkey][:-1])
        sec_diff = peak_sec[dkey][1:] - peak_sec[dkey][:-1]

        # Get the latitude break indices
        lat_break_inds[dkey] = np.where((lat_diff > lat_break)
                                   | (sec_diff > sec_break))[0] + 1

    return peak_lat, peak_sec, lat_break_inds


def fit_lines_to_peaks(peak_lat, peak_sec, lat_break, min_pnts=40):
    """Calculate linear fits to the maxima and minima for each pass.

    Parameters
    ----------
    peak_lat : dict
        Latitudes corresponding to the peak indices with keys corresponding to
        data variables
    peak_sec : dict
        Seconds from the starting time corresponding to the peak indices with
        keys corresponding to data variables
    lat_break : dict
        Indices of the latitude breaks with keys corresponding to data
        variables
    min_pnts : int
        Minimum number of points needed to calculate a linear fit (default=40)

    Returns
    -------
    lin_fit : dict
        Output from scipy.stats.linregress for each valid fit period

    """
    # Initialize the output
    lin_fit = {dkey: list() for dkey in peak_lat.keys()}

    # Cycle through all data variables
    for dkey in lin_fit.keys():
        # Cycle through each latitude break
        for i, ibreak in enumerate(lat_break[dkey]):
            # Set the start index
            istart = 0 if i == 0 else lat_break[dkey][i - 1]

            # Get the range
            tmin = peak_sec[dkey][istart:ibreak][0]
            tmax = peak_sec[dkey][istart:ibreak][-1]

            # Fit the data for this pass
            if tmin < tmax and ibreak - istart > min_pnts:
                lin_fit[dkey].append([stats.linregress(
                    peak_sec[dkey][istart:ibreak],
                    peak_lat[dkey][istart:ibreak]), tmin, tmax])

        # Set the start index for the last pass
        istart = lat_break[dkey][-1]

        # Get the range
        tmin = peak_sec[dkey][istart:][0]
        tmax = peak_sec[dkey][istart:][-1]

        # Fit the data for this pass
        if tmin < tmax and len(peak_sec[dkey]) - istart > min_pnts:
            lin_fit[dkey].append([stats.linregress(
                peak_sec[dkey][istart:], peak_lat[dkey][istart:]), tmin, tmax])

    return lin_fit
