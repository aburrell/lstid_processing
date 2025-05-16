#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
#
# DISTRIBUTION STATEMENT A: Approved for public release. Distribution is
# unlimited.
# ----------------------------------------------------------------------------
"""Data filtering routines, specifically designed to support TID analysis."""

import datetime as dt
from matplotlib import dates
import numpy as np
from scipy import interpolate


def fill_data(data, data_time=None, data_loc=None, min_val=np.nan,
              max_val=np.nan, fill_val=np.nan, method='linear'):
    """Pad instances of unspecified or bad data using grid interpolation.

    Parameters
    ----------
    data : array-like
        ND data array with potential bad values
    data_time : array-like or NoneType
        Temporal data as datetime objects, must be along axis 0 of data array
        (default=None)
    data_loc : list of array-like or None
        Location coordinate data, contained as a list, in axis order
        corresponding to their order in the data array.  For example, in a data
        array with time along the first axis and altitude along the second
        axis, this would be a list with the first and only element containing
        an array of the altitude data. Alternatively, for a data array with
        longitude along the first axis and altitude along the second axis,
        time_data would be none, and this would be a list with the longitude
        array as the first element and the altitude array as the second
        element. (default=None)
    min_val : float
        Minimum allowed value for data, applied if a number (default=NaN)
    max_val : float
        Maximum allowed value for data, applied if a number (default=NaN)
    fill_val : float
        Value used to fill in for requested points outside of the convex hull
        of the input points. This option has no effect for the ‘nearest’
        method. (default=np.nan)
    method : str
        Interpolation method, see scipy.interpolate.griddata (default='linear')

    Returns
    -------
    good_data : array-like
        ND array with no bad values

    See Also
    --------
    scipy.interpolate.griddata

    """
    # Cast all data as arrays
    good_data = np.asarray(data)

    # Organize the coordinates
    coords = list()

    if data_time is not None:
        coords.append(np.asarray(dates.date2num(data_time)))

    if data_loc is not None:
        for ldata in data_loc:
            coords.append(np.asarray(ldata))

    coords = np.asarray(coords)

    # If the extrema are not to be used, set them to the highest and lowest
    # data values
    if np.isnan(min_val):
        min_val = np.nanmin(good_data)

    if np.isnan(max_val):
        max_val = np.nanmax(good_data)

    # Select the good data points
    mgood = (~np.isnan(good_data)
             & np.greater_equal(good_data, min_val, where=~np.isnan(good_data))
             & np.less_equal(good_data, max_val, where=~np.isnan(good_data)))
    igood = np.where(mgood)
    ibad = np.where(~mgood)

    # Get the good coordinates and values
    points = np.array([cdat[igood[i]]
                       for i, cdat in enumerate(coords)],
                      dtype=object).transpose()
    xi = np.array([cdat[ibad[i]] for i, cdat in enumerate(coords)],
                  dtype=object).transpose()
    idata = interpolate.griddata(points.astype(float), good_data[igood],
                                 xi.astype(float), method=method,
                                 fill_value=fill_val)

    # Replace the bad values with interpolated values
    try:
        good_data[ibad] = idata
    except ValueError:
        # For 1D data, the output shape differs along the last axis
        good_data[ibad] = idata[:, 0]

    return good_data


def fill_time_series(data_time, data, samp_period, min_val=np.nan,
                     max_val=np.nan, method='linear', fill_val=np.nan):
    """Pad instances of unspecified or bad data using 1D interpolation.

    Parameters
    ----------
    data_time : array-like
        Temporal data as datetime objects, must be along axis 0 of data array
    data : array-like
        1D data array with potential bad values
    samp_period : float
        Sample period in minutes at which the data should be observed
    min_val : float
        Minimum allowed value for data, applied if a number (default=NaN)
    max_val : float
        Maximum allowed value for data, applied if a number (default=NaN)
    method : str
        Interpolation method, see `kind` in scipy.interpolate.interp1d
        (default='linear')
    fill_val : array-like or 'extrapolate'
        Fill value if no interpolation possible, see scipy.interpolate.interp1d
        (default=np.nan)

    Returns
    -------
    good_time : array-like
        1D array of time data without gaps
    good_data : array-like
        1D array with no bad values

    See Also
    --------
    scipy.interpolate.interp1d

    """
    # Cast all data as arrays
    xdat = np.asarray(dates.date2num(data_time))
    data = np.asarray(data)

    # If the extrema are not to be used, set them to the highest and lowest
    # data values
    if np.isnan(min_val):
        min_val = np.nanmin(data)

    if np.isnan(max_val):
        max_val = np.nanmax(data)

    # Calculate the expected temporal cadance
    good_time = [data_time[0]]
    i = 0
    while good_time[-1] < data_time[-1]:
        i += 1
        good_time.append(data_time[0] + dt.timedelta(minutes=samp_period * i))
    xgood = np.asarray(dates.date2num(good_time))

    # Select the good data points
    mgood = (~np.isnan(data)
             & np.greater_equal(data, min_val, where=~np.isnan(data))
             & np.less_equal(data, max_val, where=~np.isnan(data)))
    igood = np.where(mgood)

    # Get the good coordinates and values
    idata = interpolate.interp1d(xdat[igood], data[igood], kind=method,
                                 fill_value=fill_val)

    # Retrieve the interpolated values at the desired times
    good_data = np.array([idata(x) for x in xgood])

    return good_time, good_data
