#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
#
# DISTRIBUTION STATEMENT A: Approved for public release. Distribution is
# unlimited.
# ----------------------------------------------------------------------------
"""Data filtering routines, specifically designed to support TID analysis."""

import numpy as np
from scipy import signal
from scipy import stats

from lstid_processing.smoothing import fill_rout


def rel_data_butter(tdata, samp_period, min_period=-1, max_period=-1, axis=0,
                    relative=True):
    """Calculate (relative) differences in data using a Butterworth filter.

    Parameters
    ----------
    tdata : array-like
        Data organized consistently by time
    samp_period : int or float
        Sample period for data in seconds
    min_period : int
        Minimum period in seconds, negative for high-pass (default=-1)
    max_period : int
        Maximum period in seconds, negative for low-pass (default=-1)
    axis : int
        Axis along which data should be filtered (default=0)
    relative : bool
        If True, calculate relative differences; if False return absolute
        differences (default=True)

    Returns
    -------
    rdata : array-like
        Unitless, filtered changes in the data, of the same shape and time
        order as tdata.  Calculated by normalizing the filtered data by the
        input data

    Notes
    -----
    When both min_period and max_period are provided, a bandpass filter is used

    Raises
    ------
    ValueError
        If no valid period is supplied

    See Also
    --------
    scipy.signal.butter, scipy.signal.sosfiltfilt

    """
    # Ensure the data is array-like
    tdata = np.asarray(tdata)

    # Move from the time to the frequency domain
    samp_freq = 1.0 / samp_period

    # Determine if a low-pass, high-pass, or no filter will be applied,
    # converting from seconds to Hz
    butter_pass = None
    if min_period >= 0.0:
        freq_lim = 1.0 / min_period
        butter_pass = 'lowpass'

    if max_period >= 0.0:
        if butter_pass is None:
            freq_lim = 1.0 / max_period
            butter_pass = 'highpass'
        else:
            freq_lim = [1.0 / max_period, freq_lim]
            butter_pass = 'bandpass'

    if butter_pass is None:
        raise ValueError('must provide at least a minimum or maximum period')

    # Obtain a 12th order digital Butterworth filter
    butter_sos = signal.butter(12, freq_lim, butter_pass, output='sos',
                               fs=samp_freq)

    # Filter the data
    fdata = signal.sosfiltfilt(butter_sos, tdata, axis=axis)

    # Calculate the relative data
    rdata = fdata / tdata if relative else fdata

    return rdata


def pysat_rel_data_butter(inst, dat_key, out_key='', samp_period=None,
                          min_period=-1, max_period=-1, axis=0, nan_val=0.0,
                          relative=True):
    """Filter data within pysat as a custom function.

    Parameters
    ----------
    inst : pysat.Instrument
        pysat Instrument object
    dat_key : str
        Key to access the desired data variable in the pysat Instrument
    out_key : str
        Key to which the relative data will be saved, or an empty string to
        construct a name using the dat_key and min/max periods (default='')
    samp_period : int, float, or None
        Sample period for data in seconds or None to determine from the
        pysat Instrument time index (default=None)
    min_period : int
        Minimum period in seconds, negative for high-pass (default=-1)
    max_period : int
        Maximum period in seconds, negative for low-pass (default=-1)
    axis : int
        Axis along which data should be filtered (default=0)
    nan_val : float
        Replace NaNs in the data with this value before filtering, if they
        cannot be replaced during interpolation. (default=0.0)
    relative : bool
        If True, calculate relative differences; if False return absolute
        differences (default=True)

    Raises
    ------
    ValueError
        If no valid period is supplied, out_key already exists, or data is not
        regularly sampled

    See Also
    --------
    rel_data_butter, scipy.signal.butter, scipy.signal.sosfiltfilt

    """
    # Set the output key
    if len(out_key) == 0:
        out_key = '{:s}_rel_butter'.format(dat_key)

        if min_period >= 0.0:
            out_key = '_'.join([out_key, 'Tmin{:d}s'.format(int(min_period))])

        if max_period >= 0.0:
            out_key = '_'.join([out_key, 'Tmax{:d}s'.format(int(max_period))])

    if out_key in inst.data.keys():
        raise ValueError('relative data found in {:}'.format(out_key))

    # Determine the sample period, if needed
    if samp_period is None:
        if inst.index.freq is not None:
            samp_period = inst.index.freq.delta.total_seconds()
        else:
            # Use mode instead of min to avoid small rounding errors
            ind_diff = (inst.index[1:] - inst.index[:-1]).total_seconds()
            samp_period = stats.mode(ind_diff).mode

            # If the data is multi-modal, exit with a value error
            if len(samp_period) == 1:
                samp_period = samp_period[0]
            else:
                raise ValueError('unable to determine data sample rate')

    # If there are fill values, replace them with interpolated values
    if np.isnan(inst[dat_key]).any() and not np.isnan(inst[dat_key]).all():
        good_dat = fill_rout.fill_data(inst[dat_key].copy(), inst.index)
        good_dat[np.isnan(good_dat)] = nan_val
    else:
        good_dat = inst[dat_key]

    # Get the relative data
    rdata = rel_data_butter(good_dat, samp_period, min_period=min_period,
                            max_period=max_period, axis=axis, relative=relative)

    # Assign the data to the pysat object
    if inst.pandas_format:
        inst.data = inst.data.assign(**{out_key: rdata})
    else:
        inst.data = inst.data.assign({out_key: rdata})

    # Assign the metadata to the pysat object
    ftype = None
    tstr = ''
    if min_period >= 0.0:
        ftype = 'lowpass'
        tstr = '{:d} sec'.format(int(min_period))

    if max_period >= 0.0:
        if ftype is None:
            ftype = 'highpass'
            tstr = '{:d} sec'.format(int(max_period))
        else:
            ftype = 'bandpass'
            tstr = ' and '.join([tstr, '{:d} sec'.format(int(max_period))])

    dat_name = inst.meta[dat_key, inst.meta.labels.name]
    desc = ''.join(['Relative variations in ', dat_name, ' as specified by ',
                    'the ratio between the Butterworth filtered and ',
                    'unfiltered data. The filter used a ', ftype,
                    ' with period limit', 's' if ftype == 'bandpass' else '',
                    ' of ', tstr])
    inst.meta[out_key] = {inst.meta.labels.units: 'unitless',
                          inst.meta.labels.name: 'Filtered {:s} / {:s}'.format(
                              dat_name, dat_name),
                          inst.meta.labels.desc: desc,
                          inst.meta.labels.min_val: -np.inf,
                          inst.meta.labels.max_val: np.inf,
                          inst.meta.labels.fill_val: np.nan}

    return
