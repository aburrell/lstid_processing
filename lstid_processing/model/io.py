#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
#
# DISTRIBUTION STATEMENT A: Approved for public release. Distribution is
# unlimited.
# ----------------------------------------------------------------------------
"""Process standard SAMI3 output files."""

import datetime as dt
import numpy as np
import os
import requests
from time import sleep
import xarray as xr

from lstid_processing import logger
from lstid_processing.smoothing.filter_rout import rel_data_butter


def create_concat_files(sami3_dir, output_dir, run_name, date_list=None,
                        samp_period=300, min_period=1800, max_period=7200,
                        nl_inds=None, nf_inds=None, key_inds=None):
    """Create concatonated files from SAMI3 output f1, f2, f3 files.

    Parameters
    ----------
    sami3_dir : str
        Directory where the SAMI3 f1, f2, and f3 files exist
    output_dir : str
        Directory to which the concatonated files will be output
    run_name : str
        Run name to use to distinguish the output files
    date_list : list-like or NoneType
        List containing string specifying the YYYYDDD to concatonate. If None,
        will use ['2014084', '2014085']. (default=None)
    samp_period : int
        Sample period for model data in seconds (default=300)
    min_period : int
        Minimum period for relative variations in seconds (default=1800)
    max_period : int
        Maximum period for relative variations in seconds (default=7200)
    nl_inds : list-like or NoneType
        List of nl indices at which variations will be computed, if None
        uses [26, 26] (default=None)
    nf_inds : list-like or NoneType
        List of nl indices at which variations will be computed, if None
        uses [44, 58] (default=None)
    key_inds : list-like or NoneType
        Key to assign to the nl/nf index pairs, if None uses ['c', 'd']
        (default=None)

    Returns
    -------
    sami : xr.Dataset
        Concatonated dataset with relative and maxima included

    Notes
    -----
    Adds variations at the 26 March 2014 DMSP and C/NOFS conjugate locations.

    """
    # Set the date and indices lists
    if date_list is None:
        date_list = ['2014084', '2014085']

    if nl_inds is None:
        nl_inds = [26, 26]

    if nf_inds is None:
        nf_inds = [44, 58]

    if key_inds is None:
        key_inds = ['c', 'd']

    # Create concatonated file for multiple outputs and days
    sami_list = []
    for date_str in date_list:
        # Open the f1, f2, and f3 files
        sami1 = xr.open_dataset(
            os.path.join(sami3_dir, 'sami3_f1_{:s}.nc'.format(date_str)),
            decode_times=False)
        sami2 = xr.open_dataset(
            os.path.join(sami3_dir, 'sami3_f2_{:s}.nc'.format(date_str)),
            decode_times=False)
        sami3 = xr.open_dataset(
            os.path.join(sami3_dir, 'sami3_f3_{:s}.nc'.format(date_str)),
            decode_times=False)

        # Merge the SAMI3 files together
        sami_list.append(xr.merge([sami1, sami2, sami3]))

        # Create an array of times
        stimes = np.array([
            dt.datetime.strptime('{:d} {:d}'.format(
                int(sami_list[-1]['year'].values),
                int(sami_list[-1]['day'].values)), "%Y %j")
            + dt.timedelta(seconds=int(hrut * 3600))
            for i, hrut in enumerate(sami_list[-1]['hrut'].values)])
        sami_list[-1] = sami_list[-1].assign({'datetime': (('num_times'),
                                                           stimes)})

    # Concat all dates
    sami = xr.concat(sami_list, 'num_times')
    sami.to_netcdf(os.path.join(output_dir,
                                'sami_{:s}_f1_f2_f3_{:s}.nc'.format(
                                    run_name, '_'.join(date_list))))

    # Create file with relative variations and maximums
    stimes = np.array([
        dt.datetime.strptime('{:.0f} {:.0f}'.format(
            sami['year'].values[i], sami['day'].values[i]), "%Y %j")
        + dt.timedelta(seconds=int(hrut * 3600))
        for i, hrut in enumerate(sami['hrut'].values)])

    # Ensure the lists to follow have the right variables
    not_vars = ['year', 'day', 'hrut', 'runtime', 'glat', 'glon', 'zalt',
                'datetime']
    svars = [key for key in sami.data_vars.keys() if key not in not_vars]

    # Get the relative variations
    rel_dict = {}
    for key in svars:
        relative = True if key.find('den') == 0 else False
        rel_dict[key] = dict()
        for ind, ikey in enumerate(key_inds):
            rel_dict[key][ikey] = rel_data_butter(
                sami[key].values[:, nl_inds[ind], nf_inds[ind], :],
                samp_period, min_period=min_period, max_period=max_period,
                axis=0, relative=relative)

    # Reshape the relative variations for xarray
    rel_data_dict = {}
    for dkey in rel_dict.keys():
        for skey in rel_dict[dkey].keys():
            rkey = "_".join(['rel', dkey, skey])
            rel_data_dict[rkey] = (('num_times', 'nz'), rel_dict[dkey][skey])

    # Assign the relative data
    sami = sami.assign(rel_data_dict)

    # For the acceleration terms/special terms, also calculate variations
    # at all latitudes for the longitude slices
    rel_dict = {}
    nlind = np.unique(nl_inds)
    for key in ['u1', 'u2', 'u3']:
        relative = False
        rel_dict[key] = dict()

        for ind, nl in enumerate(nlind):
            lkey = "" if len(nlind) == 1 else "".join((key_inds[ind], "lon"))
            rel_dict[key][lkey] = list()

            for nfind in sami['nf'].values:
                # Calculate the relative variation
                rel_dict[key][lkey].append(rel_data_butter(
                    sami[key].values[:, nl, nfind, :], samp_period,
                    min_period=min_period, max_period=max_period, axis=0,
                    relative=relative))

    # Put the data in xarray format
    rel_data_dict = {}
    for dkey in rel_dict.keys():
        for skey in rel_dict[dkey].keys():
            rkey = "_".join(['rel', dkey, skey])
            if rkey[-1] == "_":
                rkey = rkey[:-1]

            rel_data_dict[rkey] = (('num_times', 'nf', 'nz'),
                                   np.array(rel_dict[dkey][skey]).swapaxes(
                                       1, 0))
    sami = sami.assign(rel_data_dict)

    # Establish which acceleration term is largest
    var_window = int(np.floor(900 / samp_period))
    max_dict = {}
    mkeys = [key for key in sami.data_vars.keys() if key.find("rel_u") == 0
             and key[-2:] not in ['_c', '_d']]
    for dkey in mkeys:
        max_dict[dkey] = list()

        for nfind in sami['nf'].values:
            # Get the rolling object for the absolute value of the relative var
            roll = abs(sami[dkey][:, nfind]).rolling(
                num_times=var_window, min_periods=var_window - 1, center=True)

            # Save the maximum value for each 15 min window
            max_dict[dkey].append(roll.max())

    # Reformat maximum data for xarray
    rel_data_dict = {}
    for dkey in max_dict.keys():
        rkey = "_".join([dkey, 'max'])
        rel_data_dict[rkey] = (('num_times', 'nf', 'nz'),
                               np.array(max_dict[dkey]).swapaxes(1, 0))

    # Assign maximum data and save the output
    sami = sami.assign(rel_data_dict)
    sami.to_netcdf(os.path.join(output_dir,
                                'sami_{:s}_f1_f2_f3_{:s}_rel_max.nc'.format(
                                    run_name, "_".join(date_list))))

    return sami


def load_concat_file(filename):
    """Load a concatonated file SAMI3 file into an xarray Dataset.

    Parameters
    ----------
    filename : str
        Filename with full path.

    Returns
    -------
    sami : xr.Dataset
        An xarray Dataset with model data and fixed times

    """
    # Time decoding doesn't work because there are multiple time variables
    sami = xr.open_dataset(filename, decode_times=False)

    # Construct an array of datetimes
    sami_times = np.array([
        dt.datetime.strptime('{:.0f} {:.0f}'.format(
            sami['year'].values[i], sami['day'].values[i]), "%Y %j")
        + dt.timedelta(seconds=int(hrut * 3600))
        for i, hrut in enumerate(sami['hrut'].values)])

    sami['datetime'] = (('num_times'), sami_times)

    return sami


def download_nrl_files(outdir, filename=None,
                       nrl_site='https://map.nrl.navy.mil/map/pub/nrl/lstid'):
    """Download SAMI3 files from the public NRL directory.

    Parameters
    ----------
    outdir : str
        Directory to which data will be saved
    filename : str or NoneType
        Desired file to download or None to get all of them (default=None)
    nrl_site : str
        URL hosting the SAMI3 files
        (default='https://map.nrl.navy.mil/map/pub/nrl/lstid')

    Returns
    -------
    model_names : list
        List of model run names with their destination directory

    Raises
    ------
    ValueError
        If `outdir` does not exist.

    """
    model_names = list()

    # Test the output directory
    if not os.path.isdir(outdir):
        raise ValueError('Bad output directory: {:}'.format(outdir))

    # Build the download filenames
    if filename is None:
        remote_names = ["fejer_sami3_rel_2014084_85.nc",
                        "model_description.txt",
                        "oneway_sami3_rel_w_hmf2_nmf2_2014084_85.nc",
                        "sami_diagh_oneway_2014084_85_rel_max.nc",
                        "sami_nohpopcoll_f1_f2_f3_2014084_85_rel_max.nc",
                        "twoway_sami3_f2_2014084_85.nc"]
    else:
        remote_names = [filename]

    # Cycle through each remote file and download it to the target directory
    for rfile in remote_names:
        # Build the URL
        url = "/".join([nrl_site, rfile])

        # Build the output filename
        ofile = os.path.join(outdir, rfile)

        # Perform the download
        try:
            with requests.get(url) as req:
                if req.status_code != 404:
                    with open(ofile, 'wb') as fin:
                        fin.write(req.content)
                    logger.info('Successfully downloaded: {:}'.format(ofile))
                    model_names.append(ofile)
                else:
                    logger.info('Remote file unavailable: {:}'.format(url))
        except requests.exceptions.RequestException as err:
            logger.info(''.join([str(err), ' - File: "', rfile,
                                 '" is not available']))

        # Pause to avoid excessive pings to the server
        sleep(0.2)

    return model_names
