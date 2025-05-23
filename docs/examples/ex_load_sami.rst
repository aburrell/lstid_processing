.. _ex-load:

Download Model Data
===================

Any SAMI3 data may be used by this package, but the files used in the study
that this code was created to support may be obtained from NRL. To download
one or all of the files from NRL, use the
:py:func:`lstid_processing.model.io.download_nrl_files` function.

::

   import lstid_processing

   # Set the output directory to be an existing directory on your computer
   outdir = "/path/to/hold/model/data"

   # Download just one file; to download all, skip the `filename` kwarg
   sami_files = lstid_processing.model.io.download_nrl_files(
       outdir, filename="oneway_sami3_rel_w_hmf2_nmf2_2014084_85.nc")



Load Model Data
===============

To load the downloaded model data file, use
:py:func:`lstid_processing.model.io.load_concat_file`.

::

  sami3 = lstid_processing.model.io.load_concat_file(sami_files[0])


This will load the data into an :py:class:`xarray.Dataset` and ensure the times
are processed correctly.
