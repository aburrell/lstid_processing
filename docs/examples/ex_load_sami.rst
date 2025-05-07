.. _ex-load:

Load Model Data
===============

The SAMI3 data included in this package is located in the directory specified
by the variable :py:var:`lstid_processing.model_run_dir`.  To load this data::

  import lstid_processing
  import os

  sami3 = lstid_processing.model.io.load_concat_file(
      os.path.join(lstid_processing.model_run_dir,
      'oneway_sami3_rel_w_hmf2_nmf2_2014084_85.nc')

This will load the data into an :py:class:`xarray.Dataset`.
