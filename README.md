Repository with tools for identifying medium and large scale Travelling
Ionospheric Disturbances (TIDs) in the Communication/Navigation Outage
Forecasting System (C/NOFS) Coupled Ion-Neutral Dynamics Investigation (CINDI)
Ion Velocity Meter (IVM) data, tools for working with SAMI3 runs, and model
runs used for a LSTID case study.

[DOI HERE] [PYPI HERE] [DOI FOR MODEL DATA]

Example
-------

To load one of the SAMI3 files into an xarray Dataset:

```
import lstid_processing
import os

sami3 = lstid_processing.model.io.load_concat_file(
    os.path.join(lstid_processing.model_run_dir,
    'oneway_sami3_rel_w_hmf2_nmf2_2014084_85.nc')
```

Notes
-----

This package and data are supplied to support the reproducibility of a
manuscript currently under review.  Frequent updates are not expected.


