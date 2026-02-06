Repository with tools for identifying medium and large scale Travelling
Ionospheric Disturbances (TIDs) in the Communication/Navigation Outage
Forecasting System (C/NOFS) Coupled Ion-Neutral Dynamics Investigation (CINDI)
Ion Velocity Meter (IVM) data, tools for working with SAMI3 runs, and a routine
to obtain model runs used for a LSTID case study are provided.

[![DOI](https://zenodo.org/badge/989170648.svg)](https://doi.org/10.5281/zenodo.15528276) [![PyPI version](https://badge.fury.io/py/lstid_processing.svg)](https://badge.fury.io/py/lstid_processing) [![Documentation Status](https://readthedocs.org/projects/lstid-processing/badge/?version=latest)](https://lstid-processing.readthedocs.io/en/latest/?badge=latest)


Example
-------

To download all of the model runs to a local directory:

```
import lstid_processing

sami3_files = lstid_processing.model.io.download_nrl_files('path/to/downloads')
```

Notes
-----

This package and data are supplied to support the reproducibility of
Burrell et al., (2026) doi:10.1029/2025JA034335. Frequent updates are not
expected.


