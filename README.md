Repository with tools for identifying medium and large scale Travelling
Ionospheric Disturbances (TIDs) in the Communication/Navigation Outage
Forecasting System (C/NOFS) Coupled Ion-Neutral Dynamics Investigation (CINDI)
Ion Velocity Meter (IVM) data, tools for working with SAMI3 runs, and a routine
to obtain model runs used for a LSTID case study are provided.

[DOI HERE] [PYPI HERE]

Example
-------

To download all of the model runs to a local directory:

```
import lstid_processing

sami3_files = lstid_processing.model.io.download_nrl_files('path/to/downloads')
```

Notes
-----

This package and data are supplied to support the reproducibility of a
manuscript currently under review.  Frequent updates are not expected.


