#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
# ----------------------------------------------------------------------------
"""Package to identify travelling ionospheric disturbances."""

import logging
from importlib import metadata

# Define a logger object to allow easier log handling
logging.raiseExceptions = False
logger = logging.getLogger('cindi_tids')

# Import the package modules and top-level classes
from cindi_tids import analysis  # noqa F401
from cindi_tids import plot_rout  # noqa F401
from cindi_tids import smoothing  # noqa F401

# Define the global variables
__version__ = metadata.version('cindi_tids')
