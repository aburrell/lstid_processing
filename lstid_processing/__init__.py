#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Full license can be found in License.md
#
# DISTRIBUTION STATEMENT A: Approved for public release. Distribution is
# unlimited.
# ----------------------------------------------------------------------------
"""Package to identify and process travelling ionospheric disturbances."""

import logging
from importlib import metadata

# Define a logger object to allow easier log handling
logging.raiseExceptions = False
logger = logging.getLogger('lstid_processing')

# Import the package modules and top-level classes
from lstid_processing import cindi  # noqa F401
from lstid_processing import model  # noqa F401
from lstid_processing import plot_rout  # noqa F401
from lstid_processing import smoothing  # noqa F401

# Define the global variables
__version__ = metadata.version('lstid_processing')
