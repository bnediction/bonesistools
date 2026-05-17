#!/usr/bin/env python

"""
Interfaces to biological databases and prior knowledge resources.

The `dbs` package provides utilities for retrieving, parsing and
standardising biological information from external resources.

Current sub-packages
--------------------
ncbi
    Gene nomenclature and synonym handling utilities.
omnipath
    Access to regulatory and interaction prior knowledge databases.
"""

from . import omnipath, ncbi
