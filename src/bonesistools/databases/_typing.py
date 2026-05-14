#!/usr/bin/env python

try:
    from mpbn import MPBooleanNetwork
    _mpbn_is_available = True
    _mpbn_import_error = None
except ImportError as error:
    _mpbn_is_available = False
    _mpbn_import_error = error
    MPBooleanNetwork = type(NotImplemented)