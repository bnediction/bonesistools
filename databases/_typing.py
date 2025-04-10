#!/usr/bin/env python

import importlib

try:
    _mpbn_is_available = importlib.util.find_spec("mpbn") is not None
except:
    _mpbn_is_available = importlib.find_loader("mpbn") is not None

if _mpbn_is_available:
    from mpbn import MPBooleanNetwork
else:
    MPBooleanNetwork = type(NotImplemented)
