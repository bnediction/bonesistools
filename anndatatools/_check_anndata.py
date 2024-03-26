#!/usr/bin/env python

import types
from anndata import AnnData

def _adata_arg_checking(
    function: types.FunctionType
):
    """Check if the first argument is an anndata instance.
    
    Parameters
    ----------
    function
        called function
    
    Returns
    -------
    `function` argument result.
    Raise an error if the first argument is not an anndata instance.
    """

    def wrapper(adata, *args, **kwargs):
        if not isinstance(adata, AnnData):
            raise TypeError(f"Argument `adata` must be of type {type(AnnData)}, not {type(adata)}")
        return function(adata, *args, **kwargs)
    
    return wrapper
