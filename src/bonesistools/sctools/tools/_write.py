#!/usr/bin/env python

from typing import Optional, Union

from pathlib import Path

import pandas as pd
from anndata import AnnData

from scipy import io, sparse

PathLike = Union[str,Path]

def to_csv(
    adata: AnnData,
    filename: PathLike,
    layer: Optional[str] = None
) -> None:

    filename = str(filename)
    if not filename.endswith(".csv"):
        filename=Path(f"{filename}.csv")

    if layer is None:
        x = adata.X
    else:
        x = adata.layers[layer]

    pd.DataFrame(x).to_csv(
        path_or_buf=filename,
        sep=","
    )

def to_mtx(
    adata: AnnData,
    filename: PathLike,
    layer: Optional[str] = None
) -> None:
    
    filename = str(filename)
    if not filename.endswith(".mtx"):
        filename=Path(f"{filename}.mtx")
    
    if layer is None:
        x = adata.X
    else:
        x = adata.layers[layer]
    
    io.mmwrite(filename, x)

def to_npz(
    adata: AnnData,
    filename: PathLike,
    layer: Optional[str] = None
) -> None:
    
    filename = str(filename)
    if not filename.endswith(".npz"):
        filename=Path(f"{filename}.npz")

    if layer is None:
        x = adata.X
    else:
        x = adata.layers[layer]
    
    sparse.save_npz(filename, x)

def to_csv_or_mtx(
    adata: AnnData,
    filename: PathLike,
    layer: Optional[str] = None,
) -> None:
    
    if layer is None:
        if not sparse.issparse(adata.X):
            to_csv(adata, filename)
        else:
            to_mtx(adata, filename)
    else:
        if not sparse.issparse(adata.layers[layer]):
            to_csv(adata, filename, layer)
        else:
            to_mtx(adata, filename, layer)

def to_csv_or_npz(
    adata: AnnData,
    filename: PathLike,
    layer: Optional[str] = None,
) -> None:
    
    if layer is None:
        if not sparse.issparse(adata.X):
            to_csv(adata, filename)
        else:
            to_npz(adata, filename)
    else:
        if not sparse.issparse(adata.layers[layer]):
            to_csv(adata, filename, layer)
        else:
            to_npz(adata, filename, layer)
