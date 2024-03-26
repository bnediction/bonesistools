#!/usr/bin/env python

from distutils.core import setup

setup(
    name="anndatatools",
    version="1.0.0",
    description="python package for handling and performing operations on annotated data (anndata)",
    author="Théo Roncalli",
    author_email="theo.roncalli@labri.fr",
    packages=["color_settings", "plot_settings", "genemarkers", "neighborgraph", "visualization"]
)