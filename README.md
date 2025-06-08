[![python](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fbnediction%2Fbonesistools%2Frefs%2Fheads%2Fmain%2Fpyproject.toml&query=%24.project.requires-python&style=flat&label=python)](https://www.python.org/)
[![PyPI](https://img.shields.io/pypi/v/bonesistools.svg)](https://pypi.org/project/bonesistools)
[![license](https://img.shields.io/pypi/l/bonesistools.svg)](https://github.com/bnediction/bonesistools/blob/main/LICENSE)

# BoNesisTools

BoNesisTools is a python package proposing bioinformatics tools for upstream and downstream analysis of [BoNesis](https://github.com/bnediction/bonesis) framework.

## Functionalities

BoNesisTools contains multiple sub-packages:
* sctools (operations on unimodal and multimodal annotated data)
* boolpy (features for Boolean algebra and Boolean functions)
* grntools (efficient features and algorithms for gene regulatory network-like graphs)
* databases (information requests to bioinformatics databases)

Note that sctools proposes a wide range of efficient features that are not available in [Scanpy](https://github.com/scverse/scanpy) and [muon](https://github.com/scverse/muon) packages.

## Installation

Install the latest released version of BoNesisTools with pip:
```sh
pip install bonesistools
```
Install with all optional dependencies:
```sh
pip install bonesistools[default]
```

for the development version, use:
```sh
git clone https://github.com/bnediction/bonesistools.git
```
or:
```sh
pip install git+https://github.com/bnediction/bonesistools.git
```

## Usage

For importing package:
```python
import bonesistools as bt
```
Each sub-package has a corresponding alias:
* sctools can be accessed with `bt.sct`, which is decomposed into three Scanpy-like parts: preprocessing (`bt.sct.pp`), tools (`bt.sct.tl`) and plotting (`bt.sct.pl`)
* boolpy can be accessed with `bt.bpy`
* grntools can be accessed with `bt.grn`
* databases can be accessed with `bt.dbs`, which is decomposed into two parts for the moment: ncbi (`bt.dbs.ncbi`) for handling gene alias-related issues and [CollecTRI](https://github.com/saezlab/CollecTRI) (`bt.dbs.collectri`) for getting a gene regulatory networks

## Bugs

Please report any bugs or ask questions [here](https://github.com/bnediction/bonesistools/issues).

## License

This package is distributed under the [CeCILL v2.1](http://www.cecill.info/index.en.html) free software license (GNU GPL compatible).
