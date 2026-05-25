[![tests](https://github.com/bnediction/bonesistools/actions/workflows/tests.yml/badge.svg)](https://github.com/bnediction/bonesistools/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/bonesistools.svg)](https://pypi.org/project/bonesistools)
[![python](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fbnediction%2Fbonesistools%2Frefs%2Fheads%2Fmain%2Fpyproject.toml&query=%24.project.requires-python&style=flat&label=python)](https://www.python.org/)
[![license](https://img.shields.io/pypi/l/bonesistools.svg)](https://github.com/bnediction/bonesistools/blob/main/LICENSE)

# BoNesisTools

`BoNesisTools` is a Python package providing utilities for Boolean modelling, regulatory influence graphs and single-cell analyses around the [BoNesis](https://github.com/bnediction/bonesis) ecosystem.

The package provides:

- Boolean algebra and partial Boolean abstractions
- Boolean network manipulation and analysis
- signed influence graph utilities
- GRN-informed Boolean predecessor inference
- single-cell and multimodal analysis helpers
- biological database interfaces

## Usage

```python
import bonesistools as bt
```

`BoNesisTools` exposes four main namespaces:

- `bt.sct` — single-cell and multimodal annotated data tools
- `bt.bpy` — Boolean modelling and graph utilities
- `bt.dbs` — biological database interfaces
- `bt.grn` — deprecated alias for `bt.bpy.ig`

---

## Single-cell tools

`bt.sct` follows a [Scanpy](https://github.com/scverse/scanpy)-like API while providing additional and complementary features for single-cell analyses.

Submodules:

- preprocessing: `bt.sct.pp`
- tools: `bt.sct.tl`
- plotting: `bt.sct.pl`
- datasets: `bt.sct.datasets`

---

## Boolean modelling utilities

`bt.bpy` provides utilities for Boolean modelling, logical abstractions and signed regulatory graphs.

Submodules:

- Boolean algebra: `bt.bpy.ba`
  - partial Boolean abstractions and hypercube representations
  - Boolean differential and predecessor inference utilities

- Boolean network: `bt.bpy.bn`
  - Boolean network manipulation and analysis
  - fixed-point computation
  - `.bnet` import/export

- influence graph: `bt.bpy.ig`
  - signed regulatory influence graphs
  - feedback circuit and SCC analysis
  - signed interaction scoring from bounded walks
  - graph compression and visualization utilities

Example:

```python
bn = bt.bpy.bn.BooleanNetwork(
    {
        "A": "B & ~C",
        "B": 1,
        "C": 0,
    }
)

graph = bn.to_influence_graph()

graph.show()
```

---

## Biological external resources

`bt.dbs` provides lightweight interfaces and utilities for biological
external resources.

Submodules:

- NCBI: `bt.dbs.ncbi`
  - gene synonym harmonization
  - gene annotation utilities

- OmniPath: `bt.dbs.omnipath`
  - DoRothEA transcription factor interactions
  - CollecTRI regulatory interaction networks

Example:

```python
genesyn = bt.dbs.ncbi.GeneSynonyms()

grn = bt.dbs.omnipath.load_collectri_grn(
    organism="mouse",
    genesyn=genesyn,
)
```

---

## Installation

Install the latest release:

```sh
pip install bonesistools
```

Install the single-cell dependencies:

```sh
pip install "bonesistools[sctools]"
```

Install all optional dependencies:

```sh
pip install "bonesistools[all]"
```

Install the development version:

```sh
git clone https://github.com/bnediction/bonesistools.git
cd bonesistools
pip install -e ".[all]"
```

or directly:

```sh
pip install git+https://github.com/bnediction/bonesistools.git
```

---

## Bugs

Please report bugs or ask questions here:

https://github.com/bnediction/bonesistools/issues

---

## License

This package is distributed under the [CeCILL v2.1](http://www.cecill.info/index.en.html) free software license (GNU GPL compatible).

This package also includes third-party data resources derived from the
NCBI Gene database (`gene_info`). NCBI places no restrictions on the
use or redistribution of these data: https://www.ncbi.nlm.nih.gov/home/about/policies/
