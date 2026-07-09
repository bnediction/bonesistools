[![tests](https://github.com/bnediction/bonesistools/actions/workflows/tests.yml/badge.svg)](https://github.com/bnediction/bonesistools/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/bonesistools.svg)](https://pypi.org/project/bonesistools)
[![python](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fbnediction%2Fbonesistools%2Frefs%2Fheads%2Fmain%2Fpyproject.toml&query=%24.project.requires-python&style=flat&label=python)](https://www.python.org/)
[![license](https://img.shields.io/pypi/l/bonesistools.svg)](https://github.com/bnediction/bonesistools/blob/main/LICENSE)

# BoNesisTools

`BoNesisTools` provides Python-implemented toolkits for upstream and downstream analyses around the [BoNesis](https://github.com/bnediction/bonesis) ecosystem.

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

---

## Single-cell tools

`bt.sct` is inspired by [Scanpy](https://github.com/scverse/scanpy) while providing additional and complementary features for single-cell analyses.

Submodules:

- preprocessing: `bt.sct.pp`
  - expression transformations, feature selection, filtering and metadata utilities

- tools: `bt.sct.tl`
  - embeddings, neighborhood graphs, clustering and differential analysis

- input/output: `bt.sct.io`
  - registered single-cell example datasets, GEO import and matrix export helpers

- plotting: `bt.sct.pl`
  - visualization helpers for embeddings, trajectories, distributions and summaries

Example:

```python
bt.sct.io.available()
bt.sct.io.info("pbmc3k")
adata = bt.sct.io.load("pbmc3k")
adata = bt.sct.io.load("nestorowa")
bt.sct.io.clear("pbmc3k")
```

---

## Boolean modelling utilities

`bt.bpy` provides utilities for Boolean modelling, logical abstractions and signed regulatory graphs.

Submodules:

- Boolean algebra: `bt.bpy.ba`
  - logical objects, configuration sets and transformations for Boolean-state reasoning

- Boolean network: `bt.bpy.bn`
  - Boolean model representation, conversion, analysis and exchange

- influence graph: `bt.bpy.ig`
  - signed regulatory graph construction, comparison, analysis and display

- input/output: `bt.bpy.io`
  - BNet, GINML, ZGINML, hypercube and influence-graph readers

Example:

```python
bn = bt.bpy.io.read_bnet("model.bnet")
graph = bn.to_influence_graph()

graph.show()
```

---

## Biological external resources

`bt.dbs` provides lightweight interfaces and utilities for biological
external resources.

Submodules:

- NCBI: `bt.dbs.ncbi`
  - gene identifier, synonym and annotation utilities

- OmniPath: `bt.dbs.omnipath`
  - regulatory interaction datasets

- HCOP: `bt.dbs.hcop`
  - orthology resources

Example:

```python
genesyn = bt.dbs.ncbi.genesyn()

grn = bt.dbs.omnipath.collectri(
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

Install the omics dependencies:

```sh
pip install "bonesistools[omics]"
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
