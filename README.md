[![tests](https://github.com/bnediction/bonesistools/actions/workflows/tests.yml/badge.svg)](https://github.com/bnediction/bonesistools/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/bonesistools.svg)](https://pypi.org/project/bonesistools)
[![python](https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fbnediction%2Fbonesistools%2Frefs%2Fheads%2Fmain%2Fpyproject.toml&query=%24.project.requires-python&style=flat&label=python)](https://www.python.org/)
[![license](https://img.shields.io/pypi/l/bonesistools.svg)](https://github.com/bnediction/bonesistools/blob/main/LICENSE)

# BoNesisTools

`BoNesisTools` is a Python package providing bioinformatics utilities for upstream and downstream analyses of the [BoNesis](https://github.com/bnediction/bonesis) framework.

## Usage

```python
import bonesistools as bt
```

`BoNesisTools` exposes four main namespaces:

- `bt.sct` — single-cell and multimodal annotated data tools
- `bt.bpy` — Boolean modelling utilities
- `bt.dbs` — biological database interfaces
- `bt.grn` — deprecated alias for `bt.bpy.ig`

`bt.sct` follows a [Scanpy](https://github.com/scverse/scanpy)-like API while providing additional and complementary features for single-cell analyses:

- preprocessing: `bt.sct.pp`
- tools: `bt.sct.tl`
- plotting: `bt.sct.pl`
- datasets: `bt.sct.datasets`

`bt.bpy` provides utilities for Boolean modelling:

- Boolean algebra and partial Boolean abstractions: `bt.bpy.ba`
- Boolean network utilities: `bt.bpy.bn`
- interaction and influence graph utilities: `bt.bpy.ig`

## Installation

Install the latest release:
```sh
pip install bonesistools
```

Install the single-cell tools dependencies:
```sh
pip install "bonesistools[sctools]"
```

Install with all extra dependencies:
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

## Bugs

Please report any bugs or ask questions [here](https://github.com/bnediction/bonesistools/issues).

## License

This package is distributed under the [CeCILL v2.1](http://www.cecill.info/index.en.html) free software license (GNU GPL compatible).
