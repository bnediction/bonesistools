# Architecture

`BoNesisTools` is organized around three main domains: Boolean modelling,
single-cell analysis utilities, and biological database resources.

## Top-level package

The top-level `bonesistools` package exposes the public namespaces used by the
library:

- `bonesistools.logic`: Boolean modelling and signed graph utilities.
- `bonesistools.omics`: single-cell and multimodal analysis helpers.
- `bonesistools.resources`: interfaces to external biological resources.
- `bonesistools.grntools`: deprecated compatibility namespace.

Internal compatibility helpers are grouped in `_compat.py`.

## Boolean modelling: `logic`

`logic` contains domain objects and algorithms for Boolean modelling.

### `boolean_algebra`

Utilities for partial Boolean abstractions and hypercube-like structures.

Main responsibilities:

- Boolean and partial Boolean values.
- Hypercubes and hypercube collections.
- Exact configuration sets over fixed Boolean components.
- Boolean algebra parsing and representation helpers.
- Structural typing utilities for Boolean abstractions.

### `boolean_network`

Boolean network objects and related model utilities.

Main responsibilities:

- Boolean network representation.
- Boolean network ensembles.
- `.bnet` export.
- Rule validation and network-level operations.

### `influence_graph`

Signed regulatory and influence graph utilities.

Main responsibilities:

- Signed directed influence graphs.
- Influence graph representation.
- Signed path and walk algorithms.
- Interaction scoring from signed walks.

### `input_output`

Readers and writers for Boolean modelling files and imported logical-model
containers.

Main responsibilities:

- `.bnet` import.
- GINML, ZGINML and SBML Level 3 Qual import, plus ZGINML export.
- Hypercube and influence-graph file readers.
- Protected executable-model snapshots with Boolean networks, influence
  graphs, named initial conditions, parameters, perturbations and metadata.

## Omics tools: `omics`

`omics` provides Scanpy-like utilities for annotated single-cell and
multimodal data objects.

### `preprocessing`

Data preparation and metadata transfer helpers.

Main responsibilities:

- Gene name handling.
- Simple preprocessing operations.
- Observation transfer between integrated and condition-specific datasets.

### `tools`

Computational utilities for single-cell analyses.

Main responsibilities:

- Matrix conversion.
- Graph construction.
- Marker and log-fold-change computation.
- K-nearest-neighbor and shared-neighbor utilities.
- Regression helpers.
- Data export helpers.

### `plotting`

Plotting utilities built around matplotlib-style APIs.

Main responsibilities:

- Boxplots.
- Composition stacked barplots.
- Embedding plots.
- Density plots.
- Graph plots.
- Shared colors and figure styling.

### `datasets`

Small packaged datasets used for examples and tests.

## Biological resources: `resources`

`resources` contains loaders and interfaces for external biological resources.

### `ncbi`

NCBI Gene-derived resources.

Main responsibilities:

- Gene synonym lookup.
- Identifier conversion.
- Packaged `gene_info`-derived resources for selected organisms.

### `omnipath`

OmniPath-derived regulatory resources.

Main responsibilities:

- CollecTRI GRN loading.
- DoRothEA GRN loading.

## Design principles

The package follows a few implementation principles:

- Keep public APIs domain-oriented and readable.
- Preserve Python 3.7 runtime compatibility.
- Use typing to clarify contracts, not to obscure the implementation.
- Use targeted casts at boundaries with dynamic scientific libraries.
- Avoid propagating `Any` into internal business logic.
- Prefer explicit runtime validation for public APIs.
- Keep deprecated compatibility namespaces until the next major API cleanup.

## Public API policy

Public objects should be exposed through the closest package `__init__.py`
when they are part of the supported API. Internal modules may stay prefixed
with `_`, but their public re-exports should remain stable once documented or
used outside their module.

Deprecated aliases are kept when practical to avoid breaking user code between
minor releases. New names should be clearer than the old ones, and deprecation
messages should point to the replacement API.

## Compatibility policy

`bonesistools` currently supports Python 3.7 at runtime. Avoid global adoption
of syntax that is not valid on Python 3.7, such as PEP 604 unions
(`A | B`) or builtin generics (`list[str]`), unless annotations are safely
deferred and the local context has been checked.

Prefer readable typing over linter-driven complexity. Use domain type aliases
for recurring concepts, keep runtime validation explicit, and use targeted
casts at boundaries with dynamic scientific libraries such as pandas, AnnData,
NetworkX, or scipy.

## Private modules

Modules prefixed with `_` are implementation details. They may contain shared
helpers, internal types, parsers, or algorithmic utilities. Public objects should
be re-exported from the corresponding package `__init__.py` when they are part
of the supported API.
