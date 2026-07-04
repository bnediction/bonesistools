# Development guidelines

## Internal helper modules

The project provides a small set of shared internal modules located at the
package root.

### `_typing.py`

Shared internal type aliases used across multiple subpackages.

Examples:

* `DataFrameAxis`
* `FileOrientation`
* `RandomStateSeed`

Only place types here when they are reused by multiple modules. Keep highly
local types close to their implementation.

AnnData-specific aliases such as `AnnDataAxis` and `AnnDataAxisWithBoth` live
in `sctools/_typing.py`.

Avoid creating aliases for short, readable type annotations in public APIs.
Prefer writing the annotation directly so users can understand the signature
from `help(...)` without looking up an alias elsewhere:

```python
legend: Union[bool, Mapping[str, Any]]
```

Only introduce a public-facing alias when the type is long enough to obscure
the signature or when the alias names a real domain concept.

### `_validation.py`

Centralized argument validation and normalization helpers.

Use `_as_*` helpers whenever a value is validated and converted to a canonical
representation.

Examples:

* `_as_positive_integer`
* `_as_non_negative_integer`
* `_as_positive_number`
* `_as_non_negative_number`
* `_as_boolean`
* `_as_callable`
* `_as_string`
* `_as_literal`
* `_as_memory_size`
* `_as_seed`
* `_as_dataframe_axis`

Avoid duplicating validation logic in public functions when a reusable helper
already exists.

AnnData-specific validation helpers such as `_as_anndata_axis` live in
`sctools/_validation.py`.

Do not centralize checks that are already handled cleanly by the underlying
library. For example, prefer using `adata.obsm[key]`, `adata.obs[key]` or
`adata.layers[key]` directly and let AnnData or pandas raise the natural
`KeyError`, unless bonesistools needs to add domain-specific behavior.

Shared single-cell statistical helpers, such as expression mean and variance
calculations, live in `sctools/_stats.py`.

### `_warnings.py`

Centralized warning and deprecation helpers.

Use this module whenever a warning follows a common project pattern.

Deprecation warnings should:

* use `FutureWarning`;
* identify the deprecated object;
* provide the recommended replacement when available;
* indicate the planned removal version;
* use an appropriate `stacklevel`.

Avoid repeating raw `warnings.warn(...)` blocks when a shared helper is
available.

### `_compat.py`

Compatibility helpers and transitional APIs.

Typical use cases include:

* deprecated aliases;
* compatibility shims;
* version-dependent behavior;
* optional dependencies.

Keep compatibility code isolated whenever possible.

---

## Validation conventions

### AnnData axes

Public AnnData APIs should use:

```python
axis="obs"
axis="var"
```

Integer axes (`0`, `1`) may be accepted temporarily for backward compatibility
but should not be introduced in new APIs.

### DataFrame axes

DataFrame-related helpers may accept:

```python
axis=0
axis=1
axis="index"
axis="columns"
```

and normalize them to a canonical representation.

### File orientation

When describing how objects are stored in tabular files, prefer:

```python
orientation="rows"
orientation="columns"
```

instead of `axis`.

---

## Dataset registry

Built-in single-cell datasets are registered in `bt.sct.datasets` through a
small Python loader registry and a JSON metadata file:

* `_DATASET_LOADERS`: maps dataset names to loader functions;
* `datasets.json`: stores public metadata only.

Keep `datasets.json` descriptive. Do not store loader functions or
implementation objects in it.

Each dataset entry should define:

* `description`: concise dataset description;
* `organism`: scientific organism name, for example `"Mus musculus"`;
* `tissue`: tissue or sample source in lower case;
* `technology`: assay or platform name;
* observation count fields, usually `cells`;
* feature count fields with explicit biological names, such as `genes` or
  `peaks`;
* `source`: source provider or accession;
* `license`: license string, or `"Not specified"` when no license is provided;
* `url`: a single canonical user-facing URL with dataset metadata, not
  necessarily the technical download URL;
* `citation`: publication or dataset citation, or `"Not specified"` when no
  citation is provided.

Use concrete feature names instead of a generic `features` field. For example,
a scRNA-seq dataset should use `genes`, while a chromatin accessibility dataset
may use `peaks`. Multimodal datasets may define several feature-count fields,
such as both `genes` and `peaks`.

Example:

```json
{
    "nestorowa": {
        "description": "Single-cell RNA sequencing dataset of mouse hematopoietic stem and progenitor cells (HSPCs).",
        "organism": "Mus musculus",
        "tissue": "bone marrow",
        "technology": "Smart-seq2",
        "cells": 1656,
        "genes": 46078,
        "source": "NCBI GEO (GSE81682)",
        "license": "Not specified",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81682",
        "citation": "Nestorowa et al. (2016). A single-cell resolution map of mouse hematopoietic stem and progenitor cell differentiation. Blood, 128(8), e20-e31."
    }
}
```

---

## Public parameter order

Public functions should use a consistent parameter hierarchy when it naturally
fits the API:

1. object being operated on;
2. hyperparameters and behavior options;
3. metric parameters;
4. where results are stored, such as `key_added`, `distances_key` or
   `connectivities_key`;
5. `seed`;
6. `n_jobs`;
7. `copy`.

Do not reorder parameters mechanically when it would create noisy diffs or
break a stable public API without a clear benefit.

---

## Plotting APIs

### Legends

Plotting functions should expose one legend parameter:

```python
legend=True
```

Use:

* `legend=False` to disable the legend;
* `legend=True` to draw it with default Matplotlib parameters;
* `legend={...}` or any other mapping to forward keyword arguments to
  `Axes.legend`.

Avoid adding separate `show_legend` parameters. During migrations, deprecated
`show_legend`/`showlegend` compatibility keywords must be ignored after
emitting a warning; `legend` remains the only authoritative parameter.

For plotting elements that can be both enabled and styled, prefer the same
polymorphic pattern:

```python
median=True
mean=False
```

Use `False` to disable the element, `True` to draw it with defaults, and a
mapping to draw it with forwarded artist properties. Avoid adding paired
`show_*` parameters when a single domain-level parameter can express both
activation and customization.

Document the forwarding target explicitly for each mapping form. For example:

* `legend={...}` is forwarded to `Axes.legend`;
* `labels={...}` is forwarded to `Axes.text`;
* `points={...}` is forwarded to `Axes.scatter`;
* `median={...}` is forwarded to `medianprops` for boxplots and applied to
  `cmedians` for violins;
* `mean={...}` is forwarded to `meanprops` for boxplots and applied to
  `cmeans` for violins.

---

## Local helpers

Prefer keeping logic inline when a private helper would be used only once and
does not clarify a genuinely complex block.

Create a private helper when it:

* is reused by multiple functions;
* isolates a coherent algorithmic step;
* improves readability of a long public function;
* provides a stable boundary around external package typing or compatibility.

Avoid adding private helpers only to give a name to two or three obvious lines
of code.

---

## Expression, representation and pairwise matrices

Use precise naming for AnnData expression matrices and observation
representations:

* `expression` is a public parameter naming an expression source:
  `.X`, `.raw.X` or `.layers[...]`;
* `representation` is a public parameter naming an observation representation
  in `.obsm[...]`, such as `"X_pca"` or `"X_umap"`;
* `pairwise` is a public parameter naming a pairwise matrix in `.obsp[...]` or
  `.varp[...]`;
* `expression_mtx` is a local variable holding the matrix-like value returned
  by `get_expression(...)`;
* `representation_mtx` is a local variable holding the matrix-like value
  returned by `get_representation(...)`;
* `pairwise_mtx` is a local variable holding the matrix-like value returned by
  `get_pairwise(...)`.

Prefer `_mtx` over `_matrix` or `_array` for these values because they may be
dense NumPy arrays or sparse matrices.

Avoid reusing public selector names for computed matrix-like values:

```python
expression_mtx = get_expression(...)
representation_mtx = get_representation(...)
pairwise_mtx = get_pairwise(...)
```

Avoid:

```python
expression = get_expression(...)
representation = get_representation(...)
pairwise = get_pairwise(...)
X = get_expression(...)
X = get_representation(...)
```

---

## Variable subsets

Use `var_subset` for generic variable selection in AnnData tools.

Prefer:

```python
pca(
    var_subset="highly_variable",
)

pca(
    var_subset=["Gata1", "Klf1", "Tal1"],
)
```

over specialized options such as:

```python
pca(
    use_highly_variable=True,
)
```

The `var_subset` argument represents a generic variable-selection mechanism:

* `var_subset=None`: use all variables;
* `var_subset="highly_variable"`: use the boolean column
  `adata.var["highly_variable"]`;
* `var_subset=["Gata1", "Klf1", "Tal1"]`: use explicit variable names.

Keep validation centralized and design the structure so future `obs_subset`
support can follow the same pattern.

---

## Metric column names

Use stable metric prefixes in AnnData annotations:

* `n_`: cardinality;
* `total` and `total_`: sum;
* `mean`: mean;
* `median`: median;
* `variance`: variance;
* `mad`: median absolute deviation;
* `pct_`: percentage;
* `log1p_`: log1p transformation.

Prefer full statistical names over abbreviations. For example, use
`variance` and `log1p_variance`, not `var` and `log1p_var`.

---

## Fitted objects

Objects that learn state during execution should distinguish between:

### Configuration

Parameters used to fit the object.

Example:

```python
params_
```

### Learned state

Objects generated during fitting.

Examples:

```python
knn_graph
shortest_path_lengths_df
cluster_counts
```

Whenever possible:

* store internal state in private attributes (`_attribute`);
* expose public read-only properties when appropriate.

---

## AnnData return documentation

For functions that mutate AnnData or MuData objects, the `Returns` section
should describe the concrete fields that are added or updated.

Prefer:

```text
If `copy=True`, returns a copy of `adata` with PCA results added.
Otherwise, updates `adata` in place and returns None.

PCA results are stored in:

- `adata.obsm[key_added]`: projected cell coordinates;
- `adata.varm["PCs"]`: principal component loadings;
- `adata.uns["pca"]`: PCA metadata.
```

Avoid only saying that the function updates the object in place. Mention
specific fields such as `adata.obs[key]`, `adata.var[key]`,
`adata.obsm[key_added]`, `adata.uns[key_added]`, `adata.layers[layer]` or
`adata.obsp[key_added]` whenever they are part of the public behavior.

---

## Docstring section order

Use the following section order for docstrings:

1. description;
2. `Parameters`;
3. `Examples`;
4. `Returns`;
5. `Raises`;
6. `Notes`;
7. `References`;

Omit sections that do not apply. Keep `Parameters` before `Examples` so the
call contract is introduced before usage snippets. Use `Notes` for secondary
implementation details, caveats, or reproducibility comments that are not
required to understand the call signature. Use `References` for papers,
specifications, external documentation or other citable sources.

Format references as compact prose entries:

```text
Author et al. (year). Title. Journal, volume(issue), pages.
```

Use `et al.` for references with three or more authors. For one or two
authors, cite the listed author names directly:

```text
Author (year). Title. Journal, volume(issue), pages.
Author and Author (year). Title. Journal, volume(issue), pages.
Author et al. (year). Title. Journal, volume(issue), pages.
```

Prefer one reference per paragraph. Wrap long references only to keep line
length readable.

Use numbered bracket references (`[1]`, `[2]`, etc.) only for specific
citations inside narrative text:

```text
This can be computed using XXX algorithm [1] or YYY algorithm [2].
```

Do not number entries in the `References` section itself. A reference that
describes the function or method as a whole should appear there without
brackets. Use bracketed references only when a specific statement in the
docstring points to a specific source. Number bracket references by order of
appearance in the narrative text.

Keep `Raises` sections for errors that users cannot trivially infer from the
signature or ordinary runtime validation:

* domain errors, such as unsupported organisms or overlapping cluster sets;
* important side-effect decisions, such as `NotImplementedError` when `.varp`
  is present and `varp="error"`;
* non-obvious required external state.

Avoid documenting routine type validation such as "`copy` is not a bool" or
"`genesyn` has an unsupported type" unless that validation is part of a
non-obvious API contract.

---

## Class organization

For classes, prefer this order when it naturally fits the implementation:

1. class attributes
2. `__*__` methods
3. public constructors / conversions
4. public mutation methods
5. public query methods
6. public analysis methods
7. visualization / export methods
8. private helper methods
9. static/class methods

Do not reorganize large classes mechanically if it makes diffs noisy.

---

## Naming conventions

Prefer concept-oriented names over implementation-oriented names.

Examples:

```python
distribution
embedding
composition
```

instead of names tied to a specific plotting backend or algorithm.

Public APIs should describe what the user manipulates rather than how the
result is computed.

---

## Typing conventions

Keep runtime compatibility with Python 3.7.

Use Python 3.7-compatible typing syntax:

```python
Optional[str]
Union[str, Path]
Dict[str, Any]
List[str]
```

Do not use PEP 604 or builtin generic syntax while Python 3.7 is supported:

```python
str | None
dict[str, Any]
list[str]
```

Use targeted `cast(...)` at boundaries with dynamic or poorly typed libraries,
including pandas, AnnData, MuData, NetworkX, scipy, matplotlib and sklearn.

Runtime validation may be kept even when static typing considers it redundant,
especially for public APIs.

Prefer domain-specific internal aliases over propagating `Any` through internal
logic.

### `copy` overloads

Public functions with the standard AnnData-style `copy` behavior should expose
overloads when the return type depends on `copy`.

Use this pattern:

```python
@overload
def function(..., copy: Literal[True] = True) -> AnnData:
    ...

@overload
def function(..., *, copy: Literal[False]) -> None:
    ...

@overload
def function(..., copy: bool = True) -> Optional[AnnData]:
    ...
```

Make the `copy=False` overload keyword-only when omitting `copy` would
otherwise overlap with the `copy=True` default overload. Decorators should be
placed on the implementation, not on overload stubs.

The same convention applies to plotting functions with `outfile`: calls without
`outfile` should be typed as returning Matplotlib objects, while calls with
`outfile` should be typed as returning None.

---

## Linting and type checking

Lint diagnostics that reveal real bugs or clearer Python should be fixed.

Avoid broad rewrites solely for linting or type-checking tools.

Use the project CI configuration when validating type changes.

Do not add global ignores. Use local ignores only with a short reason.

---

## Tests

Prefer full test files or the standard test suite over highly targeted
pytest invocations with coverage.

Tests are split by CI policy and purpose:

* `tests/regression`: ordinary correctness tests run by default and across the
  supported Python matrix;
* `tests/golden`: golden-output tests run only in a canonical environment with
  Python 3.13 and `BONESISTOOLS_RUN_GOLDEN=1`;
* `tests/reproducibility`: reproducibility tests run only in the dedicated
  reproducibility job.

For computational functions, use small synthetic inputs with known expected
results. Expected values should come from reasoning independent of the
implementation.

Golden files are frozen test resources, not tests. Store them under
`tests/golden` and update them only when an intentional algorithmic change or
bug fix changes the expected behavior. See `docs/TESTING.md` for the golden
acceptance-test update procedure.

When testing warnings or errors, prefer checking the warning/error type and the
behavior. Avoid asserting exact message strings unless the message itself is
part of the public behavior being tested.

Do not update reproducibility references unless explicitly requested. If code no
longer matches a reference, fix the implementation first.

---

## Deprecation policy

When renaming public APIs:

1. Introduce the new API.
2. Keep the previous name as a deprecated alias.
3. Emit a `FutureWarning`.
4. Remove deprecated aliases only in a major release.

Prefer long deprecation windows when practical.
