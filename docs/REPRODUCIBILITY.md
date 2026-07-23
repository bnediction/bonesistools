# Reproducibility

This document defines the reproducibility guarantees expected from
bonesistools and the procedure for reviewing reference changes. Test layout
and execution commands are documented in [TESTING.md](TESTING.md).

## Scope

Reproducibility has several distinct levels:

- seeded reproducibility controls stochastic choices for identical inputs;
- numerical reproducibility prevents supported execution environments from
  selecting different floating-point results;
- golden reproducibility compares complete scientific workflows with validated
  reference outputs;
- resource reproducibility detects changes in external scientific databases.

Passing one level does not imply the others. A fixed seed does not control CPU
instructions, BLAS kernels, reduction order, native thread scheduling or JIT
compiler optimizations. Likewise, pinned Python dependencies constrain the
software environment but do not make floating-point arithmetic independent of
the machine running it.

## Numerical reproducibility

Numerically sensitive code should be deterministic by construction. Do not
rely on CI-specific environment variables or a forced CPU target as part of the
public behavior. Instead:

1. canonicalize quantities with an exact mathematical value, such as a sample's
   distance to itself;
2. use explicit random states for stochastic algorithms;
3. control parallel execution when operation ordering affects results;
4. avoid relaxed floating-point transformations on reproducible code paths;
5. use convergence criteria that prevent solver noise from changing a
   downstream branch or initialization;
6. canonicalize mathematically equivalent representations when their arbitrary
   orientation or ordering would be amplified later.

When a numerical reference fails, locate the first divergent intermediate
result. Do not start from the final embedding or regenerate every downstream
reference. A small upstream difference can be amplified by a nonlinear method
while the downstream implementation remains correct.

Diagnostic runs may vary CPU and numerical backends, for example through
`NUMBA_CPU_NAME` or `OPENBLAS_CORETYPE`. These variables are diagnostic tools,
not required runtime settings. The implementation should produce its reference
result without forcing either variable.

### UMAP

UMAP combines several numerically sensitive stages:

1. construction of the nearest-neighbor graph;
2. estimation of local scales from neighbor distances;
3. initialization;
4. stochastic nonlinear optimization.

An epsilon-valued self-distance can change the local connectivity scale, and a
small initialization difference can lead the nonlinear optimizer to a
different embedding. KNNSC is evaluated independently from a PCA
representation so that a UMAP divergence is not reported as a second,
downstream KNNSC regression.

The reproducible UMAP path therefore:

- canonicalizes self-neighbor distances to exact zero before estimating local
  connectivity scales;
- evaluates the fuzzy-connectivity exponential for `float32` inputs through a
  fixed polynomial and accumulates probabilities in a fixed order, avoiding
  platform-specific vector math near the local-scale convergence threshold;
- canonicalizes automatically estimated curve parameters before optimization;
- uses the supplied seed throughout initialization and optimization;
- uses the canonical spectral initialization by default, converts its
  coordinates to `float32` before scaling and fixes arbitrary eigenvector
  orientations;
- prevents fused multiply-add contraction when accumulating squared Euclidean
  distances in the optimizer;
- limits native numerical thread pools according to `n_jobs`, whose default is
  one.

These controls make the seeded serial optimization path bitwise reproducible
across the processor targets covered by the test suite. They do not imply that
UMAP has a unique mathematical embedding, and external eigensolvers may still
change a spectral initialization when their numerical behavior changes.
UMAP's own
[reproducibility guide](https://umap-learn.readthedocs.io/en/latest/reproducibility.html)
explains the role of seeds and serial execution.

Spectral initialization remains available explicitly. Bonesistools fixes
arbitrary eigenvector orientations and converts the eigensolver output to
`float32` before UMAP scales it. This prevents insignificant platform-level
`float64` differences from influencing the scale supplied to the nonlinear
optimizer.

The dedicated UMAP golden decomposes this trajectory into exact checkpoints:

1. curve parameters;
2. prepared fuzzy graph;
3. raw and oriented spectral layouts;
4. noisy and normalized initialization;
5. optimizer random state;
6. embeddings at epochs 0, 1, 2, 5, 10, 25, 50, 100, 250, and 500.

The test compares these checkpoints in order and stops at the first
divergence. This diagnostic isolated two former platform-dependent operations:
SciPy's numerical fit of the UMAP curve parameters and fused multiply-add
contraction in UMAP's squared Euclidean distance. Bonesistools canonicalizes
the fitted parameters and excludes only this contraction from UMAP's otherwise
optimized distance calculation. PCA, nearest-neighbor search, spectral
initialization and random-state construction retain their original numerical
paths.

For a UMAP-related golden failure, compare outputs in this order:

1. PCA coordinates;
2. neighbor indices and distances;
3. fuzzy connectivities;
4. spectral and normalized initialization;
5. optimizer checkpoints;
6. final UMAP coordinates.

Only the first divergent stage should drive the correction.

### t-SNE

The seeded t-SNE path uses random initialization explicitly. Scikit-learn's
implicit PCA initialization depends on the installed linear algebra backend;
numerically equivalent initial coordinates can then lead the nonlinear
optimizer to different embeddings. Explicit random initialization makes the
starting coordinates a direct function of `seed` while preserving the same
distance, probability and optimization calculations.

The `n_jobs` argument also limits the native numerical thread pools used by
the Barnes-Hut optimizer. Scikit-learn otherwise applies `n_jobs` to nearest-
neighbor search but may select a different OpenMP thread count for gradient
descent. The default `n_jobs=1` therefore defines a serial seeded path for the
complete t-SNE calculation.

The t-SNE golden receives the frozen PCA reference rather than the PCA
recomputed during the same test run. PCA remains tested independently, while a
t-SNE coordinate divergence can be attributed to the t-SNE stage itself.

The dedicated t-SNE golden records exact checkpoints for:

1. the nearest-neighbor structure and squared distances;
2. conditional neighborhood probabilities;
3. the symmetrized joint-probability graph;
4. seeded random initialization and its initial gradient;
5. embeddings at iterations 1, 2, 5, 10, 25, 50, 100, 250, and 300.

Neighbor indices and every downstream checkpoint are compared exactly. Raw
squared distances use a narrow floating-point tolerance because equivalent
serial distance kernels can differ below `1e-12`; the subsequent probability
checkpoint remains exact and detects whether such noise affects t-SNE itself.

The diagnostic also verifies that its final checkpoint is exactly the result
returned by the public scikit-learn solver under the same serial settings. A
platform failure therefore identifies whether divergence first occurs during
neighbor search, probability estimation, initialization, or optimization.

## Golden references

Golden tests compare current outputs against validated references under
`tests/golden/expected/`. They are acceptance tests for important scientific
workflows, not immutable mathematical truths.

The omics golden isolates dependencies between tested stages. A stage is
collected and compared independently, then its frozen reference is supplied to
downstream stages where available. PCA therefore consumes the frozen HVG mask;
neighbors, spectral embedding, KNNSC and t-SNE consume the frozen PCA; and
Louvain and UMAP consume the frozen neighbor graph. Intentional reference
generation disables this isolation and rebuilds the complete workflow.

If a golden test fails, first determine whether the difference is:

- an unintended code regression;
- an intentional algorithmic change;
- a dependency change;
- a machine-dependent numerical divergence.

Do not update a golden file merely to make one runner pass. Correct
machine-dependent behavior at its source and verify the result in the canonical
Python 3.13 environment. Regenerate a reference only after an intentional
behavior change or a reviewed change to the canonical environment.

The update command is:

```bash
BONESISTOOLS_RUN_GOLDEN=1 python tests/golden/generate_expected.py
```

Review the generated diff and rerun the complete golden suite before accepting
it.

## External resources

Tests under `tests/reproducibility/` monitor external scientific resources whose
content can affect downstream results. They may require network access and can
fail when an upstream database changes.

The tests do not redistribute restricted upstream data. Instead, they:

1. fetch the resource at test time;
2. build the corresponding bonesistools object;
3. compute a deterministic signature from the result;
4. hash that signature with SHA-256;
5. compare the hash with a tracked `.sha256` file.

For DoRothEA, the tracked references are:

- `tests/reproducibility/dorothea_current_mouse_A.sha256`;
- `tests/reproducibility/dorothea_legacy_mouse.sha256`.

A changed signature can indicate an upstream resource change, a change in
resource processing, or an intentional reference migration. Investigate which
case applies before updating the hash.

Local JSON files under `tests/reproducibility/` are diagnostic artifacts. They
must not be modified during ordinary test updates or committed.
