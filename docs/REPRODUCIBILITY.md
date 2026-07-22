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
different embedding. In the golden omics workflow, that embedding is also
consumed by KNNSC, so one UMAP divergence can produce failures in both
reference outputs.

The reproducible UMAP path therefore:

- canonicalizes self-neighbor distances to exact zero before estimating local
  connectivity scales;
- uses the supplied seed throughout initialization and optimization;
- uses random initialization by default, avoiding platform-dependent numerical
  eigensolvers;
- uses serial, strict floating-point arithmetic for the Euclidean optimization
  path rather than relaxed `fastmath` transformations;
- limits native numerical thread pools according to `n_jobs`, whose default is
  one.

These controls select a stable representative embedding. They do not imply
that UMAP has a unique mathematical embedding. UMAP's own
[reproducibility guide](https://umap-learn.readthedocs.io/en/latest/reproducibility.html)
explains the role of seeds and serial execution; bonesistools additionally
stabilizes numerical boundaries needed for cross-runner golden tests.

Spectral initialization remains available explicitly. Bonesistools converts
its coordinates to the embedding precision and fixes arbitrary eigenvector
orientations before optimization. This reduces backend variability, but the
underlying numerical eigensolver can still produce mathematically equivalent
yet non-bitwise initial coordinates across platforms.

For a UMAP-related golden failure, compare outputs in this order:

1. PCA coordinates;
2. neighbor indices and distances;
3. fuzzy connectivities;
4. initialization;
5. final UMAP coordinates;
6. analyses that consume the embedding, such as KNNSC.

Only the first divergent stage should drive the correction.

## Golden references

Golden tests compare current outputs against validated references under
`tests/golden/expected/`. They are acceptance tests for important scientific
workflows, not immutable mathematical truths.

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
