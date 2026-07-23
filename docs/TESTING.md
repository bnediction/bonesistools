# Testing

## Test layout

Tests are separated by question:

- `tests/regression/`: is the software still correct?
- `tests/golden/`: does the current implementation match validated reference
  outputs?
- `tests/reproducibility/`: does the same input produce the same output?

Regression tests contain bug checks, invariants, edge cases, unit tests and
ordinary behavioral tests. They are the default test suite.

Reproducibility tests are intentionally separate because they can use larger
data, live resources or deterministic signatures.

Golden tests are intentionally separate because they compare against versioned
reference outputs and run in the pinned Python 3.13 environment under Linux,
macOS and Windows.

The guarantees, numerical-stability policy and reference-review procedure are
documented in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## CI layout

The GitHub CI separates static checks and three test matrices:

- `ruff`: linting;
- `pyright`: static typing with the `pyproject.toml` configuration;
- `reproducibility (Linux | macOS ARM64 | Windows x64)`: live resource
  checks against deterministic signatures on Python 3.13;
- `golden (Linux, strict)`;
- `golden (macOS ARM64, portable)`;
- `golden (Windows x64, portable)`;
- `import (Linux, 3.7 | 3.8 | 3.9)`;
- `test (Linux, 3.10 | 3.11 | 3.12 | 3.13)`;
- `test (macOS ARM64, 3.13)`;
- `test (Windows x64, 3.13)`.

Python 3.7 to 3.9 are import-only jobs because they mainly protect runtime
compatibility. Regression tests run on newer Python versions under Linux and on
Python 3.13 under macOS and Windows. Golden and reproducibility checks use
independent matrices so their failures are reported separately from behavioral
regressions.

## Reproducibility tests

Reproducibility tests are for external scientific resources whose current
content can affect downstream results. They intentionally run in a separate
Linux, macOS and Windows CI matrix because they may require network access and
can fail when an upstream database changes.

Signature construction and review policy are described in
[REPRODUCIBILITY.md](REPRODUCIBILITY.md#external-resources).

Run reproducibility tests explicitly with:

```bash
BONESISTOOLS_RUN_REPRODUCIBILITY=1 pytest tests/reproducibility
```

## Golden tests

Golden tests compare current outputs against validated reference outputs stored
under `tests/golden/expected/`. They are acceptance tests for important
scientific workflows, not immutable mathematical truths.

Before changing a reference, follow the numerical diagnosis and review policy
in [REPRODUCIBILITY.md](REPRODUCIBILITY.md#golden-references).

Run golden tests explicitly with:

```bash
BONESISTOOLS_RUN_GOLDEN=1 pytest tests/golden --golden-mode=strict
```

The strict contract is the default and requires bitwise equality with the
canonical Linux references, except for the `hvg_loess.score` values, which use
`rtol=0` and `atol=2e-15`; their mask, ranks and selected genes remain exact.
The portable contract keeps deterministic arrays exact while applying reviewed
tolerances or scientific quality invariants to numerically sensitive native
algorithms:

```bash
BONESISTOOLS_RUN_GOLDEN=1 pytest tests/golden --golden-mode=portable
```

CI selects the contract explicitly. It does not infer it from the operating
system: the Linux matrix entry is strict, while the macOS and Windows entries
are portable.

Regenerate expected outputs intentionally with:

```bash
BONESISTOOLS_RUN_GOLDEN=1 python tests/golden/generate_expected.py
```
