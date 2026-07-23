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

The GitHub CI separates fast local checks from live reproducibility checks:

- `ruff`: linting;
- `pyright`: static typing with the `pyproject.toml` configuration;
- `reproducibility`: live resource checks against deterministic signatures;
- `golden`: golden-output checks on Python 3.13;
- `import (3.7)`;
- `import (3.8)`;
- `import (3.9)`;
- `regression (3.10)`;
- `regression (3.11)`;
- `regression (3.12)`;
- `regression (3.13)`;
- `compatibility (macOS ARM64)`;
- `compatibility (Windows x64)`.

Python 3.7 to 3.9 are import-only jobs because they mainly protect runtime
compatibility. Regression tests run on newer Python versions under Linux and on
Python 3.13 under macOS and Windows. The macOS and Windows compatibility jobs
also run the complete golden suite in the pinned reference environment.

## Reproducibility tests

Reproducibility tests are for external scientific resources whose current
content can affect downstream results. They intentionally run in a separate CI
job because they may require network access and can fail when an upstream
database changes.

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
BONESISTOOLS_RUN_GOLDEN=1 pytest tests/golden
```

Regenerate expected outputs intentionally with:

```bash
BONESISTOOLS_RUN_GOLDEN=1 python tests/golden/generate_expected.py
```
