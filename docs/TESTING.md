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
reference outputs and run only in the canonical Python 3.13 CI job.

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
- `regression (3.13)`.

Python 3.7 to 3.9 are import-only jobs because they mainly protect runtime
compatibility. Full tests run on newer Python versions.

## Reproducibility tests

Reproducibility tests are for external scientific resources whose current
content can affect downstream results. They intentionally run in a separate CI
job because they may require network access and can fail when an upstream
database changes.

The tests do not redistribute restricted upstream data. Instead, they:

1. fetch the resource at test time;
2. build the corresponding bonesistools object;
3. compute a deterministic signature from the result;
4. hash that signature with SHA-256;
5. compare the hash to a tracked `.sha256` file.

For DoRothEA, the tracked files are:

- `tests/reproducibility/dorothea_current_mouse_A.sha256`;
- `tests/reproducibility/dorothea_legacy_mouse.sha256`.

These hashes are small reproducibility sentinels. If one changes, the test
failure means either:

- the upstream resource changed;
- the resource access or post-processing logic changed;
- the expected signature must be reviewed and intentionally updated.

The local JSON files in `tests/reproducibility/` are diagnostic artifacts only.
They must not be modified for normal test updates and must not be committed.

Run reproducibility tests explicitly with:

```bash
BONESISTOOLS_RUN_REPRODUCIBILITY=1 pytest tests/reproducibility
```

## Golden tests

Golden tests compare current outputs against validated reference outputs stored
under `tests/golden/expected/`. They are acceptance tests for important
scientific workflows, not immutable mathematical truths.

If a golden test fails, first determine whether the difference is an unintended
code regression, an intentional algorithmic change, or a dependency/numerical
environment change. Only intentional algorithmic changes and validated
environment changes justify regenerating the expected outputs.

Run golden tests explicitly with:

```bash
BONESISTOOLS_RUN_GOLDEN=1 pytest tests/golden
```

Regenerate expected outputs intentionally with:

```bash
BONESISTOOLS_RUN_GOLDEN=1 python tests/golden/generate_expected.py
```

This section is the golden update procedure.
