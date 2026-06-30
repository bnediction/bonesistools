# Testing

## Test layout

Tests are separated by question:

- `tests/regression/`: is the software still correct?
- `tests/reproducibility/`: does the same input produce the same output?

Regression tests contain bug checks, invariants, edge cases, unit tests and
ordinary behavioral tests. They are the default test suite.

Reproducibility tests are intentionally separate because they can use larger
data, live resources or deterministic signatures.

## CI layout

The GitHub CI separates fast local checks from live reproducibility checks:

- `ruff`: linting;
- `pyright`: static typing with the `pyproject.toml` configuration;
- `reproducibility`: live resource checks against deterministic signatures;
- `pytest (3.7, import)`;
- `pytest (3.8, import)`;
- `pytest (3.9, import)`;
- `pytest (3.10, test)`;
- `pytest (3.11, test)`;
- `pytest (3.12, test)`;
- `pytest (3.13, test)`.

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
