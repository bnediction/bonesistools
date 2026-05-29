# Testing

## Philosophy

Tests should verify expected behaviour, not only execution. A good test should
make it clear what the expected result is and why this result is mathematically,
biologically, or structurally meaningful.

Prefer small synthetic examples with known expected results over large real
datasets. Mathematical and graph-based utilities should be tested against
independently derived expectations.

## General Guidelines

- Avoid smoke tests when an expected result can be computed.
- Use small deterministic inputs.
- Prefer synthetic matrices, graphs, hypercubes and AnnData objects.
- Use fake classes to isolate external dependencies when appropriate.
- Do not call external services in unit tests.
- Keep tests fast and reproducible.
- If a bug is discovered, write a minimal failing test before changing the
  implementation.

## Expected Values

Expected values should come from reasoning independent of the implementation.
Avoid reimplementing the tested function inside the test.

Good patterns include:

- constructing a matrix from known coefficients and checking that regression
  recovers the expected residual structure;
- building a signed graph by hand and checking exact signed paths or scores;
- creating small point clouds where central and peripheral cells are visually
  and geometrically obvious;
- using tiny Boolean networks where fixed points can be enumerated by hand.

## Coverage Policy

Coverage should increase by testing meaningful behaviours, not by exercising
lines mechanically.

Prioritize missing branches when they correspond to:

- input validation;
- non-trivial numerical or graph logic;
- documented public behaviour;
- compatibility paths that users are likely to hit.

Avoid coverage-only tests for fragile or low-value branches, such as
multiprocessing fallbacks or legacy dependency fallbacks, unless the behaviour
is important enough to justify the maintenance cost.

## Warnings

Warnings are part of the tested behaviour.

- Use `pytest.warns(...)` when a warning is expected.
- Let unexpected warnings appear during normal test runs.
- Suppress warnings only locally, and only when they come from unavoidable
  third-party setup noise.
- Do not suppress warnings around the behaviour being tested.

For example, if a function is designed to handle duplicated AnnData variable
names, the test should verify that the function itself does not leak the
duplicated-name warning.

## External Resources

Unit tests should not depend on remote services, changing databases, or network
availability.

Use fake loaders, monkeypatching, tiny embedded tables, or packaged cache files
when testing database helpers. Integration tests that require external
resources should be clearly separated from unit tests.

## Static Checks

The project uses Ruff and Pyright as development checks.

Ruff diagnostics should be fixed when they reveal clearer Python or real bugs.
Pyright diagnostics should be fixed when the fix improves the code while
preserving Python 3.7 runtime compatibility.

Do not change public APIs, make annotations less precise, or add noisy
workarounds only to satisfy a checker.

## Targeted Test Runs

Targeted pytest runs are useful while developing, but they can behave
differently from the full suite when imports, plugins, or coverage are involved.

When validating a change, prefer running the full relevant test file, and use
coverage on broader runs rather than very narrow single-test invocations.

## Examples of Preferred Tests

- Boolean algebra: compare exact truth values or partial assignments.
- Boolean networks: test fixed points, predecessors and rule semantics.
- Influence graphs: test signed paths, SCCs, feedback circuits and walk scores
  on small hand-built graphs.
- KNNBS: use simple geometric configurations with visually predictable central
  and peripheral cells.
- Database helpers: use fake loaders or small embedded tables rather than live
  remote resources.
