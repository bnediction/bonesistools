#!/usr/bin/env python

from __future__ import annotations

from _workflow import EXPECTED_DIR, run_workflow, save_expected


def main() -> None:

    outputs = run_workflow()
    for path in EXPECTED_DIR.glob("*.npz"):
        path.unlink()

    save_expected(outputs)

    for path in sorted(EXPECTED_DIR.glob("*.npz")):
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
