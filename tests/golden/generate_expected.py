#!/usr/bin/env python

from __future__ import annotations

import argparse

from _boolean_workflow import (
    EXPECTED_DIR as BOOLEAN_EXPECTED_DIR,
)
from _boolean_workflow import (
    run_boolean_workflow,
)
from _boolean_workflow import (
    save_expected as save_boolean_expected,
)
from _workflow import EXPECTED_DIR, run_workflow, save_expected


def main() -> None:

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--section",
        choices=("all", "sct", "bpy"),
        default="all",
        help="Golden output section to regenerate.",
    )
    args = parser.parse_args()

    if args.section in ["all", "sct"]:
        outputs = run_workflow()

        for path in EXPECTED_DIR.glob("*.npz"):
            path.unlink()

        save_expected(outputs)

        for path in sorted(EXPECTED_DIR.glob("*.npz")):
            print(f"wrote {path}")

    if args.section in ["all", "bpy"]:
        boolean_outputs = run_boolean_workflow()

        for path in BOOLEAN_EXPECTED_DIR.glob("*.npz"):
            path.unlink()

        save_boolean_expected(boolean_outputs)

        for path in sorted(BOOLEAN_EXPECTED_DIR.glob("*.npz")):
            print(f"wrote {path}")


if __name__ == "__main__":
    main()
