#!/usr/bin/env python

from __future__ import annotations

import argparse
import sys
from pathlib import Path

GOLDEN_DIR = Path(__file__).parent
sys.path.insert(0, str(GOLDEN_DIR))

from logic._boolean_workflow import (  # noqa: E402
    EXPECTED_DIR as BOOLEAN_EXPECTED_DIR,
)
from logic._boolean_workflow import (  # noqa: E402
    run_boolean_workflow,
)
from logic._boolean_workflow import (  # noqa: E402
    save_expected as save_boolean_expected,
)
from omics._omics_workflow import (  # noqa: E402
    EXPECTED_DIR,
    run_omics_workflow,
    save_expected,
)


def main() -> None:

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--section",
        choices=("all", "omics", "logic"),
        default="all",
        help="Golden output section to regenerate.",
    )
    args = parser.parse_args()

    if args.section in ["all", "omics"]:
        outputs = run_omics_workflow()

        for path in EXPECTED_DIR.glob("*.npz"):
            path.unlink()

        save_expected(outputs)

        for path in sorted(EXPECTED_DIR.glob("*.npz")):
            print(f"wrote {path}")

    if args.section in ["all", "logic"]:
        boolean_outputs = run_boolean_workflow()

        for path in BOOLEAN_EXPECTED_DIR.glob("*.npz"):
            path.unlink()

        save_boolean_expected(boolean_outputs)

        for path in sorted(BOOLEAN_EXPECTED_DIR.glob("*.npz")):
            print(f"wrote {path}")


if __name__ == "__main__":
    main()
