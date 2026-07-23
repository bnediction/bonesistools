#!/usr/bin/env python

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path
from typing import Dict, cast

import numpy as np
from typing_extensions import Protocol


class _TSNEDiagnosticsModule(Protocol):
    EXPECTED_PATH: Path

    def run_tsne_diagnostics(self) -> Dict[str, np.ndarray]: ...

    def save_expected(
        self,
        checkpoints: Dict[str, np.ndarray],
        output_path: Path = ...,
    ) -> None: ...

GOLDEN_DIR = Path(__file__).parent
PROJECT_DIR = GOLDEN_DIR.parents[1]
sys.path.insert(0, str(PROJECT_DIR))

from tests.golden.logic import _graph_layout  # noqa: E402
from tests.golden.logic._boolean_workflow import (  # noqa: E402
    EXPECTED_DIR as BOOLEAN_EXPECTED_DIR,
)
from tests.golden.logic._boolean_workflow import (  # noqa: E402
    run_boolean_workflow,
)
from tests.golden.logic._boolean_workflow import (  # noqa: E402
    save_expected as save_boolean_expected,
)
from tests.golden.omics._omics_workflow import (  # noqa: E402
    EXPECTED_DIR,
    run_omics_workflow,
    save_expected,
)
from tests.golden.omics._umap_diagnostics import (  # noqa: E402
    EXPECTED_PATH as UMAP_DIAGNOSTIC_PATH,
)
from tests.golden.omics._umap_diagnostics import (  # noqa: E402
    run_umap_diagnostics,
)
from tests.golden.omics._umap_diagnostics import (  # noqa: E402
    save_expected as save_umap_diagnostics,
)

_tsne_diagnostics = cast(
    _TSNEDiagnosticsModule,
    importlib.import_module("tests.golden.omics._tsne_diagnostics"),
)
TSNE_DIAGNOSTIC_PATH = _tsne_diagnostics.EXPECTED_PATH


def main() -> None:

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--section",
        choices=(
            "all",
            "omics",
            "logic",
            "tsne-diagnostics",
            "umap-diagnostics",
        ),
        default="all",
        help="Golden output section to regenerate.",
    )
    args = parser.parse_args()

    if args.section in ["all", "omics"]:
        outputs = run_omics_workflow(use_expected_dependencies=False)

        for path in EXPECTED_DIR.glob("*.npz"):
            path.unlink()

        save_expected(outputs)

        for path in sorted(EXPECTED_DIR.glob("*.npz")):
            print(f"wrote {path}")

    if args.section in ["all", "omics", "umap-diagnostics"]:
        save_umap_diagnostics(run_umap_diagnostics())
        print(f"wrote {UMAP_DIAGNOSTIC_PATH}")

    if args.section in ["all", "omics", "tsne-diagnostics"]:
        _tsne_diagnostics.save_expected(
            _tsne_diagnostics.run_tsne_diagnostics()
        )
        print(f"wrote {TSNE_DIAGNOSTIC_PATH}")

    if args.section in ["all", "logic"]:
        boolean_outputs = run_boolean_workflow()
        graph_layout = _graph_layout.graphviz_input()

        for path in BOOLEAN_EXPECTED_DIR.glob("*.npz"):
            path.unlink()

        save_boolean_expected(boolean_outputs)
        _graph_layout.save_expected(graph_layout)

        for path in sorted(BOOLEAN_EXPECTED_DIR.glob("*.npz")):
            print(f"wrote {path}")

        print(f"wrote {_graph_layout.EXPECTED_PATH}")


if __name__ == "__main__":
    main()
