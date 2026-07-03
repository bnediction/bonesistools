# Golden Regression Tests

Golden tests compare the current implementation against reference outputs
generated from fixed input datasets.

Golden outputs are versioned reference outputs, not permanent mathematical
truths. They may be regenerated after manual validation when an intentional
algorithmic change or a dependency/numerical-environment change modifies the
expected behavior.

This directory is separate from `tests/regression` because golden tests may be
larger, stricter, or only valid in a canonical Python environment. The default
`pytest` command runs `tests/regression` only.

Golden tests are run in CI only with Python 3.13.

Run golden tests explicitly with:

```bash
BONESISTOOLS_RUN_GOLDEN=1 pytest tests/golden
```

## Layout

Keep frozen artifacts in this directory or a dataset-specific subdirectory:

```text
tests/golden/
    pbmc3k.h5ad
    expected/
        qc.npz
        hvg_loess.npz
        pca.npz
        neighbors.npz
```

The files are test resources, not source data loaders. Tests in this directory
should load frozen inputs and compare current outputs to frozen expected
outputs.

## Current Artifacts

### `pbmc3k.h5ad`

Compressed AnnData file generated from `bt.sct.datasets.load("pbmc3k")`.

Properties:

* shape: 2700 observations x 32738 genes;
* expression matrix: sparse CSR count matrix;
* compression: HDF5 gzip compression through `AnnData.write_h5ad`.

Generation command:

```bash
conda run --no-capture-output -n py313 python - <<'PY'
from pathlib import Path
import bonesistools as bt

output = Path("tests/golden/pbmc3k.h5ad")
adata = bt.sct.datasets.load("pbmc3k", quiet=False)
adata.write_h5ad(output, compression="gzip")
PY
```

Source:

* 10x Genomics, 3k PBMCs from a Healthy Donor;
* https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0
* license: Creative Commons Attribution 4.0 International (CC BY 4.0).

### `expected/*.npz`

Compressed NumPy archives storing reference outputs for the deterministic
PBMC3k preprocessing workflow:

* `qc.npz`: observation-level and variable-level QC columns;
* `hvg_loess.npz`: HVG mask, ranks, scores and selected feature names;
* `pca.npz`: PCA embedding, selected loadings and explained variance;
* `neighbors.npz`: sparse distance and connectivity graph arrays.

Regenerate these files intentionally with:

```bash
BONESISTOOLS_RUN_GOLDEN=1 python tests/golden/generate_expected.py
```

## Updating Golden Files

Only update golden files when an intentional algorithmic change or bug fix
changes the expected behavior.

When updating a golden file, document:

* the command or script used to regenerate it;
* the Python environment used;
* the reason the reference output changed.

Do not update golden files to hide a failing test.

See `docs/TESTING.md` for the golden acceptance-test update procedure.
