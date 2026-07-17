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
    omics/
        pbmc3k.h5ad
        expected/
            qc.npz
            hvg_loess.npz
            hvg_binning.npz
            pca.npz
            neighbors.npz
            spectral.npz
            umap.npz
            tsne.npz
    logic/
        data/
            README.md
            death_receptor_cell_fate.zginml
            sea_urchin_dorsal_ventral_axis.sbml
            synthetic_30_node_dynamics.zginml
        expected/
            aggregated_influence_graph.dot
            death_receptor_cell_fate_mp.npz
            sea_urchin_dorsal_ventral_axis.bnet
            synthetic_30_node_dynamics.npz
```

The files are test resources, not source data loaders. Tests in this directory
should load frozen inputs and compare current outputs to frozen expected
outputs.

## Current Artifacts

### `omics/pbmc3k.h5ad`

Compressed AnnData file generated from `bt.omics.io.load("pbmc3k")`.

Properties:

* shape: 2700 observations x 32738 genes;
* expression matrix: sparse CSR count matrix;
* compression: HDF5 gzip compression through `AnnData.write_h5ad`.

Generation command:

```bash
conda run --no-capture-output -n py313 python - <<'PY'
from pathlib import Path
import bonesistools as bt

output = Path("tests/golden/omics/pbmc3k.h5ad")
adata = bt.omics.io.load("pbmc3k", quiet=False)
adata.write_h5ad(output, compression="gzip")
PY
```

Source:

* 10x Genomics, 3k PBMCs from a Healthy Donor;
* https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0
* license: Creative Commons Attribution 4.0 International (CC BY 4.0).

### `omics/expected/*.npz`

Compressed NumPy archives storing reference outputs for the deterministic
PBMC3k preprocessing workflow:

* `qc.npz`: observation-level and variable-level QC columns;
* `hvg_loess.npz`: HVG mask, ranks, scores and selected feature names;
* `hvg_binning.npz`: mean-binned dispersion HVG mask, ranks, scores and
  selected feature names;
* `pca.npz`: PCA embedding, selected loadings and explained variance;
* `neighbors.npz`: sparse distance and fuzzy connectivity graph arrays;
* `spectral.npz`: spectral embedding computed from binary neighbors;
* `umap.npz`: UMAP embedding computed from fuzzy neighbors;
* `tsne.npz`: t-SNE embedding computed from the PCA representation.

Regenerate these files intentionally with:

```bash
BONESISTOOLS_RUN_GOLDEN=1 python tests/golden/generate_expected.py --section omics
```

### `logic/data/death_receptor_cell_fate.zginml`

Compressed GINML model for death-receptor cell fate decision, used as a
realistic Boolean-network fixture for most-permissive reachable-attractor
golden tests.

Properties:

* model: Calzone et al. death receptor engagement model;
* format: ZGINML;
* variables: 28 Boolean components;
* golden coverage: `most-permissive` reachable attractors from selected
  initial states.

Source:

* BioModels MODEL0912180000;
* https://www.ebi.ac.uk/biomodels/MODEL0912180000
* license: Creative Commons CC0/public domain dedication for the encoded
  model.

Citation:

Calzone, L., Tournier, L., Fourquet, S., Thieffry, D., Zhivotovsky, B.,
Barillot, E., & Zinovyev, A. (2010). Mathematical modelling of cell-fate
decision in response to death receptor engagement. PLoS computational biology,
6(3), e1000702.

### `logic/data/sea_urchin_dorsal_ventral_axis.sbml`

SBML Level 3 Qual model of the regulatory network controlling dorsal-ventral
axis specification in the sea urchin embryo. It exercises Booleanization of
components with maximum levels 2 and 3, MathML logical conditions, signed
influences, and layout metadata.

Properties:

* source components: 31;
* Booleanized components: 42;
* maximum source level: 3;
* format: SBML Level 3 Qual;
* golden coverage: rules, signed influence graph, threshold metadata, and
  layout metadata.

Source, licensing, checksums, and citations for the logical-model fixtures are
documented in `logic/data/README.md`.

### `logic/data/synthetic_30_node_dynamics.zginml`

Synthetic ZGINML Boolean model generated for BoNesisTools golden tests. It is
designed to exercise all reachable-attractor update semantics without making
the golden suite expensive.

Properties:

* variables: 30 Boolean components;
* initial state: one fully specified state;
* golden coverage: `synchronous`, `asynchronous`, `general`, and
  `most-permissive` reachable attractors.

### `logic/expected/`

Compressed NumPy archives storing reference outputs for deterministic Boolean
workflows:

* `aggregated_influence_graph.dot`: canonical Graphviz input generated from a
  fixed Boolean network ensemble, freezing node and edge insertion order
  independently of the installed Graphviz rendering engine;
* `death_receptor_cell_fate_mp.npz`: parsed rules, influence-graph shape,
  selected initial states, and most-permissive reachable attractors;
* `synthetic_30_node_dynamics.npz`: parsed rules, influence-graph shape, and
  reachable attractors under the four supported update semantics.

The readable `logic/expected/sea_urchin_dorsal_ventral_axis.bnet` reference
contains the 42 Boolean rules produced independently with BioLQM from the
bundled SBML model. Golden tests compare imported rules semantically with this
reference and separately validate the complete signed influence graph.

The BioLQM reference was generated with:

```bash
bioLQM \
    tests/golden/logic/data/sea_urchin_dorsal_ventral_axis.sbml \
    -m booleanize \
    -of bnet \
    tests/golden/logic/expected/sea_urchin_dorsal_ventral_axis.bnet
```

Regenerate only Boolean golden outputs intentionally with:

```bash
BONESISTOOLS_RUN_GOLDEN=1 python tests/golden/generate_expected.py --section logic
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
