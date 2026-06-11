#!/usr/bin/env python

import warnings

import networkx as nx
import pandas as pd
import pytest

import bonesistools as bt
from bonesistools.databases.ncbi import _genesyn

GENE_IDS = [
    "20375",
    "78688",
    "17869",
    "435318",
    "18033",
    "12394",
    "12393",
    "12305",
    "22059",
]

GENE_NAMES = [
    "PU.1",
    "B430311C09Rik",
    "Myc",
    "Gm5657",
    "NF-kappaB",
    "Aml1",
    "Aml3",
    "Ddr1",
    "Tp53",
]

NCBI_NAMES = [
    "Spi1",
    "Nol3",
    "Myc",
    "Serpina3d-ps",
    "Nfkb1",
    "Runx1",
    "Runx2",
    "Ddr1",
    "Trp53",
]


class CopyableBooleanNetworkLike(dict):
    """
    Minimal BooleanNetworkLike object with copy support but no rename method.
    """

    def copy(self):
        return type(self)(self)


@pytest.fixture(scope="module")
def mouse_genesyn():
    return bt.dbs.ncbi.GeneSynonyms(organism="mouse")


def test_mouse_gene_synonyms_convert_gene_ids_to_ncbi_names(mouse_genesyn):

    converted = mouse_genesyn.convert_sequence(
        GENE_IDS,
        input_identifier_type="gene_id",
        output_identifier_type="ncbi_name",
    )

    assert converted == NCBI_NAMES


def test_mouse_gene_synonyms_convert_gene_names_to_ncbi_names(mouse_genesyn):

    converted = mouse_genesyn.convert_sequence(
        GENE_NAMES,
        input_identifier_type="name",
        output_identifier_type="ncbi_name",
    )

    assert converted == NCBI_NAMES


def test_gene_synonyms_contains_and_find(mouse_genesyn):
    assert mouse_genesyn.contains() == []
    assert mouse_genesyn.contains("Tp53") == [True]
    assert mouse_genesyn.contains("Tp53", "not-a-gene") == [True, False]
    assert mouse_genesyn.contains("22059", identifier_type="gene_id") == [True]
    assert mouse_genesyn.contains(
        "ENSMUSG00000022346", identifier_type="ensembl_id"
    ) == [True]
    assert mouse_genesyn.contains("MGI:97250", identifier_type="MGI") == [True]
    assert mouse_genesyn.contains("bad-id", identifier_type="MGI") == [False]
    assert mouse_genesyn.find() == []
    assert mouse_genesyn.find("Tp53", "not-a-gene", "Myc") == ["Tp53", "Myc"]
    assert mouse_genesyn.find("22059", "bad-id", identifier_type="gene_id") == ["22059"]
    assert mouse_genesyn.find("bad-id", identifier_type="MGI") == []

    with pytest.raises(
        ValueError, match="invalid argument value for 'identifier_type'"
    ):
        mouse_genesyn.contains("Tp53", identifier_type="bad")

    with pytest.raises(
        ValueError, match="invalid argument value for 'identifier_type'"
    ):
        mouse_genesyn.find("Tp53", identifier_type="bad")


def test_gene_synonyms_identifier_lookups_and_database_aliases(mouse_genesyn):
    assert mouse_genesyn.get_gene_id("NF-kappaB") == "18033"
    assert mouse_genesyn.get_gene_id("ENSMUSG00000022346", "ensembl_id") == "17869"
    assert mouse_genesyn.get_gene_id("MGI:97250", "MGI") == "17869"
    assert mouse_genesyn.get_official_name("Tp53") == "Trp53"
    assert mouse_genesyn.get_ncbi_name("22059", "gene_id") == "Trp53"
    assert mouse_genesyn.get_ensembl_id("Myc") == "ENSMUSG00000022346"
    assert mouse_genesyn.get_alias_from_database("Myc", database="MGI") == "MGI:97250"
    assert (
        mouse_genesyn.conversion(
            "MGI:97250",
            input_identifier_type="MGI",
            output_identifier_type="ensembl_id",
        )
        == "ENSMUSG00000022346"
    )
    assert (
        mouse_genesyn.conversion(
            "Myc",
            output_identifier_type="MGI",
        )
        == "MGI:97250"
    )
    assert mouse_genesyn.convert_sequence(
        ["Myc"],
        output_identifier_type="MGI",
    ) == ["MGI:97250"]

    assert mouse_genesyn.get_gene_id("not-a-gene") is None
    assert mouse_genesyn.convert_sequence(
        ["not-a-gene"],
        output_identifier_type="official_name",
        keep_if_missing=False,
    ) == [None]


def test_gene_synonyms_get_mapping_returns_deepcopy():
    genesyn = object.__new__(bt.dbs.ncbi.GeneSynonyms)
    genesyn.gene_aliases_mapping = {
        "name": {
            "Myc": "17869",
        },
    }

    mapping = genesyn.get_mapping()
    mapping["name"].clear()

    assert genesyn.gene_aliases_mapping == {"name": {"Myc": "17869"}}


def test_gene_synonyms_convert_dataframe_graph_and_interactions(mouse_genesyn):
    df = pd.DataFrame(
        [[1, 2], [3, 4]],
        index=["Tp53", "Myc"],
        columns=["NF-kappaB", "unknown"],
    )
    converted_index = mouse_genesyn.convert_df(df, axis="index", copy=True)
    converted_columns = mouse_genesyn.convert_df(df, axis="columns", copy=True)

    assert converted_index.index.tolist() == ["Trp53", "Myc"]
    assert converted_columns.columns.tolist() == ["Nfkb1", "unknown"]
    assert df.index.tolist() == ["Tp53", "Myc"]
    assert mouse_genesyn(df, axis="index").index.tolist() == ["Trp53", "Myc"]

    in_place_df = df.copy()
    assert mouse_genesyn.convert_df(in_place_df, axis=1, copy=False) is None
    assert in_place_df.columns.tolist() == ["Nfkb1", "unknown"]

    graph = nx.MultiDiGraph()
    graph.add_edge("Tp53", "NF-kappaB", sign=-1)

    converted_graph = mouse_genesyn.convert_graph(graph, copy=True)
    assert set(converted_graph.nodes) == {"Trp53", "Nfkb1"}
    assert set(graph.nodes) == {"Tp53", "NF-kappaB"}
    assert set(mouse_genesyn(graph, copy=True).nodes) == {"Trp53", "Nfkb1"}

    assert mouse_genesyn.convert_graph(graph, copy=False) is None
    assert set(graph.nodes) == {"Trp53", "Nfkb1"}

    interactions = [("Tp53", "NF-kappaB", {"sign": -1})]
    assert mouse_genesyn.convert_interaction_list(interactions) == [
        ("Trp53", "Nfkb1", {"sign": -1})
    ]

    assert mouse_genesyn([]) == []
    assert mouse_genesyn(interactions) == [("Trp53", "Nfkb1", {"sign": -1})]
    assert list(mouse_genesyn(("Tp53", "unknown"))) == ["Trp53", "unknown"]


def test_gene_synonyms_standardize_wrappers_and_legacy_arguments(mouse_genesyn):
    assert mouse_genesyn.standardize_sequence(["Tp53", "NF-kappaB"]) == [
        "Trp53",
        "Nfkb1",
    ]

    bn = bt.bpy.bn.BooleanNetwork({"Tp53": "Myc", "Myc": "Tp53"})
    converted_bn = mouse_genesyn.convert_bn(bn, copy=True)
    assert isinstance(converted_bn, bt.bpy.bn.BooleanNetwork)
    assert converted_bn.rules == {"Trp53": "Myc", "Myc": "Trp53"}
    assert bn.rules == {"Tp53": "Myc", "Myc": "Tp53"}

    dispatched_bn = mouse_genesyn(bn, copy=True)
    assert isinstance(dispatched_bn, bt.bpy.bn.BooleanNetwork)
    assert dispatched_bn.rules == {"Trp53": "Myc", "Myc": "Trp53"}

    standardized_bn = mouse_genesyn.standardize_bn(bn, copy=True)
    assert isinstance(standardized_bn, bt.bpy.bn.BooleanNetwork)
    assert standardized_bn.rules == {"Trp53": "Myc", "Myc": "Trp53"}

    mapping_bn = {"Tp53": "Myc", "Myc": "Tp53"}
    converted_mapping_bn = mouse_genesyn.standardize_bn(mapping_bn, copy=True)
    assert type(converted_mapping_bn) is type(mapping_bn)
    assert converted_mapping_bn == {"Trp53": "Myc", "Myc": "Trp53"}
    assert mapping_bn == {"Tp53": "Myc", "Myc": "Tp53"}

    copyable_bn = CopyableBooleanNetworkLike({"Tp53": "Myc", "Myc": "Tp53"})
    converted_copyable_bn = mouse_genesyn.standardize_bn(copyable_bn, copy=True)
    assert type(converted_copyable_bn) is type(copyable_bn)
    assert converted_copyable_bn == {"Trp53": "Myc", "Myc": "Trp53"}
    assert isinstance(copyable_bn.copy(), CopyableBooleanNetworkLike)

    with pytest.raises(TypeError, match="requires a 'rename' method"):
        mouse_genesyn.standardize_bn(mapping_bn, copy=False)

    with pytest.raises(TypeError, match="Boolean network-like"):
        mouse_genesyn.standardize_bn({"Tp53": object()}, copy=True)

    with pytest.raises(AttributeError, match="no attribute 'get_bad'"):
        mouse_genesyn.convert_bn(bn, output_identifier_type="bad", copy=True)

    df = pd.DataFrame([[1]], index=["Tp53"], columns=["NF-kappaB"])
    assert mouse_genesyn.standardize_df(df).index.tolist() == ["Trp53"]

    graph = nx.Graph()
    graph.add_edge("Tp53", "NF-kappaB")
    assert set(mouse_genesyn.standardize_graph(graph).nodes) == {"Trp53", "Nfkb1"}

    interactions = [("Tp53", "NF-kappaB", {"sign": 1})]
    assert mouse_genesyn.standardize_interaction_list(interactions) == [
        ("Trp53", "Nfkb1", {"sign": 1})
    ]

    with pytest.warns(FutureWarning, match="'gene_type' is deprecated"):
        assert mouse_genesyn.get_gene_id("Tp53", gene_type="name") == "22059"

    with pytest.warns(FutureWarning, match="'alias_gene' is deprecated"):
        assert mouse_genesyn.convert_sequence(["Tp53"], alias_gene="gene_id") == [
            "22059"
        ]

    with pytest.warns(FutureWarning, match="'gene_type' is deprecated"):
        with pytest.raises(TypeError, match="use either 'gene_type'"):
            mouse_genesyn.get_gene_id(
                "Tp53",
                gene_type="name",
                input_identifier_type="name",
            )


def test_gene_synonyms_reset_updates_configuration_without_download(monkeypatch):
    calls = []

    def fake_download(genesyn):
        calls.append(("download", genesyn.version, genesyn.ncbi_file.name))

    def fake_initialize(genesyn, show_warnings):
        calls.append(("initialize", show_warnings))
        genesyn.show_warnings = show_warnings
        genesyn.gene_aliases_mapping = {"databases": {}}
        genesyn.databases = set()

    monkeypatch.setattr(
        bt.dbs.ncbi.GeneSynonyms,
        "_GeneSynonyms__download_gene_info",
        fake_download,
    )
    monkeypatch.setattr(
        bt.dbs.ncbi.GeneSynonyms,
        "_GeneSynonyms__initialize_mappings",
        fake_initialize,
    )

    genesyn = object.__new__(bt.dbs.ncbi.GeneSynonyms)
    genesyn.reset(
        organism="mouse",
        show_warnings=False,
    )

    assert genesyn.organism == "mouse"
    assert genesyn.version == "current"
    assert genesyn.show_warnings is False
    assert genesyn.ncbi_file.name == "mus_musculus_gene_info.tsv"
    assert calls == [
        ("download", "current", "mus_musculus_gene_info.tsv"),
        ("initialize", False),
    ]


def test_gene_synonyms_supports_current_latest_and_local_dated_versions(
    monkeypatch,
    tmp_path,
):
    calls = []
    data_dir = tmp_path / "gi"
    version_dir = data_dir / "versions" / "20240101"
    version_dir.mkdir(parents=True)

    bundled_file = data_dir / "mus_musculus_gene_info.tsv"
    version_file = version_dir / "mus_musculus_gene_info.tsv"
    bundled_file.write_text("bundled\n")
    version_file.write_text("version\n")

    monkeypatch.setattr(_genesyn, "NCBI_DIR", data_dir)
    monkeypatch.setattr(_genesyn, "NCBI_GENE_INFO_VERSION_DIR", data_dir / "versions")
    monkeypatch.setattr(
        _genesyn,
        "NCBI_GENE_INFO_FILES",
        {"mouse": bundled_file},
    )

    def fake_download(genesyn):
        calls.append(("download", genesyn.version, genesyn.ncbi_file))

    def fake_initialize(genesyn, show_warnings):
        calls.append(("initialize", show_warnings))
        genesyn.show_warnings = show_warnings
        genesyn.gene_aliases_mapping = {"databases": {}}
        genesyn.databases = set()

    monkeypatch.setattr(
        bt.dbs.ncbi.GeneSynonyms,
        "_GeneSynonyms__download_gene_info",
        fake_download,
    )
    monkeypatch.setattr(
        bt.dbs.ncbi.GeneSynonyms,
        "_GeneSynonyms__initialize_mappings",
        fake_initialize,
    )

    current = bt.dbs.ncbi.GeneSynonyms(organism="mouse", version="current")
    latest = bt.dbs.ncbi.GeneSynonyms(organism="mouse", version="latest")
    dated = bt.dbs.ncbi.GeneSynonyms(organism="mouse", version="2024-01-01")

    assert current.version == "current"
    assert current.ncbi_file == bundled_file
    assert latest.version == "latest"
    assert latest.ncbi_file == bundled_file
    assert dated.version == "20240101"
    assert dated.ncbi_file == version_file
    assert calls == [
        ("download", "current", bundled_file),
        ("initialize", False),
        ("download", "latest", bundled_file),
        ("initialize", False),
        ("download", "20240101", version_file),
        ("initialize", False),
    ]


def test_gene_synonyms_validation_errors_and_missing_warnings(
    mouse_genesyn, monkeypatch
):
    with pytest.raises(ValueError, match="invalid argument value for 'organism'"):
        bt.dbs.ncbi.GeneSynonyms(organism="rat")

    with pytest.raises(
        TypeError, match="unsupported argument type for 'show_warnings'"
    ):
        bt.dbs.ncbi.GeneSynonyms(show_warnings="yes")

    with pytest.raises(TypeError, match="unsupported argument type for 'version'"):
        bt.dbs.ncbi.GeneSynonyms(version=object())

    genesyn = mouse_genesyn
    monkeypatch.setattr(genesyn, "show_warnings", True)

    with pytest.raises(ValueError, match="invalid argument value for 'organism'"):
        genesyn.reset(organism="rat")

    with pytest.raises(
        TypeError, match="unsupported argument type for 'show_warnings'"
    ):
        genesyn.reset(show_warnings="yes")

    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")

        assert genesyn.get_gene_id("not-a-gene") is None
        assert genesyn.get_gene_id("not-an-ensembl", "ensembl_id") is None
        assert genesyn.get_gene_id("bad-id", "MGI") is None
        assert genesyn.get_ncbi_name("not-a-gene") is None
        assert genesyn.get_official_name("not-a-gene") is None
        assert genesyn.get_ensembl_id("not-a-gene") is None
        assert genesyn.get_alias_from_database("not-a-gene", database="MGI") is None
        assert (
            genesyn.get_alias_from_database(
                "12305",
                database="miRBase",
                input_identifier_type="gene_id",
            )
            is None
        )

    warning_messages = [str(warning.message) for warning in recorded]
    assert any("no correspondence" in message for message in warning_messages)
    assert any("no gene_id correspondence" in message for message in warning_messages)
    assert any(
        "no NCBI reference name correspondence" in message
        for message in warning_messages
    )
    assert any(
        "no official name correspondence" in message for message in warning_messages
    )
    assert any(
        "no Ensembl id correspondence" in message for message in warning_messages
    )
    assert any("no MGI correspondence" in message for message in warning_messages)
    assert any("no miRBase correspondence" in message for message in warning_messages)

    with pytest.raises(
        ValueError, match="invalid argument value for 'input_identifier_type'"
    ):
        genesyn.get_gene_id("Tp53", input_identifier_type="bad")

    with pytest.raises(ValueError, match="invalid argument value for 'database'"):
        genesyn.get_alias_from_database("Tp53", database="bad")

    with pytest.raises(AttributeError, match="no attribute 'get_bad'"):
        genesyn.conversion("Tp53", output_identifier_type="bad")

    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        genesyn.convert_df(pd.DataFrame([[1]]), axis="bad")

    with pytest.raises(TypeError, match="unsupported argument type for 'data'"):
        genesyn(1)
