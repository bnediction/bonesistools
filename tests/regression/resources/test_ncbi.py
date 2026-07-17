#!/usr/bin/env python

import gzip
import warnings
from pathlib import Path
from typing import Any, cast

import networkx as nx
import pandas as pd
import pytest

import bonesistools as bt
from bonesistools.resources.ncbi import _genesyn

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

NCBI_SYMBOLS = [
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
    return bt.resources.ncbi.genesyn(organism="mouse")


def test_mouse_gene_synonyms_convert_gene_ids_to_ncbi_symbols(mouse_genesyn):

    converted = mouse_genesyn.convert_sequence(
        GENE_IDS,
        input_type="gene_id",
        output_type="ncbi_symbol",
    )

    assert converted == NCBI_SYMBOLS


def test_mouse_gene_synonyms_convert_gene_names_to_ncbi_symbols(mouse_genesyn):

    converted = mouse_genesyn.convert_sequence(
        GENE_NAMES,
        input_type="name",
        output_type="ncbi_symbol",
    )

    assert converted == NCBI_SYMBOLS


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
    assert mouse_genesyn.get_symbol("Tp53") == "Trp53"
    assert mouse_genesyn.get_ncbi_symbol("22059", "gene_id") == "Trp53"
    assert mouse_genesyn.get_ensembl_id("Myc") == "ENSMUSG00000022346"
    assert mouse_genesyn.get_alias_from_database("Myc", database="MGI") == "MGI:97250"
    assert (
        mouse_genesyn.conversion(
            "22059",
            input_type="gene_id",
            output_type="gene_id",
        )
        == "22059"
    )
    assert (
        mouse_genesyn.conversion(
            "MGI:97250",
            input_type="MGI",
            output_type="ensembl_id",
        )
        == "ENSMUSG00000022346"
    )
    assert (
        mouse_genesyn.conversion(
            "Myc",
            output_type="MGI",
        )
        == "MGI:97250"
    )
    assert mouse_genesyn.convert_sequence(
        ["Myc"],
        output_type="MGI",
    ) == ["MGI:97250"]

    assert mouse_genesyn.get_gene_id("not-a-gene") is None
    assert mouse_genesyn.convert_sequence(
        ["not-a-gene"],
        output_type="symbol",
        keep_if_missing=False,
    ) == [None]


@pytest.mark.parametrize(
    ("deprecated", "replacement"),
    [
        ("get_official_name", "get_symbol"),
        ("get_ncbi_name", "get_ncbi_symbol"),
    ],
)
def test_gene_synonyms_deprecated_getters_warn(
    mouse_genesyn,
    deprecated,
    replacement,
):
    with pytest.warns(FutureWarning, match=deprecated):
        result = getattr(mouse_genesyn, deprecated)("Tp53")

    assert result == getattr(mouse_genesyn, replacement)("Tp53")


@pytest.mark.parametrize(
    ("deprecated", "replacement"),
    [
        ("official_name", "symbol"),
        ("ncbi_name", "ncbi_symbol"),
    ],
)
def test_gene_synonyms_deprecated_output_types_warn(
    mouse_genesyn,
    deprecated,
    replacement,
):
    with pytest.warns(FutureWarning, match=deprecated):
        result = mouse_genesyn.conversion(
            "Tp53",
            output_type=deprecated,
        )

    assert result == mouse_genesyn.conversion(
        "Tp53",
        output_type=replacement,
    )


@pytest.mark.parametrize(
    ("deprecated", "replacement", "gene", "value", "expected"),
    [
        ("input_identifier_type", "input_type", "22059", "gene_id", "Trp53"),
        ("output_identifier_type", "output_type", "Tp53", "gene_id", "22059"),
    ],
)
def test_gene_synonyms_deprecated_identifier_argument_names_warn(
    mouse_genesyn,
    deprecated,
    replacement,
    gene,
    value,
    expected,
):
    with pytest.warns(FutureWarning, match=f"{deprecated}.*{replacement}"):
        result = mouse_genesyn.conversion(gene, **{deprecated: value})

    assert result == expected
    assert mouse_genesyn.conversion(gene, **{replacement: value}) == expected


def test_gene_synonyms_rejects_deprecated_and_current_argument_names(
    mouse_genesyn,
):
    with pytest.warns(FutureWarning, match="input_identifier_type.*input_type"):
        with pytest.raises(TypeError, match="use either 'input_identifier_type'"):
            mouse_genesyn.conversion(
                "Tp53",
                input_type="name",
                **{"input_identifier_type": "name"},
            )


@pytest.mark.parametrize(
    ("input_type", "gene"),
    [
        ("name", "Myc"),
        ("gene_id", "17869"),
        ("ensembl_id", "ENSMUSG00000022346"),
        ("MGI", "MGI:97250"),
    ],
)
@pytest.mark.parametrize(
    ("output_type", "expected"),
    [
        ("symbol", "Myc"),
        ("ncbi_symbol", "Myc"),
        ("gene_id", "17869"),
        ("ensembl_id", "ENSMUSG00000022346"),
        ("MGI", "MGI:97250"),
    ],
)
def test_gene_synonyms_conversion_supports_all_identifier_pairs(
    mouse_genesyn,
    input_type,
    gene,
    output_type,
    expected,
):
    assert (
        mouse_genesyn.conversion(
            gene,
            input_type=input_type,
            output_type=output_type,
        )
        == expected
    )


@pytest.mark.parametrize(
    ("attribute", "value"),
    [
        ("organism", "human"),
        ("version", "latest"),
        ("ncbi_file", Path("other.tsv")),
        ("show_warnings", True),
        ("databases", frozenset()),
        ("valid_input_identifier_types", ()),
        ("valid_output_identifier_types", ()),
        ("gene_aliases_mapping", {}),
    ],
)
def test_gene_synonyms_configuration_is_read_only(
    mouse_genesyn,
    attribute,
    value,
):
    with pytest.raises(AttributeError, match="read-only"):
        setattr(mouse_genesyn, attribute, value)

    assert isinstance(mouse_genesyn.databases, frozenset)
    with pytest.raises(AttributeError):
        getattr(mouse_genesyn, "gene_aliases_mapping")


def test_gene_synonyms_parses_bundled_gene_info_and_resolves_conflicts(tmp_path):
    gene_info = tmp_path / "gene_info.tsv"
    gene_info.write_text(
        "\t".join(
            [
                "gene_id",
                "symbol",
                "ncbi_symbol",
                "synonyms",
                "dbXrefs",
                "chromosome",
                "gene_type",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "1",
                "OFF1",
                "Ref1",
                "Alias1",
                "MGI:MGI1|Ensembl:ENSMUSG00000000001",
                "1",
                "protein-coding",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "2",
                "-",
                "Ref2",
                "OFF1|Alias2",
                "MGI:MGI2",
                "2",
                "pseudo",
            ]
        )
        + "\n"
        + "\t".join(["3", "Low", "LowRef", "Shared|Future", "-", "3", "pseudo"])
        + "\n"
        + "\t".join(
            ["4", "High", "HighRef", "Shared", "MGI:MGI4", "MT", "protein-coding"]
        )
        + "\n"
        + "\t".join(
            [
                "5",
                "Future",
                "FutureRef",
                "-",
                "Ensembl:ENSG00000121410",
                "5",
                "protein-coding",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    genesyn = object.__new__(_genesyn.GeneSynonyms)
    genesyn._organism = "mouse"
    genesyn._version = "bundled"
    genesyn._ncbi_file = gene_info

    genesyn._initialize_mappings(show_warnings=False)

    assert genesyn.get_gene_id("OFF1") == "1"
    assert genesyn.get_gene_id("Ref1") == "1"
    assert genesyn.get_gene_id("Alias2") == "2"
    assert genesyn.get_gene_id("Ref2") == "2"
    assert genesyn.get_gene_id("Shared") == "4"
    assert genesyn.get_gene_id("Future") == "5"
    assert genesyn.get_gene_id("ENSG00000121410", "ensembl_id") == "5"
    assert genesyn.get_symbol("2", input_type="gene_id") == "Ref2"
    assert genesyn.get_ensembl_id("1", input_type="gene_id") == ("ENSMUSG00000000001")
    assert (
        genesyn.get_alias_from_database(
            "1",
            database="MGI",
            input_type="gene_id",
        )
        == "MGI1"
    )
    table = genesyn.to_dataframe()
    assert table.columns.tolist() == [
        "symbol",
        "ncbi_symbol",
        "ensembl_id",
        "MGI",
        "chromosome",
        "gene_type",
    ]
    assert table.index.name == "gene_id"
    assert table.index.tolist() == ["1", "2", "3", "4", "5"]
    assert table.shape == (5, 6)
    first_gene = table.loc["1"]
    assert first_gene.to_dict() == {
        "symbol": "OFF1",
        "ncbi_symbol": "Ref1",
        "ensembl_id": "ENSMUSG00000000001",
        "chromosome": "1",
        "gene_type": "protein-coding",
        "MGI": "MGI1",
    }

    table.loc["1", "symbol"] = "changed"
    assert genesyn.get_symbol("1", input_type="gene_id") == "OFF1"
    assert (
        genesyn.get_alias_from_database(
            "1",
            database="MGI",
            input_type="gene_id",
        )
        == "MGI1"
    )
    assert genesyn.get_gene_id("MGI1", input_type="MGI") == "1"
    assert genesyn.contains("MGI4", identifier_type="MGI") == [True]
    assert genesyn.find("MGI1", "missing", identifier_type="MGI") == ["MGI1"]
    assert genesyn.databases == {"MGI"}
    assert genesyn.valid_input_identifier_types == (
        "name",
        "gene_id",
        "ensembl_id",
        "MGI",
    )
    assert genesyn.valid_output_identifier_types == (
        "symbol",
        "ncbi_symbol",
        "gene_id",
        "ensembl_id",
        "MGI",
    )
    assert not Path(f"{gene_info}_cut").exists()


def test_gene_synonyms_reduces_full_gene_info_without_shell_tools(tmp_path):
    full_gene_info = tmp_path / "full_gene_info.tsv"
    bundled_gene_info = tmp_path / "bundled_gene_info.tsv.gz"
    full_gene_info.write_text(
        "\t".join(
            [
                "#tax_id",
                "GeneID",
                "Symbol",
                "LocusTag",
                "Synonyms",
                "dbXrefs",
                "chromosome",
                "map_location",
                "description",
                "type_of_gene",
                "Symbol_from_nomenclature_authority",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "10090",
                "10",
                "RefA",
                "-",
                "AliasA|AliasB",
                "Ensembl:ENSMUSG00000000010|MGI:MGI10",
                "1",
                "-",
                "-",
                "protein-coding",
                "OfficialA",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    genesyn = object.__new__(_genesyn.GeneSynonyms)
    getattr(genesyn, "_GeneSynonyms__write_bundled_gene_info_file")(
        full_gene_info,
        bundled_gene_info,
    )

    with gzip.open(bundled_gene_info, "rt", encoding="utf-8") as file:
        rows = [line.rstrip("\n").split("\t") for line in file]

    assert rows == [
        [
            "gene_id",
            "symbol",
            "ncbi_symbol",
            "synonyms",
            "dbXrefs",
            "chromosome",
            "gene_type",
        ],
        [
            "10",
            "OfficialA",
            "RefA",
            "AliasA|AliasB",
            "Ensembl:ENSMUSG00000000010|MGI:MGI10",
            "1",
            "protein-coding",
        ],
    ]
    assert not Path(f"{bundled_gene_info}.tmp").exists()


def test_gene_synonyms_parses_latest_gene_info_and_removes_temporary_files(tmp_path):
    gene_info = tmp_path / "latest.tsv"
    gene_info.write_text(
        "\t".join(
            [
                "#tax_id",
                "GeneID",
                "Symbol",
                "LocusTag",
                "Synonyms",
                "dbXrefs",
                "chromosome",
                "map_location",
                "description",
                "type_of_gene",
                "Symbol_from_nomenclature_authority",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "10090",
                "10",
                "RefA",
                "-",
                "AliasA",
                "Ensembl:ENSMUSG00000000010|MGI:MGI10",
                "1",
                "-",
                "-",
                "protein-coding",
                "OfficialA",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    gzip_file = Path(f"{gene_info}.gz")
    gzip_file.write_text("partial", encoding="utf-8")

    genesyn = object.__new__(_genesyn.GeneSynonyms)
    genesyn._organism = "mouse"
    genesyn._version = "latest"
    genesyn._ncbi_file = gene_info

    genesyn._initialize_mappings(show_warnings=False)

    assert genesyn.get_gene_id("AliasA") == "10"
    assert genesyn.get_symbol("10", input_type="gene_id") == ("OfficialA")
    assert genesyn.get_ensembl_id("10", input_type="gene_id") == ("ENSMUSG00000000010")
    assert (
        genesyn.get_alias_from_database(
            "10",
            database="MGI",
            input_type="gene_id",
        )
        == "MGI10"
    )
    assert not gene_info.exists()
    assert not gzip_file.exists()


def test_gene_synonyms_convert_dataframe_graph_and_interactions(mouse_genesyn):
    df = pd.DataFrame(
        [[1, 2], [3, 4]],
        index=["Tp53", "Myc"],
        columns=["NF-kappaB", "unknown"],
    )
    converted_index = mouse_genesyn.convert_dataframe(df, axis="index", copy=True)
    converted_columns = mouse_genesyn.convert_dataframe(
        df,
        axis="columns",
        copy=True,
    )

    assert converted_index.index.tolist() == ["Trp53", "Myc"]
    assert converted_columns.columns.tolist() == ["Nfkb1", "unknown"]
    assert df.index.tolist() == ["Tp53", "Myc"]
    assert mouse_genesyn(df, axis="index").index.tolist() == ["Trp53", "Myc"]

    in_place_df = df.copy()
    assert mouse_genesyn.convert_dataframe(in_place_df, axis=1, copy=False) is None
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

    hypercube = bt.logic.ba.Hypercube({"Tp53": 1, "NF-kappaB": 0, "unknown": "*"})
    converted_hypercube = mouse_genesyn.convert_hypercube(hypercube, copy=True)
    assert isinstance(converted_hypercube, bt.logic.ba.Hypercube)
    assert converted_hypercube == {"Trp53": 1, "Nfkb1": 0, "unknown": "*"}
    assert hypercube == {"Tp53": 1, "NF-kappaB": 0, "unknown": "*"}

    dispatched_hypercube = mouse_genesyn(hypercube)
    assert isinstance(dispatched_hypercube, bt.logic.ba.Hypercube)
    assert dispatched_hypercube == {"Trp53": 1, "Nfkb1": 0, "unknown": "*"}

    assert mouse_genesyn.convert_hypercube(hypercube, copy=False) is None
    assert hypercube == {"Trp53": 1, "Nfkb1": 0, "unknown": "*"}


def test_gene_synonyms_convert_influence_graph_preserves_opposite_signs(
    mouse_genesyn,
    monkeypatch,
):
    monkeypatch.setattr(
        _genesyn.GeneSynonyms,
        "_GeneSynonyms__conversion_function",
        lambda *args, **kwargs: lambda gene, **_kwargs: {"c": "b"}.get(gene, gene),
    )

    graph = bt.logic.ig.InfluenceGraph()
    graph.add_edge("a", "b", sign=1)
    graph.add_edge("a", "c", sign=-1)

    assert mouse_genesyn.convert_graph(graph, copy=False) is None

    assert sorted((s, t, d["sign"]) for s, t, d in graph.edges(data=True)) == [
        ("a", "b", -1),
        ("a", "b", 1),
    ]


def test_gene_synonyms_convert_influence_graph_ignores_duplicate_edges_after_relabel(
    mouse_genesyn,
    monkeypatch,
):
    monkeypatch.setattr(
        _genesyn.GeneSynonyms,
        "_GeneSynonyms__conversion_function",
        lambda *args, **kwargs: lambda gene, **_kwargs: {"c": "b"}.get(gene, gene),
    )

    graph = bt.logic.ig.InfluenceGraph()
    graph.add_edge("a", "b", sign=1)
    graph.add_edge("a", "c", sign=1)
    graph.add_edge("c", "d", sign=1)

    assert mouse_genesyn.convert_graph(graph, copy=False) is None

    assert sorted((s, t, d["sign"]) for s, t, d in graph.edges(data=True)) == [
        ("a", "b", 1),
        ("b", "d", 1),
    ]


def test_gene_synonyms_convert_aggregated_influence_graph_preserves_total_and_counts(
    mouse_genesyn,
):
    graph = bt.logic.ig.AggregatedInfluenceGraph(total=4)
    graph.add_edge("Tp53", "Myc", sign=1, count=3)

    converted = mouse_genesyn.convert_graph(graph, copy=True)

    assert isinstance(converted, bt.logic.ig.AggregatedInfluenceGraph)
    assert converted.total == 4
    assert converted.edge_count("Trp53", "Myc") == 3
    assert graph.total == 4
    assert graph.edge_count("Tp53", "Myc") == 3


def test_gene_synonyms_convert_hypercube_rejects_component_merges(
    mouse_genesyn,
    monkeypatch,
):
    monkeypatch.setattr(
        _genesyn.GeneSynonyms,
        "_GeneSynonyms__conversion_function",
        lambda *args, **kwargs: lambda gene, **_kwargs: {"a": "x", "b": "x"}[gene],
    )

    hypercube = bt.logic.ba.Hypercube({"a": 0, "b": 1})

    with pytest.raises(ValueError, match="merge hypercube components"):
        mouse_genesyn.convert_hypercube(hypercube)


def test_gene_synonyms_short_conversion_aliases_are_deprecated(mouse_genesyn):
    dataframe = pd.DataFrame([[1]], index=["Tp53"])
    with pytest.warns(FutureWarning, match="convert_df.*convert_dataframe"):
        converted_dataframe = mouse_genesyn.convert_df(dataframe)
    assert converted_dataframe.index.tolist() == ["Trp53"]

    network = bt.logic.bn.BooleanNetwork({"Tp53": "Myc", "Myc": "Tp53"})
    with pytest.warns(FutureWarning, match="convert_bn.*convert_boolean_network"):
        converted_network = mouse_genesyn.convert_bn(network, copy=True)
    assert isinstance(converted_network, bt.logic.bn.BooleanNetwork)
    assert converted_network.rules == {"Trp53": "Myc", "Myc": "Trp53"}


def test_gene_synonyms_standardize_wrappers_and_legacy_arguments(mouse_genesyn):
    def deprecated(old_name, new_name):
        return pytest.warns(FutureWarning, match=f"{old_name}.*{new_name}")

    with deprecated("standardize_sequence", "convert_sequence"):
        standardized_sequence = mouse_genesyn.standardize_sequence(
            ["Tp53", "NF-kappaB"]
        )
    assert standardized_sequence == [
        "Trp53",
        "Nfkb1",
    ]

    bn = bt.logic.bn.BooleanNetwork({"Tp53": "Myc", "Myc": "Tp53"})
    converted_bn = mouse_genesyn.convert_boolean_network(bn, copy=True)
    assert isinstance(converted_bn, bt.logic.bn.BooleanNetwork)
    assert converted_bn.rules == {"Trp53": "Myc", "Myc": "Trp53"}
    assert bn.rules == {"Tp53": "Myc", "Myc": "Tp53"}

    dispatched_bn = mouse_genesyn(bn, copy=True)
    assert isinstance(dispatched_bn, bt.logic.bn.BooleanNetwork)
    assert dispatched_bn.rules == {"Trp53": "Myc", "Myc": "Trp53"}

    with deprecated("standardize_bn", "convert_boolean_network"):
        standardized_bn = mouse_genesyn.standardize_bn(bn, copy=True)
    assert isinstance(standardized_bn, bt.logic.bn.BooleanNetwork)
    assert standardized_bn.rules == {"Trp53": "Myc", "Myc": "Trp53"}

    mapping_bn = {"Tp53": "Myc", "Myc": "Tp53"}
    with deprecated("standardize_bn", "convert_boolean_network"):
        converted_mapping_bn = mouse_genesyn.standardize_bn(mapping_bn, copy=True)
    assert type(converted_mapping_bn) is type(mapping_bn)
    assert converted_mapping_bn == {"Trp53": "Myc", "Myc": "Trp53"}
    assert mapping_bn == {"Tp53": "Myc", "Myc": "Tp53"}

    copyable_bn = CopyableBooleanNetworkLike({"Tp53": "Myc", "Myc": "Tp53"})
    with deprecated("standardize_bn", "convert_boolean_network"):
        converted_copyable_bn = mouse_genesyn.standardize_bn(copyable_bn, copy=True)
    assert type(converted_copyable_bn) is type(copyable_bn)
    assert converted_copyable_bn == {"Trp53": "Myc", "Myc": "Trp53"}
    assert isinstance(copyable_bn.copy(), CopyableBooleanNetworkLike)

    with deprecated("standardize_bn", "convert_boolean_network"):
        with pytest.raises(TypeError, match="requires a 'rename' method"):
            mouse_genesyn.standardize_bn(mapping_bn, copy=False)

    with deprecated("standardize_bn", "convert_boolean_network"):
        with pytest.raises(TypeError, match="Boolean network-like"):
            mouse_genesyn.standardize_bn({"Tp53": object()}, copy=True)

    with pytest.raises(AttributeError, match="no attribute 'get_bad'"):
        mouse_genesyn.convert_boolean_network(bn, output_type="bad", copy=True)

    df = pd.DataFrame([[1]], index=["Tp53"], columns=["NF-kappaB"])
    with deprecated("standardize_df", "convert_dataframe"):
        standardized_df = mouse_genesyn.standardize_df(df)
    assert standardized_df.index.tolist() == ["Trp53"]

    graph = nx.Graph()
    graph.add_edge("Tp53", "NF-kappaB")
    with deprecated("standardize_graph", "convert_graph"):
        standardized_graph = mouse_genesyn.standardize_graph(graph)
    assert set(standardized_graph.nodes) == {"Trp53", "Nfkb1"}

    interactions = [("Tp53", "NF-kappaB", {"sign": 1})]
    with deprecated("standardize_interaction_list", "convert_interaction_list"):
        standardized_interactions = mouse_genesyn.standardize_interaction_list(
            interactions
        )
    assert standardized_interactions == [("Trp53", "Nfkb1", {"sign": 1})]

    hypercube = bt.logic.ba.Hypercube({"Tp53": 1, "NF-kappaB": 0})
    with deprecated("standardize_hypercube", "convert_hypercube"):
        standardized_hypercube = mouse_genesyn.standardize_hypercube(hypercube)
    assert standardized_hypercube == {
        "Trp53": 1,
        "Nfkb1": 0,
    }
    assert hypercube == {"Tp53": 1, "NF-kappaB": 0}

    with pytest.warns(FutureWarning):
        assert mouse_genesyn.get_gene_id("Tp53", gene_type="name") == "22059"

    with pytest.warns(FutureWarning):
        assert mouse_genesyn.convert_sequence(["Tp53"], alias_gene="gene_id") == [
            "22059"
        ]

    with pytest.warns(FutureWarning):
        with pytest.raises(TypeError):
            mouse_genesyn.get_gene_id(
                "Tp53",
                gene_type="name",
                input_type="name",
            )


def test_gene_synonyms_reset_updates_configuration_without_download(monkeypatch):
    calls = []

    def fake_download(genesyn):
        calls.append(("download", genesyn.version, genesyn.ncbi_file.name))

    def fake_initialize(genesyn, show_warnings):
        calls.append(("initialize", show_warnings))
        genesyn._show_warnings = show_warnings

    monkeypatch.setattr(
        _genesyn.GeneSynonyms,
        "_GeneSynonyms__download_gene_info",
        fake_download,
    )
    monkeypatch.setattr(
        _genesyn.GeneSynonyms,
        "_initialize_mappings",
        fake_initialize,
    )

    genesyn = object.__new__(_genesyn.GeneSynonyms)
    genesyn.reset(
        organism="mouse",
        show_warnings=False,
    )

    assert genesyn.organism == "mouse"
    assert genesyn.version == "bundled"
    assert genesyn.show_warnings is False
    assert genesyn.ncbi_file.name == "mus_musculus_gene_info.tsv.gz"
    assert calls == [
        ("download", "bundled", "mus_musculus_gene_info.tsv.gz"),
        ("initialize", False),
    ]


def test_gene_synonyms_reset_restores_state_after_failure(monkeypatch):
    genesyn = object.__new__(_genesyn.GeneSynonyms)
    original_mapping = {"name": {"MYC": "17869"}}
    genesyn._organism = "mouse"
    genesyn._version = "bundled"
    genesyn._ncbi_file = Path("mouse.tsv")
    genesyn._show_warnings = False
    genesyn._gene_aliases_mapping = original_mapping
    genesyn._databases = frozenset({"MGI"})
    genesyn._valid_input_identifier_types = ("name", "MGI")
    genesyn._valid_output_identifier_types = ("symbol", "MGI")

    monkeypatch.setattr(
        _genesyn.GeneSynonyms,
        "_GeneSynonyms__download_gene_info",
        lambda self: None,
    )

    def fail_after_mutation(self, show_warnings):
        self._show_warnings = show_warnings
        self._gene_aliases_mapping = {"broken": {}}
        raise RuntimeError("invalid gene_info")

    monkeypatch.setattr(
        _genesyn.GeneSynonyms,
        "_initialize_mappings",
        fail_after_mutation,
    )

    with pytest.raises(RuntimeError, match="invalid gene_info"):
        genesyn.reset(organism="human", show_warnings=True)

    assert genesyn.organism == "mouse"
    assert genesyn.version == "bundled"
    assert genesyn.ncbi_file == Path("mouse.tsv")
    assert genesyn.show_warnings is False
    assert genesyn._gene_aliases_mapping is original_mapping
    assert genesyn.databases == {"MGI"}
    assert genesyn.valid_input_identifier_types == ("name", "MGI")
    assert genesyn.valid_output_identifier_types == ("symbol", "MGI")


def test_gene_synonyms_supports_bundled_latest_and_local_dated_versions(
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
        genesyn._show_warnings = show_warnings

    monkeypatch.setattr(
        _genesyn.GeneSynonyms,
        "_GeneSynonyms__download_gene_info",
        fake_download,
    )
    monkeypatch.setattr(
        _genesyn.GeneSynonyms,
        "_initialize_mappings",
        fake_initialize,
    )

    bundled = bt.resources.ncbi.genesyn(organism="mouse", version="bundled")
    latest = bt.resources.ncbi.genesyn(organism="mouse", version="latest")
    dated = bt.resources.ncbi.genesyn(organism="mouse", version="2024-01-01")

    assert bundled.version == "bundled"
    assert bundled.ncbi_file == bundled_file
    assert latest.version == "latest"
    assert latest.ncbi_file != bundled_file
    assert dated.version == "20240101"
    assert dated.ncbi_file == version_file
    assert calls[0] == ("download", "bundled", bundled_file)
    assert calls[1] == ("initialize", False)
    assert calls[2][0:2] == ("download", "latest")
    assert calls[2][2] != bundled_file
    assert calls[3:] == [
        ("initialize", False),
        ("download", "20240101", version_file),
        ("initialize", False),
    ]


def test_gene_synonyms_validation_errors_and_missing_warnings(
    mouse_genesyn, monkeypatch
):
    with pytest.raises(ValueError, match="invalid argument value for 'organism'"):
        bt.resources.ncbi.genesyn(organism="rat")

    with pytest.raises(
        TypeError, match="unsupported argument type for 'show_warnings'"
    ):
        bt.resources.ncbi.genesyn(show_warnings=cast(Any, "yes"))

    with pytest.raises(TypeError, match="unsupported argument type for 'version'"):
        bt.resources.ncbi.genesyn(version=cast(Any, object()))

    genesyn = mouse_genesyn
    monkeypatch.setattr(genesyn, "_show_warnings", True)

    with pytest.raises(ValueError, match="invalid argument value for 'organism'"):
        genesyn.reset(organism="rat")

    with pytest.raises(
        TypeError, match="unsupported argument type for 'show_warnings'"
    ):
        genesyn.reset(show_warnings=cast(Any, "yes"))

    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")

        assert genesyn.get_gene_id("not-a-gene") is None
        assert genesyn.get_gene_id("not-an-ensembl", "ensembl_id") is None
        assert genesyn.get_gene_id("bad-id", "MGI") is None
        assert genesyn.get_ncbi_symbol("not-a-gene") is None
        assert genesyn.get_symbol("not-a-gene") is None
        assert genesyn.get_ensembl_id("not-a-gene") is None
        assert genesyn.get_alias_from_database("not-a-gene", database="MGI") is None
        assert (
            genesyn.get_alias_from_database(
                "12305",
                database="miRBase",
                input_type="gene_id",
            )
            is None
        )

    warning_messages = [str(warning.message) for warning in recorded]
    assert any("no correspondence" in message for message in warning_messages)
    assert any("no gene_id correspondence" in message for message in warning_messages)
    assert any(
        "no NCBI symbol correspondence" in message for message in warning_messages
    )
    assert any("no symbol correspondence" in message for message in warning_messages)
    assert any(
        "no Ensembl id correspondence" in message for message in warning_messages
    )
    assert any("no MGI correspondence" in message for message in warning_messages)
    assert any("no miRBase correspondence" in message for message in warning_messages)

    with pytest.raises(ValueError, match="invalid argument value for 'input_type'"):
        genesyn.get_gene_id("Tp53", input_type="bad")

    with pytest.raises(ValueError, match="invalid argument value for 'database'"):
        genesyn.get_alias_from_database("Tp53", database="bad")

    with pytest.raises(AttributeError, match="no attribute 'get_bad'"):
        genesyn.conversion("Tp53", output_type="bad")

    with pytest.raises(ValueError, match="invalid argument value for 'axis'"):
        genesyn.convert_dataframe(pd.DataFrame([[1]]), axis="bad")

    with pytest.raises(TypeError, match="unsupported argument type for 'data'"):
        genesyn(cast(Any, 1))
