#!/usr/bin/env python

from typing import Any, cast

import pandas as pd
import pytest

from bonesistools.logic.influence_graph import InfluenceGraph
from bonesistools.resources.omnipath import _archive, _collectri, _dorothea


def test_resolve_interactions_archive_supports_latest_and_dates(monkeypatch):
    monkeypatch.setattr(
        _archive,
        "_list_interactions_archives",
        lambda: [
            ("20230101", "20230501", "omnipath_webservice_interactions__a.tsv.xz"),
            (
                "20230502",
                "20240101",
                "omnipath_webservice_interactions__20230502-20240101.tsv.xz",
            ),
        ],
    )

    url, version = _archive.resolve_interactions_archive("latest")

    assert version == "latest"
    assert url.endswith("omnipath_webservice_interactions__latest.tsv.gz")

    url, version = _archive.resolve_interactions_archive("2023-06-01")

    assert version == "20230601"
    assert url.endswith("omnipath_webservice_interactions__20230502-20240101.tsv.xz")


def test_resolve_interactions_archive_reports_available_dates(monkeypatch):
    monkeypatch.setattr(
        _archive,
        "_list_interactions_archives",
        lambda: [
            (
                "20230101",
                "20230501",
                "omnipath_webservice_interactions__20230101-20230501.tsv.xz",
            ),
            (
                "20230502",
                "20240101",
                "omnipath_webservice_interactions__20230502-20240101.tsv.xz",
            ),
        ],
    )

    with pytest.raises(ValueError) as error:
        _archive.resolve_interactions_archive("2024-06-01")

    assert "no OmniPath interactions archive found for version '20240601'" in str(
        error.value
    )
    assert (
        "available date ranges: 2023-01-01..2023-05-01, 2023-05-02..2024-01-01"
    ) in str(error.value)


def test_resolve_interactions_archive_rejects_current_and_default():
    for version in ["current", "default"]:
        with pytest.raises(ValueError, match="expected 'latest' or a date"):
            _archive.resolve_interactions_archive(version)


def test_dorothea_versions_lists_available_archive_ranges(monkeypatch):
    monkeypatch.setattr(
        _archive,
        "_list_interactions_archives",
        lambda: [
            (
                "20230101",
                "20230501",
                "omnipath_webservice_interactions__20230101-20230501.tsv.xz",
            ),
            (
                "20230502",
                "20230502",
                "omnipath_webservice_interactions__20230502-20230502.tsv.xz",
            ),
        ],
    )

    expected = ["latest", "2023-01-01..2023-05-01", "2023-05-02"]

    assert getattr(_dorothea.dorothea, "versions")() == expected
    assert getattr(_collectri.collectri, "versions")() == expected


def test_load_interactions_version_filters_signed_resource_and_levels(monkeypatch):
    monkeypatch.setattr(
        _archive,
        "resolve_interactions_archive",
        lambda version: ("archive.tsv", "20230601"),
    )

    def read_archive(url):
        return pd.DataFrame(
            {
                "source_genesymbol": ["TfA", "TfB", "TfC", "TfD"],
                "target_genesymbol": ["GeneA", "GeneB", "GeneC", "GeneD"],
                "dorothea": ["True", "True", "True", "False"],
                "collectri": ["False", "False", "False", "True"],
                "is_stimulation": ["True", "False", "True", "True"],
                "is_inhibition": ["False", "True", "True", "False"],
                "consensus_stimulation": ["True", "False", "True", "True"],
                "consensus_inhibition": ["False", "True", "True", "False"],
                "dorothea_level": ["A", "B", "B", "A"],
                "ncbi_tax_id_source": ["9606", "9606", "9606", "9606"],
                "ncbi_tax_id_target": ["9606", "9606", "9606", "9606"],
            }
        )

    monkeypatch.setattr(_archive, "_read_interactions_archive", read_archive)

    dorothea = _archive.load_interactions_version(
        "dorothea",
        version="2023-06-01",
        organism="human",
        levels=["B"],
    )

    assert dorothea.to_dict("records") == [
        {
            "source": "TfB",
            "target": "GeneB",
            "sign": -1,
            "version": "20230601",
            "confidence": "B",
        },
        {
            "source": "TfC",
            "target": "GeneC",
            "sign": 1,
            "version": "20230601",
            "confidence": "B",
        },
    ]

    collectri = _archive.load_interactions_version(
        "collectri",
        version="2023-06-01",
        organism="human",
    )

    assert collectri[["source", "target", "sign", "version"]].to_dict("records") == [
        {
            "source": "TfD",
            "target": "GeneD",
            "sign": 1,
            "version": "20230601",
        }
    ]


def test_load_interactions_version_translates_non_human_organisms(monkeypatch):
    monkeypatch.setattr(
        _archive,
        "resolve_interactions_archive",
        lambda version: ("archive.tsv", "20230601"),
    )
    monkeypatch.setattr(
        _archive,
        "_read_interactions_archive",
        lambda url: pd.DataFrame(
            {
                "source_genesymbol": ["TfB"],
                "target_genesymbol": ["GeneB"],
                "dorothea": ["True"],
                "is_stimulation": ["False"],
                "is_inhibition": ["True"],
                "consensus_stimulation": ["False"],
                "consensus_inhibition": ["True"],
                "dorothea_level": ["B"],
                "ncbi_tax_id_source": ["9606"],
                "ncbi_tax_id_target": ["9606"],
            }
        ),
    )

    calls = []

    class FakeOrthologs:
        def translate_df(self, net, **kwargs):
            calls.append((net.copy(), kwargs))
            translated = net.copy()
            translated["source"] = ["TfB_mouse"]
            translated["target"] = ["GeneB_mouse"]
            return translated

    def orthologs(**kwargs):
        calls.append(kwargs)
        return FakeOrthologs()

    monkeypatch.setattr(_archive, "hcop_orthologs", orthologs)

    result = _archive.load_interactions_version(
        "dorothea",
        version="2023-06-01",
        organism="mouse",
        levels=["B"],
    )

    assert calls[0] == {
        "target_organism": "mouse",
        "version": "latest",
    }
    assert calls[1][1] == {
        "columns": ["source", "target"],
        "one_to_many": 5,
    }
    assert result[["source", "target", "sign"]].to_dict("records") == [
        {"source": "TfB_mouse", "target": "GeneB_mouse", "sign": -1}
    ]


def test_load_interactions_version_uses_requested_hcop_version(monkeypatch):
    monkeypatch.setattr(
        _archive,
        "resolve_interactions_archive",
        lambda version: ("archive.tsv", "20230601"),
    )
    monkeypatch.setattr(
        _archive,
        "_read_interactions_archive",
        lambda url: pd.DataFrame(
            {
                "source_genesymbol": ["TfA"],
                "target_genesymbol": ["GeneA"],
                "dorothea": ["True"],
                "is_stimulation": ["True"],
                "is_inhibition": ["False"],
                "consensus_stimulation": ["True"],
                "consensus_inhibition": ["False"],
                "dorothea_level": ["A"],
                "ncbi_tax_id_source": ["9606"],
                "ncbi_tax_id_target": ["9606"],
            }
        ),
    )

    calls = []

    class FakeOrthologs:
        def translate_df(self, net, **kwargs):
            calls.append((net.copy(), kwargs))
            translated = net.copy()
            translated["source"] = ["TfA_mouse"]
            translated["target"] = ["GeneA_mouse"]
            return translated

    def orthologs(**kwargs):
        calls.append(kwargs)
        return FakeOrthologs()

    monkeypatch.setattr(_archive, "hcop_orthologs", orthologs)

    result = _archive.load_interactions_version(
        "dorothea",
        version="2023-06-01",
        organism="mouse",
        levels=["A"],
        hcop_version="bundled",
    )

    assert calls[0] == {
        "target_organism": "mouse",
        "version": "bundled",
    }
    assert calls[1][1] == {
        "columns": ["source", "target"],
        "one_to_many": 5,
    }
    assert result[["source", "target", "sign"]].to_dict("records") == [
        {"source": "TfA_mouse", "target": "GeneA_mouse", "sign": 1}
    ]


def test_load_interactions_version_uses_archived_mouse_dorothea(monkeypatch):
    monkeypatch.setattr(
        _archive,
        "resolve_interactions_archive",
        lambda version: ("archive.tsv", "20230601"),
    )
    monkeypatch.setattr(
        _archive,
        "_read_interactions_archive",
        lambda url: pd.DataFrame(
            {
                "source_genesymbol": ["TfA", "TfB", "TfC", "TfD"],
                "target_genesymbol": ["GeneA", "GeneB", "TfC", "GeneD"],
                "dorothea": ["True", "True", "True", "True"],
                "is_stimulation": ["False", "False", "False", "False"],
                "is_inhibition": ["False", "False", "False", "False"],
                "consensus_stimulation": ["False", "False", "False", "False"],
                "consensus_inhibition": ["False", "False", "False", "False"],
                "dorothea_level": ["A", "A", "A", "A"],
                "evidences": [
                    {
                        "positive": [
                            {
                                "dataset": "dorothea",
                                "resource": "DoRothEA",
                                "references": ["1"],
                                "via": None,
                            }
                        ],
                        "negative": [],
                        "directed": [],
                        "undirected": [],
                    },
                    {
                        "positive": [
                            {
                                "dataset": "collectri",
                                "resource": "CollecTRI",
                                "references": ["2"],
                                "via": None,
                            }
                        ],
                        "negative": [],
                        "directed": [],
                        "undirected": [],
                    },
                    {
                        "positive": [
                            {
                                "dataset": "dorothea",
                                "resource": "DoRothEA",
                                "references": ["3"],
                                "via": None,
                            }
                        ],
                        "negative": [],
                        "directed": [],
                        "undirected": [],
                    },
                    {
                        "positive": [
                            {
                                "dataset": "dorothea",
                                "resource": "DoRothEA",
                                "references": ["4"],
                                "via": None,
                            }
                        ],
                        "negative": [],
                        "directed": [],
                        "undirected": [],
                    },
                ],
                "ncbi_tax_id_source": ["10090", "10090", "10090", "9606"],
                "ncbi_tax_id_target": ["10090", "10090", "10090", "9606"],
            }
        ),
    )
    monkeypatch.setattr(
        _archive,
        "hcop_orthologs",
        lambda *args, **kwargs: pytest.fail(
            "mouse archive rows should not be translated"
        ),
    )

    result = _archive.load_interactions_version(
        "dorothea",
        version="2023-06-01",
        organism="mouse",
        levels=["A"],
        flavor="legacy",
    )

    assert result[["source", "target", "sign", "confidence"]].to_dict("records") == [
        {"source": "TfA", "target": "GeneA", "sign": 1, "confidence": "A"},
    ]


def test_load_interactions_version_deduplicates_dorothea_edges(monkeypatch):
    monkeypatch.setattr(
        _archive,
        "resolve_interactions_archive",
        lambda version: ("archive.tsv", "20230601"),
    )
    monkeypatch.setattr(
        _archive,
        "_read_interactions_archive",
        lambda url: pd.DataFrame(
            {
                "source_genesymbol": ["TfA", "TfA", "TfA"],
                "target_genesymbol": ["GeneA", "GeneA", "GeneA"],
                "dorothea": ["True", "True", "True"],
                "is_stimulation": ["True", "True", "False"],
                "is_inhibition": ["False", "False", "True"],
                "consensus_stimulation": ["True", "True", "False"],
                "consensus_inhibition": ["False", "False", "True"],
                "dorothea_level": ["A", "A", "A"],
                "ncbi_tax_id_source": ["9606", "9606", "9606"],
                "ncbi_tax_id_target": ["9606", "9606", "9606"],
            }
        ),
    )

    result = _archive.load_interactions_version(
        "dorothea",
        version="2023-06-01",
        organism="human",
        levels=["A"],
    )

    assert result[["source", "target", "sign"]].to_dict("records") == [
        {"source": "TfA", "target": "GeneA", "sign": 1},
        {"source": "TfA", "target": "GeneA", "sign": -1},
    ]


def test_deduplicate_dorothea_keeps_opposite_signs():
    dorothea = pd.DataFrame(
        {
            "source": ["TfA", "TfA", "TfA"],
            "target": ["GeneA", "GeneA", "GeneA"],
            "sign": [1, 1, -1],
            "confidence": ["A", "A", "A"],
        }
    )

    result = _archive._deduplicate_dorothea(dorothea)

    assert result[["source", "target", "sign"]].to_dict("records") == [
        {"source": "TfA", "target": "GeneA", "sign": 1},
        {"source": "TfA", "target": "GeneA", "sign": -1},
    ]


def test_deduplicate_dorothea_can_reproduce_decoupler_compatibility():
    dorothea = pd.DataFrame(
        {
            "source": ["TfA", "TfA", "TfA"],
            "target": ["GeneA", "GeneA", "GeneA"],
            "sign": [1, 1, -1],
            "confidence": ["A", "A", "A"],
        }
    )

    result = _archive._deduplicate_dorothea(dorothea, compatibility=True)

    assert result[["source", "target", "sign"]].to_dict("records") == [
        {"source": "TfA", "target": "GeneA", "sign": -1}
    ]


def test_latest_modern_dorothea_keeps_opposite_signs(monkeypatch):
    monkeypatch.setattr(
        _dorothea.pd,
        "read_csv",
        lambda *args, **kwargs: pd.DataFrame(
            {
                "source_genesymbol": ["TfA", "TfA", "TfA"],
                "target_genesymbol": ["GeneA", "GeneA", "GeneA"],
                "is_stimulation": ["True", "True", "False"],
                "is_inhibition": ["False", "False", "True"],
                "consensus_stimulation": ["True", "True", "False"],
                "dorothea_level": ["A", "A", "A"],
            }
        ),
    )

    result = _dorothea._load_latest_modern_dorothea(
        organism="human",
        levels=["A"],
        hcop_version="bundled",
    )

    assert result[["source", "target", "sign"]].to_dict("records") == [
        {"source": "TfA", "target": "GeneA", "sign": 1},
        {"source": "TfA", "target": "GeneA", "sign": -1},
    ]


def test_latest_modern_dorothea_can_reproduce_decoupler_compatibility(monkeypatch):
    monkeypatch.setattr(
        _dorothea.pd,
        "read_csv",
        lambda *args, **kwargs: pd.DataFrame(
            {
                "source_genesymbol": ["TfA", "TfA", "TfA"],
                "target_genesymbol": ["GeneA", "GeneA", "GeneA"],
                "is_stimulation": ["True", "True", "False"],
                "is_inhibition": ["False", "False", "True"],
                "consensus_stimulation": ["True", "True", "False"],
                "dorothea_level": ["A", "A", "A"],
            }
        ),
    )

    result = _dorothea._load_latest_modern_dorothea(
        organism="human",
        levels=["A"],
        hcop_version="bundled",
        compatibility=True,
    )

    assert result[["source", "target", "sign"]].to_dict("records") == [
        {"source": "TfA", "target": "GeneA", "sign": -1}
    ]


def test_dorothea_uses_omnipath_archive_version(monkeypatch):
    calls = []

    def load_version(
        resource,
        version,
        organism,
        levels=None,
        flavor="modern",
        hcop_version="latest",
        compatibility=False,
    ):
        calls.append(
            (resource, version, organism, levels, flavor, hcop_version, compatibility)
        )
        return pd.DataFrame(
            {
                "source": ["TfA", "TfB"],
                "target": ["GeneA", "GeneB"],
                "sign": [1, -1],
                "confidence": ["A", "B"],
                "version": ["20240101", "20240101"],
            }
        )

    monkeypatch.setattr(_dorothea, "load_interactions_version", load_version)

    grn = _dorothea.dorothea(
        organism=9606,
        levels=["B"],
        version="2024-01-01",
        hcop_version="bundled",
    )

    assert calls == [
        ("dorothea", "2024-01-01", 9606, ["B"], "modern", "bundled", False)
    ]
    assert isinstance(grn, InfluenceGraph)
    assert list(grn.edges(data=True)) == [
        (
            "TfB",
            "GeneB",
            {"confidence": "B", "sign": -1, "version": "20240101"},
        )
    ]


def test_dorothea_supports_legacy_flavor_and_deprecated_wrappers(monkeypatch):
    calls = []

    def load_version(
        resource,
        version,
        organism,
        levels=None,
        flavor="modern",
        hcop_version="latest",
        compatibility=False,
    ):
        calls.append(
            (resource, version, organism, levels, flavor, hcop_version, compatibility)
        )
        return pd.DataFrame(
            {
                "source": ["TfA"],
                "target": ["GeneA"],
                "sign": [1],
                "confidence": ["A"],
            }
        )

    monkeypatch.setattr(_dorothea, "load_interactions_version", load_version)

    _dorothea.dorothea(version="2024-01-01", flavor="legacy")
    with pytest.warns(FutureWarning, match="`wrapper` is deprecated"):
        _dorothea.dorothea(version="2024-01-01", wrapper="get")
    with pytest.warns(FutureWarning, match="`wrapper` is deprecated"):
        _dorothea.dorothea(version="2024-01-01", wrapper="op")

    assert calls == [
        (
            "dorothea",
            "2024-01-01",
            "mouse",
            ["A", "B", "C"],
            "legacy",
            "latest",
            False,
        ),
        (
            "dorothea",
            "2024-01-01",
            "mouse",
            ["A", "B", "C"],
            "legacy",
            "latest",
            False,
        ),
        (
            "dorothea",
            "2024-01-01",
            "mouse",
            ["A", "B", "C"],
            "modern",
            "latest",
            False,
        ),
    ]


def test_collectri_uses_omnipath_archive_version_and_genesyn(monkeypatch):
    class FakeGeneSynonyms:
        def __init__(self):
            self.calls = []

        def __call__(self, graph, **kwargs):
            self.calls.append((graph, kwargs))

    calls = []

    def load_version(resource, version, organism):
        calls.append((resource, version, organism))
        return pd.DataFrame(
            {
                "source": ["TfA"],
                "target": ["GeneA"],
                "sign": [1],
                "version": ["20240101"],
            }
        )

    genesyn = FakeGeneSynonyms()
    monkeypatch.setattr(_collectri, "load_interactions_version", load_version)
    monkeypatch.setattr(_collectri, "GeneSynonyms", FakeGeneSynonyms)

    grn = _collectri.collectri(
        organism="mouse",
        genesyn=genesyn,
        gene_identifier_type="ensembl_id",
        version="2024-01-01",
    )

    assert calls == [("collectri", "2024-01-01", "mouse")]
    assert isinstance(grn, InfluenceGraph)
    assert list(grn.edges(data=True)) == [
        ("TfA", "GeneA", {"sign": 1, "version": "20240101"})
    ]
    assert genesyn.calls == [
        (
            grn,
            {
                "input_type": "name",
                "output_type": "ensembl_id",
                "copy": False,
            },
        )
    ]


def test_omnipath_deprecated_grn_aliases_warn(monkeypatch):
    monkeypatch.setattr(
        _dorothea,
        "dorothea",
        lambda *args, **kwargs: InfluenceGraph(),
    )
    monkeypatch.setattr(
        _collectri,
        "collectri",
        lambda *args, **kwargs: InfluenceGraph(),
    )

    with pytest.warns(FutureWarning, match="load_dorothea_grn"):
        assert isinstance(_dorothea.load_dorothea_grn(), InfluenceGraph)

    with pytest.warns(FutureWarning, match="load_collectri_grn"):
        assert isinstance(_collectri.load_collectri_grn(), InfluenceGraph)


def test_omnipath_loaders_reject_invalid_arguments(monkeypatch):
    monkeypatch.setattr(
        _dorothea,
        "load_interactions_version",
        lambda *args, **kwargs: pd.DataFrame(
            {"source": ["Tf"], "target": ["Gene"], "sign": [1]}
        ),
    )
    monkeypatch.setattr(
        _collectri,
        "load_interactions_version",
        lambda *args, **kwargs: pd.DataFrame(
            {"source": ["Tf"], "target": ["Gene"], "sign": [1]}
        ),
    )

    with pytest.raises(TypeError, match="unsupported argument type for 'organism'"):
        _dorothea.dorothea(organism=cast(Any, object()))

    with pytest.raises(TypeError, match="unsupported argument type for 'genesyn'"):
        _dorothea.dorothea(genesyn=cast(Any, object()), version="latest")

    with pytest.raises(
        TypeError,
        match="unsupported argument type for 'compatibility'",
    ):
        _dorothea.dorothea(compatibility=cast(Any, "yes"))

    with pytest.raises(ValueError, match="invalid argument value for 'flavor'"):
        _dorothea.dorothea(flavor=cast(Any, "old"))

    with pytest.raises(ValueError, match="invalid argument value for 'wrapper'"):
        _dorothea.dorothea(wrapper=cast(Any, "legacy"))

    with pytest.raises(ValueError, match="conflicting DoRothEA"):
        _dorothea.dorothea(flavor="legacy", wrapper="op")

    with pytest.raises(
        TypeError,
        match="unsupported argument type for 'split_complexes'",
    ):
        _collectri.collectri(split_complexes=cast(Any, "yes"))

    with pytest.raises(TypeError, match="unsupported argument type for 'remove_pmid'"):
        _collectri.collectri(remove_pmid=cast(Any, "yes"))

    with pytest.raises(TypeError, match="unsupported argument type for 'organism'"):
        _collectri.collectri(organism=cast(Any, object()))

    with pytest.raises(TypeError, match="unsupported argument type for 'genesyn'"):
        _collectri.collectri(genesyn=cast(Any, object()), version="latest")

    with pytest.raises(ValueError, match="support human, mouse and rat"):
        _archive.load_interactions_version(
            "dorothea",
            version="2023-06-01",
            organism="zebrafish",
        )
