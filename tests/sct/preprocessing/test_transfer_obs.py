#!/usr/bin/env python

import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import bonesistools as bt

pytestmark = pytest.mark.filterwarnings("ignore:Observation names are not unique")


def _adata_with_obs(index, obs):
    return ad.AnnData(
        X=np.zeros((len(index), 1)),
        obs=pd.DataFrame(obs, index=index),
        var=pd.DataFrame(index=["g1"]),
    )


def _integrated_and_specific_adatas():
    ctrl = _adata_with_obs(
        ["cell1", "cell2"],
        {
            "local_score": [1.0, 2.0],
            "local_label": ["ctrl_cell1", "ctrl_cell2"],
        },
    )
    stim = _adata_with_obs(
        ["cell1", "cell3"],
        {
            "local_score": [10.0, 30.0],
            "local_label": ["stim_cell1", "stim_cell3"],
        },
    )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Observation names are not unique")
        integrated = _adata_with_obs(
            ["cell1", "cell2", "cell1", "cell3"],
            {
                "condition": ["ctrl", "ctrl", "stim", "stim"],
                "global_score": [101.0, 102.0, 201.0, 203.0],
                "global_label": [
                    "integrated_ctrl_cell1",
                    "integrated_ctrl_cell2",
                    "integrated_stim_cell1",
                    "integrated_stim_cell3",
                ],
            },
        )

    return integrated, [ctrl, stim]


def test_transfer_obs_to_integrated_handles_colliding_cell_ids():
    integrated, adatas = _integrated_and_specific_adatas()

    result = bt.sct.pp.transfer_obs_to_integrated(
        integrated,
        adatas,
        obs=["local_score", "local_label"],
        conditions=["ctrl", "stim"],
    )

    assert result is None
    assert integrated.obs["condition"].tolist() == ["ctrl", "ctrl", "stim", "stim"]
    assert integrated.obs["local_score"].tolist() == [1.0, 2.0, 10.0, 30.0]
    assert integrated.obs["local_label"].tolist() == [
        "ctrl_cell1",
        "ctrl_cell2",
        "stim_cell1",
        "stim_cell3",
    ]


def test_transfer_obs_to_integrated_copy_does_not_mutate_input():
    integrated, adatas = _integrated_and_specific_adatas()

    copied = bt.sct.pp.transfer_obs_to_integrated(
        integrated,
        adatas,
        obs="local_score",
        conditions=["ctrl", "stim"],
        copy=True,
    )

    assert "local_score" not in integrated.obs
    assert copied is not integrated
    assert copied.obs["local_score"].tolist() == [1.0, 2.0, 10.0, 30.0]


def test_transfer_obs_to_integrated_keeps_unmatched_integrated_cells_as_missing():
    ctrl = _adata_with_obs(
        ["cell1"],
        {"local_score": [1.0]},
    )
    stim = _adata_with_obs(
        ["cell1"],
        {"local_score": [10.0]},
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Observation names are not unique")
        integrated = _adata_with_obs(
            ["cell1", "cell4", "cell1"],
            {"condition": ["ctrl", "ctrl", "stim"]},
        )

    bt.sct.pp.transfer_obs_to_integrated(
        integrated,
        [ctrl, stim],
        obs="local_score",
        conditions=["ctrl", "stim"],
    )

    assert integrated.obs["local_score"].tolist()[0] == 1.0
    assert pd.isna(integrated.obs["local_score"].tolist()[1])
    assert integrated.obs["local_score"].tolist()[2] == 10.0


def test_transfer_obs_to_integrated_supports_custom_condition_column():
    integrated, adatas = _integrated_and_specific_adatas()
    integrated.obs = integrated.obs.rename(columns={"condition": "sample"})

    bt.sct.pp.transfer_obs_to_integrated(
        integrated,
        adatas,
        obs="local_score",
        conditions=["ctrl", "stim"],
        condition_colname="sample",
    )

    assert integrated.obs["sample"].tolist() == ["ctrl", "ctrl", "stim", "stim"]
    assert integrated.obs["local_score"].tolist() == [1.0, 2.0, 10.0, 30.0]


def test_transfer_obs_to_specific_matches_conditions():
    integrated, adatas = _integrated_and_specific_adatas()

    result = bt.sct.pp.transfer_obs_to_specific(
        integrated,
        adatas,
        obs=["global_score", "global_label"],
        conditions=["ctrl", "stim"],
    )

    assert result is None
    ctrl, stim = adatas
    assert ctrl.obs["global_score"].tolist() == [101.0, 102.0]
    assert ctrl.obs["global_label"].tolist() == [
        "integrated_ctrl_cell1",
        "integrated_ctrl_cell2",
    ]
    assert stim.obs["global_score"].tolist() == [201.0, 203.0]
    assert stim.obs["global_label"].tolist() == [
        "integrated_stim_cell1",
        "integrated_stim_cell3",
    ]


def test_transfer_obs_to_specific_copy_does_not_mutate_inputs():
    integrated, adatas = _integrated_and_specific_adatas()

    copied_adatas = bt.sct.pp.transfer_obs_to_specific(
        integrated,
        adatas,
        obs="global_score",
        conditions=["ctrl", "stim"],
        copy=True,
    )

    assert all("global_score" not in original.obs for original in adatas)
    assert copied_adatas is not adatas
    assert [copied.obs["global_score"].tolist() for copied in copied_adatas] == [
        [101.0, 102.0],
        [201.0, 203.0],
    ]


def test_transfer_obs_to_specific_keeps_unmatched_specific_cells_as_missing():
    integrated, adatas = _integrated_and_specific_adatas()
    adatas[0] = _adata_with_obs(["cell1", "cell_missing"], {})

    bt.sct.pp.transfer_obs_to_specific(
        integrated,
        adatas,
        obs="global_score",
        conditions=["ctrl", "stim"],
    )

    assert adatas[0].obs["global_score"].tolist()[0] == 101.0
    assert pd.isna(adatas[0].obs["global_score"].tolist()[1])


def test_transfer_obs_to_specific_supports_custom_condition_column():
    integrated, adatas = _integrated_and_specific_adatas()
    integrated.obs = integrated.obs.rename(columns={"condition": "sample"})

    bt.sct.pp.transfer_obs_to_specific(
        integrated,
        adatas,
        obs="global_score",
        conditions=["ctrl", "stim"],
        condition_colname="sample",
    )

    assert adatas[0].obs["global_score"].tolist() == [101.0, 102.0]
    assert adatas[1].obs["global_score"].tolist() == [201.0, 203.0]


def test_transfer_obs_to_specific_raises_when_condition_column_is_missing():
    integrated, adatas = _integrated_and_specific_adatas()

    with pytest.raises(KeyError):
        bt.sct.pp.transfer_obs_to_specific(
            integrated,
            adatas,
            obs="global_score",
            conditions=["ctrl", "stim"],
            condition_colname="missing",
        )


def test_transfer_obs_to_integrated_raises_when_requested_obs_column_is_missing():
    integrated, adatas = _integrated_and_specific_adatas()

    with pytest.raises(KeyError):
        bt.sct.pp.transfer_obs_to_integrated(
            integrated,
            adatas,
            obs="missing",
            conditions=["ctrl", "stim"],
        )
