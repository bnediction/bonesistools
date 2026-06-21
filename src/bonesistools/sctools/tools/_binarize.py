#!/usr/bin/env python

from __future__ import annotations

from typing import Any, Optional, Sequence, Union, cast

import numpy as np
import pandas as pd
from anndata import AnnData

from ..._compat import Literal
from ..._validation import _as_literal, _as_non_negative_number, _as_probability
from .._typing import VarSubset
from .._validation import _as_var_subset
from ._markers import DEA_METHODS, DEAMethod, dea
from ._stats import CORRECTION_METHODS, CorrectionMethod


class DEABinarizer:
    """
    DEA-based binarization derives Boolean abstractions directly at the group
    level by comparing each population against a reference population using
    differential expression analysis. Unlike cell-level binarization approaches
    followed by voting procedures, this strategy operates on aggregated
    expression profiles, reducing the impact of technical dropouts and focusing
    on transcriptional differences that distinguish groups. Consequently,
    genes exhibiting similar expression across all groups remain undefined,
    whereas significantly enriched or depleted genes are assigned active or
    inactive states according to the sign of their log2 fold-change.

    The model first computes a complete differential expression table during
    `fit()`. The Boolean abstraction is then derived during `binarize()` by
    applying adjusted p-value and log2 fold-change thresholds.

    Attributes
    ----------
    dea_: pandas.DataFrame
        Complete unfiltered DEA table.
    abstraction_: pandas.DataFrame
        Partial Boolean abstraction with groups as rows and features as
        columns. Values are 1, 0 or NaN.
    groups_: pandas.Index
        Groups used as macrostates.
    features_: pandas.Index
        Features used in the abstraction.
    """

    dea_: pd.DataFrame
    abstraction_: pd.DataFrame
    groups_: pd.Index
    features_: pd.Index
    method: DEAMethod
    correction: Optional[CorrectionMethod]
    alpha: Optional[float]
    min_abs_logfoldchange: float
    max_memory: Optional[Union[int, str]]
    obs_: str
    expression_: Optional[str]
    is_log_: bool

    def __init__(
        self,
        method: DEAMethod = "wilcoxon",
        correction: Optional[CorrectionMethod] = "benjamini-hochberg",
        alpha: Optional[float] = 0.05,
        min_abs_logfoldchange: float = 0.5,
        max_memory: Optional[Union[int, str]] = None,
    ) -> None:

        self.method = cast(
            DEAMethod,
            _as_literal(method, choices=DEA_METHODS, name="method"),
        )
        self.correction = cast(
            Optional[CorrectionMethod],
            _as_literal(
                correction,
                choices=CORRECTION_METHODS,
                name="correction",
                allow_none=True,
            ),
        )
        self.alpha = None if alpha is None else _as_probability(alpha, "alpha")
        self.min_abs_logfoldchange = _as_non_negative_number(
            min_abs_logfoldchange,
            "min_abs_logfoldchange",
        )
        self.max_memory = max_memory

    def fit(
        self,
        adata: AnnData,
        obs: str,
        background: Union[Literal["rest"], Sequence[Any]] = "rest",
        expression: Optional[str] = None,
        is_log: bool = False,
        var_subset: Optional[Union[str, Sequence[str]]] = None,
    ) -> "DEABinarizer":
        """
        Compute the complete DEA table.

        The stored DEA table is not filtered by `alpha` or log-fold-change.
        Filtering is applied by `binarize()`.

        Parameters
        ----------
        adata: AnnData
            Unimodal annotated data matrix.
        obs: str
            Observation column in `adata.obs` defining groups.
        background: 'rest' or sequence (default: "rest")
            Background population passed to `dea`.
        expression: str, optional
            Expression source. If None or `"X"`, use `adata.X`. If `"raw.X"`,
            use `adata.raw.X`. Otherwise, interpret as a layer key in
            `adata.layers`.
        is_log: bool (default: False)
            Whether selected expression values are already log1p-transformed.
        var_subset: str or collection of str, optional
            Variables to test and include in the abstraction.

        Returns
        -------
        DEABinarizer
            Fitted binarizer.
        """

        self.obs_ = obs
        self.expression_ = expression
        self.is_log_ = is_log
        self.dea_ = dea(
            adata,
            groupby=obs,
            groups="all",
            background=background,
            method=self.method,
            expression=expression,
            is_log=is_log,
            var_subset=var_subset,
            correction=self.correction,
            alpha=None,
            filter_logfoldchanges=None,
            max_memory=self.max_memory,
        )
        self.groups_ = pd.Index(sorted(adata.obs[obs].dropna().unique()))
        self.features_ = self._resolve_features(adata, var_subset)
        return self

    def binarize(self) -> pd.DataFrame:
        """
        Build a partial Boolean abstraction from the fitted DEA table.

        A feature is assigned:

        - 1 when significantly up-regulated;
        - 0 when significantly down-regulated;
        - NaN otherwise.

        Returns
        -------
        pandas.DataFrame
            Partial Boolean abstraction with groups as rows and features as
            columns.

        Raises
        ------
        RuntimeError
            If the binarizer has not been fitted.
        """

        if not hasattr(self, "dea_"):
            raise RuntimeError("DEABinarizer must be fitted before calling binarize().")

        dea_df = self.dea_.copy()
        if self.alpha is not None:
            dea_df = dea_df.loc[dea_df["pvals_adj"] <= self.alpha]
        dea_df = dea_df.loc[
            np.abs(dea_df["logfoldchanges"]) >= self.min_abs_logfoldchange
        ]

        abstraction = pd.DataFrame(
            data=np.nan,
            index=self.groups_,
            columns=self.features_,
        )
        for group, feature, logfoldchange in dea_df[
            ["group", "feature", "logfoldchanges"]
        ].itertuples(index=False, name=None):
            if group in abstraction.index and feature in abstraction.columns:
                abstraction.at[group, feature] = 1 if logfoldchange > 0 else 0

        self.abstraction_ = abstraction
        return abstraction

    def fit_binarize(
        self,
        adata: AnnData,
        obs: str,
        background: Union[Literal["rest"], Sequence[Any]] = "rest",
        expression: Optional[str] = None,
        is_log: bool = False,
        var_subset: Optional[Union[str, Sequence[str]]] = None,
    ) -> pd.DataFrame:
        """
        Fit the binarizer and return the Boolean abstraction.

        Parameters
        ----------
        adata: AnnData
            Unimodal annotated data matrix.
        obs: str
            Observation column in `adata.obs` defining groups.
        background: 'rest' or sequence (default: "rest")
            Background population passed to `dea`.
        expression: str, optional
            Expression source.
        is_log: bool (default: False)
            Whether selected expression values are already log1p-transformed.
        var_subset: str or collection of str, optional
            Variables to test and include in the abstraction.

        Returns
        -------
        pandas.DataFrame
            Partial Boolean abstraction.
        """

        return self.fit(
            adata=adata,
            obs=obs,
            background=background,
            expression=expression,
            is_log=is_log,
            var_subset=var_subset,
        ).binarize()

    @staticmethod
    def _resolve_features(
        adata: AnnData,
        var_subset: VarSubset,
    ) -> pd.Index:

        mask = _as_var_subset(adata, var_subset)
        if mask is None:
            return pd.Index(adata.var_names)
        return pd.Index(adata.var_names[mask])
