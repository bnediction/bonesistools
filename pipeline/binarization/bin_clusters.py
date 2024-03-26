#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

from typing import Optional, Union
from collections import namedtuple

import os, argparse
from pathlib import Path
from utils.argtype import Range
from utils.stdout import Section, disable_print

from pandas import (
    DataFrame,
    Series,
    MultiIndex,
    merge
)
import anndata as ad, anndatatools as adt

import numpy as np

from scboolseq import scBoolSeq

class Predict(object):

    def check_thresholds(
        self,
        nans_threshold: float,
        bimodal_threshold: float,
        zeroinf_threshold: float,
        unimodal_threshold: float
    ):
        if not isinstance(nans_threshold, float):
            raise ValueError("`nans_threshold` argument is not a float.")
        elif not isinstance(bimodal_threshold, float):
            raise ValueError("`bimodal_threshold` argument is not a float.")
        elif not isinstance(zeroinf_threshold, float):
            raise ValueError("`zeroinf_threshold` argument is not a float.")
        elif not isinstance(unimodal_threshold, float):
            raise ValueError("`unimodal_threshold` argument is not a float.")
        elif nans_threshold < 0 or nans_threshold > 1:
            raise ValueError(f"`nans_threshold` value ({nans_threshold}) not in range [0-1].")
        elif bimodal_threshold < 0.5 or bimodal_threshold > 1:
            raise ValueError(f"`bimodal_threshold` value ({bimodal_threshold}) not in range [0.5-1].")
        elif zeroinf_threshold < 0 or zeroinf_threshold > 0.5:
            raise ValueError(f"`zeroinf_threshold` value ({zeroinf_threshold}) not in range [0-0.5].")
        elif unimodal_threshold < 0.5 or unimodal_threshold > 1:
            raise ValueError(f"`unimodal_threshold` value ({unimodal_threshold}) not in range [0.5-1].")
        else:
            pass

    def __init__(
        self,
        nans_threshold: Optional[float] = None,
        bimodal_threshold: Optional[float] = None,
        zeroinf_threshold: Optional[float] = None,
        unimodal_threshold: Optional[float] = None
    ):
        if nans_threshold is None or bimodal_threshold is None or zeroinf_threshold is None or unimodal_threshold is None:
            pass
        else:
            self.check_thresholds(nans_threshold, bimodal_threshold, zeroinf_threshold, unimodal_threshold)
            Threshold = namedtuple("Constants", ["nans", "bimodal", "zeroinf", "unimodal"])
            self.__THRESHOLD = Threshold(
                nans_threshold,
                bimodal_threshold,
                zeroinf_threshold,
                unimodal_threshold
            )
    
    def add(
        self,
        nans_threshold: float,
        bimodal_threshold: float,
        zeroinf_threshold: float,
        unimodal_threshold: float
    ):
        if hasattr(self, f"_{self.__class__.__name__}__THRESHOLD"):
            raise AttributeError(f"'{self.__class__.__name__}' object already attribute '_{self.__class__.__name__}__THRESHOLD'")
        else:
            self.check_thresholds(nans_threshold, bimodal_threshold, zeroinf_threshold, unimodal_threshold)
            Threshold = namedtuple("Constants", ["nans", "bimodal", "zeroinf", "unimodal"])
            self.__THRESHOLD = Threshold(
                nans_threshold,
                bimodal_threshold,
                zeroinf_threshold,
                unimodal_threshold
            )
    
    def get(self):
        if hasattr(self, f"_{self.__class__.__name__}__THRESHOLD"):
            return self.__THRESHOLD
        else:
            return None
    
    def __call__(
        self,
        data: Union[Series, DataFrame],
        category: Optional[str] = None,
    ) -> Union[Series, DataFrame]:

        def boolean_prediction(self, zeros, ones, nans, category):
            denominator = zeros + ones
            total = denominator + nans
            if category=="Bimodal":
                if nans/total > self.__THRESHOLD.nans:
                    return float("nan")
                elif zeros/denominator > self.__THRESHOLD.bimodal:
                    return 0
                elif ones/denominator > self.__THRESHOLD.bimodal:
                    return 1
                else:
                    return float("nan")
            elif category=="ZeroInf":
                if nans/total > self.__THRESHOLD.nans:
                    return float("nan")
                elif ones/denominator > self.__THRESHOLD.zeroinf:
                    return 1
                else:
                    return 0
            elif category=="Unimodal":
                if nans/total > self.__THRESHOLD.nans:
                    return float("nan")
                elif zeros/denominator > self.__THRESHOLD.unimodal:
                    return 0
                elif ones/denominator > self.__THRESHOLD.unimodal:
                    return 1
            else:
                raise ValueError(f"Category argument must be `Bimodal`, `ZeroInf` or `Unimodal`, not `{category}`.")

        if isinstance(data, Series) and category is not None:
            predict_series = Series(index=data.index.get_level_values(0).unique(), name=data._name)
            if category == "Discarded":
                return predict_series
            else:
                for cluster in sorted(data.index.get_level_values(0).unique()):
                    _zeros, _ones, _nans = data[data.index.get_level_values(0) == cluster].droplevel(0)
                    _value = boolean_prediction(self, zeros=_zeros, ones=_ones, nans=_nans, category=category)
                    predict_series[cluster] = _value
                return predict_series
        elif isinstance(data, DataFrame) and category is None:
            predict_df = DataFrame(index=data.index.get_level_values(0).unique())
            for gene in data:
                predict_series = self.__call__(
                    data=data.loc[:,gene],
                    category=scbool.criteria_.loc[gene,"Category"]
                )
                predict_df = predict_df.join(predict_series)
            return predict_df
        else:
            raise ValueError(f"""`data` and `category` arguments must be either of types respectively {Series} and {str}
            or of types respectively {DataFrame} and {type(None)}, not {type(data)} and {type(category)}.""")

def cell_to_cluster_binarization(
    obs_df: DataFrame,
    columns: list,
    group: str,
    dropna: bool = False
) -> DataFrame:

    def counts(
        column_series,
        dropna
    ):
        series = column_series.value_counts(dropna=dropna).to_frame()
        series.index.set_names([column_series.keys, "value"], inplace=True)
        series.rename(columns={"count": column_series._selection}, inplace=True)
        return series

    iterables = (
        sorted(obs_df.loc[:,group].unique()),
        [float(0), float(1), np.nan]
    )
    group_df = DataFrame(index=MultiIndex.from_product(iterables, names=[group, "value"]))

    for column in columns:
        series = counts(obs_df.groupby(by=group)[column], dropna=dropna)
        group_df = group_df.join(series)

    return group_df.fillna(0).astype(int)

parser = argparse.ArgumentParser(
    prog="binarization of sc-RNAseq data",
    description="""From concatenated sc-rnaSeq data recorded in the hdf5 format (<filename>.h5ad), \
        compute cluster-related binarization based on scBoolSeq \
        method (see Magaña López et al. (2023): <https://hal.science/hal-04294917/>).""",
    usage=""""python bin_clusters.py [-h] -i <path> -c <literal> [<args>]"""
)

parser.add_argument(
    "-i", "--infile",
    dest="infile",
    type=lambda x: Path(x).resolve(),
    required=True,
    metavar="PATH",
    help="path to .h5ad file (including file)"
)

parser.add_argument(
    "-o", "--outpath",
    dest="outpath",
    type=lambda x: Path(x).resolve(),
    required=False,
    default=Path("./binarization").resolve(),
    metavar="PATH",
    help="output path (default: ./binarization)"
)

parser.add_argument(
    "-c", "--cluster",
    dest="groupby",
    type=str,
    required=True,
    nargs="+",
    metavar="LITERAL",
    help="clusters retrieving from adata.obs[`cluster`] used for cluster-related binarization"
)

parser.add_argument(
    "-e", "--exclude",
    dest="exclude",
    type=str,
    required=False,
    nargs="+",
    metavar="LITERAL",
    help="cluster names to remove for cluster-related binarization"
)

parser.add_argument(
    "-l", "--layer",
    dest="layer",
    type=str,
    required=False,
    default="log-normalize",
    metavar="LITERAL",
    help="layer used for binarization (default: `log-normalize`)"
)

parser.add_argument(
    "--hvg",
    dest="hvg",
    required=False,
    action="store_true",
    help="select the most variable genes for binarization"
)

parser.add_argument(
    "-n", "--nans-threshold",
    dest="nans_threshold",
    type=float,
    action=Range,
    min=0.,
    max=1.,
    required=False,
    default=0.3,
    help="""set binarized gene value of a cluster to nan if the proportion of nan values
    in the cluster is above `nans_threshold` (default: 0.3)"""
)

parser.add_argument(
    "-b", "--bimodal-threshold",
    dest="bimodal_threshold",
    type=float,
    action=Range,
    min=0.5,
    max=1.,
    required=False,
    default=2/3,
    help="""for a bimodal gene, set binarized gene value of a cluster to 0 (resp. 1)
    if the proportion of zero-values (resp. one-values) in the cluster is above `bimodal_threshold`
    with respect to binarized values (default: 2/3)"""
)

parser.add_argument(
    "-z", "--zeroinf-threshold",
    dest="zeroinf_threshold",
    type=float,
    action=Range,
    min=0.,
    max=0.5,
    required=False,
    default=0.3,
    help="""for a zero-inflated gene, set binarized gene value of a cluster to 1
    if the proportion of one-values in the cluster is above `zeroinf_threshold`,
    otherwise 0 (default: 0.3)"""
)

parser.add_argument(
    "-u", "--unimodal-threshold",
    dest="unimodal_threshold",
    type=float,
    action=Range,
    min=0.5,
    max=1.,
    required=False,
    default=2/3,
    help="""for a unimodal gene, set binarized gene value of a cluster to 0 (resp. 1)
    if the proportion of zero-values (resp. one-values) in the cluster is above `unimodal_threshold`
    with respect to binarized values (default: 2/3)"""
)

parser.add_argument(
    "-v", "--verbose",
    dest="verbose",
    required=False,
    action="store_true",
    help="display information about running programm"
)

args = parser.parse_args()

section = Section(verbose = args.verbose)

scbool = scBoolSeq(
    margin_quantile=0.10,
    zeroinf_binarizer="zero_or_not"
)

predict = Predict(
    args.nans_threshold,
    args.bimodal_threshold,
    args.zeroinf_threshold,
    args.unimodal_threshold
)

if not args.outpath.exists():
    os.makedirs(args.outpath)

print(f"Loading data...")

adata = ad.read_h5ad(args.infile)

adata.obs_names_make_unique()
adata.var_names_make_unique()

if args.hvg is True:
    adata = adata[:,adata.var["highly_variable"]]

gene_list = list(adata.var.index)
counts_df = adt.tl.anndata_to_dataframe(adata, layer=args.layer)

print("Data binarization...")

section("Compute estimators")
with disable_print():
    scbool.fit(counts_df)

section("Estimate boolean values by observation")
with disable_print():
    cell_df = scbool.binarize(counts_df)

section("Estimate boolean values by cluster")
cluster_d = dict()
predict_d = dict()
for _group in args.groupby:
    if len(args.groupby) > 1:
        print(f"\tcomputation for `{_group}` cluster")
    _cell_df = merge(
        cell_df,
        adata.obs.loc[:,_group],
        left_index=True,
        right_index=True,
        how="inner"
    )
    cluster_d[_group] = cell_to_cluster_binarization(
        obs_df=_cell_df,
        columns=gene_list,
        group=_group,
        dropna=False
    )
    if args.exclude:
        _index_label_to_drop = list()
        for _index_label in args.exclude:
            if _index_label in cluster_d[_group].index.get_level_values(0).unique():
                _index_label_to_drop.append(_index_label)
        if _index_label_to_drop:
            cluster_d[_group] = cluster_d[_group].drop(_index_label_to_drop)
        del _index_label_to_drop, _index_label
    predict_d[_group] = predict(cluster_d[_group])

print("Saving data...")

cell_df.to_csv(f"{args.outpath}/cell_bin.csv", sep=",", index=True)
for _group in args.groupby:
    cluster_d[_group].to_csv(f"{args.outpath}/cluster_bin_counts_{_group}.csv", sep=",", index=True)
    predict_d[_group].to_csv(f"{args.outpath}/cluster_bin_{_group}.csv", sep=",", index=True)
