# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


def preprocess(
    adata_in,
    nz_thresh=0.05,
    transform=np.arcsinh,
    coding_genes=[
        "2",
        "3",
        "5",
        "7",
        "11",
        "13",
        "17",
        "19",
        "23",
        "29",
        "31",
        "37",
        "41",
        "43",
        "47",
        "53",
        "59",
        "61",
        "67",
        "71",
        "73",
        "79",
        "83",
        "89",
        "97",
    ],
):

    # filter our data for only protein coding genes
    cols = np.array([c for c in adata_in.var.index if c in coding_genes])
    adata = adata_in[:, cols]

    # filter for genes that apoear in at least x% of cells
    gene_nz_freq = (adata.X > 0).mean(axis=0)
    adata = adata[:, cols[gene_nz_freq > nz_thresh]]
    adata.X = transform(adata.X)

    return adata.copy()


def subset(adata, days=["0"]):
    if days == "all":
        return adata
    else:
        return adata[adata.obs["day"].isin(days)].copy()


def split(adata, test_size=0.2, random_state=0):
    """
    Split anndata into train/test and return a dict of the integer indices
    (from the original anndata) in each split and a dict of each split anndata
    """
    all_inds = np.arange(len(adata))
    inds_train, inds_test = train_test_split(
        all_inds, test_size=test_size, random_state=random_state
    )
    split_inds_dict = {"train": sorted(inds_train), "test": sorted(inds_test)}
    split_inds_dict = {k: [int(i) for i in v] for k, v in split_inds_dict.items()}
    return split_inds_dict, {k: adata[v, :] for k, v in split_inds_dict.items()}


def tidy(arr):
    """Take a numpy ndarray and turn it into a tidy dataframe."""
    return pd.DataFrame([(*inds, arr[inds]) for inds in np.ndindex(arr.shape)])
