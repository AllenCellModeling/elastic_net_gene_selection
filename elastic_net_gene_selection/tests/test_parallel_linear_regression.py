# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.datasets import make_regression

from ..solvers.linear import parallel_runs
from ..utils.genes import get_thresh_lambda_df


def test_linear_works():
    X, y = make_regression(
        n_samples=1000,
        n_features=100,
        n_informative=10,
        n_targets=1,
    )

    # quantiles for y -> string day labels
    yq = pd.qcut(y, 3, labels=["0", "1", "2"])
    obs = pd.DataFrame({"day": yq})
    obs.index = obs.index.map(str)

    adata = sc.AnnData(X=X, obs=obs)

    _ = parallel_runs(
        adata,
        n_processes=4,
        n_bootstraps=32,
        X_noise=0.01,
        y_noise=0.5,
        alpha=0.9,
        lambda_path=np.geomspace(10, 0.01, num=10),
        target_col="day",
        target_map={"0": 0, "1": 1, "2": 2},
    )


def test_linear_correct():
    LAMBDAS = np.geomspace(1, 0.01, num=25)
    THRESHOLDS = np.linspace(0.01, 1, num=25)

    # make adata
    X, y, coef = make_regression(
        n_samples=1000, n_features=10, n_informative=1, n_targets=1, coef=True
    )
    n_y = 10
    y = pd.qcut(y, n_y, labels=[str(i) for i in range(n_y)])
    obs = pd.DataFrame({"target": y.astype(np.float32)})
    obs.index = obs.index.map(str)
    var = pd.DataFrame({"Gene name": [str(i) for i in range(X.shape[1])], "coef": coef})
    var.index = var.index.astype(str)
    adata = sc.AnnData(X=X, obs=obs, var=var)

    # bootstrap
    out_par = parallel_runs(
        adata,
        n_processes=1,
        n_bootstraps=10,
        lambda_path=LAMBDAS,
        target_col="target",
        target_map={},
    )

    # collate bootstrap results
    df = get_thresh_lambda_df(
        out_par,
        gene_names=adata.var["Gene name"].astype(str).values,
        thresholds=THRESHOLDS,
        lambdas=LAMBDAS,
    )

    # pick midpoint of lambdas + thresholds
    mylambda = LAMBDAS[len(LAMBDAS) // 2]
    mythresh = THRESHOLDS[len(THRESHOLDS) // 2]
    sel_row = df[(df["lambda"] == mylambda) & (df["selection threshold"] == mythresh)]
    true_gene = adata.var[adata.var["coef"] != 0]["Gene name"].item()
    sel_gene = sel_row["selected genes"].item()

    # only one gene picked / same as true gene
    assert sel_row["number of genes selected"].item() == 1
    assert sel_gene == true_gene
