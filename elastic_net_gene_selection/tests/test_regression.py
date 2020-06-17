# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.datasets import make_regression

from ..solvers.linear import parallel_runs as linear_parallel
from ..solvers.pca import parallel_runs as pca_parallel


def test_linear():
    X, y = make_regression(
        n_samples=1000, n_features=100, n_informative=10, n_targets=1,
    )

    # quantiles for y -> string day labels
    yq = pd.qcut(y, 3, labels=["0", "1", "2"])
    obs = pd.DataFrame({"day": yq})
    obs.index = obs.index.map(str)

    adata = sc.AnnData(X=X, obs=obs)

    _ = linear_parallel(
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


def test_pca():
    X, _ = make_regression(
        n_samples=1000, n_features=100, n_informative=10, n_targets=1,
    )

    adata = sc.AnnData(X=X)

    _ = pca_parallel(
        adata,
        n_processes=4,
        n_bootstraps=32,
        noise=0.001,
        n_pcs=2,
        alpha=0.9,
        lambda_path=np.geomspace(10, 0.01, num=10),
        pc_weights="scaled",
        unpenalized_genes=np.array([]),
    )
