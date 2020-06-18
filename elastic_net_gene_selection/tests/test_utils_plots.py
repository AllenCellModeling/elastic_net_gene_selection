# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.datasets import make_regression

from ..utils.plots import thresh_lambda_plot, hub_persistence_plot
from ..datasets.correlated_random_variables import HubSpokeData
from ..solvers.linear import parallel_runs as linear_parallel


def test_thresh_lambda_plot():
    X, y = make_regression(
        n_samples=1000, n_features=100, n_informative=10, n_targets=1,
    )

    # quantiles for y -> string day labels
    yq = pd.qcut(y, 3, labels=["0", "1", "2"])
    obs = pd.DataFrame({"day": yq})
    obs.index = obs.index.map(str)

    adata = sc.AnnData(X=X, obs=obs)

    boot_results = linear_parallel(
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

    _ = thresh_lambda_plot(boot_results, adata)


def test_hub_persistence_plot():
    hubspoke = HubSpokeData()
    adata = hubspoke.sample()
    adata.obs = pd.DataFrame({"day": np.random.choice(["0", "1", "2"], len(adata))})

    boot_results = linear_parallel(
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

    hub_persistence_plot(adata, boot_results)
