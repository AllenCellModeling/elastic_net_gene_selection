# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scanpy as sc

from scipy.special import expit
from sklearn.datasets import make_regression

from ..solvers.logistic import worker as logistic_serial


def test_logistic_serial():
    X, y = make_regression(n_samples=1000, n_features=10, n_informative=1, n_targets=1)
    y = expit(y / 100)

    obs = pd.DataFrame({"day": y})
    obs["day"] = pd.qcut(obs["day"], 2, labels=["0", "1"])

    obs.index = obs.index.map(str)
    var = pd.DataFrame({"Gene name": [str(i) for i in range(X.shape[1])]})
    var.index = var.index.astype(str)

    adata = sc.AnnData(X=X, obs=obs, var=var)

    _ = logistic_serial(
        np.arange(len(X)),
        adata.X,
        adata.obs["day"],
        X_noise=0.01,
        alpha=0.9,
        lambda_path=np.geomspace(1.0, 0.01, num=100),
    )
