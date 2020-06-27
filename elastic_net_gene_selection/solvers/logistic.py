# -*- coding: utf-8 -*-

import multiprocessing
from functools import partial

import numpy as np

from sklearn.preprocessing import scale
from glmnet import LogitNet


def worker(
    boot_inds,
    X,
    y,
    X_noise=0.01,
    alpha=0.9,
    lambda_path=np.geomspace(1e-3, 1e-06, num=100),
):

    X_boot = X[boot_inds, :]
    y_boot = y[boot_inds]

    X_boot = scale(
        scale(X_boot + np.random.normal(scale=X_noise * 1e-6, size=X_boot.shape))
        + np.random.normal(scale=X_noise, size=X_boot.shape)
    )

    m = LogitNet(alpha=alpha, lambda_path=lambda_path, fit_intercept=False,)
    m.fit(X_boot, y_boot)

    lambdas_enet = m.lambda_path_
    coefs_enet = m.coef_path_.squeeze()

    return {
        "beta": coefs_enet != 0,
        "lambda_path": lambdas_enet,
    }


def parallel_runs(
    adata,
    n_processes=10,
    n_bootstraps=1000,
    X_noise=0.01,
    alpha=0.9,
    lambda_path=np.geomspace(1e-3, 1e-06, num=100),
    target_col="timepoint_coarse",
):

    boot_inds = [
        np.random.choice(len(adata.X), size=len(adata.X)) for _ in range(n_bootstraps)
    ]

    X = adata.X
    y = adata.obs[target_col].values

    worker_partial = partial(
        worker, X=X, y=y, X_noise=X_noise, alpha=alpha, lambda_path=lambda_path,
    )

    pool = multiprocessing.Pool(processes=n_processes)
    result_list = pool.map(worker_partial, boot_inds)
    pool.close()
    pool.join()

    return result_list
