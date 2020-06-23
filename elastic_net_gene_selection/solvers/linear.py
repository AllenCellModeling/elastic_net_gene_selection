# -*- coding: utf-8 -*-

import multiprocessing
from functools import partial

import numpy as np

from sklearn.preprocessing import scale
from sklearn.linear_model import enet_path


def worker(
    boot_inds,
    X,
    y,
    X_noise=0.01,
    y_noise=0.5,
    alpha=0.9,
    lambda_path=np.geomspace(1.0, 0.01, num=100),
):

    X_boot = X[boot_inds, :]
    y_boot = y[boot_inds]

    X_boot = scale(
        scale(X_boot + np.random.normal(scale=X_noise * 1e-6, size=X_boot.shape))
        + np.random.normal(scale=X_noise, size=X_boot.shape)
    )
    y_boot = scale(
        scale(y_boot + np.random.normal(scale=y_noise * 1e-6, size=len(y_boot)))
        + np.random.normal(scale=y_noise, size=len(y_boot))
    )

    lambdas_enet, coefs_enet, _ = enet_path(
        X_boot, y_boot, l1_ratio=alpha, alphas=lambda_path, fit_intercept=False
    )

    return {
        "beta": coefs_enet != 0,
        "lambda_path": lambdas_enet,
    }


def parallel_runs(
    adata,
    n_processes=10,
    n_bootstraps=1000,
    X_noise=0.01,
    y_noise=0.5,
    alpha=0.9,
    lambda_path=np.geomspace(10, 0.01, num=10),
    target_col="day",
    target_map={"0": 0, "1": 1, "2": 2},
):

    boot_inds = [
        np.random.choice(len(adata.X), size=len(adata.X)) for _ in range(n_bootstraps)
    ]

    X = adata.X
    y = np.array([target_map[d] for d in adata.obs[target_col]])

    worker_partial = partial(
        worker,
        X=X,
        y=y,
        X_noise=X_noise,
        y_noise=y_noise,
        alpha=alpha,
        lambda_path=lambda_path,
    )

    pool = multiprocessing.Pool(processes=n_processes)
    result_list = pool.map(worker_partial, boot_inds)
    pool.close()
    pool.join()

    return result_list
