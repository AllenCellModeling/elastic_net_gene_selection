# -*- coding: utf-8 -*-

import multiprocessing
from functools import partial

import numpy as np

from sklearn.preprocessing import scale
from sklearn.decomposition import PCA

import glmnet_python  # noqa: F401
from glmnet import glmnet

from .utils import get_selected_betas


def enet_pca(
    X,
    n_components=10,
    pc_weights="scaled",
    alpha=0.9,
    penalty_factor=None,
    lambda_path=None,
    **kwargs
):

    # pca and tranform X into top n_components pc coords
    pca = PCA(n_components=n_components, svd_solver="randomized")
    pca.fit(X)
    X_pca = pca.transform(X)

    # weight the pcs or not
    if pc_weights == "variance_explained":
        pass
    elif pc_weights == "scaled":
        X_pca = scale(X_pca)
    else:
        X_pca = pc_weights * scale(X_pca)

    # default enet kwargs
    enet_kwargs = {
        "x": X.copy(),
        "y": X_pca.copy(),
        "alpha": alpha,
        "family": "mgaussian",
    }

    # maybe add unevenly weighted penalties
    if penalty_factor is not None:
        enet_kwargs["penalty_factor"] = penalty_factor

    # maybe specify lambda path
    if lambda_path is not None:
        enet_kwargs["lambdau"] = lambda_path

    # run the regression over the lambda path
    mfit = glmnet(**enet_kwargs)

    # return lambda path and which features are selected at those lambdas
    return {"beta": get_selected_betas(mfit), "lambda_path": mfit["lambdau"]}


def worker(
    boot_inds,
    adata,
    noise=0.001,
    n_pcs=5,
    alpha=0.9,
    lambda_path=np.geomspace(10, 0.01, num=10),
    pc_weights="scaled",
    unpenalized_genes=np.array([]),
):
    X_boot = adata.X[boot_inds, :]
    X_boot = scale(
        scale(X_boot + np.random.normal(scale=noise * 1e-6, size=X_boot.shape))
        + np.random.normal(scale=noise, size=X_boot.shape)
    )
    enet_kwargs = dict(
        n_components=n_pcs,
        alpha=alpha,
        lambda_path=lambda_path,
        pc_weights=pc_weights,
        penalty_factor=np.array(
            [
                gene not in set(unpenalized_genes)
                for gene in [g.split("_")[0] for g in adata.var.index.values]
            ]
        ).astype(float),
    )
    return enet_pca(X_boot, **enet_kwargs)


def parallel_runs(
    adata,
    n_processes=10,
    n_bootstraps=1000,
    noise=0.001,
    n_pcs=5,
    alpha=0.9,
    lambda_path=np.geomspace(10, 0.01, num=10),
    pc_weights="scaled",
    unpenalized_genes=np.array([]),
):

    boot_inds = [
        np.random.choice(len(adata.X), size=len(adata.X)) for _ in range(n_bootstraps)
    ]

    worker_partial = partial(
        worker,
        adata=adata,
        noise=noise,
        n_pcs=n_pcs,
        alpha=alpha,
        lambda_path=lambda_path,
        pc_weights=pc_weights,
        unpenalized_genes=unpenalized_genes,
    )

    pool = multiprocessing.Pool(processes=n_processes)
    result_list = pool.map(worker_partial, boot_inds)
    pool.close()
    pool.join()

    return result_list
