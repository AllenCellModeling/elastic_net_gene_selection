# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

from ..utils.genes import (
    get_gene_set,
    get_selected_genes,
)
from ..datasets.correlated_random_variables import HubSpokeData
from ..solvers.linear import parallel_runs as linear_parallel


def test_get_gene_set(L=100):
    df = pd.DataFrame(
        {
            "Number of genes": np.random.choice([2, 3, 5, 7, 11], L),
            "Gene name": np.random.choice(["a", "b", "c"], L),
        }
    )
    _ = get_gene_set(df, num_genes=10)


def test_get_selected_genes():
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
        lambda_path=np.geomspace(10, 0.01, num=100),
        target_col="day",
        target_map={"0": 0, "1": 1, "2": 2},
    )

    _ = get_selected_genes(
        boot_results,
        adata,
        lambda_index=50,
        selection_threshold_index=75,
        thresholds=np.linspace(0.01, 1, num=100),
    )
