# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

from ..utils.data import preprocess, subset, split, tidy
from ..datasets.correlated_random_variables import HubSpokeData


def test_preprocess():
    hubspoke = HubSpokeData()
    adata = hubspoke.sample()
    _ = preprocess(adata)


def test_subset():
    hubspoke = HubSpokeData()
    adata = hubspoke.sample()
    adata.obs = pd.DataFrame({"day": np.random.choice(["0", "1", "2"], len(adata))})
    for days in ["all", ["0"]]:
        _ = subset(adata, days=days)


def test_split():
    hubspoke = HubSpokeData()
    adata = hubspoke.sample()
    _ = split(adata)


def test_tidy():
    hubspoke = HubSpokeData()
    adata = hubspoke.sample()
    tidy(adata.X)
