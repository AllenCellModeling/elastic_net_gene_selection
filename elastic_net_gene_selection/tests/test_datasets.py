# -*- coding: utf-8 -*-

from ..datasets.correlated_random_variables import (
    HubSpokeData,
    CorrelatedNormal,
    random_corr_mat,
)


def test_HubSpokeData():
    hubspoke = HubSpokeData()
    _ = hubspoke.sample()


def test_CorrelatedNormal():
    sigma = random_corr_mat()
    corrnorm = CorrelatedNormal(sigma)
    _ = corrnorm.sample()
