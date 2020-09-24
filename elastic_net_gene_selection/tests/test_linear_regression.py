# -*- coding: utf-8 -*-

import numpy as np
from sklearn.datasets import make_regression

from ..solvers.linear import worker as linear_serial


def test_linear_serial():
    X, y = make_regression(
        n_samples=1000,
        n_features=100,
        n_informative=10,
        n_targets=1,
    )

    _ = linear_serial(
        np.arange(len(X)),
        X,
        y,
        X_noise=0.01,
        y_noise=0.5,
        alpha=0.9,
        lambda_path=np.geomspace(1.0, 0.01, num=100),
    )
