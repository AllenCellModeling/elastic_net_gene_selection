import numpy as np
import pandas as pd
import scanpy.api as sc


def random_corr_mat(D=10, beta=1):
    """Generate random valid correlation matrix of dimension D.
    Smaller beta gives larger off diagonal correlations (beta > 0)."""

    P = np.zeros([D, D])
    S = np.eye(D)

    for k in range(0, D - 1):
        for i in range(k + 1, D):
            P[k, i] = 2 * np.random.beta(beta, beta) - 1
            p = P[k, i]
            for l in reversed(range(k)):
                p = (
                    p * np.sqrt((1 - P[l, i] ** 2) * (1 - P[l, k] ** 2))
                    + P[l, i] * P[l, k]
                )
            S[k, i] = S[i, k] = p

    p = np.random.permutation(D)
    for i in range(D):
        S[:, i] = S[p, i]
    for i in range(D):
        S[i, :] = S[i, p]
    return S


def hub_spoke_corr_mat(D=50, groups=5, v=0.3, u=0.1):
    """Port of data generation code from Wasserman's Huge package."""
    from statsmodels.stats.moment_helpers import cov2corr

    G = D // groups  # group size
    Theta = np.zeros([D, D])

    for g in range(groups):
        for i in range(G):
            Theta[g * G, g * G + i] = Theta[g * G + i, g * G] = 1

    Theta[np.diag_indices(D)] = 0
    Omega = Theta * v
    Omega[np.diag_indices(D)] = np.abs(np.min(np.linalg.eigvals(Omega))) + 0.1 + u
    Sigma = cov2corr(np.linalg.inv(Omega))
    Omega = np.linalg.inv(Sigma)

    return Omega, Sigma


class HubSpokeData:
    def __init__(
        self,
        n_groups=5,
        group_size=10,
        n_singeltons=50,
        diagonal_weight=1 / np.e,
        off_diagonal_weight=1,
    ):

        D_small = n_groups * group_size
        _, Sigma_small = hub_spoke_corr_mat(
            D=D_small, groups=n_groups, v=off_diagonal_weight, u=diagonal_weight
        )
        D = n_groups * group_size + n_singeltons
        Sigma = np.eye(D, D)
        Sigma[:D_small, :D_small] = Sigma_small

        var = pd.DataFrame(
            {
                "Type": (["hub"] + ["spoke"] * (group_size - 1)) * n_groups
                + ["singleton"] * n_singeltons,
                "Group": [
                    i
                    for s in ([g + 1] * group_size for g in range(n_groups))
                    for i in s
                ]
                + [0] * n_singeltons,
            }
        )

        self.Sigma = Sigma
        self.var = var

    def sample(self, n_samples=1000):
        X = np.random.multivariate_normal(
            mean=np.zeros(len(self.Sigma)), cov=self.Sigma, size=n_samples
        )
        return sc.AnnData(X=X, var=self.var)


class CorrelatedNormal:
    def __init__(self, Sigma):

        self.Sigma = Sigma

    def sample(self, n_samples=1000):
        X = np.random.multivariate_normal(
            mean=np.zeros(len(self.Sigma)), cov=self.Sigma, size=n_samples
        )
        return sc.AnnData(X=X)
