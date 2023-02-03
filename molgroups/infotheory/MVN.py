import numpy as np


def cov_entropy(C):
    """
    Entropy estimate from covariance matrix C
    """
    absdet = np.abs(np.linalg.det(C))
    if absdet > 0:
        return 0.5 * (len(C) * np.log2(2 * np.pi * np.e) + np.log2(absdet))
    else:
        return None


class MVNEntropy(object):
    def __init__(self, x=None, alpha=0.05, max_points=1000, cov=None):
        # max points was 1000
        # compute Mardia test coefficient
        if cov is not None:
            self.C = cov
        else:
            n, p = x.shape  # num points, num dimensions
            self.mu = np.mean(x, axis=0)
            self.C = np.cov(x.T, bias=True) if p > 1 else np.array([[np.var(x.T, ddof=1)]])

    def entropy(self):
        return cov_entropy(self.C)

    def marginal_entropy(self, independent_pars):
        MC = self.C.copy()
        MC = np.delete(MC, independent_pars, 0)
        MC = np.delete(MC, independent_pars, 1)
        return cov_entropy(MC)
