import numpy
from sklearn.mixture import BayesianGaussianMixture as GMM


class GMMEntropy(object):
    def __init__(self, x):
        self.n_components = int(5 * numpy.sqrt(x.shape[1]))
        self.gmmpredictor = GMM(n_components=self.n_components, max_iter=10000)
        self.gmmpredictor.fit(x)

    def entropy(self, n_est):
        sample, _ = self.gmmpredictor.sample(n_est)
        ll = self.gmmpredictor.score_samples(sample)
        h = -numpy.mean(ll) / numpy.log(2)
        return h

    def score_samples(self, x):
        return self.gmmpredictor.score_samples(x)