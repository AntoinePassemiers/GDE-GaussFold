# -*- coding: utf-8 -*-
# restraint.py
# author : Antoine Passemiers

import numpy as np
import scipy.spatial


class Model:

    def __init__(self, L):
        self.L = L
        self.mu = np.full((self.L, self.L), np.nan, dtype=np.float)
        self.sigma = np.full((self.L, self.L), np.nan, dtype=np.float)

    def add_restraint(self, i, j, mu, sigma):
        self.mu[i, j] = self.mu[j, i] = mu
        self.sigma[i, j] = self.sigma[j, i] = sigma

    def evaluate(self, coords):
        indices = np.triu_indices(self.L, -1)
        distances = scipy.spatial.distance.cdist(coords, coords, metric='euclidean')
        distances, mu, sigma = distances[indices], self.mu[indices], self.sigma[indices]

        indices = ~np.isnan(mu)
        distances, mu, sigma = distances[indices], mu[indices], sigma[indices]

        logp = ((distances - mu) / sigma) ** 2.
        return -0.5 * logp.sum()
