# -*- coding: utf-8 -*-
# all_atom_model.py
# author : Antoine Passemiers

import numpy as np
import scipy.spatial


class AllAtomModel:

    def __init__(self, n_atoms, beta=1.):
        self.n_atoms = n_atoms
        self.beta = beta
        self.triu_indices = np.triu_indices(self.n_atoms, k=-1)

        # Parameters for modelling gaussian restraints
        self.mu = np.full((self.n_atoms, self.n_atoms), np.nan, dtype=np.float)
        self.sigma = np.full((self.n_atoms, self.n_atoms), np.nan, dtype=np.float)

        # Parameters for modelling uniform restraints
        self.lb = np.full((self.n_atoms, self.n_atoms), np.nan, dtype=np.float)

    def add_gaussian_restraint(self, i, j, mu, sigma):
        self.mu[i, j] = self.mu[j, i] = mu
        self.sigma[i, j] = self.sigma[j, i] = sigma

    def add_uniform_restraint(self, i, j, lb):
        self.lb[i, j] = self.lb[j, i] = True

    def evaluate(self, coords):
        distances = scipy.spatial.distance.cdist(coords, coords, metric='euclidean')
        distances = distances[self.triu_indices]
        mu, sigma = self.mu[self.triu_indices], self.sigma[self.triu_indices]
        lb = self.lb[self.triu_indices]

        indices_g = ~np.isnan(mu)
        distances_g, mu, sigma = distances[indices_g], mu[indices_g], sigma[indices_g]

        indices_u = ~np.isnan(lb)
        distances_u, lb = distances[indices_u], lb[indices_u]

        logp = ((distances_g - mu) / sigma) ** 2.
        penalties = (distances_u < lb)
        return -0.5 * logp.sum() + self.beta * penalties.sum()
