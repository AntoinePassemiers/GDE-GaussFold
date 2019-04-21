# -*- coding: utf-8 -*-
# restraint.py
# author : Antoine Passemiers

import numpy as np
import scipy.spatial


class Model:
    """Gaussian-restrained 3D model of a protein.

    Attributes:
        L (int): Number of residues in the protein.
        mu (np.ndarray): Array of shape (L, L) where element (i, j)
            is the average expected distance between residues i and j.
        sigma (np.ndarray): Array of shape (L, L) where element (i, j)
            is the standard deviation of expected distance between
            residues i and j.
    """

    def __init__(self, L):
        self.L = L
        self.mu = np.full((self.L, self.L), np.nan, dtype=np.float)
        self.sigma = np.full((self.L, self.L), np.nan, dtype=np.float)

    def add_restraint(self, i, j, mu, sigma):
        """Adds Gaussian restraint to the model.

        Parameters:
            i (int): Identifier of first residue
            j (int): Identifier of second residue
            mu (float): Average expected distance
            sigma (float): Standard deviation of expected distance
        """
        self.mu[i, j] = self.mu[j, i] = mu
        self.sigma[i, j] = self.sigma[j, i] = sigma

    def evaluate(self, coords):
        """Computes log-likelihood given the Gaussian parameters `mu` and `sigma`.

        Parameters:
            coords (np.ndarray): Array of shape (L, 3) where ith sub-array represents
                the coordinates of residue i in three-dimensional space.

        Returns:
            float: Log-likelihood of the coordinates given the Gaussian parameters.
        """
        indices = np.triu_indices(self.L, -1)
        distances = scipy.spatial.distance.cdist(coords, coords, metric='euclidean')
        distances, mu, sigma = distances[indices], self.mu[indices], self.sigma[indices]

        indices = ~np.isnan(mu)
        distances, mu, sigma = distances[indices], mu[indices], sigma[indices]

        logp = ((distances - mu) / sigma) ** 2.
        return -0.5 * logp.sum()
