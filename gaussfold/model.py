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
        weights (np.ndarray): Array of shape (L, L) where element (i, j)
            is the weight of restraint (i, j) in the log-likelihood.
        weighted (bool): Whether restraints are weighted
    """

    def __init__(self, L):
        self.L = L
        self.mu = np.full((self.L, self.L), np.nan, dtype=np.float)
        self.sigma = np.full((self.L, self.L), np.nan, dtype=np.float)
        self.weights = np.ones((self.L, self.L), dtype=np.float)
        self.triu_indices = np.triu_indices(self.L, k=-1)
        self.tril_indices = np.tril_indices(self.L, k=0)
        self.weighted = False

    def add_restraint(self, i, j, mu, sigma, weight=1.):
        """Adds Gaussian restraint to the model.

        Parameters:
            i (int): Identifier of first residue
            j (int): Identifier of second residue
            mu (float): Average expected distance
            sigma (float): Standard deviation of expected distance
            weight (float): Restraint weight in the log-likelihood
        """
        self.mu[i, j] = self.mu[j, i] = mu
        self.sigma[i, j] = self.sigma[j, i] = sigma
        self.weights[i, j] = self.weights[j, i] = weight
        if weight != 1.:
            self.weighted = True

    def evaluate(self, coords):
        """Computes log-likelihood given the Gaussian parameters `mu` and `sigma`.

        Parameters:
            coords (np.ndarray): Array of shape (L, 3) where ith sub-array represents
                the coordinates of residue i in three-dimensional space.

        Returns:
            float: Log-likelihood of the coordinates given the Gaussian parameters.
        """
        distances = scipy.spatial.distance.cdist(coords, coords, metric='euclidean')
        distances = distances[self.triu_indices]
        mu, sigma = self.mu[self.triu_indices], self.sigma[self.triu_indices]

        indices = ~np.isnan(mu)
        distances, mu, sigma = distances[indices], mu[indices], sigma[indices]

        logp = ((distances - mu) / sigma) ** 2.
        if self.weighted:
            weights = self.weights[self.triu_indices]
            logp *= weights[indices]
        return -0.5 * logp.sum()

    def gradient(self, coords):
        """Computes gradient of log-likelihood given the Gaussian parameters `mu` and `sigma`,
        with respect to 3D coordinates.

        Parameters:
            coords (np.ndarray): Array of shape (L, 3) where ith sub-array represents
                the coordinates of residue i in three-dimensional space.

        Returns:
            np.ndarray: Array of shape (L, 3) representing the log-likelihood gradient
                with respect to 3D coordinates.
        """
        L = coords.shape[0]
        D = scipy.spatial.distance.cdist(coords, coords, metric='euclidean')
        M, S = self.mu, self.sigma


        X = np.tile(coords, (L, 1, 1))
        delta = X - np.transpose(X, (1, 0, 2))

        F = ((D - M) / (D * S ** 2.))[..., np.newaxis]
        F[np.isnan(F)] = 0
        grad = 4. * delta * F
        grad[self.tril_indices] = 0
        grad[np.isnan(self.mu)] = 0

        if self.weighted:
            grad *= weights[..., np.newaxis]
        grad = np.nan_to_num(grad)

        grad = grad.sum(axis=0)
        return grad
