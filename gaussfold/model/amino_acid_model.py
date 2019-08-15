# -*- coding: utf-8 -*-
# amino_acid_model.py
# author : Antoine Passemiers

from gaussfold.constraints.gaussian_constraint import GaussianConstraint

import numpy as np
import scipy.spatial


class AminoAcidModel:
    """Gaussian-restrained 3D model of a protein.

    Attributes:
        mu (np.ndarray): Array of shape (L, L) where element (i, j)
            is the average expected distance between residues i and j.
        sigma (np.ndarray): Array of shape (L, L) where element (i, j)
            is the standard deviation of expected distance between
            residues i and j.
        weights (np.ndarray): Array of shape (L, L) where element (i, j)
            is the weight of restraint (i, j) in the log-likelihood.
        weighted (bool): Whether restraints are weighted
    """

    def __init__(self, weighted=False):
        self._constraints = list()
        self._atom_to_id = dict()
        self._id_to_atom = dict()
        self._initialized = False
        self._weighted = weighted

    def _initialize_matrices(self, n_atoms):
        self._n_atoms = n_atoms
        self._mu = np.full((n_atoms, n_atoms), np.nan, dtype=np.float)
        self._sigma = np.full((n_atoms, n_atoms), np.nan, dtype=np.float)
        self._weights = np.ones((n_atoms, n_atoms), dtype=np.float)
        self._triu_indices = np.triu_indices(n_atoms, k=-1)
        self._tril_indices = np.tril_indices(n_atoms, k=0)
        self._distances = np.empty((n_atoms, n_atoms), dtype=np.float)

    def add_constraint(self, constraint):
        if constraint not in self._constraints:
            self._constraints.append(constraint)
        self._initialized = False

    def initialize(self):
        atoms = set()
        atoms.add(GaussianConstraint.__CENTER_OF_MASS__)
        for constraint in self._constraints:
            for atom in constraint.atoms():
                atoms.add(atom)
        atoms = list(atoms)
        self._atom_to_id = { atom: i for i, atom in enumerate(atoms) }
        self._id_to_atom = { i: atom for i, atom in enumerate(atoms) }

        self._initialize_matrices(len(atoms))

        for constraint in self._constraints:
            atom_a, atom_b = constraint.atoms()
            i = self._atom_to_id[atom_a]
            j = self._atom_to_id[atom_b]
            mu = constraint.mu()
            sigma = constraint.sigma()
            weight = constraint.weight()
            self._add_restraint(i, j, mu, sigma, weight=weight)

        self._initialized = True
        self._weighted = False # TODO
        return self

    def set_coords(self, coords):
        for i in range(len(coords)):
            atom = self._id_to_atom[i]
            atom.set_coords(*coords[i])

    def get_coords(self):
        assert(self._initialized)
        coords = np.empty((self._n_atoms, 3), dtype=np.float)
        for i in range(self._n_atoms):
            coords[i, :] = self._id_to_atom[i].get_coords()
        return np.nan_to_num(coords)

    def _add_restraint(self, i, j, mu, sigma, weight=1.):
        """Adds Gaussian restraint to the model.

        Parameters:
            i (int): Identifier of first residue
            j (int): Identifier of second residue
            mu (float): Average expected distance
            sigma (float): Standard deviation of expected distance
            weight (float): Restraint weight in the log-likelihood
        """
        self._mu[i, j] = self._mu[j, i] = mu
        self._sigma[i, j] = self._sigma[j, i] = sigma
        self._weights[i, j] = self._weights[j, i] = weight
        if weight != 1.:
            self._weighted = True

    def evaluate(self, coords):
        """Computes log-likelihood given the Gaussian parameters `mu` and `sigma`.

        Parameters:
            coords (np.ndarray): Array of shape (L, 3) where ith sub-array represents
                the coordinates of residue i in three-dimensional space.

        Returns:
            float: Log-likelihood of the coordinates given the Gaussian parameters.
        """
        scipy.spatial.distance.cdist(coords, coords, metric='euclidean', out=self._distances)
        distances = self._distances[self._tril_indices]
        mu, sigma = self._mu[self._tril_indices], self._sigma[self._tril_indices]

        indices = ~np.isnan(mu)
        distances, mu, sigma = distances[indices], mu[indices], sigma[indices]

        logp = ((distances - mu) / sigma) ** 2.
        if self._weighted:
            weights = self._weights[self._tril_indices]
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
        scipy.spatial.distance.cdist(coords, coords, metric='euclidean', out=self._distances)
        D = self._distances
        M, S = self._mu, self._sigma

        X = np.tile(coords, (L, 1, 1))
        delta = X - np.transpose(X, (1, 0, 2))

        F = ((D - M) / (D * S ** 2.))[..., np.newaxis]
        F[np.isnan(F)] = 0
        grad = 4. * delta * F
        grad[self._tril_indices] = 0
        grad[np.isnan(self._mu)] = 0

        if self._weighted:
            weights = self._weights[self._tril_indices]
            grad *= weights[..., np.newaxis]
        grad = np.nan_to_num(grad)

        grad = grad.sum(axis=0)
        return grad
