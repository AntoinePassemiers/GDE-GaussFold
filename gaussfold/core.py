# -*- coding: utf-8 -*-
# core.py: GDE-GaussFold core algorithm
# author : Antoine Passemiers

from gaussfold.corrector import DeviationCorrector
from gaussfold.graph import Graph
from gaussfold.metrics import tm_score
from gaussfold.model import Model
from gaussfold.optimizer import Optimizer

import numpy as np
import random
from sklearn.manifold import MDS


class GaussFold:
    """GDE-GaussFold base class.

    Attributes:
        sep (int): Minimum sequence separation, excluding
            backbone contacts.
        n_runs (int): Number of times to run Multi-Dimensional
            Scaling algorithm from scratch.
        max_n_iter (int): Maximum number of iterations of
            Multi-Dimensional Scaling algorithm.
        eps (float): Convergence threshold of Multi-Dimensional
            Scaling algorithm.
    """

    def __init__(self, sep=1, n_init_sols=20, n_runs=1, max_n_iter=300, eps=1e-3):
        self.sep = sep
        self.n_init_sols = n_init_sols
        self.n_runs = n_runs
        self.max_n_iter = max_n_iter
        self.eps = eps
        self._model = None
        self._optimizer = None

    def run(self, cmap, ssp, verbose=True):
        """Runs GDE-GaussFold algorithm.

        Parameters:
            cmap (:obj:`np.ndarray`): Array of shape (L, L) representing
                predicted contact probabilities, where L is the number of
                residues in the sequence.
            ssp (:obj:`np.ndarray`): Array of shape (L,) representing
                3-state secondary structure prediction.
                0 stands for 'H', 1 for 'E' and 2 for 'C'.
            verbose (bool): Whether to display messages in stdout.

        Returns:
            :obj:`np.ndarray`: Array of shape (L, 3) representing the
                protein in the 3D space.
        """
        L = len(cmap)

        # Set diagonal to zeros
        cmap = np.asarray(cmap)
        cmap[np.isnan(cmap)] = 0.
        np.fill_diagonal(cmap, 0)

        # Choose threshold such that exactly 4.5*L contacts
        # are obtained
        L = len(cmap)
        proba = cmap[np.triu_indices(L, -self.sep)]
        proba.sort()
        n_top = int(np.round(4.5 * L))
        threshold = proba[-n_top]

        A = (cmap > threshold)
        G = Graph(A)
        gds = G.distances()
        missing = np.zeros(L, dtype=np.bool)
        if not G.is_connected():
            missing = (gds == 0).all(axis=0)
            missing_idx = np.where(missing)[0]
            connected_idx = np.random.choice(
                    np.arange(L), size=len(missing_idx), replace=False)
            A[missing_idx, connected_idx] = 1
            A[connected_idx, missing_idx] = 1
            G = Graph(A)

            """
            gds[missing_idx, :] = G.distances()[missing_idx, :]
            gds[:, missing_idx] = gds[missing_idx, :].T
            missing_norm = gds[missing, :].mean()
            nonmissing_norm = gds[~missing, :].mean()
            gds = gds.astype(np.float)
            gds[missing_idx, :] *= (nonmissing_norm / missing_norm)
            gds[:, missing_idx] = gds[missing_idx, :].T
            """
            

        # Compute confidence indexes
        #weights = cmap - threshold
        #weights[weights < 0.] = 0.
        #weights /= np.max(weights)
        weights = np.ones((L, L), dtype=np.float)
        weights[missing, :] = 0.
        weights[:, missing] = 0.

        # Apply theoretical linear correspondence between graph
        # distance and Angstroms distance based on statistical
        # observations on euclidean distances found for a graph
        # distance of 1.
        # The empirical value of 4.846 could be used as well.
        distances = gds * 5.72

        if verbose:
            print('Apply Multi-Dimensional Scaling algorithm')

        init_solutions = list()
        for k in range(self.n_init_sols):
            # Multi-dimensional scaling to obtain approximate
            # 3D coordinates
            embedding = MDS(
                    n_components=3,
                    metric=True,
                    n_init=self.n_runs,
                    max_iter=self.max_n_iter,
                    eps=self.eps,
                    n_jobs=None,
                    random_state=None,
                    dissimilarity='precomputed')
            X_transformed = embedding.fit_transform(distances)

            # Apply correction on pairs of adjacent residues
            # based on known C_alpha-C_alpha (or C_beta-C_beta) distance
            if verbose:
                print('Apply deviation correction')
            corrector = DeviationCorrector(len(distances))
            try:
                X_transformed = corrector.fit_transform(X_transformed)
            except np.linalg.linalg.LinAlgError:
                if verbose:
                    print('[Warning] Invalid value encountered in deviation corrector')
            init_solutions.append(X_transformed)

        # Align all initial solution in 3D space
        if len(init_solutions) > 1:
            for i in range(1, len(init_solutions)):
                init_solutions[i] = tm_score(
                    init_solutions[i], init_solutions[0], return_coords=True)[1]

        # Create Gaussian model if not set by the user
        if not isinstance(self._model, Model):
            if verbose:
                print('Model not set by user. Creating model from scratch...')
            self._model = self.create_model(gds, ssp, weights)

        # Create optimizer if not set by the user.
        # Use default hyper-parameters.
        if not isinstance(self._optimizer, Optimizer):
            if verbose:
                print('Optimizer not set by user. Using default parameters.')
            self._optimizer = Optimizer()
        
        # Run optimizer on the Gaussian model
        best_coords = self._optimizer.run(
                init_solutions, self._model, verbose=verbose)
        return best_coords

    def create_model(self, gds, ssp, weights):
        """Creates a Gaussian model for the protein.

        Parameters:
            gds (:obj:`np.ndarray`): Matrix of graph distance
                between each pair of residues in the protein.
            ssp (:obj:`np.ndarray`): Array of shape (L,) representing
                3-state secondary structure prediction.
                0 stands for 'H', 1 for 'E' and 2 for 'C'.
            weights (:obj:`np.ndarray`): Matrix of restraint weights

        Returns:
            :obj:`gaussold.Model`: Gaussian model with
                restraints defined for each pair of residues.
        """
        # Cut predicted secondary structure
        # into contiguous segments
        L = len(ssp)
        segment_ids = [0]
        for i in range(1, L):
            if ssp[i] != ssp[i-1]:
                segment_ids.append(segment_ids[i-1] + 1)
            else:
                segment_ids.append(segment_ids[i-1])

        # Instantiate an empty model. Restraints have
        # to be defined for each pair of residues before
        # optimizing this model.
        model = Model(L)

        # Add backbone restraints:
        # Restraints based on average
        # C-alpha - C-alpha distance for residues
        # with a sequence separation of 1.
        for i in range(L - 1):
            model.add_restraint(i, i + 1, 3.82, 0.39, weight=1.)
        
        # Add restraints based on graph distances:
        # either contacts or non-contacts.
        for i in range(L):
            for j in range(max(0, i-self.sep)):
                if gds[i, j] > 0:
                    mu = gds[i, j] * 5.72
                    sigma = gds[i, j] * 1.34
                    model.add_restraint(i, j, mu, sigma, weight=weights[i, j])

        # Add restraints based on contacts in predicted
        # secondary structures
        for i in range(L):
            for j in range(max(0, i-self.sep)):
                if gds[i, j] == 1: # Contact
                    sep = np.abs(i - j)
                    if segment_ids[i] == segment_ids[j]:
                        if ssp[i] == 0: # Helix
                            if sep == 1:
                                model.add_restraint(i, j, 3.82, 0.35, weight=weights[i, j])
                            elif sep == 2:
                                model.add_restraint(i, j, 5.50, 0.52, weight=weights[i, j])
                            elif sep == 3:
                                model.add_restraint(i, j, 5.33, 0.93, weight=weights[i, j])
                            elif sep == 4:
                                model.add_restraint(i, j, 6.42, 1.04, weight=weights[i, j])
                        elif ssp[i] == 1: # Beta strand
                            if sep == 1:
                                model.add_restraint(i, j, 3.80, 0.28, weight=weights[i, j])
                            elif sep == 2:
                                model.add_restraint(i, j, 6.66, 0.30, weight=weights[i, j])
                    elif (ssp[i] == 0 and ssp[j] == 1) or (ssp[i] == 1 and ssp[j] == 0):
                        if sep >= 4: # Alpha/beta contact
                            model.add_restraint(i, j, 6.05, 0.95, weight=weights[i, j])
                    elif (ssp[i] == 0 and ssp[j] == 2) or (ssp[i] == 2 and ssp[j] == 0):
                        if sep >= 4: # Helix/coil contact
                            model.add_restraint(i, j, 6.60, 0.92, weight=weights[i, j])
                    elif (ssp[i] == 1 and ssp[j] == 2) or (ssp[i] == 2 and ssp[j] == 1):
                        if sep >= 4: # Beta/coil contact
                            model.add_restraint(i, j, 6.44, 1.00, weight=weights[i, j])
        return model

    @property
    def model(self):
        return self._model
    
    @model.setter
    def model(self, model):
        self._model = model

    @property
    def optimizer(self):
        return self._optimizer
    
    @optimizer.setter
    def optimizer(self, optimizer):
        self._optimizer = optimizer
