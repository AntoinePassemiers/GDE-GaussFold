# -*- coding: utf-8 -*-
# core.py: GDE-GaussFold core algorithm
# author : Antoine Passemiers

from gaussfold.corrector import DeviationCorrector
from gaussfold.model import Model
from gaussfold.optimizer import Optimizer

import numpy as np
import random
from sklearn.manifold import MDS
import networkx as nx
from networkx.algorithms.shortest_paths.generic import shortest_path_length


class GaussFold:
    """GDE-GaussFold base class.

    Attributes:
        n_runs (int): Number of times to run Multi-Dimensional
            Scaling algorithm from scratch.
        max_n_iter (int): Maximum number of iterations of
            Multi-Dimensional Scaling algorithm.
        eps (float): Convergence threshold of Multi-Dimensional
            Scaling algorithm.
    """

    def __init__(self, n_runs=1, max_n_iter=300, eps=1e-3):
        """
        """
        self.n_runs = n_runs
        self.max_n_iter = max_n_iter
        self.eps = eps
        self.model = None
        self.optimizer = None

    def run(self, cmap, ssp):

        L = len(cmap)

        # Set diagonal to ones
        cmap = np.asarray(cmap)
        np.fill_diagonal(cmap, 1)
        # TODO: add more ones if not connected graph

        # Choose threshold such that exactly 4.5*L contacts
        # are obtained
        proba = cmap[np.triu_indices(L, -6)]
        np.sort(proba)
        threshold = proba[-int(np.round(4.5 * L))]

        # Create grah from adjacency matrix
        cmap = np.nan_to_num(cmap)
        A = (cmap > threshold)
        G = nx.from_numpy_array(A, parallel_edges=False)
        # TODO: cost of edges in fuzzy case?

        # Compute graph distance map
        path_lengths = shortest_path_length(G)
        gds = np.zeros((L, L), dtype=np.int)
        for i, i_lengths in path_lengths:
            for j in i_lengths.keys():
                gds[i, j] = i_lengths[j]
                gds[j, i] = gds[i, j]

        # Apply theoretical linear correspondence between graph
        # distance and Angstroms distance based on statistical
        # observations on euclidean distances found for a graph
        # distance of 1.
        # The empirical value of 4.846 could be used as well.
        distances = gds * 5.72

        print('MDS')

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

        print('Apply deviation correction')

        # Apply correction on pairs of adjacent residues
        # based on known C_alpha-C_alpha (or C_beta-C_beta) distance
        corrector = DeviationCorrector(L)
        X_transformed = corrector.fit_transform(X_transformed)

        if not isinstance(self.model, Model):
            self.model = self.create_model(gds, ssp)

        if not isinstance(self.optimizer, Optimizer):
            self.optimizer = Optimizer(X_transformed, self.model.evaluate)
        best_coords = self.optimizer.run(n_iter=100000)
        return best_coords

    def create_model(self, gds, ssp):
        L = len(ssp)
        segment_ids = [0]
        for i in range(1, L):
            if ssp[i] != ssp[i-1]:
                segment_ids.append(segment_ids[i-1] + 1)
            else:
                segment_ids.append(segment_ids[i-1])

        model = Model(L)
        
        for i in range(L):
            for j in range(max(0, i-1)):
                mu = gds[i, j] * 5.72
                sigma = gds[i, j] * 1.34
                model.add_restraint(i, j, mu, sigma)

        for i in range(L - 1):
            model.add_restraint(i, i + 1, 3.82, 0.39)

        for i in range(L):
            for j in range(i):
                sep = np.abs(i - j)
                if segment_ids[i] == segment_ids[j]:
                    if gds[i, j] == 1: # Intra-segment contact
                        if ssp[i] == 0: # Helix
                            if sep == 1:
                                model.add_restraint(i, j, 3.82, 0.35)
                            elif sep == 2:
                                model.add_restraint(i, j, 5.50, 0.52)
                            elif sep == 3:
                                model.add_restraint(i, j, 5.33, 0.93)
                            elif sep == 4:
                                model.add_restraint(i, j, 6.42, 1.04)
                        elif ssp[i] == 1: # Beta strand
                            if sep == 1:
                                model.add_restraint(i, j, 3.80, 0.28)
                            elif sep == 2:
                                model.add_restraint(i, j, 6.66, 0.30)
                elif (ssp[i] == 0 and ssp[j] == 1) or (ssp[i] == 1 and ssp[j] == 0):
                    if gds[i, j] == 1 and sep >= 4: # Alpha/beta contact
                        model.add_restraint(i, j, 6.05, 0.95)
                elif (ssp[i] == 0 and ssp[j] == 2) or (ssp[i] == 2 and ssp[j] == 0):
                    if gds[i, j] == 1 and sep >= 4: # Helix/coil contact
                        model.add_restraint(i, j, 6.60, 0.92)
                elif (ssp[i] == 1 and ssp[j] == 2) or (ssp[i] == 2 and ssp[j] == 1):
                    if gds[i, j] == 1 and sep >= 4: # Beta/coil contact
                        model.add_restraint(i, j, 6.44, 1.00)
        return model
