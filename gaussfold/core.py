# -*- coding: utf-8 -*-
# core.py: GDE-GaussFold core algorithm
# author : Antoine Passemiers

from gaussfold.aa import Glycine, Cysteine
from gaussfold.chain.chain import Chain
from gaussfold.constraints import *
from gaussfold.corrector import DeviationCorrector
from gaussfold.graph import Graph
from gaussfold.model.amino_acid_model import AminoAcidModel
from gaussfold.model.all_atom_model import AllAtomModel
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

    def __init__(self, sep=1, n_runs=1, max_n_iter=300, eps=1e-3):
        self.sep = sep
        self.n_runs = n_runs
        self.max_n_iter = max_n_iter
        self.eps = eps
        self._model = None
        self._optimizer = None

    def run(self, cmap, ssp, acc, seq, verbose=True):
        """Runs GDE-GaussFold algorithm.

        Parameters:
            cmap (:obj:`np.ndarray`): Array of shape (L, L) representing
                predicted contact probabilities, where L is the number of
                residues in the sequence.
            ssp (:obj:`np.ndarray`): Array of shape (L,) representing
                3-state secondary structure prediction.
                0 stands for 'H', 1 for 'E' and 2 for 'C'.
            acc (:obj:`np.ndarray`): Array of shape (L,) representing
                3-state solvent accessibility prediction.
                0 stands for 'buried', 1 for 'medium' and 2 for 'exposed'.
            seq (str): Protein primary structure.
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
        n_top = int(np.round(2.5 * L))
        threshold = proba[-n_top]

        A = (cmap > threshold)

        G = Graph(A)
        gds = G.distances()
        missing = np.zeros(L, dtype=np.bool)
        
        """
        if not G.is_connected():
            missing = (gds == 0).all(axis=0)
            missing_idx = np.where(missing)[0]
            connected_idx = np.random.choice(
                    np.arange(L), size=len(missing_idx), replace=False)
            A[missing_idx, connected_idx] = 1
            A[connected_idx, missing_idx] = 1
            G = Graph(A)
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

        # Graph distances above 14 are statistically impossible
        gds = np.minimum(gds, 14)

        """
        while (gds == 14).any():
            values = np.copy(cmap)
            values[gds < 14] = 0
            i, j = np.unravel_index(values.argmax(), values.shape)
            A[i, j] = 1
            G = Graph(A)
            gds = G.distances()
            gds = np.minimum(gds, 14)
        """

        import matplotlib.pyplot as plt
        plt.imshow(gds)
        plt.show()
        #import sys; sys.exit(0)

        # Apply theoretical linear correspondence between graph
        # distance and Angstroms distance based on statistical
        # observations on euclidean distances found for a graph
        # distance of 1.
        # The empirical value of 4.846 could be used as well.
        distances = gds * 5.72

        if verbose:
            print('Apply Multi-Dimensional Scaling algorithm')

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
        except ValueError:
            if verbose:
                print('[Warning] Invalid value encountered in deviation corrector')
        initial_coords = X_transformed

        # Create Gaussian model if not set by the user
        if not isinstance(self._model, AminoAcidModel):
            if verbose:
                print('Model not set by user. Creating model from scratch...')
            chain = Chain.from_string(seq, c='CA')
            self._model = self.create_model(chain, cmap, gds, ssp, acc, weights)
        for i in range(len(chain)):
            chain[i].ref().set_coords(*initial_coords[i])

        # Create optimizer if not set by the user.
        # Use default hyper-parameters.
        if not isinstance(self._optimizer, Optimizer):
            if verbose:
                print('Optimizer not set by user. Using default parameters.')
            self._optimizer = Optimizer()
        
        # Run optimizer on the Gaussian model
        self._optimizer.run(self._model, verbose=verbose)
        best_coords = np.empty((len(chain), 3), dtype=np.float)
        for i in range(len(chain)):
            best_coords[i, :] = chain[i].ref().get_coords()
        return best_coords

    def create_model(self, chain, cmap, gds, ssp, acc, weights):
        """Creates a Gaussian model for the protein.

        Parameters:
            gds (:obj:`np.ndarray`): Matrix of graph distance
                between each pair of residues in the protein.
            ssp (:obj:`np.ndarray`): Array of shape (L,) representing
                3-state secondary structure prediction.
                0 stands for 'H', 1 for 'E' and 2 for 'C'.
            weights (:obj:`np.ndarray`): Matrix of restraint weights

        Returns:
            :obj:`gaussold.model.AminoAcidModel`: Gaussian model with
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
        model = AminoAcidModel(L)

        # Repulsion constraints
        for i in range(L):
            for j in range(max(0, i-self.sep)):
                model.add_constraint(Repulsion(chain[i].ref(), chain[j].ref()))

        # Surface accessibility
        for i in range(L):
            if acc[i] == 0:
                model.add_constraint(Interior(chain[i].ref()))
            #elif acc[i] == 2:
            #    model.add_constraint(Exterior(chain[i].ref()))

        # Disulfide bonds
        for constraint in self.make_disulfide_bonds(cmap, chain):
            model.add_constraint(constraint)

        # Add backbone restraints:
        # Restraints based on average
        # C-alpha - C-alpha distance for residues
        # with a sequence separation of 1.
        for i in range(L - 1):
            model.add_constraint(Adjacent(chain[i].ref(), chain[i + 1].ref(), 1))

        # Next adjacent restraints
        for i in range(L - 2):
            model.add_constraint(Adjacent(chain[i].ref(), chain[i + 2].ref(), 2))
        for i in range(L - 3):
            model.add_constraint(Adjacent(chain[i].ref(), chain[i + 3].ref(), 3))

        # Regular contacts
        for i in range(L):
            for j in range(max(0, i-self.sep)):
                if gds[i, j] == 1:
                    model.add_constraint(DistanceRestraint(chain[i].ref(), chain[j].ref(), 3.82, 0.35))

        # Add restraints based on contacts in predicted
        # secondary structures
        for i in range(L):
            for j in range(0, i):
                sep = np.abs(i - j)
                if segment_ids[i] == segment_ids[j]:
                    if ssp[i] == 0: # Helix
                        if sep == 2:
                            model.add_constraint(DistanceRestraint(
                                    chain[i].ref(), chain[j].ref(), 5.48, 0.14, weight=weights[i, j]))
                        elif sep == 3:
                            model.add_constraint(DistanceRestraint(
                                    chain[i].ref(), chain[j].ref(), 5.20, 0.14, weight=weights[i, j]))
                        elif sep == 4:
                            model.add_constraint(DistanceRestraint(
                                    chain[i].ref(), chain[j].ref(), 6.28, 0.26, weight=weights[i, j]))
                        elif sep == 5:
                            model.add_constraint(DistanceRestraint(
                                    chain[i].ref(), chain[j].ref(), 8.75, 0.26))
                    elif ssp[i] == 1: # Beta strand
                        if sep == 2:
                            model.add_constraint(DistanceRestraint(
                                    chain[i].ref(), chain[j].ref(), 6.74, 0.28, weight=weights[i, j]))
                        elif sep == 3:
                            model.add_constraint(DistanceRestraint(
                                    chain[i].ref(), chain[j].ref(), 10.10, 0.32, weight=weights[i, j]))
                        elif sep == 4:
                            model.add_constraint(DistanceRestraint(
                                    chain[i].ref(), chain[j].ref(), 13.30, 1.41, weight=weights[i, j]))

        for i in range(L):
            for j in range(i):
                if gds[i, j] == 1: # Contact
                    if segment_ids[i] == segment_ids[j]:
                        sep = np.abs(i - j)
                        if ssp[i] == 1 and ssp[j] == 1:
                            if sep >= 4:
                                model.add_constraint(DistanceRestraint(
                                        chain[i].ref(), chain[j].ref(), 4.54, 0.32))
                        elif (ssp[i] == 0 and ssp[j] == 1) or (ssp[i] == 1 and ssp[j] == 0):
                            if sep >= 4: # Alpha/beta contact
                                model.add_constraint(DistanceRestraint(
                                        chain[i].ref(), chain[j].ref(), 6.05, 0.95, weight=weights[i, j]))
                        elif (ssp[i] == 0 and ssp[j] == 2) or (ssp[i] == 2 and ssp[j] == 0):
                            if sep >= 4: # Helix/coil contact
                                model.add_constraint(DistanceRestraint(
                                        chain[i].ref(), chain[j].ref(), 6.60, 0.92, weight=weights[i, j]))
                        elif (ssp[i] == 1 and ssp[j] == 2) or (ssp[i] == 2 and ssp[j] == 1):
                            if sep >= 4: # Beta/coil contact
                                model.add_constraint(DistanceRestraint(
                                        chain[i].ref(), chain[j].ref(), 6.44, 1.00, weight=weights[i, j]))
        return model.initialize()

    def make_disulfide_bonds(self, cmap, chain):
        cysteine_ids = [i for i, amino_acid in enumerate(chain) if isinstance(amino_acid, Cysteine)]
        n_cysteines = len(cysteine_ids)
        id_pairs = list()
        for i in cysteine_ids:
            for j in cysteine_ids:
                if i < j:
                    id_pairs.append((i, j, cmap[i, j]))

        paired_ids = list()
        constraints = list()
        while n_cysteines >= 2:
            idx = np.argmax([score for _, _, score in id_pairs])
            i, j, _ = id_pairs[idx]
            paired_ids += [i, j]
            print('[Model] Disulfide bond created between residues %i and %i' % (i, j))
            constraints.append(DisulfideBond(chain[i].ref(), chain[j].ref()))

            new_id_pairs = list()
            for i, j, score in id_pairs:
                if i not in paired_ids and j not in paired_ids:
                    new_id_pairs.append((i, j, score))
            id_pairs = new_id_pairs

            n_cysteines -= 2
        return constraints

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
