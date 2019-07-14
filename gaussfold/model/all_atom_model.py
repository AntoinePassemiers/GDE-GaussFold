# -*- coding: utf-8 -*-
# all_atom_model.py
# author : Antoine Passemiers

import numpy as np
import scipy.spatial


class AllAtomModel:

    def __init__(self, chain):
        self.chain = chain
        self.n_atoms = len(self.chain.atoms())
        self.n_residues = len(self.chain)
        self.triu_indices = np.triu_indices(self.n_atoms, k=-1)
        self.distances = np.empty((self.n_atoms, self.n_atoms), dtype=np.float)

        # Parameters for modelling distance restraints
        self.lb = np.full((self.n_atoms, self.n_atoms), -np.inf, dtype=np.float)
        self.ub = np.full((self.n_atoms, self.n_atoms), np.inf, dtype=np.float)

        # Parameters for modelling angle restraints
        self.phi_lb = np.full(self.n_residues, -np.inf, dtype=np.float)
        self.phi_ub = np.full(self.n_residues, np.inf, dtype=np.float)
        self.psi_lb = np.full(self.n_residues, -np.inf, dtype=np.float)
        self.psi_ub = np.full(self.n_residues, np.inf, dtype=np.float)

    def add_distance_restraint(self, i, j, lb=None, ub=None):
        if lb is not None and lb > self.lb[i, j]:
            self.lb[i, j] = self.lb[j, i] = lb
        if ub is not None and ub < self.ub[i, j]:
            self.ub[i, j] = self.ub[j, i] = ub

    def add_angle_restraint(self, i, phi_lb=None, phi_ub=None, psi_lb=None, psi_ub=None):
        if phi_lb is not None and phi_lb > self.phi_lb[i]:
            self.phi_lb[i] = phi_lb
        if phi_ub is not None and phi_ub < self.phi_ub[i]:
            self.phi_ub[i] = phi_ub
        if psi_lb is not None and psi_lb > self.psi_lb[i]:
            self.psi_lb[i] = psi_lb
        if psi_ub is not None and psi_ub < self.psi_ub[i]:
            self.psi_ub[i] = psi_ub

    def evaluate(self, coords):
        self.chain.set_atoms_coords(coords)
        
        scipy.spatial.distance.cdist(coords, coords, metric='euclidean', out=self.distances)
        distances = self.distances[self.triu_indices]
        lb = self.lb[self.triu_indices]
        ub = self.ub[self.triu_indices]
        cost = (distances < lb).sum()
        cost += (distances > ub).sum()

        phi, psi = self.chain.dihedral_angles()
        cost += (phi < self.phi_lb).sum()
        cost += (phi > self.phi_ub).sum()
        cost += (psi < self.psi_lb).sum()
        cost += (psi > self.psi_ub).sum()

        return -cost
