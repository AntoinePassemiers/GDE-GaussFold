# -*- coding: utf-8 -*-
# distance_restraint.py
# author : Antoine Passemiers

from gaussfold.constraints.gaussian_constraint import GaussianConstraint


class DistanceRestraint(GaussianConstraint):

    def __init__(self, atom_a, atom_b, mu, sigma, **kwargs):
        GaussianConstraint.__init__(self, atom_a, atom_b, **kwargs)
        self._mu = mu
        self._sigma = sigma

    def mu(self):
        return self._mu

    def sigma(self):
        return self._sigma
