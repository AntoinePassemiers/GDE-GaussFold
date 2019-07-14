# -*- coding: utf-8 -*-
# gaussian_constraint.py
# author : Antoine Passemiers

from gaussfold.constraints.gaussian_constraint import GaussianConstraint


class Repulsion(GaussianConstraint):

    __MU__ = 20.0
    __SIGMA__ = 10.95

    def __init__(self, atom_a, atom_b, **kwargs):
        GaussianConstraint.__init__(self, atom_a, atom_b, **kwargs)

    def mu(self):
        return Repulsion.__MU__

    def sigma(self):
        return Repulsion.__SIGMA__
