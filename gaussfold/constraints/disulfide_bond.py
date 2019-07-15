# -*- coding: utf-8 -*-
# disulfide_bond.py
# author : Antoine Passemiers

from gaussfold.constraints.gaussian_constraint import GaussianConstraint


class DisulfideBond(GaussianConstraint):

    __MU__ = 5.50
    __SIGMA__ = 0.45

    def __init__(self, atom_a, atom_b, **kwargs):
        GaussianConstraint.__init__(self, atom_a, atom_b, **kwargs)

    def mu(self):
        return DisulfideBond.__MU__

    def sigma(self):
        return DisulfideBond.__SIGMA__
