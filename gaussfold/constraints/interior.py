# -*- coding: utf-8 -*-
# interior.py
# author : Antoine Passemiers

from gaussfold.constraints.gaussian_constraint import GaussianConstraint


class Interior(GaussianConstraint):

    __MU__ = 5.0
    __SIGMA__ = 3.16

    def __init__(self, atom, **kwargs):
        center = GaussianConstraint.__CENTER_OF_MASS__
        GaussianConstraint.__init__(self, atom, center, **kwargs)

    def mu(self):
        return Interior.__MU__

    def sigma(self):
        return Interior.__SIGMA__
