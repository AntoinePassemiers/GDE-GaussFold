# -*- coding: utf-8 -*-
# exterior.py
# author : Antoine Passemiers

from gaussfold.constraints.gaussian_constraint import GaussianConstraint


class Exterior(GaussianConstraint):

    __MU__ = 12.5
    __SIGMA__ = 3.87

    def __init__(self, atom, **kwargs):
        center = GaussianConstraint.__CENTER_OF_MASS__
        GaussianConstraint.__init__(self, atom, center, **kwargs)

    def mu(self):
        return Exterior.__MU__

    def sigma(self):
        return Exterior.__SIGMA__
