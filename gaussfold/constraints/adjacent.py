# -*- coding: utf-8 -*-
# adjacent.py
# author : Antoine Passemiers

from gaussfold.constraints.gaussian_constraint import GaussianConstraint


class Adjacent(GaussianConstraint):

    def __init__(self, atom_a, atom_b, sep, **kwargs):
        GaussianConstraint.__init__(self, atom_a, atom_b, **kwargs)
        self._sep = sep
        if sep == 1:
            self._mu, self._sigma = 3.81, 0.39
        elif sep == 2:
            self._mu, self._sigma = 5.20, 0.55
        elif sep == 3:
            self._mu, self._sigma = 7.00, 0.71
        else:
            pass # TODO: exception

    def mu(self):
        return self._mu

    def sigma(self):
        return self._sigma
