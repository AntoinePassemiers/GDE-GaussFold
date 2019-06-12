# -*- coding: utf-8 -*-
# bond.py: Atom bonds
# author : Antoine Passemiers


class Bond:

    def __init__(self, atom1, atom2, order=1):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        self.distance, self.sigma = self.compute_distance(
                atom1, atom2, order=order)

        self.atom1.add_bonded_atom(atom2)
        self.atom2.add_bonded_atom(atom1)

    def compute_distance(self, atom1, atom2, order=1):
        a1, a2 = atom1.element, atom2.element
        if ord(a2) < ord(a1):
            a1, a2 = a2, a1
            atom1, atom2 = atom2, atom1
        elements = [a1, a2]
        distance, sigma = 1.3, 0.01
        if elements == ['C', 'N']:
            if atom1.name == 'CA' and atom2.name == 'N':
                distance, sigma = 1.32, 0.01
            else:
                distance, sigma = 1.46, 0.01
        elif elements == ['C', 'C']:
            distance, sigma = 1.49, 0.01
        elif elements == ['C', 'O']:
            if order == 1:
                distance, sigma = 1.33, 0.01
            else:
                distance, sigma = 1.22, 0.01
        # TODO: peptide bonds, etc.
        return distance, sigma
