# -*- coding: utf-8 -*-
# bond.py: Atom bonds
# author : Antoine Passemiers


class Bond:

    def __init__(self, atom1, atom2, order=1):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        self.distance = self.compute_distance(atom1, atom2, order=order)

        self.atom1.add_bonded_atom(atom2)
        self.atom2.add_bonded_atom(atom1)

    def compute_distance(self, atom1, atom2, order=1):
        a1, a2 = atom1.element, atom2.element
        elements = set([a1, a2])
        if elements == set(['N', 'C']):
            distance = 1.45
        elif elements == set(['C', 'C']):
            distance = 1.49
        elif elements == set(['C', 'O']):
            if order == 1:
                distance = 1.33
            else:
                distance = 1.22
        # TODO: peptide bonds, etc.
        return distance
