# -*- coding: utf-8 -*-
# ramachandran.py: Dihedral angles and Ramachandran diagrams
# author : Antoine Passemiers

import numpy as np


def dihedral_angle(p1, p2, p3, p4):
    
    # Compute vectors from coordinates
    v1 = p2 - p1
    v1 /= np.linalg.norm(v1)
    v2 = p3 - p2
    v2 /= np.linalg.norm(v2)
    v3 = p4 - p3
    v3 /= np.linalg.norm(v3)

    # Normal vector of the plane containing v1 and v2
    nv1 = np.cross(v1, v2)

    # Normal vector of the plane containing v2 and v3
    nv2 = np.cross(v2, v3)

    # Vector that is orthogonal to nv1 and v2
    o1 = np.cross(nv1, v2)

    # Coordinates of nv2 in the frame (nv1, v2, o1)
    x = np.dot(nv1, nv2)
    y = np.dot(o1, nv2)

    return np.atan2(y, x)
