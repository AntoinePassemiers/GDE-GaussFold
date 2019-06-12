# -*- coding: utf-8 -*-
# corrector.py: CA-CA distance deviation correction
# author : Antoine Passemiers

import numpy as np
from scipy.optimize import minimize_scalar


class DeviationCorrector:
    """Applies corrections to 3D coordinates of residues that do
    not have a distance of 3.8 Angstroms with adjacent residues.
    This distance corresponds to the average C-alpha - C-alpha
    distance.

    Attributes:
        L (int): Number of residues in the protein.
        n_polynomials (int): Number of polynomial interpolations
            to be created.
        polynomials (list): List of polynomial interpolations.
            Each interpolation is a parabola constructed
            from 3 adjacent residues.
        inter_residue
    """

    CA_CA_DISTANCE = 3.82

    def __init__(self, L):
        """Constructs a deviation corrector.

        Parameters:
            L (int): Number of residues in the protein.
        """
        self.L = L
        self.n_polynomials = int(np.floor((self.L - 2.) / 2.))
        self.polynomials = [None] * self.n_polynomials

    def fit_transform(self, coords):
        """Applies corrections to coordinates of adjacent residues,
        based on average C-alpha - C-alpha distance.

        Parameters:
            coords (np.ndarray): Array of shape (L, 3) representing
                the initial coordinates.

        Returns;
            np.ndarray: Array of shape (L, 3) representing the
                refined coordinates.
        """
        coords = np.copy(coords)
        old_coords = np.copy(coords)
        k = 0
        while True:
            offset = (k % 2) + 1
            for i in range(offset, offset + self.n_polynomials * 2, 2):

                # Find the plane defined by the 3 points.
                # Points are assumed to be not aligned.
                P1, P2, P3 = coords[i-1, :], coords[i, :], coords[i+1, :]
                AB = P2 - P1
                AC = P3 - P1
                pi = np.cross(AB, AC)
                d = -np.dot(P1, pi)
                nu = pi / np.linalg.norm(pi)
                zeta = -d / np.linalg.norm(pi)

                # Compute rotation and inverse rotation matrices
                O = np.eye(3) + np.asarray(
                    [[0, 0, -nu[0]],
                     [0, 0, -nu[1]],
                     [nu[0], nu[1], 0]]) - \
                    (1. / (1. + nu[2])) * np.asarray(
                    [[nu[0]**2., nu[0]*nu[1], 0],
                     [nu[0]*nu[1], nu[1]**2., 0],
                     [0, 0, nu[0]**2. + nu[1]**2.]])
                Oinv = np.linalg.inv(O)

                # Define projection and inverse projection
                project = lambda p: (np.dot(O, p) - np.asarray([0, 0, zeta]))[:2]
                project_back = lambda p: np.dot(Oinv, np.asarray([p[0], p[1], zeta]))

                # Project points onto a plane, apply corrections and project 
                # them back in 3D
                XY = np.asarray([project(p) for p in coords[i-1:i+2, :]])
                coords[i, :] = project_back(self.correct_triple(XY))

            # Check for convergence
            tau = np.linalg.norm(old_coords - coords)
            if tau < 1e-3 or k >= 12:
                break
            old_coords[:, :] = coords[:, :]
            k += 1

        return coords

    def correct_triple(self, coords):
        """Apply correction on second point such that
        the C-alpha - C-alpha distances are minimized
        with neighbouring points.

        Parameters:
            coords (:obj:`np.ndarray`): Array of shape (3, 2)
                containing the 2D coordinates of the 3 points.

        Returns:
            :obj:`np.ndarray`: Array of shape (2,) representing
                the new 2D coordinates of the second point.
        """
        assert(coords.shape == (3, 2))
        
        # Create a parabolic interpolation
        coefs = np.polyfit(coords[:, 0], coords[:, 1], 2)
        poly = lambda x: coefs[0] * x ** 2. + coefs[1] * x + coefs[2]

        # Define the objective function as the sum of the deviation of
        # the distance between first and second points from the C-alpha -
        # C-alpha distance, and the deviation of the distance between
        # second and third points from the C-alpha - C-alpha distance.
        def objective(x_i):
            """
            Parameters:
                x_i (float): Position of the second point on the X-axis.
                    Position on Y-axis is then determined by the parabolic
                    interpolation.
            """
            y_i = poly(x_i)
            left_distance = np.sqrt(
                (x_i - coords[0, 0]) ** 2. + (y_i - coords[0, 1]) ** 2.)
            left_cost = np.abs(left_distance - DeviationCorrector.CA_CA_DISTANCE)
            right_distance = np.sqrt(
                (x_i - coords[2, 0]) ** 2. + (y_i - coords[2, 1]) ** 2.)
            right_cost = np.abs(right_distance - DeviationCorrector.CA_CA_DISTANCE)
            return left_cost + right_cost

        # Optimize the objective function and find a new position
        # for the second point
        bounds = (coords[0, 0], coords[2, 0])
        result = minimize_scalar(objective, bounds=bounds)
        x_hat = result.x
        y_hat = poly(x_hat)
        return np.asarray([x_hat, y_hat])
