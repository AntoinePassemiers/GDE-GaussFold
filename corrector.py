# -*- coding: utf-8 -*-
# corrector.py: CA-CA and CB-CB distance deviation correction
# author : Antoine Passemiers

import numpy as np
from scipy.optimize import minimize_scalar


class DeviationCorrector:

    CA_CA_DISTANCE = 3.8

    def __init__(self, L):
        self.L = L
        self.n_polynomials = int(np.floor((self.L - 2.) / 2.))
        self.polynomials = [None] * self.n_polynomials
        self.inter_residue_distance = DeviationCorrector.CA_CA_DISTANCE

    def fit_transform(self, coords):
        coords = np.copy(coords)
        old_coords = np.copy(coords)
        k = 0
        while True:
            offset = (k % 2) + 1
            for i in range(offset, offset + self.n_polynomials * 2, 2):
                P1, P2, P3 = coords[i-1, :], coords[i, :], coords[i+1, :]
                AB = P2 - P1
                AC = P3 - P1
                pi = np.cross(AB, AC)
                d = -np.dot(P1, pi)
                nu = pi / np.linalg.norm(pi)
                zeta = -d / np.linalg.norm(pi)

                O = np.eye(3) + np.asarray(
                    [[0, 0, -nu[0]],
                     [0, 0, -nu[1]],
                     [nu[0], nu[1], 0]]) - \
                    (1. / (1. + nu[2])) * np.asarray(
                    [[nu[0]**2., nu[0]*nu[1], 0],
                     [nu[0]*nu[1], nu[1]**2., 0],
                     [0, 0, nu[0]**2. + nu[1]**2.]])
                Oinv = np.linalg.inv(O)

                project = lambda p: (np.dot(O, p) - np.asarray([0, 0, zeta]))[:2]
                project_back = lambda p: np.dot(Oinv, np.asarray([p[0], p[1], zeta]))

                XY = np.asarray([project(p) for p in coords[i-1:i+2, :]])
                coords[i, :] = project_back(self.correct_triple(XY))
            tau = np.linalg.norm(old_coords - coords)
            if tau < 1e-3 or k >= 12:
                break
            old_coords[:, :] = coords[:, :]
            k += 1

        return coords

    def correct_triple(self, coords):
        assert(coords.shape == (3, 2))
        coefs = np.polyfit(coords[:, 0], coords[:, 1], 2)
        poly = lambda x: coefs[0] * x ** 2. + coefs[1] * x + coefs[2]

        def objective(x_i):
            y_i = poly(x_i)
            left_distance = np.sqrt(
                (x_i - coords[0, 0]) ** 2. + (y_i - coords[0, 1]) ** 2.)
            left_cost = np.abs(left_distance - self.inter_residue_distance)
            right_distance = np.sqrt(
                (x_i - coords[2, 0]) ** 2. + (y_i - coords[2, 1]) ** 2.)
            right_cost = np.abs(right_distance - self.inter_residue_distance)
            return left_cost + right_cost

        bounds = (coords[0, 0], coords[2, 0])
        result = minimize_scalar(objective, bounds=bounds)
        x_hat = result.x
        y_hat = poly(x_hat)
        return np.asarray([x_hat, y_hat])
