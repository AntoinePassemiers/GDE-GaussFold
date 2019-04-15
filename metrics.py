# -*- coding: utf-8 -*-
# metrics.py
# author : Antoine Passemiers

import numpy as np
from scipy.optimize import minimize


def transform(coords, a, b, c, phi, psi, theta, swapx, swapy, swapz):
    Rx = np.asarray([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
    Ry = np.asarray([[np.cos(psi), 0, np.sin(psi)], [0, 1, 0], [-np.sin(psi), 0, np.cos(psi)]])
    Rz = np.asarray([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    A = np.dot(Rx, np.dot(Ry, Rz))
    coords_projected = np.dot(A, coords.T).T + np.asarray([a, b, c])
    if swapx > 0.5:
        coords_projected[:, 0] = -coords_projected[:, 0]
    if swapy > 0.5:
        coords_projected[:, 1] = -coords_projected[:, 1]
    if swapz > 0.5:
        coords_projected[:, 2] = -coords_projected[:, 2]
    return coords_projected


def __tm_score(coords_predicted, coords_target):
    L = len(coords_predicted)
    d_0 = 1.24 * np.cbrt(L - 15.) - 1.8
    distances = np.sqrt(((coords_predicted - coords_target) ** 2.).sum(axis=1))
    return (1. / (1. + (distances / d_0) ** 2.)).sum() / float(L)

def tm_score(coords_predicted, coords_target):
    coords_predicted = coords_predicted - np.mean(coords_predicted, axis=0)
    coords_target = coords_target - np.mean(coords_target, axis=0)

    def objective(x):
        coords_projected = transform(coords_predicted, *x)
        score = __tm_score(coords_projected, coords_target)
        return -score

    scores = list()
    d = 1000
    for swapx in [0, 1]:
        for swapy in [0, 1]:
            for swapz in [0, 1]:
                bounds = [(-d, d), (-d, d), (-d, d), (-2*np.pi, 2*np.pi), (-2*np.pi, 2*np.pi), (-2*np.pi, 2*np.pi), (0, 1), (0, 1), (0, 1)]
                a = 0
                b = 0
                c = 0
                phi = np.random.uniform(-2.*np.pi, 2.*np.pi)
                psi = np.random.uniform(-2.*np.pi, 2.*np.pi)
                theta = np.random.uniform(-2.*np.pi, 2.*np.pi)
                res = minimize(objective, [a, b, c, phi, psi, theta, swapx, swapy, swapz], bounds=bounds)
                score = -objective(res.x)
                scores.append(score)
    print(np.max(scores))
    return np.max(scores)
