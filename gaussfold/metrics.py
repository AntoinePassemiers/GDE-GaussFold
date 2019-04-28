# -*- coding: utf-8 -*-
# metrics.py: Evaluation metrics for 3D models
# author : Antoine Passemiers

import numpy as np
from scipy.optimize import minimize


OPTIMIZATION_BOUNDS = [
    (  -1000.,   1000.), # Translation on X-axis (in angstroms)
    (  -1000.,   1000.), # Translation on Y-axis (in angstroms)
    (  -1000.,   1000.), # Translation on Z-axis (in angstroms)
    (-2*np.pi, 2*np.pi), # Rotation around X-axis (in radians)
    (-2*np.pi, 2*np.pi), # Rotation around Y-axis (in radians)
    (-2*np.pi, 2*np.pi), # Rotation around Z-axis (in radians)
    (       0,       1), # Symmetry with respect to X (boolean)
    (       0,       1), # Symmetry with respect to Y (boolean)
    (       0,       1)  # Symmetry with respect to Z (boolean)
]


def transform(coords, a, b, c, phi, psi, theta, swapx, swapy, swapz):
    """Applies symmetry, rotation and translation to 3D points.

    Parameters: 
        coords (:obj:`np.ndarray`): Array of shape (n_points, 3).
        a (float): Translation on the X-axis.
        b (float): Translation on the Y-axis.
        c (float): Translation on the Z-axis.
        phi (float): Rotation angle (in radians) around X-axis.
        psi (float): Rotation angle (in radians) around Y-axis.
        theta (float): Rotation angle (in radians) around Z-axis.
        swapx (bool): Whether to apply symmetry w.r.t. X-axis.
        swapy (bool): Whether to apply symmetry w.r.t. Y-axis.
        swapz (bool): Whether to apply symmetry w.r.t. Z-axis.

    Returns:
        :obj:`np.ndarray`: Array of shape (n_points, 3) representing
            the transformed points.
    """
    # Rotation matrix around X-axis
    Rx = np.asarray(
            [[1,           0,            0],
             [0, np.cos(phi), -np.sin(phi)],
             [0, np.sin(phi),  np.cos(phi)]])

    # Rotation matrix around Y-axis
    Ry = np.asarray(
            [[ np.cos(psi), 0, np.sin(psi)],
             [           0, 1,           0],
             [-np.sin(psi), 0, np.cos(psi)]])

    # Rotation matrix around Z-axis
    Rz = np.asarray(
            [[np.cos(theta), -np.sin(theta), 0],
             [np.sin(theta),  np.cos(theta), 0],
             [            0,              0, 1]])

    # Complete rotation matrix
    A = np.dot(Rx, np.dot(Ry, Rz))

    # Apply projection
    coords_projected = np.dot(A, coords.T).T + np.asarray([a, b, c])
    if swapx > 0.5:
        coords_projected[:, 0] = -coords_projected[:, 0]
    if swapy > 0.5:
        coords_projected[:, 1] = -coords_projected[:, 1]
    if swapz > 0.5:
        coords_projected[:, 2] = -coords_projected[:, 2]
    return coords_projected


def projection_invariant(method='Maximize'):
    """Decorator for making evaluation metrics
    symmetry, rotation and translation-invariant in 3D space.

    Parameters:
        maximize (bool, optional): Whether to maximize the
            evaluation metric.

    Returns:
        callable: The same function but transformation-invariant.
    """

    assert(method in ['Maximize', 'Minimize'])
    maximize = (method == 'Maximize')

    def decorator(func):
        """
        Parameters:
            func (callable): Function taking as first argument an
                array of shape (n_points, 3) representing predicted
                points and as second argument an array of shape
                (n_points, 3) representing ground-truth points.
        """

        def new_func(coords_predicted, coords_target, return_coords=False):
            
            # Remove invalid indices
            valid_indices = np.asarray(
                [i for i, x in enumerate(coords_target) if x is not None])
            coords_predicted = np.asarray(coords_predicted)[valid_indices]
            coords_target = np.asarray([coords_target[i] for i in valid_indices])

            coords_predicted = coords_predicted - np.mean(coords_predicted, axis=0)
            coords_target = coords_target - np.mean(coords_target, axis=0)

            def objective(x):
                coords_projected = transform(coords_predicted, *x)
                score = func(coords_projected, coords_target)
                return score if not maximize else -score

            scores = list()
            sols = list()
            for swapx in [0, 1]:
                for swapy in [0, 1]:
                    for swapz in [0, 1]:
                        
                        # Start at the origin
                        a = b = c = 0
                        
                        # Start with random angles
                        phi = np.random.uniform(-2.*np.pi, 2.*np.pi)
                        psi = np.random.uniform(-2.*np.pi, 2.*np.pi)
                        theta = np.random.uniform(-2.*np.pi, 2.*np.pi)
                        
                        # Optimize
                        res = minimize(
                                objective,
                                [a, b, c, phi, psi, theta, swapx, swapy, swapz],
                                bounds=OPTIMIZATION_BOUNDS)
                        sols.append(res.x)
                        score = objective(res.x)
                        scores.append(score)
            scores = np.asarray(scores)
            best_score = np.max(-scores) if maximize else np.min(scores)
            if return_coords:
                best_idx = np.argmax(-scores) if maximize else np.argmin(scores)
                x = sols[best_idx]
                coords_projected = transform(coords_predicted, *x)
                return best_score, coords_projected
            else:
                return best_score
        new_func.__name__ = func.__name__
        return new_func
    return decorator


@projection_invariant('Maximize')
def tm_score(coords_predicted, coords_target):
    """Computes template-modelling score.

    Parameters:
        coords_predicted (:obj:`np.ndarray`): Array of shape
            (n_points, 3) representing predicted points.
        coords_target (:obj:`np.ndarray`): Array of shape
            (n_points, 3) representing ground-truth points.

    Returns:
        float: Template-modelling score.
    """
    L = len(coords_predicted)
    d_0 = 1.24 * np.cbrt(L - 15.) - 1.8
    distances = np.sqrt(((coords_predicted - coords_target) ** 2.).sum(axis=1))
    return (1. / (1. + (distances / d_0) ** 2.)).sum() / float(L)


@projection_invariant('Minimize')
def rmsd(coords_predicted, coords_target):
    """Computes root-mean-square deviation.

    Parameters:
        coords_predicted (:obj:`np.ndarray`): Array of shape
            (n_points, 3) representing predicted points.
        coords_target (:obj:`np.ndarray`): Array of shape
            (n_points, 3) representing ground-truth points.

    Returns:
        float: Root-mean-square deviation.
    """
    return np.sqrt(((coords_predicted - coords_target) ** 2.).sum(axis=1).mean())
