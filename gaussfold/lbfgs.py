# -*- coding: utf-8 -*-
# lbfgs.py: Limited-memory BFGS algorithm
# author : Antoine Passemiers

import scipy.optimize


def lbfgs(initial_solution, model, verbose=True):
    """Run L-BFGS on an initial solution,
    with given objective function.

    Parameters:
        initial_solution (:obj:`np.ndarray`): Array of shape (L, 3)
            representing the initial solution, where L is the number
            of residues in the protein.
        model (:obj:`gaussfold.Model`): Gaussian model.
        verbose (bool): Whether to display messages in stdout.

    Returns:
        :obj:`np.ndarray`: Locally optimal solution.
    """
    L = initial_solution.shape[0]

    # Define objective function
    def obj(x):
        coords = x.reshape(L, 3)
        return -model.evaluate(coords)

    # Define gradient function
    def jac(x):
        coords = x.reshape(L, 3)
        grad = model.gradient(coords)
        return grad.flatten()

    # Define callback function
    def callback(x):
        if verbose:
            print(obj(x))

    # Solve the optimization problem
    x0 = initial_solution.flatten()
    res = scipy.optimize.minimize(
        obj, x0, jac=jac, method='L-BFGS-B', callback=callback)
    return res.x.reshape(L, 3)
