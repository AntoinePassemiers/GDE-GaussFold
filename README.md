# Graph Distance Evolutionary Gaussian Folding


GDE-GaussFold uses Gaussian restraints to model distances between pairs
of residues and creates the three-dimensional representation that
maximizes log-likelihood. Optimization is done heuristically with
a non-stationary genetic algorithm.

```python
from gaussfold import GaussFold

# cmap is an array of shape (L, L) of
# predicted probabilities.
cmap = ...

# ssp is an array of shape (L,) representing
# 3-state secondary structure prediction.
# 0 stands for 'H', 1 for 'E' and 2 for 'C'.
ssp = ...

coords_predicted = GaussFold().run(cmap, ssp)
```

The library provides projection-invariant evaluation metrics
for the predicted 3D models:

```python
from gaussfold import tm_score, rmsd

# coords_target is an array of shape (L, 3)
# representing the 3D coordinates of the protein
# native structure. If at least one residue is unknown,
# then coords_target is a list of triples where
# missing residues are replaced by None.
coords_target = ...

print(tm_score(coords_predicted, coords_target))
print(rmsd(coords_predicted, coords_target))
```

A complete example is available in example/ folder.


### Installation

GDE-GaussFold can be installed with the following command:

```
python setup.py install
```

### Dependencies

* Numpy
* Scipy
* NetworkX
* Scikit-learn

Optional dependencies (but required in order to use the PDBParser):

* minineedle (https://github.com/scastlara/minineedle)
* miniseq (https://github.com/scastlara/miniseq)