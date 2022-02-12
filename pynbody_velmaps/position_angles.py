import numpy as np
from pafit import fit_kinematic_pa


def infer_coordinates(velmap):
    x = np.arange(velmap.shape[1])
    y = np.arange(velmap.shape[0])

    x -= len(x) // 2
    y -= len(y) // 2
    x, y = np.meshgrid(x, y)
    return x, y[::-1,:]


def calc_pa(velmap, **kwargs):
    x, y = infer_coordinates(velmap)

    pa = fit_kinematic_pa.fit_kinematic_pa(x, y, velmap, plot=False, **kwargs)
    return np.asarray(pa)