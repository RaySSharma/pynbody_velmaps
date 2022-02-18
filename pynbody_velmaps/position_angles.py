import numpy as np
from pafit import fit_kinematic_pa


def infer_coordinates(vel_map):
    x = np.linspace(-vel_map.image_width_kpc/2, vel_map.image_width_kpc/2, vel_map.npixels)
    x, y = np.meshgrid(x, x)
    return x, y[::-1,:]


def calc_pa(vel_map, **kwargs):
    x, y = infer_coordinates(vel_map)
    pa = fit_kinematic_pa.fit_kinematic_pa(x, y, vel_map.vel_map, plot=False, **kwargs)
    return pa