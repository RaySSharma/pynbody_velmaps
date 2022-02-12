import argparse
import pathlib

import numpy as np
import pynbody
import matplotlib.pyplot as plt
from astropy.cosmology import Planck13

from pynbody_velmaps import generate, plot, load


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("simulation", help="Path to simulation.", type=str)
    parser.add_argument(
        "redshift", help="Redshift at which to generate image.", type=float
    )
    parser.add_argument(
        "image_width",
        help="Physical size of output image, in kpc.",
        type=float,
    )
    return parser.parse_args()


def calc_half_mass_radius(gal, ndim=2):
    bins = np.arange(0.1, gal["r"].max(), 0.1)
    pro = pynbody.analysis.profile.Profile(gal, ndim=ndim, bins=bins)
    half_mass = 0.5 * pro["mass_enc"].max()
    half_mass_r = pro["rbins"][abs(pro["mass_enc"] - half_mass).argmin()]
    return half_mass_r

def plot_star_map(filename, redshift, image_width=20, orientation="sideon", ax=None):
    halo_families = load.load_halo_families(filename, orientation=orientation)
    half_mass_r = calc_half_mass_radius(halo_families["star"])
    aperture_radius = 1.5*half_mass_r
    bh_xy = halo_families['bh']['pos'][:,:2]

    vel_map = generate.VelocityMap(
        halo_families["star"],
        z=redshift,
        cosmo=Planck13,
        image_width_kpc=image_width,
        aperture_kpc=aperture_radius,
        pixel_scale_arcsec=0.5,
        fwhm_arcsec=2.5,
    )
    ax = plot.plot_map(
        vel_map,
        cmap="PuOr",
        vmin=None,
        vmax=None,
        bh_xy=bh_xy,
        scalebar_size=5,
        aperture=aperture_radius,
        ax=ax
    )
    return ax

def plot_gas_map(filename, redshift, image_width=20, orientation="sideon", ax=None):
    halo_families = load.load_halo_families(filename, orientation=orientation)
    half_mass_r = calc_half_mass_radius(halo_families["star"])
    aperture_radius = 1.5*half_mass_r
    bh_xy = halo_families['bh']['pos'][:,:2]

    vel_map = generate.VelocityMap(
        halo_families["gas"],
        z=redshift,
        cosmo=Planck13,
        image_width_kpc=image_width,
        aperture_kpc=aperture_radius,
        pixel_scale_arcsec=0.5,
        fwhm_arcsec=2.5,
    )
    ax = plot.plot_map(
        vel_map,
        cmap="RdBu",
        vmin=None,
        vmax=None,
        bh_xy=bh_xy,
        scalebar_size=5,
        aperture=aperture_radius,
        ax=ax
    )
    return ax

if __name__ == "__main__":
    args = parse_args()
    filename = pathlib.Path(args.simulation)
    if args.redshift == 0:
        redshift = 0.03
    else:
        redshift = args.redshift
    image_width = args.image_width

    fig, ax = plt.subplots(1, 2, figsize=(8,4))
    plot_star_map(filename, redshift, image_width, orientation="sideon", ax=ax[0])
    plot_gas_map(filename, redshift, image_width, orientation="sideon", ax=ax[1])

    fig.tight_layout()
    fig.savefig("/home/ray/Desktop/r454_velmap.png")

