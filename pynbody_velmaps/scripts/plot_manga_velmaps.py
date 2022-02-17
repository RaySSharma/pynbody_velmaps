import argparse
import pathlib

import numpy as np
import pynbody
import matplotlib.pyplot as plt
from astropy.cosmology import Planck13

from pynbody_velmaps import generate, plot, load, position_angles


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "simulation", help="Path to simulation.", type=pathlib.Path
    )
    parser.add_argument(
        "redshift", help="Redshift at which to generate image.", type=float
    )
    parser.add_argument(
        "image_width",
        help="Physical size of output image, in kpc.",
        type=float,
    )
    parser.add_argument(
        "outfile",
        nargs="?",
        help="Output image name",
        type=pathlib.Path,
        default=None,
    )
    return parser.parse_args()


def calc_half_mass_radius(gal, ndim=2):
    bins = np.arange(0.5, gal["r"].max(), 0.1)
    pro = pynbody.analysis.profile.Profile(gal, ndim=ndim, bins=bins)
    half_mass = 0.5 * pro["mass_enc"].max()
    half_mass_r = pro["rbins"][abs(pro["mass_enc"] - half_mass).argmin()]
    return half_mass_r


def plot_star_map(filename, redshift, image_width=20, orientation="sideon", ax=None):
    import matplotlib.font_manager as fm

    halo_families = load.load_halo_families(filename, orientation=orientation)
    half_mass_r = calc_half_mass_radius(halo_families["star"])
    bh_xy = halo_families["bh"]["pos"][:, :2]

    vel_map = generate.VelocityMap(
        halo_families["star"],
        z=redshift,
        cosmo=Planck13,
        image_width_kpc=image_width,
        aperture_kpc=2.5 * half_mass_r,
        pixel_scale_arcsec=0.5,
        fwhm_arcsec=2.5,
    )
    ax = plot.plot_map(
        vel_map,
        cmap="PuOr",
        vmin=None,
        vmax=None,
        ax=ax,
    )
    plot.plot_bh(bh_xy / vel_map.kpc_per_pixel, ax)
    plot.plot_aperture(1.5 * half_mass_r / vel_map.kpc_per_pixel, ax)
    plot.plot_scalebar(5, vel_map.kpc_per_pixel, ax, fontproperties = fm.FontProperties(size=14, family="monospace"))
    #ax.tick_params(axis="both", which="both", width=0, labelsize=0)
    return vel_map, ax


def plot_gas_map(filename, redshift, image_width=20, orientation="sideon", ax=None):
    import matplotlib.font_manager as fm

    halo_families = load.load_halo_families(filename, orientation=orientation)
    half_mass_r = calc_half_mass_radius(halo_families["star"])
    bh_xy = halo_families["bh"]["pos"][:, :2]

    vel_map = generate.VelocityMap(
        halo_families["gas"],
        z=redshift,
        cosmo=Planck13,
        image_width_kpc=image_width,
        aperture_kpc=2.5 * half_mass_r,
        pixel_scale_arcsec=0.5,
        fwhm_arcsec=2.5,
    )
    ax = plot.plot_map(
        vel_map,
        cmap="RdBu",
        vmin=None,
        vmax=None,
        scalebar_size=5,
        ax=ax,
    )
    plot.plot_bh(bh_xy / vel_map.kpc_per_pixel, ax)
    plot.plot_aperture(1.5 * half_mass_r / vel_map.kpc_per_pixel, ax)
    plot.plot_scalebar(5, vel_map.kpc_per_pixel, ax, fontproperties = fm.FontProperties(size=14, family="monospace"))
    ax.tick_params(axis="both", which="both", width=0, labelsize=0)
    return vel_map, ax


if __name__ == "__main__":
    """ Example usage: >python plot_manga_velmaps.py /path/to/simulation.tipsy 0.05 20 /path/to/output/image.png
    """
    args = parse_args()
    filename = args.simulation.as_posix()
    out_filename = args.outfile
    if args.redshift == 0:
        redshift = 0.03
    else:
        redshift = args.redshift
    image_width = args.image_width

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    star_map, _ = plot_star_map(
        filename, redshift, image_width, orientation="sideon", ax=ax[0]
    )
    gas_map, _ = plot_gas_map(
        filename, redshift, image_width, orientation="sideon", ax=ax[1]
    )

    star_pa = position_angles.calc_pa(star_map.data.data * star_map.data.mask)
    gas_pa = position_angles.calc_pa(gas_map.data.data * gas_map.data.mask)

    plot.plot_pa(star_map, star_pa, ax[0])
    plot.plot_pa(gas_map, gas_pa, ax[1])

    if out_filename is not None:
        fig.tight_layout()
        fig.savefig(out_filename, bbox_inches="tight")