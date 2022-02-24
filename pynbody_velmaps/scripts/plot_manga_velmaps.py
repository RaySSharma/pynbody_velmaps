import argparse
import pathlib

import numpy as np
import pynbody
import matplotlib.pyplot as plt
from astropy.cosmology import Planck13

from pynbody_velmaps import generate, plot, load, position_angles


def calc_half_mass_radius(gal, ndim=2):
    bins = np.arange(0.5, gal["r"].max(), 0.1)
    pro = pynbody.analysis.profile.Profile(gal, ndim=ndim, bins=bins)
    half_mass = 0.5 * pro["mass_enc"].max()
    half_mass_r = pro["rbins"][abs(pro["mass_enc"] - half_mass).argmin()]
    return half_mass_r


def plot_manga_map(
    filename,
    redshift,
    particles="star",
    weights="mass",
    image_width=20,
    orientation="sideon",
    ax=None,
    vmin=None,
    vmax=None,
    cmap="PuOr",
    title=None,
    ax_labels=True,
    show_cbar=True,
):
    halo_families = load.load_halo_families(filename, orientation=orientation)

    if isinstance(weights, str):
        halo_families[particles]["map_weights"] = halo_families[particles][weights]
    elif callable(weights):
        halo_families[particles]["map_weights"] = weights(halo_families[particles])
    else:
        halo_families[particles]["map_weights"] = np.ones(len(halo_families[particles]))

    halo_families[particles]["map_qty"] = (
        halo_families[particles]["vz"] * halo_families[particles]["map_weights"]
    )

    half_mass_r = calc_half_mass_radius(halo_families["star"])
    bh_xy = halo_families["bh"]["pos"][:, :2]

    vel_map = generate.VelocityMap(
        halo_families[particles],
        z=redshift,
        cosmo=Planck13,
        image_width_kpc=image_width,
        aperture_kpc=1.5 * half_mass_r,
        pixel_scale_arcsec=0.5,
        fwhm_arcsec=None,
        halo_max_size=2.5 * half_mass_r,
    )
    ax = plot.plot_map(
        vel_map,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        ax=ax,
        title=title,
        ax_labels=ax_labels,
        show_cbar=show_cbar,
    )
    plot.plot_bh(bh_xy, ax)
    plot.plot_aperture(1.5 * half_mass_r, ax)
    plot.plot_scalebar(5, ax, size_vertical=0.1, pad=0.5, sep=10)
    plot.plot_pa(vel_map, ax)
    return vel_map, ax


if __name__ == "__main__":
    """Example usage: >python plot_manga_velmaps.py /path/to/simulation.tipsy 0.05 20 /path/to/output/image.png"""

    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument("simulation", help="Path to simulation.", type=pathlib.Path)
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

    def assign_args():
        args = parse_args()
        filename = args.simulation.as_posix()
        out_filename = args.outfile
        if args.redshift == 0:
            redshift = 0.03
        else:
            redshift = args.redshift
        image_width = args.image_width
        return filename, redshift, image_width, out_filename

    filename, redshift, image_width, out_filename = assign_args()

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))

    def rho_sq(gas):
        return gas["rho"] ** 2

    star_map, _ = plot_manga_map(
        filename,
        redshift,
        particles="star",
        weights="mass",
        image_width=image_width,
        orientation="sideon",
        ax=ax[0],
        cmap="PuOr",
    )
    gas_map, _ = plot_manga_map(
        filename,
        redshift,
        particles="gas",
        weights=rho_sq,
        image_width=image_width,
        orientation="sideon",
        ax=ax[1],
        cmap="RdBu",
    )

    ax[0].tick_params(axis="both", which="both", width=0, labelsize=0)
    ax[1].tick_params(axis="both", which="both", width=0, labelsize=0)

    if out_filename is not None:
        fig.tight_layout()
        fig.savefig(out_filename, bbox_inches="tight")