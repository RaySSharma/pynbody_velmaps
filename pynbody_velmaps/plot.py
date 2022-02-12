import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from . import position_angles


def plot_map(
    vel_map,
    cmap="PuOr",
    vmin=None,
    vmax=None,
    scalebar_size=5,
    ax=None
):
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(4, 4))


    kpc_per_pixel = vel_map.image_width_kpc / vel_map.npixels
    width = vel_map.data.shape[0]
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap(cmap)

    ax.imshow(
        vel_map.data,
        extent=(-width / 2, width / 2, -width / 2, width / 2),
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        norm=norm,
    )

    fontprops = fm.FontProperties(size=14, family='monospace')
    scalebar = AnchoredSizeBar(
        ax.transData,
        scalebar_size / kpc_per_pixel,
        "{} kpc".format(scalebar_size),
        "lower right",
        pad=0.5,
        color="black",
        frameon=False,
        size_vertical=0.5,
        sep=10,
        fontproperties=fontprops
    )

    ax.add_artist(scalebar)
    ax.tick_params(axis="both", which="both", width=0, labelsize=0)
    return ax

def plot_aperture(aperture, ax, cen=(0, 0)):
    aperture_outline = plt.Circle(
        cen, aperture, color="k", fill=False, linestyle='--', lw=3
    )
    ax.add_patch(aperture_outline)
    return ax

def plot_bh(bh_xy, ax):
    for (bh_x, bh_y) in bh_xy:
        ax.plot(bh_x, bh_y, color='k', marker='o', ls='none')
    return ax

def plot_pa(velmap, pa, ax):
    import numpy as np
    angBest, angErr, vSyst = pa
    x, y = position_angles.infer_coordinates(velmap)
    rad = np.sqrt(np.max(x**2 + y**2))
    ang = [0,np.pi] + np.radians(angBest)
    ax.plot(rad*np.cos(ang), rad*np.sin(ang), 'k--', linewidth=3) # Zero-velocity line
    ax.plot(-rad*np.sin(ang), rad*np.cos(ang), color="limegreen", linewidth=3) # Major axis PA