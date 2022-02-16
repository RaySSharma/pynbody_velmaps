import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from . import position_angles


def plot_map(vel_map, cmap="PuOr", vmin=None, vmax=None, scalebar_size=5, ax=None):
    """Plot velocity map with reasonable defaults.

    Args:
        vel_map (array-like): Pixel-by-pixel velocity map.
        cmap (str, optional): Colormap of output image. Defaults to "PuOr".
        vmin (float, optional): Minimum (toward) velocity. Setting to None will automatically calculate limit. Defaults to None.
        vmax (float, optional): Maximum (away) velocity. Setting to None will automatically calculate limit. Defaults to None.
        scalebar_size (float, optional): Physical size of scale bar in kpc. Defaults to 5.
        ax (matplotlib.axes.Axes, optional): Matplotlib axes object in which to plot velocity map. Setting to None will generate a new Axes object. Defaults to None.

    Returns:
        matplotlib.axes.Axes: Matplotlib axes object used to plot.
    """
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

    fontprops = fm.FontProperties(size=14, family="monospace")
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
        fontproperties=fontprops,
    )

    ax.add_artist(scalebar)
    ax.tick_params(axis="both", which="both", width=0, labelsize=0)
    return ax


def plot_aperture(aperture, ax, cen=(0, 0)):
    """Plot a dashed circle outlinining the location of the aperture.

    Args:
        aperture (float): Radius of aperture in pixels.
        ax (matplotlib.axes.Axes): Matplotlib axes object.
        cen (tuple, optional): Location of aperture center, in pixels. Defaults to (0, 0).

    Returns:
        matplotlib.axes.Axes: Matplotlib axes object.
    """
    aperture_outline = plt.Circle(
        cen, aperture, color="k", fill=False, linestyle="--", lw=3
    )
    ax.add_patch(aperture_outline)
    return ax


def plot_bh(bh_xy, ax):
    """Plot locations of black holes with black dots.

    Args:
        bh_xy (array-like): Length N array containing (x,y) coordinates of black holes, in pixel coordinates.
        ax (matplotlib.axes.Axes): Matplotlib axes object.

    Returns:
        matplotlib.axes.Axes: Matplotlib axes object
    """
    for (bh_x, bh_y) in bh_xy:
        ax.plot(bh_x, bh_y, color="k", marker="o", ls="none")
    return ax


def plot_pa(velmap, pa, ax):
    """Plots a line for the position angle fit from pafit

    Args:
        velmap (array-like): Pixel-by-pixel velocity map.
        pa (tuple): Position angle tuple output by pafit.fit_kinematic_pa()
        ax (matplotlib.axes.Axes): Matplotlib axes object.
    """
    import numpy as np

    angBest, angErr, vSyst = pa
    x, y = position_angles.infer_coordinates(velmap)
    rad = np.sqrt(np.max(x ** 2 + y ** 2))
    ang = [0, np.pi] + np.radians(angBest)
    ax.plot(
        rad * np.cos(ang), rad * np.sin(ang), "k--", linewidth=3
    )  # Zero-velocity line
    ax.plot(
        -rad * np.sin(ang), rad * np.cos(ang), color="limegreen", linewidth=3
    )  # Major axis PA
    return ax