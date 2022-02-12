import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


def plot_map(
    vel_map,
    cmap="PuOr",
    vmin=None,
    vmax=None,
    bh_xy=None,
    aperture=None,
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

    if aperture is not None:
        aperture_outline = plt.Circle(
            (0, 0), aperture / kpc_per_pixel, color="k", fill=False, linestyle='--', lw=3
        )
        ax.add_patch(aperture_outline)

    if bh_xy is not None:
        for (bh_x, bh_y) in bh_xy:
            ax.plot(bh_x / kpc_per_pixel, bh_y / kpc_per_pixel, color='k', marker='o', ls='none')

    ax.add_artist(scalebar)
    ax.tick_params(axis="both", which="both", width=0, labelsize=0)
    return ax
