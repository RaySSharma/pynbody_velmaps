# pynbody_velmaps
Generate stellar/gas line-of-sight velocity maps and position angles using pynbody.

## Dependencies

    numpy >= 1.19.2
    scipy >= 1.5.2
    matplotlib >= 3.0.0
    astropy >= 4.0.0
    pynbody >= 1.0.2
    pafit >= 2.0.7

## Installation

To install *pynbody_velmaps*:

    > git clone git@github.com:RaySSharma/pynbody_velmaps.git
    > cd pynbody_velmaps
    > pip install -e .

## Usage

An example script can be found at `scripts/plot_manga_velmaps.py`, which generates velocity maps with MaNGA-style observing characteristics (Koudmani+2021).

To make line-of-sight stellar/gas velocity maps in MaNGA style:

    from pynbody_velmaps.scripts.plot_manga_velmaps import *

    stellar_map, stellar_pa, ax = plot_manga_map(filename, redshift, "star", image_width, orientation="sideon")

    gas_map, gas_pa, ax = plot_gas_map(filename, redshift,  "gas", image_width, orientation="sideon")

where `filename` is the simulation path, `redshift` is the desired observing redshift, and `image_width` is the width of the image in kpc. `stellar_map` and `gas_map` contain the pixel-by-pixel values of the velocity maps. `stellar_pa` and `gas_pa` contain the position angle tuples (angBest, angErr, velSys) from *pafit*. Setting z=0 maps instead to z=0.03, the mean redshift of the primary MaNGA sample.

To manually measure position angles:

    from pynbody_velmaps.position_angles import *

    stellar_pa = calc_pa(star_map)
    
    gas_pa = calc_pa(gas_map)
