import pynbody


def stars_faceon(h, **kwargs):
    """

    Reposition and rotate the simulation containing the halo h to see
    h's disk face on.

    Given a simulation and a subview of that simulation (probably the
    halo of interest), this routine centers the simulation and rotates
    it so that the disk lies in the x-y plane. This gives a face-on
    view for SPH images, for instance.

    """

    return stars_sideon(h, pynbody.analysis.angmom.calc_faceon_matrix, **kwargs)


def stars_sideon(
    h,
    vec_to_xform=pynbody.analysis.angmom.calc_sideon_matrix,
    cen_size="1 kpc",
    disk_size="5 kpc",
    cen=None,
    vcen=None,
    move_all=True,
    **kwargs
):
    """

    Reposition and rotate the simulation containing the halo h to see
    h's disk edge on.

    Given a simulation and a subview of that simulation (probably the
    halo of interest), this routine centers the simulation and rotates
    it so that the disk lies in the x-z plane. This gives a side-on
    view for SPH images, for instance.

    """
    from pynbody import transformation, filt
    from pynbody.analysis import halo

    global config

    if move_all:
        top = h.ancestor
    else:
        top = h

    # Top is the top-level view of the simulation, which will be
    # transformed

    if cen is None:
        # or h['pos'][h['phi'].argmin()]
        cen = halo.center(h, retcen=True, **kwargs)

    tx = transformation.inverse_translate(top, cen)

    if vcen is None:
        vcen = halo.vel_center(h, retcen=True, cen_size=cen_size)

    tx = transformation.inverse_v_translate(tx, vcen)

    # Use stars from inner 10kpc to calculate angular momentum vector
    if len(h.s) > 0:
        cen = h.s[filt.Sphere(disk_size)]
    else:
        cen = h[filt.Sphere(disk_size)]

    trans = vec_to_xform(pynbody.analysis.angmom.ang_mom_vec(cen))

    tx = transformation.transform(tx, trans)

    return tx


def load_halo_families(filename, orientation="sideon"):
    """Load all particles families from the halo.

    Args:
        filename (str): Full path to simulation file.
        orientation (str, optional): Orientation of halo, feeds into pynbody. Defaults to "sideon".

    Returns:
        dict: Dictionary containing particles families.
    """
    halo = pynbody.load(filename)
    halo.physical_units()
    rotate = {
        "sideon": stars_sideon,
        "faceon": stars_faceon,
    }
    rotate[orientation](halo, mode="pot")
    return load_particles(halo)


def load_particles(sim):
    families = sim.families()
    family_names = [family.name for family in families]
    family_particles = [sim[family] for family in families]
    family_dict = dict(zip(family_names, family_particles))
    if "bh" not in family_dict.keys():
        try:
            family_dict["bh"] = sim.s[sim.s["tform"] < 0]
        except KeyError:
            print("Cannot find BH")
    return family_dict
