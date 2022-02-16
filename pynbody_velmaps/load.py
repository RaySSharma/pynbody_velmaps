import pynbody


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
        "sideon": pynbody.analysis.angmom.sideon,
        "faceon": pynbody.analysis.angmom.faceon,
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
