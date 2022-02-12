import pynbody


def load_halo_families(filename, orientation="sideon"):
    halo = pynbody.load(filename.as_posix())
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
