import warnings as _warnings

import ProtoCaller as _PC

def charmmWrapper(params, filename, molecule_type, charge=None, *args,
                  **kwargs):
    # TODO: cofactors
    ff_path = _PC.SHAREDDIR + "/charmm-parameters/"
    if molecule_type == "protein":
        suffix = params.protein_ff.split("charmm")[-1]
        force_fields = [
            ff_path + "top_all36_prot.rtf",
            ff_path + "par_all{}_prot.prm".format(suffix),
            ff_path + "toppar_water_ions.str",
        ]
        return PDBReader(filename, *args, **kwargs) + force_fields
    elif molecule_type == "complex_anion":
        _warnings.warn("CHARMM parametrisation failed: polyatomic anions not "
                       "supported")
        return
    elif molecule_type == "complex_cation":
        _warnings.warn("CHARMM parametrisation failed: transition metals not "
                       "supported")
        return
    elif molecule_type == "cofactor":
        _warnings.warn("CHARMM parametrisation failed: cofactors not "
                       "supported")
        return
    elif molecule_type == "ligand":
        force_fields = [
            ff_path + "top_all36_cgenff.rtf",
            ff_path + "par_all36_cgenff.prm",
        ]
        return ligandReader(filename, *args, **kwargs) + force_fields
    else:
        raise ValueError("Value %s for molecule_type not supported " %
                         molecule_type)



