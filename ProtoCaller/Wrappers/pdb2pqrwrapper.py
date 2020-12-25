# TODO:
# 1. check proteins with modified AA's
# 2. check numbering of ignored residues
import os as _os
import warnings as _warnings

from ProtoCaller.IO.PDB import PDB as _PDB
from ProtoCaller.Utils.runexternal import runExternal as _runExternal

__all__ = ["pdb2pqrTransform"]


def pdb2pqrTransform(filename, pdb2pqr_executable="pdb2pqr30", **kwargs):
    """
    A thin wrapper around PDB2PQR which protonates an input PDB file.

    Parameters
    ----------
    filename : str
        Name of input file.
    pdb2pqr_executable : str
        Name or path to the pdb2pqr executable.
    kwargs:
        Keyword arguments to be passed on to PDB2PQR. All supported arguments can be found at e.g. pdb2pqr30 --help.

    Returns
    -------
    filename : str
        Absolute path to the protonated file.
    """

    default_kwargs = {
        'ff': 'AMBER',
        'ffout': 'AMBER',
        'drop-water': True,
        'keep-chain': True,
        'with-ph': 7,
        'titration-state-method': 'propka',
        'log-level': "INFO",
    }

    default_kwargs = {**default_kwargs, **kwargs}
    filename_output = _os.path.splitext(filename)[0] + "_pdb2pqr.pdb"
    command = f"{pdb2pqr_executable} {filename} {filename_output}"
    for k, v in default_kwargs.items():
        if v is not False:
            command += f" --{k}"
        if v not in [True, False]:
            command += f" {v}"
    _runExternal(command, procname="PDB2PQR")

    return fixPdb2pqrPDB(filename_output, filename, filename_output)


def fixPdb2pqrPDB(filename_modified, filename_original, filename_output=None):
    """
    Used to regenerate some data lost by PDB2PQR.

    Parameters
    ----------
    filename_modified : str
        Name of the modified PDB file.
    filename_original : str
        Name of the original PDB file.
    filename_output : str
        Name of the fixed output PDB file.

    Returns
    -------
    filename_output : str
        The absolute path to the fixed output PDB file.
    """
    pdb_original = _PDB(filename_original)
    pqr_modified = _PDB(filename_modified)

    for res_orig in pdb_original.filter("type=='amino_acid'"):
        filter = "chainID=='{}'&resSeq=={}&iCode=='{}'".format(
            res_orig.chainID, res_orig.resSeq, res_orig.iCode)
        res_mod = pqr_modified.filter(filter)[0]
        res_orig.clear()
        res_orig.__init__(res_mod)
        for atom in res_orig:
            for attr in ["occupancy", "tempFactor", "element", "charge"]:
                setattr(atom, attr, "")
    pdb_original.missing_atoms = []

    for bond in pdb_original.disulfide_bonds:
        for res in bond:
            if res.resName != "CYX":
                _warnings.warn("Disulfide bond found at a residue not labelled "
                               "as CYX ({} {}{}{}). Please check your "
                               "structure. Correcting residue to CYX...".format(
                    res.resName, res.chainID, res.resSeq, res.iCode))
                res.resName = "CYX"

    pdb_original.reNumberAtoms()
    if filename_output is None:
        filename_output = _os.path.splitext(pdb_original.filename)[0] + \
                          "_modified.pdb"

    return pdb_original.writePDB(filename_output)
