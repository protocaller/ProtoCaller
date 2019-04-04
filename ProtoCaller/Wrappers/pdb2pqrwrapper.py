# TODO:
# 1. check proteins with modified AA's
# 2. check numbering of ignored residues
import os as _os

import ProtoCaller as _PC
from ProtoCaller.IO.PDB import PDB as _PDB
import ProtoCaller.Utils.runexternal as _runexternal

__all__ = ["pdb2pqrTransform"]


# all options can be found here:
# https://apbs-pdb2pqr.readthedocs.io/en/latest/pdb2pqr/invoking.html
def pdb2pqrTransform(filename, **kwargs):
    """
    A thin wrapper around PDB2PQR which protonates an input PDB file.

    Parameters
    ----------
    filename : str
        Name of input file.
    kwargs:
        Keyword arguments to be passed on to PDB2PQR. All options can be found here:
        https://apbs-pdb2pqr.readthedocs.io/en/latest/pdb2pqr/invoking.html

    Returns
    -------
    filename : str
        Absolute path to the protonated file.
    """
    if not _PC.PDB2PQREXE:
        raise OSError("Cannot find PDB2PQR binary for the current operating "
                      "system.")
    filename_output = _os.path.splitext(filename)[0] + "_pdb2pqr.pdb"

    default_kwargs = {
        'chain': True,
        'ff': 'amber',
        'ffout': 'amber',
        'verbose': True,
        'drop-water': True,
    }
    for key, val in kwargs.items():
        default_kwargs[key] = val

    runstring = "\"{}\" \"{}\" \"{}\"".format(_PC.PDB2PQREXE, filename,
                                              filename_output)
    for key, val in default_kwargs.items():
        if val is True:
            runstring += " --%s" % key
        elif val is not False:
            runstring += " --%s=%s" % (key, val)
    _runexternal.runExternal(runstring, procname="PDB2PQR")

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

    pdb_original.reNumberAtoms()
    if filename_output is None:
        filename_output = _os.path.splitext(pdb_original.filename)[0] + \
                          "_modified.pdb"

    return pdb_original.writePDB(filename_output)
