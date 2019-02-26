# TODO:
# 1. check proteins with modified AA's
# 2. check numbering of ignored residues

import ProtoCaller as _PC
import ProtoCaller.IO.PDB as _PDB
import ProtoCaller.Utils.runexternal as _runexternal

import os as _os

__all__ = ["PDB2PQRtransform"]


# all options can be found here:
# https://apbs-pdb2pqr.readthedocs.io/en/latest/pdb2pqr/invoking.html
def PDB2PQRtransform(filename_input, filename_output=None, **kwargs):
    """
    A thin wrapper around PDB2PQR which protonates an input PDB file.

    Parameters
    ----------
    filename_input : str
        Name of input file.
    filename_output : str
        Name of output file. Default: 'inputfile_pdb2pqr.pdb'.
    kwargs:
        Keyword arguments to be passed on to PDB2PQR. All options can be found here:
        https://apbs-pdb2pqr.readthedocs.io/en/latest/pdb2pqr/invoking.html

    Returns
    -------
    filename : str
        Absolute path to the protonated file.
    """
    if not _PC.PDB2PQREXE:
        raise OSError("Cannot find PDB2PQR binary for the current operating system.")
    if filename_output is None:
        filename_output = _os.path.splitext(filename_input)[0] + "_pdb2pqr.pqr"

    default_kwargs = {'chain': True, 'ff': 'amber', 'ffout': 'amber', 'verbose': True, 'drop-water': True}
    for key, val in kwargs.items():
        default_kwargs[key] = val

    runstring = " ".join([_PC.PDB2PQREXE, filename_input, filename_output])
    for key, val in default_kwargs.items():
        if val is True:
            runstring += " --%s" % key
        elif val is not False:
            runstring += " --%s=%s" % (key, val)
    _runexternal.runExternal(runstring, procname="PDB2PQR")

    return fixPDB2PQRPQR(_PDB.PDB(filename_input), filename_output)


def fixPDB2PQRPQR(pdb_obj, pqr_file, filebase=None):
    """
    Used to regenerate some PDB data lost by PDB2PQR.

    Parameters
    ----------
    pdb_obj : ProtoCaller.IO.PDB.PDB
        An object of the original unmodified PDB file.
    pqr_file : str
        The path to the PDB2PQR output PDB file.
    filebase: str
        The base name of the fixed output PDB file.


    Returns
    -------
    filename : str
        The absolute path to the fixed output PDB file.

    """
    if filebase is None:
        filebase = _os.path.splitext(pdb_obj.filename)[0] + "_pdb2pqr"
    filename = filebase + ".pdb"

    def updateResidue():
        nonlocal residue, pqr_obj

        protonated_residue = pqr_obj.filter("resSeq==%d&iCode=='%s'" % (residue.resSeq, residue.iCode))[0]
        residue.clear()
        residue.__init__(protonated_residue)
        for atom in residue:
            for attr in ["occupancy", "tempFactor", "element", "charge"]:
                setattr(atom, attr, "")

    pqr_obj = _PDB.PDB(pqr_file)

    pdb_obj.missing_atoms = []
    for i, chain in enumerate(pdb_obj):
        if chain.type == "chain":
            for j, residue in enumerate(chain):
                if residue.type == "amino_acid":
                    updateResidue()

    pdb_obj.reNumberResidues()
    pdb_obj.reNumberAtoms()
    pdb_obj.writePDB(filename)

    return _os.path.abspath(filename)
