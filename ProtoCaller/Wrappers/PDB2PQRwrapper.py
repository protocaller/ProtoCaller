# TODO:
# 1. check proteins with modified AA's
# 2. check numbering of ignored residues

import ProtoCaller as _PC
import ProtoCaller.IO.PDB as _PDB
import ProtoCaller.Utils.runexternal as _runexternal

import os as _os


# all options can be found here:
# https://apbs-pdb2pqr.readthedocs.io/en/latest/pdb2pqr/invoking.html
def PDB2PQRtransform(filename_input, filename_output=None, **kwargs):
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

    pdb_obj._missing_atoms = []
    for i, chain in enumerate(pdb_obj):
        if chain.type == "chain":
            for j, residue in enumerate(chain):
                if residue.type == "amino_acid":
                    updateResidue()

    pdb_obj.reNumberResidues()
    pdb_obj.reNumberAtoms()
    pdb_obj.writePDB(filename)

    return _os.path.abspath(filename)
