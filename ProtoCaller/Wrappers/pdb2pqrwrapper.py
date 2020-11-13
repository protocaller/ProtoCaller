# TODO:
# 1. check proteins with modified AA's
# 2. check numbering of ignored residues
from contextlib import redirect_stdout as _redirect_stdout
import logging as _logging
import os as _os
import sys as _sys
import warnings as _warnings

from pdb2pqr.main import runPDB2PQR as _runPDB2PQR
from pdb2pqr.src.pdb import readPDB as _readPDB
import propka.lib as _lib

from ProtoCaller.IO.PDB import PDB as _PDB

__all__ = ["pdb2pqrTransform"]


def pdb2pqrTransform(filename, **kwargs):
    """
    A thin wrapper around PDB2PQR which protonates an input PDB file.

    Parameters
    ----------
    filename : str
        Name of input file.
    kwargs:
        Keyword arguments to be passed on to PDB2PQR. Supported arguments are:

        ph:
            The desired ph of the system (float).
        verbose:
            Whether all information is output to the logger.
        ph_calc_method:
            pKa calculation method (None, "propka","propka31","pdb2pka").
        ph_calc_options:
            optionParser like option object for propka30.
        neutraln:
            Make the N-terminus of this protein neutral
        neutralc:
            Make the C-terminus of this protein neutral
        assign_only:
            Only assign charges and radii - do not add atoms, debump, or optimize.
        drop_water:
            Remove water molecules from output.
        debump:
            When 1, debump heavy atoms (int).
        opt:
            When 1, run hydrogen optimization (int).
        include_old_header:
            Include most of the PDB header in output.
        holdlist:
            A list of residues not to be optimized, as [(resid, chain, icode)].

    Returns
    -------
    filename : str
        Absolute path to the protonated file.
    """

    default_kwargs = {
        'chain': True,
        'ff': 'amber',
        'ffout': 'amber',
        'verbose': True,
        'drop_water': True,
        'ph': 7,
        'ph_calc_method': 'propka31',
        'ph_calc_options': _lib.loadOptions('--quiet')[0]
    }

    default_kwargs = {**default_kwargs, **kwargs}

    def write(msg):
        msg = msg.strip()
        if msg in ["", "*"] or "100%" in msg or "progress" in msg:
            return None
        if msg[-1] == ":":
            msg = msg[:-1] + "..."
        return _logging.info(msg)

    _logging.basicConfig(stream=_sys.stdout)
    _logging.write = write
    with _redirect_stdout(_logging):
        pdb = _readPDB(open(filename))[0]
        pdb = _runPDB2PQR(pdb, **default_kwargs)

    filename_output = _os.path.splitext(filename)[0] + "_pdb2pqr.pdb"
    with open(filename_output, "w") as f:
        for line in pdb['lines']:
            f.write(line)

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
