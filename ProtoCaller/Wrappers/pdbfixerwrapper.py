import copy as _copy
import os as _os
import warnings as _warnings

import pdbfixer as _pdbfix
from simtk.openmm.app import PDBFile as _PDBFile

from ProtoCaller.IO.PDB import PDB as _PDB


def pdbfixerTransform(filename, add_missing_residues, add_missing_atoms):
    """
    Adds missing residues and/or missing atoms to a PDB file.

    Parameters
    ----------
    filename : str
        Name of the input PDB file.
    add_missing_residues : bool
        Whether to add missing residues.
    add_missing_atoms : bool
        Whether to add missing atoms.

    Returns
    -------
    filename_output : str
        Absolute path to the modified file.
    """
    if not add_missing_atoms and not add_missing_residues:
        return _os.path.abspath(filename)

    fix = _pdbfix.PDBFixer(filename=filename)

    if add_missing_residues:
        fix.findMissingResidues()
    else:
        fix.missingResidues = []

    if add_missing_atoms:
        fix.findMissingAtoms()
    else:
        fix.missingAtoms = []
        fix.missingTerminals = []

    fix.addMissingAtoms()

    filename_output = _os.path.splitext(filename)[0] + "_pdbfixer.pdb"
    _PDBFile.writeFile(fix.topology, fix.positions, open(filename_output, "w"))

    return fixPDBFixerPDB(filename_output, filename, add_missing_residues,
                          add_missing_atoms, filename_output)


def fixPDBFixerPDB(filename_modified, filename_original, add_missing_residues,
                   add_missing_atoms, filename_output=None):
    """
    Used to regenerate some data lost by PDBFixer.

    Parameters
    ----------
    filename_modified : str
        Name of the modified PDB file.
    filename_original : str
        Name of the original PDB file.
    add_missing_residues : bool
        Whether to add missing residues.
    add_missing_atoms : bool
        Whether to add missing atoms.
    filename_output : str
        Name of the fixed output PDB file.

    Returns
    -------
    filename_output : str
        The absolute path to the fixed output PDB file.
    """
    pdb_original = _PDB(filename_original)
    pdb_modified = _PDB(filename_modified)

    all_res_mod = pdb_modified.totalResidueList()
    all_res_orig = pdb_original.totalResidueList()
    if len(all_res_orig) != len(all_res_mod):
        raise ValueError("Mismatch between original number of residues ({}) "
                         "and number of residues output by PDBFixer ({}).".
                         format(len(all_res_orig), len(all_res_mod)))

    if add_missing_residues:
        for miss_res in pdb_original.missing_residues:
            fixed_res = _copy.copy(all_res_mod[all_res_orig.index(miss_res)])
            fixed_res.chainID = miss_res.chainID
            fixed_res.resSeq = miss_res.resSeq
            if fixed_res.resName != miss_res.resName:
                _warnings.warn("Mismatch between original residue name ({}) "
                               "and the residue name output by Modeller ({}) "
                               "in chain {}, residue number {}.".format(
                                miss_res.resName, fixed_res.resName,
                                miss_res.chainID, miss_res.resSeq))
            chain = pdb_original.filter("chainID=='{}'".format(
                miss_res.chainID), type="chains")[0]

            if miss_res > chain[-1]:
                chain.append(fixed_res)
            else:
                for i, res in enumerate(chain):
                    if res > miss_res:
                        chain.insert(i, fixed_res)
                        break
        pdb_original.missing_residues = []

    if add_missing_atoms:
        for missing_atom in pdb_original.missing_atoms:
            filter = "chainID=='{}'&resSeq=={}&iCode=='{}'".format(
                missing_atom.chainID, missing_atom.resSeq, missing_atom.iCode)
            res = pdb_original.filter(filter)[0]
            fixed_res = all_res_mod[all_res_orig.index(res)]
            fixed_res.chainID = missing_atom.chainID
            fixed_res.resSeq = missing_atom.resSeq
            res.__init__(fixed_res)
        pdb_original.missing_atoms = []

    pdb_original.reNumberAtoms()
    pdb_original.reNumberResidues()
    if filename_output is None:
        filename_output = _os.path.splitext(pdb_original.filename)[0] + \
                          "_modified.pdb"
    pdb_original.writePDB(filename_output)

    return _os.path.abspath(filename_output)
