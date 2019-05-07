from collections import Counter as _Counter
import copy as _copy
import os as _os

from rdkit import Chem as _Chem
from rdkit.Chem import AllChem as _AllChem
from rdkit.Chem import rdForceFieldHelpers as _FF
from rdkit.Chem import rdFMCS as _FMCS
from rdkit.Chem import rdmolops as _rdmolops
from rdkit.Geometry import rdGeometry as _Geom
import parmed as _pmd

import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Utils.runexternal as _runexternal
import ProtoCaller.Utils.stdio as _stdio
from . import babelwrapper as _babel


def openFileAsRdkit(filename, **kwargs):
    """
    Opens an input file and returns an RDKit molecule.

    Parameters
    ----------
    filename:
        The name of the input file.
    kwargs
        Keyword arguments to be supplied to the more specialised RDKit functions.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The file loaded as an RDKit Mol object.
    """
    extension = filename.split(".")[-1].lower()

    func_dict = {
        "mol": _Chem.MolFromMolFile,
        "mol2": _Chem.MolFromMol2File,
        "pdb": _Chem.MolFromPDBFile,
    }

    if extension in func_dict.keys():
        mol = func_dict[extension](filename, **kwargs)
    elif extension == "sdf":
        mol = _Chem.SDMolSupplier(filename, **kwargs)[0]
    else:
        raise TypeError("Unrecognised input extension for RDKit: {}".format(extension))

    if mol is None:
        raise ValueError("Invalid file format for file: {}".format(filename))
    return mol


@_stdio.stdout_stderr()
def openInChIAsRdkit(inchi_str, removeHs=True, **kwargs):
    """
    Converts an InChI string into an RDKit Mol object. This function is a simple wrapper of rdkit.Chem.MolFromInchi.

    Parameters
    ----------
    inchi_str : str
        Input InChI string.
    removeHs : bool
        Whether to remove hydrogens from the created molecule.
    kwargs
        Keyword arguments to be passed to rdkit.Chem.MolFromInchi.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The InChI loaded as an RDKit Mol object.
    """
    mol = _Chem.MolFromInchi(inchi_str, removeHs=removeHs, **kwargs)
    if mol is None:
        raise ValueError("InChI string not recognised: %s" % inchi_str)
    return mol


@_stdio.stdout_stderr()
def openSmilesAsRdkit(smiles_str, **kwargs):
    """
    Converts a SMILES string into an RDKit Mol object. This function is a simple wrapper of rdkit.Chem.MolFromSmiles.

    Parameters
    ----------
    smiles_str : str
        Input SMILES string.
    kwargs
        Keyword arguments to be passed to rdkit.Chem.MolFromSmiles.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The SMILES loaded as an RDKit Mol object.
    """
    mol = _Chem.MolFromSmiles(smiles_str, **kwargs)
    if mol is None:
        raise ValueError("SMILES string not recognised: %s" % smiles_str)
    return mol


def openAsRdkit(val, minimise=None, **kwargs):
    """
    A general wrapper which can convert a variety of representations for a molecule into an rdkit.Chem.rdchem.Mol object.

    Parameters
    ----------
    val : str
        Input value - SMILES, InChI strings or a filename.
    minimise : bool or None
        Whether to perform a GAFF minimisation using OpenBabel. None means minimisation for molecules initialised from
        strings and no minimisation for molecules initialised from files.
    kwargs
        Keyword arguments to be passed to the more specialsied RDKit functions.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The input string opened as an RDKit Mol object.
    """
    if isinstance(val, str) and len(val.split(".")) == 1:
        if minimise is None: minimise = True
        flist = [openSmilesAsRdkit, openInChIAsRdkit]
        for i, f in enumerate(flist):
            try:
                mol = f(val, **kwargs)
                _AllChem.EmbedMolecule(mol, useRandomCoords=True)
                break
            except:
                if i == len(flist) - 1:
                    raise ValueError("String not recognised as a valid SMILES or InChI input")
        if minimise:
            # make sure we preserve the chirality of the input string and minimise with obminimize / GAFF
            with _fileio.Dir("Temp", temp=True):
                smiles = _Chem.MolToSmiles(mol)
                with open("molecule.smi", "w") as f:
                    f.write(smiles)
                _babel.babelTransform("molecule.smi", output_extension="sdf", generate_3D_coords=True)
                _runexternal.runExternal("obminimize -sd -c 1e-6 -n 10000 -ff GAFF -osdf molecule.sdf > minimised.sdf")
                mol = openFileAsRdkit("minimised.sdf")
    else:
        if minimise is None: minimise = False
        try:
            mol = openFileAsRdkit(val, **kwargs)
        except:
            raise ValueError("File is not in a valid format")
        # minimise the molecule using obminimize / GAFF
        if minimise:
            with _fileio.Dir("Temp", temp=True):
                saveFromRdkit(mol, "molecule.sdf")
                _runexternal.runExternal("obminimize -sd -c 1e-6 -n 10000 -ff GAFF -osdf molecule.sdf > minimised.sdf")
                mol = openFileAsRdkit("minimised.sdf")

    return mol


def saveFromRdkit(mol, filename, **kwargs):
    """
    Saves an RDKit Mol object to a file.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input molecule.
    filename : str
        The name of the output file.
    kwargs
        Keyword arguments to be passed to more specialised RDKit functions.

    Returns
    -------
    filename : str
        The absolute path to the written file.
    """
    if _os.path.exists(filename):
        _os.remove(filename)
    extension = filename.split(".")[-1]
    if extension.lower() == "sdf":
        writer = _Chem.SDWriter(filename, **kwargs)
        writer.write(mol)
        writer.close()
    elif extension.lower() == "mol":
        _Chem.MolToMolFile(mol, filename=filename, **kwargs)
    else:
        tempfilename = _os.path.splitext(filename)[0] + ".mol"
        _Chem.MolToMolFile(mol, filename=tempfilename, **kwargs)
        if extension == "mol2":
            _babel.babelTransform(tempfilename, output_extension=extension, pH=None)
        else:
            _pmd.load_file(tempfilename)[0].save(filename)

    return _os.path.abspath(filename)


def translateMolecule(mol, vector):
    """
    Translates an input molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be translated.
    vector : (float, float, float)
        A 3D vector.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The translated molecule.
    """
    mol_conf = mol.GetConformer(-1)
    for i in range(mol.GetNumAtoms()):
        new_atom_position = mol_conf.GetAtomPosition(i) + _Geom.Point3D(*vector)
        mol_conf.SetAtomPosition(i, new_atom_position)
    return mol


def getMCSMap(ref, mol, match="any", always_maximum=False, matchChiralTag=True, **kwargs):
    """
    Generates the Maximum Common Substructure (MCS) mapping between two molecules.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    match : str
        One of "any" (matches any pairs of atoms) and "elements" (matches only the same atoms).
    always_maximum : bool
        Always give the substructure that has the maximum common number of atoms - good for single
       topology mapping. False is good for overall structure alignment. Note that in some more
       complex molecules there might be an MCS that satisfies our criteria and is bigger. Here we
       rely on there being one obvious MCS - something that should be the case in free energy
       calculations anyway. Note that this option does not apply to chiral symmetric molecules.
    matchChiralTag : bool
        Take care of stereochemistry. Here we don't use rdkit's default behaviour because it is not
       good for our purposes. We instead opt for a maximum common substructure and then pruning all
       atoms that will not obey the stereochemistry. While this option does not work perfectly for
       symmetric chiral molecules, it should work rather well in the general case. Again, the
       results are best when the two molecules are similar in structure which should be the case
       in FE calculations.
    kwargs
        Keyword arguments passed on to rdkit.Chem.rdFMCS.FindMCS.

    Returns
    -------

    mcs : [tuple]
        A list of tuples corresponding to the atom index matches between the reference and the other molecule.
    """
    def matchAndReturnMatches(*args, **kwargs):
        mcs_string = _FMCS.FindMCS(*args, **kwargs).smartsString
        mcs = _Chem.MolFromSmarts(mcs_string)
        return mcs, list(zip(ref.GetSubstructMatch(mcs), mol.GetSubstructMatch(mcs)))

    def transformIndices(indices_new, indices_to_delete):
        # take into account the fact that indices change when atoms are deleted
        return [x + sum(x >= y for y in sorted(indices_to_delete)) for x in sorted(indices_new)]

    def fixIndicesToDelete(indices_to_delete, mcs, matches):
        # accept a molecule and indices to delete and return the indices from the old molecule corresponding to the
        # biggest fragment in the new one
        mcs_edit = _Chem.EditableMol(_copy.deepcopy(mcs))
        for index in reversed(indices_to_delete):
            mcs_edit.RemoveAtom(index)
        indices_new = max(_rdmolops.GetMolFrags(mcs_edit.GetMol(), asMols=False, sanitizeFrags=False), key=len)
        return [x for i, x in enumerate(matches) if i in transformIndices(indices_new, indices_to_delete)]

    def fixChiralIndices(indices, mcs, matches):
        # delete all specified indices and create fragments recursively so that each one of them contains at least one
        # of the indices
        fragments = {frozenset(i for i in range(mcs.GetNumAtoms()))}
        for index in indices:
            fragments_new = set()
            for fragment in fragments:
                # every fragment is stored as a set of indices. Here we recreate this fragment from the MCS
                # then we fragment this fragment further by deleting the index from it
                if index not in fragment:
                    fragments_new.add(fragment)
                    continue
                indices_to_delete = [x for x in range(mcs.GetNumAtoms()) if x == index or x not in fragment]
                mcs_edit = _Chem.EditableMol(_copy.deepcopy(mcs))
                for index_to_delete in reversed(indices_to_delete):
                    mcs_edit.RemoveAtom(index_to_delete)
                # the new fragments resulting from the fragment have new indices that need to be translated back
                fragment_new = {frozenset(i for i in transformIndices(x, indices_to_delete) + [index])
                                for x in _rdmolops.GetMolFrags(mcs_edit.GetMol(), asMols=False, sanitizeFrags=False)}
                fragments_new |= fragment_new
            fragments = fragments_new
        return [x for i, x in enumerate(matches) if i in max(fragments, key=len)]

    def fixMCS(mcs, matches):
        # check whether there are any improperly assigned rings (e.g. 5-membered to 6-membered) and flag the atoms
        # also check if there is any chirality mismatch and flag the atoms if this is the case
        indices_to_delete, indices_chiral = [], []
        for i_mcs, (i_ref, i_mol) in enumerate(matches):
            rings_ref, rings_mol = getRings(i_ref, ref), getRings(i_mol, mol)
            is_wrong_mapping = not isSublist(rings_ref, rings_mol) and not isSublist(rings_mol, rings_ref)
            is_different_chirality = matchChiralTag and i_ref in chiral_ref.keys() and i_mol in chiral_mol.keys() and \
                                     chiral_ref[i_ref] != chiral_mol[i_mol]
            if is_wrong_mapping:
                indices_to_delete += [i_mcs]
            elif is_different_chirality:
                indices_chiral += [i_mcs]
        # delete invalid atoms
        if indices_to_delete:
            matches = fixIndicesToDelete(indices_to_delete, mcs, matches)
        if indices_chiral:
            matches = fixChiralIndices(indices_chiral, mcs, matches)
        return matches

    match = match.strip().lower()
    if match == "any":
        comp_method = _FMCS.AtomCompare.CompareAny
    elif match == "elements":
        comp_method = _FMCS.AtomCompare.CompareElements
    else:
        raise ValueError("'match' arguments needs to be either 'any' or 'elements'")

    kwargs = {
        "atomCompare": comp_method,
        "bondCompare": _FMCS.BondCompare.CompareAny,
        "ringMatchesRingOnly": False,
        **kwargs,
    }

    # get the chiral centres
    chiral_ref = dict(_Chem.FindMolChiralCenters(ref))
    chiral_mol = dict(_Chem.FindMolChiralCenters(mol))

    # find loose MCS. matchChiralTag is always False because of bad combined functionality of completeRingsOnly False
    # and matchChiralTag True. We try to fix any wrong results later on in the code.
    if always_maximum:
        mcs, matches = matchAndReturnMatches([ref, mol], completeRingsOnly=False, matchChiralTag=False, **kwargs)
        matches = fixMCS(mcs, matches)

    # now get the strict mapping and discard the loose mapping if it doesn't perform better
    mcs_strict, matches_strict = matchAndReturnMatches([ref, mol], completeRingsOnly=True, matchChiralTag=False, **kwargs)
    if matchChiralTag:
        matches_strict = fixMCS(mcs_strict, matches_strict)
    if not always_maximum or len(matches) < len(matches_strict):
        return matches_strict

    return matches


def alignTwoMolecules(ref, mol, n_min=-1, mcs=None, **kwargs):
    """
    Aligns two molecules based on an input MCS.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    n_min : int
        Minimum number of force field minimisation iterations. -1 is no limit.
    mcs : [tuple] or None
        The maximum common substucture. None means the one generated from getMCSMap.
    kwargs
        Keyword arguments passed onto getMCSMap if mcs is None.

    Returns
    -------
        mol : rdkit.Chem.rdchem.Mol
            The aligned molecule.
        mcs : [tuple]
            The maximum common substructure.
    """
    if mcs is None:
        mcs = getMCSMap(ref, mol, **kwargs)

    ref_conf = ref.GetConformer(-1)
    mol_conf = mol.GetConformer(-1)

    ff = _FF.UFFGetMoleculeForceField(mol)
    for i_ref, i_mol in mcs:
        mol_conf.SetAtomPosition(i_mol, ref_conf.GetAtomPosition(i_ref))
        ff.AddFixedPoint(i_mol)

    if mol_conf.GetNumAtoms() != len(mcs):
        ff.Initialize()
        more = ff.Minimize()
        while more and n_min:
            more = ff.Minimize()
            # n_min = -1 means infinite minimisation
            if n_min != -1:
                n_min -= 1

    return mol, mcs


def getRings(atom_idx, molecule):
    """
    Returns an array with the sorted ring sizes of all rings the target atom is contained in.

    Parameters
    ----------
    atom_idx : int
        Index of the target atom.
    molecule : rdkit.Chem.rdchem.Mol
        The input molecule.

    Returns
    -------
    rings : list
        A sorted list of all rings the atom is a part of.
    """
    rings = molecule.GetRingInfo().AtomRings()
    return sorted([len(ring) for ring in rings if atom_idx in ring])


def isSublist(l1, l2):
    """
    Returns whether a list is a sublist of another, accounting for repeated elements.

    Parameters
    ----------
    l1 : list
        The sublist.
    l2 : list
        The main list.

    Returns
    -------
    value : bool:
        Whether l1 is contained in l2.
    """
    c1, c2 = _Counter(l1), _Counter(l2)
    for k, n in c1.items():
        if n > c2[k]:
            return False
    return True