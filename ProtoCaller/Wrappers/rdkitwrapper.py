import copy as _copy
import os as _os

from rdkit import Chem as _Chem
from rdkit.Chem import AllChem as _AllChem
from rdkit.Chem import rdForceFieldHelpers as _FF
from rdkit.Chem import rdmolops as _rdmolops
from rdkit.Geometry import rdGeometry as _Geom
import parmed as _pmd

import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Utils.runexternal as _runexternal
import ProtoCaller.Utils.stdio as _stdio
with _stdio.stdout_stderr_cls():
    # Here we use the deprecated function because it is more robust TODO: fix
    from rdkit.Chem import MCS as _MCS
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


def getMCSMap(ref, mol, atomCompare="any", bondCompare="any", **kwargs):
    """
    Generates the Maximum Common Substructure (MCS) mapping between two molecules.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    atomCompare : str
        One of "any" (matches any pairs of atoms) and "elements" (matches only the same atoms).
    bondCompare : str
        One of "any" (matches any bonds) and "elements" (matches only the same bonds).
    kwargs
        Additional keyword arguments passed on to rdkit.Chem.MCS.FindMCS.

    Returns
    -------

    mcs : [[tuple]]
        A list of lists tuples corresponding to the atom index matches between the reference and the other molecule.
    """
    def matchAndReturnMatches(*args, **kwargs):
        # returns the MCS of two molecules as given by RDKit
        mcs_string = _MCS.FindMCS(*args, **kwargs).smarts
        if mcs_string is None:
            return None, [], []
        mcs = _Chem.MolFromSmarts(mcs_string)
        # fixes some issues with RDKit concerning SMARTS
        mcs = _Chem.MolFromSmarts(_Chem.MolToSmiles(mcs, canonical=False, allBondsExplicit=True, allHsExplicit=True))
        return mcs, ref.GetSubstructMatches(mcs), mol.GetSubstructMatches(mcs)

    def carbonify(mol):
        # converts all atoms in the molecule to carbons
        mol_c = _copy.deepcopy(mol)
        for atom in mol_c.GetAtoms():
            atom.SetAtomicNum(6)
        _Chem.AssignStereochemistry(mol_c, cleanIt=True, force=True)
        return mol_c

    def generateFragment(mol, indices_to_delete, getAsFrags=False):
        # generates a molecule from another given indices to delete
        mol_frag_edit = _Chem.EditableMol(_copy.deepcopy(mol))
        for idx in reversed(sorted(indices_to_delete)):
            mol_frag_edit.RemoveAtom(idx)
        mol_frag = mol_frag_edit.GetMol()
        if getAsFrags:
            frags = _rdmolops.GetMolFrags(mol_frag, asMols=False, sanitizeFrags=False)
            return [transformIndices(f, indices_to_delete) for f in frags]
        else:
            mol_frag.UpdatePropertyCache()
            return mol_frag

    def transformIndices(indices_new, indices_to_delete):
        # take into account the fact that indices change when atoms are deleted
        if not indices_new:
            return []
        if not indices_to_delete:
            return indices_new
        remaining_indices = [x for x in range(max(indices_to_delete) + max(indices_new) + 2) if x not in indices_to_delete]
        return [remaining_indices[x] for x in indices_new]

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

    def fixMCS(ref, mol, match_ref, match_mol, **kwargs):
        kwargs = {**kwargs, 'ringMatchesRingOnly': True}
        # make sure that correct mapping is obtained between the two molecules given an MCS
        indices_to_delete_ref = [x for x in range(ref.GetNumAtoms()) if x not in match_ref]
        indices_to_delete_mol = [x for x in range(mol.GetNumAtoms()) if x not in match_mol]
        mcs_ref = generateFragment(ref, indices_to_delete_ref)
        mcs_mol = generateFragment(mol, indices_to_delete_mol)
        mcs, _, _ = matchAndReturnMatches([mcs_ref, mcs_mol], completeRingsOnly=True, **kwargs)
        if mcs is None:
            return [[]]
        matches = [list(zip(x, y)) for x in ref.GetSubstructMatches(mcs) for y in mol.GetSubstructMatches(mcs)]

        matches_new = []

        # check if both atoms are chiral and flag the atoms if this is the case
        for match in matches:
            _Chem.AssignStereochemistry(ref, cleanIt=True, force=True)
            _Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            chiral_ref = dict(_Chem.FindMolChiralCenters(ref))
            chiral_mol = dict(_Chem.FindMolChiralCenters(mol))
            indices_chiral = [i_mcs for i_mcs, (i_ref, i_mol) in enumerate(match)
                              if i_ref in chiral_ref.keys() and i_mol in chiral_mol.keys()]

            # delete invalid atoms
            if indices_chiral:
                # take care of R/S changing with different atom types
                ref_c = carbonify(ref)
                mol_c = carbonify(mol)
                chiral_ref_c = dict(_Chem.FindMolChiralCenters(ref_c))
                chiral_mol_c = dict(_Chem.FindMolChiralCenters(mol_c))

                indices_chiral_final = []
                for idx in indices_chiral:
                    i_ref, i_mol = match[idx]
                    if chiral_ref_c[i_ref] != chiral_mol_c[i_mol]:
                        indices_chiral_final += [idx]

                if indices_chiral_final:
                    match = fixChiralIndices(indices_chiral_final, mcs, match)

            matches_new += [match]

        return matches_new

    len_ref = ref.GetNumAtoms()
    len_mol = mol.GetNumAtoms()
    if not all([len_ref, len_mol]):
        return [[]]
    elif any([len_ref == 1, len_mol == 1]):
        return [[(x, y)] for x in range(len_ref) for y in range(len_mol)]

    if "iterate" in kwargs.keys():
        iterate = True
        kwargs.pop("iterate")
    else:
        iterate = False

    kwargs = {**kwargs, 'atomCompare': atomCompare, 'bondCompare': bondCompare}

    # get initial pruned MCS
    matches = []
    mcs, matches_ref, matches_mol = matchAndReturnMatches([ref, mol], completeRingsOnly=False, **kwargs)
    if mcs is None:
        return [[]]
    for match_ref in matches_ref:
        for match_mol in matches_mol:
            matches += fixMCS(ref, mol, match_ref, match_mol, **kwargs)

    # only keep the biggest MCS's
    max_len = max([len(match) for match in matches])
    matches = {frozenset(match) for match in matches if len(match) == max_len}

    # split the remainder of the molecule and recursively search for more MCS's
    matches_new = _copy.copy(matches)
    for match in matches:
        if match == frozenset():
            continue
        match_ref, match_mol = zip(*match)
        frags_ref = generateFragment(ref, match_ref, getAsFrags=True)
        frags_mol = generateFragment(mol, match_mol, getAsFrags=True)
        for frag_ref in frags_ref:
            indices_to_delete_ref = [x for x in range(ref.GetNumAtoms()) if
                                     x not in frag_ref]
            for frag_mol in frags_mol:
                indices_to_delete_mol = [x for x in range(mol.GetNumAtoms()) if
                                         x not in frag_mol]
                mol_new = generateFragment(mol, indices_to_delete_mol)
                ref_new = generateFragment(ref, indices_to_delete_ref)
                match_new_total = getMCSMap(ref_new, mol_new, iterate=True,
                                            **kwargs)
                if match_new_total == {frozenset()}:
                    continue
                for match_new in match_new_total:
                    match_new_ref, match_new_mol = zip(*match_new)
                    match_new_ref = transformIndices(match_new_ref,
                                                     indices_to_delete_ref)
                    match_new_mol = transformIndices(match_new_mol,
                                                     indices_to_delete_mol)

                    # only include if mappings are unique
                    #if not any([len(set(match_ref).intersection(match_new_ref)),
                    #            len(set(match_mol).intersection(match_new_mol))]):
                    # TODO: deal with repeated indices
                    matches_new |= {m | frozenset(zip(match_new_ref, match_new_mol)) for m in matches}

    # don't prune disconnected MCS's just yet
    if iterate:
        return matches_new if len(matches_new) else {frozenset()}

    # prune all of the resulting MCS combinations so that they are connected
    matches_final = [[]]
    for match_new in matches_new:
        if match_new == frozenset():
            continue
        match_new_ref, match_new_mol = zip(*match_new)
        indices_to_delete_ref = [x for x in range(ref.GetNumAtoms()) if
                                 x not in match_new_ref]
        indices_to_delete_mol = [x for x in range(mol.GetNumAtoms()) if
                                 x not in match_new_mol]
        frags_ref = generateFragment(ref, indices_to_delete_ref, getAsFrags=True)
        frags_mol = generateFragment(mol, indices_to_delete_mol, getAsFrags=True)
        len_frags_ref = len(frags_ref)
        len_frags_mol = len(frags_mol)

        # only keep MCS if it is connected
        # TODO: deal with multiple fragments
        if len_frags_ref < 2 and len_frags_mol < 2:
            n_atoms = len(match_new)
            if n_atoms == 0:
                matches_final += [[]]
            elif n_atoms == 1:
                matches_final += [(x, y) for x in range(len_ref) for y in range(len_mol)]
            else:
                matches_final += fixMCS(ref, mol, match_new_ref, match_new_mol, **kwargs)

    # only keep largest unique MCS's
    max_len = max([len(match) for match in matches_final])
    matches = {frozenset(match) for match in matches_final if len(match) == max_len}
    matches = [list(match) for match in matches]

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
        # TODO: fix
        mcs = getMCSMap(ref, mol, **kwargs)[0]

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
