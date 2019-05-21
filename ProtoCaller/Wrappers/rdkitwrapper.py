import collections as _collections
import copy as _copy
import functools as _functools
import operator as _operator
import os as _os

import numpy as _np
from rdkit import Chem as _Chem
from rdkit.Chem import AllChem as _AllChem
from rdkit.Chem import rdForceFieldHelpers as _FF
from rdkit.Chem import rdmolops as _rdmolops
from rdkit.Chem import rdMolTransforms as _Transforms
from rdkit.Geometry import rdGeometry as _Geom
import parmed as _pmd
from scipy.optimize import minimize as _minimize

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


def getMCSMap(ref, mol, atomCompare="any", bondCompare="any", maxRecursions=3,
              **kwargs):
    """
    Generates the Maximum Common Substructure (MCS) mapping between two molecules.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    atomCompare : str
        One of "any" (matches any pairs of atoms) and "elements" (matches only
        the same atoms).
    bondCompare : str
        One of "any" (matches any bonds) and "elements" (matches only the same
        bonds).
    kwargs
        Additional keyword arguments passed on to rdkit.Chem.MCS.FindMCS.

    Returns
    -------
    mcs : {frozenset([tuple])}
        A set of frozensets of tuples corresponding to the atom index matches
        between the reference and the other molecule.
    """
    if maxRecursions == 0:
        return {frozenset()}

    len_ref = ref.GetNumAtoms()
    len_mol = mol.GetNumAtoms()
    if not all([len_ref, len_mol]):
        return {frozenset()}
    elif any([len_ref == 1, len_mol == 1]):
        return {frozenset([(x, y)] )for x in range(len_ref)
                for y in range(len_mol)}

    if "iterate" in kwargs.keys():
        iterate = True
        kwargs.pop("iterate")
    else:
        iterate = False

    kwargs = {**kwargs, 'atomCompare': atomCompare, 'bondCompare': bondCompare}

    # get initial pruned MCS
    mcs, matches_ref, matches_mol = _matchAndReturnMatches([ref, mol], **kwargs)
    if mcs is None:
        return {frozenset()}

    # get initial fixed matches
    matches = set()
    for match_ref in matches_ref:
        for match_mol in matches_mol:
            matches |= _fixMCS(ref, mol, match_ref, match_mol, **kwargs)

    # only keep the biggest MCS's
    max_len = max([len(match) for match in matches])
    matches = {match for match in matches if len(match) == max_len}

    # split the remainder of the molecule and recursively search for more MCS's
    matches_new = _copy.copy(matches)
    for match in matches:
        if match == frozenset():
            continue
        match_ref, match_mol = zip(*match)
        frags_ref = _generateFragment(ref, match_ref, getAsFrags=True)
        frags_mol = _generateFragment(mol, match_mol, getAsFrags=True)
        for frag_ref in frags_ref:
            indices_to_delete_ref = [x for x in range(ref.GetNumAtoms()) if
                                     x not in frag_ref]
            for frag_mol in frags_mol:
                indices_to_delete_mol = [x for x in range(mol.GetNumAtoms()) if
                                         x not in frag_mol]
                mol_new = _generateFragment(mol, indices_to_delete_mol)
                ref_new = _generateFragment(ref, indices_to_delete_ref)
                match_new_total = getMCSMap(ref_new, mol_new, iterate=True,
                                            maxRecursions=maxRecursions-1,
                                            **kwargs)
                if match_new_total == {frozenset()}:
                    continue
                for match_new in match_new_total:
                    match_new_ref, match_new_mol = zip(*match_new)
                    match_new_ref = _transformIndices(match_new_ref,
                                                     indices_to_delete_ref)
                    match_new_mol = _transformIndices(match_new_mol,
                                                     indices_to_delete_mol)
                    match_new = frozenset(zip(match_new_ref, match_new_mol))
                    matches_new |= {m | match_new for m in matches}

    # don't prune disconnected MCS's just yet
    if iterate:
        return matches_new if len(matches_new) else {frozenset()}

    # prune all of the resulting MCS combinations so that they are connected
    matches_final = set()
    max_len = 0
    #matches_new = _eliminateSubsets(matches_new)
    for match_new in matches_new:
        if match_new == frozenset() or len(match_new) < max_len:
            continue
        match_new_ref, match_new_mol = zip(*match_new)
        indices_to_delete_ref = [x for x in range(ref.GetNumAtoms()) if
                                 x not in match_new_ref]
        indices_to_delete_mol = [x for x in range(mol.GetNumAtoms()) if
                                 x not in match_new_mol]
        frags_ref = _generateFragment(ref, indices_to_delete_ref, getAsFrags=True)
        frags_mol = _generateFragment(mol, indices_to_delete_mol, getAsFrags=True)

        # only keep MCS if it is connected
        for frag_ref in frags_ref:
            for frag_mol in frags_mol:
                match_new_frag = [(x, y) for x in frag_ref for y in frag_mol
                                  if (x, y) in match_new]
                n_atoms = len(match_new_frag)

                if n_atoms >= max_len:
                    if n_atoms == 0:
                        matches_final |= {frozenset()}
                    elif n_atoms == 1:
                        max_len = 1
                        matches_final |= {frozenset([x, y] )for x in range(len_ref) for y
                                          in range(len_mol)}
                    else:
                        fixed_mcss = _fixMCS(ref, mol, match_new_ref,
                                            match_new_mol, min_len=max_len-1,
                                             **kwargs)
                        max_len = max(max_len, len(max(fixed_mcss, key=len)))
                        matches_final |= {frozenset(x) for x in fixed_mcss}

    # remove all subsets and proceed to
    matches_final |= set(matches)
    matches_final = _eliminateSubsets(matches_final)

    # only keep largest unique MCS's
    matches = {frozenset(match) for match in matches_final if len(match) == max_len}
    matches = [list(match) for match in matches]

    return matches


def alignTwoMolecules(ref, mol, n_min=-1, mcs=None, mcs_parameters=None,
                      minimiser_parameters=None):
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
    if mcs_parameters is None:
        mcs_parameters = {}
    if minimiser_parameters is None:
        minimiser_parameters = {}

    if mcs is None:
        mcss = getMCSMap(ref, mol, **mcs_parameters)
        # only choose the maps that give maximum atom similarity
        scores = [getMatchingAtomScore(ref,mol, mcs) for mcs in mcss]
        scores_max = max(scores)
        mcss = [mcs for mcs, score in zip(mcss, scores) if score == scores_max]
    else:
        mcss = [mcs]

    # if we have multiple equivalent MCS's we pick the one with best RMSD
    mol_final, rmsd_final, mcs_final = None, None, None
    for mcs in mcss:
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

        mol, rmsd = minimiseRMSD(ref, mol, frozen_atoms=list(zip(*mcs))[1],
                                 **minimiser_parameters)

        if rmsd_final is None or rmsd < rmsd_final:
            mol_final, rmsd_final, mcs_final = mol, rmsd, mcs

    return mol_final, mcs_final


def minimiseRMSD(ref, mol, frozen_atoms=None, confId1=-1, confId2=-1,
                 algorithm=_minimize, **kwargs):
    if algorithm is None:
        return _copy.deepcopy(mol), getRMSD(ref, mol, confId1=confId1,
                                            confId2=confId2)
    if frozen_atoms is None:
        frozen_atoms = []

    mol = _copy.deepcopy(mol)
    mol_conf = mol.GetConformer(confId2)

    # find rotatable bonds outside of the frozen atom domain
    rbonds = _Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    rbond_indices = mol.GetSubstructMatches(rbonds)
    rbond_indices = [(i1, i2) for i1, i2 in rbond_indices
                     if not (i1 in frozen_atoms and i2 in frozen_atoms)]

    all_bonds = mol.GetBonds()
    all_bonds_indices = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx()) for x in all_bonds]

    dihedrals = []
    initial_values = []
    for i1, i2 in rbond_indices:
        neighbours_1 = [x for x in all_bonds_indices if i1 in x and i2 not in x]
        neighbours_2 = [x for x in all_bonds_indices if i2 in x and i1 not in x]

        # only rotate if the bond belongs to a dihedral
        if len(neighbours_1) and len(neighbours_2):
            dihedral = [list(set(neighbours_1[0]) - {i1})[0], i1,
                        i2, list(set(neighbours_2[0]) - {i2})[0]]
            mol_temp = _Chem.EditableMol(_copy.deepcopy(mol))
            mol_temp.RemoveBond(i1, i2)
            frags = _Chem.GetMolFrags(mol_temp.GetMol(),
                                      asMols=False, sanitizeFrags=False)

            for frag in frags:
                if i2 in frag and set(frozen_atoms).intersection(frag):
                    dihedral = list(reversed(dihedral))

            dihedrals += [dihedral]
            initial_values += [_Transforms.GetDihedralRad(mol_conf, *dihedral)]

    def rotateDihedral(*args):
        nonlocal dihedrals, mol, confId1

        mol_temp = _copy.deepcopy(mol)
        mol_conf_temp = mol_temp.GetConformer(confId1)
        for dihedral, arg in zip(dihedrals, list(args[0])):
            _Transforms.SetDihedralRad(mol_conf_temp, *dihedral, arg)

        return getRMSD(ref, mol_temp, confId1=confId1, confId2=confId2)

    # rotate dihedrals until RMSD minimisation
    if len(initial_values):
        result = algorithm(rotateDihedral, _np.asarray(initial_values),
                           **kwargs)
        optimal_angles, final_rmsd = list(result.x), result.fun
        for dihedral, optimal_angle in zip(dihedrals, optimal_angles):
            _Transforms.SetDihedralRad(mol_conf, *dihedral, optimal_angle)
        return mol, final_rmsd
    else:
        return mol, getRMSD(ref, mol, confId1=confId1, confId2=confId2)


def getRMSD(ref, mol, confId1=-1, confId2=-1):
    ref_conf = ref.GetConformer(confId1)
    mol_conf = mol.GetConformer(confId2)

    rmsd = 0
    for atom2 in mol.GetAtoms():
        rmsds = []
        for atom1 in ref.GetAtoms():
            x1 = ref_conf.GetAtomPosition(atom1.GetIdx())
            x2 = mol_conf.GetAtomPosition(atom2.GetIdx())
            rmsds += [(x1 - x2).LengthSq()]
        rmsd += min(rmsds)

    return rmsd


def getMatchingAtomScore(mol1, mol2, matches):
    matching_atoms = 0
    for i1, i2 in matches:
        n1 = mol1.GetAtomWithIdx(i1).GetAtomicNum()
        n2 = mol2.GetAtomWithIdx(i2).GetAtomicNum()
        matching_atoms += (n1 == n2)
    return matching_atoms


def _carbonify(mol):
    """
    Converts all atoms in a molecule to carbons.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input molecule.

    Returns
    -------
    mol_c : rdkit.Chem.rdchem.Mol
        The output molecule.
    """
    mol_c = _copy.deepcopy(mol)
    for atom in mol_c.GetAtoms():
        atom.SetAtomicNum(6)
    _Chem.AssignStereochemistry(mol_c, cleanIt=True, force=True)
    return mol_c


def _eliminateSubsets(sets):
    """Return a list of the elements of `sequence_of_sets`, removing all
    elements that are subsets of other elements.  Assumes that each
    element is a set or frozenset and that no element is repeated."""

    def is_power_of_two(n):
        """Returns True iff n is a power of two.  Assumes n > 0."""
        return (n & (n - 1)) == 0

    # The code below does not handle the case of a sequence containing
    # only the empty set, so let's just handle all easy cases now.
    if len(sets) <= 1:
        return list(sets)
    # We need an indexable sequence so that we can use a bitmap to
    # represent each set.
    if not isinstance(sets, _collections.Sequence):
        sets = list(sets)
    # For each element, construct the list of all sets containing that
    # element.
    sets_containing_element = {}
    for i, s in enumerate(sets):
        for element in s:
            try:
                sets_containing_element[element] |= 1 << i
            except KeyError:
                sets_containing_element[element] = 1 << i
    # For each set, if the intersection of all of the lists in which it is
    # contained has length != 1, this set can be eliminated.
    out = {s for s in sets
           if s and is_power_of_two(_functools.reduce(
               _operator.and_, (sets_containing_element[x] for x in s)))}
    return out


def _breakMismatchingBonds(ref, mol, match_ref_strict, match_mol_strict,
                           match_ref_loose, match_mol_loose):
    """
    Breaks all bonds bordering on an MCS if they are a subset of another MCS.
    Helpful at breaking rings in order to utilise the strict MCS algorithm.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    match_ref_strict : [int]
        The indices of the strictly matched reference.
    match_mol_strict : [int]
        The indices of the strictly matched molecule.
    match_ref_loose : [int]
        The indices of the loosely matched reference.
    match_mol_loose : [int]
        The indices of the loosely matched molecule.

    Returns
    -------

    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    successful : bool
        Whether any bonds were broken
    """
    swap = False
    # get all bonds
    bonds_ref = [frozenset([x.GetBeginAtomIdx(), x.GetEndAtomIdx()])
                 for x in ref.GetBonds()]
    bonds_mol = [frozenset([x.GetBeginAtomIdx(), x.GetEndAtomIdx()])
                 for x in mol.GetBonds()]

    # only keep bordering bonds that are part of the loose MCS

    bonds_ref_temp = {x for x in bonds_ref if len(x | set(match_ref_strict))
                 == len(match_ref_strict) + 1 }
    bonds_mol_temp = {x for x in bonds_mol if len(x | set(match_mol_strict))
                 == len(match_mol_strict) + 1 }

    bonds_ref = {x for x in bonds_ref if len(x | set(match_ref_strict))
                 == len(match_ref_strict) + 1 and x.issubset(match_ref_loose)}
    bonds_mol = {x for x in bonds_mol if len(x | set(match_mol_strict))
                 == len(match_mol_strict) + 1 and x.issubset(match_mol_loose)}

    # map the molecule indices to the reference indices and vice versa
    bonds_ref_new = set()
    for atom1, atom2 in bonds_ref:
        try:
            b1 = match_mol_loose[match_ref_loose.index(atom1)]
            b2 = match_mol_loose[match_ref_loose.index(atom2)]
        except (IndexError, ValueError):
            continue
        bonds_ref_new |= {frozenset([b1, b2])}

    bonds_mol_new = set()
    for atom1, atom2 in bonds_mol:
        try:
            b1 = match_ref_loose[match_mol_loose.index(atom1)]
            b2 = match_ref_loose[match_mol_loose.index(atom2)]
        except (IndexError, ValueError):
            continue
        bonds_mol_new |= {frozenset([b1, b2])}

    # swap the molecules if the reference is a subset of the molecule
    if bonds_ref == bonds_mol_new:
        return ref, mol, False
    elif bonds_ref.issubset(bonds_mol_new):
        bonds_to_break = bonds_mol - bonds_ref_new
        swap = not swap
        ref, mol = mol, ref
    elif bonds_mol_new.issubset(bonds_ref):
        bonds_to_break = bonds_ref - bonds_mol_new
    else:
        return ref, mol, False

    # break all bonds
    ref_copy = _Chem.EditableMol(_copy.deepcopy(ref))
    for bond in bonds_to_break:
        ref_copy.RemoveBond(*bond)
    ref_copy = ref_copy.GetMol()
    ref_copy.UpdatePropertyCache()

    if swap:
        return mol, ref_copy, True
    else:
        return ref_copy, mol, True


def _fixMCS(ref, mol, match_ref, match_mol, min_len=0, **kwargs):
    kwargs = {x: y for x, y in kwargs.items()
              if x not in ["ringMatchesRingOnly", "completeRingsOnly"]}

    # create a makeshift static memoisation variable
    if "memo" not in _fixMCS.__dict__:
        _fixMCS.memo = {}

    memo_key = (_Chem.MolToSmiles(ref), _Chem.MolToSmiles(mol),
                frozenset(match_ref), frozenset(match_mol))
    if memo_key in _fixMCS.memo.keys():
        return _fixMCS.memo[memo_key]

    # make sure that correct mapping is obtained between the two molecules given an MCS
    indices_to_delete_ref = [x for x in range(ref.GetNumAtoms())
                             if x not in match_ref]
    indices_to_delete_mol = [x for x in range(mol.GetNumAtoms())
                             if x not in match_mol]
    mcs_ref = _generateFragment(ref, indices_to_delete_ref)
    mcs_mol = _generateFragment(mol, indices_to_delete_mol)
    mcs_strict, _, _ = _matchAndReturnMatches([mcs_ref, mcs_mol],
                                             ringMatchesRingOnly=True,
                                             completeRingsOnly=True,
                                             **kwargs)

    if mcs_strict is None:
        return {frozenset()}
    matches_strict = {frozenset(zip(x, y))
                      for x in set(ref.GetSubstructMatches(mcs_strict))
                      for y in set(mol.GetSubstructMatches(mcs_strict))}

    mcs_loose, _, _ = _matchAndReturnMatches([mcs_ref, mcs_mol],
                                            ringMatchesRingOnly=False,
                                            completeRingsOnly=False,
                                            **kwargs)
    matches_loose = {frozenset(zip(x, y))
                     for x in set(ref.GetSubstructMatches(mcs_loose))
                     for y in set(mol.GetSubstructMatches(mcs_loose))}
    matches = _copy.deepcopy(matches_strict)

    # break the recursion early if the matches are not big enough
    max_len_strict = len(max(matches_strict, key=len))
    max_len_loose = len(max(matches_loose, key=len))
    if max_len_strict <= min_len:
        return {frozenset()}

    # recursively break bonds to maximise strict MCS
    if max_len_strict < max_len_loose:
        for match_strict in matches_strict:
            for match_loose in matches_loose:
                if match_strict != match_loose:
                    ref_broken, mol_broken, flag = _breakMismatchingBonds(
                        ref, mol, *zip(*match_strict), *zip(*match_loose))
                    #print(_Chem.MolToSmiles(ref), _Chem.MolToSmiles(mol))
                    #(_Chem.MolToSmiles(mcs_ref), _Chem.MolToSmiles(mcs_mol))
                    if flag:
                        max_len = max([len(x) for x in matches])
                        matches |= _fixMCS(ref_broken, mol_broken,
                                           match_ref, match_mol,
                                           min_len=max_len, **kwargs)
                        #print(_Chem.MolToSmiles(ref_broken), _Chem.MolToSmiles(mol_broken))


    # check if both atoms are chiral and flag the atoms if this is the case
    matches_final = set()
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
            ref_c = _carbonify(ref)
            mol_c = _carbonify(mol)
            chiral_ref_c = dict(_Chem.FindMolChiralCenters(ref_c))
            chiral_mol_c = dict(_Chem.FindMolChiralCenters(mol_c))

            indices_chiral_final = []
            for idx in indices_chiral:
                i_ref, i_mol = match[idx]
                if chiral_ref_c[i_ref] != chiral_mol_c[i_mol]:
                    indices_chiral_final += [idx]

            if indices_chiral_final:
                match = _fixChiralIndices(indices_chiral_final, mcs_strict,
                                          match)

        matches_final |= {match}

    _fixMCS.memo[memo_key] = matches_final

    return matches_final

def _matchAndReturnMatches(*args, **kwargs):
    # returns the MCS of two molecules as given by RDKit
    mcs_string = _MCS.FindMCS(*args, **kwargs).smarts
    if mcs_string is None:
        return None, [], []
    mcs = _Chem.MolFromSmarts(mcs_string)
    # fixes some issues with RDKit concerning SMARTS
    mcs = _Chem.MolFromSmarts(_Chem.MolToSmiles(mcs, canonical=False,
                                                allBondsExplicit=True,
                                                allHsExplicit=True))
    return tuple([mcs] + [x.GetSubstructMatches(mcs) for x in args[0]])


def _generateFragment(mol, indices_to_delete, getAsFrags=False):
    # generates a molecule from another given indices to delete
    mol_frag_edit = _Chem.EditableMol(_copy.deepcopy(mol))
    for idx in reversed(sorted(set(indices_to_delete))):
        mol_frag_edit.RemoveAtom(idx)
    mol_frag = mol_frag_edit.GetMol()
    if getAsFrags:
        frags = _rdmolops.GetMolFrags(mol_frag, asMols=False, sanitizeFrags=False)
        return [_transformIndices(f, indices_to_delete) for f in frags]
    else:
        mol_frag.UpdatePropertyCache()
        return mol_frag


def _fixChiralIndices(indices, mcs, matches):
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
            fragment_new = {frozenset(i for i in _transformIndices(x, indices_to_delete) + [index])
                            for x in _rdmolops.GetMolFrags(mcs_edit.GetMol(), asMols=False, sanitizeFrags=False)}
            fragments_new |= fragment_new
        fragments = fragments_new
    return {x for i, x in enumerate(matches) if i in max(fragments, key=len)}


def _transformIndices(current_indices, deleted_indices):
    """
    Generates the old indices of the molecule given the new and the deleted
    indices.

    Parameters
    ----------
    current_indices : [int]
        Current indices that are to be transformed.
    deleted_indices : [int]
        Deleted indices from the previous molecule.

    Returns
    -------
    transformed_indices : [int]
        Transformed indices.
    """
    if not current_indices:
        return []
    if not deleted_indices:
        return current_indices
    max_len = max(deleted_indices) + max(current_indices) + 2
    remaining_indices = [x for x in range(max_len) if x not in deleted_indices]
    return [remaining_indices[x] for x in current_indices]