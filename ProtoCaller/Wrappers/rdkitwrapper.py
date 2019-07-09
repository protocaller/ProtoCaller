import itertools as _it
import copy as _copy
import os as _os
import sys as _sys

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
        Keyword arguments to be supplied to the more specialised RDKit
        functions.

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
        raise TypeError("Unrecognised input extension for RDKit: {}".format(
            extension))

    if mol is None:
        raise ValueError("Invalid file format for file: {}".format(filename))
    return mol


@_stdio.stdout_stderr()
def openInChIAsRdkit(inchi_str, removeHs=True, **kwargs):
    """
    Converts an InChI string into an RDKit Mol object. This function is a
    simple wrapper of rdkit.Chem.MolFromInchi.

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
    Converts a SMILES string into an RDKit Mol object. This function is a
    simple wrapper of rdkit.Chem.MolFromSmiles.

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
    A general wrapper which can convert a variety of representations for a
    molecule into an rdkit.Chem.rdchem.Mol object.

    Parameters
    ----------
    val : str
        Input value - SMILES, InChI strings or a filename.
    minimise : bool or None
        Whether to perform a GAFF minimisation using OpenBabel. None means
        minimisation for molecules initialised from
        strings and no minimisation for molecules initialised from files.
    kwargs
        Keyword arguments to be passed to the more specialsied RDKit functions.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The input string opened as an RDKit Mol object.
    """
    if isinstance(val, str) and len(val.split(".")) == 1:
        if minimise is None:
            minimise = True
        flist = [openSmilesAsRdkit, openInChIAsRdkit]
        for i, f in enumerate(flist):
            try:
                mol = f(val, **kwargs)
                _AllChem.EmbedMolecule(mol, useRandomCoords=True)
                break
            except:
                if i == len(flist) - 1:
                    raise ValueError("String not recognised as a valid SMILES "
                                     "or InChI input")
        if minimise:
            # make sure we preserve the chirality of the input string and
            # minimise with obminimize / GAFF
            with _fileio.Dir("Temp", temp=True):
                smiles = _Chem.MolToSmiles(mol)
                with open("molecule.smi", "w") as f:
                    f.write(smiles)
                _babel.babelTransform("molecule.smi", output_extension="sdf",
                                      generate_3D_coords=True)
                _runexternal.runExternal("obminimize -sd -c 1e-6 -n 10000 -ff "
                                         "GAFF -osdf molecule.sdf > "
                                         "minimised.sdf")
                mol = openFileAsRdkit("minimised.sdf")
    else:
        if minimise is None:
            minimise = False
        try:
            mol = openFileAsRdkit(val, **kwargs)
        except:
            raise ValueError("File is not in a valid format")
        # minimise the molecule using obminimize / GAFF
        if minimise:
            with _fileio.Dir("Temp", temp=True):
                saveFromRdkit(mol, "molecule.sdf")
                _runexternal.runExternal("obminimize -sd -c 1e-6 -n 10000 -ff "
                                         "GAFF -osdf molecule.sdf > "
                                         "minimised.sdf")
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
            _babel.babelTransform(tempfilename, output_extension=extension,
                                  pH=None)
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
        new_position = mol_conf.GetAtomPosition(i) + _Geom.Point3D(*vector)
        mol_conf.SetAtomPosition(i, new_position)
    return mol


def getMCSMap(ref, mol, atomCompare="any", bondCompare="any", **kwargs):
    """
    Generates the Maximum Common Substructure (MCS) mapping between two
    molecules. This algorithm calls the getFixedMCS() function which improves
    on RDKit's default MCS algorithm. It does so recursively until it is
    certain that the MCS has definitely been found.

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
        Additional keyword arguments passed on to getFixedMCS() or to
        rdkit.Chem.MCS.FindMCS().

    Returns
    -------
    mcs : [[tuple]]]
        A list of lists of tuples corresponding to the atom index matches
        between the reference and the other molecule.
    """
    kwargs_default = {
        "maximize": "atoms",
        "timeout": 60,
    }

    kwargs = {
        **kwargs_default,
        **kwargs,
        'atomCompare': atomCompare,
        'bondCompare': bondCompare,
        'completeRingsOnly': True,
        'ringMatchesRingOnly': True,
    }

    matches = getFixedMCS(ref, mol, [i for i in range(ref.GetNumAtoms())],
                          [i for i in range(mol.GetNumAtoms())], **kwargs)
    if matches == {frozenset()}:
        return []
    max_len_final = 0
    while True:
        # get the first of or all other smaller initial substructures
        master_set = set().union(*matches)
        max_len = max([len(x) for x in matches])
        mapped_atoms_ref, mapped_atoms_mol = list(zip(*master_set))
        unmapped_atoms_ref = ref.GetNumAtoms() - len(mapped_atoms_ref)
        unmapped_atoms_mol = mol.GetNumAtoms() - len(mapped_atoms_mol)

        # only keep iterating if it is theoretically possible to get a bigger
        # structure than what we have
        if max_len <= max_len_final or max_len > unmapped_atoms_mol or \
                max_len > unmapped_atoms_ref:
            break
        else:
            frags_ref = _generateFragment(ref, mapped_atoms_ref,
                                          getAsFrags=True)
            frags_ref = [x for x in frags_ref if len(x) >= max_len]
            frags_mol = _generateFragment(mol, mapped_atoms_mol,
                                          getAsFrags=True)
            frags_mol = [x for x in frags_mol if len(x) >= max_len]
            for frag_ref in frags_ref:
                for frag_mol in frags_mol:
                    matches |= getFixedMCS(ref, mol, frag_ref, frag_mol,
                                           **kwargs)
            max_len_final = max_len

    # only keep largest unique MCS's
    matches = _onlyKeepLongest(matches)[0]
    matches = [list(match) for match in matches]

    return matches


def getFixedMCS(ref, mol, match_ref, match_mol, break_recursively=True,
                valid_mcs=False, two_way_matching=True, **kwargs):
    """
    A helper function to getMCSMap which expands on RDKit's strict MCS
    algorithm. The first addition is recursive bond breaking to find a
    connected MCS where aliphatic atoms can be mapped to ring atoms. The second
    addition is the detection of atoms with different stereochemistry and
    finding the biggest fragment that doesn't contain them. This algorithm
    still doesn't support the functionality of mapping two very similar atoms
    of different stereochemistry. This might be added in the future.

    Parameters
    ----------
    ref : rdkit.Chem.Mol
        The reference molecule.
    mol : rdkit.Chem.Mol
        The molecule to be aligned.
    match_ref : list
        A list of indices corresponding to the MCS for the reference molecule.
    match_mol : list
        A list of indices corresponding to the MCS for the molecule to be
        aligned.
    break_recursively : bool
        Whether to use the recursive bond breaking algorithm to improve on
        RDKit's functionality
    valid_mcs : bool
        Whether the input MCS is ordered, unique and valid.
    two_way_matching : bool
        This algorithm matches certain aliphatic mol atoms to certain ring
        ref atoms by default. This parameter determines if the opposite is also
        true. A value of True is good for matching Ligand B to ligand A and
        a value of False is good for matching Ligand A to the binding pocket
        ligand.
    kwargs :
        Keyword arguments to be supplied to _matchAndReturnMatches()

    Returns
    -------
    matches : {frozenset([tuple])}
        A set of frozensets corresponding to the largest unique MCS's given the
        input.
    """
    kwargs = {x: y for x, y in kwargs.items()}
    kwargs["ringMatchesRingOnly"] = True
    kwargs["completeRingsOnly"] = True

    # make sure that correct mapping is obtained between the two molecules
    # given an MCS
    if not valid_mcs:
        indices_to_delete_ref = [x for x in range(ref.GetNumAtoms())
                                 if x not in match_ref]
        indices_to_delete_mol = [x for x in range(mol.GetNumAtoms())
                                 if x not in match_mol]
        mcs_ref = _generateFragment(ref, indices_to_delete_ref)
        mcs_mol = _generateFragment(mol, indices_to_delete_mol)
        mcs, _, _ = _matchAndReturnMatches([mcs_ref, mcs_mol], **kwargs)

        if mcs is None:
            return {frozenset()}
        matches = {frozenset(zip(x, y))
                   for x in set(ref.GetSubstructMatches(mcs))
                   for y in set(mol.GetSubstructMatches(mcs))}
    else:
        matches = {frozenset(zip(match_ref, match_mol))}

    if break_recursively:
        # recursively break bonds to maximise MCS
        matches_new = set()
        matches_rec_prev = {x: {x} for x in matches}
        matches = list(sorted(matches, key=len))
        explored_sets = set()
        memo_dict = {(_Chem.MolToSmiles(ref, allHsExplicit=True),
                      _Chem.MolToSmiles(mol, allHsExplicit=True)) : matches}
        for match in matches:
            if match in explored_sets:
                continue
            while True:
                matches_rec = set()
                for submatch in matches_rec_prev[match]:
                    refs_broken, mols_broken = \
                        _breakMismatchingBonds(ref, mol, *list(zip(*submatch)))
                    pairs_broken = [(x, mol) for x in refs_broken]
                    if two_way_matching:
                        pairs_broken += [(ref, x) for x in mols_broken]

                    # get recursive matches
                    for ref_broken, mol_broken in pairs_broken:
                        key_ref = _Chem.MolToSmiles(ref_broken,
                                                    allHsExplicit=True)
                        key_mol = _Chem.MolToSmiles(mol_broken,
                                                    allHsExplicit=True)

                        # memoise the MCSs
                        if (key_ref, key_mol) in memo_dict.keys():
                            matches_rec |= {x for x in memo_dict[(key_ref, key_mol)]
                                            if _haveCommonElements(x, submatch)}
                            continue

                        mcs_broken, matches_ref_broken, matches_mol_broken = \
                            _matchAndReturnMatches([ref_broken, mol_broken],
                                                   **kwargs)
                        matches_broken = {frozenset(zip(x, y))
                                          for x in matches_ref_broken
                                          for y in matches_mol_broken}
                        # only keep the matches that are connected to the
                        # original one
                        memo_dict[(key_ref, key_mol)] = matches_broken
                        matches_rec |= {x for x in matches_broken
                                        if _haveCommonElements(x, submatch)}

                # combine the new MCS's in an optimal way
                matches_rec_max, max_len_rec = \
                    _optimalMergedSets(*matches_rec_prev[match], *matches_rec,
                                       seed=None)
                # don't duplicate computational effort for already explored
                # sets
                matches_rec_max -= explored_sets
                explored_sets |= matches_rec | matches_rec_max

                # break the recursion if the new MCS is not larger
                if matches_rec_max and max_len_rec > len(max(
                        matches_rec_prev[match], key=len)):
                    matches_rec_prev[match] = matches_rec_max
                else:
                    matches_new |= matches_rec_prev[match]
                    break

        matches = matches_new

    # break mismatching neighbours next to a double / amide / ester bond
    matches_new = set()
    for match in matches:
        results = _breakEZBonds(ref, mol, *list(zip(*match)),
                                ignore_single_bonds=not two_way_matching)

        if not len(results[0]):
            # matching is already correct
            matches_new |= {match}
            continue

        match_new = set()
        for ref_broken, mol_broken, mcs_broken, del_ref, del_mol in zip(*results):
            mcs_broken = set(mcs_broken)
            _, matches_ref_broken, matches_mol_broken = \
                _matchAndReturnMatches([ref_broken, mol_broken],
                                       **kwargs)
            matches_broken = {frozenset(zip(_transformIndices(x, del_ref),
                                            _transformIndices(y, del_mol)))
                              for x in matches_ref_broken
                              for y in matches_mol_broken}
            # only keep the matches that are connected to the
            # original one
            match_new |= {x for x in matches_broken
                          if _haveCommonElements(x, mcs_broken)}

            matches_new |= _optimalMergedSets(*match_new, seed=mcs_broken)[0]
    matches = matches_new

    # here we deal with mismatching atoms of different chirality
    matches_final = set()
    # take care of R/S changing with different atom types
    _Chem.AssignStereochemistry(ref, cleanIt=True, force=True)
    _Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    ref_c = _carbonify(ref)
    mol_c = _carbonify(mol)
    chiral_ref_c = dict(_Chem.FindMolChiralCenters(ref_c))
    chiral_mol_c = dict(_Chem.FindMolChiralCenters(mol_c))
    EZ_ref_c = getEZStereochemistry(ref)
    EZ_mol_c = getEZStereochemistry(mol)

    # check if both atoms are chiral and flag the atoms if this is the case
    for match in matches:
        if match == frozenset():
            continue
        match_rev_dict = {mol_atom: ref_atom for ref_atom, mol_atom in match}

        # only deal with one of the molecules
        chiral_ref_indices = {x[0] for x in match
                              if x[0] in chiral_ref_c.keys()
                              and x[1] in chiral_mol_c.keys()
                              and chiral_ref_c[x[0]] != chiral_mol_c[x[1]]}
        mismatching_bonds = []

        if not two_way_matching:
            for bond_mol in EZ_mol_c.keys():
                if set(bond_mol).issubset(set(list(match_rev_dict.keys()))):
                    bond_ref = [match_rev_dict[x] for x in bond_mol]
                    if frozenset(bond_ref) not in EZ_ref_c.keys():
                        mismatching_bonds += [frozenset(bond_ref)]
            chiral_ref_indices |= set().union(*mismatching_bonds)

        # delete invalid atoms
        frags_ref_total = {frozenset(list(zip(*match))[0])}
        for i_ref in chiral_ref_indices:
            frags_ref = set()
            for frag_prev in frags_ref_total:
                # ignore if the fragment doesn't contain the chiral atom
                if i_ref not in frag_prev:
                    frags_ref |= {frag_prev}
                    continue
                # only include indices that are not a result of ring breaking
                indices_to_delete = {i for i in range(ref.GetNumAtoms())} - \
                                    frag_prev | {i_ref}
                frags_ref_temp = _generateFragment(ref, indices_to_delete,
                                                   getAsFrags=True)
                frags_ref |= {frozenset(set(x) | {i_ref})
                              for x in frags_ref_temp if i_ref not in x}
            frags_ref_total = frags_ref

        frags_ref_total = {x.union(*[y for y in mismatching_bonds if y & x])
                           for x in frags_ref_total}
        frags_ref_total = _onlyKeepLongest(frags_ref_total)[0]
        for frag_ref in frags_ref_total:
            matches_final |= {frozenset([x for x in match
                                         if x[0] in frag_ref])}

    matches_final, max_len_final = _onlyKeepLongest(matches_final)
    if not matches_final:
        return {frozenset()}

    return matches_final


def alignTwoMolecules(ref, mol, n_min=-1, two_way_matching=True, mcs=None,
                      mcs_parameters=None, minimiser_parameters=None,
                      minimise_score=True):
    """
    Aligns two molecules based on an input MCS. The algorithm uses atom
    freezing of the common core and force field minimisation of the rest.
    Additional minimisation using minimiseAlignmentScore() is then performed.

    If there is a choice between several equally long MCS's the ones with
    the highest atom-atom matches are first selected (e.g. C->C mapping trumps
    C->H mapping). After that the first structure with the lowest MSD is
    selected.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    two_way_matching : bool
        Whether to treat ref and mol equally in terms of matching.
    n_min : int
        Minimum number of force field minimisation iterations. -1 is no limit.
    mcs : [tuple] or None
        The maximum common substucture. None means the one generated from
        getMCSMap.
    mcs_parameters : dict
        A dictionary of the parameters to be passed on to getMCSMap().
    minimiser_parameters : dict
        A dictionary of the parameters to be passed on to
        minimiseAlignmentScore().

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The aligned molecule.
    mcs : [tuple]
        The maximum common substructure.
    """
    if mcs_parameters is None:
        mcs_parameters = {}
    mcs_parameters = {
        **mcs_parameters,
        "two_way_matching": two_way_matching,
    }
    if minimiser_parameters is None:
        minimiser_parameters = {}

    if mcs is None:
        mcss = getMCSMap(ref, mol, **mcs_parameters)
        # only choose the maps that give maximum atom similarity
        scores = [getMatchingAtomScore(ref, mol, mcs) for mcs in mcss]
        scores_max = max(scores)
        mcss = [mcs for mcs, score in zip(mcss, scores) if score == scores_max]
    else:
        mcss = [mcs]

    # if we have multiple equivalent MCS's we pick the one with best score
    mol_final, score_final, mcs_final = None, None, None
    for mcs in mcss:
        ref_conf = ref.GetConformer(-1)

        # get the non MCS dihedrals and store their initial values
        dihedrals_ini = _nonMCSDihedrals(mol, mcs=mcs)

        # MMFF gives better conformations but can occasionally fail
        ff_to_use = "mmff"
        while True:
            # translate the molecule randomly to prevent minimisation failure
            try:
                vec = _np.random.uniform(-1, 1, size=3)
                mol = translateMolecule(mol, tuple(vec))
                mol_conf = mol.GetConformer(-1)

                if ff_to_use == "mmff":
                    ff = _FF.MMFFGetMoleculeForceField(
                        mol, _FF.MMFFGetMoleculeProperties(mol), confId=0)
                else:
                    ff = _FF.UFFGetMoleculeForceField(mol)

                for i_ref, i_mol in mcs:
                    mol_conf.SetAtomPosition(i_mol,
                                             ref_conf.GetAtomPosition(i_ref))
                    ff.AddFixedPoint(i_mol)

                if mol_conf.GetNumAtoms() != len(mcs):
                    ff.Initialize()
                    more = ff.Minimize()
                    while more and n_min:
                        more = ff.Minimize()
                        # n_min = -1 means infinite minimisation
                        if n_min != -1:
                            n_min -= 1

                # restore the dihedral angles to their initial values
                for dihedral, val in dihedrals_ini.items():
                    _Transforms.SetDihedralRad(mol_conf, *dihedral, val)
                break
            except:
                # sometimes MMFF fails and in this case we use UFF
                ff_to_use = "uff"
                continue

        if not minimise_score:
            return mol, mcs

        mol_new, score = minimiseAlignmentScore(ref, mol, mcs=mcs,
                                                **minimiser_parameters)

        if score_final is None or score < score_final:
            mol_final = mol_new
            score_final = score
            mcs_final = mcs

    return mol_final, mcs_final


def getAlignmentScore(ref, mol, mcs=None, confId1=-1, confId2=-1):
    """
    Returns the alignment score between two molecules. The way this is done is
    by calculating all possible distances between atom i of ref and atom j of
    mol and computing the sum of their squares. All atoms from mol that are
    within 1 Angstrom of atom j contribute inverse square distances to the
    score in order to prevent unfavourable clashes.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    mcs : [(tuple)] or None
        The maximum common substructure of the two molecules.
    confId1 : int
        The conformer number of ref.
    confId2 : int
        The conformer number of mol.

    Returns
    -------
    alignment_score : float
        The alignment score of the two conformers.
    """
    if not mcs is None:
        frozen_atoms_ref, frozen_atoms_mol = list(zip(*mcs))
    else:
        frozen_atoms_ref, frozen_atoms_mol = [], []

    ref_conf = ref.GetConformer(confId1)
    mol_conf = mol.GetConformer(confId2)

    alignment_score = 0
    mol_indices = [x for x in range(mol.GetNumAtoms())
                   if x not in frozen_atoms_mol]
    ref_indices = [x for x in range(ref.GetNumAtoms())
                   if x not in frozen_atoms_ref]
    for mol_idx in mol_indices:
        for ref_idx in ref_indices:
            x1 = ref_conf.GetAtomPosition(ref_idx)
            x2 = mol_conf.GetAtomPosition(mol_idx)
            alignment_score += (x1 - x2).LengthSq()

        for mol_idx2 in range(ref.GetNumAtoms()):
            x1 = mol_conf.GetAtomPosition(mol_idx)
            x2 = mol_conf.GetAtomPosition(mol_idx2)
            len_sq = (x1 - x2).LengthSq()
            if len_sq > 1:
                continue
            elif len_sq < 2 / _sys.maxsize:
                return _sys.maxsize
            else:
                alignment_score += 2 / len_sq**6

    return alignment_score


def minimiseAlignmentScore(ref, mol, mcs=None, confId1=-1, confId2=-1,
                           minimisation_algorithm=_minimize,
                           scoring_algorithm=getAlignmentScore, **kwargs):
    """
    Rotates all of the rotatable dihedrals outside of the MCS  until a balance
    between good alignment to the reference and minimal clashing with other
    atoms in the same molecule is achieved.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    mcs : [(tuple)] or None
        The maximum common substructure of the two molecules.
    confId1 : int
        The conformer number of ref.
    confId2 : int
        The conformer number of mol.
    minimisation_algorithm : function
        A scipy algorithm to be used for minimisation. Default is regular
        local minimisation using scipy.optimize.minimize. Custom algorithms
        can also be used as long as they obey the same format as the ones in
        scipy.optimize.
    kwargs:
        Keyword arguments to be supplied to the minimiser algorithm.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        The aligned molecule.
    alignment_score : float
        The alignment score after alignment obtained from getAlignmentScore().
    """
    if mcs is None:
        mcs = []
        frozen_atoms = []
    else:
        frozen_atoms = list(zip(*mcs))[1]

    if minimisation_algorithm is None:
        return _copy.deepcopy(mol), scoring_algorithm(ref, mol, mcs=mcs,
                                                      confId1=confId1,
                                                      confId2=confId2)

    mol = _copy.deepcopy(mol)
    mol_conf = mol.GetConformer(confId2)

    all_bonds = mol.GetBonds()
    all_bonds_indices = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx()) for x in
                         all_bonds]

    # a dictionary which shows how many frozen atoms each atom is bonded to
    n_frozen_neighbours = {}
    unfrozen_neighbours = {}
    for atom in range(mol.GetNumAtoms()):
        atom_neighbours = set().union(*[set(x) for x in all_bonds_indices
                                        if atom in x]) - {atom}
        unfrozen_neighbours[atom] = atom_neighbours - set(frozen_atoms)
        n_frozen_neighbours[atom] = len(atom_neighbours) - \
                                    len(atom_neighbours - set(frozen_atoms))

    # find rotatable bonds at the edge of and outside of the frozen atom domain
    amide = "[NX3]!@[C]=[OX1]"
    ester = "[OX2]!@[C]=[OX1]"
    smarts = "[!$({0})&!$({1})&!D1&!$(*#*)]-&!@[!$({0})&!$({1})&!D1&!$(*#*)]"
    rbonds = _Chem.MolFromSmarts(smarts.format(amide, ester))
    rbond_indices = mol.GetSubstructMatches(rbonds)
    rbond_indices = [{i1, i2} for i1, i2 in rbond_indices
                     if n_frozen_neighbours[i1] <= 1
                     or n_frozen_neighbours[i2] <= 1]

    # get non MCS dihedrals corresponding to non-rotatable bonds
    dihedrals = _nonMCSDihedrals(mol, mcs, confId=confId2)
    dihedrals = {k: v for k, v in dihedrals.items()
                 if set(k[1:3]) in rbond_indices}

    def rotateDihedral(*args):
        nonlocal dihedrals, mol, confId1

        mol_temp = _copy.deepcopy(mol)
        mol_conf_temp = mol_temp.GetConformer(confId1)
        for dihedral, arg in zip(dihedrals.keys(), list(args[0])):
            _Transforms.SetDihedralRad(mol_conf_temp, *dihedral, arg)

        return scoring_algorithm(ref, mol_temp, mcs=mcs, confId1=confId1,
                                 confId2=confId2)

    # rotate dihedrals until MSD minimisation
    if len(dihedrals):
        ini_vals = _np.array(list(dihedrals.values()))
        result = minimisation_algorithm(rotateDihedral, ini_vals, **kwargs)
        optimal_angles, final_alignment_score = list(result.x), result.fun
        for dihedral, optimal_angle in zip(dihedrals.keys(), optimal_angles):
            _Transforms.SetDihedralRad(mol_conf, *dihedral, optimal_angle)
        return mol, final_alignment_score
    else:
        return mol, scoring_algorithm(ref, mol, mcs=mcs, confId1=confId1,
                                      confId2=confId2)


def getEZStereochemistry(mol, carbonify=True, confId=-1, extra_bonds=None):
    """
    Determines the stereochemistry of double bonds, esters and amides. Note
    that the corresponding E/Z labels might not correspond to the IUPAC E/Z
    labels. None corresponds to a symmetric bond.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input molecule.
    carbonify : bool
        Whether to ignore atom type when determining the stereochemistry.
    confId : int
        The conformer ID to be used for stereochemistry determination.
    extra_bonds : [(int, int)]
        Extra bonds to be treated as double.

    Returns
    -------
    bonds : dict(frozenset(int, int), str)
        The relevant double bonds and their labels.
    """
    if carbonify:
        mol_c = _carbonify(mol)
    else:
        mol_c = _copy.deepcopy(mol)

    if extra_bonds is None:
        extra_bonds = []

    all_bonds = [x for x in mol.GetBonds()]
    all_bonds_c = [x for x in mol_c.GetBonds()]
    all_bond_indices = [{x.GetBeginAtomIdx(), x.GetEndAtomIdx()}
                        for x in all_bonds]

    amide = "[N;X2,X3;R0]!@[C;R0]=[OX1]"
    ester = "[OX2;R0]!@[C;R0]=[OX1]"
    double = "[$([*;X3;R0])&!$(*=*!-*)]=&!@[$([*;X2,X3;R0])&!$(*=*!-*)]"

    matches_amides = list(mol.GetSubstructMatches(_Chem.MolFromSmarts(amide)))
    matches_esters = list(mol.GetSubstructMatches(_Chem.MolFromSmarts(ester)))
    matches_double = list(mol.GetSubstructMatches(_Chem.MolFromSmarts(double)))
    matches_all = {frozenset(x[:2]) for x in matches_amides + matches_esters +
                   matches_double + extra_bonds}

    bond_mapper = {
        _Chem.BondStereo.STEREOE: "E",
        _Chem.BondStereo.STEREOZ: "Z",
        _Chem.BondStereo.STEREONONE: None,
    }

    bonds = {}

    # generate E/Z information
    for bond in matches_all:
        idx = all_bond_indices.index(bond)
        all_bonds_c[idx].SetBondType(_Chem.BondType.DOUBLE)
        _Chem.DetectBondStereochemistry(mol_c, confId=confId)
        _Chem.AssignStereochemistry(mol_c, cleanIt=True, force=True)
        bonds[bond] = bond_mapper[mol_c.GetBondWithIdx(idx).GetStereo()]
        all_bonds_c[idx].SetBondType(_Chem.BondType.SINGLE)

    return bonds


def getMatchingAtomScore(mol1, mol2, matches):
    """
    Returns the matching atom score between two molecules. This is simply
    the number of atoms that are mapped onto the same element from the other
    molecule.

    Parameters
    ----------
    mol1 : rdkit.Chem.rdchem.Mol
        The first molecule.
    mol2 : rdkit.Chem.rdchem.Mol
        The second molecule.
    matches : [tuple]
        The maximum common substructure matches.

    Returns
    -------
    matching_atoms : float
        The matching atom score of the two molecules.
    """
    matching_atoms = 0
    for i1, i2 in matches:
        n1 = mol1.GetAtomWithIdx(i1).GetAtomicNum()
        n2 = mol2.GetAtomWithIdx(i2).GetAtomicNum()
        matching_atoms += (n1 == n2)
    return matching_atoms


def _areCompatibleSets(set1, set2):
    """
    Determines whether two sets of tuples are compatible.

    Parameters
    ----------
    set1 : set
        The first input set.
    set2 : set
        The second input set.

    Returns
    -------
    are_compatible : bool
        Whether the sets have common tuple elements.
    """
    set12 = set1.union(set2)
    tuple1, tuple2 = zip(*set12)
    if len(set(tuple1)) == len(set(tuple2)) == len(set12):
        return True
    else:
        return False


def _breakEZBonds(ref, mol, match_ref, match_mol, ignore_single_bonds=False):
    """
    Detects double bonds with mismatching E/Z stereochemistry and breaks all
    adjacent bonds one by one so that the MCS algorithm is forced to obey the
    stereochemistry.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    match_ref : [int]
        The indices of the matched reference.
    match_mol : [int]
        The indices of the matched molecule.
    ignore_single_bonds : bool
        Whether to ignore mappings of double bonds onto single bonds.

    Returns
    -------
    refs : [rdkit.Chem.rdchem.Mol]
        The reference molecules with broken bonds.
    mols : [rdkit.Chem.rdchem.Mol]
        The target molecules with broken bonds.
    mcss : [[tuple]]
        The new pruned MCSs - one per ref-mol pair
    ids_ref : [[int]]
        The deleted atoms from the original ref resulting in each fragment.
    ids_mol : [[int]]
        The deleted atoms from the original mol resulting in each fragment.
    """
    match_rev_map = {x: y for x, y in zip(match_mol, match_ref)}

    # get E/Z information of the original molecule
    EZ_mol = getEZStereochemistry(mol)

    # get MCS as a fragment and obtain its E/Z information
    indices_to_delete_ref = [x for x in range(ref.GetNumAtoms())
                             if x not in match_ref]
    indices_to_delete_mol = [x for x in range(mol.GetNumAtoms())
                             if x not in match_mol]
    mcs_ref = _generateFragment(ref, indices_to_delete_ref)
    mcs_mol = _generateFragment(mol, indices_to_delete_mol)
    # initialise the ring info, otherwise an error is thrown
    _Chem.GetSymmSSSR(mcs_ref)
    _Chem.GetSymmSSSR(mcs_mol)

    EZ_ref_mcs = getEZStereochemistry(mcs_ref)
    EZ_ref_mcs = {frozenset(_transformIndices(k, indices_to_delete_ref)): v
                  for k, v in EZ_ref_mcs.items()}
    EZ_mol_mcs = getEZStereochemistry(mcs_mol)
    EZ_mol_mcs = {frozenset(_transformIndices(k, indices_to_delete_mol)): v
                  for k, v in EZ_mol_mcs.items()}
    # get rid of symmetric bonds with asymmetric mappings
    EZ_mol_mcs = {key: val for key, val in EZ_mol_mcs.items()
                  if frozenset(key) in EZ_mol.keys()
                  and EZ_mol[key] is not None}

    # get all bonds
    bonds_ref = {(x.GetBeginAtomIdx(), x.GetEndAtomIdx()): x
                 for x in ref.GetBonds()}
    bonds_mol = {(x.GetBeginAtomIdx(), x.GetEndAtomIdx()): x
                 for x in mol.GetBonds()}
    bonds_mcs_ref = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx())
                     for x in mcs_ref.GetBonds()]
    bonds_mcs_mol = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx())
                     for x in mcs_mol.GetBonds()]
    # transform fragment indices to master indices
    bonds_mcs_ref = [tuple(_transformIndices(x, indices_to_delete_ref))
                     for x in bonds_mcs_ref]
    bonds_mcs_mol = [tuple(_transformIndices(x, indices_to_delete_mol))
                     for x in bonds_mcs_mol]

    # bonds that we break in the first instance (tetrahedral bonds) that are
    # not in the MCS
    bonds_to_break_ref = []
    # pairs of bonds that we break independently of one another
    bonds_to_break_pairs = {}

    # populate the above empty containers using bond information from original
    # molecules (will need to transform the indices back later)
    for bond_mol, chirality_mol in EZ_mol_mcs.items():
        # make sure we have an ordered bond
        bond_mol = [key for key, val in bonds_mol.items()
                    if set(key) == set(bond_mol)][0]
        # only consider bonds that are inside the MCS
        if set(bond_mol).issubset(match_mol):
            bond_ref = tuple([match_ref[match_mol.index(x)] for x in bond_mol])
            bonds_to_break_pairs[(bond_ref, bond_mol)] = []

            # get all adjacent bonds for atom 0 and atom 1 for both mols
            neighbours = {"ref": {}, "mol": {}}
            neighbours["mol"][0] = [k for k, v in bonds_mol.items()
                                    if bond_mol[0] in k
                                    and bond_mol[1] not in k
                                    and set(k).issubset(match_mol)]
            neighbours["mol"][1] = [k for k, v in bonds_mol.items()
                                    if bond_mol[1] in k
                                    and bond_mol[0] not in k
                                    and set(k).issubset(match_mol)]
            neighbours["ref"][0] = [(match_rev_map[x[0]], match_rev_map[x[1]])
                                    for x in neighbours["mol"][0]]
            neighbours["ref"][1] = [(match_rev_map[x[0]], match_rev_map[x[1]])
                                    for x in neighbours["mol"][1]]

            # no bond breakage if less than 2 adjacent bonds
            single_bonds = True
            for i in ["ref", "mol"]:
                for j in [0, 1]:
                    neighbours[i][j] += [None] * (2 - len(neighbours[i][j]))
                    # don't break in the case we have a symmetric ester
                    single_bonds = single_bonds and None not in neighbours[i][j]

            # planar bonds mapped onto tetrahedral bonds
            if frozenset(bond_ref) not in EZ_ref_mcs.keys():
                if ignore_single_bonds:
                    del bonds_to_break_pairs[(bond_ref, bond_mol)]
                    continue

                # get neighbouring bonds outside of MCS
                breakbonds_ref = {frozenset(x) for x in bonds_ref.keys()
                                  if len(set(x) & set(bond_ref)) == 1}
                breakbonds_ref -= {frozenset(x) for x in
                                   neighbours["ref"][0] + neighbours["ref"][1]
                                   if x is not None}

                # add to unconditional break list
                bonds_to_break_ref += list(breakbonds_ref)

                # regenerate stereochemical information
                extra_bonds = [_revTransformIndices(list(bond_ref),
                                                    indices_to_delete_ref)]
                EZ_ref_mcs_new = getEZStereochemistry(mcs_ref,
                                                      extra_bonds=extra_bonds)
                EZ_ref_mcs_new = {frozenset(
                    _transformIndices(k, indices_to_delete_ref)): v
                    for k, v in EZ_ref_mcs_new.items()}
                EZ_ref_mcs = {**EZ_ref_mcs, **EZ_ref_mcs_new}

            # asymmetric planar bonds mapped onto symmetric planar bonds
            if single_bonds and EZ_ref_mcs[frozenset(bond_ref)] is None:
                symm_bonds = []
                bond_atoms_ref = neighbours["ref"][0] + neighbours["ref"][1]
                bond_atoms_mol = neighbours["mol"][0] + neighbours["mol"][1]

                # determine which bound atoms are symmetric
                for b_ref, b_mol in zip(bond_atoms_ref, bond_atoms_mol):
                    idx_ref, = set(b_ref) - set(bond_ref)
                    idx_mol, = set(b_mol) - set(bond_mol)

                    # transform the indices to MCS indices
                    idx_ref = _revTransformIndices([idx_ref],
                                                   indices_to_delete_ref)[0]
                    idx_mol = _revTransformIndices([idx_mol],
                                                   indices_to_delete_mol)[0]

                    # determine pro-E and pro-Z atoms
                    mcs_ref_c = _carbonify(mcs_ref)
                    mcs_mol_c = _carbonify(mcs_mol)
                    mcs_ref_c.GetAtomWithIdx(idx_ref).SetAtomicNum(14)
                    mcs_mol_c.GetAtomWithIdx(idx_mol).SetAtomicNum(14)

                    idx_bond_ref, = [i for i, x in enumerate(bonds_mcs_ref)
                                     if set(bond_ref) == set(x)]
                    mcs_ref_c.GetBondWithIdx(idx_bond_ref).SetBondType(
                        _Chem.BondType.DOUBLE)
                    idx_bond_mol, = [i for i, x in enumerate(bonds_mcs_mol)
                                     if set(bond_mol) == set(x)]
                    mcs_mol_c.GetBondWithIdx(idx_bond_mol).SetBondType(
                        _Chem.BondType.DOUBLE)

                    t_bond_ref = _revTransformIndices(bond_ref, indices_to_delete_ref)
                    EZ_mcs_ref_c = getEZStereochemistry(mcs_ref_c,
                        carbonify=False, extra_bonds=[t_bond_ref])
                    t_bond_mol = _revTransformIndices(bond_mol, indices_to_delete_mol)
                    EZ_mcs_mol_c = getEZStereochemistry(mcs_mol_c,
                        carbonify=False, extra_bonds=[t_bond_mol])

                    # these are the symmetric atoms
                    EZ_curbond_ref = EZ_mcs_ref_c[frozenset(t_bond_ref)]
                    if EZ_curbond_ref is not None:
                        EZ_curbond_mol = EZ_mcs_mol_c[frozenset(t_bond_mol)]
                        if EZ_curbond_ref != EZ_curbond_mol:
                            symm_bonds += [[b_ref, b_mol]]

                # swap the substituents if we have a mismatch
                if symm_bonds:
                    bonds_to_break_pairs[(bond_ref, bond_mol)] += \
                        [(symm_bonds[0][0], symm_bonds[1][1]),
                         (symm_bonds[1][0], symm_bonds[0][1])]
            # E mapped to Z and vice-versa
            elif EZ_ref_mcs[frozenset(bond_ref)] != chirality_mol:
                # swap the substituents and break bonds on both atoms
                for k2 in [0, 1]:
                    for k3 in [0, 1]:
                        bond1 = neighbours["ref"][k2][k3]
                        bond2 = neighbours["mol"][k2][k3 - 1]
                        bonds_to_break_pairs[(bond_ref, bond_mol)] += [(bond1,
                                                                        bond2)]
            else:
                # all mappings are fine
                del bonds_to_break_pairs[(bond_ref, bond_mol)]

    if not bonds_to_break_pairs:
        return [], [], [], [], []

    # now we finally break the bonds
    refs_prev, mols_prev = [_copy.deepcopy(ref)], [_copy.deepcopy(mol)]
    ids_ref_prev, ids_mol_prev = [[]], [[]]

    for (dbond_ref_orig, dbond_mol_orig), v in bonds_to_break_pairs.items():
        refs, mols, ids_ref, ids_mol = [], [], [], []
        # iterate over previously created fragments
        for ref_prev, mol_prev, id_ref_prev, id_mol_prev in \
                zip(refs_prev, mols_prev, ids_ref_prev, ids_mol_prev):
            # continue if bond is not a part of the fragment
            if len(set(dbond_ref_orig) & set(id_ref_prev)) or \
                    len(set(dbond_mol_orig) & set(id_mol_prev)):
                continue
            # transform the double bond representation back to the original
            # indices
            dbond_ref = _revTransformIndices(dbond_ref_orig, id_ref_prev)
            dbond_mol = _revTransformIndices(dbond_mol_orig, id_mol_prev)

            for breakbond_ref, breakbond_mol in v:
                if breakbond_ref is not None:
                    if len(set(breakbond_ref) & set(id_ref_prev)):
                        continue
                    breakbond_ref = _revTransformIndices(breakbond_ref,
                                                         id_ref_prev)
                if breakbond_mol is not None:
                    if len(set(breakbond_mol) & set(id_mol_prev)):
                        continue
                    breakbond_mol = _revTransformIndices(breakbond_mol,
                                                         id_mol_prev)

                # break the bonds in the reference molecule
                ref_copy = _Chem.EditableMol(_copy.deepcopy(ref_prev))
                for bond in [breakbond_ref] + bonds_to_break_ref:
                    if bond:
                        ref_copy.RemoveBond(*bond)
                ref_copy = ref_copy.GetMol()
                ref_frags = _rdmolops.GetMolFrags(ref_copy,
                                                  sanitizeFrags=False)
                deleted_frags_ref = list(_it.chain(
                    *[list(x) for x in ref_frags
                      if not set(dbond_ref).issubset(x)]))
                ids_ref += [_transformIndices(deleted_frags_ref, id_ref_prev)
                            + id_ref_prev]
                refs += [_generateFragment(ref_copy, deleted_frags_ref)]

                # break the bonds in the target molecule
                mol_copy = _Chem.EditableMol(_copy.deepcopy(mol_prev))
                if breakbond_mol:
                    mol_copy.RemoveBond(*breakbond_mol)
                mol_copy = mol_copy.GetMol()
                mol_frags = _rdmolops.GetMolFrags(mol_copy,
                                                  sanitizeFrags=False)
                deleted_frags_mol = list(_it.chain(
                    *[list(x) for x in mol_frags
                      if not set(dbond_mol).issubset(x)]))
                ids_mol += [_transformIndices(deleted_frags_mol, id_mol_prev)
                            + id_mol_prev]
                mols += [_generateFragment(mol_copy, deleted_frags_mol)]

        bonds_to_break_ref = []

        # update old values
        refs_prev, mols_prev = refs, mols
        ids_ref_prev, ids_mol_prev = ids_ref, ids_mol

    # get new pruned MCS for every broken molecule
    mcss = [[(x, y) for x, y in zip(match_ref, match_mol)
             if x not in r and y not in m]
             for r, m in zip(ids_ref_prev, ids_mol_prev)]

    return refs_prev, mols_prev, mcss, ids_ref_prev, ids_mol_prev

def _breakMismatchingBonds(ref, mol, match_ref, match_mol):
    """
    Breaks all bonds bordering on an MCS if they are a subset of another MCS.
    Helpful at breaking rings in order to utilise the strict MCS algorithm.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    match_ref : [int]
        The indices of the matched reference.
    match_mol : [int]
        The indices of the matched molecule.

    Returns
    -------
    ref : [rdkit.Chem.rdchem.Mol]
        The reference molecules with broken bonds.
    mol : [rdkit.Chem.rdchem.Mol]
        The molecules with broken bonds to be aligned.
    """
    # get all bonds
    bonds_ref = [frozenset([x.GetBeginAtomIdx(), x.GetEndAtomIdx()])
                 for x in ref.GetBonds()]
    bonds_mol = [frozenset([x.GetBeginAtomIdx(), x.GetEndAtomIdx()])
                 for x in mol.GetBonds()]

    # only keep bordering bonds
    bonds_ref = {x for x in bonds_ref if len(x | set(match_ref))
                 == len(match_ref) + 1}
    bonds_mol = {x for x in bonds_mol if len(x | set(match_mol))
                 == len(match_mol) + 1}

    # break the bonds if they don't fragment the molecule
    refs_broken = []
    for bond in bonds_ref:
        ref_copy = _Chem.EditableMol(_copy.deepcopy(ref))
        ref_copy.RemoveBond(*bond)
        ref_copy = ref_copy.GetMol()
        if len(_rdmolops.GetMolFrags(ref_copy, sanitizeFrags=False)) == 1:
            ref_copy.UpdatePropertyCache()
            refs_broken += [ref_copy]

    mols_broken = []
    for bond in bonds_mol:
        mol_copy = _Chem.EditableMol(_copy.deepcopy(mol))
        mol_copy.RemoveBond(*bond)
        mol_copy = mol_copy.GetMol()
        if len(_rdmolops.GetMolFrags(mol_copy, sanitizeFrags=False)) == 1:
            mol_copy.UpdatePropertyCache()
            mols_broken += [mol_copy]

    return refs_broken, mols_broken


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
    for bond in mol_c.GetBonds():
        bond.SetBondType(_Chem.rdchem.BondType.SINGLE)
    _Chem.AssignStereochemistry(mol_c, cleanIt=True, force=True)
    return mol_c


def _generateFragment(mol, indices_to_delete, getAsFrags=False):
    """
    Generates a molecule from another given indices of atoms to be deleted.
    It can either return the new molecule or the indices of the resulting
    fragments.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        The input molecule.
    indices_to_delete : list
        A list of the atom indices to be deleted.
    getAsFrags : bool
        False returns an rdkit.Chem.Mol object, True returns a list of tuples.

    Returns
    -------
    mol_frag :
        Either an rdkit.Chem.Mol object or a list of tuples corresponding to
        the new molecule.
    """
    mol_frag_edit = _Chem.EditableMol(_copy.deepcopy(mol))
    for idx in reversed(sorted(set(indices_to_delete))):
        mol_frag_edit.RemoveAtom(idx)
    mol_frag = mol_frag_edit.GetMol()
    if getAsFrags:
        frags = _rdmolops.GetMolFrags(mol_frag, sanitizeFrags=False)
        return [_transformIndices(f, indices_to_delete) for f in frags]
    else:
        mol_frag.UpdatePropertyCache()
        return mol_frag


def _haveCommonElements(set1, set2):
    """
    Determines whether two sets of tuples have common tuple elements. They may
    not necessarily have common elements.

    Parameters
    ----------
    set1 : set
        The first input set.
    set2 : set
        The second input set.

    Returns
    -------
    have_common_elements : bool
        Whether the sets have common tuple elements.
    """
    if not set1 or not set2:
        return set()
    tuple11, tuple12 = list(zip(*set1))
    tuple21, tuple22 = list(zip(*set2))
    return set(tuple11).intersection(tuple21) and \
           set(tuple12).intersection(tuple22)


def _matchAndReturnMatches(*args, **kwargs):
    """
    A light wrapper around RDKit's FindMCS function. Here we use the deprecated
    version of the function because it supports strict matching of rings to
    rings of the same size.

    Parameters
    ----------
    args:
        Positional arguments to be given to FindMCS. The first argument is
        assumed to be input molecules.
    kwargs:
        Keyword arguments to be given to FindMCS.

    Returns
    -------
    transformed_indices : [int]
        Transformed indices.
    """
    mcs_string = _MCS.FindMCS(*args, **kwargs).smarts
    if mcs_string is None:
        return None, [], []
    mcs = _Chem.MolFromSmarts(mcs_string)
    # fixes some issues with RDKit concerning SMARTS
    mcs = _Chem.MolFromSmarts(_Chem.MolToSmiles(mcs, canonical=False,
                                                allBondsExplicit=True,
                                                allHsExplicit=True))
    return tuple([mcs] + [set(x.GetSubstructMatches(mcs)) for x in args[0]])


def _nonMCSDihedrals(mol, mcs=None, confId=-1):
    """
    Retrieves all non-ring dihedrals which are not a part of the MCS. In case
    they are bordering dihedrals, the ordering is done so that directly
    feeding the resulting dihedrals into RDKit's SetDihedralRad() does not move
    the common core.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    mcs : [(tuple)] or None
        The maximum common substructure of the two molecules.
    confId : int
        The conformer number of mol.

    Returns
    -------
    dihedrals : dict((int, int, int, int), float)
        All of the non-ring dihedrals outside of the MCS and their values in
        radians.
    """
    if mcs is None:
        frozen_atoms = []
    else:
        frozen_atoms = list(zip(*mcs))[1]

    mol = _copy.deepcopy(mol)
    mol_conf = mol.GetConformer(confId)

    # only rotate if the bond belongs to a dihedral
    pattern = _Chem.MolFromSmarts("[*]~[*;X2,X3,X4]~&!@[*;X2,X3,X4]~[*]")
    dihedrals_ini =  mol.GetSubstructMatches(pattern)
    dihedrals_unique = {frozenset(x[1:3]): x for x in reversed(dihedrals_ini)}

    dihedrals = {}

    for i0, i1, i2, i3 in dihedrals_unique.values():
        mol_temp = _Chem.EditableMol(_copy.deepcopy(mol))
        mol_temp.RemoveBond(i1, i2)
        frags = _Chem.GetMolFrags(mol_temp.GetMol(), sanitizeFrags=False)

        for frag in frags:
            n_frozen = len(set(frozen_atoms) & set(frag))
            if n_frozen <= 1:
                if i1 in frag:
                    dihedral = (i3, i2, i1, i0)
                else:
                    dihedral = (i0, i1, i2, i3)

                dihedrals[dihedral] = _Transforms.GetDihedralRad(mol_conf,
                                                                 *dihedral)
                break

    return dihedrals


def _onlyKeepLongest(input_set):
    """
    Only keeps the longest elements of the input set.

    Parameters
    ----------
    input_set : set
        The input set to be pruned.

    Returns
    -------
    pruned_set : set
        The resulting pruned set.
    max_len : float
        The length of the largest element(s).
    """
    if not input_set:
        return set(), 0
    max_len = max([len(x) for x in input_set])
    return {x for x in input_set if len(x) == max_len}, max_len


def _optimalMergedSets(*sets, seed=None):
    """
    Generates the optimal merged set from a sequence of sets of tuples. If
    two sets have common tuple elements but not common tuples they are deemed
    incompatible. The algorithm finds the longest combination(s) of these sets.

    Parameters
    ----------
    seed : set
        A set that will definitely be present in the returned result.
    sets : set
        The sets to be added

    Returns
    -------
    sets_final : {frozenset([tuple])}
        A set containing all the optimal set(s).
    max_len_final:
        The length of the optimal set(s).
    """
    def increaseByOne(inpsets):
        outpsets = set()
        for inpset in inpsets:
            all_sets = [compatible_dict[i] for i in inpset]
            available_sets = set(all_sets[0]).intersection(*all_sets[1:])
            available_sets -= inpset
            for x in available_sets:
                outpsets |= {frozenset([*inpset, x])}
        return outpsets

    def getLen(number_set):
        return len(set().union(*[sets[i] for i in number_set]))

    def getSet(number_set):
        return set().union(*[sets[i] for i in number_set])

    if len(sets) <= 1:
        return _onlyKeepLongest(set([frozenset(x) for x in sets]))

    if seed is None:
        seed = set()
    sets = [seed, *sets]

    # get compatible pairs of sets in terms of their position
    compatible_sets = []
    for i in range(len(sets)):
        set1 = sets[i]
        for j in range(i):
            set2 = sets[j]
            if _areCompatibleSets(set1, set2):
                compatible_sets += [{i, j}]

    # turn compatible_sets into a dictionary as well
    compatible_dict = {}
    for i in range(len(sets)):
        compatible_dict[i] = {i}.union(*[x for x in compatible_sets if i in x])

    # build up larger subsets that obey all the rules
    # we do this by translating the sets into numbers corresponding to set
    # positions
    sets_incl = {frozenset([0])}
    sets_final = set()
    max_len_final = 0
    for size in range(len(sets)):
        max_len = max([getLen(x) for x in sets_incl]) if sets_incl else 0

        if max_len > max_len_final:
            max_len_final = max_len
            sets_final = {x for x in sets_incl if getLen(x) == max_len}
        elif max_len == max_len_final:
            sets_final |= {x for x in sets_incl if getLen(x) == max_len}
        else:
            break

        if size != len(sets) - 1:
            sets_incl = increaseByOne(sets_incl)

    # translate back into the actual sets
    sets_final = {frozenset(getSet(x)) for x in sets_final}

    return sets_final, max_len_final


def _revTransformIndices(prev_indices, deleted_indices):
    """
    Generates the new indices of the molecule given the old and the deleted
    indices.

    Parameters
    ----------
    prev_indices : [int]
        Previous indices that are to be transformed.
    deleted_indices : [int]
        Deleted indices from the previous molecule.

    Returns
    -------
    transformed_indices : [int]
        Transformed indices.
    """
    if not prev_indices:
        return []
    if not deleted_indices:
        return prev_indices
    if set(prev_indices) & set(deleted_indices):
        raise ValueError("Cannot transform deleted indices.")
    return [x - sum(x >= y for y in sorted(deleted_indices))
            for x in sorted(prev_indices)]


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
