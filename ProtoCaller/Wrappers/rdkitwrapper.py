import copy as _copy
import os as _os
import sys as _sys

import numpy as _np
from rdkit import Chem as _Chem
from rdkit.Chem import AllChem as _AllChem
from rdkit.Chem import rdForceFieldHelpers as _FF
from rdkit.Chem import rdmolops as _rdmolops
from rdkit.Chem import rdMolTransforms as _Transforms
from rdkit.Chem import TorsionFingerprints as _Torsions
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
    kwargs = {
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
        matches_rec_prev = {x: {x} for x in matches}
        for match in matches:
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
                        mcs_broken, matches_ref_broken, matches_mol_broken = \
                            _matchAndReturnMatches([ref_broken, mol_broken],
                                                   **kwargs)
                        matches_broken = {frozenset(zip(x, y))
                                          for x in matches_ref_broken
                                          for y in matches_mol_broken}
                        # only keep the matches that are connected to the
                        # original one
                        matches_rec |= {x for x in matches_broken
                                        if _haveCommonElements(x, submatch)}

                # combine the new MCS's in an optimal way
                matches_rec, max_len_rec = \
                    _optimalMergedSets(*matches_rec_prev[match], *matches_rec)

                # break the recursion if the new MCS is not larger
                if max_len_rec > len(max(matches_rec_prev[match], key=len)):
                    matches_rec_prev[match] = matches_rec
                else:
                    break

        matches = set().union(*[matches_rec_prev[x] for x in matches])

    # here we deal with mismatching atoms of different chirality
    matches_final = set()
    # take care of R/S changing with different atom types
    _Chem.AssignStereochemistry(ref, cleanIt=True, force=True)
    _Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    ref_c = _carbonify(ref)
    mol_c = _carbonify(mol)
    chiral_ref_c = dict(_Chem.FindMolChiralCenters(ref_c))
    chiral_mol_c = dict(_Chem.FindMolChiralCenters(mol_c))

    # check if both atoms are chiral and flag the atoms if this is the case
    for match in matches:
        if match == frozenset():
            continue
        # only deal with one of the molecules
        chiral_ref_indices = [x[0] for x in match
                              if x[0] in chiral_ref_c.keys()
                              and x[1] in chiral_mol_c.keys()
                              and chiral_ref_c[x[0]] != chiral_mol_c[x[1]]]

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

        frags_ref_total = _onlyKeepLongest(frags_ref_total)[0]
        for frag_ref in frags_ref_total:
            matches_final |= {frozenset([x for x in match
                                         if x[0] in frag_ref])}

    matches_final, max_len_final = _onlyKeepLongest(matches_final)
    if not matches_final:
        return {frozenset()}

    return matches_final


def alignTwoMolecules(ref, mol, n_min=-1, two_way_matching=True, mcs=None,
                      mcs_parameters=None, minimiser_parameters=None):
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
        A dictionary of the parameters to be apssed on to
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

    # if we have multiple equivalent MCS's we pick the one with best MSD
    mol_final, score_final, mcs_final = None, None, None
    for mcs in mcss:
        ref_conf = ref.GetConformer(-1)

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
                break
            except:
                # sometimes MMFF fails and in this case we use UFF
                ff_to_use = "uff"
                continue

        frozen_atoms = list(zip(*mcs))[1]
        mol, score = minimiseAlignmentScore(ref, mol,
                                            frozen_atoms=frozen_atoms,
                                            **minimiser_parameters)

        if score_final is None or score < score_final:
            mol_final = mol
            score_final = score
            mcs_final = mcs

    return mol_final, mcs_final


def minimiseAlignmentScore(ref, mol, frozen_atoms=None, confId1=-1, confId2=-1,
                           minimisation_algorithm=_minimize, **kwargs):
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
    frozen_atoms : list
        The list of indices of the atoms in mol not to be moved.
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
    mcs : {frozenset([tuple])}
        A set of frozensets of tuples corresponding to the atom index matches
        between the reference and the other molecule.
    alignment_score : float
        The alignment score after alignment obtained from getAlignmentScore()
    """
    if minimisation_algorithm is None:
        return _copy.deepcopy(mol), getAlignmentScore(ref, mol,
                                                      confId1=confId1,
                                                      confId2=confId2)
    if frozen_atoms is None:
        frozen_atoms = []

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
    smarts = "[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]"
    rbonds = _Chem.MolFromSmarts(smarts)
    rbond_indices = mol.GetSubstructMatches(rbonds)
    rbond_indices = [{i1, i2} for i1, i2 in rbond_indices
                     if n_frozen_neighbours[i1] <= 1
                     or n_frozen_neighbours[i2] <= 1]

    # only rotate if the bond belongs to a dihedral
    dihedrals_ini = [x[0][0] for x in _Torsions.CalculateTorsionLists(mol)[0]]
    dihedrals_ini = [x for x in dihedrals_ini if set(x[1:3]) in rbond_indices]

    dihedrals = []
    initial_values = []

    for i0, i1, i2, i3 in dihedrals_ini:
        mol_temp = _Chem.EditableMol(_copy.deepcopy(mol))
        mol_temp.RemoveBond(i1, i2)
        frags = _Chem.GetMolFrags(mol_temp.GetMol(),
                                  asMols=False, sanitizeFrags=False)

        for frag in frags:
            n_frozen = len(set(frozen_atoms) & set(frag))
            if n_frozen > 1:
                if i3 in frag:
                    dihedral = (i3, i2, i1, i0)
                else:
                    dihedral = (i0, i1, i2, i3)
                break

        dihedrals += [dihedral]
        initial_values += [_Transforms.GetDihedralRad(mol_conf, *dihedral)]

    def rotateDihedral(*args):
        nonlocal dihedrals, mol, confId1

        mol_temp = _copy.deepcopy(mol)
        mol_conf_temp = mol_temp.GetConformer(confId1)
        for dihedral, arg in zip(dihedrals, list(args[0])):
            _Transforms.SetDihedralRad(mol_conf_temp, *dihedral, arg)

        return getAlignmentScore(ref, mol_temp, confId1=confId1,
                                 confId2=confId2)

    # rotate dihedrals until MSD minimisation
    if len(initial_values):
        result = minimisation_algorithm(rotateDihedral,
                                        _np.asarray(initial_values), **kwargs)
        optimal_angles, final_alignment_score = list(result.x), result.fun
        for dihedral, optimal_angle in zip(dihedrals, optimal_angles):
            _Transforms.SetDihedralRad(mol_conf, *dihedral, optimal_angle)
        return mol, final_alignment_score
    else:
        return mol, getAlignmentScore(ref, mol, confId1=confId1,
                                      confId2=confId2)


def getAlignmentScore(ref, mol, confId1=-1, confId2=-1):
    """
    Returns the alignment score between two molecules. The way this is done is
    by calculating all possible distances between atom i of ref and atom j of
    mol and computing the sum of their squares. All atoms from mol that are
    within 1 Anstrom of atom j contribute inverse square distances to the score
    in order to prevent unfavourable clashes.

    Parameters
    ----------
    ref : rdkit.Chem.rdchem.Mol
        The reference molecule.
    mol : rdkit.Chem.rdchem.Mol
        The molecule to be aligned.
    confId1 : int
        The conformer number of ref.
    confId2 : int
        The conformer number of mol.

    Returns
    -------
    alignment_score : float
        The alignment score of the two conformers.
    """
    ref_conf = ref.GetConformer(confId1)
    mol_conf = mol.GetConformer(confId2)

    alignment_score = 0
    for mol_idx in range(mol.GetNumAtoms()):
        for ref_idx in range(ref.GetNumAtoms()):
            x1 = ref_conf.GetAtomPosition(ref_idx)
            x2 = mol_conf.GetAtomPosition(mol_idx)
            alignment_score += (x1 - x2).LengthSq()

        for mol_idx2 in range(mol_idx):
            x1 = mol_conf.GetAtomPosition(mol_idx)
            x2 = mol_conf.GetAtomPosition(mol_idx2)
            len_sq = (x1 - x2).LengthSq()
            if len_sq > 1:
                continue
            elif len_sq < 2 / _sys.maxsize:
                return _sys.maxsize
            else:
                alignment_score += 2 / len_sq

    return alignment_score


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
        if len(_rdmolops.GetMolFrags(ref_copy, asMols=False,
                                     sanitizeFrags=False)) == 1:
            ref_copy.UpdatePropertyCache()
            refs_broken += [ref_copy]

    mols_broken = []
    for bond in bonds_mol:
        mol_copy = _Chem.EditableMol(_copy.deepcopy(mol))
        mol_copy.RemoveBond(*bond)
        mol_copy = mol_copy.GetMol()
        if len(_rdmolops.GetMolFrags(mol_copy, asMols=False,
                                     sanitizeFrags=False)) == 1:
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
        frags = _rdmolops.GetMolFrags(mol_frag, asMols=False,
                                      sanitizeFrags=False)
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
        return {}, 0
    max_len = max([len(x) for x in input_set])
    return {x for x in input_set if len(x) == max_len}, max_len


def _optimalMergedSets(*sets):
    """
    Generates the optimal merged set from a sequence of sets of tuples. If
    two sets have common tuple elements but not common tuples they are deemed
    incompatible. The algorithm finds the longest combination(s) of these sets.

    Parameters
    ----------
    sets : set
        The sets to be added

    Returns
    -------
    sets_final : {frozenset([tuple])}
        A set containing all the optimal set(s).
    max_len_final:
        The length of the optimal set(s).
    """
    def increaseByOne(inpsets, max_len):
        outpsets = set()
        max_len_set = {i for i in range(max_len)}
        for inpset in inpsets:
            for x in max_len_set - inpset:
                outpsets |= {frozenset([*inpset, x])}
        return outpsets

    def isCompatibleSet(inpset):
        return not any([x.issubset(inpset) for x in incompatible_sets])

    def getLen(number_set):
        return len(set().union(*[sets[i] for i in number_set]))

    def getSet(number_set):
        return set().union(*[sets[i] for i in number_set])

    if len(sets) <= 1:
        return set([frozenset(x) for x in sets]), len(sets)

    # get maximum possible lengths of the sets
    tuples1 = set().union(*[list(zip(*x))[0] for x in sets])
    tuples2 = set().union(*[list(zip(*x))[1] for x in sets])
    lengths = [len(x) for x in sets]
    max_lens = [sum(sorted(lengths[:i+1], reverse=True))
                for i in range(len(sets))]
    max_lens = [min(x, len(tuples1), len(tuples2)) for x in max_lens]

    # get incompatible pairs of sets in terms of their position
    incompatible_sets = []
    for i in range(len(sets)):
        set1 = sets[i]
        for j in range(i):
            set2 = sets[j]
            set12 = set1.union(set2)
            tuple1, tuple2 = zip(*set12)
            if not (len(set(tuple1)) == len(set(tuple2)) == len(set12)):
                incompatible_sets += [{i, j}]

    # build up larger subsets that obey all the rules
    # we do this by translating the sets into numbers corresponding to set
    # positions
    sets_incl = {frozenset([i]) for i in range(len(sets))}
    sets_final = set()
    max_len_final = 0
    for size in range(len(sets)):
        max_len = max([getLen(x) for x in sets_incl]) if sets_incl else 0
        if max_len:
            if max_len > max_len_final:
                max_len_final = max_len
                sets_final = {x for x in sets_incl if getLen(x) == max_len}
            elif max_len == max_len_final:
                sets_final |= {x for x in sets_incl if getLen(x) == max_len}
            else:
                break

        if size == len(sets) - 1 or max_len == max_lens[size] \
                or max_len < max_len_final:
            break

        sets_incl = increaseByOne(sets_incl, len(sets))
        sets_incl = {x for x in sets_incl if isCompatibleSet(x)}

    # translate back into the actual sets
    sets_final = {frozenset(getSet(x)) for x in sets_final}

    return sets_final, max_len_final


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
