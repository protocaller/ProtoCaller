import os as _os

from rdkit import Chem as _Chem
from rdkit.Chem import AllChem as _AllChem
from rdkit.Chem import rdForceFieldHelpers as _FF
from rdkit.Chem import rdFMCS as _FMCS
from rdkit.Chem import rdmolops as _rdmolops
from rdkit.Geometry import rdGeometry as _Geom
import parmed as _pmd

import ProtoCaller.Utils.stdio as _stdio
from . import babelwrapper as _babel


def openFileAsRdkit(filename, **kwargs):
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
    mol = _Chem.MolFromInchi(inchi_str, removeHs=removeHs, **kwargs)
    if mol is None:
        raise ValueError("InChI string not recognised: %s" % inchi_str)
    return mol


@_stdio.stdout_stderr()
def openSmilesAsRdkit(smiles_str, **kwargs):
    mol = _Chem.MolFromSmiles(smiles_str, **kwargs)
    if mol is None:
        raise ValueError("SMILES string not recognised: %s" % smiles_str)
    return mol


def openAsRdkit(val, **kwargs):
    if isinstance(val, str) and len(val.split(".")) == 1:
        flist = [openSmilesAsRdkit, openInChIAsRdkit]
        for i, f in enumerate(flist):
            try:
                mol = f(val, **kwargs)
                _AllChem.EmbedMolecule(mol, _AllChem.ETKDG())
                break
            except:
                if i == len(flist) - 1:
                    raise ValueError("String not recognised as a valid SMILES or InChI input")
    else:
        try:
            mol = openFileAsRdkit(val, **kwargs)
        except:
            raise ValueError("File is not in a valid format")

    return mol


def saveFromRdkit(mol, filename, **kwargs):
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
    mol_conf = mol.GetConformer(-1)
    for i in range(mol.GetNumAtoms()):
        new_atom_position = mol_conf.GetAtomPosition(i) + _Geom.Point3D(*vector)
        mol_conf.SetAtomPosition(i, new_atom_position)
    return mol


def getMCSMap(ref, mol, match="any", always_maximum=True):
    """
    :param mol: molecule to be aligned
    :param ref: reference molecule
    :param always_maximum: always give the substructure that has the maximum common number of atoms - good for single
                           topology mapping. False is good for overall structure alignment. Note that in some more
                           complex molecules there might be an MCS that satisfies our criteria and is bigger. Here we
                           rely on there being one obvious MCS - something that should be the case in free energy
                           calculations anyway
    :param match: "any" matches any pairs of atoms - useful for single-topology;
                  "elements" matches only the same atoms - useful for dual-topology
    :return: the aligned molecule
    """
    match = match.strip().lower()
    if match == "any":
        comp_method = _FMCS.AtomCompare.CompareAny
    elif match == "elements":
        comp_method = _FMCS.AtomCompare.CompareElements
    else:
        raise ValueError("'match' arguments needs to be either 'any' or 'elements'")

    # find loose MCS
    mcs_string = _FMCS.FindMCS([ref, mol], atomCompare=comp_method, bondCompare=_FMCS.BondCompare.CompareAny,
                               ringMatchesRingOnly=False, completeRingsOnly=False).smartsString
    mcs = _Chem.MolFromSmarts(mcs_string)
    matches = list(zip(ref.GetSubstructMatch(mcs), mol.GetSubstructMatch(mcs)))

    # check whether there are any improperly assigned rings (e.g. 5-membered to 6-membered) and flag the atoms
    indices_to_delete = []
    for i_mcs, (i_ref, i_mol) in enumerate(matches):
        rings_ref, rings_mol = getRings(i_ref, ref), getRings(i_mol, mol)
        if not set(rings_ref).issubset(rings_mol) and not set(rings_mol).issubset(rings_ref):
            indices_to_delete += [i_mcs]

    # delete these atoms from the MCS and get the biggest fragment
    if indices_to_delete:
        mcs = _Chem.EditableMol(mcs)
        for index in reversed(indices_to_delete):
            mcs.RemoveAtom(index)

        indices_new = max(_rdmolops.GetMolFrags(mcs.GetMol(), asMols=False, sanitizeFrags=False), key=len)
        matches = [x for i, x in enumerate(matches) if i in indices_new]

    # now get the strict mapping and discard the loose mapping if it doesn't perform better
    if always_maximum:
        mcs_string_strict = _FMCS.FindMCS([ref, mol], atomCompare=comp_method, bondCompare=_FMCS.BondCompare.CompareAny,
                                          ringMatchesRingOnly=False, completeRingsOnly=True).smartsString
        mcs_strict = _Chem.MolFromSmarts(mcs_string_strict)
        matches_strict = list(zip(ref.GetSubstructMatch(mcs_strict), mol.GetSubstructMatch(mcs_strict)))
        matches = matches if len(matches) > len(matches_strict) else matches_strict

    return matches


def alignTwoMolecules(ref, mol, n_min=-1, mcs=None, **kwargs):
    if mcs is None:
        mcs = getMCSMap(ref, mol, **kwargs)

    ref_conf = ref.GetConformer(-1)
    mol_conf = mol.GetConformer(-1)

    ff = _FF.MMFFGetMoleculeForceField(mol, _FF.MMFFGetMoleculeProperties(mol), confId=0)
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


# returns an array with the sorted ring sizes of all rings the target atom is contained in
def getRings(atom_idx, molecule):
    rings = molecule.GetRingInfo().AtomRings()
    return sorted([len(ring) for ring in rings if atom_idx in ring])
