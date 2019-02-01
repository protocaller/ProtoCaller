import os as _os

from rdkit import Chem as _Chem
from rdkit.Chem import AllChem as _AllChem
from rdkit.Chem import rdForceFieldHelpers as _FF
from rdkit.Chem import rdFMCS as _FMCS
from rdkit.Geometry import rdGeometry as _Geom
import parmed as _pmd

import ProtoCaller.Utils.stdio as _stdio
from . import babelwrapper as _babel

def openFileAsRdkit(filename, **kwargs):
    extension = filename.split(".")[-1].lower()

    func_dict = {
        "mol": _Chem.MolFromMolFile,
        "mol2": _Chem.MolFromMol2File,
        "pdb":_Chem.MolFromPDBFile,
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

def alignTwoMolecules(ref, mol, n_min=-1, match="any"):
    """
    :param mol: molecule to be aligned
    :param ref: reference molecule
    :param n_min: number of minimisation iterations before minimisation stops. Default is -1 (infinite)
    :param match: "any" matches any pairs of atoms - useful for single-topology;
                  "elements" matches only the same atoms - useful for dual-topology
    :return: the aligned molecule
    """
    if match == "any":
        comp_method = _FMCS.AtomCompare.CompareAny
    elif match == "elements":
        comp_method = _FMCS.AtomCompare.CompareElements
    else:
        raise ValueError("'match' arguments needs to be either 'any' or 'elements'")

    #now we match aliphatic chains to rings, but not rings to other rings with different size
    #this functionality is not at the moment present in RDKit, so we implement it manually
    #since ring breaking is not good practice, we throw an error if we could match two different rings to each other
    #TODO maybe fix this behaviour later on?
    mcs_string_strict = _FMCS.FindMCS([ref, mol], atomCompare=comp_method, bondCompare=_FMCS.BondCompare.CompareAny,
                                       ringMatchesRingOnly=False, completeRingsOnly=True).smartsString
    mcs_string_loose = _FMCS.FindMCS([ref, mol], atomCompare=comp_method, bondCompare=_FMCS.BondCompare.CompareAny,
                                      ringMatchesRingOnly=False, completeRingsOnly=False).smartsString

    mcs_strict = _Chem.MolFromSmarts(mcs_string_strict)
    mcs_loose = _Chem.MolFromSmarts(mcs_string_loose)

    match1_strict, match2_strict = ref.GetSubstructMatch(mcs_strict), mol.GetSubstructMatch(mcs_strict)
    match1_loose, match2_loose = ref.GetSubstructMatch(mcs_loose), mol.GetSubstructMatch(mcs_loose)
    
    #find extra pairs of atoms found with loose matching
    extra_pairs = list(set(zip(match1_loose, match2_loose)) - set(zip(match1_strict, match2_strict)))

    #If any extra pairs contain matched atoms contained in different rings,
    #then this is not a complete ring and we throw an error.
    #Otherwise, we pick the loose match
    for i_ref, i_mol in extra_pairs:
        rings_ref, rings_mol = getRings(i_ref, ref), getRings(i_mol, mol)
        if [] not in [rings_ref, rings_mol] and rings_ref != rings_mol:
            raise ValueError("MCS cannot match molecules with rings of different sizes")

    ref_conf = ref.GetConformer(-1)
    mol_conf = mol.GetConformer(-1)

    ff = _FF.MMFFGetMoleculeForceField(mol, _FF.MMFFGetMoleculeProperties(mol), confId=0)
    for i_ref, i_mol in zip(match1_loose, match2_loose):
        mol_conf.SetAtomPosition(i_mol, ref_conf.GetAtomPosition(i_ref))
        ff.AddFixedPoint(i_mol)

    if mol_conf.GetNumAtoms() != len(match1_loose):
        ff.Initialize()
        more = ff.Minimize()
        while more and n_min:
            more = ff.Minimize()
            #n_min = -1 means infinite minimisation
            if n_min != -1:
                n_min -= 1

    return mol, list(zip(match1_loose, match2_loose))

#returns an array with the sorted ring sizes of all rings the target atom is contained in
def getRings(atom_idx, molecule):
    rings = molecule.GetRingInfo().AtomRings()
    return sorted([len(ring) for ring in rings if atom_idx in ring])