from rdkit import Chem as _Chem
from rdkit.Chem import AllChem as _AllChem
from rdkit.Chem import rdMolAlign as _MolAlign
from rdkit.Chem import rdForceFieldHelpers as _FF
from rdkit.Chem import rdFMCS as _FMCS

def openFileAsRdkit(filename, **kwargs):
    extension = filename.split(".")[-1]

    if extension == "mol2":
        mol = _Chem.MolFromMol2File(filename, **kwargs)
    else:
        exec("mol = _Chem.MolFrom%sFile(filename, **kwargs)" % extension.upper(), locals(), globals())

    if (mol is None):
        raise ValueError("File format not recognised: %s" % filename)
    return mol

def openInChIAsRdkit(inchi_str, removeHs=True, **kwargs):
    mol = _Chem.MolFromInchi(inchi_str, removeHs=removeHs, **kwargs)
    if(mol is None):
        raise ValueError("InChI string not recognised: %s" % inchi_str)
    return mol

def openSmilesAsRdkit(smiles_str, **kwargs):
    mol = _Chem.MolFromSmiles(smiles_str, **kwargs)
    if (mol is None):
        raise ValueError("SMILES string not recognised: %s" % smiles_str)
    return mol

def openAsRdkit(val, **kwargs):
    if len(val.split(".")) == 1:
        try:
            mol = openSmilesAsRdkit(val, **kwargs)
            _AllChem.EmbedMolecule(mol, _AllChem.ETKDG())
        except:
            try:
                mol = openInChIAsRdkit(val, **kwargs)
                _AllChem.EmbedMolecule(mol, _AllChem.ETKDG())
            except:
                raise ValueError("String not recognised as a valid SMILES or InChI input")
    else:
        mol = openFileAsRdkit(val, **kwargs)
    return mol

def matchTwoMolecules(mol, ref, n_min=-1):
    mcs = _Chem.MolFromSmarts(_FMCS.FindMCS([ref, mol], atomCompare=_FMCS.AtomCompare.CompareAny,
                                            ringMatchesRingOnly=True, completeRingsOnly=True).smartsString)
    match1 = ref.GetSubstructMatch(mcs)
    match2 = mol.GetSubstructMatch(mcs)

    ref_conf = ref.GetConformer(-1)
    mol_conf = mol.GetConformer(-1)

    ff = _FF.MMFFGetMoleculeForceField(mol, _FF.MMFFGetMoleculeProperties(mol), confId=0)
    for i_ref, i_mol in zip(match1, match2):
        mol_conf.SetAtomPosition(i_mol, ref_conf.GetAtomPosition(i_ref))
        ff.AddFixedPoint(i_mol)

    ff.Initialize()
    more = ff.Minimize()
    while more and n_min:
        more = ff.Minimize()
        #n_min = -1 means infinite minimisation
        if n_min != -1:
            n_min -= 1
    return mol