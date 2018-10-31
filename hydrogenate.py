from rdkit import Chem as _Chem
from rdkit.Chem import AllChem as _AllChem
from rdkit.Chem import rdmolops as _MolOps
from rdkit.Chem import rdMolAlign as _MolAlign
from rdkit.Chem import rdmolfiles as _MolFiles
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

def openUnknownAsRdkit(val, **kwargs):
    if len(val.split(".")) == 1:
        try:
            mol = openSmilesAsRdkit(val, **kwargs)
            _AllChem.EmbedMolecule(mol, _AllChem.ETKDG())
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

def matchTwoMolecules(mol, ref, **kwargs):
    mol1 = openUnknownAsRdkit(mol, **kwargs)
    mol2 = openUnknownAsRdkit(ref, **kwargs)
    #print(_MolFiles.MolToMolBlock(mol1, kekulize=False))
    #_AllChem.EmbedMolecule(mol2, _AllChem.ETKDG())

    #find the maximum common substructure of both molecules
    _MolOps.RemoveHs(mol1)
    _MolOps.RemoveHs(mol2)
    mcs = _Chem.MolFromSmarts(_FMCS.FindMCS([mol1, mol2]).smartsString)
    match1 = mol1.GetSubstructMatch(mcs)
    match2 = mol2.GetSubstructMatch(mcs)
    core = _AllChem.DeleteSubstructs(_AllChem.ReplaceSidechains(_Chem.RemoveHs(mol2), mcs), _Chem.MolFromSmiles('*'))
    core.UpdatePropertyCache()
    #core = _Chem.PathToSubmol(mol2, mcs)
    #core.UpdatePropertyCache()
    mcs.AddConformer(core.GetConformer())
    mmmm = mcs.GetSubstructMatch(mcs)
    #matchhhh = mol1.GetSubstructMatch(mcs)
    #return mol1
    GetFF = lambda x, confId=-1: _AllChem.MMFFGetMoleculeForceField(x, _AllChem.MMFFGetMoleculeProperties(x),
                                                                   confId=confId)
    #align the first molecule to the second one
    #_AllChem.ConstrainedEmbed(mol1, mcs, getForceField=GetFF)
    _MolAlign.AlignMol(mol1, mcs, atomMap=list(zip(mmmm, match1)))

    #return hydrogenated molecule with correct geometry
    print(_MolFiles.MolToMolBlock(mol1, kekulize=False))
    print(_MolFiles.MolToMolBlock(mcs, kekulize=False))
    print(mol1.GetProp('EmbedRMS'))
    return _MolOps.AddHs(mol1)