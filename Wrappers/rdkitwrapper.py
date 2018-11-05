import tempfile as _tempfile

from rdkit import Chem as _Chem
from rdkit.Chem import AllChem as _AllChem
from rdkit.Chem import rdForceFieldHelpers as _FF
from rdkit.Chem import rdFMCS as _FMCS
import parmed as _pmd

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
    if isinstance(val, str) and len(val.split(".")) == 1:
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
        try:
            mol = openFileAsRdkit(val, **kwargs)
        except:
            raise ValueError("File is not in a valid format")

    return mol

def saveFromRdkit(mol, filename):
    extension = filename.split(".")[-1]
    if extension.lower() == "sdf":
        writer = _Chem.SDWriter(filename)
        writer.save(mol)
        writer.close()
    elif extension.lower() == "mol":
        _Chem.MolToMolFile(mol, filename=filename)
    else:
        tempfilename = _tempfile.NamedTemporaryFile(suffix=".mol").name
        _Chem.MolToMolFile(mol, filename=tempfilename)
        _pmd.load_file(tempfilename)[0].save(filename)

    return filename

def alignTwoMolecules(mol, ref, n_min=-1, match="any"):
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

    mcs = _Chem.MolFromSmarts(_FMCS.FindMCS([ref, mol], atomCompare=comp_method,
                                            ringMatchesRingOnly=True, completeRingsOnly=True).smartsString)
    match1 = ref.GetSubstructMatch(mcs)
    match2 = mol.GetSubstructMatch(mcs)

    ref_conf = ref.GetConformer(-1)
    mol_conf = mol.GetConformer(-1)

    ff = _FF.MMFFGetMoleculeForceField(mol, _FF.MMFFGetMoleculeProperties(mol), confId=0)
    for i_ref, i_mol in zip(match1, match2):
        mol_conf.SetAtomPosition(i_mol, ref_conf.GetAtomPosition(i_ref))
        ff.AddFixedPoint(i_mol)

    if mol_conf.GetNumAtoms() != len(match1):
        ff.Initialize()
        more = ff.Minimize()
        while more and n_min:
            more = ff.Minimize()
            #n_min = -1 means infinite minimisation
            if n_min != -1:
                n_min -= 1
    return mol