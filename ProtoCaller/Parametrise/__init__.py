import ProtoCaller as _PC

if _PC.BIOSIMSPACE:
    import BioSimSpace as _BSS

import ProtoCaller.Wrappers.parmedwrapper as _pmdwrap
from . import amber


class Params:
    """
    A simple class which contains force field names for different molecule types.

    Parameters
    ----------
    protein_ff: str
        Initialises the protein force field.
    ligand_ff: str
        Initialises the ligand force field.
    water_ff: str
        Initialises the water / simple ion force field.
    """
    def __init__(self, protein_ff=_PC.AMBERDEFAULTPROTEINFF, ligand_ff=_PC.AMBERDEFAULTLIGANDFF,
                 water_ff=_PC.AMBERDEFAULTWATERFF):
        self.protein_ff = protein_ff
        self.ligand_ff = ligand_ff
        self.water_ff = water_ff

    @property
    def protein_ff(self):
        """The protein force field."""
        return self._protein_ff

    @protein_ff.setter
    def protein_ff(self, val):
        val = val.strip()
        if val not in _PC.PROTEINFFS:
            raise ValueError("Value %s not supported. Supported values: " % val, _PC.PROTEINFFS)
        self._protein_ff = val

    @property
    def ligand_ff(self):
        """The ligand force field."""
        return self._ligand_ff

    @ligand_ff.setter
    def ligand_ff(self, val):
        val = val.strip()
        if val not in _PC.LIGANDFFS:
            raise ValueError("Value %s not supported. Supported values: " % val, _PC.LIGANDFFS)
        self._ligand_ff = val

    @property
    def water_ff(self):
        """The water / simple ion force field."""
        return self._water_ff

    @water_ff.setter
    def water_ff(self, val):
        val = val.strip()
        if val not in _PC.WATERFFS:
            raise ValueError("Value %s not supported. Supported values: " % val, _PC.WATERFFS)
        self._water_ff = val


def parametriseAndLoadPmd(params, *args, **kwargs):
    """
    Calls parametriseFile and loads the files into ParmEd.

    Parameters
    ----------
    params : ProtoCaller.Parametrise.Params
        Force field parameters.
    args
        Positional arguments to be passed to parametriseFile.
    kwargs
        Keyword arguments to be passed to parametriseFile.

    Returns
    -------
    system : parmed.structure.Structure
        The loaded parametrised files as a ParmEd Structure.
    """
    files = parametriseFile(params, *args, **kwargs)
    return _pmdwrap.openFilesAsParmed(files)


def parametriseAndLoadBSS(params, *args, **kwargs):
    """
    Calls parametriseFile and loads the files into BioSimSpace.

    Parameters
    ----------
    params : ProtoCaller.Parametrise.Params
        Force field parameters.
    args
        Positional arguments to be passed to parametriseFile.
    kwargs
        Keyword arguments to be passed to parametriseFile.

    Returns
    -------
    system : BioSimSpace.system
        The loaded parametrised files as a BioSimSpace System.
    """
    if _PC.BIOSIMSPACE:
        files = parametriseFile(params, *args, **kwargs)
        return _BSS.IO.readMolecules(files)
    else:
        raise ImportError("BioSimSpace module cannot be imported")


def parametriseFile(params, filename, molecule_type, fix_charge=True, id=None, *args, **kwargs):
    """
    Parametrises an input file according to a force field.

    Parameters
    ----------
    params : ProtoCaller.Parametrise.Params
        Force field parameters.
    filename : str
        Name of the input file.
    molecule_type : str
        The type of the molecule. One of: "protein", "ligand", "cofactor", "water", "simple_anion", "complex_anion",
        "simple_cation", "complex_cation".
    fix_charge : bool
        Whether to enforce integer net charge.
    id : str
        The name of the molecule. Default: equal to molecule_type.
    args
        Positional arguments to be passed to the relevant wrapper.
    kwargs
        Keyword arguments to be passed to the relevant wrapper.

    Returns
    -------
    files : [str]
        The output parametrised file(s). If there is more than one file, the first one is always the topology file and
        the second one - the coordinate file.
    """
    if id is None: id = molecule_type
    files = []

    if molecule_type == "protein":
        if params.protein_ff in _PC.AMBERPROTEINFFS:
            files = amber.amberWrapper(params, filename, molecule_type, id, *args, **kwargs)
    elif molecule_type in ["ligand", "cofactor"]:
        if params.ligand_ff in _PC.AMBERLIGANDFFS:
            files = amber.amberWrapper(params, filename, molecule_type, id, *args, **kwargs)
    elif molecule_type in ["water", "simple_anion", "complex_anion", "simple_cation", "complex_cation"]:
        if params.water_ff in _PC.AMBERWATERFFS:
            files = amber.amberWrapper(params, filename, molecule_type, id, *args, **kwargs)
    else:
        raise ValueError("Invalid argument: %s. Argument must be one of: protein, ligand, cofactor, simple_anion, "
                         "complex_anion, simple_cation, complex_cation or water." % molecule_type)

    if files:
        return _pmdwrap.fixCharge(files) if fix_charge else files
    else:
        raise ValueError("No force fields available for the input system")
