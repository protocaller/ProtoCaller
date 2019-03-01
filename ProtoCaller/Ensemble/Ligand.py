import os as _os

import rdkit.Chem.rdchem as _rdchem
import rdkit.Chem.rdmolfiles as _rdmolfiles
import rdkit.Chem.rdmolops as _rdmolops

import ProtoCaller.Parametrise as _parametrise
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Wrappers.babelwrapper as _babel
import ProtoCaller.Wrappers.rdkitwrapper as _rdkit


class Ligand:
    """
    The main class responsible for handling ligands and cofactors in ProtoCaller.

    Parameters
    ----------
    input: str or rdkit.Chem.rdChem.Mol
        Initialises the ligand from a SMILES string, InChI string, file or an RDKit Mol object.
    parametrised_files: [str]
        Initialises the parametrised file(s).
    name: str or None
        Initialises the ligand name.
    protonated : bool
        Whether the input corresponds to the protonated or deprotonated state of the molecule.
    minimise : bool or None
        Initialises minimise.
    workdir : str, optional
        Initialises workdir.

    Attributes
    ----------
    name : str
        The ligand name. Default: e.g. "ligand1" for the first instantiated Ligand object.
    workdir : ProtoCaller.Utils.fileio.Dir
        The working directory of the ligand.
    minimise : bool
        Whether to perform a GAFF minimisation using OpenBabel. None means minimisation for molecules initialised from
        strings and no minimisation for molecules initialised from files.
    """
    _counter = 1

    def __init__(self, input, parametrised_files=None, name=None, protonated=False, minimise=None, workdir="."):
        self.name = name
        self.workdir = _fileio.Dir(workdir)
        self.minimise = minimise
        # always set protonated to False and if a valid protonated file is given it is automatically set to True
        self.protonated_filename = None

        with self.workdir:
            if isinstance(input, str):
                if _os.path.exists(input):
                    if protonated:
                        self.protonated_filename = input
                    else:
                        self.molecule = _rdkit.openAsRdkit(input, minimise=minimise)
                else:
                    self.string = input
            elif isinstance(input, _rdchem.Mol):
                self.molecule = input
            else:
                raise TypeError("Need a SMILES, InChI string, filename or an RDKit object as an input")
            self.parametrised_files = parametrised_files
            self.minimise = False

    def __hash__(self):
        return hash(self.name)

    @property
    def name(self):
        """The name of the ligand"""
        return self._name

    @name.setter
    def name(self, val):
        if val is None:
            self._name = "ligand%d" % self._counter
            Ligand._counter += 1
        else:
            self._name = val

    @property
    def string(self):
        """The SMILES string of the ligand."""
        return self._string

    @string.setter
    def string(self, val):
        if not _os.path.exists(val):
            self._molecule = _rdkit.openAsRdkit(val, minimise=self.minimise)
            self._string = _rdmolfiles.MolToSmiles(self._molecule)
        else:
            raise ValueError("Need a SMILES or InChI string instead of a filename")

    @property
    def molecule(self):
        """The rdkit.Chem.rdChem.Mol object of the ligand."""
        return self._molecule

    @molecule.setter
    def molecule(self, val):
        if isinstance(val, _rdchem.Mol):
            self._molecule = val
            self._string = _rdmolfiles.MolToSmiles(self._molecule)
        else:
            raise TypeError("Need an object of type RDKit Mol")

    @property
    def protonated(self):
        """Whether the ligand has been protonated."""
        return self._protonated

    @property
    def protonated_filename(self):
        """The absolute path to the protonated file."""
        return self._protonated_filename

    @protonated_filename.setter
    def protonated_filename(self, val):
        with self.workdir:
            if val is None:
                self._protonated_filename = None
                self._protonated = False
            else:
                self._protonated_filename = _fileio.checkFileExists(val)
                self._molecule = _rdkit.openAsRdkit(self._protonated_filename, removeHs=False, minimise=self.minimise)
                self._string = _rdmolfiles.MolToSmiles(self.molecule)
                self._protonated = True

    @property
    def parametrised(self):
        """Whether the ligand has been parametrised."""
        return self._parametrised

    @property
    def parametrised_files(self):
        """The absolute path(s) to the parametrised file(s)."""
        return self._parametrised_files

    @parametrised_files.setter
    def parametrised_files(self, val):
        with self.workdir:
            if val is None:
                self._parametrised_files = None
                self._parametrised = False
            else:
                for i, file in enumerate(val):
                    val[i] = _fileio.checkFileExists(file)
                self._parametrised_files = val
                self._parametrised = True

    def protonate(self, reprotonate=False, babel_parameters=None, rdkit_parameters=None):
        """
        Protonates the ligand using OpenBabel.

        Parameters
        ----------
        reprotonate : bool
            Whether to reprotonate an already protonated ligand.
        babel_parameters : dict
            Keyword arguments to be passed to ProtoCaller.Wrappers.babelwrapper
        rdkit_parameters : dict
            Keyword arguments to be passed to ProtoCaller.Wrappers.rdkitwrapper
        """
        with self.workdir:
            if babel_parameters is None: babel_parameters = {}
            babel_parameters = {"pH": 7.0, **babel_parameters}
            if rdkit_parameters is None: rdkit_parameters = {}

            if self.protonated and not reprotonate:
                print("Ligand %s is already protonated." % self.name)
            else:
                filename_temp = _rdkit.saveFromRdkit(self.molecule, filename="%s.mol" % self.name)
                # here we use SDF because of parser differences between OpenBabel and RDKit concerning mol and mol2 files
                self.protonated_filename = _babel.babelTransform(filename_temp, "sdf", **babel_parameters)
                _os.remove(filename_temp)
                self.molecule = _rdkit.openFileAsRdkit(self.protonated_filename, removeHs=False, **rdkit_parameters)

    def parametrise(self, params=None, molecule_type="ligand", id=None, reparametrise=False):
        """
        Parametrises the ligand using ProtoCaller.Parametrise.

        Parameters
        ----------
        params : ProtoCaller.Parametrise.Params
            Force field parameters.
        molecule_type : str
            The type of the molecule. One of: "ligand" and "cofactor".
        id : str
            The name of the molecule. Default: equal to the ligand name.
        reparametrise : bool
            Whether to reparametrise an already parametrised ligand.
        """
        with self.workdir:
            print("Parametrising ligand %s..." % self.name)
            if self._parametrised and not reparametrise:
                print("Ligand %s is already parametrised." % self.name)
                return

            if not self.protonated:
                print("Cannot parametrise unprotonated ligand. Protonating first with default parameters...")
                self.protonate()

            if params is None: params = _parametrise.Params()
            # we convert the protonated file into a pdb so that antechamber can read it
            filename = _babel.babelTransform(self.protonated_filename, "pdb")
            if id is None: id = self.name

            charge = _rdmolops.GetFormalCharge(self.molecule)
            self.parametrised_files = _parametrise.parametriseFile(params=params, filename=filename,
                                                                   molecule_type=molecule_type, id=id, charge=charge)
