import logging as _logging
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
        Ignored when parametrised files are passed.
    minimise : bool or None
        Initialises minimise. Ignored when parametrised files are passed.
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
            if parametrised_files:
                self.molecule = _rdkit.openAsRdkit(parametrised_files, removeHs=False, minimise=False,
                                                   template=input)
                self.protonated_filename = _rdkit.saveFromRdkit(self.molecule, "{}.pdb".format(name))
            elif isinstance(input, str):
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
            self._molecule = _rdmolops.RemoveHs(self._molecule, updateExplicitCount=True)
            self._molecule = _rdmolops.AddHs(self._molecule, addCoords=True)
            self._string = _rdmolfiles.MolToSmiles(self._molecule, allHsExplicit=True)
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
            self._string = _rdmolfiles.MolToSmiles(self._molecule, allHsExplicit=True)
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
                self._string = _rdmolfiles.MolToSmiles(self.molecule, allHsExplicit=True)
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
                val = _fileio.checkFileExists(val)
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
                _logging.info("Ligand %s is already protonated." % self.name)
            else:
                # Here we are careful about reading bond orders generated by Open Babel
                filename_temp = _rdkit.saveFromRdkit(self.molecule, filename="%s_protonated.sdf" % self.name)
                protonated_filename = _babel.babelTransform(filename_temp, "mol", **babel_parameters)
                molecule = _rdkit.openAsRdkit(protonated_filename, removeHs=False, sanitize=False)
                molecule = _rdkit.AssignBondOrdersFromTemplate(self.molecule, molecule, assign_charge=False)
                _os.remove(protonated_filename)
                self.protonated_filename = _rdkit.saveFromRdkit(molecule, filename_temp)

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
            if self._parametrised and not reparametrise:
                _logging.debug("Ligand %s is already parametrised." % self.name)
                return

            _logging.info("Parametrising ligand %s..." % self.name)
            if not self.protonated:
                _logging.warning("Cannot parametrise unprotonated ligand. Protonating first with default parameters...")
                self.protonate()

            if params is None:
                params = _parametrise.Params()

            # we convert the protonated file into a pdb so that antechamber can read it
            filename = _babel.babelTransform(self.protonated_filename, "pdb")
            if id is None: id = self.name

            charge = _rdmolops.GetFormalCharge(self.molecule)
            self.parametrised_files = _parametrise.parametriseFile(params=params, filename=filename,
                                                                   molecule_type=molecule_type, id=id, charge=charge)
