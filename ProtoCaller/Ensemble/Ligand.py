import os as _os

import rdkit.Chem.rdmolfiles as _rdmolfiles
import rdkit.Chem.rdchem as _rdchem

import ProtoCaller.Parametrise as _parametrise
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Wrappers.rdkitwrapper as _rdkit
import ProtoCaller.Wrappers.babelwrapper as _babel


class Ligand:
    _counter = 1

    def __init__(self, input, parametrised_files=None, name=None, protonated=False):
        self.name = name
        if isinstance(input, str):
            if _os.path.exists(input):
                if protonated:
                    self.protonated_filename = input
                else:
                    self.molecule = _rdkit.openAsRdkit(input)
            else:
                protonated = False
                self.string = input
        elif isinstance(input, _rdchem.Mol):
            self.molecule = input
        self.parametrised_files = parametrised_files
        self.protonated = protonated

    def __hash__(self):
        return self.string

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, val):
        if val is None:
            self._name = "ligand%d" % self._counter
        else:
            self._name = val
        Ligand._counter += 1

    @property
    def string(self):
        return self._string

    @string.setter
    def string(self, val):
        if not _os.path.exists(val):
            self._molecule = _rdkit.openAsRdkit(val)
            self._string = _rdmolfiles.MolToSmiles(self._molecule)
        else:
            raise ValueError("Need a SMILES or InChI string instead of a filename")

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, val):
        if isinstance(val, _rdchem.Mol):
            self._molecule = val
            self._string = _rdmolfiles.MolToSmiles(self._molecule)
        else:
            raise TypeError("Need an object of type RDKit Mol")

    @property
    def protonated_filename(self):
        return self._protonated_filename

    @protonated_filename.setter
    def protonated_filename(self, val):
        if val is None:
            self.protonated = None
            self._protonated_filename = None
        else:
            self._protonated_filename = _fileio.checkFileExists(val)
            self._molecule = _rdkit.openAsRdkit(self._protonated_filename, removeHs=False)
            self._string = _rdmolfiles.MolToSmiles(self.molecule)

    @property
    def parametrised_files(self):
        return self._parametrised_files

    @parametrised_files.setter
    def parametrised_files(self, val):
        if val is None:
            self._parametrised_files = None
            self.parametrised = False
        else:
            for i, file in enumerate(val):
                val[i] = _fileio.checkFileExists(file)
            self._parametrised_files = val
            self.parametrised = True

    def protonate(self, reprotonate=False, babel_parameters=None, rdkit_parameters=None):
        if babel_parameters is None: babel_parameters = {}
        if rdkit_parameters is None: rdkit_parameters = {}

        if self.protonated and not reprotonate:
            print("Ligand %s is already protonated." % self.name)
        else:
            filename_temp = _rdkit.saveFromRdkit(self.molecule, filename="%s.mol" % self.name)
            self.protonated_filename = _babel.babelTransform(filename_temp, "mol2", **babel_parameters)
            _os.remove(filename_temp)
            self.molecule = _rdkit.openFileAsRdkit(self.protonated_filename, removeHs=False, **rdkit_parameters)
            self.protonated = True

    def parametrise(self, params, filename=None, molecule_type="ligand", id=None, reparametrise=False):
        if filename is None: filename=self.protonated_filename
        if id is None: id = self.name
        print("Parametrising ligand %s..." % self.name)
        if self.parametrised and not reparametrise:
            print("Ligand %s is already parametrised." % self.name)
            return

        if not self.protonated:
            print("Cannot parametrise unprotonated ligand. Protonating first with default parameters...")
            self.protonate()

        self.parametrised_files = _parametrise.parametriseFile(params=params, filename=filename,
                                                               molecule_type=molecule_type, id=id)
        self.parametrised = True
