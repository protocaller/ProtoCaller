import copy as _copy
import os as _os

# TODO support the native Morph class when it becomes viable
import BioSimSpace as _BSS

from .Ligand import Ligand as _Ligand
import ProtoCaller.Wrappers.rdkitwrapper as _rdkit


class Perturbation:
    _counter = 1

    def __init__(self, ligand1, ligand2, name=None):
        try:
            assert isinstance(ligand1, _Ligand) and isinstance(ligand2, _Ligand)
        except AssertionError:
            raise TypeError("Input ligands need to be of type Ligand")

        self._ligand1 = ligand1
        self._ligand2 = ligand2
        self._ligand1_molecule = None
        self._ligand2_molecule = None
        self._ligand1_coords = None
        self._ligand2_coords = None
        self._morph = None
        self.name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, val):
        if val is True:
            self._name = "Pair %d" % self._counter
            Perturbation._counter += 1
        elif val is None:
            self._name = "{}~{}".format(self._ligand1.name, self._ligand2.name)
        else:
            self._name = val

    @property
    def parametrised_files1(self):
        topol = None if not self._ligand1.parametrised else self._ligand1.parametrised_files[0]
        if self._ligand1_coords is not None:
            coords = self._ligand1_coords
        else:
            coords = None if not self._ligand1.parametrised else self._ligand1.parametrised_files[1]

        if topol is None or coords is None:
            return None
        else:
            return [topol, coords]

    @property
    def parametrised_files2(self):
        topol = None if not self._ligand2.parametrised else self._ligand2.parametrised_files[0]
        if self._ligand2_coords is not None:
            coords = self._ligand2_coords
        else:
            coords = None if not self._ligand2.parametrised else self._ligand2.parametrised_files[1]

        if topol is None or coords is None:
            return None
        else:
            return [topol, coords]

    @property
    def morph(self):
        return self._morph

    def alignToReference(self, ref, output_filename=None, **kwargs):
        if output_filename is None:
            if not self._ligand1.parametrised:
                self._ligand1.parametrise()
            output_filename = "{}_aligned_{}.{}".format(self._ligand1.name, self.name,
                                                        self._ligand1.parametrised_files[1].split(".")[-1])
        if "always_maximum" not in kwargs.keys(): kwargs["always_maximum"] = False

        ref_temp = _copy.deepcopy(ref.molecule)
        lig_temp = _copy.deepcopy(self._ligand1.molecule)
        self._ligand1_molecule, mcs = _rdkit.alignTwoMolecules(ref_temp, lig_temp, **kwargs)
        self._ligand1_coords = _os.path.abspath(_rdkit.saveFromRdkit(self._ligand1_molecule, output_filename))

        return mcs

    def alignToEachOther(self, output_filename=None, **kwargs):
        if output_filename is None:
            if not self._ligand1.parametrised:
                self._ligand1.parametrise()
            if not self._ligand2.parametrised:
                self._ligand2.parametrise()
            output_filename = "{}_aligned_{}.{}".format(self._ligand2.name, self.name,
                                                        self._ligand2.parametrised_files[1].split(".")[-1])
        if "always_maximum" not in kwargs.keys(): kwargs["always_maximum"] = True

        lig1_temp = self._ligand1_molecule if self._ligand1_molecule is not None else self._ligand1.molecule
        lig2_temp = _copy.deepcopy(self._ligand2.molecule)
        self._ligand2_molecule, mcs = _rdkit.alignTwoMolecules(lig1_temp, lig2_temp, **kwargs)
        self._ligand2_coords = _os.path.abspath(_rdkit.saveFromRdkit(self._ligand2_molecule, output_filename))

        # load the shifted ligands into BioSimSpace
        ligand1_BSS = _BSS.IO.readMolecules(self.parametrised_files1).getMolecules()[0]
        ligand2_BSS = _BSS.IO.readMolecules(self.parametrised_files2).getMolecules()[0]

        # translate RDKit mapping into BioSimSpace mapping
        mapping = {}
        indices_1 = [atom.index() for atom in ligand1_BSS._sire_molecule.edit().atoms()]
        indices_2 = [atom.index() for atom in ligand2_BSS._sire_molecule.edit().atoms()]
        for idx1, idx2 in mcs:
            mapping[indices_1[idx1]] = indices_2[idx2]

        # finally create morph
        self._morph = _BSS.Align.merge(ligand1_BSS, ligand2_BSS, mapping=mapping)

        return self._morph, mcs

    # default alignment method wrapping the two alignment stages
    def alignAndCreateMorph(self, ref):
        self.alignToReference(ref)
        return self.alignToEachOther()
