import copy as _copy
import logging as _logging
import os as _os

# TODO support the native Morph class when it becomes viable
import BioSimSpace as _BSS

from .Ligand import Ligand as _Ligand
import ProtoCaller.Wrappers.rdkitwrapper as _rdkit


class Perturbation:
    """
    The main class responsible for single-topology perturbations in ProtoCaller.

    Parameters
    ----------
    ligand1 : ProtoCaller.Ensemble.Ligand.Ligand
        Initialises ligand1.
    ligand2 : ProtoCaller.Ensemble.Ligand.Ligand
        Initialises ligand1.
    name : str or None, optional
        Initialises name.

    Attributes
    ----------
    ligand1 : ProtoCaller.Ensemble.Ligand.Ligand
        The first ligand.
    ligand2 : ProtoCaller.Ensemble.Ligand.Ligand
        The second ligand.
    current_ref : ProtoCaller.Ensemble.Ligand.Ligand
        The current reference ligand.
    """
    _counter = 1

    def __init__(self, ligand1, ligand2, name=None):
        if not isinstance(ligand1, _Ligand) and not isinstance(ligand2, _Ligand):
            raise TypeError("Input ligands need to be of type Ligand")

        self.ligand1 = ligand1
        self.ligand2 = ligand2
        self.current_ref = None
        self._ligand1_molecule = {}
        self._ligand2_molecule = {}
        self._ligand1_coords = {}
        self._ligand2_coords = {}
        self._morph = {}
        self.name = name

    @property
    def name(self):
        """The name of the perturbation."""
        return self._name

    @name.setter
    def name(self, val):
        if val is True:
            self._name = "Pair %d" % self._counter
            Perturbation._counter += 1
        elif val is None:
            self._name = "{}~{}".format(self.ligand1.name, self.ligand2.name)
        else:
            self._name = val

    def get_parametrised_files1(self, ref):
        """
        Returns parametrised files for the first ligand corresponding to a reference ligand.

        Parameters
        ----------
        ref : ProtoCaller.Ensemble.Ligand.Ligand
            The reference ligand.

        Returns
        -------
        files : [str, str]
            The parametrised files for the first ligand aligned to the reference ligand. The first file is a
            topology file and the second file is a coordinate file.
        """
        topol = None if not self.ligand1.parametrised else self.ligand1.parametrised_files[0]
        if self._ligand1_coords is not None:
            coords = self._ligand1_coords[ref]
        else:
            coords = None if not self.ligand1.parametrised else self.ligand1.parametrised_files[1]

        if topol is None or coords is None:
            return None
        else:
            return [topol, coords]

    @property
    def parametrised_files1(self):
        """The parametrised files for the first ligand corresponding to the current reference ligand."""
        return self.get_parametrised_files1(ref=self.current_ref)

    def get_parametrised_files2(self, ref):
        """
        Returns parametrised files for the second ligand corresponding to a reference ligand.

        Parameters
        ----------
        ref : ProtoCaller.Ensemble.Ligand.Ligand
            The reference ligand.

        Returns
        -------
        files : [str, str]
            The parametrised files for the first second aligned to the reference ligand. The first file is a
            topology file and the second file is a coordinate file.
        """
        topol = None if not self.ligand2.parametrised else self.ligand2.parametrised_files[0]
        if self._ligand2_coords is not None:
            coords = self._ligand2_coords[ref]
        else:
            coords = None if not self.ligand2.parametrised else self.ligand2.parametrised_files[1]

        if topol is None or coords is None:
            return None
        else:
            return [topol, coords]

    @property
    def parametrised_files2(self):
        """The parametrised files for the second ligand corresponding to the current reference ligand."""
        return self.get_parametrised_files2(ref=self.current_ref)

    def get_morph(self, ref):
        """
        Returns the mixed topology object corresponding to a reference ligand.

        Parameters
        ----------
        ref : ProtoCaller.Ensemble.Ligand.Ligand
            The reference ligand.


        Returns
        -------
        morph : BioSimSpace.System
            The mixed topology object aligned to the reference ligand.
        """
        return self._morph[ref]

    @property
    def morph(self):
        """The mixed topology object corresponding to the current reference ligand."""
        return self._morph[self.current_ref]

    def isAlignedTo(self, ref):
        """
        Returns whether there is a mixed topology file that has been aligned to the reference ligand.

        Parameters
        ----------
        ref : ProtoCaller.Ensemble.Ligand.Ligand
            The reference ligand.

        Returns
        -------
        exists : bool
            Whether the mixed topology file has already been aligned to the reference ligand.
        """
        return ref in self._morph.keys()

    @property
    def isAligned(self):
        """Whether there is a mixed topology file that has been aligned to the current reference ligand."""
        return self.isAlignedTo(self.current_ref)

    def alignToReference(self, ref, output_filename=None, realign=False, **kwargs):
        """
        Aligns the first ligand to a reference ligand.

        Parameters
        ----------
        ref : ProtoCaller.Ensemble.Ligand.Ligand
            The reference ligand.
        output_filename : str
            The name of the output file. None uses the default value.
        realign : bool
            Whether to realign an already aligned ligand.
        kwargs
            Keyword arguments to pass to ProtoCaller.Wrappers.rdkitwrapper.alignTwoMolecules.

        Returns
        -------
        mcs : [tuple]
            The maximum common substructure of the two molecules.
        """
        self.current_ref = ref
        if ref in self._ligand1_molecule.keys() and ref in self._ligand1_coords.keys() and not realign:
            _logging.info("Morph %s is already aligned to a reference" % self.name)
            return

        if output_filename is None:
            self.ligand1.parametrise(reparametrise=False)
            output_filename = "{}_aligned_{}.{}".format(self.ligand1.name, self.name,
                                                        self.ligand1.parametrised_files[1].split(".")[-1])
        if "two_way_matching" not in kwargs.keys():
            kwargs["two_way_matching"] = False

        lig_temp = _copy.deepcopy(self.ligand1.molecule)
        self._ligand1_molecule[ref], mcs = _rdkit.alignTwoMolecules(ref.molecule, lig_temp, **kwargs)
        self._ligand1_coords[ref] = _rdkit.saveFromRdkit(self._ligand1_molecule[ref], output_filename)

        return mcs

    def alignToEachOther(self, output_filename=None, realign=False, **kwargs):
        """
        Aligns the second ligand to the first ligand.

        Parameters
        ----------
        output_filename : str
            The name of the output file. None uses the default value.
        realign : bool
            Whether to realign an already aligned ligand.
        kwargs
            Keyword arguments to pass to ProtoCaller.Wrappers.rdkitwrapper.alignTwoMolecules.

        Returns
        -------
        morph : BioSimSpace.System
            The mixed topology object.
        mcs : [tuple]
            The maximum common substructure of the two molecules.
        """
        if self.current_ref in self._morph.keys() and not realign:
            _logging.info("Morph %s is already aligned to a reference" % self.name)
            return

        if output_filename is None:
            self.ligand1.parametrise(reparametrise=False)
            self.ligand2.parametrise(reparametrise=False)
            output_filename = "{}_aligned_{}.{}".format(self.ligand2.name, self.name,
                                                        self.ligand2.parametrised_files[1].split(".")[-1])
        if "two_way_matching" not in kwargs.keys():
            kwargs["two_way_matching"] = True

        if self.current_ref is None:
            self._ligand1_molecule[None] = self.ligand1.molecule
            self._ligand1_coords[None] = self.ligand1.parametrised_files[1]

        lig1_temp = self._ligand1_molecule[self.current_ref]
        lig2_temp = _copy.deepcopy(self.ligand2.molecule)
        self._ligand2_molecule[self.current_ref], mcs = _rdkit.alignTwoMolecules(lig1_temp, lig2_temp, **kwargs)
        self._ligand2_coords[self.current_ref] = _os.path.abspath(
            _rdkit.saveFromRdkit(self._ligand2_molecule[self.current_ref], output_filename))

        # load the shifted ligands into BioSimSpace
        ligand1_BSS = _BSS.IO.readMolecules(self.parametrised_files1).getMolecules()[0]
        ligand2_BSS = _BSS.IO.readMolecules(self.parametrised_files2).getMolecules()[0]

        # translate RDKit mapping into BioSimSpace mapping
        mapping = {}
        indices_1 = [atom.index() for atom in ligand1_BSS._sire_object.edit().atoms()]
        indices_2 = [atom.index() for atom in ligand2_BSS._sire_object.edit().atoms()]
        for idx1, idx2 in mcs:
            mapping[indices_1[idx1]] = indices_2[idx2]

        # finally create morph
        self._morph[self.current_ref] = _BSS.Align.merge(
            ligand1_BSS, ligand2_BSS, mapping=mapping, allow_ring_breaking=True)

        return self._morph[self.current_ref], mcs

    def alignAndCreateMorph(self, ref):
        """
        Default alignment method wrapping the two alignment stages - aligning ligand1 to a reference and aligning and
        mixing ligand2 to ligand1.

        Parameters
        ----------
        ref : ProtoCaller.Ensemble.Ligand.Ligand
            The reference ligand.

        Returns
        -------
        morph : BioSimSpace.System
            The mixed topology object object.
        mcs : [tuple]
            The maximum common substructure of the two molecules.
        """
        self.alignToReference(ref)
        return self.alignToEachOther()
