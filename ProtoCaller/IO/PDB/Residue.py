from Bio.Data.IUPACData import protein_letters_3to1 as _3_to_1

import ProtoCaller.Utils.ConditionalList as _CondList
from . import Missing as _Missing


class Residue(_Missing.MissingResidue, _CondList.ConditionalList):
    """
    Represents a single residue in a PDB file.

    Parameters
    ----------
    atoms : [ProtoCaller.IO.PDB.Atom.Atom]
        A list of atoms.

    Attributes
    ----------
    resName : str
        Residue name.
    chainID : str
        Chain identifier.
    resSeq : int
        Residue sequence number.
    iCode : str
        Code for insertion of residues.
    """
    # properties which have to be conserved within the whole residue
    _common_properties = ["chainID", "resName", "resSeq", "iCode"]

    def __init__(self, atoms=None):
        if atoms is None:
            atoms = []
        _CondList.ConditionalList.__init__(self, atoms, self._checkAtom)

    def __getattr__(self, item):
        if item in self._common_properties:
            if not len(self):
                return None
            return getattr(self[0], item)
        else:
            raise AttributeError("Invalid attribute: {}. Needs to be one of: {}".format(item, self._common_properties))

    def __setattr__(self, key, value):
        if key in self._common_properties:
            for atom in self:
                setattr(atom, key, value)
            _Missing.MissingResidue.__setattr__(self, key, value)
        else:
            _CondList.ConditionalList.__setattr__(self, key, value)

    def __str__(self):
        return "".join([atom.__str__() for atom in self])

    def __repr__(self):
        return "<Residue of {} atoms>".format(len(self))

    @property
    def numberOfAtoms(self):
        """int: Returns the number of atoms."""
        return len(self)

    @property
    def sequence(self):
        """str: Returns the single-character code of the amino acid residue."""
        try:
            code = _3_to_1[self.resName.capitalize()]
        except:
            code = "-"
        return code

    def reNumberAtoms(self, start=1, custom_serials=None):
        """
        Renumbers the atoms.

        Parameters
        ----------
        start : int, optional
            Initial index.
        custom_serials : [int] or None, optional
            A list of custom indices. If None, atoms are renumbered sequentially.
        """
        if custom_serials is None:
            serials = [start + i for i in range(start, start + len(self))]
        elif len(custom_serials) == len(self):
            serials = custom_serials
        else:
            raise ValueError("Custom number list does not match number of atoms")

        for i, atom in zip(serials, self):
            atom.serial = i

    def purgeAtoms(self, atoms, mode):
        """
        Removes or keeps the selected atoms.

        Parameters
        ----------
        atoms : [ProtoCaller.IO.PDB.Atom.Atom]
            Atoms to be removed or kept
        mode : str
            One of "keep" and "discard" - either keeps only the selection or discards the latter.
        """
        assert mode in ["keep", "discard"]
        atoms_to_remove = atoms if mode == "discard" else [x for x in self if x not in atoms]
        for atom in atoms_to_remove:
            self.remove(atom)

    def _checkAtom(self, atom):
        """Checks whether the atom has the same chainID, resName, resSeq and iCode as the current object."""
        try:
            for prop in self._common_properties:
                if getattr(self, prop) is not None:
                    assert (getattr(self, prop) == getattr(atom, prop))
        except AttributeError:
            raise TypeError("Need to pass a valid object with chainID, resSeq, resName and iCode attributes")
        except AssertionError:
            raise ValueError("Input atom(s) from an incompatible residue")
