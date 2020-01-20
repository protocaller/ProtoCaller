import ProtoCaller.Utils.ConditionalList as _CondList
from . import Residue as _Residue


class Chain(_Residue, _CondList.ConditionalList):
    """
    Represents a single Chain in a PDB file.

    Parameters
    ----------
    atoms : [ProtoCaller.IO.PDB.Residue.Residue]
        A list of residues.

    Attributes
    ----------
    chainID : str
        Chain identifier.
    """
    # properties which have to be conserved within the whole chain
    _common_properties = ["chainID"]

    def __init__(self, residues=None):
        if residues is None:
            residues = []
        _CondList.ConditionalList.__init__(self, residues, self._checkResidue)

    def __eq__(self, other):
        return self.sameChain(other)

    def __hash__(self):
        return id(self)

    def __repr__(self):
        return "<Chain of {} residues>".format(len(self))

    @property
    def sequence(self):
        """str: Returns the single-character sequence of the amino acids in
                the chain. Missing amino acids are not taken into account."""
        try:
            seq = "".join([residue.sequence for residue in self
                           if residue.type in ["amino_acid",
                                               "amino_acid_modified"]])
        except:
            seq = ""
        return seq

    @property
    def type(self):
        """str: Returns the type of the chain: "chain" if it has amino acid
                residues and "molecules" otherwise."""
        if not len(self):
            return None
        else:
            for residue in self:
                if residue.type in ["amino_acid"]:
                    return "chain"
            return "molecules"

    @property
    def numberOfAtoms(self):
        return sum([residue.numberOfAtoms for residue in self])

    @property
    def numberOfResidues(self):
        """int: Returns the number of Residues."""
        return len(self)

    def reNumberAtoms(self, start=1, custom_serials=None):
        if custom_serials is None:
            custom_serials = [start + i for i in range(self.numberOfAtoms)]
        if not len(custom_serials) == self.numberOfAtoms:
            raise ValueError("Custom number of atoms does not match chain number of atoms")

        i = 0
        for residue in self:
            residue.reNumberAtoms(custom_serials=custom_serials[i:i + residue.numberOfAtoms])
            i += residue.numberOfAtoms

    def reNumberResidues(self, start=1, custom_resSeqs=None, custom_iCodes=None):
        """
        Renumbers the residues.

        Parameters
        ----------
        start : int, optional
            Initial index.
        custom_resSeqs : [int] or None, optional
            A list of custom resSeqs. If None, residues are renumbered sequentially.
        custom_iCodes : [str] or None, optional
            A list of custom iCodes. If None, residues are renumbered sequentially.
        """
        if custom_resSeqs is None:
            custom_resSeqs = [start + i for i in range(len(self))]
        if custom_iCodes is None:
            custom_iCodes = [" "] * len(custom_resSeqs)
        if not len(custom_resSeqs) == len(custom_iCodes) == len(self):
            raise ValueError("Custom number of residues does not match chain number of residues")

        for resSeq, iCode, residue in zip(custom_resSeqs, custom_iCodes, self):
            residue.resSeq, residue.iCode = resSeq, iCode

    def purgeAtoms(self, atoms, mode):
        for residue in self:
            residue.purgeAtoms(atoms, mode)
        self.purgeEmpty()

    def purgeEmpty(self):
        """Removes empty elements."""
        empty_residues = [res for res in self if not len(res)]
        for empty_residue in empty_residues:
            self.remove(empty_residue)

    def purgeResidues(self, residues, mode):
        """
        Removes or keeps the selected residues.

        Parameters
        ----------
        atoms : [ProtoCaller.IO.PDB.Atom.Atom]
            Residues to be removed or kept
        mode : str
            One of "keep" and "discard" - either keeps only the selection or discards the latter.
        """
        assert mode in ["keep", "discard"]
        residues_to_remove = residues if mode == "discard" else [x for x in self if x not in residues]
        for residue in residues_to_remove:
            self.remove(residue)

    def _checkResidue(self, residue):
        """Checks whether the residue has the same chainID as the current object."""
        try:
            for prop in self._common_properties:
                if getattr(self, prop) is not None:
                    assert (getattr(self, prop) == getattr(residue, prop))
        except AttributeError:
            raise TypeError("Need to pass a valid object with a chainID attribute")
        except AssertionError:
            raise ValueError("Input residue(s) from an incompatible chain")
