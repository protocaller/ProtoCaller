import ProtoCaller.Utils.ConditionalList as _CondList
from . import Residue as _Residue


class Chain(_Residue, _CondList.ConditionalList):
    # properties which have to be conserved within the whole chain
    _common_properties = ["chainID"]

    def __init__(self, residues=None):
        if residues is None:
            residues = []
        _CondList.ConditionalList.__init__(self, residues, self._checkResidue)

    def __eq__(self, other):
        return self.sameChain(other)

    def __repr__(self):
        return "<Chain of {} residues>".format(len(self))

    @property
    def type(self):
        if not len(self):
            return None
        else:
            for residue in self:
                if residue.type in ["amino_acid"]:
                    return "chain"
            return "molecules"

    def numberOfAtoms(self):
        return sum([residue.numberOfAtoms() for residue in self])

    def numberOfResidues(self):
        return len(self)

    def reNumberAtoms(self, start=1, custom_serials=None):
        if custom_serials is None:
            custom_serials = [start + i for i in range(self.numberOfAtoms())]
        if not len(custom_serials) == self.numberOfAtoms():
            raise ValueError("Custom number of atoms does not match chain number of atoms")

        i = 0
        for residue in self:
            residue.reNumberAtoms(custom_serials=custom_serials[i:i + residue.numberOfAtoms()])
            i += residue.numberOfAtoms()

    def reNumberResidues(self, start=1, custom_resSeqs=None, custom_iCodes=None):
        if custom_resSeqs is None:
            custom_resSeqs = [start + i for i in range(len(self))]
        if custom_iCodes is None:
            custom_iCodes = [" "] * len(custom_resSeqs)
        if not len(custom_resSeqs) == len(custom_iCodes) == len(self):
            raise ValueError("Custom number of residues does not match chain number of residues")

        for resSeq, iCode, residue in zip(custom_resSeqs, custom_iCodes, self):
            residue.resSeq, residue.iCode = resSeq, iCode

    def purgeAtoms(self, atoms):
        for residue in self:
            residue.purgeAtoms(atoms)
        self.purgeEmpty()

    def purgeEmpty(self):
        empty_residues = [res for res in self if not len(res)]
        for empty_residue in empty_residues:
            self.remove(empty_residue)

    def purgeResidues(self, residues):
        for residue in residues:
            self.remove(residue)

    def _checkResidue(self, residue):
        # checks whether the residue has the same chainID as the current object
        try:
            for prop in self._common_properties:
                if getattr(self, prop) is not None:
                    assert (getattr(self, prop) == getattr(residue, prop))
        except AttributeError:
            raise TypeError("Need to pass a valid object with a chainID attribute")
        except AssertionError:
            raise ValueError("Input residue(s) from an incompatible chain")
