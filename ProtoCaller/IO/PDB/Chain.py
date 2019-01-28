import ProtoCaller.Utils.ConditionalList as _CondList
from . import Residue as _Residue


class Chain(_Residue.Residue, _CondList.ConditionalList):
    # properties which have to be conserved within the whole chain
    _common_properties = ["chainID"]
    # methods partially inherited from Residue
    _inherited_methods = ["__copy__", "__deepcopy__", "__getattr__", "__hash__", "__setattr__", "__str__", "sameChain",
                          "sameResidue"]

    for method in _inherited_methods:
        vars()[method] = getattr(_Residue.Residue, method)

    def __init__(self, residues=None):
        if residues is None:
            residues = []
        _CondList.ConditionalList.__init__(self, residues, self._checkResidue)

    def __eq__(self, other):
        return self.sameChain(other)

    def __repr__(self):
        return "<Chain of {} residues>".format(self.__len__())

    @property
    def type(self):
        if not self.__len__():
            return None
        else:
            for residue in self.__iter__():
                if residue.type in ["amino_acid"]:
                    return "chain"
            return "molecules"

    def numberOfAtoms(self):
        return sum([residue.numberOfAtoms() for residue in self.__iter__()])

    def numberOfResidues(self):
        return self.__len__()

    def reNumberAtoms(self, start=1, custom_serials=None):
        if custom_serials is None:
            custom_serials = [start + i for i in range(self.numberOfAtoms())]
        if not len(custom_serials) == self.numberOfAtoms():
            raise ValueError("Custom number of atoms does not match chain number of atoms")

        i = 0
        for residue in self.__iter__():
            residue.reNumberAtoms(custom_serials=custom_serials[i:i + residue.numberOfAtoms()])
            i += residue.numberOfAtoms()

    def reNumberResidues(self, start=1, custom_resSeqs=None, custom_iCodes=None):
        if custom_resSeqs is None:
            custom_resSeqs = [start + i for i in range(self.__len__())]
        if custom_iCodes is None:
            custom_iCodes = [" "] * len(custom_resSeqs)
        if not len(custom_resSeqs) == len(custom_iCodes) == self.__len__():
            raise ValueError("Custom number of residues does not match chain number of residues")

        for resSeq, iCode, residue in zip(custom_resSeqs, custom_iCodes, self.__iter__()):
            residue.resSeq, residue.iCode = resSeq, iCode

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
