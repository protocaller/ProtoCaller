import copy as _copy

import ProtoCaller as _PC

from . import _Helper_Mixin


class MissingResidue(_Helper_Mixin.HelperMixin):
    def __init__(self, resName, chainID, resSeq, iCode=" "):
        self.resName = resName
        self.chainID = chainID
        self.resSeq = resSeq
        self.iCode = iCode

    def __copy__(self):
        cls = self.__class__
        obj = cls.__new__(cls)
        obj.__dict__.update(self.__dict__)
        return obj

    def __deepcopy__(self, memo):
        cls = self.__class__
        obj = cls.__new__(cls)
        obj.__dict__.update(self.__dict__)
        memo[id(self)] = obj
        for k, v in self.__dict__.items():
            setattr(obj, k, _copy.deepcopy(v, memo))
        return obj

    def __lt__(self, other):
        if not isinstance(other, MissingResidue):
            raise TypeError("Argument needs to be of type Residue")
        else:
            str1 = "{:1.1}{:5d}{:1.1}".format(self.chainID, self.resSeq, self.iCode)
            str2 = "{:1.1}{:5d}{:1.1}".format(other.chainID, other.resSeq, other.iCode)
            return True if str1 < str2 else False

    def __gt__(self, other):
        return not self.__lt__(other)

    def __eq__(self, other):
        return self.sameResidue(other)

    def __getattr__(self, item):
        return getattr(self, "_" + item)

    def __hash__(self):
        return id(self)

    def __setattr__(self, key, value):
        value = value.strip()

        # checks
        if key == "chainID":
            if not isinstance(value, str) or len(value) != 1:
                raise ValueError("ChainID must be a single character")
            value = value.upper()
        elif key == "resSeq":
            value = int(float(value))
        elif key == "iCode":
            if not isinstance(value, str) or len(value) > 1:
                raise ValueError("iCode must be a single character")
            value = value if value != "" else " "
        elif key != "resName":
            raise AttributeError("No attribute {}".format(key))

        # setter
        super(MissingResidue, self).__setattr__(self, "_" + key, value)

    def __str__(self):
        string = "REMARK 465     {:3.3} {:1.1} {:5d}{:1.1}\n".format(self.resName, self.chainID, self.resSeq,
                                                                     self.iCode)
        return string

    @property
    def type(self):
        return _PC.RESIDUETYPE(self.resName)


class MissingAtoms(MissingResidue):
    def __init__(self, resName, chainID, resSeq, iCode=" ", atoms=None):
        super().__init__(resName=resName, chainID=chainID, resSeq=resSeq, iCode=iCode)
        self._atomlist = atoms

    def __str__(self):
        string = "REMARK 470     {:3.3} {:1.1}{:4d}{:1.1}    ".format(self.resName, self.chainID, self.resSeq,
                                                                      self.iCode)
        for atom in self._atomlist:
            string += "{:<5.5}".format(atom)
        return string + "\n"
