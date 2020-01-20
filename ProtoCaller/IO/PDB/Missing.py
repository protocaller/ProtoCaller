import ProtoCaller as _PC

from . import _Helper_Mixin


class MissingResidue(_Helper_Mixin.HelperMixin):
    """
    Represents a missing residue in a PDB file.

    Parameters
    ----------
    resName : str
        Initialises resName.
    chainID : str
        Initialises chainID.
    resSeq : str, int
        Initialises resSeq.
    iCode : str
        Initialises iCode.

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
    def __init__(self, resName, chainID, resSeq, iCode=" "):
        self.resName = resName
        self.chainID = chainID
        self.resSeq = resSeq
        self.iCode = iCode

    def __lt__(self, other):
        if not isinstance(other, MissingResidue):
            raise TypeError("Argument needs to be of type Residue")
        else:
            try:
                str1 = "{:1.1}{:5d}{:1.1}".format(self.chainID, self.resSeq, self.iCode)
                str2 = "{:1.1}{:5d}{:1.1}".format(other.chainID, other.resSeq, other.iCode)
                return True if str1 < str2 else False
            except TypeError:
                raise TypeError("Cannot compare an empty residue")

    def __gt__(self, other):
        return not self.__lt__(other)

    def __eq__(self, other):
        return self.sameResidue(other)

    def __getattr__(self, item):
        return getattr(self, "_" + item)

    def __hash__(self):
        return id(self)

    def __setattr__(self, key, value):
        if isinstance(value, str): value = value.strip()

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
        super(MissingResidue, self).__setattr__("_" + key, value)

    def __str__(self):
        string = "REMARK 465     {:3.3} {:1.1} {:5d}{:1.1}\n".format(
            self.resName, self.chainID, self.resSeq, self.iCode)
        return string

    @property
    def sequence(self):
        """str: Returns the single-character code of the amino acid residue. For
                missing residues this is always '-'."""
        return "-"

    @property
    def type(self):
        """str: Returns the type of the residue, determined by
               ProtoCaller.RESIDUETYPE."""
        return _PC.RESIDUETYPE(self.resName)


class MissingAtoms(MissingResidue, list):
    """
    Represents all missing atoms corresponding to a single residue in a PDB file.

    Parameters
    ----------
    resName : str
        Initialises resName.
    chainID : str
        Initialises chainID.
    resSeq : str, int
        Initialises resSeq.
    iCode : str
        Initialises iCode.
    atoms : [ProtoCaller.IO.PDB.Atom.Atom]
        Initialises the missing atoms.
    """
    def __init__(self, resName, chainID, resSeq, iCode=" ", atoms=None):
        MissingResidue.__init__(self, resName=resName, chainID=chainID, resSeq=resSeq, iCode=iCode)
        list.__init__(self, atoms)

    def __str__(self):
        string = "REMARK 470     {:3.3} {:1.1}{:4d}{:1.1}    ".format(
            self.resName, self.chainID, self.resSeq, self.iCode)
        for atom in self:
            string += "{:<5.5}".format(atom)
        return string + "\n"
