class Atom:
    _common_properties = ["type", "serial", "name", "altloc", "resName", "chainID", "resSeq", "iCode", "x", "y", "z",
                          "occupancy", "tempFactor", "element", "charge"]

    def __init__(self, pdb_line):
        self.type = pdb_line[:6]
        self.serial = pdb_line[6:11]
        self.name = pdb_line[12:16]
        self.altloc = pdb_line[16]
        self.resName = pdb_line[17:20]
        self.chainID = pdb_line[21]
        self.resSeq = pdb_line[22:26]
        self.iCode = pdb_line[26]
        self.x = pdb_line[30:38]
        self.y = pdb_line[38:46]
        self.z = pdb_line[46:54]
        self.occupancy = pdb_line[54:60]
        self.tempFactor = pdb_line[60:66]
        self.element = pdb_line[76:78]
        self.charge = pdb_line[78:80]

    def __getattr__(self, item):
        return getattr(self, "_" + item)

    def __str__(self):
        string = "{:<6.6}{:5d} {:<4.4}{:>1.1}{:>3.3} {:>1.1}{:4d}{:>1.1}   {:>8.8}{:>8.8}{:>8.8}{:>6.6}{:>6.6}" \
                 "          {:>2.2}{:>2.2}".format(self.type, self.serial, self.name, self.altloc, self.resName,
                                                   self.chainID, self.resSeq, self.iCode, self.x, self.y, self.z,
                                                   self.occupancy, self.tempFactor, self.element, self.charge)
        return string.strip() + "\n"

    def __setattr__(self, key, value):
        if key not in self._common_properties:
            raise ValueError("Invalid attribute {}. Attributes need to be one of {}".format(key,
                                                                                            self._common_properties))
        value = value.strip()

        # checks
        if key == "chainID":
            if not isinstance(value, str) or len(value) != 1:
                raise ValueError("ChainID must be a single character")
            value = value.upper()
        elif key in ["serial", "resSeq"]:
            value = int(float(value))
        elif key in ["x", "y", "z"]:
            value = float(value)
        elif key == "iCode":
            if not isinstance(value, str) or len(value) > 1:
                raise ValueError("iCode must be a single character")
            value = value if value != "" else " "
        elif key == "type":
            value = value.upper()
            if value not in ["ATOM", "HETATM"]:
                raise ValueError("Atom type needs to be either ATOM or HETATM")

        # setter
        super(Atom, self).__setattr__("_" + key, value)

    def __repr__(self):
        return "<Atom>"

    def sameChain(self, obj):
        from . import Chain as _Chain
        try:
            return all([getattr(self, prop) == getattr(obj, prop) for prop in _Chain.Chain._common_properties])
        except AttributeError:
            print("Need to pass a valid object with {} attributes".format(_Chain.Chain._common_properties))

    def sameResidue(self, obj):
        from . import Residue as _Residue
        try:
            return all([getattr(self, prop) == getattr(obj, prop) for prop in _Residue.Residue._common_properties])
        except AttributeError:
            print("Need to pass a valid object with {} attributes".format(_Residue.Residue._common_properties))
