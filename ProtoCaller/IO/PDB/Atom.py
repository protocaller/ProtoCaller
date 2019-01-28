from . import _Helper_Mixin


class Atom(_Helper_Mixin.HelperMixin):
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
        if item not in self._common_properties:
            raise ValueError("Invalid attribute {}. Attributes need to be one of {}".format(item,
                                                                                            self._common_properties))
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
        if key in ["serial", "resSeq"]:
            value = int(float(value))
        elif key in ["x", "y", "z"]:
            value = float(value)
        elif key in ["iCode", "chainID"]:
            if not isinstance(value, str) or len(value) > 1:
                raise ValueError("{} must be a single character".format(key))
            value = value.upper() if value != "" else " "
        elif key == "type":
            value = value.upper()
            if value not in ["ATOM", "HETATM"]:
                raise ValueError("Atom type needs to be either ATOM or HETATM")

        # setter
        super(Atom, self).__setattr__("_" + key, value)

    def __repr__(self):
        return "<Atom of type {}>".format(self.type)
