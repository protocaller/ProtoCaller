import const as _const

import copy as _copy
import os as _os
import re as _re
import subprocess as _subprocess
import sys as _sys
import tempfile as _tempfile

import BioSimSpace as _BSS
import parmed as _pmd
import pypdb as _pypdb

class Atom:
    def __init__(self, pdb_line):
        self._type = pdb_line[:6].strip()
        self._serial = int(float(pdb_line[6:11]))
        self._name = pdb_line[12:16].strip()
        self._altLoc = pdb_line[16]
        self._resName = pdb_line[17:20]
        self._chainID = pdb_line[21]
        self._resSeq = pdb_line[22:26]
        self._iCode = pdb_line[26]
        self._x = pdb_line[30:38]
        self._y = pdb_line[38:46]
        self._z = pdb_line[46:54]
        self._occupancy = pdb_line[54:60]
        self._tempFactor = pdb_line[60:66]
        self._element = pdb_line[76:78]
        self._charge = pdb_line[78:80]

    def __str__(self):
        string = "{:<6.6}{:5d}  {:<3.3}{:>1.1}{:>3.3} {:>1.1}{:4d}{:>1.1}" \
                 "   {:>8.8}{:>8.8}{:>8.8}{:>6.6}{:>6.6}          {:>2.2}{:>2.2}".format(self._type,
                 self._serial, self._name, self._altLoc, self._resName, self._chainID, self._resSeq, self._iCode,
                 self._x, self._y, self._z, self._occupancy, self._tempFactor, self._element, self._charge)
        return string.strip() + "\n"

    def __repr__(self):
        return self.__str__()

    @property
    def _type(self):
        return self.__type

    @_type.setter
    def _type(self, val):
        if (val not in ["ATOM", "HETATM"]):
            raise ValueError("Cannot initialise Atom type from non-atom PDB line")
        self.__type = val

    @property
    def _serial(self):
        return self.__serial

    @_serial.setter
    def _serial(self, val):
        self.__serial = self._convert_to_int(val)

    @property
    def _chainID(self):
        return self.__chainID

    @_chainID.setter
    def _chainID(self, val):
        if not isinstance(val, str) and not len(val) == 0:
            raise ValueError("ChainID must be a single character")
        else:
            self.__chainID = val

    @property
    def _resSeq(self):
        return self.__resSeq

    @_resSeq.setter
    def _resSeq(self, val):
        self.__resSeq = self._convert_to_int(val)

    @classmethod
    def _convert_to_int(self, val):
        if (isinstance(val, int)):
            return val
        else:
            return int(self._convert_to_float(val))

    @classmethod
    def _convert_to_float(self, val):
        if(isinstance(val, float)):
            return val
        else:
            try:
                return (float(val))
            except:
                raise ValueError("Not a valid value, needs to be a number or a number string")

    def sameResidue(self, obj):
        try:
            return all([self._chainID == obj._chainID, self._resSeq == obj._resSeq, self._iCode == obj._iCode])
        except:
            print("Need to pass a valid object with _chainID, _resSeq and _iCode attributes")

class MissingResidue:
    def __init__(self, resName, chainID, resSeq, iCode=" "):
        self._resName = resName
        self._chainID = chainID
        self._resSeq = Atom._convert_to_int(resSeq)
        self._iCode = iCode

    def __lt__(self, other):
        if not isinstance(other, MissingResidue):
            raise TypeError("Argument needs to be of type Residue")
        else:
            str1 = "{:1.1}{:5d}{:1.1}".format(self._chainID, self._resSeq, self._iCode)
            str2 = "{:1.1}{:5d}{:1.1}".format(other._chainID, other._resSeq, other._iCode)
            if str1 < str2:
                return True
            else:
                return False

    def __gt__(self, other):
        return not self.__lt__(other)

    @property
    def _type(self):
        return MissingResidue.classify(self._resName)

    @staticmethod
    def classify(resName):
        resName = resName.strip().upper()
        if resName in _const.WATERALIAS:
            return "water"
        elif resName in _const.CATIONLIST:
            return "cation"
        elif resName in _const.ANIONLIST:
            return "anion"
        elif resName in _const.AMINOACIDLIST:
            return "amino acid"
        elif resName in _const.MISCLIST:
            return "misc"
        else:
            return "ligand"

class Residue(MissingResidue):
    def __init__(self):
        self.__atomList = []
        self._chainID = None
        self._resSeq = None
        self._iCode = None

    def __iter__(self):
        for atom in self.__atomList:
            yield atom

    def __getitem__(self, key):
        return self.__atomList[key]

    def __len__(self):
        return len(self.__atomList)

    def __str__(self):
        return "".join([atom.__str__() for atom in self.__atomList])

    def __repr__(self):
        return "<Residue of {} atoms>".format(len(self.__atomList))

    @property
    def _chainID(self):
        return self.__chainID

    @_chainID.setter
    def _chainID(self, value):
        for atom in self.__atomList:
            atom._chainID = value
        self.__chainID = value

    @property
    def _resSeq(self):
        return self.__resSeq

    @_resSeq.setter
    def _resSeq(self, value):
        for atom in self.__atomList:
            atom._resSeq = value
        self.__resSeq = value

    @property
    def _iCode(self):
        return self.__iCode

    @_iCode.setter
    def _iCode(self, value):
        for atom in self.__atomList:
            atom._iCode = value
        self.__iCode = value

    def addAtoms(self, *atoms):
        for atom in atoms:
            if not isinstance(atom, Atom):
                raise TypeError("Input variables must be of type Atom")

            if self.__atomList == []:
                self.__atomList.append(atom)
                self._resSeq = atom._resSeq
                self._chainID = atom._chainID
                self._resName = atom._resName
                self._iCode = atom._iCode

            elif self._resSeq == atom._resSeq and self._chainID == atom._chainID and self._iCode == atom._iCode:
                self.__atomList.append(atom)

            else:
                raise ValueError("The following atom did not have the same residue number as the rest of the atoms: %s", atom)

    def numberOfAtoms(self):
        return len(self.__atomList)

    def removeAtoms(self, *atoms):
        self.__atomList = [atom for atom in self.__atomList if atom in atoms]

    def reAssign(self, chainID=None, resName=None, iCode=None):
        for name, attribute in zip(["chainID", "resName", "iCode"], [chainID, resName, iCode]):
            if(attribute is not None):
                setattr(self, name, attribute)

    def reNumberAtoms(self, start=1, custom_serials=None):
        if custom_serials is None:
            serials = [start + i for i in range(start, start + len(self.__atomList))]
        elif len(custom_serials) == len(self.__atomList):
            serials = custom_serials
        else:
            raise ValueError("Custom number list does not match number of atoms")

        for i, atom in zip(serials, self.__atomList):
            atom._serial = i

class Chain:
    def __init__(self):
        self.__residueList = []
        self._chainID = None

    def __iter__(self):
        for residue in self.__residueList:
            yield residue

    def __getitem__(self, key):
        return self.__residueList[key]

    def __len__(self):
        return len(self.__residueList)

    def __str__(self):
        return "".join([residue.__str__() for residue in self.__residueList])

    def __repr__(self):
        return "<Chain of {} residues>".format(len(self.__residueList))

    @property
    def _chainID(self):
        return self.__chainID

    @_chainID.setter
    def _chainID(self, value):
        for residue in self.__residueList:
            residue._chainID = value
        self.__chainID = value

    @property
    def _type(self):
        return self.classify()

    def classify(self):
        if not len(self.__residueList):
            return None
        else:
            for residue in self.__residueList:
                if residue._type in ["amino acid"]:
                    return "chain"
            return "molecules"

    def addResidues(self, *residues):
        for residue in residues:
            if not isinstance(residue, Residue):
                raise TypeError("Input variables must be of type Residue")

            if self.__residueList == []:
                self.__residueList.append(residue)
                self._chainID = residue._chainID

            elif self._chainID == residue._chainID:
                self.__residueList.append(residue)

            else:
                raise ValueError("The following residue did not have the same chain ID as the rest of the residues: %s", residue)

    def insertResidue(self, position, residue):
        if not isinstance(residue, Residue):
            raise TypeError("Element needs to be residue")
        self.__residueList.insert(position, residue)

    def numberOfAtoms(self):
        n = 0
        for residue in self.__residueList:
            n += residue.numberOfAtoms()
        return n

    def removeResidues(self, *residues):
        self.__residueList = [residue for residue in self.__residueList if residue in residues]

    def reNumberResidues(self, start=1, custom_resSeqs=None, custom_iCodes=None):
        if custom_resSeqs is None:
            for i, residue in enumerate(self.__residueList):
                residue.reAssign(seqRes = i + start, iCode = " ")
        else:
            if len(custom_resSeqs) != len(self.__residueList) and (custom_iCodes is not None and
                                                                   len(custom_iCodes) != len(custom_resSeqs)):
                raise ValueError("Length of parameters and chain do not match")
            else:
                if custom_iCodes is None:
                    for resSeq, residue in zip(custom_resSeqs, self.__residueList):
                        residue.reAssign(resSeq=resSeq)
                else:
                    for resSeq, iCode, residue in zip(custom_resSeqs, custom_iCodes, self.__residueList):
                        residue.reAssign(resSeq=resSeq, iCode=iCode)

    def reNumberAtoms(self, start=1, custom_serials=None):
        if custom_serials is None:
            serials = [start + i for i in range(self.numberOfAtoms())]
        elif len(custom_serials) == self.numberOfAtoms():
            serials = custom_serials
        else:
            raise ValueError("Custom number list does not match number of atoms")

        i = 0
        for residue in self.__residueList:
            residue.reNumberAtoms(custom_serials=serials[i:i+residue.numberOfAtoms()])
            i += residue.numberOfAtoms()


    def reLabel(self, chainID="A"):
        for residue in self.__residueList:
            residue.reAssign(chainID=chainID)

    @property
    def _chainID(self):
        return self.__chainID

    @_chainID.setter
    def _chainID(self, val):
        for residue in self.__residueList:
            residue._chainID = val
        if len(self.__residueList):
            self.__chainID = self.__residueList[0]._chainID
        else:
            self.__chainID = None

class PDB:
    def __init__(self, filename):
        self.readPDB(filename)

    def __getitem__(self, key):
        return self.__chainList[key]

    def __iter__(self):
        for chain in self.__chainList:
            yield chain

    def readPDB(self, filename):
        self.filename = filename
        self.__chainList = []
        pdb_string = open(filename).readlines()

        curr_atom, curr_res, curr_chain = None, Residue(), Chain()

        #first populate with all atoms
        for i, line in enumerate(pdb_string):
            if line[:4] == "ATOM" or line[:6] == "HETATM" or line[:3] == "TER" or i == len(pdb_string) - 1:
                if line[:4] == "ATOM" or line[:6] == "HETATM":
                    curr_atom = Atom(line)
                if len(curr_res) and (line[:3] == "TER" or i == len(pdb_string) - 1 or
                                      not curr_atom.sameResidue(curr_res)):
                    curr_chain.addResidues(curr_res)
                    curr_res = Residue()
                if len(curr_chain) and (line[:3] == "TER" or i == len(pdb_string) - 1 or
                                        curr_atom._chainID != curr_chain._chainID):
                    self.addChains(curr_chain)
                    curr_chain = Chain()
                if line[:4] == "ATOM" or line[:6] == "HETATM":
                    curr_res.addAtoms(curr_atom)

        self._disulfide_bonds = []
        self._site_residues = []
        self._modified_residues = []
        self._missing_residues = []

        #then re-read and create links to other parts of the PDB
        for line in pdb_string:
            if line[:6] == "SSBOND":
                chainIDs = [line[15], line[29]]
                resSeqs = [int(float(line[17:21])), int(float(line[31:35]))]
                iCodes = [line[21], line[35]]
                self._disulfide_bonds += [self.filterResidues(includedict={"_chainID" : chainIDs, "_resSeq" : resSeqs,
                                                                           "_iCode" : iCodes})]
            elif line[:4] == "SITE":
                chainIDs, resSeqs, iCodes = [], [], []
                for i in range(22, len(line.strip()), 11):
                    chainIDs += [line[i]]
                    resSeqs += [int(float(line[i+1:i+5]))]
                    iCodes += [line[i+5]]
                self._site_residues += self.filterResidues(includedict={"_chainID" : chainIDs, "_resSeq" : resSeqs,
                                                                        "_iCode": iCodes})
            elif line[:6] == "MODRES":
                self._modified_residues += self.filterResidues(includedict={"_chainID" : [line[16]],
                                                                            "_resSeq" : [int(float(line[18:22]))],
                                                                            "_iCode" : [line[22]]})
            else:
                match = _re.findall(r"^REMARK 465\s*([\w]{3})\s*([\w])\s*([\d]*)\s*$", line)
                if len(match) != 0:
                    self._missing_residues += [MissingResidue(*match[0])]

    def writePDB(self, filename=None):
        if filename is None: filename = self.filename
        with open(filename, "w") as file:
            for i, residue in enumerate(self._modified_residues):
                file.write("MODRES {:>4d} {:>3.3} {:1.1} {:>4d}{:1.1}\n".format(i + 1, residue._resName, residue._chainID,
                                                                                residue._resSeq, residue._iCode))
            for i, pair in enumerate(self._disulfide_bonds):
                file.write("SSBOND{:>4d} CYS {:1.1} {:>4d}{:1.1}   CYS {:>1.1} {:>4d}{:1.1}\n".format(i + 1, pair[0]._chainID,
                                    pair[0]._resSeq, pair[0]._iCode, pair[1]._chainID, pair[1]._resSeq, pair[1]._iCode))
            for i in range(0, len(self._site_residues), 4):
                str_list = ["SITE              "]
                for residue in self._site_residues[i:min(i + 4, len(self._site_residues))]:
                    str_list += "{:>3.3} {:1.1}{:>4d}{:1.1} ".format(residue._resName, residue._chainID, residue._resSeq,
                                                                    residue._iCode)
                file.write("".join(str_list) + "\n")
            for chain in self.__chainList:
                file.write(chain.__str__())
                if chain._type == "chain":
                    file.write("TER\n")
        return filename

    def writeHetatms(self, filebase=None):
        if filebase == None: filebase = self.filename.split(".")[0]
        filenames = []
        moltypes = []
        for chain in self.__chainList:
            if chain._type == "chain":
                continue
            for residue in chain:
                filename = "%s_%s.pdb" % (filebase, residue._type)
                with open(filename, "a") as file:
                    file.write(residue.__str__())
                if filename not in filenames:
                    filenames += [filename]
                    moltypes += [residue._type]
        return filenames, moltypes

    #dictionaries whose keys are Residue variables and whose values are lists of permitted values
    def filterResidues(self, includedict=None, excludedict=None, includejoint=True, excludejoint=False,
                       includepolicy="AND", excludepolicy="OR"):
        if includedict == None: includedict = {}
        if excludedict == None: excludedict = {}
        if not isinstance(includedict, dict) or not isinstance(excludedict, dict):
            raise TypeError("Inputs must be dictionaries")
        for key1, values1 in includedict.items():
            for key2, values2 in excludedict.items():
                if not isinstance(values1, list) or not isinstance(values2, list):
                    raise TypeError("Both dictionaries should have lists as values")
                if key1 == key2:
                    raise ValueError("Cannot have the same keys in both dictionaries")
        if includejoint == True:
            length = None
            for key, values in includedict.items():
                if length == None:
                    length = len(values)
                elif len(values) != length:
                    raise ValueError("Cannot have joint criteria with different sizes")
        if excludejoint == True:
            length = None
            for key, values in excludedict.items():
                if length == None:
                    length = len(values)
                elif len(values) != length:
                    raise ValueError("Cannot have joint criteria with different sizes")

        reslist = []
        for i, chain in enumerate(self.__chainList):
            for j, residue in enumerate(self.__chainList[i]):
                include = True
                if not includejoint:
                    for key, values in includedict.items():
                        if getattr(residue, key) in values and includepolicy == "OR":
                            include = True
                            break
                        elif getattr(residue, key) not in values and includepolicy == "AND":
                            include = False
                            break
                else:
                    keys = list(includedict.keys())
                    for argset in zip(*includedict.values()):
                        include = True
                        for key, arg in zip(keys, argset):
                            if getattr(residue, key) != arg:
                                include = False
                                break
                        if include:
                            break
                if not excludejoint:
                    for key, values in excludedict.items():
                        if getattr(residue, key) in values and excludepolicy == "OR":
                            include = False
                            break
                        elif getattr(residue, key) not in values and excludepolicy == "AND":
                            include = True
                            break
                else:
                    keys = list(excludedict.keys())
                    for argset in zip(*excludedict.values()):
                        include = True
                        for key, arg in zip(keys, argset):
                            if getattr(residue, key) != arg:
                                include = False
                                break
                        if not include:
                            break
                if include:
                    reslist += [self.__chainList[i][j]]

        return reslist

    def purgeResidues(self, residues):
        for i, reslist in enumerate([self._missing_residues, self._modified_residues, self._site_residues]):
            for j in reversed(range(len(reslist))):
                if reslist[j] not in residues:
                    del reslist[j]

        for i in reversed(range(len(self._disulfide_bonds))):
            if self._disulfide_bonds[i][0] not in residues or self._disulfide_bonds[i][1] not in residues:
                del self._disulfide_bonds[i]

        for i in reversed(range(len(self.__chainList))):
            self.__chainList[i].removeResidues(*residues)
            if not len(self.__chainList[i]):
                del self.__chainList[i]

    def addChains(self, *chains):
        for chain in chains:
            if not isinstance(chain, Chain):
                raise TypeError("Input variables must be of type Chain")
            else:
                self.__chainList += [chain]

    def numberOfAtoms(self):
        n = 0
        for chain in self.__chainList:
            n += chain.numberOfAtoms()
        return n

    def reNumberResidues(self, start=1, custom_resSeqs=None, custom_iCodes=None):
        if custom_resSeqs is None:
            index = 1
            for chain in self.__chainList:
                chain.reNumberResidues(start=index)
                index += len(chain)
        else:
            length = 0
            for chain in self.__chainList:
                length += len(chain)

            if len(custom_resSeqs) != length and (custom_iCodes is not None and len(custom_iCodes) != len(custom_resSeqs)):
                raise ValueError("Length of parameters and chain do not match")
            else:
                if custom_iCodes is None:
                    index = 0
                    for chain in self.__chainList:
                        chain.reAssign(resSeq=custom_resSeqs[index:index+len(chain)])
                else:
                    index = 0
                    for chain in self.__chainList:
                        chain.reAssign(resSeq=custom_resSeqs[index:index + len(chain)],
                                       iCode=custom_iCodes[index:index + len(chain)])

    def reNumberAtoms(self, start=1, custom_serials=None):
        if custom_serials is None:
            serials = [start + i for i in range(self.numberOfAtoms())]
        elif len(custom_serials) == self.numberOfAtoms():
            serials = custom_serials
        else:
            raise ValueError("Custom number list does not match number of atoms")

        i = 0
        for chain in self.__chainList:
            chain.reNumberAtoms(custom_serials=serials[i:i+chain.numberOfAtoms()])
            i += chain.numberOfAtoms()


    def totalResidueList(self, sorted=True):
        total_res = self._missing_residues + self.filterResidues(includedict={"_type" : ["amino acid"]})
        if(sorted):
            PDB.sortResidueList(total_res)
        return total_res


    @staticmethod
    def sortResidueList(residuelist):
        if not isinstance(residuelist, list):
            raise TypeError("Need to input a list of Residues")

        sortingfunc = lambda res: "{}{:5d}{}".format(res._chainID, res._resSeq, res._iCode)
        residuelist.sort(key=sortingfunc)

def downloadPDB(id, filename=None):
    if type(id) is not str:
        raise TypeError("'id' must be of type 'str' with length 4")
    else:
        id = id.replace(" ", "").upper()
        if len(id) != 4:
            raise TypeError("'id' must be of type 'str' with length 4")

    if(filename == None):
        filename = id + ".pdb"

    try:
        pdb_string = _pypdb.get_pdb_file(id, filetype="pdb", compression=False)
    except:
        raise ValueError("Invalid PDB ID: '%s'" % id)

    with open("%s.pdb" % id, "w") as file:
        file.write(pdb_string)

    return filename