import re as _re

import pypdb as _pypdb

from . import Atom as _Atom
from . import Missing as _Missing
from . import Residue as _Residue
from . import Chain as _Chain

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
        if pdb_string[-1][:3] != "END": pdb_string += ["END" + " "*80]

        curr_atom, curr_res, curr_chain = None, _Residue.Residue(), _Chain.Chain()

        #first populate with all atoms
        for i, line in enumerate(pdb_string):
            if line[:4] == "ATOM" or line[:6] == "HETATM" or line[:3] in ["TER", "END"] or i == len(pdb_string) - 1:
                if line[:4] == "ATOM" or line[:6] == "HETATM":
                    curr_atom = _Atom.Atom(line)
                if len(curr_res) and (line[:3] in ["TER", "END"] or i == len(pdb_string) - 1 or
                                      not curr_atom.sameResidue(curr_res)):
                    curr_chain.append(curr_res)
                    curr_res = _Residue.Residue()
                if len(curr_chain) and (line[:3] in ["TER", "END"] or i == len(pdb_string) - 1 or
                                        curr_atom.chainID != curr_chain.chainID):
                    self.addChains(curr_chain)
                    curr_chain = _Chain.Chain()
                if line[:4] == "ATOM" or line[:6] == "HETATM":
                    curr_res.addAtoms(curr_atom)

        self._disulfide_bonds = []
        self._site_residues = []
        self._modified_residues = []
        self._missing_residues = []
        self._missing_atoms = []

        #then re-read and create links to other parts of the PDB
        for line in pdb_string:
            if line[:6] == "SSBOND":
                chainIDs = [line[15], line[29]]
                resSeqs = [int(float(line[17:21])), int(float(line[31:35]))]
                iCodes = [line[21], line[35]]
                masks = []
                for chainID, resSeq, iCode in zip(chainIDs, resSeqs, iCodes):
                    masks += ["chainID=='{}'&resSeq=={}&iCode=='{}'".format(chainID, resSeq, iCode)]
                mask = "|".join(masks)
                self._disulfide_bonds += [self.filter(mask)]
            elif line[:4] == "SITE":
                for i in range(22, len(line.strip()), 11):
                    chainID, resSeq, iCode = line[i], int(float(line[i+1:i+5])), line[i+5]
                    mask = "chainID=='{}'&resSeq=={}&iCode=='{}'".format(chainID, resSeq, iCode)
                    self._site_residues += self.filter(mask)
            elif line[:6] == "MODRES":
                chainID, resSeq, iCode = line[16], int(float(line[18:22])), line[22]
                mask = "chainID=='{}'&resSeq=={}&iCode=='{}'".format(chainID, resSeq, iCode)
                self._modified_residues += self.filter(mask)
            else:
                match = _re.findall(r"^REMARK 465\s*([\w]{3})\s*([\w])\s*([\d]+)([\D|\S]?)\s*$", line)
                if len(match) != 0:
                    self._missing_residues += [_Missing.MissingResidue(*match[0])]
                    continue

                match = _re.findall(r"^REMARK 470\s*([\w]{3})\s*([\w])\s*([\d]+)([\D|\S]?)\s*", line)
                if len(match) != 0:
                    args = match[0]
                    arr = line[21:].split()
                    self._missing_atoms += [_Missing.MissingAtoms(*args, atoms=arr)]

    def writePDB(self, filename=None):
        if filename is None: filename = self.filename
        with open(filename, "w") as file:
            for residue in self._missing_residues + self._missing_atoms:
                file.write(residue.__str__())
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
            file.write("END")
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
                mode = "w" if filename not in filenames else "a"
                with open(filename, mode) as file:
                    file.write(residue.__str__())
                if filename not in filenames:
                    filenames += [filename]
                    moltypes += [residue._type]
        return filenames, moltypes

    def filter(self, mask, filter="residues"):
        def add(elem):
            nonlocal i, condition, all_elems, curr_list
            if i == 0:
                all_elems += [elem]
            if eval("elem" + "." + condition, globals(), vars()):
                curr_list += [elem]
        all_elems = []
        all_sets = []
        expr = r"[^\(\)&|]+"
        conditions = [x for x in _re.findall(expr, mask) if x not in [None, ""]]

        for i, condition in enumerate(conditions):
            curr_list = []
            for chain in self.__chainList:
                if filter != "chains":
                    for residue in chain:
                        if filter != "residues":
                            for atom in residue:
                                add(atom)
                        else:
                            add(residue)
                else:
                    add(chain)

            all_sets += [frozenset(curr_list)]
            mask = mask.replace(condition, "all_sets[%d]" % i)

        total_set = eval(mask, globals(), locals())
        return [x for x in all_elems if x in total_set]

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
            if not isinstance(chain, _Chain.Chain):
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
        total_res = self._missing_residues + self.filter("_type=='amino_acid'")
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
