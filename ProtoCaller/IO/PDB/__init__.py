import os as _os
import re as _re

import pypdb as _pypdb

import ProtoCaller.Utils.ConditionalList as _CondList
from . import _Helper_Mixin
from . import Atom as _Atom
from . import Missing as _Missing
from . import Residue as _Residue
from . import Chain as _Chain


class PDB(_Chain.Chain, _CondList.ConditionalList):
    # properties which have to be conserved within the whole PDB
    _common_properties = []

    def __init__(self, input):
        _CondList.ConditionalList.__init__(self, [])
        self.readPDB(input)

    @property
    def type(self):
        return "protein"

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = _os.path.abspath(value)
        if not _os.path.exists(self._filename):
            try:
                self._filename = self.downloadPDB(value)
            except:
                raise ValueError("Input '{}' is not a valid filename or PDB code".format(value))

    def readPDB(self, filename):
        self.filename = filename
        pdb_string = open(self.filename).readlines()
        if pdb_string[-1][:3] != "END": pdb_string += ["END" + " " * 80]

        curr_atom, curr_res, curr_chain = None, _Residue.Residue(), _Chain.Chain()

        # first populate with all atoms
        for i, line in enumerate(pdb_string):
            is_atom = True if line[:4] == "ATOM" or line[:6] == "HETATM" else False
            is_end = True if line[:3] in ["TER", "END"] or i == len(pdb_string) - 1 else False

            if any([is_atom, is_end]):
                if is_atom:
                    curr_atom = _Atom.Atom(line)
                if len(curr_res) and any([is_end, not curr_atom.sameResidue(curr_res)]):
                    curr_chain.append(curr_res)
                    curr_res = _Residue.Residue()
                if len(curr_chain) and any([is_end, not curr_atom.sameChain(curr_chain)]):
                    self.append(curr_chain)
                    curr_chain = _Chain.Chain()
                if is_atom:
                    curr_res.append(curr_atom)

        self._disulfide_bonds = []
        self._site_residues = []
        self._modified_residues = []
        self._missing_residues = []
        self._missing_atoms = []

        # then re-read and create links to other parts of the PDB
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
                    chainID, resSeq, iCode = line[i], int(float(line[i + 1:i + 5])), line[i + 5]
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
                file.write("MODRES {:>4d} {:>3.3} {:1.1} {:>4d}{:1.1}\n".format(i + 1, residue.resName, residue.chainID,
                                                                                residue.resSeq, residue.iCode))
            for i, pair in enumerate(self._disulfide_bonds):
                template = "SSBOND{:>4d} CYS {:1.1} {:>4d}{:1.1}   CYS {:>1.1} {:>4d}{:1.1}\n"
                file.write(template.format(i + 1, pair[0].chainID, pair[0].resSeq, pair[0].iCode,
                                           pair[1].chainID, pair[1].resSeq, pair[1].iCode))
            for i in range(0, len(self._site_residues), 4):
                str_list = ["SITE              "]
                for residue in self._site_residues[i:min(i + 4, len(self._site_residues))]:
                    str_list += "{:>3.3} {:1.1}{:>4d}{:1.1} ".format(residue.resName, residue.chainID, residue.resSeq,
                                                                     residue.iCode)
                file.write("".join(str_list) + "\n")
            for chain in self:
                file.write(chain.__str__())
                if chain.type == "chain":
                    file.write("TER\n")
            file.write("END")
        return filename

    def writeHetatms(self, filebase=None):
        if filebase is None: filebase = _os.path.splitext(self.filename)[0]
        filenames = []
        moltypes = []
        for chain in self.filter("type!='chain'", type="chains"):
            for residue in chain:
                filename = _os.path.abspath("%s_%s.pdb" % (filebase, residue.type))
                mode = "w" if filename not in filenames else "a"
                with open(filename, mode) as file:
                    file.write(residue.__str__())
                if filename not in filenames:
                    filenames += [filename]
                    moltypes += [residue.type]
        return filenames, moltypes

    def filter(self, mask, type="residues"):
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
            for chain in self:
                if type != "chains":
                    for residue in chain:
                        if type != "residues":
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

    def numberOfResidues(self):
        return sum([chain.numberOfResidues() for chain in self])

    def numberOfChains(self):
        return len(self)

    def reNumberResidues(self, start=1, custom_resSeqs=None, custom_iCodes=None):
        if custom_resSeqs is None:
            custom_resSeqs = [start + i for i in range(self.numberOfResidues())]
        if custom_iCodes is None:
            custom_iCodes = [" "] * len(custom_resSeqs)
        if not len(custom_resSeqs) == len(custom_iCodes) == self.numberOfResidues():
            raise ValueError("Custom number of residues does not match chain number of residues")

        i = 0
        for chain in self:
            n = chain.numberOfResidues()
            chain.reNumberResidues(custom_resSeqs=custom_resSeqs[i:i + n], custom_iCodes=custom_iCodes[i:i + n])
            i += n

    def purgeAtoms(self, atoms):
        raise NotImplementedError("Inherited function not implemented for object of type 'PDB'")

    def purgeEmpty(self):
        for chain in self:
            chain.purgeEmpty()

        _Chain.Chain.purgeEmpty(self)

    def purgeResidues(self, residues):
        for name in ["_missing_residues", "_modified_residues", "_site_residues"]:
            setattr(self, name, [item for item in getattr(self, name) if item not in residues])

        for i in reversed(range(len(self._disulfide_bonds))):
            if self._disulfide_bonds[i][0] not in residues or self._disulfide_bonds[i][1] not in residues:
                del self._disulfide_bonds[i]

        for chain in self:
            chain.purgeResidues(residues)
        self.purgeEmpty()

    def totalResidueList(self, sort=True):
        total_res = self._missing_residues + self.filter("type=='amino_acid'")
        if sort:
            PDB.sortResidueList(total_res)
        return total_res

    @staticmethod
    def sortResidueList(residuelist):
        if not isinstance(residuelist, list):
            raise TypeError("Need to input a list of Residues")

        sortingfunc = lambda res: "{}{:5d}{}".format(res.chainID, res.resSeq, res.iCode)
        residuelist.sort(key=sortingfunc)

    @staticmethod
    def downloadPDB(code, filename=None):
        code = code.strip().upper()
        try:
            pdb_string = _pypdb.get_pdb_file(code, filetype="pdb", compression=False)
        except:
            raise ValueError("Invalid PDB code: '%s'" % code)
        if filename is None:
            filename = code + ".pdb"

        with open(filename, "w") as file:
            file.write(pdb_string)

        return _os.path.abspath(filename)
