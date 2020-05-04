import os as _os
import re as _re

import ProtoCaller.Utils.ConditionalList as _CondList
from . import _Helper_Mixin
from .Atom import *
from .Missing import *
from .Residue import *
from .Chain import *


class PDB(Chain, _CondList.ConditionalList):
    """
    Represents a whole PDB file.

    Parameters
    ----------
    input : str
        Name of the input PDB file.

    Attributes
    ----------
    """
    # properties which have to be conserved within the whole PDB
    _common_properties = []

    def __init__(self, input):
        _CondList.ConditionalList.__init__(self, [])
        self.readPDB(input)

    def __repr__(self):
        return "<PDB of {} chains>".format(len(self))

    @property
    def type(self):
        """str: Returns "protein". A potentially useful placeholder for more
                functionality in the future."""
        return "protein"

    @property
    def filename(self):
        """str: Name of the file corresponding to the PDB object."""
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = _os.path.abspath(value)
        if not _os.path.exists(self._filename):
            raise ValueError("Input '{}' is not a valid filename or PDB code".format(value))

    def readPDB(self, filename):
        """Reads the input PDB."""
        self.filename = filename
        pdb_string = open(self.filename).readlines()
        if pdb_string[-1][:3] != "END": pdb_string += ["END" + " " * 80]

        curr_atom, curr_res, curr_chain = None, Residue(), Chain()

        # first populate with all atoms
        for i, line in enumerate(pdb_string):
            is_atom = True if line[:4] == "ATOM" or line[:6] == "HETATM" else False
            is_end = True if line[:3] in ["TER", "END"] or i == len(pdb_string) - 1 else False

            if any([is_atom, is_end]):
                if is_atom:
                    curr_atom = Atom(line)
                if len(curr_res) and any([is_end, not curr_atom.sameResidue(curr_res)]):
                    curr_chain.append(curr_res)
                    curr_res = Residue()
                if len(curr_chain) and any([is_end, not curr_atom.sameChain(curr_chain)]):
                    self.append(curr_chain)
                    curr_chain = Chain()
                if is_atom:
                    curr_res.append(curr_atom)

        self.disulfide_bonds = []
        self.site_residues = []
        self.modified_residues = []
        self.missing_residues = []
        self.missing_atoms = []

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
                self.disulfide_bonds += [self.filter(mask)]
            elif line[:4] == "SITE":
                for i in range(22, len(line.strip()), 11):
                    chainID, resSeq, iCode = line[i], int(float(line[i + 1:i + 5])), line[i + 5]
                    mask = "chainID=='{}'&resSeq=={}&iCode=='{}'".format(chainID, resSeq, iCode)
                    self.site_residues += self.filter(mask)
            elif line[:6] == "MODRES":
                chainID, resSeq, iCode = line[16], int(float(line[18:22])), line[22]
                mask = "chainID=='{}'&resSeq=={}&iCode=='{}'".format(chainID, resSeq, iCode)
                self.modified_residues += self.filter(mask)
            else:
                match = _re.findall(r"^REMARK 465\s*([\w]{3})\s*([\w])\s*([-]?[\d]+)([\D|\S]?)\s*$", line)
                if len(match) != 0:
                    self.missing_residues += [MissingResidue(*match[0])]
                    continue

                match = _re.findall(r"^REMARK 470\s*([\w]{3})\s*([\w])\s*([\d]+)([\D|\S]?)\s*", line)
                if len(match) != 0:
                    args = match[0]
                    arr = line[21:].split()
                    self.missing_atoms += [MissingAtoms(*args, atoms=arr)]

    def writePDB(self, filename=None):
        """
        Writes the object as a PDB file.

        Parameters
        ----------
        filename : str or None
            Name of the output file. Default: original filename.

        Returns
        -------
        filename : str
            The absolute path to the written PDB file.
        """
        if filename is None: filename = self.filename
        with open(filename, "w") as file:
            # write SEQRES
            total_residue_list = self.totalResidueList()
            chainIDs = [x.chainID for x in total_residue_list]
            chain_lengths = {x: chainIDs.count(x) for x in set(chainIDs)}
            line, curr_chainID, i, j = None, None, None, None
            for residue in total_residue_list:
                if residue.chainID != curr_chainID:
                    if line is not None:
                        file.write(line + "\n")
                    curr_chainID = residue.chainID
                    i, j = 1, 1
                    line = "SEQRES {:>3d} {:1.1}  {:>3d}  ".format(i, curr_chainID, chain_lengths[curr_chainID])
                if j == 14:
                    file.write(line + "\n")
                    i += 1
                    j = 1
                    line = "SEQRES {:>3d} {:1.1}  {:>3d}  ".format(i, curr_chainID, chain_lengths[curr_chainID])
                line += "{} ".format(residue.resName)
                j += 1
            if j != 1:
                file.write(line + "\n")

            # write missing atoms and residues
            for residue in self.missing_residues + self.missing_atoms:
                file.write(residue.__str__())

            # write modified residues
            for i, residue in enumerate(self.modified_residues):
                file.write("MODRES {:>4d} {:>3.3} {:1.1} {:>4d}{:1.1}\n".format(i + 1, residue.resName, residue.chainID,
                                                                                residue.resSeq, residue.iCode))

            #write disulfide bonds
            for i, pair in enumerate(self.disulfide_bonds):
                template = "SSBOND{:>4d} CYS {:1.1} {:>4d}{:1.1}   CYS {:>1.1} {:>4d}{:1.1}\n"
                file.write(template.format(i + 1, pair[0].chainID, pair[0].resSeq, pair[0].iCode,
                                           pair[1].chainID, pair[1].resSeq, pair[1].iCode))

            # write site residues
            for i in range(0, len(self.site_residues), 4):
                str_list = ["SITE              "]
                for residue in self.site_residues[i:min(i + 4, len(self.site_residues))]:
                    str_list += "{:>3.3} {:1.1}{:>4d}{:1.1} ".format(residue.resName, residue.chainID, residue.resSeq,
                                                                     residue.iCode)
                file.write("".join(str_list) + "\n")

            # write all residues
            for chain in self:
                file.write(chain.__str__())
                if chain.type == "chain":
                    file.write("TER\n")
            file.write("END")

        return _os.path.abspath(filename)

    def writeHetatms(self, filebase=None):
        """
        Partitions every type of HETATM molecules into a separate PDB file.

        Parameters
        ----------
        filebase : str or None
            Base name of the output files. Default: the original PDB name.

        Returns
        -------
        filenames : [str]
            The absolute paths of all output files.
        moltypes : [str]
            Their corresponding moltypes.
        """
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
        """
        Returns a selection based on a mask. E.g. mask="chainID='A'" returns all elements which belong to chain A.

        Parameters
        ----------
        mask : str
            A conditional expression which is to be evaluated.
        type : str
            The type of objects to be returned. One of "chains", "residues" and "atoms".

        Returns
        -------
        total_list : list
            All elements satisfying the given conditions.
        """
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
            mask = _re.sub(r"(\A|{0})({1})(\Z|{0})".format("[\(\)&|]+", _re.escape(condition)),
                           r"\1all_sets[%d]\3" % i, mask)

        total_set = eval(mask, globals(), locals())
        return [x for x in all_elems if x in total_set]

    @property
    def numberOfResidues(self):
        """int: Returns the totel number of Residues and MissingResidues."""
        return sum([chain.numberOfResidues for chain in self]) + len(self.missing_residues)

    @property
    def numberOfChains(self):
        """int: Returns the number of chains."""
        return len(self)

    def reNumberResidues(self, start=1, custom_resSeqs=None, custom_iCodes=None):
        if custom_resSeqs is None:
            custom_resSeqs = [start + i for i in range(self.numberOfResidues)]
        if custom_iCodes is None:
            custom_iCodes = [" "] * len(custom_resSeqs)
        if not len(custom_resSeqs) == len(custom_iCodes) == self.numberOfResidues:
            raise ValueError("Custom number of residues does not match chain number of residues")

        all_residues = self.totalResidueList()
        number_dict = {}
        for resSeq, iCode, residue in zip(custom_resSeqs, custom_iCodes, all_residues):
            number_dict[(residue.resSeq, residue.iCode)] = (resSeq, iCode)
            residue.resSeq, residue.iCode = resSeq, iCode
        for atom in self.missing_atoms:
            atom.resSeq, atom.iCode = number_dict[(atom.resSeq, atom.iCode)]

    def purgeAtoms(self, atoms):
        raise NotImplementedError("Inherited function not implemented for object of type 'PDB'")

    def purgeEmpty(self):
        for chain in self:
            chain.purgeEmpty()

        Chain.purgeEmpty(self)

    def purgeResidues(self, residues, mode):
        assert mode in ["keep", "discard"]
        for name in ["missing_residues", "modified_residues", "site_residues"]:
            if mode == "keep":
                setattr(self, name, [item for item in getattr(self, name) if item in residues])
            else:
                setattr(self, name, [item for item in getattr(self, name) if item not in residues])

        for i in reversed(range(len(self.disulfide_bonds))):
            if mode == "keep":
                if self.disulfide_bonds[i][0] not in residues or self.disulfide_bonds[i][1] not in residues:
                    del self.disulfide_bonds[i]
            else:
                if self.disulfide_bonds[i][0] in residues or self.disulfide_bonds[i][1] in residues:
                    del self.disulfide_bonds[i]

        for chain in self:
            chain.purgeResidues([residue for residue in residues if residue in chain], mode=mode)
        self.purgeEmpty()

    @property
    def sequence(self):
        """str: Returns the single-character sequence of all amino acids in
                the protein."""
        try:
            residues = self.totalResidueList()
            seq = {}
            for res in residues:
                if res.chainID not in seq.keys():
                    seq[res.chainID] = []
                seq[res.chainID] += [res.sequence]
            seq = "/".join(["".join(seq[k]) for k in sorted(seq.keys())])
        except:
            seq = ""
        return seq

    def totalResidueList(self, sort=True):
        """
        Returns only missing and non-missing residues with no HETATM molecules.

        Parameters
        ----------
        sort : bool
            Whether to sort the residues by number.

        Returns
        -------
        total_res : [ProtoCaller.IO.PDB.Residue.Residue]
            All residues in the PDB file.
        """
        total_res = self.missing_residues + self.filter("type in ['amino_acid', 'amino_acid_modified']")
        if sort:
            PDB.sortResidueList(total_res)
        return total_res

    @staticmethod
    def sortResidueList(residuelist):
        """A helper method which sorts the residues in the list by number."""
        if not isinstance(residuelist, list):
            raise TypeError("Need to input a list of Residues")

        sortingfunc = lambda res: (res.chainID, res.resSeq, res.iCode)
        residuelist.sort(key=sortingfunc)
