# TODO:
# 1. See how modified protein residues are affected
# 2. Renumber residues - sometimes PDB files provide negative numbers - might crash other pieces of software

import copy as _copy
import os as _os
import re as _re

import modeller as _modeller
import modeller.automodel as _modellerautomodel

import ProtoCaller.IO as _IO


# taken from https://salilab.org/modeller/manual/node36.html
class MyLoop(_modellerautomodel.loopmodel):
    def __init__(self, filename=None, *args, **kwargs):
        super(MyLoop, self).__init__(*args, **kwargs)
        if filename is not None:
            self.pdb = _IO.PDB.PDB(filename)
            self.total_residue_list = self.pdb.totalResidueList()
            self.transformed_total_residue_list = self._transform_id(self.total_residue_list)
            self.transformed_missing_residue_list = [self.transformed_total_residue_list[i]
                                                     for i, res in enumerate(self.total_residue_list)
                                                     if type(res) == _IO.PDB.MissingResidue]
        else:
            raise ValueError("Need to specify input PDB filename")

        self.selection = _modeller.selection()

    def select_loop_atoms(self):
        for residue in self.transformed_missing_residue_list:
            try:
                self.selection.add(self.residues[residue])
            except KeyError:
                try:
                    self.selection.add(self.residues[residue[:-2]])
                except KeyError:
                    print("Error while selecting missing residue: %s. "
                          "Residue either missing from FASTA sequence or there is a problem with Modeller." % residue)
        return self.selection

    @staticmethod
    def _transform_id(reslist):
        modified_id_list = []
        for i, res in enumerate(reslist):
            modified_id_list += ["{}:{}".format(i + 1, res.chainID)]
        return modified_id_list


def modellerTransform(filename_pdb, filename_fasta, add_missing_atoms):
    env = _modeller.environ()
    pdb_code = _os.path.splitext(_os.path.basename(filename_pdb))[0]

    # convert FASTA to PIR
    filename_pir = FASTA2PIR(filename_fasta)

    m = MyLoop(env=env,
               alnfile=filename_pir,
               knowns="PROT",
               filename=filename_pdb,
               sequence=pdb_code.upper(),  # code of the target
               loop_assess_methods=_modeller.automodel.assess.DOPE)

    m.auto_align()  # get an automatic alignment
    m.make()

    return fixModellerPDB(m, add_missing_atoms=add_missing_atoms, filename_output=pdb_code + "_modeller.pdb")


def FASTA2PIR(filename_input):
    file = open(filename_input).readlines()
    pdb_code = ""
    new_file = []

    for line in file:
        if not pdb_code:
            match = _re.search(r">([\w]{4})", line)
            if match:
                pdb_code = match.group(1)
            new_file = [">P1;%s\nsequence:::::::::\n" % pdb_code]
        else:
            if _re.search(r">[\w]{4}", line):
                new_file += "/"
            else:
                new_file += line

    filename_output = pdb_code + ".pir"
    with open(filename_output, "w") as file:
        file.write(">P1;PROT\nstructureX:%s:FIRST:@ END::::::\n*\n\n" % pdb_code)
        for i, line in enumerate(new_file):
            if i != len(new_file) - 1:
                file.write(line)
            else:
                file.write(line[:-1])
        file.write("*")

    return filename_output


def fixModellerPDB(model, add_missing_atoms, filename_output=None):
    def getAndFixResidue(idx, missing_residue):
        nonlocal pdb_modified
        modified_residue = _copy.deepcopy(pdb_modified.filter("resSeq==%d" % idx)[0])
        modified_residue.chainID = missing_residue.chainID
        modified_residue.resSeq = missing_residue.resSeq
        modified_residue.iCode = missing_residue.iCode
        return modified_residue

    def equalChainIDs(res_orig, res_mod):
        return True if (res_orig.chainID == "A" and res_mod.chainID == " ") or res_orig.sameChain(res_mod) else False

    try:
        filename_modified = model.loop.outputs[0]['name']
    except:
        print("Modeller failed to create file. Please see the log.")
        return -1

    pdb_modified = _IO.PDB.PDB(filename_modified)
    for missing_residue in model.pdb._missing_residues:
        index = 1
        breakloops = False
        for i, chain in enumerate(model.pdb):
            if chain.type == "chain":
                for j, residue in enumerate(chain):
                    if j == 0 and equalChainIDs(missing_residue, residue) and missing_residue < residue:
                        chain.insert(j, getAndFixResidue(index, missing_residue))
                        breakloops = True
                    elif j == len(chain) - 1 and equalChainIDs(missing_residue, residue) and residue < missing_residue:
                        chain.insert(j + 1, getAndFixResidue(index, missing_residue))
                        breakloops = True
                    elif 0 < j < len(chain) and chain[j - 1] < missing_residue < chain[j]:
                        chain.insert(j, getAndFixResidue(index, missing_residue))
                        breakloops = True
                    if residue.type == "amino_acid":
                        index += 1
                    if breakloops:
                        break
            if breakloops:
                break
    model.pdb._missing_residues = []

    if add_missing_atoms:
        for missing_atom in model.pdb._missing_atoms:
            index = 1
            breakloops = False
            for i, chain in enumerate(model.pdb):
                if chain.type == "chain":
                    for j, residue in enumerate(chain):
                        if residue == missing_atom:
                            residue_new = getAndFixResidue(index, missing_atom)
                            residue.clear()
                            residue.__init__(residue_new)
                        if residue.type == "amino acid":
                            index += 1
                        if breakloops:
                            break
                if breakloops:
                    break
        model.pdb._missing_atoms = []

    model.pdb.reNumberAtoms()
    if filename_output is None:
        filename_output = _os.path.splitext(model.pdb.filename)[0] + "_modified.pdb"
    model.pdb.writePDB(filename_output)

    return _os.path.abspath(filename_output)
