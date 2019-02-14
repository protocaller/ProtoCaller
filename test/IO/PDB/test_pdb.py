import ProtoCaller as PC
import ProtoCaller.IO.PDB as PDB
from ProtoCaller.Utils.fileio import Dir

import copy
import tempfile


def test_read_write_1bji():
    with Dir(PC.TESTDIR + "/shared"):
        obj = PDB.PDB("1bji.pdb")
        filestr = open("1bji.pdb").readlines()
        filestr = [line.strip() for line in filestr
                   if line[:3] in ["TER", "END"] or line[:4] == "ATOM" or line[:6] == "HETATM"]

        assert obj.numberOfAtoms == 3398
        assert len(obj.missing_atoms) == 0
        assert len(obj.missing_residues) == 0
        assert len(obj.modified_residues) == 3
        assert len(obj.disulfide_bonds) == 9
        assert len(obj.site_residues) == 74
        assert len(obj[0]) == 388
        assert len(obj[1]) == 202

        new_file = tempfile.NamedTemporaryFile()
        obj.writePDB(new_file.name)

        filestr_new = open(new_file.name).readlines()
        # filter out only the relevant sections
        filestr_new = [line for line in filestr_new
                       if line[:3] in ["TER", "END"] or line[:4] == "ATOM" or line[:6] == "HETATM"]
        # fix a known formatting difference with the original PDB
        filestr_new = [line[:12] + " " + line[12:16] + line[17:] for line in filestr_new]
        filestr_new = [line.strip() for line in filestr_new]

        obj_new = PDB.PDB(new_file.name)
        assert obj_new.numberOfAtoms == 3398
        assert len(obj_new.missing_atoms) == 0
        assert len(obj_new.missing_residues) == 0
        assert len(obj.modified_residues) == 3
        assert len(obj_new.disulfide_bonds) == 9
        assert len(obj_new.site_residues) == 74
        assert len(obj_new[0]) == 388
        assert len(obj_new[1]) == 202

        # the ultimate test: rewriting the PDB should result in the same PDB
        assert filestr_new == filestr


def test_copy_1bji():
    with Dir(PC.TESTDIR + "/shared"):
        obj = PDB.PDB("1bji.pdb")
        obj2 = copy.copy(obj)

        assert obj2.numberOfAtoms == 3398
        assert len(obj2.missing_atoms) == 0
        assert len(obj2.missing_residues) == 0
        assert len(obj2.modified_residues) == 3
        assert len(obj2.disulfide_bonds) == 9
        assert len(obj2.site_residues) == 74
        assert len(obj2[0]) == 388
        assert len(obj2[1]) == 202

        obj2.purgeResidues(obj2.filter("type=='amino_acid'"), mode="keep")
        assert obj.numberOfAtoms == obj2.numberOfAtoms == 3067


def test_deepcopy_1bji():
    with Dir(PC.TESTDIR + "/shared"):
        obj = PDB.PDB("1bji.pdb")
        obj2 = copy.deepcopy(obj)

        assert obj2.numberOfAtoms == 3398
        assert len(obj2.missing_atoms) == 0
        assert len(obj2.missing_residues) == 0
        assert len(obj2.modified_residues) == 3
        assert len(obj2.disulfide_bonds) == 9
        assert len(obj2.site_residues) == 74
        assert len(obj2[0]) == 388
        assert len(obj2[1]) == 202

        obj2.purgeResidues(obj2.filter("type=='amino_acid'"), mode="keep")
        assert 3067 == obj2.numberOfAtoms != obj.numberOfAtoms
