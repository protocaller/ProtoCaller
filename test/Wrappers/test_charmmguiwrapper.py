import ProtoCaller as PC
import ProtoCaller.IO.PDB as PDB
from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Wrappers.charmmguiwrapper import charmmguiTransform

from shutil import copyfile


def test_model_3zg0():
    with Dir(PC.TESTDIR + "/shared"):
        with Dir("temp", temp=True):
            copyfile("../3ZG0.pdb", "3ZG0.pdb")

            obj = PDB.PDB("3ZG0.pdb")
            obj_mod = PDB.PDB(charmmguiTransform("3ZG0.pdb"))

            assert len(obj.missing_atoms) == len(obj_mod.missing_atoms) == 0
            assert len(obj.missing_residues) == 7
            assert len(obj_mod.missing_residues) == 0
            assert len(obj.disulfide_bonds) == len(obj_mod.disulfide_bonds) == 0
            assert len(obj.modified_residues) == len(obj_mod.modified_residues) == 0
            assert len(obj.site_residues) == len(obj_mod.site_residues) == 120
            assert len(obj_mod[0]) == 642
            assert len(obj_mod[1]) == 642
            assert len(obj_mod[2]) == 9
            assert len(obj_mod[3]) == 9
            assert len(obj_mod[4]) == 239
            assert len(obj_mod[5]) == 171
