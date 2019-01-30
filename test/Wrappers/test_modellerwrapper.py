import ProtoCaller as PC
import ProtoCaller.IO.PDB as PDB
from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Wrappers.modellerwrapper import modellerTransform

from shutil import copyfile


def test_protonate_3zg0():
    with Dir(PC.TESTDIR + "/shared"):
        with Dir("temp", temp=True):
            copyfile("../3ZG0.pdb", "3ZG0.pdb")
            copyfile("../3ZG0.fasta", "3ZG0.fasta")

            obj = PDB.PDB("3ZG0.pdb")
            obj_modeller = PDB.PDB(modellerTransform("3ZG0.pdb", "3ZG0.fasta", add_missing_atoms=False))

            assert len(obj._missing_atoms) == len(obj_modeller._missing_atoms) == 0
            assert len(obj._missing_residues) == 7
            assert len(obj_modeller._missing_residues) == 0
            assert len(obj._disulfide_bonds) == len(obj_modeller._disulfide_bonds) == 0
            assert len(obj._modified_residues) == len(obj_modeller._modified_residues) == 0
            assert len(obj._site_residues) == len(obj_modeller._site_residues) == 120
            assert len(obj_modeller[0]) == 642
            assert len(obj_modeller[1]) == 642
            assert len(obj_modeller[2]) == 9
            assert len(obj_modeller[3]) == 9
            assert len(obj_modeller[4]) == 239
            assert len(obj_modeller[5]) == 171