import ProtoCaller as PC
import ProtoCaller.IO.PDB as PDB
from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Wrappers.PDB2PQRwrapper import PDB2PQRtransform

# residues with their correct protonated atom length
size_dict = {
                "ALA": 10,
                "ARG": 24,
                "ASH": 13,
                "ASN": 14,
                "ASP": 12,
                "CYS": 11,
                "CYX": 10,
                "GLH": 16,
                "GLN": 17,
                "GLU": 15,
                "GLY": 7,
                "HID": 17,
                "HIE": 17,
                "HIS": 17,
                "HIP": 18,
                "ILE": 19,
                "LEU": 19,
                "LYN": 21,
                "LYS": 22,
                "MET": 17,
                "PHE": 20,
                "PRO": 14,
                "SER": 11,
                "THR": 14,
                "TRP": 24,
                "TYR": 21,
                "VAL": 16,
            }


def test_protonate_1bji():
    with Dir(PC.TESTDIR + "/shared"):
        with Dir("temp", temp=True):
            file_pdb2pqr = PDB2PQRtransform("../1bji.pdb", filename_output="1bji_pdb2pqr.pdb")
            obj_protonated = PDB.PDB(file_pdb2pqr)

            assert len(obj_protonated._missing_atoms) == 0
            assert len(obj_protonated._missing_residues) == 0
            assert len(obj_protonated._disulfide_bonds) == 9
            assert len(obj_protonated._site_residues) == 74
            assert len(obj_protonated[0]) == 388
            assert len(obj_protonated[1]) == 202

            for i, residue in enumerate(obj_protonated[0]):
                # don't deal with terminal amino acids
                if i == 0 or i == len(obj_protonated[0]) - 1:
                    continue
                if residue.resName in size_dict.keys():
                    assert len(residue) == size_dict[residue.resName]
