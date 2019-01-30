import ProtoCaller.IO.PDB as PDB

import tempfile


def test_read_write_1bji():
    obj = PDB.PDB("1bji.pdb")
    filestr = open("1bji.pdb").readlines()
    filestr = [line.strip() for line in filestr
               if line[:3] in ["TER", "END"] or line[:4] == "ATOM" or line[:6] == "HETATM"]

    assert len(obj._missing_atoms) == 0
    assert len(obj._missing_residues) == 0
    assert len(obj._disulfide_bonds) == 9
    assert len(obj._site_residues) == 74
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
    assert len(obj_new._missing_atoms) == 0
    assert len(obj_new._missing_residues) == 0
    assert len(obj_new._disulfide_bonds) == 9
    assert len(obj_new._site_residues) == 74
    assert len(obj_new[0]) == 388
    assert len(obj_new[1]) == 202

    # the ultimate test: rewriting the PDB should result in the same PDB
    assert filestr_new == filestr
