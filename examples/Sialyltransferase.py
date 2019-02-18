import os, sys
# TODO: need to append the folder where ProtoCaller is installed. This is going to be fixed in a later version
sys.path.append(os.path.expanduser("~/ProtoCaller"))
# add an alternative default version for GROMACS. Otherwise, use bash default
# os.environ["GROMACSHOME"] = os.path.expanduser("~/gromacs-2018.4")

from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Ensemble import Ligand, Protein, Ensemble

with Dir("Sialyltransferase", overwrite=True):
    # create a protein from its PDB code and the residue number of the ligand we are going to use for mapping
    protein = Protein("2WNB", ligand_ref="1344")

    # change the selenomethionine residues to methionine
    for residue in protein.pdb_obj.modified_residues:
        residue.resName = "MET"
        for atom in residue:
            if atom.name == "SE":
                atom.name = "SD"
                atom.element = "S"

    # delete any atoms with altLoc B
    for chain in protein.pdb_obj:
        for residue in chain:
            atoms_to_purge = [x for x in residue if x.altLoc == "B"]
            residue.purgeAtoms(atoms_to_purge, "discard")
    protein.pdb_obj.writePDB()

    # create two ligands from SMILES strings
    lig1 = Ligand("O([P@]([O-])(=O)CC[C@@]1(C(=O)N)C[C@@H]([C@H]([C@H]([C@@H]([C@@H](CO)O)O)O1)NC(=O)C)O)C"
                  "[C@H]1O[C@@H](n2c(=O)nc(N)cc2)[C@H](O)[C@@H]1O", workdir="Ligands")
    lig2 = Ligand("[P@]([O-])(=O)(OC[C@H]1O[C@@H](n2ccc(nc2=O)N)[C@H](O)[C@@H]1O)CC[C@]1(CO)O[C@H]([C@@H]"
                  "([C@H](C1)O)NC(=O)C)[C@H](O)[C@H](O)CO", workdir="Ligands")
    # create the morphs from the ligands
    morphs = [[lig1, lig2], [lig2, lig1]]

    # create a system from the protein and the morphs and set up some default settings regarding system preparation
    system = Ensemble("GROMACS", protein=protein, morphs=morphs, box_length_complex=7, ligand_ff="gaff2",
                      workdir=protein.workdir.path)
    # only keep the reference ligand and keep all crystallographic waters
    system.protein.filter(ligands=None, waters="all")
    # add missing residues with Modeller and add missing atoms and protonate with PDB2PQR
    protein.prepare(add_missing_residues="modeller", add_missing_atoms="pdb2pqr")
    # prepare the complex and solvated leg starting structures and save them as GROMACS input
    # parametrisation here will be very slow. Be patient
    system.prepareComplexes()
