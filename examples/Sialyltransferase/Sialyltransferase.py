# import os
# add an alternative default version for GROMACS. Otherwise, use bash default
# os.environ["GROMACSHOME"] = os.path.expanduser("~/gromacs-2018.4")
import logging
logging.basicConfig(level=logging.INFO)

from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Ensemble import Ligand, Protein, Ensemble

with Dir("Sialyltransferase", overwrite=False):
    # create a protein and ligands using coustom files
    lig_ref = Ligand("../../lig.sdf", protonated=True,
                                   workdir="Ligands", name="lig_ref")
    lig1 = Ligand("../../in1.sdf", protonated=True,
                                   workdir="Ligands", name="lig1")
    lig2 = Ligand("../../in2.sdf", protonated=True,
                                   workdir="Ligands", name="lig2")
    protein = Protein("2WNB", pdb_file="../../pro.pdb", ligand_ref=lig_ref)
    protein.pdb_obj[0].chainID = "A"

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


    # create the morphs from the ligands
    morphs = [[lig1, lig2], [lig2, lig1]]

    # create a system from the protein and the morphs and set up some default
    # settings regarding system preparation
    system = Ensemble("GROMACS", protein=protein, morphs=morphs,
                      box_length_complex=7, ligand_ff="gaff2",
                      workdir=protein.workdir.path)
    # only keep the reference ligand and keep all crystallographic waters
    system.protein.filter(ligands=None, waters="all")
    # add missing residues with Modeller and add missing atoms and protonate
    # with PDB2PQR
    system.protein.prepare(add_missing_residues="modeller",
                           add_missing_atoms="pdb2pqr")
    # prepare the complex and solvated leg starting structures and save them as
    # GROMACS input. Parametrisation here will be very slow. Be patient
    system.prepareComplexes()
