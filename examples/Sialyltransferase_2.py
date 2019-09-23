# import os
# add an alternative default version for GROMACS. Otherwise, use bash default
# os.environ["GROMACSHOME"] = os.path.expanduser("~/gromacs-2018.4")
import logging
logging.basicConfig(level=logging.INFO)

from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Ensemble import Ligand, Protein, Ensemble

with Dir("Sialyltransferase", overwrite=False):
    # create a protein and ligands using custom files
    lig_ref = Ligand("lig_ref.sdf", protonated=True, workdir="Ligands",
                     name="lig_ref")
    lig1 = Ligand("lig1.sdf", protonated=True, workdir="Ligands", name="lig1")
    lig2 = Ligand("lig2.sdf", protonated=True, workdir="Ligands", name="lig2")
    protein = Protein("2WNB", pdb_file="protein.pdb", ligand_ref=lig_ref)

    # create the morphs from the ligands
    morphs = [[lig1, lig2], [lig2, lig1]]

    # create a system from the protein and the morphs and set up some default
    # settings regarding system preparation
    system = Ensemble("GROMACS", protein=protein, morphs=morphs,
                      box_length_complex=7, ligand_ff="gaff2",
                      workdir=protein.workdir.path)
    # only keep the reference ligand
    system.protein.filter(ligands=None)
    # we still need PDB2PQR to give proper naming to residues
    system.protein.prepare(protonate_proteins="pdb2pqr")
    # prepare the complex and solvated leg starting structures and save them as
    # GROMACS input. Parametrisation here will be very slow. Be patient
    system.prepareComplexes()
