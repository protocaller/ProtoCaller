import os, sys
# TODO: need to append the folder where ProtoCaller is installed. This is going to be fixed in a later version
sys.path.append("/home/ms2m18/")
# add an alternative default version for GROMACS. Otherwise, use bash default
# os.environ["GROMACSHOME"] = "/home/ms2m18/gromacs-2018.4"

from ProtoCaller.Utils.fileio import Dir
import ProtoCaller.Ensemble as Ensemble
from ProtoCaller.Ensemble import Ligand, Protein

with Dir("T4-lysozyme", overwrite=True):
    # create a protein from its PDB code and the residue number of the ligand we are going to use for mapping
    protein = Protein("181L", ligand_ref="400")
    # create two ligands from SMILES strings and name them
    benzol = Ligand("C1=CC=CC=C1", name="benzol", workdir="Ligands")
    o_xylene = Ligand("CC1=CC=CC=C1C", name="o-xylene", workdir="Ligands")
    # create the morphs from the ligands
    morphs = [[benzol, o_xylene], [o_xylene, benzol]]

    # create a system from the protein and the morphs and set up some default settings regarding system preparation
    system = Ensemble.Ensemble("GROMACS", protein=protein, morphs=morphs, box_length_complex=7, ligand_ff="gaff2",
                               workdir=protein.workdir.path)
    # only keep the reference ligand and keep all crystallographic waters
    system.protein.filter(ligands=None, waters="all")
    # protonate the protein using PDB2PQR and add missing atoms / residues if needed
    system.protein.prepare()
    # prepare the complex and solvated leg starting structures and save them as GROMACS input
    system.prepareComplexes()