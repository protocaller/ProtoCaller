# import os
# add an alternative default version for GROMACS. Otherwise, use bash default
# os.environ["GROMACSHOME"] = os.path.expanduser("~/gromacs-2018.4")
import logging
logging.basicConfig(level=logging.INFO)

from BioSimSpace._SireWrappers import System

from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Ensemble import Ligand, Perturbation
from ProtoCaller.Solvate import solvate
from ProtoCaller.IO.GROMACS import saveAsGromacs
from ProtoCaller.Parametrise import Params
from ProtoCaller.Wrappers.biosimspacewrapper import resize

with Dir("Relative_Solvation", overwrite=True):
    with Dir("Ligands"):
        # create two ligands from SMILES strings and name them
        toluene = Ligand("CC1=CC=CC=C1C", name="toluene", workdir="Ligands")
        benzene = Ligand("C1=CC=CC=C1", name="benzene", workdir="Ligands")

        # protonate parametrise the ligands at pH 7 with GAFF2
        for ligand in [toluene, benzene]:
            ligand.protonate(babel_parameters={"pH": 7.0})
            ligand.parametrise(Params(ligand_ff="gaff2"))

    with Dir("toluene~benzene"):
        # create the perturbation from the ligands
        perturbation = Perturbation(toluene, benzene)

        # we create the mixed topology object by using the toluene as a reference ligand
        perturbation.alignAndCreateMorph(toluene)

        # get the morph as a BioSimSpace molecule
        morph_vac = perturbation.morph
        # promote the molecule to a BioSimSpace system
        morph_vac_system = System(morph_vac)
        # give a dummy box length of e.g. 10 nm
        morph_vac_system = resize(morph_vac_system, 10)
        # save the vacuum leg as .gro and .top files
        saveAsGromacs("perturbation_vac", morph_vac_system)

        # solvate and save the solvated leg
        morph_soln = solvate(morph_vac, Params(water_ff="tip3p"), box_length=4, ion_conc=0.154)
        saveAsGromacs("perturbation_soln", morph_soln)
