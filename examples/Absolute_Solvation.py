# import os
# add an alternative default version for GROMACS. Otherwise, use bash default
# os.environ["GROMACSHOME"] = os.path.expanduser("~/gromacs-2018.4")
import logging
logging.basicConfig(level=logging.INFO)

from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Ensemble import Ligand
from ProtoCaller.Solvate import solvate
from ProtoCaller.IO.GROMACS import saveAsGromacs
from ProtoCaller.Parametrise import Params
from ProtoCaller.Wrappers.parmedwrapper import openFilesAsParmed, resize

with Dir("Absolute_Solvation", overwrite=True):
    with Dir("Ligands"):
        # create two ligands from SMILES strings and name them
        toluene = Ligand("CC1=CC=CC=C1C", name="toluene", workdir="Ligands")
        benzene = Ligand("C1=CC=CC=C1", name="benzene", workdir="Ligands")

        # protonate parametrise the ligands at pH 7 with GAFF2
        for ligand in [toluene, benzene]:
            ligand.protonate(babel_parameters={"pH": 7.0})
            ligand.parametrise(Params(ligand_ff="gaff2"))

    for ligand in [toluene, benzene]:
        with Dir(ligand.name):
            # open the parametrised ligand as a ParmEd object
            lig_vac = openFilesAsParmed(ligand.parametrised_files)
            # give a dummy box length of e.g. 10 nm
            lig_vac = resize(lig_vac, 10)
            # save the vacuum leg as .gro and .top files
            saveAsGromacs(ligand.name + "_vac", lig_vac)

            # solvate and save the solvated leg
            lig_solv = solvate(lig_vac, Params(water_ff="tip3p"), box_length=4, ion_conc=0.154)
            saveAsGromacs(ligand.name + "_solv", lig_solv)
