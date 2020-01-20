import logging
logging.basicConfig(level=logging.INFO)

import glob

from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Ensemble import Ensemble, Ligand
from ProtoCaller.Parametrise import Params

params = Params(ligand_ff="GAFF2")

with Dir("FXA", overwrite=False):
    # create a ligand dictionary and populate it from the Ligands folder.
    # this block makes sure we don't reparametrise if we rerun the script.
    ligands = {}
    with Dir("Ligands") as ligdir:
        lignames = [x.split(".")[0] for x in glob.glob("*.sdf")]
        for ligname in lignames:
            sdf = "{}.sdf".format(ligname)
            prmtop = "{}.prmtop".format(ligname)
            inpcrd = "{}.inpcrd".format(ligname)

            if len(glob.glob(prmtop)) and len(glob.glob(inpcrd)):
                parametrised_files = [prmtop, inpcrd]
            else :
                parametrised_files = None

            # here we initialise the ligand with custom structures.
            ligands[ligname] = Ligand(sdf, protonated=True, minimise=False,
                                      parametrised_files=parametrised_files,
                                      name=ligname, workdir=ligdir.path)

            if parametrised_files is None:
                ligands[ligname].parametrise(params=params)

    # all of the perturbations
    morphs = [
        [ligands["Me_H"], ligands["H_H"]],
        [ligands["Me_Me"], ligands["H_H"]],
        [ligands["Me_Me"], ligands["Me_H"]],
        [ligands["iPr_H"], ligands["Et_H"]],
        [ligands["iPr_H"], ligands["Me_H"]],
        [ligands["Et_Me"], ligands["Me_H"]],
        [ligands["Et_Me"], ligands["Et_H"]],
        [ligands["Et_Et"], ligands["Et_Me"]],
    ]

    print("Generating input files for 2P3T")
    # here ProtoCaller detects the reference ligand automatically
    # (there is only one possible molecule in the PDB file which can be classed
    # as a reference)
    system = Ensemble(protein="2P3T", morphs=morphs, box_length_complex=8,
                      workdir="2P3T")
    # we keep the crystal waters and the calcium ion
    system.protein.filter(ligands=None, simple_anions=None, simple_cations="all",
                          complex_anions=None, waters="all")
    # we have to use Modeller to treat the missing atoms and residues
    system.protein.prepare(add_missing_residues="modeller",
                           add_missing_atoms="modeller")
    # we shrink the dummy bonds to half the initial length
    system.prepareComplexes(scale_dummy_bonds=0.5)
