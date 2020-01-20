import logging
logging.basicConfig(level=logging.INFO)

import glob

from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Ensemble import Ensemble, Protein, Ligand
from ProtoCaller.Parametrise import Params

params = Params(ligand_ff="GAFF2")

with Dir("DHFR", overwrite=False):
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

            # there is no point in protonating and parametrising the reference
            if parametrised_files is None and sdf != "reference.sdf":
                ligands[ligname].parametrise(params=params)

    # all of the perturbations
    morphs = [
        [ligands["Me_H"], ligands["Me_Me"]],
        [ligands["Me_Me"], ligands["H_H"]],
        [ligands["Me_H"], ligands["H_H"]],
        [ligands["Cl_Cl"], ligands["Cl_H"]],
        [ligands["Cl_H"], ligands["H_H"]],
        [ligands["Me_Me"], ligands["Cl_H"]],
        [ligands["Me_H"], ligands["Cl_Cl"]],
        [ligands["NO2"], ligands["H_H"]],
        [ligands["Br"], ligands["H_H"]],
        [ligands["NO2"], ligands["Br"]],
        [ligands["OMe"], ligands["H_H"]],
        [ligands["OMe"], ligands["Br"]],
    ]

    # now we generate all input structures for 5HPB
    print("Generating input files for 5HPB...")
    # here we initialise the Ensemble class which is responsible for handling
    # all perturbations
    system = Ensemble(protein="5HPB", morphs=morphs, box_length_complex=8,
                      workdir="5HPB", ligand_ref="202")
    # we discard everything apart from the crystal waters
    system.protein.filter(ligands=None, simple_anions=None, simple_cations=None,
                          complex_anions=None, waters="all")
    system.protein.prepare()
    # save all merged structures
    system.prepareComplexes()

    print("Generating input files for: 6DAV with altLoc A")
    # here we supply a custom ligand orientation as a reference
    protein_A = Protein("6DAV", ligand_ref=ligands["reference"],
                        workdir="6DAV_A")
    # delete any atoms with altLoc B
    for chain in protein_A.pdb_obj:
        for residue in chain:
            atoms_to_purge = [x for x in residue if x.altLoc == "B"]
            residue.purgeAtoms(atoms_to_purge, "discard")
    protein_A.pdb_obj.writePDB()

    system_A = Ensemble(protein=protein_A, morphs=morphs,
                        box_length_complex=8, workdir="6DAV_A")
    # here we only keep chain A and the crystal waters that belong to it
    system_A.protein.filter(ligands=None, simple_anions=None,
                            simple_cations=None, waters="chains",
                            chains=["A"])
    system_A.protein.prepare()
    system_A.prepareComplexes()

    # repeat but for altLoc B
    print("Generating input files for: 6DAV with altLoc B")
    protein_B = Protein("6DAV", ligand_ref=ligands["reference"],
                        workdir="6DAV_B")
    for chain in protein_B.pdb_obj:
        for residue in chain:
            atoms_to_purge = [x for x in residue if x.altLoc == "A"]
            residue.purgeAtoms(atoms_to_purge, "discard")
    protein_B.pdb_obj.writePDB()

    system_B = Ensemble(protein=protein_B, morphs=morphs,
                        box_length_complex=8, workdir="6DAV_B")
    system_B.protein.filter(ligands=None, simple_anions=None,
                            simple_cations=None, waters="chains",
                            chains=["A"])
    system_B.protein.prepare()
    system_B.prepareComplexes()
