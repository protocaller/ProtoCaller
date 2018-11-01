import const as _const

import os as _os

import BioSimSpace as _BSS

import ensemble as _ensemble
import subdir as _subdir
import Protocol as _protocol

prot = _protocol.Protocol(integrator="stochastic", thermostat="berendsen")
prot.write("GROMACS")

exit()
proteins = ["1BJI", "1A0O", "1AZ5"]
ligands = ["C1(=CC=CC=C1)CN(CC)CC", "C1(=CC=CC=C1)CN(CC)C"]
proteins = ["1BJI"]
complexes = [[*proteins, *ligands]]

with _subdir.Subdir("Temp", overwrite=True):
    protein = "1BJI"
    morphs = [["C1(=CC=CC=C1)CN(CC)CC", "C1(=CC=CC=C1)CN(CC)C"]]
    ensemble = _ensemble.Ensemble("GROMACS", protein=protein, morphs=morphs, ligand_id="479G", box_length=6)
    ensemble.filterPDB(ligands=None, waters="site")
    ensemble.preparePDB(add_missing_residues="")
    ensemble.parametrisePDB()
    ensemble.prepareComplexes()