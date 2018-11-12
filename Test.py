import ensemble as _ensemble
import subdir as _subdir


'''
import Simulation as _Simulation
import os as _os
_os.chdir("Temp/1BJI/Pair 1")
run = _Simulation.GMXSerialRuns("1BJI_run", "complex_final.gro", "complex_final.top", coulomb_lambdas = [0.0, 0.5, 1.0, 1.0, 1.0, 1.0], vdw_lambdas = [0.0, 0.0, 0.0, 0.5, 0.8, 1.0])
run.addProtocol(name="Minimisation", use_preset="minimisation")
run.addProtocol(name="Equilibration_NVT", use_preset="equilibration_nvt")
run.addProtocol(name="Equilibration_NPT", use_preset="equilibration_npt")
run.addProtocol(name="Production", use_preset="production", n_steps = 100000)
run.runSimulations()

exit()
'''

#proteins = ["1BJI", "1A0O", "1AZ5"]

#e.g. Neuraminidase
with _subdir.Subdir("Temp", overwrite=True):
    protein = "1BJI"
    morphs = [["C1(=CC=CC=C1)CN(CC)CC", "C1(=CC=CC=C1)CN(CC)C"]]
    ligand_id = "479G"
    ensemble = _ensemble.Ensemble("GROMACS", protein=protein, morphs=morphs, ligand_id=ligand_id, box_length=6)
    ensemble.filterPDB(ligands=None, waters="site")
    ensemble.preparePDB(add_missing_residues=None)
    ensemble.parametrisePDB()
    ensemble.parametriseLigands()
    ensemble.prepareComplexes()
