import const as _const

import os as _os

import BioSimSpace as _BSS

import ensemble as _ensemble
import subdir as _subdir

"""
_os.chdir("./Temp/1BJI/Pair 1")
#system = _BSS.IO.readMolecules(["complex_mod.top", "complex_mod.gro"])

ligand1_prep = _BSS.IO.readMolecules(["../ligand_1.inpcrd", "../ligand_1.prmtop"]).getMolecules()[0]
ligand2_prep = _BSS.IO.readMolecules(["../ligand_2.inpcrd", "../ligand_2.prmtop"]).getMolecules()[0]
ligand_ref = _BSS.IO.readMolecules("../1BJI_DPC_A_479G_split.mol2").getMolecules()[0]
complex_template = _BSS.IO.readMolecules(["../complex_template.gro", "../complex_template.top"])

ligand1_prep = _BSS.Align.rmsdAlign(ligand1_prep, ligand_ref)
ligand2_prep = _BSS.Align.rmsdAlign(ligand2_prep, ligand1_prep)
morph = _BSS.Align.merge(ligand1_prep, ligand2_prep)

_BSS.IO.saveMolecules("morph", morph, "Gro87")

complex = complex_template + morph

print("Solvating...")
while True:
    try:
        complex = _BSS.Solvent.solvate("tip3p", molecule=complex, box=3 * [8 * _BSS.Units.Length.nanometer], ion_conc=0.154)
        break
    except ValueError:
        print("Err")
        break

print("Saving as GROMACS...")
_BSS.IO.saveMolecules("complex", complex, "Gro87")
_BSS.IO.saveMolecules("complex", complex, "GroTop")
import parmed as _pmd
file = _pmd.load_file("complex.grotop", xyz="complex.gro87")
try:
    _os.remove("complex.gro")
except:
    pass
file.save("complex.gro")
print("Saving as AMBER...")
_BSS.IO.saveMolecules("complex", complex, "PRM7,RST7")
exit()
_BSS.Solvent.solvate("tip3p", molecule=system, box=3 * [12 * _BSS.Units.Length.nanometer], ion_conc=0.154)
exit()
"""
proteins = ["1BJI", "1A0O", "1AZ5"]
ligands = ["C1(=CC=CC=C1)CN(CC)CC", "C1(=CC=CC=C1)CN(CC)C"]
proteins = ["1BJI"]
complexes = [[*proteins, *ligands]]

with _subdir.Subdir("Temp1", overwrite=True):
    protein = "1BJI"
    morphs = [["C1(=CC=CC=C1)CN(CC)CC", "C1(=CC=CC=C1)CN(CC)C"]]
    ensemble = _ensemble.Ensemble("GROMACS", protein=protein, morphs=morphs, ligand_id="479G", box_length=6)
    ensemble.filterPDB(ligands=None, waters="site")
    ensemble.preparePDB(add_missing_residues="")
    ensemble.parametrisePDB()
    ensemble.prepareComplexes()