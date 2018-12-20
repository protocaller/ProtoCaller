import os as _os
import tempfile as _tempfile

import BioSimSpace as _BSS
import parmed as _pmd
import Sire.Maths as _SireMaths
import Sire.Vol as _SireVol

import ProtoCaller as _PC
import ProtoCaller.Parametrise as _parametrise
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Utils.runexternal as _runexternal
import ProtoCaller.Utils.stdio as _stdio

def solvate(complex, params, box_length=8, shell=0, neutralise=True, ion_conc=0.154, centre=True, work_dir=None,
            filebase="complex"):
    if work_dir is None:
        work_dir = _os.path.basename(_tempfile.TemporaryDirectory().name)
        temp = True
    else:
        temp = False
    with _fileio.Dir(dirname=work_dir, temp=temp):
        """CENTERING"""

        box_length = box_length
        if centre:
            box = complex._getAABox()
            min_coords, max_coords = box.minCoords(), box.maxCoords()
            centre = (min_coords + max_coords) / 2
            difference = max_coords - min_coords
            if any([x / 10 > box_length for x in difference]):
                box_length = int(max(difference / 10)) + 1
                print("Insufficient input box size. Changing to a box length of %d nm..." % box_length)
                complex.translate(tuple(-centre + _SireMaths.Vector(3 * (box_length * 5,))))
                complex._sire_system.setProperty("space",
                                                 _SireVol.PeriodicBox(_SireMaths.Vector(3 * (box_length * 10,))))

        """SOLVATING"""

        _BSS.IO.saveMolecules(filebase, complex, "Gro87,GroTop")
        files = [filebase + ".gro", filebase + ".top"]
        _os.rename(filebase + ".gro87", files[0])
        _os.rename(filebase + ".grotop", files[1])

        new_gro = filebase + "_solvated.gro"
        command = "{0} solvate -shell {1} -box {2} {2} {2} -cp {3} -o {4}".format(_PC.GROMACSEXE, shell,
                                                                                  box_length, files[0], new_gro)

        try:
            """CREATE UNPARAMETRISED WATERS AND LOAD THEM INTO MEMORY"""

            _runexternal.runExternal(command, procname="gmx solvate")
            complex_solvated = _pmd.load_file(new_gro)
            waters = complex_solvated[":SOL"]

            """PREPARE WATERS FOR TLEAP AND PARAMETRISE"""

            for residue in waters.residues:
                residue.name = "WAT"
            for i, atom in enumerate(waters.atoms):
                if "H" in atom.name:
                    atom.name = atom.name[0] + atom.name[2]
                elif "O" in atom.name:
                    atom.name = "O"

            waters.save(filebase + "_waters.pdb")
            waters_prep = _parametrise.parametriseAndLoadPmd(params, filebase + "_waters.pdb", "water")
            waters_prep.box = _stdio.ignore_warnings(_pmd.load_file)(files[1], xyz=files[0], skip_bonds=True).box
            for residue in waters_prep.residues:
                residue.name = "SOL"
            waters_prep_filenames = [filebase + "_waters.gro", filebase + "_waters.top"]
            for filename in waters_prep_filenames:
                waters_prep.save(filename)

        except:
            print("Solvation failed. Returning original system...")
            return complex

        """ADD IONS"""

        if any([neutralise, ion_conc, shell]):
            """WRITE MDP FILE"""

            with open("ions.mdp", "w") as file:
                file.write("; Neighbour searching\n")
                file.write("cutoff-scheme           = Verlet\n")
                file.write("rlist                   = 1.1\n")
                file.write("pbc                     = xyz\n")
                file.write("verlet-buffer-tolerance = -1\n")
                file.write("\n; Electrostatics\n")
                file.write("coulombtype             = cut-off\n")
                file.write("\n; VdW\n")
                file.write("rvdw                    = 1.0\n")

            charge = round(complex.charge().magnitude()) if neutralise else 0
            n_Na, n_Cl = [int((box_length * 10 ** -8) ** 3 * 6.022 * 10 ** 23 * ion_conc) - abs(charge) // 2] * 2
            if neutralise:
                if charge < 0:
                    n_Na -= charge
                else:
                    n_Cl += charge

            try:
                ions_prep_filenames = [filebase + "_ions.gro", filebase + "_ions.top"]
                command = "{0} grompp -f ions.mdp -c {1} -p {2} -o {3}_solvated.tpr".format(
                    _PC.GROMACSEXE, *waters_prep_filenames, filebase)
                _runexternal.runExternal(command, procname="gmx grompp")

                command = "{{ echo 2; }} | {0} genion -s {1}_solvated.tpr -o {2} -nn {3} -np {4}".format(
                    _PC.GROMACSEXE, filebase, ions_prep_filenames[0], n_Cl, n_Na)
                _runexternal.runExternal(command, procname="gmx genion")

                """PREPARE WATERS FOR TLEAP AND PARAMETRISE"""
                ions = _pmd.load_file(ions_prep_filenames[0])
                for residue in ions.residues:
                    if residue.name == "SOL":
                        residue.name = "WAT"

                ions.save(filebase + "_ions.pdb")
                _os.remove(ions_prep_filenames[0])
                ions_prep = _parametrise.parametriseAndLoadPmd(params, filebase + "_ions.pdb", "water")

                for filename in ions_prep_filenames:
                    ions_prep.save(filename)
            except:
                print("Ion addition failed. Returning solvated system...")
                return complex + _BSS.IO.readMolecules(waters_prep_filenames)
        return complex + _BSS.IO.readMolecules(ions_prep_filenames)