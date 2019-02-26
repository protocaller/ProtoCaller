import ProtoCaller as _PC

import os as _os
import tempfile as _tempfile

if _PC.BIOSIMSPACE:
    import BioSimSpace as _BSS
    import ProtoCaller.Wrappers.BioSimSpacewrapper as _BSSwrap
import parmed as _pmd

import ProtoCaller.Parametrise as _parametrise
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Utils.runexternal as _runexternal
import ProtoCaller.Utils.stdio as _stdio
import ProtoCaller.Wrappers.parmedwrapper as _pmdwrap


def solvate(complex, params, box_length=8, shell=0, neutralise=True, ion_conc=0.154, centre=True, work_dir=None,
            filebase="complex"):
    """
    Uses gmx solvate and gmx genion to solvate the system and (optionally) add NaCl ions. This function preserves the
    crystal water molecules.

    Parameters
    ----------
    complex : BioSimSpace.System or parmed.structure.Structure
        The input unsolvated system.
    params : ProtoCaller.Parametrise.Params
        The input force field parameters.
    box_length : float
        Size of the box in nm. Cubic shape is assumed.
    shell : float
        Places a layer of water of the specified thickness in nm around the solute.
    neutralise : bool
        Whether to add counterions to neutralise the system.
    ion_conc : float
        Ion concentration of NaCl in mol/L.
    centre : bool
        Whether to centre the system.
    work_dir : str
        Work directory. Default: current directory.
    filebase : str
        Output base name of the file.

    Returns
    -------
    complex : BioSimSpace.System or parmed.structure.Structure
        The solvated system.
    """
    if isinstance(complex, _pmd.Structure):
        centrefunc = _pmdwrap.centre
        chargefunc = lambda x: round(sum([atom.charge for atom in x.atoms]))
        readfunc = _pmdwrap.openFilesAsParmed
    elif _PC.BIOSIMSPACE and isinstance(complex, (_BSS._SireWrappers._molecule.Molecule,
                                                  _BSS._SireWrappers._system.System)):
        if isinstance(complex, _BSS._SireWrappers._molecule.Molecule):
            complex = complex.toSystem()
        centrefunc = _BSSwrap.centre
        chargefunc = lambda x: round(x.charge().magnitude())
        readfunc = _BSS.IO.readMolecules
    else:
        raise TypeError("Cannot solvate object of type %s" % type(complex))

    if work_dir is None:
        work_dir = _os.path.basename(_tempfile.TemporaryDirectory().name)
        temp = True
    else:
        temp = False

    with _fileio.Dir(dirname=work_dir, temp=temp):
        # center
        if centre:
            complex, box_length, _ = centrefunc(complex, box_length)

        # solvate with gmx solvate and load unparametrised waters into memory
        files = _PC.IO.GROMACS.saveAsGromacs(filebase, complex)
        new_gro = filebase + "_solvated.gro"
        command = "{0} solvate -shell {1} -box {2} {2} {2} -cp \"{3}\" -o \"{4}\"".format(
            _PC.GROMACSEXE, shell, box_length, files[0], new_gro)
        _runexternal.runExternal(command, procname="gmx solvate")
        complex_solvated = _pmd.load_file(new_gro)
        waters = complex_solvated[":SOL"]

        # prepare waters for tleap and parametrise
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

        # add ions
        if any([neutralise, ion_conc, shell]):
            # write an MDP file
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

            # neutralise if needed
            charge = chargefunc(complex) if neutralise else 0
            n_Na, n_Cl = [int((box_length * 10 ** -8) ** 3 * 6.022 * 10 ** 23 * ion_conc) - abs(charge) // 2] * 2
            if neutralise:
                if charge < 0:
                    n_Na -= charge
                else:
                    n_Cl += charge

            # add ions with gmx genion
            ions_prep_filenames = [filebase + "_ions.gro", filebase + "_ions.top"]
            command = "{0} grompp -f ions.mdp -c {1} -p {2} -o \"{3}_solvated.tpr\"".format(
                _PC.GROMACSEXE, *waters_prep_filenames, filebase)
            _runexternal.runExternal(command, procname="gmx grompp")

            command = "{{ echo 2; }} | {0} genion -s \"{1}_solvated.tpr\" -o \"{2}\" -nn {3} -np {4}".format(
                _PC.GROMACSEXE, filebase, ions_prep_filenames[0], n_Cl, n_Na)
            _runexternal.runExternal(command, procname="gmx genion")

            # prepare waters for tleap and parametrise
            ions = _pmd.load_file(ions_prep_filenames[0])
            for residue in ions.residues:
                if residue.name == "SOL":
                    residue.name = "WAT"

            ions.save(filebase + "_ions.pdb")
            _os.remove(ions_prep_filenames[0])
            ions_prep = _parametrise.parametriseAndLoadPmd(params, filebase + "_ions.pdb", "water")

            for filename in ions_prep_filenames:
                ions_prep.save(filename)

            return complex + readfunc(ions_prep_filenames)

        else:
            complex + readfunc(waters_prep_filenames)
