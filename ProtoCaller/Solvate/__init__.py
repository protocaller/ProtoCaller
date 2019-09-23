import ProtoCaller as _PC

from collections.abc import Iterable as _Iterable
import fileinput as _fileinput
import os as _os
import tempfile as _tempfile

if _PC.BIOSIMSPACE:
    import BioSimSpace as _BSS
    import ProtoCaller.Wrappers.biosimspacewrapper as _BSSwrap
import parmed as _pmd

import ProtoCaller.Parametrise as _parametrise
import ProtoCaller.Protocol as _protocol
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Utils.runexternal as _runexternal
import ProtoCaller.Wrappers.parmedwrapper as _pmdwrap


def solvate(complex, params=None, box_length=8, shell=0, neutralise=True, ion_conc=0.154, centre=True,
            work_dir=None, filebase="complex"):
    """
    Uses gmx solvate and gmx genion to solvate the system and (optionally) add NaCl ions. This function preserves the
    crystal water molecules.

    Parameters
    ----------
    complex : BioSimSpace.System or parmed.structure.Structure
        The input unsolvated system.
    params : ProtoCaller.Parametrise.Params
        The input force field parameters.
    box_length : float, iterable
        Size of the box in nm.
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
    if params is None:
        params = _parametrise.Params()

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

    if not isinstance(box_length, _Iterable):
        box_length = 3 * [box_length]

    if work_dir is None:
        work_dir = _os.path.basename(_tempfile.TemporaryDirectory().name)
        temp = True
    else:
        temp = False

    with _fileio.Dir(dirname=work_dir, temp=temp):
        # centre
        if centre:
            complex, box_length, _ = centrefunc(complex, box_length)

        # solvate with gmx solvate and load unparametrised waters into memory
        files = _PC.IO.GROMACS.saveAsGromacs(filebase, complex)
        # reloading the complex fixes some problems with ParmEd
        if isinstance(complex, _pmd.Structure):
            complex = _pmdwrap.openFilesAsParmed(files)
        new_gro = filebase + "_solvated.gro"
        command = "{0} solvate -shell {1} -box {2[0]} {2[1]} {2[2]} -cp \"{3}\" -o \"{4}\"".format(
            _PC.GROMACSEXE, shell, box_length, files[1], new_gro)
        if params.water_points == 4:
            command += " -cs tip4p.gro"
        _runexternal.runExternal(command, procname="gmx solvate")
        complex_solvated = _pmd.load_file(new_gro, skip_bonds=True)
        waters = complex_solvated[":SOL"]

        # prepare waters for tleap and parametrise
        for residue in waters.residues:
            residue.name = "WAT"
        for i, atom in enumerate(waters.atoms):
            if "H" in atom.name:
                atom.name = atom.name[0] + atom.name[2]
            elif "O" in atom.name:
                atom.name = "O"
            else:
                atom.name = "EPW"

        # here we only parametrise a single water molecule in order to gain performance
        waters_prep_filenames = [filebase + "_waters.top", filebase + "_waters.gro"]
        waters[":1"].save(filebase + "_single_wat.pdb")
        water = _parametrise.parametriseAndLoadPmd(params, filebase + "_single_wat.pdb", "water")
        _pmdwrap.saveFilesFromParmed(waters, [waters_prep_filenames[1]], combine="all")
        _pmdwrap.saveFilesFromParmed(water, [waters_prep_filenames[0]])
        for line in _fileinput.input(waters_prep_filenames[0], inplace=True):
            line_new = line.split()
            if len(line_new) == 2 and line_new == ["WAT", "1"]:
                line = line.replace("1", "{}".format(len(waters.positions) // params.water_points))
            print(line, end="")
        waters_prep = _pmdwrap.openFilesAsParmed(waters_prep_filenames)

        waters_prep.box = _pmd.load_file(files[1], skip_bonds=True).box
        if any([neutralise, ion_conc, shell]):
            for residue in waters_prep.residues:
                residue.name = "SOL"
        _pmdwrap.saveFilesFromParmed(waters_prep, waters_prep_filenames)

        # add ions
        if any([neutralise, ion_conc, shell]):
            # write an MDP file
            _protocol.Protocol(use_preset="vacuum").write("GROMACS", "ions")

            # neutralise if needed
            charge = chargefunc(complex) if neutralise else 0
            volume = box_length[0] * box_length[1] * box_length[2] * 10 ** -24
            n_Na, n_Cl = [int(volume * 6.022 * 10 ** 23 * ion_conc) - abs(charge) // 2] * 2
            if neutralise:
                if charge < 0:
                    n_Na -= charge
                else:
                    n_Cl += charge

            # add ions with gmx genion
            ions_prep_filenames = [filebase + "_ions.top", filebase + "_ions.gro"]
            command = "{0} grompp -f ions.mdp -p {1} -c {2} -o \"{3}_solvated.tpr\"".format(
                _PC.GROMACSEXE, *waters_prep_filenames, filebase)
            _runexternal.runExternal(command, procname="gmx grompp")

            command = "{{ echo 2; }} | {0} genion -s \"{1}_solvated.tpr\" -o \"{2}\" -nn {3} -np {4}".format(
                _PC.GROMACSEXE, filebase, ions_prep_filenames[1], n_Cl, n_Na)
            _runexternal.runExternal(command, procname="gmx genion")

            # prepare waters for tleap
            ions = _pmd.load_file(ions_prep_filenames[1], skip_bonds=True)
            for residue in ions.residues:
                if residue.name == "SOL":
                    residue.name = "WAT"

            # here we only parametrise single ions to gain performance
            ion = ions[":WAT"][":1"] + ions[":NA"][":1"] + ions[":CL"][":1"]
            max_len = len(ion.residues)
            ion.save(filebase + "_single_ion.pdb")
            ion = _parametrise.parametriseAndLoadPmd(params, filebase + "_single_ion.pdb", "water")
            _pmdwrap.saveFilesFromParmed(ions, [ions_prep_filenames[1]], combine="all")
            _pmdwrap.saveFilesFromParmed(ion, [ions_prep_filenames[0]])
            mol_dict = {}
            for line in _fileinput.input(ions_prep_filenames[0], inplace=True):
                line_new = line.split()
                if len(line_new) == 2 and line_new[0] in ["WAT", "NA", "CL"] and line_new[1] == "1":
                    n_mols = len(ions[":{}".format(line_new[0])].positions)
                    if line_new[0] == "WAT":
                        n_mols //= params.water_points
                    line = line.replace("1", "{}".format(n_mols))
                    mol_dict[line_new[0]] = line
                    # preserve the order of water, sodium and chloride
                    if len(mol_dict) == max_len:
                        for x in ["WAT", "NA", "CL"]:
                            if x in mol_dict.keys():
                                print(mol_dict[x], end="")
                        mol_dict = {}
                else:
                    print(line, end="")

            return complex + readfunc(ions_prep_filenames)
        else:
            complex + readfunc(waters_prep_filenames)
