import ProtoCaller as _PC
if not _PC.BIOSIMSPACE:
    raise ImportError("BioSimSpace module cannot be imported")

import BioSimSpace as _BSS
import Sire.MM as _SireMM
import Sire.Maths as _SireMaths
import Sire.Mol as _SireMol
import Sire.Vol as _SireVol


def centre(system, box_length):
    """
    Centres the system given a box length (cubic shape is assumed).

    Parameters
    ----------
    system : BioSimSpace.System
        The input system to be centred.
    box_length : float
        Length of the cubic box.

    Returns
    -------
    system : BioSimspace.System
        The centred BioSimSpace system.
    box_length : float
        The new box length. Only different from the input value if the box is too small for the system.
    translation_vec: Sire.Maths.vector
        The centering translation vector.
    """
    box = system._getAABox()
    min_coords, max_coords = box.minCoords(), box.maxCoords()
    centre = (min_coords + max_coords) / 2
    difference = max_coords - min_coords
    if any([x / 10 > box_length for x in difference]):
        box_length = int(max(difference / 10)) + 1
        print("Insufficient input box size. Changing to a box length of %d nm..." % box_length)
    translation_vec = -centre + _SireMaths.Vector(3 * (box_length * 5,))
    system.translate(tuple(translation_vec))
    system = resize(system, box_length)

    return system, box_length, translation_vec


def resize(system, box_length):
    """
    Changes the box size of the system or adds one if there is no box. Only valid for cubic boxes.

    Parameters
    ----------
    system : BioSimSpace.System
        The input system to be resized.
    box_length : float
        Length of the cubic box.

    Returns
    -------
    system : BioSimSpace.System
        The resized system.
    """
    system._sire_system.setProperty("space", _SireVol.PeriodicBox(_SireMaths.Vector(3 * (box_length * 10,))))
    return system


def rescaleSystemParams(system, scale, includelist=None, excludelist=None, neutralise=True):
    """
    Rescales charge, Lennard-Jones and dihedral parameters for REST(2).

    Parameters
    ----------
    system : BioSimSpace.System
        Input system to be rescaled.
    scale : float
        Multiplication factor for the charge, Lennard-Jones and dihedral parameters. 0 < scale < 1.
    includelist : [str]
        Molecule names to be included in the rescaling. Default: None.
    excludelist : [str]
        Molecule names to be excluded in the rescaling. Default: ["WAT"].
    neutralise : bool
        Whether to rescale a minimum number of Na+ or Cl- ions so that the rescaled part of the system is neutralised.

    Returns
    -------
    system_new : BioSimSpace.System
        The rescaled system.
    """
    if excludelist is None and includelist is None:
        excludelist = ["WAT"]
    elif excludelist is not None and includelist is not None:
        raise ValueError("Only one of includelist and excludelist can be set")

    if neutralise:
        if excludelist is not None:
            total_charge = round(sum([mol.charge().magnitude() for mol in system.getMolecules()
                                      if mol._sire_molecule.name().value() not in excludelist]))
        elif includelist is not None:
            total_charge = round(sum([mol.charge().magnitude() for mol in system.getMolecules()
                                      if mol._sire_molecule.name().value() in includelist]))

    mols_mod = []

    for mol in system.getMolecules():
        if neutralise and mol._sire_molecule.name().value().lower() == "na" and total_charge < 0:
            total_charge += 1
        elif neutralise and mol._sire_molecule.name().value().lower() == "cl" and total_charge > 0:
            total_charge -= 1
        elif excludelist is not None:
            if mol._sire_molecule.name().value() in excludelist:
                mols_mod += [mol]
                continue
        elif includelist is not None:
            if mol._sire_molecule.name().value() not in includelist:
                mols_mod += [mol]
                continue

        mol_edit = mol._sire_molecule.edit()

        if scale != 1:
            for prop in ["dihedral", "dihedral0", "dihedral1", "improper", "improper0", "improper1"]:
                try:
                    potentials = mol_edit.property(prop).potentials()
                    prop_new = _SireMM.FourAtomFunctions(mol_edit.info())
                    for potential in potentials:
                        prop_new.set(potential.atom0(), potential.atom1(), potential.atom2(), potential.atom3(),
                                     potential.function() * scale)
                    mol_edit = mol_edit.setProperty(prop, prop_new).molecule()
                except:
                    continue

        for atom in mol_edit.atoms():
            idx = atom.index()

            for prop in ["atomtype", "atomtype0", "atomtype1"]:
                try:
                    atomtype_new = atom.property(prop) + "_"
                    mol_edit = mol_edit.atom(idx).setProperty(prop, atomtype_new).molecule()
                except:
                    continue

            if scale != 1:
                for prop in ["LJ", "LJ0", "LJ1"]:
                    try:
                        LJ = atom.property(prop)
                        LJ = LJ.fromSigmaAndEpsilon(LJ.sigma(), LJ.epsilon() * scale)
                        mol_edit = mol_edit.atom(idx).setProperty(prop, LJ).molecule()
                    except:
                        continue

                for prop in ["charge", "charge0", "charge1"]:
                    try:
                        charge = atom.property(prop)
                        charge = charge * scale
                        mol_edit = mol_edit.atom(idx).setProperty(prop, charge).molecule()
                    except:
                        continue

        if mol._sire_molecule.name().value().lower() in ["na", "cl"]:
            resname = mol_edit.residue().name()
            resname_new = _SireMol.ResName(resname.value() + "_")
            mol_edit = mol_edit.residue(resname).rename(resname_new)

        mol._sire_molecule = mol_edit.commit()
        mols_mod += [mol]

    system_new = _BSS._SireWrappers._system.System(mols_mod)
    box = system._sire_system.property("space")
    system_new._sire_system.setProperty("space", box)

    return system_new
