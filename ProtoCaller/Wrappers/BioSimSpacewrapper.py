import ProtoCaller as _PC
if not _PC.BIOSIMSPACE:
    raise ImportError("BioSimSpace module cannot be imported")

import BioSimSpace as _BSS
import Sire.MM as _SireMM
import Sire.Mol as _SireMol

def rescaleSystemParams(system, scale, includelist=None, excludelist=None, neutralise=True):
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