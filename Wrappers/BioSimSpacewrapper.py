import BioSimSpace as _BSS

def rescaleSystemParams(system, scale):
    if scale == 1:
        return system
    mols_mod = []

    for mol in system.getMolecules():
        mol_edit = mol._sire_molecule.edit()
        for atom in mol_edit.atoms():
            idx = atom.index()
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
        mol._sire_molecule = mol_edit.commit()
        mols_mod += [mol]

    system_new = _BSS._SireWrappers._system.System(mols_mod)
    box = system._sire_system.property("space")
    system_new._sire_system.setProperty("space", box)

    return system_new