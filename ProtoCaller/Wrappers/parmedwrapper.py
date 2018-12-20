import os as _os

import parmed as _pmd

def openFilesAsParmed(filearr, fix_dihedrals=True):
    if len(filearr) == 1:
        mol = _pmd.load_file(filearr[0])
    else:
        try:
            mol = _pmd.load_file(filearr[0], xyz=filearr[1])
        except:
            try:
                mol = _pmd.load_file(filearr[1], xyz=filearr[0])
            except:
                raise OSError("There was an error while reading the input files.")

    if fix_dihedrals:
        for i, d_i in enumerate(mol.dihedrals):
            for j, d_j in enumerate(mol.dihedrals):
                if i < j:
                    continue
                if (d_i.atom1 == d_j.atom1 and d_i.atom2 == d_j.atom2 and
                    d_i.atom3 == d_j.atom3 and d_i.atom4 == d_j.atom4):
                    d_i.funct, d_j.funct = 9, 9

    return mol

def saveFilesFromParmed(system, filearr, overwrite=True):
    saved_files = []
    for filename in filearr:
        if _os.path.isfile(filename):
            if not overwrite:
                continue
            else:
                _os.remove(filename)
        system.save(filename)
        saved_files += [filename]
    return saved_files

def fixCharge(filearr):
    system = openFilesAsParmed(filearr)
    charges = [atom.charge for atom in system]
    total_charge, n_atoms = sum(charges), len(charges)
    d_charge = (round(total_charge) - total_charge) / n_atoms

    if not d_charge:
        return filearr

    for atom in system:
        atom.charge += d_charge

    return saveFilesFromParmed(system, filearr)