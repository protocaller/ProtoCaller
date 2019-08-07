import os as _os

import numpy as _np
import parmed as _pmd


def openFilesAsParmed(filelist, **kwargs):
    """
    Opens an input file list and returns a ParmEd Structure.

    Parameters
    ----------
    filename : str
        The name of the input file.
    kwargs
        Keyword arguments to be supplied to parmed.load_file.

    Returns
    -------
    mol : parmed.structure.Structure
        The file loaded as a ParmEd Structure object.
    """
    if len(filelist) == 1:
        mol = _pmd.load_file(filelist[0], **kwargs)
    else:
        try:
            mol = _pmd.load_file(filelist[0], xyz=filelist[1], **kwargs)
        except:
            try:
                mol = _pmd.load_file(filelist[1], xyz=filelist[0], **kwargs)
            except:
                raise OSError("There was an error while reading the input files.")

    return mol


def saveFilesFromParmed(system, filelist, overwrite=True, **kwargs):
    """
    Saves a ParmEd Structure object to a file.

    Parameters
    ----------
    system : parmed.structure.Structure
        The input system.
    filelist : [str]
        The names of the output files.
    overwrite : bool
        Whether to overwrite existing files.
    kwargs :
        Additional arguments to be pased on to save()

    Returns
    -------
    saved_files : [str]
        The absolute path to the written files.
    """
    saved_files = []
    for filename in filelist:
        if _os.path.isfile(filename):
            if not overwrite:
                continue
            else:
                _os.remove(filename)
        system.save(filename, **kwargs)
        saved_files += [filename]
    return saved_files


def fixCharge(filelist):
    """
    Makes sure that the net charge is always an integer up to a float precision.

    Parameters
    ----------
    filelist: [str]
        An array of input files.

    Returns
    -------
    filelist: [str]
        An array of fixed output files.

    """
    system = openFilesAsParmed(filelist)
    charges = [atom.charge for atom in system]
    total_charge, n_atoms = sum(charges), len(charges)
    d_charge = (round(total_charge) - total_charge) / n_atoms

    if not d_charge:
        return filelist

    for atom in system:
        atom.charge += d_charge

    return saveFilesFromParmed(system, filelist)


def centre(system, box_length):
    """
    Centres the system given a box length (cubic shape is assumed).

    Parameters
    ----------
    system : parmed.structure.Structure
        The input system to be centred.
    box_length : float
        Length of the cubic box.

    Returns
    -------
    system : parmed.structure.Structure
        The centred ParmEd system.
    box_length : float
        The new box length. Only different from the input value if the box is too small for the system.
    translation_vec: numpy.array
        The centering translation vector.
    """
    coords = system.coordinates
    min_coords, max_coords = _np.amin(coords, axis=0), _np.amax(coords, axis=0)
    centre = (min_coords + max_coords) / 2
    difference = max_coords - min_coords
    if any([x / 10 > box_length for x in difference]):
        box_length = int(max(difference / 10)) + 1
        print("Insufficient input box size. Changing to a box length of %d nm..." % box_length)
    translation_vec = -centre + _np.asarray(3 * [5 * box_length])
    for atom in system.atoms:
        atom.xx += translation_vec[0]
        atom.xy += translation_vec[1]
        atom.xz += translation_vec[2]
    system = resize(system, box_length)

    return system, box_length, translation_vec


def resize(system, box_length):
    """
    Changes the box size of the system or adds one if there is no box. Only valid for cubic boxes.

    Parameters
    ----------
    system : parmed.structure.Structure
        The input system to be resized.
    box_length : float
        Length of the cubic box.

    Returns
    -------
    system : parmed.structure.Structure
        The resized system.
    """
    system.box = 3 * [10 * box_length] + 3 * [90]
    system.box[0] = system.box[1] = system.box[2] = 10 * box_length
    return system
