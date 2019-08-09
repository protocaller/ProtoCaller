from collections.abc import Iterable as _Iterable
import os as _os
import warnings as _warnings

import numpy as _np
import parmed as _pmd


def openFilesAsParmed(filelist, **kwargs):
    """
    Opens an input file list and returns a ParmEd Structure.

    Parameters
    ----------
    filename : str, list
        The name of the input file.
    kwargs
        Keyword arguments to be supplied to parmed.load_file.

    Returns
    -------
    mol : parmed.structure.Structure
        The file loaded as a ParmEd Structure object.
    """
    if isinstance(filelist, str):
        filelist = [filelist]
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
    filelist : str, [str]
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
    if isinstance(filelist, str):
        filelist = [filelist]
        single_file = True
    else:
        single_file = False

    for filename in filelist:
        if _os.path.isfile(filename):
            if not overwrite:
                continue
            else:
                _os.remove(filename)
        system.save(filename, **kwargs)
        saved_files += [filename]

    return saved_files[0] if single_file else saved_files


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
    Centres the system given a box length.

    Parameters
    ----------
    system : parmed.structure.Structure
        The input system to be centred.
    box_length : float, iterable
        Length of the box.

    Returns
    -------
    system : parmed.structure.Structure
        The centred ParmEd system.
    box_length : tuple
        The new box length. Only different from the input value if the box is too small for the system.
    translation_vec: numpy.ndarray
        The centering translation vector.
    """
    if not isinstance(box_length, _Iterable):
        cubic = True
        box_length = 3 * [box_length]
    else:
        cubic = (len(set(box_length)) == 1)

    coords = system.coordinates
    min_coords, max_coords = _np.amin(coords, axis=0), _np.amax(coords, axis=0)
    centre = (min_coords + max_coords) / 2
    difference = max_coords - min_coords

    box_length_new = []
    for x, y in zip(difference, box_length):
        box_length_new += [y] if x < 10 * y else [x / 10 + 1]
    if cubic:
        box_length_new = 3 * [max(box_length_new)]
    if box_length_new != box_length:
        _warnings.warn("Insufficient input box size. Changing to a box size of ({:.3f}, {:.3f}, {:.3f}) nm...".format(
            *box_length_new))

    translation_vec = -centre + _np.asarray([5 * x for x in box_length_new])
    for atom in system.atoms:
        atom.xx += translation_vec[0]
        atom.xy += translation_vec[1]
        atom.xz += translation_vec[2]
    system = resize(system, box_length_new)

    return system, tuple(box_length_new), translation_vec


def resize(system, box_length):
    """
    Changes the box size of the system or adds one if there is no box.

    Parameters
    ----------
    system : parmed.structure.Structure
        The input system to be resized.
    box_length : float, iterable
        Length of the box.

    Returns
    -------
    system : parmed.structure.Structure
        The resized system.
    """
    if not isinstance(box_length, _Iterable):
        box_length = 3 * [box_length]
    box_length = [10 * x for x in box_length]
    system.box = box_length + 3 * [90]
    return system
