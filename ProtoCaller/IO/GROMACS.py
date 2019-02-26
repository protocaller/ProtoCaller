import ProtoCaller as _PC

import os as _os
import shutil as _shutil
import tempfile as _tempfile

if _PC.BIOSIMSPACE:
    import BioSimSpace as _BSS
import parmed as _pmd

import ProtoCaller.Utils.fileio as _fileio

def saveAsGromacs(filebase, system):
    """
    Saves the input object to a GRO and a TOP file.

    Parameters
    ----------
    filebase : str
        Base of the output filenames.
    system : BioSimSpace.System or parmed.structure.Structure
        Input object to be saved.

    Returns
    -------
    files : [str, str]
        A list containing the absolute paths of the TOP and GRO files.
    """
    name = _os.path.basename(_tempfile.TemporaryDirectory().name)
    with _fileio.Dir(dirname=name, temp=True) as dir:
        if "BioSimSpace" in str(type(system)):
            _BSS.IO.saveMolecules(filebase, system, "Gro87,GroTop")
            _pmd.load_file(filebase + ".gro87").save(filebase + ".gro", combine="all")
            _shutil.move(filebase + ".gro", _os.path.join(dir.workdirname, filebase + ".gro"))
            _shutil.move(filebase + ".grotop", _os.path.join(dir.workdirname, filebase + ".top"))
        elif isinstance(system, (_pmd.Structure, _PC.Morph.Morph)):
            system.save(filebase + ".gro")
            system.save(filebase + ".top")
            _shutil.move(filebase + ".gro", _os.path.join(dir.workdirname, filebase + ".gro"))
            _shutil.move(filebase + ".top", _os.path.join(dir.workdirname, filebase + ".top"))
        else:
            raise TypeError("Passed object to save as GROMACS not recognised. Please pass BioSimSpace or ParmEd "
                            "objects only.")
    return [_os.path.abspath(filebase + ".gro"), _os.path.abspath(filebase + ".top")]