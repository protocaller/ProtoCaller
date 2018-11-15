import shutil as _shutil
import tempfile as _tempfile

import BioSimSpace as _BSS
import parmed as _pmd

import ProtoCaller.Utils.fileio as _fileio

def saveAsGromacs(filebase, system):
    tmp_dir = _tempfile.TemporaryDirectory()
    with _fileio.Dir(dirname=tmp_dir.name) as dir:
        _BSS.IO.saveMolecules(filebase, system, "Gro87,GroTop")
        _pmd.load_file(filebase + ".gro87").save(filebase + ".gro")
        _shutil.move(filebase + ".gro", "%s/%s" % (dir.initialdirname, filebase + ".gro"))
        _shutil.move(filebase + ".grotop", "%s/%s" % (dir.initialdirname, filebase + ".top"))
    return [filebase + ".gro", filebase + ".top"]