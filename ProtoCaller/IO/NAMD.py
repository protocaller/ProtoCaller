import ProtoCaller as _PC

import os as _os

if _PC.BIOSIMSPACE:
    import BioSimSpace as _BSS
import parmed as _pmd


def saveAsNamd(filebase, system):
    """
        Saves the input object to a PDB and a PSF file.

        Parameters
        ----------
        filebase : str
            Base of the output filenames.
        system : BioSimSpace.System or parmed.structure.Structure
            Input object to be saved.

        Returns
        -------
        files : [str, str]
            A list containing the absolute paths of the PSF and PDB files.
        """
    if "BioSimSpace" in str(type(system)):
        _BSS.IO.saveMolecules(filebase, system, "PDB,PSF")
    elif isinstance(system, (_pmd.Structure, _PC.Morph.Morph)):
        filenames = [filebase + ".psf", filebase + ".pdb"]
        for filename in filenames:
            if _os.path.exists(filename):
                _os.remove(filename)
            system.save(filename)
    else:
        raise TypeError(
            "Passed object to save as NAMD not recognised. Please pass BioSimSpace or ParmEd "
            "objects only.")
    return [_os.path.abspath(filebase + ".psf"),
            _os.path.abspath(filebase + ".pdb")]