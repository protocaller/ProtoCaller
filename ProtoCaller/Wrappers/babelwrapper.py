import os as _os

from ProtoCaller.Utils import runexternal as _runexternal


def babelTransform(input_filename, output_extension="mol2", pH=None, generate_3D_coords=False):
    """
    A light wrapper of some of OpenBabel's file conversion capabilities.

    Parameters
    ----------
    input_filename : str
        Name of the input file.
    output_extension : str
        Type of output extension.
    pH : float or None
        Hydrogenate the molecule for a certain pH. None means no hydrogenation.
    generate_3D_coords : bool
        Whether to generate 3D coordinates for the molecule.

    Returns
    -------
    filename : str
        The absolute path to the output file.
    """
    if input_filename is None:
        return None

    input_filebase, input_extension = _os.path.splitext(input_filename)[0], input_filename.split(".")[-1]
    output_filename = _os.path.abspath(_os.path.basename(input_filebase) + "." + output_extension)
    command = "obabel -i{0} \"{1}\" -o{2} -O \"{3}\"".format(input_extension, input_filename, output_extension,
                                                             output_filename)
    if pH is not None: command += " -p {}".format(pH)
    if generate_3D_coords: command += " --gen3d"
    _runexternal.runExternal(command, procname="OpenBabel")

    return _os.path.abspath(output_filename)
