import os as _os

from ProtoCaller.Utils import runexternal as _runexternal


def babelTransform(input_filename, output_extension="mol2", pH=None):
    if input_filename is None:
        return None

    input_filebase, input_extension = _os.path.splitext(input_filename)[0], input_filename.split(".")[-1]
    output_filename = _os.path.abspath(_os.path.basename(input_filebase) + "." + output_extension)
    command = "obabel -i {0} \"{1}.{0}\" -o \"{2}\" -O {3}".format(input_extension, input_filebase, output_extension,
                                                                   output_filename)
    if pH is not None: command += " -p {}".format(pH)
    _runexternal.runExternal(command, procname="OpenBabel")

    return _os.path.abspath(output_filename)
