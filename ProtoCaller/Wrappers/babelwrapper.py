from ProtoCaller.Utils import runexternal as _runexternal

def babelTransform(input_filename, output_extension="mol2", pH=7.0):
    if input_filename is None:
        return None

    input_filebase, input_extension = input_filename.split(".")
    command = "obabel -i {0} {1}.{0} -o {2} -O {1}.{2} -p {3}".format(input_extension, input_filebase,
                                                                      output_extension, pH)
    _runexternal.runExternal(command, procname="OpenBabel")
    output_filename = input_filebase + "." + output_extension

    return output_filename
