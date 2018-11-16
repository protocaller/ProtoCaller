'''TO DO:
1. check proteins with modified AA's
2. check numbering of ignored residues'''
import ProtoCaller as _PC
from ProtoCaller.Utils import runexternal as _runexternal

#all options can be found here:
#https://apbs-pdb2pqr.readthedocs.io/en/latest/pdb2pqr/invoking.html
def PDB2PQRtransform(filename_input, filename_output=None, **kwargs):
    if not _PC.PDB2PQREXE:
        raise OSError("Cannot find PDB2PQR binary for the current operating system.")
    if filename_output == None:
        filename_output = filename_input[:-4] + "_pdb2pqr.pqr"
    default_kwargs = {'chain' : True, 'ff' : 'amber', 'ffout' : 'amber', 'verbose' : True, 'drop-water' : True}
    for key, val in kwargs.items():
        default_kwargs[key] = val

    runstring = " ".join([_PC.PDB2PQREXE, filename_input, filename_output])
    for key, val in default_kwargs.items():
        if val is True:
            runstring += " --%s" % key
        elif val is not False:
            runstring += " --%s=%s" % (key, val)
    _runexternal.runExternal(runstring, procname="PDB2PQR")

    return filename_output
