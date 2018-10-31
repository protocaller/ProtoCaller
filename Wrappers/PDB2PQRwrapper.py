'''TO DO:
1. check proteins with modified AA's
2. check numbering of ignored residues'''
import const as _const
from Wrappers import runexternal as _runexternal


#all options can be found here:
#https://apbs-pdb2pqr.readthedocs.io/en/latest/pdb2pqr/invoking.html
def PDB2PQRtransform(filename_input, filename_output=None, **kwargs):
    if filename_output == None:
        filename_output = filename_input[:-4] + "_pdb2pqr.pqr"
    default_kwargs = {'chain' : True, 'ff' : 'amber', 'ffout' : 'amber', 'verbose' : True, 'drop-water' : True}
    for key, val in kwargs.items():
        default_kwargs[key] = val

    runstring = "%s/pdb2pqr %s %s" % (_const.PDB2PQRDIR, filename_input, filename_output)
    for key, val in default_kwargs.items():
        if val is True:
            runstring += " --%s" % key
        elif val is not False:
            runstring += " --%s=%s" % (key, val)
    _runexternal.runExternal(runstring, procname="PDB2PQR")

    return filename_output