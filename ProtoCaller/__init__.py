import glob as _glob
import os as _os
import platform as _platform
import re as _re
import subprocess as _subprocess
import sys as _sys
import warnings as _warnings

HOMEDIR = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.append("%s/shared" % HOMEDIR)
OS = _platform.system()

def searchForPath(exe_name=None, var_name=None):
    if exe_name is not None:
        p = _subprocess.Popen(["which " + exe_name], stdout=_subprocess.PIPE, stderr=_subprocess.PIPE, shell=True)
        out, err = p.communicate()
        out = out.strip().decode()

    if exe_name is None or err:
        if var_name in _os.environ.keys():
            out = _os.environ[var_name]
        else:
            out = None
    return out

#setting MPI path variable
MPIEXE = searchForPath("mpirun", "MPIEXE")
if not MPIEXE:
    _warnings.warn("Cannot find MPI binary for the current operating system.")

#setting GROMACS path variables
GROMACSHOME = searchForPath(var_name="GROMACSHOME")
if not GROMACSHOME:
    GROMACSMPIEXE = searchForPath("gmx_mpi", "GROMACSMPIEXE")
    GROMACSEXE = searchForPath("gmx", "GROMACSEXE")
    if not GROMACSMPIEXE and not GROMACSEXE:
        _warnings.warn("GROMACS executable not found. Some functionality will be affected.")
    else:
        GROMACSMPIEXE = GROMACSMPIEXE if GROMACSMPIEXE else GROMACSEXE
        GROMACSEXE = GROMACSEXE if GROMACSEXE else GROMACSMPIEXE
        GROMACSHOME = _re.match(r"([\S]*)/bin/gmx", GROMACSEXE).group(1)
else:
    GROMACSEXE = _os.path.join(GROMACSHOME, "bin/gmx")
    GROMACSMPIEXE = _os.path.join(GROMACSHOME, "bin/gmx_mpi")
    if not _os.path.isfile(GROMACSEXE) and not _os.path.isfile(GROMACSMPIEXE):
        print("GROMACSHOME not a valid variable. Some functionality will be affected")
    else:
        GROMACSMPIEXE = GROMACSMPIEXE if _os.path.isfile(GROMACSMPIEXE) else GROMACSEXE
        GROMACSEXE = GROMACSEXE if _os.path.isfile(GROMACSEXE) else GROMACSMPIEXE

#setting PDB2PQR path variable
PDB2PQREXE = searchForPath("%s/shared/pdb2pqr-%s/pdb2pqr" % (HOMEDIR, OS), "PDB2PQREXE")
if not PDB2PQREXE:
    _warnings.warn("Cannot find PDB2PQR binary for the current operating system.")

_os.environ["SIRE_SILENT_PHONEHOME"] = "1"
#AMBERHOME = _os.environ["AMBERHOME"]

#identifiers for different chemical groups
ENGINES = ["AMBER", "GROMACS"]

def RESIDUETYPE(res):
    res = res.upper()
    if res in ["HOH", "WAT", "H2O", "SOL"]:
        return "water"
    elif res in ["CL", "BR", "F", "I"]:
        return "simple anion"
    elif res in ["SO4", "PO4", "CO3"]:
        return "complex anion"
    elif res in ["MG", "NA", "CA", "K"]:
        return "simple cation"
    elif res in ["FE", "CU", "NI", "CO", "MN", "CD"]:
        return "complex cation"
    elif res in ["ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", "HID", "HIE",
                 "HIS", "HIP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]:
        return "amino acid"
    elif res in ["ATP", "ADP", "GTP", "GDP", "FMN", "FAD", "HEM", "HEME", "NAD", "NAI", "NAP", "NDP"]:
        return "cofactor"
    else:
        return "ligand"

AMBERPROTEINFFS = ["ff14SB", "ff99SB", " ff15ipq", "fb15", "ff03.r1", "ff03ua"]
AMBERLIGANDFFS = ["gaff2", "gaff"]
AMBERWATERFFS = ["tip3p", "tip4p"]
AMBERFFS = AMBERPROTEINFFS + AMBERLIGANDFFS + AMBERWATERFFS
AMBEROLDFFS = ["ff99SB"]
AMBERDEFAULTPROTEINFF = "ff14SB"
AMBERDEFAULTLIGANDFF = "gaff"
AMBERDEFAULTWATERFF = "tip3p"

PROTEINFFS = AMBERPROTEINFFS
LIGANDFFS = AMBERLIGANDFFS
WATERFFS = AMBERWATERFFS

with _warnings.catch_warnings():
    _warnings.filterwarnings("ignore")
    from . import Ensemble
    from . import IO
    from . import Parametrise
    from . import Protocol
    from . import Simulation
    from . import Utils