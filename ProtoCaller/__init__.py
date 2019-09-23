import os as _os
import platform as _platform
import re as _re
import subprocess as _subprocess
import sys as _sys
import warnings as _warnings

HOMEDIR = _os.path.dirname(_os.path.abspath(__file__))
TESTDIR = _os.path.abspath(_os.path.dirname(HOMEDIR) + "/test")
SHAREDDIR = _os.path.abspath(HOMEDIR + "/shared")

_sys.path.append(SHAREDDIR)
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
        _warnings.warn("GROMACSHOME not a valid variable. Some functionality will be affected")
    else:
        GROMACSMPIEXE = GROMACSMPIEXE if _os.path.isfile(GROMACSMPIEXE) else GROMACSEXE
        GROMACSEXE = GROMACSEXE if _os.path.isfile(GROMACSEXE) else GROMACSMPIEXE


_os.environ["SIRE_SILENT_PHONEHOME"] = "1"
#AMBERHOME = _os.environ["AMBERHOME"]

ENGINES = ["GROMACS"]

#identifiers for different chemical groups
def RESIDUETYPE(res):
    namesdict = {
        "water": WATERNAMES,
        "simple_anion": SIMPLEANIONNAMES,
        "complex_anion": COMPLEXANIONNAMES,
        "simple_cation": SIMPLECATIONNAMES,
        "complex_cation": COMPLEXCATIONNAMES,
        "amino_acid": AMINOACIDNAMES,
        "amino_acid_modified": MODIFIEDAMINOACIDNAMES,
        "cofactor": COFACTORNAMES,
    }

    res = res.upper().strip()
    for key, value in namesdict.items():
        if res in value:
            return key
    return "ligand"

WATERNAMES = ["HOH", "WAT", "H2O", "SOL"]
SIMPLEANIONNAMES = ["CL", "BR", "F", "I"]
COMPLEXANIONNAMES = ["SO4", "PO4", "CO3"]
SIMPLECATIONNAMES = ["MG", "NA", "CA", "K"]
COMPLEXCATIONNAMES = ["FE", "CU", "NI", "CO", "MN", "CD"]
AMINOACIDNAMES = ["ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", "HID", "HIE",
                 "HIS", "HIP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
MODIFIEDAMINOACIDNAMES = ["MSE"]
COFACTORNAMES = ["ATP", "ADP", "GTP", "GDP", "FMN", "FAD", "HEM", "HEME", "NAD", "NAI", "NAP", "NDP"]

AMBERPROTEINFFS = ["ff14SB", "ff14SBonlysc", "ff99SB", "ff15ipq", "ff15ipq-vac", "fb15", "ff03.r1"]
AMBERLIGANDFFS = ["gaff2", "gaff"]
AMBERWATERFFS = ["tip3p", "tip4pew", "opc", "spce", "spceb", "fb3", "fb4"]
AMBERFFS = AMBERPROTEINFFS + AMBERLIGANDFFS + AMBERWATERFFS
AMBEROLDFFS = ["ff99SB"]
AMBERDEFAULTPROTEINFF = "ff14SB"
AMBERDEFAULTLIGANDFF = "gaff"
AMBERDEFAULTWATERFF = "tip3p"

PROTEINFFS = AMBERPROTEINFFS
LIGANDFFS = AMBERLIGANDFFS
WATERFFS = AMBERWATERFFS

try:
    import BioSimSpace as _BSS
    BIOSIMSPACE = True
except:
    BIOSIMSPACE = False

try:
    import Sire as _Sire
    SIRE = True
except:
    SIRE = False

try:
    import modeller as modeller
    MODELLER = True
except:
    MODELLER = False

with _warnings.catch_warnings():
    _warnings.filterwarnings("ignore")
    from . import Ensemble
    from . import IO
    from . import Morph
    from . import Parametrise
    from . import Protocol
    from . import Simulation
    from . import Solvate
    from . import Utils
    from . import Wrappers