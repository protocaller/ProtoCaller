import os as _os
import platform as _platform
import re as _re
import subprocess as _subprocess
import sys as _sys

HOMEDIR = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.append("%s/shared" % HOMEDIR)
OS = _platform.system()
GROMACSEXE = None
if OS == "Linux":
    if not "GROMACSHOME" in _os.environ:
        if any(sys in _platform.linux_distribution()[0].lower() for sys in ["red hat", "centos"]):
            p = _subprocess.Popen(["module load gromacs\n which gmx"], stdout=_subprocess.PIPE, stderr=_subprocess.PIPE,
                                  shell=True)
            out, err = p.communicate()
            if err:
                p = _subprocess.Popen(["which gmx"], stdout=_subprocess.PIPE, stderr=_subprocess.PIPE, shell=True)
                out, err = p.communicate()
                if err:
                    p = _subprocess.Popen(["which gmx_mpi"], stdout=_subprocess.PIPE, stderr=_subprocess.PIPE, shell=True)
                    out, err = p.communicate()
                    if err:
                        print("GROMACS cannot be found on this system. This will affect some functionality.")
                    else:
                        _os.environ['GROMACSHOME'] = _re.match(r"([\S]*)/bin/gmx", out.decode()).group(1)
                        GROMACSEXE = _os.environ['GROMACSHOME'] + "/bin/gmx_mpi"
                else:
                    _os.environ['GROMACSHOME'] = _re.match(r"([\S]*)/bin/gmx", out.decode()).group(1)
            if out:
                GROMACSHOME = _re.match(r"([\S]*)/bin/gmx", out.decode()).group(1)
                _os.environ["GROMACSHOME"] = GROMACSHOME

if "GROMACSHOME" in _os.environ and GROMACSEXE is None:
    GROMACSEXE = _os.environ["GROMACSHOME"] + "/bin/gmx"
PDB2PQRDIR = "%s/shared/pdb2pqr-%s" % (HOMEDIR, OS)
assert _os.path.isdir(PDB2PQRDIR), "Cannot find pdb2pqr binary for the current operating system"
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

from . import Ensemble
from . import IO
from . import Parametrise
from . import Protocol
from . import Simulation
from . import Utils
