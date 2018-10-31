import os as _os
import platform as _platform
import re as _re
import subprocess as _subprocess
import sys as _sys

HOMEDIR = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.append("%s/shared" % HOMEDIR)
OS = _platform.system()
if OS == "Linux":
    if any(sys in _platform.linux_distribution()[0].lower() for sys in ["red hat", "centos"]):
        p = _subprocess.Popen(["module load gromacs\n which gmx"], stdout=_subprocess.PIPE, stderr=_subprocess.PIPE, shell=True)
        out, err = p.communicate()
        if err:
            print("Module 'gromacs' not found on a Red Hat / CentOS system. Continuing without modules...")
        else:
            GROMACSHOME = _re.match(r"([\S]*)/bin/gmx", out.decode()).group(1)
            _os.environ["GROMACSHOME"] = GROMACSHOME
PDB2PQRDIR = "%s/shared/pdb2pqr-%s" % (HOMEDIR, OS)
assert _os.path.isdir(PDB2PQRDIR), "Cannot find pdb2pqr binary for the current operating system"
_os.environ["SIRE_SILENT_PHONEHOME"] = "1"
#WORKDIR = os.path.dirname(os.path.abspath(sys.argv[1]))
#_os.environ["AMBERHOME"] = HOMEDIR + "/amber16"
#AMBERHOME = _os.environ["AMBERHOME"]

#identifiers for different chemical groups
ENGINES = ["AMBER", "GROMACS"]
WATERALIAS = ["HOH", "WAT", "H2O", "SOL"]
ANIONLIST = ["CL", "BR", "F", "I"]
CATIONLIST = ["MG", "NA", "CA", "K", "FE", "CU", "NI", "CO", "MN", "CD"]
AMINOACIDLIST = ["ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", "HID", "HIE", "HIS", "HIP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
MISCLIST = ["SO4", "PO4", "CO3"]

AMBERPROTEINFFS = ["ff14SB", "ff99SB", " ff15ipq", "fb15", "ff03.r1", "ff03ua"]
AMBERLIGANDFFS = ["gaff2", "gaff"]
AMBERWATERFFS = ["tip3p", "tip4p"]
AMBEROLDFFS = ["ff99SB"]
AMBERDEFAULTPROTEINFF = "ff14SB"
AMBERDEFAULTLIGANDFF = "gaff"
AMBERDEFAULTWATERFF = "tip3p"

PROTEINFFS = AMBERPROTEINFFS
LIGANDFFS = AMBERLIGANDFFS
WATERFFS = AMBERWATERFFS

VERBOSE = True