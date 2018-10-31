import os as _os
import platform as _platform
import sys as _sys

HOMEDIR = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.append("%s/shared" % HOMEDIR)
OS = _platform.system()
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