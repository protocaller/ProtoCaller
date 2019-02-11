import glob as _glob
import os as _os
import warnings as _warnings

import ProtoCaller as _PC
import ProtoCaller.Utils.runexternal as _runexternal

def amberWrapper(params, filename, molecule_type, id=None, *args, **kwargs):
    if id is None: id = molecule_type
    force_fields, files, param_files = [], [filename], []

    if molecule_type == "protein":
        force_fields = [params.protein_ff, params.water_ff]
    elif molecule_type in ["water", "simple_anion", "simple_cation"]:
        force_fields = [params.water_ff]
    elif molecule_type == "complex anion":
        _warnings.warn("AMBER parametrisation failed: polyatomic anions not supported")
        return
    elif molecule_type == "complex cation":
        _warnings.warn("AMBER parametrisation failed: transition metals not supported")
        return
    elif molecule_type == "cofactor":
        force_fields = [params.ligand_ff]
        param_files = _glob.glob("%s/shared/amber-parameters/cofactors/GDP.*" % _PC.HOMEDIR)
    elif molecule_type == "ligand":
        force_fields = [params.ligand_ff]
        files = [runAntechamber(params.ligand_ff, filename)]
        param_files = [runParmchk(params.ligand_ff, files[0])]
    else:
        raise ValueError("Value %s for molecule_type not supported " % molecule_type)

    return runTleap(force_fields=force_fields, files=files, param_files=param_files, id=id, *args, **kwargs)

def runAntechamber(force_field, file, output_ext="mol2"):
    input_base, input_ext = _os.path.splitext(file)[0], file.split(".")[-1]
    output_name = "%s_antechamber.%s" % (input_base, output_ext)

    commandstr = "antechamber -i '%s' -fi %s -o '%s' -fo %s -c bcc -at %s -s 2" % (
        file, input_ext, output_name, output_ext, force_field)
    _runexternal.runExternal(commandstr, procname="antechamber")

    return _os.path.abspath(output_name)

def runParmchk(force_field, file):
    input_base, input_ext = _os.path.splitext(file)[0], file.split(".")[-1]
    commandstr = "parmchk2 -i '{0}.{1}' -f {1} -o '{0}.frcmod' -s %s".format(input_base, input_ext, force_field)
    _runexternal.runExternal(commandstr, procname="parmchk2")

    return _os.path.abspath(input_base + ".frcmod")

def runTleap(force_fields=None, files=None, param_files=None, id=None, disulfide_bonds=None):
    if id is None: id = "molecule"
    if disulfide_bonds is None: disulfide_bonds = []

    filename_tleap = "tleap_script_%s.in" % id
    filenames = [id + ".prmtop", id + ".inpcrd"]

    with open(filename_tleap, "w+") as out:
        for force_field in force_fields:
            out.write("source \"%s\"\n" % returnFFPath(force_field))
        if param_files is not None:
            for param_file in param_files:
                ext = param_file.split(".")[-1].lower()
                if ext in ["frcmod", "frcfld"]:
                    command = "loadamberparams"
                elif ext in ["off", "lib"]:
                    command = "loadoff"
                elif ext == "prep":
                    command = "loadamberprep"
                else:
                    _warnings.warn("%s is not a valid parameter file. Skipping value..." % param_file)
                    continue
                out.write("%s \"%s\"\n" % (command, param_file))
        for file in files:
            ext = file.split(".")[-1]
            if ext == "pqr": ext = "pdb"
            out.write("MOL = load%s \"%s\"\n" % (ext, file))
        for disulfide_bond in disulfide_bonds:
            out.write("bond MOL.%d.SG MOL.%d.SG\n" % (disulfide_bond[0].resSeq, disulfide_bond[1].resSeq))
        out.write("check MOL\n")
        out.write("saveAmberParm MOL {0} {1}\n".format(*filenames))
        out.write("quit\n")

    _runexternal.runExternal("tleap -f %s" % filename_tleap, procname="tleap")
    return [_os.path.abspath(f) for f in filenames]

def returnFFPath(name):
    if(name in _PC.AMBEROLDFFS):
        return "oldff/leaprc." + name
    elif name in _PC.AMBERPROTEINFFS:
        return "leaprc.protein.%s" % name
    elif name in _PC.AMBERWATERFFS:
        return "leaprc.water.%s" % name
    else:
        return "leaprc." + name
