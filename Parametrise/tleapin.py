import const as _const

from Wrappers import runexternal as _runexternal

"""
Note that all scripts must return the filenames they create as a list of strings.
"""

def RunTleapScript(filename):
    _runexternal.runExternal("tleap -f %s" % filename, procname="tleap")

def RunAntechamber(params, input_filebase, input_extension, output_filebase=None, output_extension="mol2"):
    if output_filebase is None:
        output_filebase = input_filebase + "_antechamber"
    input_filename = input_filebase + "." + input_extension
    output_filename = output_filebase + "." + output_extension
    commandstr = "antechamber -i %s -fi %s -o %s -fo %s -c -at %s bcc -s 2" % (
        input_filename, input_extension, output_filename, output_extension, params.ligand_ff)

    _runexternal.runExternal(commandstr, procname="antechamber")

def RunParmchk(params, input_filebase, input_extension):
    commandstr = "parmchk2 -i {0}.{1} -f {1} -o {0}.frcmod -s %s".format(input_filebase, input_extension, params.ligand_ff)
    _runexternal.runExternal(commandstr, procname="parmchk2")

#generates tleap input for amber protein parametrisation
def WriteProteinIN(params, input_filebase="protein", input_extension="pdb", output_filename="tleap_protein.in"):
    with open(output_filename, "w") as file:
        file.write("source %s\n" % returnFFPath(params.protein_ff, "protein"))
        file.write("source %s\n" % returnFFPath(params.water_ff, "water"))
        file.write("PRO = loadpdb %s.%s\n" % (input_filebase, input_extension))
        file.write("saveAmberParm PRO {0}.prmtop {0}.inpcrd\n".format(input_filebase))
        file.write("quit\n")
    return [input_filebase + ".prmtop", input_filebase + ".inpcrd"]

#generates tleap input for amber ligand parametrisation
def WriteLigandIN(params, input_filebase="ligand", input_extension="mol2", output_filename="tleap_ligand.in"):
    with open(output_filename, "w") as file:
        file.write("source %s\n" % returnFFPath(params.protein_ff, "protein"))
        file.write("source %s\n" % returnFFPath(params.ligand_ff, "ligand"))
        file.write("LIG = load{0} {1}.{0}\n".format(input_extension, input_filebase))
        file.write("check LIG\n")
        file.write("loadamberparams {0}.frcmod\n".format(input_filebase))
        file.write("saveoff LIG LIG.lib\n")
        file.write("saveamberparm LIG {0}.prmtop {0}.inpcrd\n".format(input_filebase))
        file.write("quit\n")
    return [input_filebase + ".prmtop", input_filebase + ".inpcrd"]

#generates tleap input for amber solvent parametrisation
def WriteSolventIN(params, input_filebase="solvent", input_extension="pdb", output_filename="tleap_solvent.in"):
    with open(output_filename, "w") as file:
        file.write("source %s\n" % returnFFPath(params.water_ff, "water"))
        file.write("SOL = load{0} {1}.{0}\n".format(input_extension, input_filebase))
        file.write("saveamberparm SOL {0}.prmtop {0}.inpcrd\n".format(input_filebase))
        file.write("quit\n")
    return [input_filebase + ".prmtop", input_filebase + ".inpcrd"]

def WriteAnionIN(params, input_filebase="anion", input_extension="pdb", output_filename="tleap_anion.in"):
    return WriteSolventIN(params=params, input_filebase=input_filebase, input_extension=input_extension,
                          output_filename=output_filename)

def WriteCationIN(params, input_filebase="cation", input_extension="pdb", output_filename="tleap_cation.in"):
    return WriteSolventIN(params=params, input_filebase=input_filebase, input_extension=input_extension,
                          output_filename=output_filename)

def returnFFPath(name, type):
    if(name in _const.AMBEROLDFFS):
        return "oldff/leaprc." + name
    elif(type in ["protein", "water"]):
        return "leaprc.%s.%s" % (type, name)
    else:
        return "leaprc." + name