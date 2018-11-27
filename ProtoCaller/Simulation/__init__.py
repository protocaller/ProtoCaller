import copy as _copy
import glob as _glob
import os as _os

import ProtoCaller as _PC
import ProtoCaller.Protocol as _Protocol
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Utils.runexternal as _runexternal

class GMXSingleRun:
    def __init__(self, name, gro_file, top_file, lambda_index, work_dir=None, replica_top_files=None):
        if work_dir is None:
            work_dir = _os.getcwd()
        self.files = {}
        self.files["gro"] = _os.path.abspath(gro_file)
        self.files["top"] = _os.path.abspath(top_file)
        self._workdir = _fileio.Subdir("%s/%s/Lambda_%d" % (work_dir, name, lambda_index))
        self.lambda_index = lambda_index
        self.replica_tops = []
        if replica_top_files not in [None, []]:
            self.replica_tops = ["%s/%s" % (_os.getcwd(), file) for file in replica_top_files]

    def runSimulation(self, name, protocol, replex=None, n_cores=None, **replica_params):
        print("Running %s..." % name)
        with self._workdir:
            with _fileio.Subdir(name, overwrite=True):
                protocol.current_lambda = self.lambda_index
                self.files["mdp"] = _os.path.abspath(protocol.write(engine="GROMACS"))

                filebase = "%s_%d" % (name, self.lambda_index)
                grompp_args = {
                    "-f" : "mdp",
                    "-c" : "gro",
                    "-p" : "top",
                    "-t" : "cpt",
                }

                if replex is None or _PC.MPIEXE is None:
                    grompp_command = "%s grompp -maxwarn 10 -o input.tpr" % _PC.GROMACSEXE
                    for grompp_arg, filetype in grompp_args.items():
                        if filetype in self.files.keys():
                            grompp_command += " %s '%s'" % (grompp_arg, self.files[filetype])
                    _runexternal.runExternal(grompp_command, procname="gmx grompp")

                    mdrun_command = "%s mdrun -s input.tpr -deffnm %s" % (_PC.GROMACSEXE, filebase)
                else:
                    protocol_temp = _copy.copy(protocol)
                    n_replicas = len(self.replica_tops) + 1
                    for key, value in replica_params.items():
                        if isinstance(value, list) and len(value) != n_replicas:
                            raise ValueError("Size of parameter list not the same as the number of replicas.")
                        elif not isinstance(value, list):
                            replica_params[key] = [value] * n_replicas
                    for i in range(n_replicas):
                        grompp_command = "%s grompp -maxwarn 10 -o input%d.tpr" % (_PC.GROMACSEXE, i)
                        for grompp_arg, filetype in grompp_args.items():
                            if filetype in self.files.keys():
                                if filetype == "top":
                                    grompp_command += " %s '%s'" % (grompp_arg, ([self.files["top"]] + self.replica_tops)[i])
                                elif filetype == "mdp":
                                    for key, value in replica_params.items():
                                        protocol_temp.__setattr__(key, value[i])
                                    filename_mdp = _os.path.join(_os.getcwd(), protocol.write("GROMACS", "input%d" % i))
                                    grompp_command += " %s '%s'" % (grompp_arg, filename_mdp)
                                else:
                                    grompp_command += " %s '%s'" % (grompp_arg, self.files[filetype])
                        _runexternal.runExternal(grompp_command, procname="gmx grompp")
                    open("plumed.dat", "a").close()
                    n_cores = n_replicas if not n_cores else n_cores
                    # disable dynamic load balancing due to GROMACS bug
                    mdrun_command = "%s -np %d %s mdrun -dlb no -plumed plumed.dat -multi %d -s input.tpr -deffnm %s " \
                                    "-replex %d -hrex" % (_PC.MPIEXE, n_cores, _PC.GROMACSMPIEXE, n_replicas, filebase,
                                                          replex)

                _runexternal.runExternal(mdrun_command, procname="gmx mdrun")

                output_files = _glob.glob("%s/%s.*" % (_os.getcwd(), filebase))
                self.files = {"top" : self.files["top"],}
                for output_file in output_files:
                    ext = output_file.split(".")[-1].lower()
                    if ext == "tpr": ext = "gro"
                    self.files[ext] = output_file

class GMXSerialRuns:
    def __init__(self, name, gro_file, top_file, work_dir=None, replica_top_files=None, **lambda_dict):
        for key, arr in lambda_dict.items():
            if len(arr):
                self.lambda_size = len(arr)
            elif len(arr) != self.lambda_size:
                raise ValueError("Need lists of the same size")

        self.lambda_dict = lambda_dict
        self.gmx_run_list = [GMXSingleRun(name, gro_file, top_file, i, work_dir, replica_top_files)
                             for i in range(self.lambda_size)]
        self.protocols = []

    def addProtocol(self, name, use_preset=None, run_options=None, **kwargs):
        run_options = run_options if run_options is not None else {}
        self.protocols += [(name,
                            _Protocol.Protocol(use_preset=use_preset, **kwargs, **self.lambda_dict),
                            run_options)]

    def runSimulations(self):
        for i, gmx_run in enumerate(self.gmx_run_list):
            print("Running simulation %d..." % (i + 1))
            for name, protocol, kwargs in self.protocols:
                try:
                    gmx_run.runSimulation(name, protocol, **kwargs)
                except:
                    "An error occurred. Check the log. Continuing..."

class GMX_REST_FEP_Runs():
    def __init__(self, name, gro_file, top_files, work_dir=None, **lambda_dict):
        for key, arr in lambda_dict.items():
            if len(arr) and len(arr) != len(top_files):
                raise ValueError("Need lists of the same size")
        self.lambda_size = len(top_files)
        if work_dir is None:
            work_dir = _os.getcwd()

        self.files = [{} for i in range(self.lambda_size)]
        for f, top_file in zip(self.files, top_files):
            f["gro"] = _os.path.abspath(gro_file)
            f["top"] = _os.path.abspath(top_file)
        self._workdir = _fileio.Subdir("%s/%s/" % (work_dir, name))
        self.lambda_dict = lambda_dict
        self.protocols = []

    def runSimulation(self, name, parallel, use_preset=None, replex=500, n_cores=None, **protocol_params):
        protocol = _Protocol.Protocol(use_preset=use_preset, **protocol_params, **self.lambda_dict)

        print("Running %s..." % name)
        with self._workdir:
            with _fileio.Subdir(name, overwrite=True):
                for i in range(self.lambda_size):
                    filebase = "%s_%d" % (name, i)
                    protocol.init_lambda_state = i
                    self.files[i]["mdp"] = _os.path.abspath(protocol.write(engine="GROMACS", filebase=filebase))
                    grompp_args = {
                        "-f" : "mdp",
                        "-c" : "gro",
                        "-p" : "top",
                        "-t" : "cpt",
                    }

                    grompp_command = "%s grompp -maxwarn 10 -o %s.tpr" % (_PC.GROMACSEXE, filebase)
                    for grompp_arg, filetype in grompp_args.items():
                        if filetype in self.files[i].keys():
                            grompp_command += " %s '%s'" % (grompp_arg, self.files[i][filetype])
                    _runexternal.runExternal(grompp_command, procname="gmx grompp")

                    if not parallel:
                        mdrun_command = "{0} mdrun -s {1}.tpr -deffnm {1}".format(_PC.GROMACSEXE, filebase)
                        _runexternal.runExternal(mdrun_command, procname="gmx mdrun")

                        self.files[i] = {"top" : self.files[i]["top"],}
                        output_files = _glob.glob("%s.*" % filebase)
                        for output_file in output_files:
                            ext = output_file.split(".")[-1].lower()
                            if ext == "tpr": ext = "gro"
                            self.files[i][ext] = _os.path.abspath(output_file)

                if parallel:
                    open("plumed.dat", "a").close()
                    filebase = name + "_"
                    n_cores = self.lambda_size if not n_cores else n_cores
                    # disable dynamic load balancing due to GROMACS bug
                    mdrun_command = "{0} -np {1} {2} mdrun -dlb no -plumed plumed.dat -multi {3} -s {4}.tpr -deffnm {4} " \
                                    "-replex {5} -hrex".format(_PC.MPIEXE, n_cores, _PC.GROMACSMPIEXE, self.lambda_size,
                                                               filebase, replex)

                    _runexternal.runExternal(mdrun_command, procname="gmx mdrun")

def generateMBARInput(input_files, trajectory_files):
    if len(input_files) != len(trajectory_files):
        raise TypeError("Need to have a trajectory file for each input file.")

    xvg_files = []

    for i, tpr in enumerate(input_files):
        for j, trr in enumerate(trajectory_files):
            filebase = "Energy_%d_%d" % (i, j)
            mdrun_command = "%s mdrun -s %s -rerun %s -deffnm %s" % (_PC.GROMACSEXE, tpr, trr, filebase)
            try:
                _runexternal.runExternal(mdrun_command, procname="gmx mdrun")
            except:
                print("Mdrun failed for input file %d and trajectory %d. You might want to check your files.")
                continue
            xvg_files += [filebase + ".xvg"]

    return xvg_files