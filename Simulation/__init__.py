import const as _const

import glob as _glob
import os as _os

import Protocol as _Protocol
import subdir as _subdir
import Wrappers.runexternal as _runexternal

class GMXSingleRun:
    def __init__(self, name, gro_file, top_file, lambda_index, work_dir=None):
        if work_dir is None:
            work_dir = _os.getcwd()
        self.files = {}
        self.files["gro"] = "%s/%s" % (_os.getcwd(), gro_file)
        self.files["top"] = "%s/%s" % (_os.getcwd(), top_file)
        self._subdir = _subdir.Subdir("%s/%s/Lambda_%d" % (work_dir, name, lambda_index))
        self.lambda_index = lambda_index

    def runSimulation(self, name, protocol):
        with self._subdir:
            with _subdir.Subdir(name, overwrite=True):
                protocol.current_lambda = self.lambda_index
                self.files["mdp"] = "%s/%s" % (_os.getcwd(), protocol.write(engine="GROMACS"))

                filebase = "%s_%d" % (name, self.lambda_index)
                grompp_args = {
                    "-f" : "mdp",
                    "-c" : "gro",
                    "-p" : "top",
                    "-t" : "cpt",
                }
                grompp_command = "%s grompp -maxwarn 10 -o input.tpr" % _const.GROMACSEXE
                for grompp_arg, filetype in grompp_args.items():
                    if filetype in self.files.keys():
                        grompp_command += " %s '%s'" % (grompp_arg, self.files[filetype])
                _runexternal.runExternal(grompp_command, procname="gmx grompp")

                mdrun_command = "%s mdrun -s input.tpr -deffnm %s" % (_const.GROMACSEXE, filebase)
                _runexternal.runExternal(mdrun_command, procname="gmx mdrun")

                output_files = _glob.glob("%s/%s.*" % (_os.getcwd(), filebase))
                self.files = {"top" : self.files["top"],}
                for output_file in output_files:
                    ext = output_file.split(".")[-1].lower()
                    if ext == "tpr": ext = "gro"
                    self.files[ext] = output_file

class GMXSerialRuns:
    def __init__(self, gro_file, top_file, name, work_dir=None, **lambda_dict):
        for key, arr in lambda_dict.items():
            if len(arr):
                self.lambda_size = len(arr)
            elif len(arr) != self.lambda_size:
                raise ValueError("Need lists of the same size")

        self.lambda_dict = lambda_dict
        self.gmx_run_list = [GMXSingleRun(gro_file, top_file, name, i, work_dir) for i in range(self.lambda_size)]
        self.protocols = []

    def addProtocol(self, name, use_preset=None, **kwargs):
        self.protocols += [(name, _Protocol.Protocol(use_preset=use_preset, **kwargs, **self.lambda_dict))]

    def runSimulations(self):
        for gmx_run in self.gmx_run_list:
            for name, protocol in self.protocols:
                try:
                    gmx_run.runSimulation(name, protocol)
                except:
                    "An error occurred. Check the log. Continuing..."