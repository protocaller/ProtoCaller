import copy as _copy
import glob as _glob
import os as _os

import MDAnalysis as _MDAnalysis
import numpy as _np
import pymbar as _pymbar

import ProtoCaller as _PC
import ProtoCaller.Protocol as _Protocol
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Utils.runexternal as _runexternal


class GMXSingleRun:
    """
    A wrapper for a single GROMACS run.

    Parameters
    ----------
    name : str
        The name of the run.
    gro_file : str
        Path to the input GRO file.
    top_file : str
        Path to the input TOP file.
    lambda_index : int
        Initialises lambda_index.
    work_dir : str, optional
        Initialises the run in a custom directory.
    replica_top_files : [str], optional
        Initialises replica_tops.

    Attributes
    ----------
    files : dict
        A dictionary with extensions as keys and absolute paths to files as values.
    lambda_index : int
        The index of the current lambda value.
    replica_tops : [str]
        A list of absolute paths to additional replica topology files for REST(2).
    """
    def __init__(self, name, gro_file, top_file, lambda_index, work_dir=None, replica_top_files=None):
        if work_dir is None:
            work_dir = _os.getcwd()
        self.files = {}
        self.files["gro"] = _os.path.abspath(gro_file)
        self.files["top"] = _os.path.abspath(top_file)
        self._workdir = _fileio.Dir("%s/%s/Lambda_%d" % (work_dir, name, lambda_index))
        self.lambda_index = lambda_index
        self.replica_tops = []
        if replica_top_files not in [None, []]:
            self.replica_tops = ["%s/%s" % (_os.getcwd(), file) for file in replica_top_files]

    def runSimulation(self, name, protocol, replex=None, n_cores=None, **replica_params):
        """
        Runs a single GROMACS simulation.

        Parameters
        ----------
        name : str
            The name of this simulation.
        protocol : Protocaller.Protocol.Protocol
            The input MD protocol.
        replex : int or None, optional
            Attempts replica exchange after replex number of steps. None means no replica exchange.
        n_cores : int or None, optional
            Number of cores used. Default: the same as the number of replicas.
        replica_params
            Keyword arguments specifying parameters which are to be set for each replica.
        """
        print("Running %s..." % name)
        with self._workdir:
            with _fileio.Dir(name, overwrite=True):
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
                    if ext == "tpr": continue
                    self.files[ext] = output_file


class GMXSerialRuns:
    """
    A wrapper for a several consecutive GROMACS runs.

    Parameters
    ----------
    name : str
        The name of the run.
    gro_file : str
        Path to the input GRO file.
    top_file : str
        Path to the input TOP file.
    work_dir : str, optional
        Initialises the run in a custom directory.
    replica_top_files : [str], optional
        Initialises replica_tops.
    lambda_dict
        Initialises lambda_dict

    Attributes
    ----------
    lambda_dict : dict
        A dictionary of lists of lambda values. All of them must be of the same length.
    gmx_run_list : [ProtoCaller.Simulation.GMXSingleRun]
        Contains an instance of GMXSingleRun for every lambda value.
    protocols : [(str, ProtoCaller.Protocol.Protocol, dict)]
        A list of tuples corresponding to the protocol name, the Protocol object itself and any extra run options.
    """
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
        """
        Adds a protocol for execution.

        Parameters
        ----------
        name : str
            Protocol name.
        use_preset : str, None
            Which default preset to use. One of: "minimisation", "equilibration_nvt", "equilibration_npt", "production"
            "vacuum".
        run_options : dict, optional
            Any extra run options to be forwarded to GMXSingleRun.
        kwargs
            Keyword arguments to be passed on to Protocol.__init__.
        """
        run_options = run_options if run_options is not None else {}
        self.protocols += [(name,
                            _Protocol.Protocol(use_preset=use_preset, **kwargs, **self.lambda_dict),
                            run_options)]

    def runSimulations(self):
        """Runs all protocols consecutively for each lambda value."""
        for i, gmx_run in enumerate(self.gmx_run_list):
            print("Running simulation %d..." % (i + 1))
            for name, protocol, kwargs in self.protocols:
                try:
                    gmx_run.runSimulation(name, protocol, **kwargs)
                except:
                    "An error occurred. Check the log. Continuing..."


class GMX_REST_FEP_Runs:
    """
    A wrapper for performing a REST2 run in GROMACS.

    Parameters
    ----------
    name : str
        The name of the run.
    gro_file : str
        Path to the input GRO file.
    top_filse : [str]
        Paths to the input scaled TOP files.
    work_dir : str, optional
        Initialises the run in a custom directory.
    lambda_dict
        Initialises lambda_dict

    Attributes
    ----------
    files : dict
        A dictionary with extensions as keys and absolute paths to files as values.
    lambda_dict : dict
        A dictionary of lists of lambda values. All of them must be of the same length.
    protocols : [ProtoCaller.Protocol.Protocol]
        A list of Protocol objects corresponding to the already run protocols.
    mbar_data : numpy.ndarray
        Data which is to be passed to pymbar.
    """
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
        self._workdir = _fileio.Dir("%s/%s" % (work_dir, name))
        self.lambda_dict = lambda_dict
        self.protocols = []
        self.mbar_data = []

    def runSimulation(self, name, parallel=True, use_mpi=False, use_preset=None, replex=None, n_cores=None, n_nodes=1,
                      **protocol_params):
        """
        Runs a REST2 simulation in GROMACS.

        Parameters
        ----------
        name : str
            The name of the simulation.
        parallel : bool, optional
            Whether to run the simulations in parallel, using -multi.
        use_mpi : bool, optional
            Whether to use ProtoCaller.GROMACSEXE or GROMACSMPIEXE.
        use_preset : str, None
            Which default preset to use. One of: "minimisation", "equilibration_nvt", "equilibration_npt", "production"
            "vacuum".
        replex : int or None, optional
            Attempts replica exchange after replex number of steps. None means no replica exchange.
        n_cores : int or None, optional
            Number of cores used. Default: the same as the number of replicas.
        n_nodes : int, optional
            Number of nodes used.
        protocol_params
            Keyword arguments passed to ProtoCaller.Protocol.Protocol
        """
        if n_cores is None: n_cores = self.lambda_size
        ppn = n_cores // n_nodes
        if n_nodes > 1: use_mpi = True
        if replex is not None: parallel=True

        protocol = _Protocol.Protocol(use_preset=use_preset, **protocol_params, **self.lambda_dict)
        self.protocols += [protocol]

        print("Running %s..." % name)
        with self._workdir:
            with _fileio.Dir(name, overwrite=True):
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
                        if not use_mpi:
                            mdrun_command = "{0} mdrun -s {1}.tpr -deffnm {1}".format(_PC.GROMACSEXE, filebase)
                        else:
                            # disable dynamic load balancing due to GROMACS bug
                            mdrun_command = "{0} -np {1} --map-by ppr:{2}:node {3} " \
                                            "mdrun -dlb no -s {4}.tpr -deffnm {4}".format( _PC.MPIEXE, n_cores, ppn,
                                                                                           _PC.GROMACSMPIEXE, filebase)
                        _runexternal.runExternal(mdrun_command, procname="gmx mdrun")

                        self.files[i] = {"top" : self.files[i]["top"],}
                        output_files = _glob.glob("%s.*" % filebase)
                        for output_file in output_files:
                            ext = output_file.split(".")[-1].lower()
                            if ext == "tpr": continue
                            self.files[i][ext] = _os.path.abspath(output_file)

                if parallel:
                    open("plumed.dat", "a").close()
                    filebase = name + "_"
                    # disable dynamic load balancing due to GROMACS bug
                    mdrun_command = "{0} -np {1} --map-by ppr:{2}:node {3} mdrun -dlb no -multi {4} -s {5}.tpr " \
                                    "-deffnm {5} ".format(_PC.MPIEXE, n_cores, ppn, _PC.GROMACSMPIEXE, self.lambda_size,
                                                          filebase)
                    if replex is not None:
                        mdrun_command += " -plumed plumed.dat -replex {} -hrex".format(replex)

                    _runexternal.runExternal(mdrun_command, procname="gmx mdrun")

                    for i in range(self.lambda_size):
                        filebase = "%s_%d" % (name, i)
                        self.files[i] = {"top": self.files[i]["top"], }
                        output_files = _glob.glob("%s.*" % filebase)
                        for output_file in output_files:
                            ext = output_file.split(".")[-1].lower()
                            if ext == "tpr": continue
                            self.files[i][ext] = _os.path.abspath(output_file)

    def generateMBARData(self, n_cores=None, n_nodes=1, cont=True):
        """
        Generates input data for MBAR using gmx mdrun -rerun

        Parameters
        ----------
         n_cores : int or None, optional
            Number of cores used. Default: the same as the number of replicas.
        n_nodes : int, optional
            Number of nodes used.
        cont : bool, optional
            Whether to overwrite existing files or to continue from the last one.

        Returns
        -------
        mbar_data : numpy.ndarray
            Data which is to be passed to pymbar.
        """
        # GROMACS might generate a bunch of warnings when we apply a non-dummy Hamiltonian to a trajectory with dummies
        # due to clashes and backup too many structures. Here we suppress this backing up
        gmx_suppress_dump = _os.environ["GMX_SUPPRESS_DUMP"] if "GMX_SUPPRESS_DUMP" in _os.environ.keys() else None
        _os.environ["GMX_SUPPRESS_DUMP"] = "1"
        self.mbar_data = []

        print("Generating Energy files...")
        with self._workdir:
            with _fileio.Dir("MBAR"):
                for i, file in enumerate(self.files):
                    self.mbar_data += [[]]
                    gro, top = file["gro"], file["top"]
                    for j, file in enumerate(self.files):
                        run = True
                        trr = file["trr"]
                        filebase = "Energy_%d_%d" % (i, j)
                        # we only overwrite the last generated energies, because they might be incomplete
                        if cont:
                            filebase_next = "Energy_%d_%d" % (i + (j + 1) // len(self.files) , (j + 1) % len(self.files))
                            files_current = _glob.glob("%s/%s.xvg" % (_os.getcwd(), filebase))
                            files_next = _glob.glob("%s/%s.xvg" % (_os.getcwd(), filebase_next))
                            if len(files_current) and len(files_next):
                                run = False
                        if run:
                            #use the last protocol or create a default one
                            if len(self.protocols):
                                protocol = self.protocols[-1]
                            else:
                                protocol = _Protocol.Protocol(use_preset="default", **self.lambda_dict)
                            protocol.init_lambda_state = i
                            protocol.skip_positions = 0
                            protocol.skip_velocities = 0
                            protocol.skip_forces = 0
                            protocol.write_derivatives = False
                            protocol.random_velocities = False
                            protocol.constraint = "no"
                            protocol.__setattr__("continuation", "yes")
                            protocol.__setattr__("dhdl-print-energy", "potential")
                            protocol.__setattr__("calc-lambda-neighbors", "0")
                            protocol.__setattr__("calc-lambda-neighbors", "0")

                            mdp = protocol.write("GROMACS", filebase=filebase)
                            tpr = filebase + ".tpr"

                            grompp_command = "%s grompp -maxwarn 10 -f '%s' -c '%s' -p '%s' -o '%s'" % (
                                _PC.GROMACSEXE, mdp, gro, top, tpr)
                            _runexternal.runExternal(grompp_command, procname="gmx grompp")

                            if n_cores is None:
                                mdrun_command = "%s mdrun -s '%s' -rerun '%s' -deffnm '%s'" % (_PC.GROMACSEXE, tpr, trr,
                                                                                               filebase)
                            else:
                                ppn = n_cores // n_nodes
                                mdrun_command = "{0} -np {1} --map-by ppr:{2}:node {3} mdrun -s '{4}' -rerun '{5}' " \
                                                "-deffnm {6}".format(_PC.MPIEXE, n_cores, ppn, _PC.GROMACSMPIEXE, tpr,
                                                                     trr, filebase)


                            _runexternal.runExternal(mdrun_command, procname="gmx mdrun")

                        self.mbar_data[i] += list(_MDAnalysis.auxiliary.XVG.XVGReader(filebase + ".xvg").
                                                  _auxdata_values[:, 1])

        # restore the original environment variable
        del _os.environ["GMX_SUPPRESS_DUMP"]
        if gmx_suppress_dump is not None: _os.environ["GMX_SUPPRESS_DUMP"] = gmx_suppress_dump

        return self.mbar_data

    def runMBAR(self, n_points_to_ignore=1, temperature=298):
        """
        Runs pymbar and outputs to stdout.

        Parameters
        ----------
        n_points_to_ignore : int, optional
            Number of initial points which are to be excluded in the calculation.
        temperature : float, optional
            The temperature at which the simulations were originally ran.
        """
        RT = 0.008314 * temperature
        if not self.mbar_data:
            raise ValueError("No MBAR data to analyse")

        if n_points_to_ignore:
            mbar_data = []
            for row in self.mbar_data:
                mbar_data += [[x for i, x in enumerate(row)
                               if i % (len(self.mbar_data[0]) / len(self.mbar_data)) >= n_points_to_ignore]]
        else:
            mbar_data = self.mbar_data

        N = len(mbar_data[0])
        size = len(mbar_data)
        N_k = [N / size] * size
        u_kln = _np.reshape(mbar_data, (size, N)) / RT
        mbar = _pymbar.mbar.MBAR(u_kn=u_kln, N_k=N_k, verbose=True, initialize="BAR", maximum_iterations=10000)
        f, df, _ = mbar.getFreeEnergyDifferences()

        print()
        print("Free Energy change from A to B in kJ/mol: ", f[0][size - 1] * RT)
        print("Uncertainty in kJ/mol: ", df[0][size - 1] * RT)
        print("Free Energy change from A to B in kcal/mol: ", f[0][size - 1] * RT / 4.184)
        print("Uncertainty in kcal/mol: ", df[0][size - 1] * RT / 4.184)
