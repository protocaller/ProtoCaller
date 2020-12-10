import glob as _glob
import logging as _logging
import os as _os
import subprocess as _subprocess

import MDAnalysis as _MDAnalysis
import numpy as _np
import pymbar as _pymbar

import ProtoCaller as _PC
import ProtoCaller.Protocol as _Protocol
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Utils.runexternal as _runexternal


class RunGMX:
    """
    A wrapper for performing a run in GROMACS.

    Parameters
    ----------
    name : str
        The name of the run.
    gro_files : str
        Path(s) to the input GRO file(s).
    top_files : [str] or str
        Path(s) to the input TOP file(s).
    work_dir : str, optional
        Initialises the run in a custom directory.
    lambda_dict
        Initialises lambda_dict.

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
    def __init__(self, name, gro_files, top_files, work_dir=None, **lambda_dict):
        try:
            len_arr, = {len(arr) for arr in lambda_dict.values() if len(arr)}
        except ValueError:
            raise ValueError("Need lists of the same non-zero size")

        if not isinstance(gro_files, list):
            gro_files = [gro_files] * len_arr
        else:
            assert len(gro_files) == len_arr, "Number of GRO files does not match number of lambda values"

        if not isinstance(top_files, list):
            top_files = [top_files] * len_arr
        else:
            assert len(top_files) == len_arr, "Number of TOP files does not match number of lambda values"

        if work_dir is None:
            work_dir = _os.getcwd()

        self.files = [{} for _ in range(len_arr)]
        for f, gro_file, top_file in zip(self.files, gro_files, top_files):
            f["gro"] = _os.path.abspath(gro_file)
            f["top"] = _os.path.abspath(top_file)
        self._workdir = _fileio.Dir(f"{work_dir}/{name}")
        self.lambda_dict = lambda_dict
        self.protocols = []
        self.mbar_data = []

    @property
    def lambda_size(self):
        """int: Returns the number of lambda windows."""
        return len(self.files)

    @staticmethod
    def _dict_to_arguments(argument_str, argument_dict, i):
        for k, v in argument_dict.items():
            argument_str += f" -{k}"
            if v is not None:
                argument_str += f" {str(v).format(i)}"
        return argument_str

    @staticmethod
    def _gmx_version(executable):
        version = _subprocess.check_output(f"{executable} -version", shell=True).decode()
        version = int(next(x.split(".")[0] for x in version.split("\n")[0].split(" ") if x.split(".")[0].isnumeric()))
        return version

    def _update_files(self, i, filebase):
        self.files[i] = {"top": self.files[i]["top"]}
        output_files = _glob.glob(f"{filebase}.*")
        for output_file in output_files:
            ext = output_file.split(".")[-1].lower()
            if ext != "tpr":
                self.files[i][ext] = _os.path.abspath(output_file)

    def runSimulation(self, name, multi=False, multidir=False, single_lambda=None, use_mpi=False, mdrun_mpi=False,
                      gpu_id=None, use_preset=None, replex=None, plumed_file=None, n_cores_per_process=None, n_nodes=1,
                      n_processes=None, dlb=False, gmx_kwargs=None, **protocol_params):
        """
        Runs a simulation in GROMACS.

        Parameters
        ----------
        name : str
            The name of the simulation.
        multi : bool, optional
            Whether to run the simulations in parallel, using -multi.
        multidir : bool, optional
            Whether to run the simulations in parallel, using -multidir. Overrides multi.
        single_lambda : None, int, optional
            An integer runs a simulation at a single lambda value and overrides multi and replex.
            None runs all lambda values.
        use_mpi : bool, optional
            Whether to use ProtoCaller.GROMACSEXE or GROMACSMPIEXE.
        mdrun_mpi : bool, optional
            Whether to call mdrun as "mdrun" or as "mdrun_mpi".
        gpu_id : str, optional
            Which GPU id to be used, if applicable.
        use_preset : str, Protocaller.Protocol.Protocol, None
            Which default preset to use. One of: "minimisation", "equilibration_nvt", "equilibration_npt", "production"
            "vacuum". You can alternatively pass a custom protocol here, in which case **protocol_params will not be
            used.
        replex : int or None, optional
            Attempts replica exchange after replex number of steps using PLUMED. None means no replica exchange.
            Overrides use_mpi and multi if not None.
        plumed_file : str, optional
            If replex is True, it uses this file as a configuration for PLUMED, otherwise an empty file is used.
        n_cores_per_process : int or None, optional
            Number of cores used per process. Default: let GROMACS decide.
        n_nodes : int, optional
            Number of physical nodes used.
        n_processes : int, optional
            Number of processes. Default: the same as the number of simulation.
        dlb : bool
            Whether to enable dynamic load balancing. Default is False due to some possible instabilities with
            gmx_mpi mdrun in some cases.
        gmx_kwargs : dict
            Additional arguments to be passed to mdrun. The keys of the dictionary need to be the name of the option,
            e.g. "cpt" for checkpoint interval, while the values need to be the value of the option if it permits one
            or None if it doesn't. If the values contain "{}" while the user is running the lambda windows in serial,
            this will be replaced by the lambda number.
        protocol_params
            Keyword arguments passed to ProtoCaller.Protocol.Protocol.
        """
        # perform some checks and initialise some default values
        if n_processes is None:
            n_processes = 1 if single_lambda else self.lambda_size
        if n_nodes > 1 or multi or multidir:
            use_mpi = True
        if single_lambda is not None:
            replex, multi, multidir = None, False, False
        if replex is not None:
            use_mpi, multi = True, True

        if mdrun_mpi:
            MDRUNEXE = _os.path.dirname(_PC.GROMACSMPIEXE) + "/mdrun_mpi"
        elif use_mpi:
            MDRUNEXE = _PC.GROMACSMPIEXE + " mdrun"
        else:
            MDRUNEXE = _PC.GROMACSEXE + " mdrun"

        if use_mpi or mdrun_mpi:
            base_mdrun_command = f"{_PC.MPIEXE} -np {n_processes} --map-by ppr:{n_processes // n_nodes}:node {MDRUNEXE}"
        else:
            base_mdrun_command = MDRUNEXE

        gmx_kwargs = {} if gmx_kwargs is None else gmx_kwargs
        gmx_kwargs["dlb"] = "yes" if dlb else "no"
        if n_cores_per_process is not None:
            gmx_kwargs["ntomp"] = n_cores_per_process
        if gpu_id is not None:
            gmx_kwargs["gpu_id"] = gpu_id

        # this is due to incompatibility between GROMACS 2019+ and multi
        multidir = True if self._gmx_version(_PC.GROMACSMPIEXE) >= 2019 and multi else multidir

        # initialise the protocol
        if isinstance(use_preset, _Protocol.Protocol):
            protocol = use_preset
        else:
            protocol = _Protocol.Protocol(use_preset=use_preset, **protocol_params, **self.lambda_dict)
        self.protocols += [protocol]

        # run the simulations
        _logging.info(f"Running {name}...")
        with self._workdir:
            with _fileio.Dir(name, overwrite=False):
                # run single lambda if needed
                if single_lambda is not None:
                    it = [single_lambda]
                else:
                    it = range(self.lambda_size)

                # call grompp for each lambda
                for i in it:
                    filebase = name if multidir else f"{name}_{i}"
                    dirname = f"Lambda_{i}" if multidir else "."
                    with _fileio.Dir(dirname):
                        protocol.init_lambda_state = i
                        self.files[i]["mdp"] = _os.path.abspath(protocol.write(engine="GROMACS", filebase=filebase))
                        grompp_args = {
                            "-f": "mdp",
                            "-c": "gro",
                            "-p": "top",
                            "-t": "cpt",
                        }

                        grompp_command = f"{_PC.GROMACSEXE} grompp -maxwarn 10 -o {filebase}.tpr"
                        for grompp_arg, filetype in grompp_args.items():
                            if filetype in self.files[i].keys():
                                grompp_command += f" {grompp_arg} '{self.files[i][filetype]}'"
                        _runexternal.runExternal(grompp_command, procname="gmx grompp")

                # call the simulations consecutively for each lambda if multi is not specified
                if not (multi or multidir):
                    for i in it:
                        filebase = f"{name}_{i}"
                        gmx_kwargs["s"] = f"'{filebase}.tpr'"
                        gmx_kwargs["deffnm"] = f"'{filebase}'"

                        # run GROMACS
                        mdrun_command = self._dict_to_arguments(base_mdrun_command, gmx_kwargs, i)
                        _runexternal.runExternal(mdrun_command, procname="gmx mdrun")

                        # update files after the run
                        self._update_files(i, filebase)
                # alternatively, run all simulations in parallel
                else:
                    filebase = name if multidir else name + "_"
                    gmx_kwargs["s"] = f"'{filebase}.tpr'"
                    gmx_kwargs["deffnm"] = f"'{filebase}'"
                    if multidir:
                        gmx_kwargs["multidir"] = " ".join(f"Lambda_{i}" for i in it)
                    else:
                        gmx_kwargs["multi"] = f"{self.lambda_size}"

                    # run replica exchange with PLUMED if needed
                    if replex is not None:
                        if plumed_file is None:
                            plumed_file = "plumed.dat"
                            open(plumed_file, "a").close()
                        gmx_kwargs["plumed"] = f"'{plumed_file}'"
                        gmx_kwargs["replex"] = replex
                        gmx_kwargs["hrex"] = None

                    # run GROMACS
                    mdrun_command = self._dict_to_arguments(base_mdrun_command, gmx_kwargs, i)
                    _runexternal.runExternal(mdrun_command, procname="gmx mdrun")

                    # update files after the run
                    for i in it:
                        filebase = name if multidir else f"{name}_{i}"
                        dirname = f"Lambda_{i}" if multidir else "."
                        with _fileio.Dir(dirname):
                            self._update_files(i, filebase)

    def generateMBARData(self, n_cores=None, n_nodes=1, cont=True):
        """
        Manually generates input data for MBAR using gmx mdrun -rerun using N
        topologies and N trajectories generated from N simulations (N squared
        energy evaluations are needed). Doesn't work with constrained
        simulations. For most cases one should instead use the XVG files
        readily generated by GROMACS during the simulation.

        Parameters
        ----------
         n_cores : int or None, optional
            Number of cores used. Default: the same as the number of replicas.
        n_nodes : int, optional
            Number of nodes used.
        cont : bool, optional
            Whether to overwrite existing files or to continue from the last
            one.

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

        _logging.info("Generating Energy files...")
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
                            filebase_next = "Energy_%d_%d" % (i + (j + 1) // len(self.files), (j + 1) % len(self.files))
                            files_current = _glob.glob("%s/%s.xvg" % (_os.getcwd(), filebase))
                            files_next = _glob.glob("%s/%s.xvg" % (_os.getcwd(), filebase_next))
                            if len(files_current) and len(files_next):
                                run = False
                        if run:
                            # use the last protocol or create a default one
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
        if gmx_suppress_dump is not None:
            _os.environ["GMX_SUPPRESS_DUMP"] = gmx_suppress_dump

        return self.mbar_data

    def runMBAR(self, n_points_to_ignore=1, temperature=298):
        """
        Runs pymbar and outputs to the logger. Note that this option is only
        valid after running
        :func:`~ProtoCaller.Simulation.RunGMX.generateMBARData`.
        For regular data output by GROMACS one should
        use external scripts such as `alchemical_analysis\
         <https://github.com/MobleyLab/alchemical-analysis>`_.

        Parameters
        ----------
        n_points_to_ignore : int, optional
            Number of initial points which are to be excluded in the calculation.
        temperature : float, optional
            The temperature at which the simulations were originally run.

        Returns
        -------
        f : numpy.ndarray
            The free energy matrix output by pymbar.
        df : numpy.ndarray
            The statistical uncertainty matrix output by pymbar.
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

        _logging.info("Free Energy change from A to B in kJ/mol: ", f[0][size - 1] * RT)
        _logging.info("Uncertainty in kJ/mol: ", df[0][size - 1] * RT)
        _logging.info("Free Energy change from A to B in kcal/mol: ", f[0][size - 1] * RT / 4.184)
        _logging.info("Uncertainty in kcal/mol: ", df[0][size - 1] * RT / 4.184)

        return f, df
