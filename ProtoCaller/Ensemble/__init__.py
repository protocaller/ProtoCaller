import os as _os
import tempfile as _tempfile
import warnings as _warnings

import ProtoCaller as _PC
from .Ligand import Ligand
from .Perturbation import Perturbation
from .PerturbationList import PerturbationList
from .Protein import Protein
import ProtoCaller.IO as _IO
import ProtoCaller.Parametrise as _parametrise
import ProtoCaller.Solvate as _solvate
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Wrappers.biosimspacewrapper as _BSSwrap
import ProtoCaller.Utils.stdio as _stdio


class Ensemble:
    """
    The main class responsible for handling batches of complexes in ProtoCaller.

    Parameters
    ----------
    engine : str
        Initialises the engine.
    box_length_complex : float
        Initialises box_length_complex.
    box_length_morph : float
        Initialises box_length_morph.
    ion_conc : float
        Initialises ion_conc.
    shell : float
        Initialises shell.
    neutralise : bool
        Initialises neutralise.
    centre : bool
        Initialises centre.
    protein_ff : str
        The protein force field. Used to initialise params.
    ligand_ff : str
        The ligand force field. Used to initialise params.
    water_ff : str
        The water / simple ion force field. Used to initialise params.
    protein : ProtoCaller.Ensemble.Protein.Protein
        Initialises protein.
    ligand_ref : str or ProtoCaller.Ensemble.Ligand.Ligand or None or False
        Initialises protein.ligand_ref.
    morphs : ProtoCaller.Ensemble.PerturbationList.PerturbationList or [ProtoCaller.Ensemble.Perturbation.Perturbation] \
             or [[ProtoCaller.Ensemble.Ligand.Ligand, ProtoCaller.Ensemble.Ligand.Ligand]]
        Initialises morphs.
    workdir : str
        Initialises workdir.

    Attributes
    ----------
    workdir : ProtoCaller.Utils.fileio.Dir
        The working directory of the ensemble.
    engine : str
        The engine for which the input files are created. One of: "GROMACS".
    box_length_complex : float
        Size of the solvated complex box in nm. Cubic shape is assumed.
    box_length_morph : float
        Size of the solvated morph box in nm. Cubic shape is assumed.
    shell : float
        Places a layer of water of the specified thickness in nm around the solute.
    neutralise : bool
        Whether to add counterions to neutralise the system.
    ion_conc : float
        Ion concentration of NaCl in mol/L.
    centre : bool
        Whether to centre the system.
    params : ProtoCaller.Parametrise.Params
        Force field parameters.
    protein : ProtoCaller.Ensemble.Protein
        The common protein.
    morphs : ProtoCaller.Ensemble.PerturbationList.PerturbationList
         A list of all perturbations.
    systems_prep : dict
        The prepared systems.
    """
    _common_properties = ["engine", "box_length_complex", "box_length_morph", "ion_conc", "shell", "protein",
                          "ligand_ref", "morphs"]

    def __init__(self, engine="GROMACS", box_length_complex=9, box_length_morph=4, ion_conc=0.154, shell=0,
                 neutralise=True, centre=True, protein_ff=_PC.AMBERDEFAULTPROTEINFF, ligand_ff=_PC.AMBERDEFAULTLIGANDFF,
                 water_ff=_PC.AMBERDEFAULTWATERFF, protein=None, ligand_ref=None, morphs=None, workdir=None):
        # work directory
        self.workdir = _fileio.Dir(workdir) if workdir else _fileio.Dir(".")

        # writing options
        self.engine = engine

        # solvation options
        self.box_length_complex = box_length_complex
        self.box_length_morph = box_length_morph
        self.ion_conc = ion_conc
        self.shell = shell
        self.neutralise = neutralise
        self.centre = centre

        # parametrisation options
        self.params = _parametrise.Params(protein_ff=protein_ff, ligand_ff=ligand_ff, water_ff=water_ff)

        # molecular systems
        if not isinstance(protein, Protein):
            self.protein = Protein(protein, ligand_ref=ligand_ref)
        else:
            self.protein = protein
        self.morphs = morphs
        self._complex_template = None
        self.systems_prep = {}

    def __getattr__(self, item):
        if item not in self._common_properties:
            return super().__getattribute__(item)
        else:
            return super().__getattribute__("_" + item)

    def __setattr__(self, key, value):
        if key not in self._common_properties:
            super(Ensemble, self).__setattr__(key, value)

        if isinstance(value, str): value = value.strip()

        # checks
        if key == "engine":
            value = value.upper()
            if value not in _PC.ENGINES:
                raise ValueError("Value {} not supported. Supported values: {}".format(value, _PC.ENGINES))
        elif key in ["box_length_complex", "box_length_morph", "ion_conc", "shell"]:
            value = float(value)
        elif key == "morphs":
            value = PerturbationList(value) if not isinstance(value, PerturbationList) else value
        elif key == "protein":
            if value is not None:
                value = Protein(value) if not isinstance(value, Protein) else value
            self.systems_prep = {}

        # setter
        super(Ensemble, self).__setattr__("_" + key, value)

    def prepareComplexes(self, replica_temps=None, intermediate_files=False, store_complexes=False, output_files=True):
        """
        Batch prepares all complexes with an option to output files for REST(2).

        Parameters
        ----------
        replica_temps : [float] or None
            A list of replica temperatures. Everything is normalised with respect to the lowest temperature. None
            means only output the normal files.
        intermediate_files : bool
            Whether to store all intermediate files.
        store_complexes : bool
            Whether to store the final complexes as a dictionary of BioSimSpace System objects.
        output_files : bool
            Whether to write output files immediately or later via saveSystems
        """
        # make sure the proteins / ligands are parametrised before proceeding
        with self.workdir:
            _stdio.stdout_stderr()(self.protein.parametrise)(params=self.params, reparametrise=False)
            for morph in self.morphs:
                _stdio.stdout_stderr()(morph.ligand1.parametrise)(params=self.params, reparametrise=False)
                _stdio.stdout_stderr()(morph.ligand2.parametrise)(params=self.params, reparametrise=False)

            # take care of the replicas if there are any
            if replica_temps is None or len(replica_temps) == 1:
                scales = [1]
            else:
                sorted_list = sorted(replica_temps)
                if sorted_list != replica_temps:
                    _warnings.warn("Input replica temperatures were not in ascending order. Sorting...")
                    replica_temps = sorted_list
                scales = [replica_temps[0] / elem for elem in replica_temps]

            for i, morph in enumerate(self.morphs):
                if intermediate_files:
                    curdir = _fileio.Dir(morph.name, overwrite=True)
                else:
                    name = _os.path.basename(_tempfile.TemporaryDirectory().name)
                    curdir = _fileio.Dir(name, overwrite=True, temp=True)

                with curdir:
                    print("Creating morph %s..." % morph.name)
                    _stdio.stdout_stderr()(morph.alignAndCreateMorph)(self.protein.ligand_ref)

                    complexes = [self.protein.complex_template + morph.get_morph(self.protein.ligand_ref)]

                    # solvate and save the prepared complex and morph with the appropriate box size
                    print("Solvating...")
                    complexes = [_solvate.solvate(complexes[0], self.params, box_length=self.box_length_complex,
                                                  shell=self.shell, neutralise=self.neutralise, ion_conc=self.ion_conc,
                                                  centre=self.centre, work_dir=curdir.path, filebase="complex")]
                    morph_sol = _solvate.solvate(morph.get_morph(self.protein.ligand_ref), self.params,
                                                 box_length=self.box_length_morph, shell=self.shell,
                                                 neutralise=self.neutralise, ion_conc=self.ion_conc, centre=self.centre,
                                                 work_dir=curdir.path, filebase="morph")

                    # rescale complexes for replica exchange if needed
                    if len(scales) != 1:
                        print("Creating replicas...")
                        complexes = [_BSSwrap.rescaleSystemParams(complexes[0], scale, includelist=["Merged_Molecule"])
                                     for scale in scales]

                    if store_complexes:
                        self.systems_prep[morph.name] = (morph_sol, complexes)

                if output_files:
                    self.saveSystems({morph.name: (morph_sol, complexes)})

    def saveSystems(self, systems=None):
        """
        Saves all systems.

        Parameters
        ----------
        systems : dict or None
            Input systems to be saved. None means use self.systems_prep.
        """
        if systems is None: systems = self.systems_prep
        # TODO support other engines
        print("Saving solvated complexes as {}...".format(self.engine))
        with self.workdir:
            if self.engine == "GROMACS":
                savefunc = _IO.GROMACS.saveAsGromacs
            elif self.engine == "NAMD":
                savefunc = _IO.NAMD.saveAsNamd
            else:
                raise ValueError("Value {} not supported. Supported values: "
                                 "{}".format(self.engine, _PC.ENGINES))

            for name, (morph, complexes) in systems.items():
                with _fileio.Dir(name):
                    savefunc("morph", morph)
                    if len(complexes) == 1:
                        savefunc("complex_final", complexes[0])
                    else:
                        for j, complex in enumerate(complexes):
                            savefunc("complex_final%d" % j, complex)
