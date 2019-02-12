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
import ProtoCaller.Wrappers.BioSimSpacewrapper as _BSSwrap
import ProtoCaller.Utils.stdio as _stdio


class Ensemble:
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
        self.protein = protein
        self.ligand_ref = ligand_ref
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
                raise ValueError("Value %s not supported. Supported values: " % value, _PC.ENGINES)
        elif key in ["box_length_complex", "box_length_morph", "ion_conc", "shell"]:
            value = float(value)
        elif key == "morphs":
            value = PerturbationList(value) if not isinstance(value, PerturbationList) else value
        elif key == "protein":
            if value is not None:
                value = Protein(value) if not isinstance(value, Protein) else value
            self.systems_prep = {}
        elif key == "ligand_ref":
            if value is not None:
                if self.protein is None:
                    raise ValueError("Cannot assign ligand_id without a protein")
                else:
                    self.protein.ligand_ref = value

        # setter
        super(Ensemble, self).__setattr__("_" + key, value)

    def prepareComplexes(self, replica_temps=None, intermediate_files=False, store_complexes=False, output_files=True):
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
        if systems is None: systems = self.systems_prep
        # TODO support other engines
        print("Saving solvated complexes as GROMACS...")
        with self.workdir:
            for name, (morph, complexes) in systems.items():
                with _fileio.Dir(name):
                    _IO.GROMACS.saveAsGromacs("morph", morph)

                    if len(complexes) == 1:
                        _IO.GROMACS.saveAsGromacs("complex_final", complexes[0])
                    else:
                        for j, complex in enumerate(complexes):
                            _IO.GROMACS.saveAsGromacs("complex_final%d" % j, complex)
