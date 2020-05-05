import ProtoCaller as _PC
if not _PC.BIOSIMSPACE:
    raise ImportError("BioSimSpace module cannot be imported")

from collections.abc import Iterable as _Iterable
import logging as _logging
import os as _os
import tempfile as _tempfile
import warnings as _warnings

import BioSimSpace as _BSS
import rdkit.Chem as _Chem

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


class Ensemble:
    """
    The main class responsible for handling batches of complexes in ProtoCaller.

    Parameters
    ----------
    engine : str
        Initialises the engine.
    box_length_complex : float, iterable
        Initialises box_length_complex.
    box_length_morph : float, iterable
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
    box_length_complex : float, iterable
        Size of the solvated complex box in nm. Cubic shape is assumed.
    box_length_morph : float, iterable
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
                raise ValueError("Value %s not supported. Supported values: " % value, _PC.ENGINES)
        elif key in ["ion_conc", "shell"]:
            value = float(value)
        elif key in ["box_length_complex", "box_length_morph"]:
            if isinstance(value, _Iterable):
                value = list(value)
            else:
                value = float(value)
        elif key == "morphs":
            value = PerturbationList(value) if not isinstance(value, PerturbationList) else value
        elif key == "protein":
            if value is not None:
                value = Protein(value) if not isinstance(value, Protein) else value
            self.systems_prep = {}

        # setter
        super(Ensemble, self).__setattr__("_" + key, value)

    def prepareComplexes(self, replica_temps=None, scale_dummy_bonds=1,
                         dummy_bond_smarts="[*]~[*]", intermediate_files=False,
                         store_complexes=False, output_files=True):
        """
        Batch prepares all complexes with an option to output files for REST(2).

        Parameters
        ----------
        replica_temps : [float] or None
            A list of replica temperatures. Everything is normalised with respect to the lowest temperature. None
            means only output the normal files.
        scale_dummy_bonds : float
            Sets the dummy bond length distance as a fraction of the real bond length distance.
        dummy_bond_smarts : str
            SMARTS string which indicates which dummy bonds are to be affected by scale_dummy_bonds.
        intermediate_files : bool
            Whether to store all intermediate files.
        store_complexes : bool
            Whether to store the final complexes as a dictionary of BioSimSpace System objects.
        output_files : bool
            Whether to write output files immediately or later via saveSystems.
        """
        # make sure the proteins / ligands are parametrised before proceeding
        with self.workdir:
            self.protein.parametrise(params=self.params, reparametrise=False)
            for morph in self.morphs:
                morph.ligand1.parametrise(params=self.params, reparametrise=False)
                morph.ligand2.parametrise(params=self.params, reparametrise=False)

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
                    _logging.info("Creating morph %s..." % morph.name)
                    morph_BSS, mcs =  morph.alignAndCreateMorph(self.protein.ligand_ref)
                    morph_BSS = _BSS._SireWrappers.System(morph_BSS)
                    box = self.protein.complex_template._sire_object.property(
                        "space")
                    morph_BSS._sire_object.setProperty("space", box)

                    # here we scale the equilibrium bond lengths if needed
                    if scale_dummy_bonds != 1:
                        n1 = morph.ligand1.molecule.GetNumAtoms()
                        n2 = morph.ligand2.molecule.GetNumAtoms()
                        inv_map = {y: x for x, y in mcs}
                        du2 = [i for i in range(n2) if i not in inv_map.keys()]
                        inv_map = {**{x: n1 + y
                                      for x, y in zip(du2, range(n2 - n1))},
                                   **inv_map}
                        mcs_smarts = _Chem.MolFromSmarts(dummy_bond_smarts)
                        matches_lig1 = morph.ligand1.molecule.\
                            GetSubstructMatches(mcs_smarts)
                        matches_lig2 = morph.ligand2.molecule.\
                            GetSubstructMatches(mcs_smarts)

                        # here we take care of the index transformation
                        matches_lig2 = [(inv_map[x[0]], inv_map[x[1]])
                                        for x in matches_lig2
                                        if set(x).issubset(inv_map.keys())]
                        matches_total = [*matches_lig1, *matches_lig2]
                        morph_BSS = _BSSwrap.rescaleBondedDummies(
                            morph_BSS, scale_dummy_bonds,
                            {"Merged_Molecule": matches_total}
                        )

                    complexes = [self.protein.complex_template + morph_BSS]

                    # solvate and save the prepared complex and morph with the appropriate box size
                    _logging.info("Solvating...")
                    complexes = [_solvate.solvate(complexes[0], self.params, box_length=self.box_length_complex,
                                                  shell=self.shell, neutralise=self.neutralise, ion_conc=self.ion_conc,
                                                  centre=self.centre, work_dir=curdir.path, filebase="complex")]
                    morph_sol = _solvate.solvate(morph_BSS, self.params,
                                                 box_length=self.box_length_morph, shell=self.shell,
                                                 neutralise=self.neutralise, ion_conc=self.ion_conc, centre=self.centre,
                                                 work_dir=curdir.path, filebase="morph")

                    # rescale complexes for replica exchange if needed
                    if len(scales) != 1:
                        _logging.info("Creating replicas...")
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
        _logging.info("Saving solvated complexes as GROMACS...")
        with self.workdir:
            for name, (morph, complexes) in systems.items():
                with _fileio.Dir(name):
                    _IO.GROMACS.saveAsGromacs("morph", morph)

                    if len(complexes) == 1:
                        _IO.GROMACS.saveAsGromacs("complex_final", complexes[0])
                    else:
                        for j, complex in enumerate(complexes):
                            _IO.GROMACS.saveAsGromacs("complex_final%d" % j, complex)
