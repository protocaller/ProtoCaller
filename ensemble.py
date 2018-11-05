import const as _const

import os as _os
import re as _re
import tempfile as _tempfile
import warnings as _warnings

import BioSimSpace as _BSS
import numpy as _np
from rdkit import Chem as _Chem
import parmed as _pmd
import Sire.Maths as _SireMaths
import Sire.Vol as _SireVol

import IO as _IO
import pdbconnect as _pdbconnect
import Parametrise as _parametrise
import subdir as _subdir
from Wrappers import babelwrapper as _babel
from Wrappers import modellerwrapper as _modeller
from Wrappers import PDB2PQRwrapper as _PDB2PQR
from Wrappers import protosswrapper as _protoss
from Wrappers import rdkitwrapper as _rdkit
from Wrappers import runexternal as _runexternal

class Ensemble:
    def __init__(self, engine, box_length=12, ion_conc=0.154, shell=0, neutralise=True, centre=True,
                 protein_ff=_const.AMBERDEFAULTPROTEINFF, ligand_ff=_const.AMBERDEFAULTLIGANDFF,
                 water_ff=_const.AMBERDEFAULTWATERFF, protein=None, ligand_id=None, morphs=None):
        self.engine = engine
        self.box_length = box_length
        self.ion_conc = ion_conc
        self.shell = shell
        self.neutralise = neutralise
        self.centre = centre
        self.params = _parametrise.Params(protein_ff=protein_ff, ligand_ff=ligand_ff, water_ff=water_ff)
        self.protein = protein
        self.ligand_id = ligand_id
        self._ligands = {}
        self.morphs = morphs
        self._complex_template = None
        self.complexes = []

    @property
    def engine(self):
        return self._engine

    @engine.setter
    def engine(self, val):
        val = val.strip().upper()
        if val not in _const.ENGINES:
            raise ValueError("Value %s not supported. Supported values: " % val, _const.ENGINES)
        self._engine = val

    @property
    def box_length(self):
        return self._box_length

    @box_length.setter
    def box_length(self, val):
        self._box_length = int(val)

    @property
    def ion_conc(self):
        return self._ion_conc

    @ion_conc.setter
    def ion_conc(self, val):
        self._ion_conc = float(val)

    @property
    def protein(self):
        return self._protein

    @protein.setter
    def protein(self, val):
        # set protein ID and directory name
        self._protein = Ensemble._checkProtein(val)
        self._subdir = _subdir.Subdir(val, overwrite=True)

        with self._subdir:
            #download protein PDB
            self._protein_file = _IO.PDB.downloadPDB(id=self.protein)
            self._protein_obj = _IO.PDB.PDB(self._protein_file)

            # download SDF files and split them
            self._downloader = _pdbconnect.PDBDownloader(id=self.protein)
            self._ligand_files_pdb = self._downloader.downloadLigandSDF()
            self._ligand_files_pdb = _IO.SDF.splitSDFs(self._ligand_files_pdb)

        self._subdir.overwrite = False

        # remove ligands from PDB file and non-ligands from SDF files
        self.filterPDB(chains="all", waters="all", anions="all", cations="all", ligands="all")

    @property
    def ligand_id(self):
        return self._ligand_id

    @ligand_id.setter
    def ligand_id(self, id):
        for i, filename in enumerate(self._ligand_files_pdb):
            _, _, _, resSeq = _re.search(r"^([\w]*)_([\w]*)_([\w])_([A-Z0-9]*)", filename).groups()
            if resSeq == id:
                self._ligand_id = filename
                del self._ligand_files_pdb[i]
                return
        raise ValueError("Molecule ID not found: %s" % id)

    @property
    def morphs(self):
        return self._morphs

    @morphs.setter
    def morphs(self, vals):
        self._morphs = []
        if vals is None: return
        for val in vals:
            if len(val) != 2:
                _warnings.warn("Value %s not recognised as a container of ligand1 and ligand2. "
                               "Skipping value..." % val)
                continue
            ligand1, ligand2 = val

            try:
                ligand1 = _rdkit.openAsRdkit(ligand1)
                ligand2 = _rdkit.openAsRdkit(ligand2)
            except:
                _warnings.warn("One or more of these values is not valid: %s. Skipping values..." % val)
                continue

            if ligand1 not in self._ligands.keys():
                self._ligands[ligand1] = []
            if ligand2 not in self._ligands.keys():
                self._ligands[ligand2] = []

            self._morphs += [[ligand1, ligand2]]

    def filterPDB(self, chains="all", waters="site", ligands="chain", anions="chain", cations="chain",
                  include_mols=None, exclude_mols=None):
        """
        :type chains: "all" OR a list of characters for chains to keep
        :param waters: which waters to keep
        :type waters: "all" OR "chain" (only waters belonging to a chain) OR "site" OR None;
                      other molecules can be added via include_mols
        :param ligands: which ligands / cofactors to keep
        :type ligands: "all" OR "chain" OR "site" OR None
        :param anions: which anions to keep
        :type anions: "all" OR "chain" OR "site" OR None
        :param cations: which cations to keep
        :type cations: "all" OR "chain" OR "site" OR None
        :param include_mols: include residue name; overrides previous filters
        :type include_mols: list of strings
        :param exclude_mols: exclude residue name; overrides previous filters
        :type exclude_mols: list of strings
        :return: nothing; rewrites PDB file
        :rtype: void
        """
        if include_mols is None: include_mols = []
        if exclude_mols is None: exclude_mols = []

        with self._subdir:
            # filter ligands
            ligand_files = []
            for filename in self._ligand_files_pdb:
                _, resname, chainID, resSeq = _re.search(r"^([\w]*)_([\w]*)_([\w])_([A-Z0-9]*)", filename).groups()
                if _IO.PDB.Residue.classify(resname) != "ligand":
                    continue
                if ligands == "all" or (ligands == "chain" and chains == "all"):
                    ligand_files += [filename]
                    continue
                elif ligands == "chain" and chainID in chains:
                    ligand_files += [filename]
                    continue
                elif resSeq in include_mols:
                    ligand_files += [filename]
                    continue
                elif resSeq in exclude_mols:
                    continue
                elif ligands == "site":
                    resSeq, iCode = Ensemble._residTransform(resSeq)
                    for residue in self._protein_obj._site_residues:
                        if residue._resSeq == resSeq and residue._iCode == iCode:
                            ligand_files += [filename]
            self._ligand_files_pdb = ligand_files

            # filter residues / molecules in protein

            # filter by chain
            filter = []
            if chains == "all":
                filter += self._protein_obj.filterResidues(includedict={"_type": ["amino acid"]})
            else:
                filter += self._protein_obj.filterResidues(includedict={"_type": ["amino acid"], "_chainID": chains})

            # filter by waters / anions / cations
            for parameter, name in zip([waters, anions, cations], ["water", "anion", "cation"]):
                if parameter == "all" or (parameter == "chain" and chains == "all"):
                    filter += self._protein_obj.filterResidues(includedict={"_type": [name]})
                elif parameter == "chain":
                    filter += self._protein_obj.filterResidues(includedict={"_type": [name], "_chainID": chains},
                                                      includejoint=False)
                elif parameter == "site":
                    if chains == "all":
                        filter_temp = self._protein_obj.filterResidues(includedict={"_type": [name]})
                    else:
                        filter_temp = self._protein_obj.filterResidues(includedict={"_type": [name], "_chainID": chains},
                                                              includejoint=False)
                    filter += [res for res in filter_temp if res in self._protein_obj._site_residues]

            # include extra molecules / residues
            for include_mol in include_mols:
                resSeq, iCode = Ensemble._residTransform(include_mol)
                residue = self._protein_obj.filterResidues(includedict={"_resSeq": [resSeq], "_iCode": iCode})[0]
                if residue._type != "ligand":
                    filter += [residue]

            # exclude extra molecules / residues
            excl_filter = []

            for exclude_mol in exclude_mols:
                resSeq, iCode = Ensemble._residTransform(exclude_mol)
                residue = self._protein_obj.filterResidues(includedict={"_resSeq": [resSeq], "_iCode": iCode})[0]
                if residue._type != "ligand":
                    excl_filter += [residue]

            filter = list(set(filter) - set(excl_filter))
            self._protein_obj.purgeResidues(filter)
            self._protein_obj.writePDB(self._protein_file)

    def preparePDB(self, add_missing_residues="modeller", add_missing_atoms="pdb2pqr", protonate_proteins="pdb2pqr",
                   protonate_ligands="babel"):
        with self._subdir:
            add_missing_residues = add_missing_residues.strip().lower()
            add_missing_atoms = add_missing_atoms.strip().lower()
            protonate_proteins = protonate_proteins.strip().lower()
            protonate_ligands = protonate_ligands.strip().lower()

            if protonate_proteins != "protoss" and protonate_ligands == "protoss":
                _warnings.warn("Warning: Protoss cannot be individually called on a ligand. Changing protein protonation "
                               "method to Protoss...")
                protonate_proteins = "protoss"
            if add_missing_atoms == "pdb2pqr" and protonate_proteins != "pdb2pqr":
                print("Warning: Cannot run PDB2PQR without protonation. This will be fixed in a later version. "
                      "Changing protein protonation method to PDB2PQR...")

            if add_missing_residues.lower() == "modeller":
                filename_fasta = self._downloader.downloadFASTA()
                self._protein_file = _modeller.modellerTransform(self._protein_file, filename_fasta)

            if add_missing_atoms.lower() == "pdb2pqr":
                self._protein_file = _PDB2PQR.PDB2PQRtransform(self._protein_file)

            if protonate_proteins.lower() == "protoss":
                if protonate_ligands.lower() != "protoss":
                    ligand_file = None
                else:
                    ligand_file = _IO.SDF.mergeSDFs(self._ligand_files_pdb)
                self._protein_file, ligand_file, _ = _protoss.protossTransform(self._protein_file, ligand_file)
                if protonate_ligands.lower() == "protoss":
                    self._ligand_files_pdb = _IO.SDF.splitSDFs([ligand_file])

            if protonate_ligands.lower() == "babel":
                self._ligand_files_pdb = [_babel.babelTransform(f) for f in self._ligand_files_pdb if f is not None]

    def parametrisePDB(self):
        with self._subdir:
            print("Parametrising original crystal system...")
            #extract non-protein residues from pdb file and save them as separate pdb files
            hetatm_files, hetatm_types = self._protein_obj.writeHetatms()
            non_protein_residues = self._protein_obj.filterResidues(excludedict={"_type": ["water", "cation", "anion"]})
            self._protein_obj.purgeResidues(non_protein_residues)
            self._protein_obj.writePDB()

            # create a merged parmed object with all parametrised molecules and save to top/gro which is then read by
            # BioSimSpace
            # we can't use a direct BioSimSpace object because it throws an error if the file contains single ions
            system = _parametrise.parametriseAndLoadPmd(params=self.params, input_filename=self._protein_file,
                                                        molecule_type="protein")

            self._ligand_id = _babel.babelTransform(self._ligand_id, "mol2")
            self._ligand_ref = _rdkit.openAsRdkit(self._ligand_id, removeHs=False)

            for filename, type in zip(self._ligand_files_pdb + hetatm_files,
                                      ["ligand"] * len(self._ligand_files_pdb) + hetatm_types):
                system += _parametrise.parametriseAndLoadPmd(params=self.params, input_filename=filename,
                                                             molecule_type=type)

            if self.centre:
                coords = system.coordinates
                c_min, c_max = _np.amin(coords, axis=0), _np.amax(coords, axis=0)
                centre = (c_min + c_max) / 2
                box_length = max(c_max - c_min) / 10
                if box_length > self.box_length:
                    self.box_length = int(box_length) + 1
                    print("Insufficient input box size. Changing to a box length of %d nm..." % self.box_length)
                    centre -= _np.asarray(3 * [5 * box_length])
                for atom in system.atoms:
                    atom.xx -= centre[0]
                    atom.xy -= centre[1]
                    atom.xz -= centre[2]

                _rdkit.translateMolecule(self._ligand_ref, -centre)
                _rdkit.saveFromRdkit(self._ligand_ref, self._ligand_id)

            system.box = [10 * self.box_length, 10 * self.box_length, 10 * self.box_length, 90, 90, 90]
            system.save("complex_template.top")
            system.save("complex_template.gro")
            self._complex_template = _BSS.IO.readMolecules(["complex_template.top", "complex_template.gro"])

    def parametriseLigands(self):
        print("Parametrising all ligands. This will take a while...")
        with self._subdir:
            for i, ligand in enumerate(self._ligands.keys()):
                filename_temp = _rdkit.saveFromRdkit(ligand, filename="ligand_%d.mol" % (i + 1))
                filename_babel = _babel.babelTransform(filename_temp, "mol2")
                self._ligands[ligand] += [filename_babel]
                self._ligands[ligand] += [_parametrise.parametriseFile(params=self.params, input_filename=filename_babel,
                                                                       molecule_type="ligand")]

    def prepareComplexes(self):
        with self._subdir:
            for i, (ligand1, ligand2) in enumerate(self.morphs):
                with _subdir.Subdir("Pair %d" % (i + 1)) as curdir:
                    print("Creating morph %d..." % (i + 1))
                    #getting parametrised ligand filenames
                    ligand1_prep = ["../" + i for i in self._ligands[ligand1][1]]
                    ligand2_prep = ["../" + i for i in self._ligands[ligand2][1]]

                    #loading protonated ligands into rdkit
                    ligand1_H = _rdkit.openAsRdkit("../" + self._ligands[ligand1][0], removeHs=False)
                    ligand2_H = _rdkit.openAsRdkit("../" + self._ligands[ligand2][0], removeHs=False)

                    #aligning protonated ligands
                    ligand1_H = _rdkit.alignTwoMolecules(ligand1_H, self._ligand_ref)
                    ligand2_H = _rdkit.alignTwoMolecules(ligand2_H, ligand1_H)

                    #replacing coordinates of parametrised ligands with new coordinates while keeping the topology
                    filenames = ["morph%d_1.inpcrd" % (i + 1), "morph%d_2.inpcrd" % (i + 1)]
                    ligand1_prep[1] = _rdkit.saveFromRdkit(ligand1_H, filenames[0])
                    ligand2_prep[1] = _rdkit.saveFromRdkit(ligand2_H, filenames[1])

                    #loading the shifted ligands into BioSimSpace
                    ligand1_BSS = _BSS.IO.readMolecules(ligand1_prep).getMolecules()[0]
                    ligand2_BSS = _BSS.IO.readMolecules(ligand2_prep).getMolecules()[0]
                    mapping = _BSS.Align.matchAtoms(ligand1_BSS, ligand2_BSS)

                    #merging the two ligands into a morph and adding it to the protein template
                    morph = _BSS.Align.merge(ligand1_BSS, ligand2_BSS, mapping=mapping)
                    _BSS.IO.saveMolecules("morph", morph, "Gro87,GroTop")
                    complex = self._complex_template + morph

                    print("Solvating...")
                    complex = self._solvate(complex, work_dir=curdir.path)
                    self.complexes += [complex]

                    print("Saving solvated complex as GROMACS...")
                    _BSS.IO.saveMolecules("complex_final", complex, "Gro87,GroTop")
                    _pmd.load_file("complex_final.gro87").save("complex_final.gro")
                    _os.rename("complex_final.grotop", "complex_final.top")

    def _solvate(self, complex, work_dir=None):
        if work_dir is None:
            tmp_dir = _tempfile.TemporaryDirectory()
            work_dir = tmp_dir.name
        with _subdir.Subdir(dirname=work_dir):
            """CENTERING"""

            box_length = self.box_length
            if self.centre:
                box = complex._getAABox()
                min_coords, max_coords = box.minCoords(), box.maxCoords()
                centre = (min_coords + max_coords) / 2
                difference = max_coords - min_coords
                if any([x / 10 > self.box_length for x in difference]):
                    box_length = int(max(difference / 10)) + 1
                    print("Insufficient input box size. Changing to a box length of %d nm..." % box_length)
                    complex.translate(tuple(-centre + _SireMaths.Vector(3 * (box_length * 5, ))))
                    complex._sire_system.setProperty("space",
                                                     _SireVol.PeriodicBox(_SireMaths.Vector(3 * (box_length * 10,))))

            """SOLVATING"""

            _BSS.IO.saveMolecules("complex", complex, "Gro87,GroTop")
            files = ["complex.gro", "complex.top"]
            _os.rename("complex.gro87", files[0])
            _os.rename("complex.grotop", files[1])

            new_gro = "complex_solvated.gro"
            command = "{0} solvate -shell {1} -box {2} {2} {2} -cp {3} -o {4}".format(_const.GROMACSEXE, self.shell,
                                                                                      box_length, files[0], new_gro)

            try:
                """CREATE UNPARAMETRISED WATERS AND LOAD THEM INTO MEMORY"""

                _runexternal.runExternal(command, procname="gmx solvate")
                complex_solvated = _pmd.load_file(new_gro)
                waters = complex_solvated[":SOL"]

                """PREPARE WATERS FOR TLEAP AND PARAMETRISE"""

                for residue in waters.residues:
                    residue.name = "WAT"
                for i, atom in enumerate(waters.atoms):
                    if "H" in atom.name:
                        atom.name = atom.name[0] + atom.name[2]
                    elif "O" in atom.name:
                        atom.name = "O"

                waters.save("waters.pdb")
                waters_prep = _parametrise.parametriseAndLoadPmd(self.params, "waters.pdb", "water")
                waters_prep.box = _pmd.load_file(files[1], xyz=files[0], skip_bonds=True).box
                for residue in waters_prep.residues:
                    residue.name = "SOL"
                waters_prep_filenames = ["waters.gro", "waters.top"]
                for filename in waters_prep_filenames:
                    waters_prep.save(filename)

            except:
                print("Solvation failed. Returning original system...")
                return complex

            """ADD IONS"""

            if any([self.neutralise, self.ion_conc, self.shell]):
                """WRITE MDP FILE"""

                with open("ions.mdp", "w") as file:
                    file.write("; Neighbour searching\n")
                    file.write("cutoff-scheme           = Verlet\n")
                    file.write("rlist                   = 1.1\n")
                    file.write("pbc                     = xyz\n")
                    file.write("verlet-buffer-tolerance = -1\n")
                    file.write("\n; Electrostatics\n")
                    file.write("coulombtype             = cut-off\n")
                    file.write("\n; VdW\n")
                    file.write("rvdw                    = 1.0\n")

                n_Na, n_Cl = [int((box_length * 10**-8)** 3 * 6.022 * 10**23 * self.ion_conc)] * 2
                if self.neutralise:
                    charge = round(complex.charge().magnitude())
                    if charge > 0:
                        n_Na += charge
                    else:
                        n_Cl += charge

                try:
                    ions_prep_filenames = ["ions.gro", "ions.top"]
                    command = "{0} grompp -f ions.mdp -c {1} -p {2} -o complex_solvated.tpr".format(
                        _const.GROMACSEXE, *waters_prep_filenames)
                    _runexternal.runExternal(command, procname="gmx grompp")

                    command = "{{ echo 2; }} | {0} genion -s complex_solvated.tpr -o {1} -nn {2} -np {3}".format(
                        _const.GROMACSEXE, ions_prep_filenames[0], n_Cl, n_Na)
                    _runexternal.runExternal(command, procname="gmx genion")

                    _pmd.load_file(ions_prep_filenames[0])[":CL,NA"].save("ions.pdb")
                    _os.remove(ions_prep_filenames[0])
                    ions_prep = _parametrise.parametriseAndLoadPmd(self.params, "ions.pdb", "water")
                    for filename in ions_prep_filenames:
                        (ions_prep + waters_prep).save(filename)
                except:
                    print("Ion addition failed. Returning solvated system...")
                    return complex + _BSS.IO.readMolecules(waters_prep_filenames)
            return complex + _BSS.IO.readMolecules(ions_prep_filenames)

    @staticmethod
    def _checkProtein(val):
        val = val.strip().upper()
        if not val.isalnum():
            raise ValueError("Value %s not a valid PDB code. Skipping value..." % val)
        return val

    @staticmethod
    def _residTransform(id):
        id = id.strip()
        if id[-1].isalpha():
            iCode = id[-1]
            resSeq = int(float(id[:-1]))
        else:
            iCode = " "
            resSeq = int(float(id))
        return resSeq, iCode