import const as _const

import re as _re
import warnings as _warnings

import BioSimSpace as _BSS
from rdkit import Chem as _Chem
from rdkit.Chem import AllChem as _rdkitAllChem

import IO as _IO
import pdbconnect as _pdbconnect
import Parametrise as _parametrise
import subdir as _subdir
from Wrappers import babelwrapper as _babel
from Wrappers import modellerwrapper as _modeller
from Wrappers import PDB2PQRwrapper as _PDB2PQR
from Wrappers import protosswrapper as _protoss

class Ensemble:
    def __init__(self, engine, box_length=12, ion_conc=0.154, protein_ff=_const.AMBERDEFAULTPROTEINFF,
                 ligand_ff=_const.AMBERDEFAULTLIGANDFF, water_ff=_const.AMBERDEFAULTWATERFF, protein=None,
                 ligand_id=None, morphs=None):
        self.engine = engine
        self.box_length = box_length
        self.ion_conc = ion_conc
        self.params = _parametrise.Params(protein_ff=protein_ff, ligand_ff=ligand_ff, water_ff=water_ff)
        self.protein = protein
        self.ligand_id = ligand_id
        self._ligands = []
        self._ligand_prep = []
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
                ligand1 = Ensemble._checkLigand(ligand1)
                ligand2 = Ensemble._checkLigand(ligand2)
            except:
                _warnings.warn("One or more of these values is not valid: %s. Skipping values..." % val)
                continue

            if ligand1 not in self._ligands:
                self._ligands += [ligand1]
            if ligand2 not in self._ligands:
                self._ligands += [ligand2]

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

            for filename, type in zip(self._ligand_files_pdb + hetatm_files,
                                      ["ligand"] * len(self._ligand_files_pdb) + hetatm_types):
                system += _parametrise.parametriseAndLoadPmd(params=self.params, input_filename=filename,
                                                             molecule_type=type)

            system.save("complex_template.top")
            system.save("complex_template.gro")
            self._complex_template = _BSS.IO.readMolecules(["complex_template.top", "complex_template.gro"])

    def _parametriseLigands(self):
        with self._subdir:
            self._ligand_prep = []
            for i, ligand in enumerate(self._ligands):
                filename = "ligand_%d.sdf" % (i + 1)
                writer = _Chem.SDWriter(filename)
                writer.write(ligand)
                writer.close()
                filename = _babel.babelTransform(filename, "mol2")
                self._ligand_prep += [_parametrise.parametriseAndLoadBSS(params=self.params, input_filename=filename,
                                                                         molecule_type="ligand").getMolecules()[0]]

    def prepareComplexes(self):
        box = 3 * [self.box_length * _BSS.Units.Length.nanometer]
        print("Parametrising all ligands. This will take a while...")
        self._parametriseLigands()

        with self._subdir:
            print("Parametrising original crystal ligand...")
            ligand_ref = _babel.babelTransform(self.ligand_id, "mol2")
            ligand_ref = _BSS.IO.readMolecules(ligand_ref).getMolecules()[0]

            for i, (ligand1, ligand2) in enumerate(self.morphs):
                with _subdir.Subdir("Pair %d" % (i + 1)):
                    print("Creating morph %d..." % (i + 1))
                    ligand1_prep = self._ligand_prep[self._ligands.index(ligand1)]
                    ligand2_prep = self._ligand_prep[self._ligands.index(ligand2)]

                    ligand1_prep = _BSS.Align.rmsdAlign(ligand1_prep, ligand_ref)
                    morph = _BSS.Align.merge(ligand1_prep, ligand2_prep)

                    complex = self._complex_template + morph

                    print("Solvating...")
                    while True:
                        try:
                            complex = _BSS.Solvent.solvate(self.params.water_ff, molecule=complex, box=box,
                                                           ion_conc=self.ion_conc)
                            break
                        except ValueError:
                            print("Box size of %d nm is too small, attempting box size of %d nm" %
                                  (self.box_length, self.box_length + 1))
                            self.box_length += 1
                            box = 3 * [self.box_length * _BSS.Units.Length.nanometer]
                            continue

                    self.complexes += [complex]

                    print("Saving as GROMACS...")
                    _BSS.IO.saveMolecules("complex", complex, "Gro87")
                    _BSS.IO.saveMolecules("complex", complex, "GroTop")
                    print("Saving as AMBER...")
                    _BSS.IO.saveMolecules("complex", complex, "PRM7,RST7")

    @staticmethod
    def _checkProtein(val):
        val = val.strip().upper()
        if not val.isalnum():
            raise ValueError("Value %s not a valid PDB code. Skipping value..." % val)
        return val

    @staticmethod
    def _checkLigand(val):
        if len(val.split(".")) > 1:
            mol = _Chem.MolFromMolFile(val)
        else:
            try:
                mol = _Chem.MolFromSmiles(val)
            except:
                try:
                    mol = _Chem.MolFromInchi(val)
                except:
                    raise ValueError("Value %s not a valid file, SMILES string, or InChI identifier.")
            _rdkitAllChem.EmbedMolecule(mol, _rdkitAllChem.ETKDG())
        return mol

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