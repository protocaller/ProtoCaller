import os as _os
import re as _re
import tempfile as _tempfile
import warnings as _warnings

import BioSimSpace as _BSS
import numpy as _np
import parmed as _pmd
import Sire.Maths as _SireMaths
import Sire.Vol as _SireVol

import ProtoCaller as _PC
import ProtoCaller.IO as _IO
import ProtoCaller.Parametrise as _parametrise
import ProtoCaller.Utils.pdbconnect as _pdbconnect
import ProtoCaller.Utils.runexternal as _runexternal
import ProtoCaller.Utils.stdio as _stdio
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Wrappers.babelwrapper as _babel
import ProtoCaller.Wrappers.BioSimSpacewrapper as _BSSwrap
import ProtoCaller.Wrappers.modellerwrapper as _modeller
import ProtoCaller.Wrappers.PDB2PQRwrapper as _PDB2PQR
import ProtoCaller.Wrappers.protosswrapper as _protoss
import ProtoCaller.Wrappers.rdkitwrapper as _rdkit

class Ensemble:
    def __init__(self, engine, box_length=12, ion_conc=0.154, shell=0, neutralise=True, centre=True,
                 protein_ff=_PC.AMBERDEFAULTPROTEINFF, ligand_ff=_PC.AMBERDEFAULTLIGANDFF,
                 water_ff=_PC.AMBERDEFAULTWATERFF, protein=None, ligand_id=None, morphs=None):
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
        self.systems_prep = {}

    @property
    def engine(self):
        return self._engine

    @engine.setter
    def engine(self, val):
        val = val.strip().upper()
        if val not in _PC.ENGINES:
            raise ValueError("Value %s not supported. Supported values: " % val, _PC.ENGINES)
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
        self._workdir = _fileio.Subdir(val)

        with self._workdir:
            #download protein PDB
            self._protein_file = _IO.PDB.downloadPDB(id=self.protein)
            self._protein_obj = _IO.PDB.PDB(self._protein_file)

            # download SDF files and split them
            self._downloader = _pdbconnect.PDBDownloader(id=self.protein)
            self._ligand_files_pdb = self._downloader.downloadLigandSDF()
            self._ligand_files_pdb = _IO.SDF.splitSDFs(self._ligand_files_pdb)

        # remove ligands from PDB file and non-ligands from SDF files
        self.filterPDB(chains="all", waters="all", simple_anions="all", complex_anions="all", simple_cations="all",
                       complex_cations="all", ligands="all", cofactors="all")

    @property
    def ligand_id(self):
        return self._ligand_id

    @ligand_id.setter
    def ligand_id(self, id):
        if id is None and len(self._ligand_files_pdb) == 1:
            self._ligand_id = self._ligand_files_pdb[0]
            self._ligand_files_pdb = []
            return

        for i, filename in enumerate(self._ligand_files_pdb):
            _, _, _, resSeq = _re.search(r"^([\w]+)_([\w]+)_([\w])_([A-Z0-9]+)", filename).groups()
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
                mol1 = _rdkit.openAsRdkit(ligand1)
                mol2 = _rdkit.openAsRdkit(ligand2)
            except:
                _warnings.warn("One or more of these values is not valid: %s. Skipping values..." % val)
                continue

            if ligand1 not in self._ligands.keys():
                self._ligands[ligand1] = [mol1]
            if ligand2 not in self._ligands.keys():
                self._ligands[ligand2] = [mol2]

            self._morphs += [[ligand1, ligand2]]

    def filterPDB(self, missing_residues="middle", chains="all", waters="site", ligands="chain", cofactors="chain",
                  simple_anions="chain", complex_anions="chain", simple_cations="chain", complex_cations="chain",
                  include_mols=None, exclude_mols=None):
        """
        :type missing_residues: "all" or "middle"
        :param missing_residues: which missing residues to keep
        :type chains: "all" OR a list of characters for chains to keep
        :type everything inbetween: "all" OR "chain" (only molecules belonging to a chain) OR "site" OR None;
                                    other molecules can be added via include_mols
        :param include_mols: include residue name; overrides previous filters
        :type include_mols: list of strings
        :param exclude_mols: exclude residue name; overrides previous filters
        :type exclude_mols: list of strings
        :return: nothing; rewrites PDB file
        :rtype: void
        """
        if include_mols is None: include_mols = []
        if exclude_mols is None: exclude_mols = []

        with self._workdir:
            # filter ligands
            ligand_files = []
            for filename in self._ligand_files_pdb:
                _, resname, chainID, resSeq = _re.search(r"^([\w]+)_([\w]+)_([\w])_([A-Z0-9]+)", filename).groups()
                if _PC.RESIDUETYPE(resname) != "ligand":
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
            mask = "_type=='amino acid'"
            if chains != "all":
                mask += "&_chainID in %s" % str(chains)
            filter += self._protein_obj.filter(mask)

            # filter missing residues
            if missing_residues == "middle":
                missing_residue_list = self._protein_obj.totalResidueList()
                for i in range(2):
                    missing_residue_list.reverse()
                    current_chain = None
                    for j in reversed(range(0, len(missing_residue_list))):
                        res = missing_residue_list[j]
                        if type(res) == _IO.PDB.MissingResidue and current_chain != res._chainID:
                            del missing_residue_list[j]
                        else:
                            current_chain = res._chainID
                missing_residue_list = [x for x in missing_residue_list if type(x) == _IO.PDB.MissingResidue]
            filter += missing_residue_list

            # filter by waters / anions / cations
            for param, name in zip([waters, simple_anions, complex_anions, simple_cations, complex_cations, cofactors],
                ["water", "simple_anion", "complex_anion", "simple_cation", "complex_cation", "cofactor"]):
                if param == "all" or (param == "chain" and chains == "all"):
                    filter += self._protein_obj.filter("_type=='%s'" % name)
                elif param == "chain":
                    filter += self._protein_obj.filter("_type=='%s'|_chainID in %s" % (name, str(chains)))
                elif param == "site":
                    if chains == "all":
                        filter_temp = self._protein_obj.filter("_type=='%s'" % name)
                    else:
                        filter_temp = self._protein_obj.filter("_type=='%s'|_chainID in %s" % (name, str(chains)))
                    filter += [res for res in filter_temp if res in self._protein_obj._site_residues]

            # include extra molecules / residues
            for include_mol in include_mols:
                residue = self._protein_obj.filter(include_mol + "&_type!='ligand'")
                filter += [residue]

            # exclude extra molecules / residues
            excl_filter = []
            for exclude_mol in exclude_mols:
                residue = self._protein_obj.filter(exclude_mol + "&_type!='ligand'")
                excl_filter += [residue]

            filter = list(set(filter) - set(excl_filter))
            self._protein_obj.purgeResidues(filter)
            self._protein_obj.writePDB(self._protein_file)

    def loadParametrisedLigands(self, ligand_dict):
        for ligand, filearr in ligand_dict.items():
            if len(filearr) != 2:
                _warnings.warn("Need to pass a protonated structure file and a list of parametrised files. "
                               "Skipping value...")
            else:
                try:
                    mol = _rdkit.openAsRdkit(ligand)
                    self._ligands[ligand] = [mol, *filearr]
                    self._ligands[ligand][1] = _os.path.abspath(filearr[0])
                    self._ligands[ligand][2] = [_os.path.abspath(x) for x in filearr[1]]

                    if any([not _os.path.exists(f) for f in [self._ligands[ligand][1], self._ligands[ligand][2][0],
                                                             self._ligands[ligand][2][1]]]):
                        raise ValueError()
                except:
                    _warnings.warn("An error occured. Skipping ligand...")
                    self._ligands[ligand] = [mol]

    def preparePDB(self, add_missing_residues="modeller", add_missing_atoms="pdb2pqr", protonate_proteins="pdb2pqr",
                   protonate_ligands="babel"):
        with self._workdir:
            add_missing_residues = add_missing_residues.strip().lower() if add_missing_residues is not None else ""
            add_missing_atoms = add_missing_atoms.strip().lower() if add_missing_atoms is not None else ""
            protonate_proteins = protonate_proteins.strip().lower() if protonate_proteins is not None else ""
            protonate_ligands = protonate_ligands.strip().lower() if protonate_ligands is not None else ""

            if protonate_proteins != "protoss" and protonate_ligands == "protoss":
                _warnings.warn("Warning: Protoss cannot be individually called on a ligand. Changing protein protonation "
                               "method to Protoss...")
                protonate_proteins = "protoss"
            if add_missing_atoms == "pdb2pqr" and protonate_proteins != "pdb2pqr":
                print("Warning: Cannot run PDB2PQR without protonation. This will be fixed in a later version. "
                      "Changing protein protonation method to PDB2PQR...")

            if len(self._protein_obj._missing_residues) and not add_missing_atoms == "modeller":
                if add_missing_residues == "modeller":
                    atoms = True if add_missing_atoms == "modeller" else False
                    filename_fasta = self._downloader.downloadFASTA()
                    self._protein_file = _modeller.modellerTransform(self._protein_file, filename_fasta, atoms)
                    self._protein_obj = _IO.PDB.PDB(self._protein_file)

            if add_missing_atoms == "pdb2pqr":
                self._protein_file = _PDB2PQR.PDB2PQRtransform(self._protein_file)
                self._protein_obj = _IO.PDB.PDB(self._protein_file)

            if protonate_proteins == "protoss":
                if protonate_ligands != "protoss":
                    ligand_file = None
                else:
                    ligand_file = _IO.SDF.mergeSDFs(self._ligand_files_pdb)
                self._protein_file, ligand_file, _ = _protoss.protossTransform(self._protein_file, ligand_file)
                if protonate_ligands == "protoss":
                    self._ligand_files_pdb = _IO.SDF.splitSDFs([ligand_file])

            if protonate_ligands == "babel":
                self._ligand_files_pdb = [_babel.babelTransform(f) for f in self._ligand_files_pdb if f is not None]

    def parametrisePDB(self):
        with self._workdir:
            print("Parametrising original crystal system...")
            #extract non-protein residues from pdb file and save them as separate pdb files
            hetatm_files, hetatm_types = self._protein_obj.writeHetatms()
            non_protein_residues = self._protein_obj.filter("_type not in ['water', 'simple cation', 'simple anion']")
            self._protein_obj.purgeResidues(non_protein_residues)
            self._protein_obj.writePDB()

            # create a merged parmed object with all parametrised molecules and save to top/gro which is then read by
            # BioSimSpace
            # we can't use a direct BioSimSpace object because it throws an error if the file contains single ions
            system = _parametrise.parametriseAndLoadPmd(params=self.params, filename=self._protein_file,
                molecule_type="protein", disulfide_bonds=self._protein_obj._disulfide_bonds)

            self._ligand_id = _babel.babelTransform(self._ligand_id, "mol2")
            self._ligand_ref = _rdkit.openAsRdkit(self._ligand_id, removeHs=False)

            for filename, type in zip(self._ligand_files_pdb + hetatm_files,
                                      ["ligand"] * len(self._ligand_files_pdb) + hetatm_types):
                system += _parametrise.parametriseAndLoadPmd(params=self.params, filename=filename,
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
            _IO.GROMACS.saveAsGromacs("complex_template", system)
            self._complex_template = _BSS.IO.readMolecules(["complex_template.top", "complex_template.gro"])

    def parametriseLigands(self):
        print("Parametrising all ligands. This will take a while...")
        with self._workdir:
            for i, ligand in enumerate(self._ligands.keys()):
                if len(self._ligands[ligand]) == 3: continue
                filename_temp = _rdkit.saveFromRdkit(self._ligands[ligand][0], filename="ligand_%d.mol" % (i + 1))
                filename_babel = _babel.babelTransform(filename_temp, "mol2")
                self._ligands[ligand] += [filename_babel]
                self._ligands[ligand] += [_parametrise.parametriseFile(params=self.params, filename=filename_babel,
                                                                       molecule_type="ligand", id="ligand_%d" % (i + 1))]

    def prepareComplexes(self, replica_temps=None, intermediate_files=False):
        if replica_temps is None or len(replica_temps) == 1:
            scales = [1]
        else:
            sorted_list = sorted(replica_temps)
            if sorted_list != replica_temps:
                _warnings.warn("Input replica temperatures were not in ascending order. Sorting...")
                replica_temps = sorted_list
            scales = [replica_temps[0] / elem for elem in replica_temps]

        with self._workdir:
            for i, (ligand1, ligand2) in enumerate(self.morphs):
                if intermediate_files:
                    curdir =_fileio.Dir("Pair %d" % (i + 1), overwrite=True)
                else:
                    curdir = _fileio.Dir(_tempfile.TemporaryDirectory().name, overwrite=True)
                with curdir:
                    print("Creating morph %d..." % (i + 1))
                    #getting parametrised ligand filenames
                    ligand1_prep = self._ligands[ligand1][2]
                    ligand2_prep = self._ligands[ligand2][2]

                    #loading protonated ligands into rdkit
                    ligand1_H = _rdkit.openAsRdkit(self._ligands[ligand1][1], removeHs=False)
                    ligand2_H = _rdkit.openAsRdkit(self._ligands[ligand2][1], removeHs=False)

                    #aligning protonated ligands
                    ligand1_H, _ = _rdkit.alignTwoMolecules(self._ligand_ref, ligand1_H)
                    ligand2_H, mcs = _rdkit.alignTwoMolecules(ligand1_H, ligand2_H)

                    #replacing coordinates of parametrised ligands with new coordinates while keeping the topology
                    filenames = ["morph%d_1.inpcrd" % (i + 1), "morph%d_2.inpcrd" % (i + 1)]
                    ligand1_prep[1] = _os.path.abspath(_rdkit.saveFromRdkit(ligand1_H, filenames[0]))
                    ligand2_prep[1] = _os.path.abspath(_rdkit.saveFromRdkit(ligand2_H, filenames[1]))

                    #loading the shifted ligands into BioSimSpace
                    ligand1_BSS = _BSS.IO.readMolecules(ligand1_prep).getMolecules()[0]
                    ligand2_BSS = _BSS.IO.readMolecules(ligand2_prep).getMolecules()[0]

                    #translating RDKit mapping into BioSimSpace mapping
                    mapping = {}
                    indices_1 = [atom.index() for atom in ligand1_BSS._sire_molecule.edit().atoms()]
                    indices_2 = [atom.index() for atom in ligand2_BSS._sire_molecule.edit().atoms()]
                    for idx1, idx2 in mcs:
                        mapping[indices_1[idx1]] = indices_2[idx2]

                    #merging the two ligands into a morph and adding it to the protein template
                    morph = _BSS.Align.merge(ligand1_BSS, ligand2_BSS, mapping=mapping)
                    complexes = [self._complex_template + morph]

                    #solvating and saving the prepared complex and morph with the appropriate box size
                    print("Solvating...")
                    complexes = [self._solvate(complexes[0], work_dir=curdir.path)]
                    box = complexes[0]._sire_system.property("space")
                    morph = morph.toSystem()
                    morph._sire_system.setProperty("space", box)

                    #rescaling complexes for replica exchange
                    if len(scales) != 1:
                        print("Creating replicas...")
                        complexes = [_BSSwrap.rescaleSystemParams(complexes[0], scale, includelist=["Merged_Molecule"])
                                     for scale in scales]

                    self.systems_prep[tuple(self.morphs[i])] = (morph, complexes)

    def saveSystems(self):
        print("Saving solvated complexes as GROMACS...")
        with self._workdir:
            for i, pair in enumerate(self.morphs):
                with _fileio.Subdir("Pair %d" % (i + 1)):
                    try:
                        morph, complexes = self.systems_prep[tuple(pair)]
                    except:
                        continue
                    _IO.GROMACS.saveAsGromacs("morph", morph)

                    if len(complexes) == 1:
                        _IO.GROMACS.saveAsGromacs("complex_final", complexes[0])
                    else:
                        for j, complex in enumerate(complexes):
                            _IO.GROMACS.saveAsGromacs("complex_final%d" % j, complex)

    def _solvate(self, complex, work_dir=None):
        if work_dir is None:
            tmp_dir = _tempfile.TemporaryDirectory()
            work_dir = tmp_dir.name
        with _fileio.Dir(dirname=work_dir):
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
            command = "{0} solvate -shell {1} -box {2} {2} {2} -cp {3} -o {4}".format(_PC.GROMACSEXE, self.shell,
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
                waters_prep.box = _stdio.ignore_warnings(_pmd.load_file)(files[1], xyz=files[0], skip_bonds=True).box
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

                charge = round(complex.charge().magnitude()) if self.neutralise else 0
                n_Na, n_Cl = [int((box_length * 10**-8)**3 * 6.022 * 10**23 * self.ion_conc) - abs(charge) // 2] * 2
                if self.neutralise:
                    if charge < 0:
                        n_Na -= charge
                    else:
                        n_Cl += charge

                try:
                    ions_prep_filenames = ["ions.gro", "ions.top"]
                    command = "{0} grompp -f ions.mdp -c {1} -p {2} -o complex_solvated.tpr".format(
                        _PC.GROMACSEXE, *waters_prep_filenames)
                    _runexternal.runExternal(command, procname="gmx grompp")

                    command = "{{ echo 2; }} | {0} genion -s complex_solvated.tpr -o {1} -nn {2} -np {3}".format(
                        _PC.GROMACSEXE, ions_prep_filenames[0], n_Cl, n_Na)
                    _runexternal.runExternal(command, procname="gmx genion")

                    """PREPARE WATERS FOR TLEAP AND PARAMETRISE"""
                    ions = _pmd.load_file(ions_prep_filenames[0])
                    for residue in ions.residues:
                        if residue.name == "SOL":
                            residue.name = "WAT"

                    ions.save("ions.pdb")
                    _os.remove(ions_prep_filenames[0])
                    ions_prep = _parametrise.parametriseAndLoadPmd(self.params, "ions.pdb", "water")

                    for filename in ions_prep_filenames:
                        ions_prep.save(filename)
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
