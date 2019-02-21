import re as _re
import warnings as _warnings

import ProtoCaller as _PC

if _PC.BIOSIMSPACE:
    import BioSimSpace as _BSS

from ProtoCaller.Ensemble import Ligand
import ProtoCaller.IO.PDB as _PDB
import ProtoCaller.IO.SDF as _SDF
import ProtoCaller.IO.GROMACS as _GROMACS
import ProtoCaller.Parametrise as _parametrise
import ProtoCaller.Utils.fileio as _fileio
import ProtoCaller.Utils.pdbconnect as _pdbconnect
if _PC.MODELLER:
    import ProtoCaller.Wrappers.modellerwrapper as _modeller
import ProtoCaller.Wrappers.parmedwrapper as _pmdwrap
import ProtoCaller.Wrappers.PDB2PQRwrapper as _PDB2PQR


class Protein:
    def __init__(self, code, pdb_file=None, ligand_ref=None, fasta_file=None, complex_template=None, name=None,
                 workdir=None):
        self.workdir = _fileio.Dir(workdir) if workdir else _fileio.Dir(code)
        with self.workdir:
            self.name = name if name is not None else code
            self._downloader = _pdbconnect.PDBDownloader(code)
            self.pdb = pdb_file
            self.fasta = fasta_file
            self.complex_template = complex_template
            ligands_original = _SDF.splitSDFs(self._downloader.getLigands())
            self.ligands = [Ligand(x, name=x, workdir=".", minimise=False) for x in ligands_original]
            self.cofactors = []
            # remove ligands from PDB file and non-ligands from SDF files
            self.filter(chains="all", waters="all", simple_anions="all", complex_anions="all", simple_cations="all",
                        complex_cations="all", ligands="all", cofactors="all")
            self.ligand_ref = ligand_ref

    @property
    def ligand_ref(self):
        return self._ligand_ref

    @ligand_ref.setter
    def ligand_ref(self, input):
        if hasattr(self, "_ligand_ref") and self._ligand_ref is not None:
            raise ValueError("Cannot reassign an already assigned reference ligand. Please create a new class instance")
        elif input is None and len(self.ligands) == 1:
            self._ligand_ref = self.ligands[0]
            self.ligands = []
        elif input is False:
            self._ligand_ref = None
        elif isinstance(input, Ligand):
            self._ligand_ref = input
        elif isinstance(input, str):
            for i, ligand in enumerate(self.ligands):
                _, _, _, resSeq = _re.search(r"^([\w]+)_([\w]+)_([\w])_([A-Z0-9]+)", ligand.name).groups()
                if resSeq == input:
                    self._ligand_ref = ligand
                    del self.ligands[i]
                    return
            raise ValueError("Molecule ID not found: %s" % input)
        else:
            try:
                self._ligand_ref = Ligand(input, name=self.name + "_ref", workdir=".", minimise=False)
            except:
                raise TypeError("Unrecognised type of input. Need either a Ligand or a PDB ID")

    @property
    def pdb(self):
        with self.workdir:
            if self._pdb is None:
                return self._downloader.getPDB()
            return self._pdb

    @pdb.setter
    def pdb(self, value):
        if value is not None:
            value = _fileio.checkFileExists(value)
        self._pdb = value
        self._pdb_obj = _PDB.PDB(self.pdb)

    @property
    def pdb_obj(self):
        return self._pdb_obj

    @property
    def fasta(self):
        with self.workdir:
            if self._fasta is None:
                return self._downloader.getFASTA()
            return self._fasta

    @fasta.setter
    def fasta(self, value):
        if value is not None:
            value = _fileio.checkFileExists(value)
        self._fasta = value

    def filter(self, missing_residues="middle", chains="all", waters="site", ligands="chain", cofactors="chain",
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

        with self.workdir:
            # filter ligands / cofactors
            temp_dict = {"ligand": [], "cofactor": []}
            for ligand in self.ligands:
                filename = ligand.name
                for param, name in zip([ligands, cofactors], ["ligand", "cofactor"]):
                    _, resname, chainID, resSeq = _re.search(r"^([\w]+)_([\w]+)_([\w])_([A-Z0-9]+)", filename).groups()
                    if _PC.RESIDUETYPE(resname) != name or resSeq in exclude_mols:
                        continue
                    if any([param == "all", param == "chain" and chains == "all",
                            param == "chain" and chainID in chains, resSeq in include_mols]):
                        temp_dict[name] += [filename]
                    elif param == "site":
                        resSeq, iCode = self._residTransform(resSeq)
                        for residue in self._pdb_obj.site_residues:
                            if residue.resSeq == resSeq and residue.iCode == iCode:
                                temp_dict[name] += [filename]
            self.ligands = [Ligand(x, name=x, workdir=".", minimise=False) for x in temp_dict["ligand"]]
            self.cofactors = [Ligand(x, name=x, workdir=".", minimise=False) for x in temp_dict["cofactor"]]

            # filter residues / molecules in protein

            # filter by chain
            filter = []
            mask = "type in ['amino_acid', 'amino_acid_modified']"
            if chains != "all":
                mask += "&chainID in %s" % str(chains)
            filter += self._pdb_obj.filter(mask)

            # filter missing residues
            if missing_residues == "middle":
                missing_residue_list = self._pdb_obj.totalResidueList()
                for i in range(2):
                    missing_residue_list.reverse()
                    current_chain = None
                    for j in reversed(range(0, len(missing_residue_list))):
                        res = missing_residue_list[j]
                        if type(res) == _PDB.Missing.MissingResidue and current_chain != res.chainID:
                            del missing_residue_list[j]
                        else:
                            current_chain = res.chainID
                missing_residue_list = [x for x in missing_residue_list if type(x) == _PDB.MissingResidue]
                filter += missing_residue_list

            # filter by waters / anions / cations
            for param, name in zip([waters, simple_anions, complex_anions, simple_cations, complex_cations],
                                   ["water", "simple_anion", "complex_anion", "simple_cation", "complex_cation"]):
                if param == "all" or (param == "chain" and chains == "all"):
                    filter += self._pdb_obj.filter("type=='%s'" % name)
                elif param == "chain":
                    filter += self._pdb_obj.filter("type=='%s'|chainID in %s" % (name, str(chains)))
                elif param == "site":
                    if chains == "all":
                        filter_temp = self._pdb_obj.filter("type=='%s'" % name)
                    else:
                        filter_temp = self._pdb_obj.filter("type=='%s'|chainID in %s" % (name, str(chains)))
                    filter += [res for res in filter_temp if res in self._pdb_obj.site_residues]

            # include extra molecules / residues
            for include_mol in include_mols:
                residue = self._pdb_obj.filter(include_mol + "&type!='ligand'")
                filter += residue

            # exclude extra molecules / residues
            excl_filter = []
            for exclude_mol in exclude_mols:
                residue = self._pdb_obj.filter(exclude_mol + "&type!='ligand'")
                excl_filter += residue

            filter = list(set(filter) - set(excl_filter))
            self._pdb_obj.purgeResidues(filter, "keep")
            self._pdb_obj.writePDB(self.pdb)

    def prepare(self, add_missing_residues="modeller", add_missing_atoms="pdb2pqr", protonate_proteins="pdb2pqr",
                protonate_ligands="babel"):
        with self.workdir:
            add_missing_residues = add_missing_residues.strip().lower() if add_missing_residues is not None else ""
            add_missing_atoms = add_missing_atoms.strip().lower() if add_missing_atoms is not None else ""
            protonate_proteins = protonate_proteins.strip().lower() if protonate_proteins is not None else ""
            protonate_ligands = protonate_ligands.strip().lower() if protonate_ligands is not None else ""

            if protonate_proteins != "protoss" and protonate_ligands == "protoss":
                _warnings.warn("Protoss cannot be individually called on a ligand. Changing protein protonation "
                               "method to Protoss...")
                protonate_proteins = "protoss"
            if add_missing_atoms == "pdb2pqr" and protonate_proteins != "pdb2pqr":
                _warnings.warn("Cannot currently run PDB2PQR without protonation. "
                               "This will be fixed in a later version. "
                               "Changing protein protonation method to PDB2PQR...")

            if len(self._pdb_obj.missing_residues):
                if add_missing_residues == "modeller" and _PC.MODELLER:
                    atoms = True if add_missing_atoms == "modeller" else False
                    filename_fasta = self._downloader.getFASTA()
                    self.pdb = _modeller.modellerTransform(self.pdb, filename_fasta, atoms)
                else:
                    _warnings.warn("Protein has missing residues. Please check your PDB file or choose a valid "
                                   "automatisation protocol")

            if add_missing_atoms == "pdb2pqr":
                self.pdb = _PDB2PQR.PDB2PQRtransform(self.pdb)

            if protonate_proteins == "protoss":
                # TODO: fix
                raise ValueError("Protoss support is not yet fully functional")
                # this code is obsolete and needs to be fixed
                # ligand_file = _SDF.mergeSDFs(self.ligands) if protonate_ligands == "protoss" else None
                # self.pdb, ligand_file, _ = _protoss.protossTransform(self.pdb, ligand_file)
                # if protonate_ligands == "protoss":
                #     self.ligands = _SDF.splitSDFs([ligand_file])

            if len(self.ligands + self.cofactors):
                if protonate_ligands == "babel":
                    for ligand in self.ligands + self.cofactors:
                        if not ligand.protonated:
                            ligand.protonate()
                else:
                    _warnings.warn("Need to protonate all relevant ligands / cofactors before any parametrisation")

    def parametrise(self, params, reparametrise=False):
        if self.complex_template is not None and not reparametrise:
            print("Protein complex template %s is already parametrised." % self.name)
            return

        with self.workdir:
            print("Parametrising original crystal system...")
            # extract non-protein residues from pdb file and save them as separate pdb files
            hetatm_files, hetatm_types = self._pdb_obj.writeHetatms()
            non_protein_residues = self._pdb_obj.filter("type not in ['water', 'simple_cation', 'simple_anion']")
            self._pdb_obj.purgeResidues(non_protein_residues, "keep")
            self._pdb_obj.writePDB(self.pdb)

            # create a merged parmed object with all parametrised molecules and save to top/gro which is then read by
            # BioSimSpace
            # we can't use a direct BioSimSpace object because it throws an error if the file contains single ions
            system = _parametrise.parametriseAndLoadPmd(params=params, filename=self.pdb, molecule_type="protein",
                                                        disulfide_bonds=self._pdb_obj.disulfide_bonds)

            # add ligands to the system
            for ligand in self.ligands:
                ligand.parametrise(params, reparametrise=reparametrise, molecule_type="ligand")
                system += _pmdwrap.openFilesAsParmed(ligand.parametrised_files)

            # add cofactors to the system
            for cofactor in self.cofactors:
                cofactor.parametrise(params, reparametrise=reparametrise, molecule_type="cofactor")
                system += _pmdwrap.openFilesAsParmed(cofactor.parametrised_files)

            # add all other HETATMS to the system
            for filename, type in zip(hetatm_files, hetatm_types):
                system += _parametrise.parametriseAndLoadPmd(params=params, filename=filename, molecule_type=type)

            _GROMACS.saveAsGromacs("complex_template", system)
            if _PC.BIOSIMSPACE:
                self.complex_template = _BSS.IO.readMolecules(["complex_template.top", "complex_template.gro"])
            else:
                self.complex_template = system

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
