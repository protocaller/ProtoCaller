import re as _re
import warnings as _warnings

from Bio import SeqIO as _SeqIO

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
import ProtoCaller.Wrappers.charmmguiwrapper as _charmmwrap
if _PC.MODELLER:
    import ProtoCaller.Wrappers.modellerwrapper as _modeller
import ProtoCaller.Wrappers.parmedwrapper as _pmdwrap
import ProtoCaller.Wrappers.pdbfixerwrapper as _pdbfix
import ProtoCaller.Wrappers.pdb2pqrwrapper as _PDB2PQR


class Protein:
    """
    The main class responsible for handling proteins in ProtoCaller.

    Parameters
    ----------
    code : str
        The PDB code of the protein.
    pdb_file : str, optional
        Path to custom PDB file.
    ligand_ref : str or ProtoCaller.Ensemble.Ligand.Ligand or None or False
        Initialises ligand_ref from a "resSeq/iCode" (e.g. "400G") or from a Ligand or from a file or automatically
        (None) or enforces no reference ligand (False)
    fasta_file : str, optional
        Path to custom FASTA file.
    complex_template : BioSimSpace.System
        Initialises the complex from an already created complex_template
    name  optional : str
        Initialises the protein name.
    workdir : str, optional
        Initialises workdir.

    Attributes
    ----------
    workdir : ProtoCaller.Utils.fileio.Dir
        The working directory of the protein.
    name : str
        The protein name. Default is the PDB code thereof.
    """
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
            self.filter(missing_residues="all", chains="all", waters="all",
                        simple_anions="all", complex_anions="all",
                        simple_cations="all", complex_cations="all",
                        ligands="all", cofactors="all")
            self.ligand_ref = ligand_ref

    @property
    def ligand_ref(self):
        """ProtoCaller.Ensemble.Ligand.Ligand: The reference ligand."""
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
        """str: The absolute path to the PDB file for the protein."""
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
        self._checkfasta()

    @property
    def pdb_obj(self):
        """ProtoCaller.IO.PDB.PDB: The object corresponding to the PDB file."""
        return self._pdb_obj

    @property
    def fasta(self):
        """str: The absolute path to the FASTA file for the protein."""
        with self.workdir:
            if self._fasta is None:
                return self._downloader.getFASTA()
            return self._fasta

    @fasta.setter
    def fasta(self, value):
        if value is not None:
            value = _fileio.checkFileExists(value)
        self._fasta = value
        self._checkfasta()

    def filter(self, missing_residues="middle", chains="all", waters="site", ligands="chain", cofactors="chain",
               simple_anions="chain", complex_anions="chain", simple_cations="chain", complex_cations="chain",
               include_mols=None, exclude_mols=None):
        """
        Conditionally removes certain molecules from the PDB object and rewrites the PDB file.

        Parameters
        ----------
        missing_residues : str
            One of "all" or "middle". Determines whether to only add missing residues when they are non-terminal.
        chains : str
            One of "all" or an iterable of characters for chains to keep.
        waters : str or None
            One of "all", "chain" (only keep molecules belonging to a chain), "site" (only keep if they are mentioned in
            the PDB SITE directive) or None (no molecules are included).
        ligands : str or None
            One of "all", "chain" (only keep molecules belonging to a chain), "site" (only keep if they are mentioned in
            the PDB SITE directive) or None (no molecules are included).
        cofactors : str or None
            One of "all", "chain" (only keep molecules belonging to a chain), "site" (only keep if they are mentioned in
            the PDB SITE directive) or None (no molecules are included).
        simple_anions : str or None
            One of "all", "chain" (only keep molecules belonging to a chain), "site" (only keep if they are mentioned in
            the PDB SITE directive) or None (no molecules are included).
        complex_anions : str or None
            One of "all", "chain" (only keep molecules belonging to a chain), "site" (only keep if they are mentioned in
            the PDB SITE directive) or None (no molecules are included).
        simple_cations : str or None
            One of "all", "chain" (only keep molecules belonging to a chain), "site" (only keep if they are mentioned in
            the PDB SITE directive) or None (no molecules are included).
        complex_cations : str or None
            One of "all", "chain" (only keep molecules belonging to a chain), "site" (only keep if they are mentioned in
            the PDB SITE directive) or None (no molecules are included).
        include_mols : [str]
            A list of strings which specify molecules that should be included. Overrides
            previous filters.
        exclude_mols : [str]
            A list of strings which specify molecules that should be excluded. Overrides
            previous filters.
        """
        if include_mols is None: include_mols = []
        if exclude_mols is None: exclude_mols = []

        with self.workdir:
            # filter ligands / cofactors
            temp_dict = {"ligand": [], "cofactor": []}
            for molecule in self.ligands + self.cofactors:
                filename = molecule.name
                for param, name in zip([ligands, cofactors], ["ligand", "cofactor"]):
                    # turn the ligand into a pseudo-residue
                    _, resname, chainID, resSeq_iCode = _re.search(r"^([\w]+)_([\w]+)_([\w])_([A-Z0-9]+)", filename).groups()
                    resSeq, iCode = self._residTransform(resSeq_iCode)

                    # filter
                    if _PC.RESIDUETYPE(resname) == name and resSeq_iCode not in exclude_mols:
                        if any([param == "all",
                                param == "chain" and chains == "all",
                                param == "chain" and chainID in chains,
                                resSeq_iCode in include_mols]):
                            temp_dict[name] += [molecule]
                        elif param == "site":
                            for residue in self._pdb_obj.site_residues:
                                if residue.resSeq == resSeq and residue.iCode == iCode:
                                    temp_dict[name] += [molecule]
            self.ligands = temp_dict["ligand"]
            self.cofactors = temp_dict["cofactor"]

            # filter residues / molecules in protein

            # filter by chain
            filter = []
            mask = "type in ['amino_acid', 'amino_acid_modified']"
            if chains != "all":
                mask += "&chainID in %s" % str(chains)
            filter += self._pdb_obj.filter(mask)

            # filter missing residues
            if missing_residues == "middle":
                fasta = next(_SeqIO.parse(open(self.fasta), 'fasta'))
                seq = fasta.seq.tomutable()
                missing_residue_list = self._pdb_obj.totalResidueList()
                for i in range(2):
                    missing_residue_list.reverse()
                    seq.reverse()
                    current_chain = None
                    for j in reversed(range(0, len(missing_residue_list))):
                        res = missing_residue_list[j]
                        if type(res) is _PDB.Missing.MissingResidue and current_chain != res.chainID:
                            del missing_residue_list[j]
                            del seq[j]
                        else:
                            current_chain = res.chainID
                fasta.seq = seq
                _SeqIO.write(fasta, self.fasta, "fasta")
                missing_residue_list = [x for x in missing_residue_list if type(x) == _PDB.MissingResidue]
                filter += missing_residue_list
            else:
                filter += self._pdb_obj.missing_residues

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
                resSeq, iCode = self._residTransform(include_mol)
                filter_str = "resSeq=={}&iCode=='{}'&type not in ['ligand', " \
                             "'cofactor']".format(resSeq, iCode)
                residue = self._pdb_obj.filter(filter_str)
                if not len(residue):
                    _warnings.warn("Could not find residue {}.".format(include_mol))
                filter += residue

            # exclude extra molecules / residues
            excl_filter = []
            for exclude_mol in exclude_mols:
                resSeq, iCode = self._residTransform(exclude_mol)
                filter_str = "resSeq=={}&iCode=='{}'&type not in ['ligand', " \
                             "'cofactor']".format(resSeq, iCode)
                residue = self._pdb_obj.filter(filter_str)
                if not len(residue):
                    _warnings.warn(
                        "Could not find residue {}.".format(exclude_mol))
                excl_filter += residue

            filter = list(set(filter) - set(excl_filter))
            self._pdb_obj.purgeResidues(filter, "keep")
            self._pdb_obj.writePDB(self.pdb)

    def prepare(self, add_missing_residues="pdbfixer",
                add_missing_atoms="pdb2pqr", protonate_proteins="pdb2pqr",
                protonate_ligands="babel", missing_residues_options=None,
                missing_atom_options=None, protonate_proteins_options=None,
                protonate_ligands_options=None):
        """
        Adds missing residues / atoms to the protein and protonates it and the
        relevant ligands.

        Parameters
        ----------
        add_missing_residues : str or None
            How to add missing residues. One of "modeller", "charmm-gui",
            "pdbfixer" and None (no addition).
        add_missing_atoms : str or None
            How to add missing atoms. One of "modeller", "pdb2pqr", "pdbfixer"
            and None (no addition).
        protonate_proteins : str or None
            How to protonate the protein. One of "pdb2pqr" and
            None (no protonation).
        protonate_ligands : str or None
            How to protonate the related ligands / cofactors. One of "babel"
            and None (no protonation).
        missing_residues_options : dict
            Keyword arguments to pass on to the relevant wrapper responsible
            for the addition of missing protein residues.
        missing_atom_options : dict
            Keyword arguments to pass on to the relevant wrapper responsible
            for the addition of missing heavy atoms.
        protonate_proteins_options : dict
            Keyword arguments to pass on to the relevant wrapper responsible
            for protein protonation.
        protonate_ligands_options : dict
            Keyword arguments to pass on to the relevant wrapper responsible
            for ligand protonation.
        """
        with self.workdir:
            add_missing_residues = add_missing_residues.strip().lower() \
                if add_missing_residues is not None else ""
            add_missing_atoms = add_missing_atoms.strip().lower() \
                if add_missing_atoms is not None else ""
            protonate_proteins = protonate_proteins.strip().lower() \
                if protonate_proteins is not None else ""
            protonate_ligands = protonate_ligands.strip().lower() \
                if protonate_ligands is not None else ""

            if missing_residues_options is None:
                missing_residues_options = {}
            if missing_atom_options is None:
                missing_atom_options = {}
            if protonate_proteins_options is None:
                protonate_proteins_options = {}
            if protonate_ligands_options is None:
                protonate_ligands_options = {}

            # compatibility checks
            if add_missing_atoms == "pdb2pqr" and \
                    protonate_proteins != "pdb2pqr":
                _warnings.warn("Cannot currently run PDB2PQR without "
                               "protonation. This will be fixed in a later "
                               "version. Changing protein protonation method "
                               "to PDB2PQR...")
            if len(self._pdb_obj.missing_residues) and \
                    add_missing_residues == "modeller" and not _PC.MODELLER:
                _warnings.warn("Invalid Modeller license. Switching to CHARMM"
                               "-GUI for the addition of missing residues...")
                add_missing_residues = "charmm-gui"
            if len(self._pdb_obj.missing_atoms) and \
                    add_missing_atoms == "modeller" and not _PC.MODELLER:
                _warnings.warn("Invalid Modeller license. Switching to "
                               "PDB2PQR...")
                add_missing_atoms = "pdb2pqr"

            # add missing residues
            if len(self._pdb_obj.missing_residues):
                kwargs = missing_residues_options
                if add_missing_residues == "modeller":
                    if add_missing_atoms == "modeller":
                        atoms = True
                        kwargs = {**kwargs, **missing_atom_options}
                    else:
                        atoms = False
                    fasta = self._downloader.getFASTA()
                    self.pdb = _modeller.modellerTransform(self.pdb, fasta,
                                                           atoms, **kwargs)
                elif add_missing_residues == "charmm-gui":
                    self.pdb = _charmmwrap.charmmguiTransform(self.pdb,
                                                              **kwargs)
                elif add_missing_residues == "pdbfixer":
                    atoms = True if add_missing_atoms == "pdbfixer" else False
                    self.pdb = _pdbfix.pdbfixerTransform(self.pdb, True,
                                                         atoms, **kwargs)
                else:
                    _warnings.warn("Protein has missing residues. Please check"
                                   " your PDB file or choose a valid "
                                   "automation protocol")

            # add missing atoms
            if len(self.pdb_obj.missing_atoms):
                kwargs = missing_atom_options
                if add_missing_atoms == "modeller" and \
                        (not len(self._pdb_obj.missing_residues) or
                         add_missing_residues != "modeller"):
                    _warnings.warn("Cannot currently add missing atoms with "
                                   "Modeller when there are no missing "
                                   "residues. Switching to PDB2PQR...")
                    add_missing_atoms = "pdb2pqr"
                elif add_missing_atoms == "pdbfixer" and \
                         add_missing_residues != "pdbfixer":
                    self.pdb = _pdbfix.pdbfixerTransform(self.pdb, False, True,
                                                         **kwargs)

            # protonate proteins
            if "pdb2pqr" in [add_missing_atoms, protonate_proteins]:
                kwargs = {}
                if add_missing_atoms == "pdb2pqr":
                    kwargs = missing_atom_options
                if protonate_proteins == "pdb2pqr":
                    kwargs = {**kwargs, **protonate_proteins_options}
                self.pdb = _PDB2PQR.pdb2pqrTransform(self.pdb, **kwargs)

            # protonate ligands
            if len(self.ligands + self.cofactors):
                kwargs = protonate_ligands_options
                if protonate_ligands == "babel":
                    for ligand in self.ligands + self.cofactors:
                        if not ligand.protonated:
                            ligand.protonate(**kwargs)
                    if not self.ligand_ref.protonated:
                        self.ligand_ref.protonate(**kwargs)
                else:
                    _warnings.warn("Need to protonate all relevant ligands / "
                                   "cofactors before any parametrisation")

    def parametrise(self, params, reparametrise=False):
        """
        Parametrises the whole protein system.

        Parameters
        ----------
        params : ProtoCaller.Parametrise.Params
            Force field parameters.
        reparametrise : bool
            Whether to reparametrise an already parametrised complex.
        """
        if self.complex_template is not None and not reparametrise:
            print("Protein complex template %s is already parametrised." % self.name)
            return

        with self.workdir:
            print("Parametrising original crystal system...")
            # extract non-protein residues from pdb file and save them as separate pdb files
            hetatm_files, hetatm_types = self._pdb_obj.writeHetatms()
            filter = "type in ['amino_acid', 'amino_acid_modified']"
            non_protein_residues = self._pdb_obj.filter(filter)
            self._pdb_obj.purgeResidues(non_protein_residues, "keep")
            self._pdb_obj.reNumberResidues()
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
                id = _re.search(r"^([\w]+)_([\w]+)_([\w])_([A-Z0-9]+)", cofactor.name).group(2)
                cofactor.parametrise(params, reparametrise=reparametrise, molecule_type="cofactor", id=id)
                system += _pmdwrap.openFilesAsParmed(cofactor.parametrised_files)

            # add all other HETATMS to the system
            for filename, type in zip(hetatm_files, hetatm_types):
                system += _parametrise.parametriseAndLoadPmd(params=params, filename=filename, molecule_type=type)

            _GROMACS.saveAsGromacs("complex_template", system)
            if _PC.BIOSIMSPACE:
                self.complex_template = _BSS.IO.readMolecules(["complex_template.top", "complex_template.gro"])
            else:
                self.complex_template = system

    def _checkfasta(self):
        if hasattr(self, "_fasta") and hasattr(self, "_pdb_obj"):
            seqlen = len(next(_SeqIO.parse(open(self.fasta), 'fasta')).seq)
            reslen = len(self._pdb_obj.totalResidueList())
            if seqlen != reslen:
                _warnings.warn("Length of FASTA sequence ({}) does not match "
                               "the length of PDB sequence ({}). Please check "
                               "your input files.".format(seqlen, reslen))

    @staticmethod
    def _residTransform(id):
        """
        Transforms a string from e.g. "400G" to (400, "G")

        Parameters
        ----------
        id : str

        Returns
        -------
        resSeq : int
        iCode : str
        """
        id = id.strip()
        if id[-1].isalpha():
            iCode = id[-1]
            resSeq = int(float(id[:-1]))
        else:
            iCode = " "
            resSeq = int(float(id))
        return resSeq, iCode
