import logging as _logging
import os as _os
import requests as _requests

import ProtoCaller as _PC
import ProtoCaller.IO.PDB as _PDB


class PDBDownloader:
    """
    A class which acts as a client to Protein Data Bank.

    Parameters
    ----------
    code : str
        Initialises code.

    Attributes
    ----------
    code : str
        The PDB code of the protein.
    """
    def __init__(self, code):
        self.code = code

    @property
    def code(self):
        return self._code

    @code.setter
    def code(self, val):
        self._fasta = None
        self._ligands = []
        self._pdb = None
        try:
            assert isinstance(val, str) and len(val) == 4
            self._code = val.upper()
        except AssertionError:
            raise ValueError("Need to pass a valid PDB code")

    def getFASTA(self):
        """
        Downloads the FASTA file from the Protein Data Bank.

        Returns
        -------
        pdb : str
            Returns the absolute path to the FASTA file downloaded from the Protein Data Bank.
        """
        if self._fasta:
            return self._fasta
        fasta_url = "https://www.rcsb.org/fasta/entry/{}".format(self.code.upper())
        fasta_filename = self._code + ".fasta"

        try:
            r = _requests.get(fasta_url)
            with open(fasta_filename, "wb") as f:
                f.write(r.content)
            self._fasta = _os.path.abspath(fasta_filename)
        except _requests.HTTPError:
            _logging.warning("Could not download file: %s from %s" % (fasta_filename, fasta_url))
            return -1

        return self._fasta

    def getLigands(self, ligands="all"):
        """
        Retrieves the ligands as SDF files from the Protein Data Bank.

        Parameters
        ----------
        ligands : str, [str]
            Either a list of Ligand ID's to keep or "all", in which case all ligands are downloaded.

        Returns
        -------
        ligands : [str]
            A list of the absolute paths of all ligands.
        """
        if self._ligands:
            return self._ligands

        pdb_obj = _PDB.PDB(self.getPDB())
        matches = []
        full_ligand_names = []
        urls = []
        for chain in pdb_obj:
            for residue in chain:
                if isinstance(residue, _PDB.Residue):
                    type = _PC.RESIDUETYPE(residue.resName)
                    if any([x in type for x in ["cofactor", "ligand"]]):
                        matches += [residue]
                        full_ligand_names += ["{}_{}_{}_{}{}_NO_H".format(
                            self.code.lower(), residue.resName, residue.chainID, residue.resSeq, residue.iCode.strip())]
                        url = f"https://models.rcsb.org/v1/{self.code.lower()}/ligand?auth_asym_id={residue.chainID}&" \
                              f"label_comp_id={residue.resName}&auth_seq_id={residue.resSeq}&encoding=sdf&" \
                              f"copy_all_categories=false"
                        if len(residue.iCode.strip()):
                            url += f"&pdbx_PDB_ins_code={residue.iCode}"
                        urls += [url]

        filename_list = []
        _logging.info("Downloading ligand files from the Protein Data Bank...")
        for url, ligname in zip(urls, full_ligand_names):
            if ligands == "all" or ligname in ligands:
                try:
                    r = _requests.get(url)
                    content = r.content.decode()
                    if "error" in content:
                        raise _requests.HTTPError
                    content = content.split("\n")
                    content[0] = ligname
                    with open(ligname + ".sdf", "w") as f:
                        f.write("\n".join(content))
                    filename_list += [_os.path.abspath(ligname + ".sdf")]
                except _requests.HTTPError:
                    _logging.warning("Could not download file: %s.sdf from %s" % (ligname, url))

        self._ligands = filename_list
        return self._ligands

    def getPDB(self):
        """
        Downloads the PDB file from the Protein Data Bank.

        Returns
        -------
        pdb : str
            Returns the absolute path to the PDB file downloaded from the Protein Data Bank.
        """
        if self._pdb:
            return self._pdb

        pdb_url = "https://files.rcsb.org/download/{}.pdb".format(self.code.upper())
        pdb_filename = self._code + ".pdb"
        try:
            r = _requests.get(pdb_url)
            with open(pdb_filename, "wb") as f:
                f.write(r.content)
            self._pdb = _os.path.abspath(pdb_filename)
        except _requests.HTTPError:
            _logging.warning("Could not download file: %s from %s" % (pdb_filename, pdb_url))
            return -1

        return self._pdb
