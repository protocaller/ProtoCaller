import logging as _logging
import os as _os
import requests as _requests

import pypdb as _pypdb

from ProtoCaller.IO.PDB import PDB as _PDB


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
        fasta_url = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList={}&" \
                    "compressionType=uncompressed".format(self.code.lower())
        fasta_filename = self._code + ".fasta"

        try:
            r = _requests.get(fasta_url)
            with open(fasta_filename, "wb") as f:
                f.write(r.content)
            self._fasta = _os.path.abspath(fasta_filename)
        except:
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
        ligand_metadata = _pypdb.get_ligands(self.code)
        try:
            ligand_metadata = ligand_metadata["ligandInfo"]["ligand"]
        except KeyError:
            return []

        if not isinstance(ligand_metadata, list):
            ligand_names = [ligand_metadata["@chemicalID"]]
        else:
            ligand_names = [x["@chemicalID"] for x in ligand_metadata]

        pdb_obj = _PDB(self.getPDB())
        matches = pdb_obj.filter("|".join(["resName=='{}'".format(x) for x in ligand_names]))
        full_ligand_names = ["{}_{}_{}_{}_NO_H.sdf".format(
            self.code.lower(), x.resName, x.chainID, x.resSeq) for x in matches]
        if ligands != "all":
            full_ligand_names = list(set(full_ligand_names) & set(ligands))

        filename_list = []
        _logging.info("Downloading ligand files from the Protein Data Bank...")
        for ligname in full_ligand_names:
            file_url = "https://files.rcsb.org/cci/view/" + ligname
            try:
                r = _requests.get(file_url)
                with open(ligname + ".sdf", "wb") as f:
                    f.write(r.content)
                filename_list += [_os.path.abspath(ligname + ".sdf")]
            except:
                _logging.warning("Could not download file: %s.sdf from %s" % (ligname, file_url))

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
        try:
            pdb_string = _pypdb.get_pdb_file(self._code, filetype="pdb", compression=False)
        except:
            raise ValueError("Invalid PDB code: '%s'" % self._code)

        filename = _os.path.abspath(self._code + ".pdb")
        with open(filename, "w") as file:
            file.write(pdb_string)

        self._pdb = filename
        return self._pdb
