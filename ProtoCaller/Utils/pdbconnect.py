import os as _os
import re as _re
import requests as _requests

from ProtoCaller.shared import pypdb as _pypdb


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
        self._html = []
        self._fasta = None
        self._ligands = []
        self._pdb = None
        try:
            assert isinstance(val, str) and len(val) == 4
            self._code = val.upper()
        except AssertionError:
            raise ValueError("Need to pass a valid PDB code")

    def _download_html(self):
        if self._html: return
        print("Attempting to access Protein Data Bank... ", end="")
        self._html = _requests.get("https://www.rcsb.org/structure/" + self.code).text
        print("OK") if self._html else "FAILED"

    def getFASTA(self):
        """
        Downloads the FASTA file from the Protein Data Bank.

        Returns
        -------
        pdb : str
            Returns the absolute path to the FASTA file downloaded from the Protein Data Bank.
        """
        self._download_html()
        if self._fasta: return self._fasta
        fasta_url = _re.search(r'"/pdb/[\S]*downloadFasta[\S]*"', self._html).group(0)
        fasta_filename = self._code + ".fasta"

        try:
            fasta_url = "https://www.rcsb.org" + fasta_url[1:-1].replace("&amp;", "&")
            r = _requests.get(fasta_url)
            with open(fasta_filename, "wb") as f:
                f.write(r.content)
            self._fasta = _os.path.abspath(fasta_filename)
        except:
            print("Could not download file: %s from %s" % (fasta_filename, fasta_url))
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
        self._download_html()
        if self._ligands: return self._ligands
        raw_file_url_list = _re.findall(r'"/pdb/download/downloadLigandFiles.do\?ligandIdList=[^"]*"', self._html)
        if ligands != "all":
            for ligand in ligands:
                raw_file_url_list = [x for x in raw_file_url_list if "ligandIdList=%s&" % ligand in x]

        filename_list = []
        for raw_file_url in raw_file_url_list:
            ligname = _re.search(r'ligandIdList=([\w\d]*)', raw_file_url).group(1)
            file_url = "https://www.rcsb.org" + raw_file_url[1:-1].replace("&amp;", "&")
            try:
                r = _requests.get(file_url)
                with open(ligname + ".sdf", "wb") as f:
                    f.write(r.content)
                filename_list += [_os.path.abspath(ligname + ".sdf")]
            except:
                print("Could not download file: %s.sdf from %s" % (ligname, file_url))
        print("Finished downloading.")

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
        if self._pdb: return self._pdb
        try:
            pdb_string = _pypdb.get_pdb_file(self._code, filetype="pdb", compression=False)
        except:
            raise ValueError("Invalid PDB code: '%s'" % self._code)

        filename = _os.path.abspath(self._code + ".pdb")
        with open(filename, "w") as file:
            file.write(pdb_string)

        self._pdb = filename
        return self._pdb
