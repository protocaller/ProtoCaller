import os as _os
import re as _re
import urllib as _urllib

from ProtoCaller.shared import pypdb as _pypdb


class PDBDownloader:
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
        self._html = _urllib.request.urlopen("https://www.rcsb.org/structure/" + self.code).read().decode('utf-8')
        print("OK") if self._html else "FAILED"

    def getFASTA(self):
        self._download_html()
        if self._fasta: return self._fasta
        fasta_url = _re.search(r'"/pdb/[\S]*downloadFasta[\S]*"', self._html).group(0)
        fasta_filename = self._code + ".fasta"

        try:
            fasta_url = "https://www.rcsb.org" + fasta_url[1:-1].replace("&amp;", "&")
            _urllib.request.urlretrieve(fasta_url, fasta_filename)
            self._fasta = _os.path.abspath(fasta_filename)
        except:
            print("Could not download file: %s from %s" % (fasta_filename, fasta_url))
            return -1

        return self._fasta

    def getLigands(self, ligands="all"):
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
                _urllib.request.urlretrieve(file_url, ligname + ".sdf")
                filename_list += [_os.path.abspath(ligname + ".sdf")]
            except:
                print("Could not download file: %s.sdf from %s" % (ligname, file_url))
        print("Finished downloading.")

        self._ligands = filename_list
        return self._ligands

    def getPDB(self):
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
