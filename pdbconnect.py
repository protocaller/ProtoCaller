import const as _const

import re as _re
import urllib as _urllib

class PDBDownloader:
    def __init__(self, id):
        self.id = id
        self.download_html()

    def download_html(self):
        print("Attempting to access Protein Data Bank... ", end="")
        self.html = _urllib.request.urlopen("https://www.rcsb.org/structure/" + self.id.upper()).read().decode('utf-8')
        if (not self.html):
            print("FAILED")
        else:
            print("OK")

    def downloadLigandSDF(self, ligands="all"):
        raw_file_url_list = _re.findall(r'"/pdb/download/downloadLigandFiles.do\?ligandIdList=[^"]*"', self.html)
        if (ligands != "all"):
            for ligand in ligands:
                raw_file_url_list = [x for x in raw_file_url_list if "ligandIdList=%s&" % ligand in x]

        filename_list = []
        for raw_file_url in raw_file_url_list:
            ligname = _re.search(r'ligandIdList=([\w\d]*)', raw_file_url).group(1)
            file_url = "https://www.rcsb.org" + raw_file_url[1:-1].replace("&amp;", "&")
            try:
                _urllib.request.urlretrieve(file_url, ligname + ".sdf")
                filename_list += [ligname + ".sdf"]
            except:
                print("Could not download file: %s.sdf from %s" % (ligname, file_url))
                return -1

        print("Finished downloading.")

        if filename_list == []:
            filename_list = None
        return filename_list

    def downloadFASTA(self):
        fasta_url = _re.search(r'"/pdb/[\S]*downloadFasta[\S]*"', self.html).group(0)
        if(fasta_url):
            fasta_url = "https://www.rcsb.org" + fasta_url[1:-1].replace("&amp;", "&")

        fasta_filename = self.id + ".fasta"
        try:
            _urllib.request.urlretrieve(fasta_url, fasta_filename)
        except:
            print("Could not download file: %s from %s" % (fasta_filename, fasta_url))
            return -1
        return fasta_filename