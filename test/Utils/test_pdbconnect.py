import ProtoCaller as PC
from ProtoCaller.Utils.pdbconnect import PDBDownloader
from ProtoCaller.Utils.fileio import Dir

def test_pdbconnect():
    with Dir(PC.TESTDIR + "/Utils/temp", temp=True):
        downloader = PDBDownloader("1BJI")
        downloader.getPDB()
        downloader.getFASTA()
        assert len(downloader.getLigands()) == 11