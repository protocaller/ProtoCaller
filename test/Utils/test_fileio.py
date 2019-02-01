import ProtoCaller as PC
from ProtoCaller.Utils.fileio import Dir

import os


def test_Dir():
    orig_dir = os.getcwd()
    os.chdir(PC.TESTDIR + "/Utils")

    target_dirname = os.path.join(PC.TESTDIR, "Utils", "temp")
    f = Dir("temp", temp=True, purge_immediately=False)
    with f:
        assert os.getcwd() == target_dirname
        with f:
            assert os.getcwd() == target_dirname
        assert os.getcwd() == target_dirname
    assert os.path.exists(target_dirname)
    assert os.getcwd() == os.path.join(PC.TESTDIR, "Utils")

    target_dirname2 = os.path.join(PC.TESTDIR, "Utils", "temp2")
    with Dir("temp2", temp=True, purge_immediately=True):
        pass
    assert not os.path.exists(target_dirname2)

    os.chdir(orig_dir)
