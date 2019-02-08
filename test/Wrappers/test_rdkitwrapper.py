from ProtoCaller.Wrappers.rdkitwrapper import *


def test_MCS():
    ref = openAsRdkit("C1CNC(CC1)C2CC(CCC2)C3CCCC3")
    mol = openAsRdkit("C1CCC(CC1)C2OC(CC2)C3CCNC3")
    assert len(getMCSMap(ref, mol, match="any", always_maximum=False)) == 6
    assert len(getMCSMap(ref, mol, match="any", always_maximum=True)) == 11

    ref = openAsRdkit("C1CCC(CC1)CCCC2CCCC2")
    mol = openAsRdkit("C1CCC(CC1)C2OC(CC2)C3CCNC3")
    assert len(getMCSMap(ref, mol)) == 14
