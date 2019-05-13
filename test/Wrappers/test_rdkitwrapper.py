from ProtoCaller.Wrappers.rdkitwrapper import *


def test_MCS():
    ref = openAsRdkit("C1CNC(CC1)C2CC(CCC2)C3CCCC3")
    mol = openAsRdkit("C1CCC(CC1)C2OC(CC2)C3CCNC3")
    assert len(getMCSMap(ref, mol)[0]) == 11

    ref = openAsRdkit("C1CCC(CC1)CCCC2CCCC2")
    mol = openAsRdkit("C1CCC(CC1)C2OC(CC2)C3CCNC3")
    assert len(getMCSMap(ref, mol)[0]) == 14

    ref = openAsRdkit("C1CCCC1")
    mol = openAsRdkit("C1CCC1")
    assert len(getMCSMap(ref, mol)[0]) == 0

    # test some stereospecific mapping
    ref = openAsRdkit("C1CC[C@@H](C2CCCC2)[C@@H](C3CCC3)C1")

    mol = openAsRdkit("C1CC[C@@H](C2CCCC2)[C@H](C3CCC3)C1")
    assert len(getMCSMap(ref, mol)[0]) == 11
    mol = openAsRdkit("C1CC[C@H](C2CCCC2)[C@@H](C3CCC3)C1")
    assert len(getMCSMap(ref, mol)[0]) == 10
    mol = openAsRdkit("C1CC[C@H](C2CCCC2)[C@H](C3CCC3)C1")
    assert len(getMCSMap(ref, mol)[0]) == 6
    mol = openAsRdkit("C1CC[C@@H](C2CCCC2)[C@@H](C3CCC3)C1")
    assert len(getMCSMap(ref, mol)[0]) == 15
