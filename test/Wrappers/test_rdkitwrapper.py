from ProtoCaller.Wrappers.rdkitwrapper import *


def test_MCS():
    # test some extended mapping
    ref = openAsRdkit("C1CNC(CC1)C2CC(CCC2)C3CCCC3")
    mol = openAsRdkit("C1CCC(CC1)C2OC(CC2)C3CCNC3")
    assert len(getMCSMap(ref, mol)[0]) == 11

    ref = openAsRdkit("C1CCC(CC1)CCCC2CCCC2")
    mol = openAsRdkit("C1CCC(CC1)C2OC(CC2)C3CCNC3")
    assert len(getMCSMap(ref, mol)[0]) == 14

    ref = openAsRdkit("C1CCCC1")
    mol = openAsRdkit("C1CCC1")
    assert len(getMCSMap(ref, mol)) == 0

    ref = openAsRdkit("c1Nc2c(c1)Cc1c2cccc1")
    mol = openAsRdkit("c1ccc(c(c1)C)c1Nccc1")
    assert len(getMCSMap(ref, mol)[0]) == 12

    ref = openAsRdkit("c1ccc2c(c1)Cc1c2cccc1")

    mol = openAsRdkit("c1ccc(c(c1)C)c1ccccc1")
    assert len(getMCSMap(ref, mol)[0]) == 13
    mol = openAsRdkit("c1ccc(c(c1)C)c1c(cccc1)C")
    assert len(getMCSMap(ref, mol)[0]) == 13

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
    mol = openAsRdkit("C3CC[C@@H](C1CCCC1)[C@@](C2CCC2)(C3)S")
    assert len(getMCSMap(ref, mol)[0]) == 15
