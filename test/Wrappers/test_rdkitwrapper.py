from ProtoCaller.Wrappers.rdkitwrapper import *


def test_MCS():
    ref = openAsRdkit("C1CNC(CC1)C2CC(CCC2)C3CCCC3")
    mol = openAsRdkit("C1CCC(CC1)C2OC(CC2)C3CCNC3")
    assert len(getMCSMap(ref, mol, match="any", always_maximum=False)) == 6
    assert len(getMCSMap(ref, mol, match="any", always_maximum=True)) == 11

    ref = openAsRdkit("C1CCC(CC1)CCCC2CCCC2")
    mol = openAsRdkit("C1CCC(CC1)C2OC(CC2)C3CCNC3")
    assert len(getMCSMap(ref, mol, always_maximum=True)) == 14

    # test some stereospecific mapping
    ref = openAsRdkit("C1CC[C@@H](C2CCCC2)[C@@H](C3CCC3)C1")

    mol = openAsRdkit("C1CC[C@@H](C2CCCC2)[C@H](C3CCC3)C1")
    assert len(getMCSMap(ref, mol, match="any", always_maximum=True, matchChiralTag=True)) == 11
    mol = openAsRdkit("C1CC[C@H](C2CCCC2)[C@@H](C3CCC3)C1")
    assert len(getMCSMap(ref, mol, match="any", always_maximum=True, matchChiralTag=True)) == 10
    mol = openAsRdkit("C1CC[C@H](C2CCCC2)[C@H](C3CCC3)C1")
    assert len(getMCSMap(ref, mol, match="any", always_maximum=True, matchChiralTag=True)) == 6
    mol = openAsRdkit("C1CC[C@@H](C2CCCC2)[C@@H](C3CCC3)C1")
    assert len(getMCSMap(ref, mol, match="any", always_maximum=True, matchChiralTag=True)) == 15
