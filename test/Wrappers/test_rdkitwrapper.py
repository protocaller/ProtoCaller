from pytest import approx

import ProtoCaller as PC
from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Wrappers.rdkitwrapper import *


def test_MCS_extended():
    # test some extended mapping
    ref = openAsRdkit("C1CNC(CC1)C2CC(CCC2)C3CCCC3", minimise=False)
    mol = openAsRdkit("C1CCC(CC1)C2OC(CC2)C3CCNC3", minimise=False)
    assert len(getMCSMap(ref, mol)[0]) == 11

    ref = openAsRdkit("C1CCC(CC1)CCCC2CCCC2", minimise=False)
    mol = openAsRdkit("C1CCC(CC1)C2OC(CC2)C3CCNC3", minimise=False)
    assert len(getMCSMap(ref, mol)[0]) == 14

    ref = openAsRdkit("C1CCCC1", minimise=False)
    mol = openAsRdkit("C1CCC1", minimise=False)
    assert len(getMCSMap(ref, mol)) == 0

    ref = openAsRdkit("c1Nc2c(c1)Cc1c2cccc1", minimise=False)
    mol = openAsRdkit("c1ccc(c(c1)C)c1Nccc1", minimise=False)
    assert len(getMCSMap(ref, mol)[0]) == 12
    assert len(getMCSMap(mol, ref, two_way_matching=False)[0]) == 6

    ref = openAsRdkit("c1ccc2c(c1)Cc1c2cccc1", minimise=False)

    mol = openAsRdkit("c1ccc(c(c1)C)c1ccccc1", minimise=False)
    assert len(getMCSMap(ref, mol)[0]) == 13
    mol = openAsRdkit("c1ccc(c(c1)C)c1c(cccc1)C", minimise=False)
    assert len(getMCSMap(ref, mol)[0]) == 13

    ref = openAsRdkit("c1cc2c(cc1)c1c([nH]2)c2c([nH]1)c1c([nH]2)cc[nH]1", minimise=False)
    mol = openAsRdkit("c1cc2c(cc1)c1c(c3c([nH]1)c1c([nH]3)cc[nH]1)cc2", minimise=False)
    assert len(getMCSMap(ref, mol)[0]) == 15

    ref = openAsRdkit("c1ccc2c(c1)c1c([nH]2)[nH]c2c1[nH]cc2", minimise=False)
    mol = openAsRdkit("c1cccc(c1)c1c[nH]c2c1[nH]cc2", minimise=False)
    assert len(getMCSMap(ref, mol)[0]) == 14


def test_MCS_EZ():
    with Dir(PC.TESTDIR + "/shared"):
        # test mapping of an asymmetric double bond onto symmetric one
        ref = openFileAsRdkit("EZ_ref1.mol2", removeHs=False)
        mol = openFileAsRdkit("EZ_mol1.mol2", removeHs=False)
        results = getMCSMap(ref, mol)
        assert all([len(x) == 38 for x in results])
        assert all([{(11, 12), (12, 11)}.issubset(x) for x in results])

        # test mapping between two esters of different E/Z conformations
        ref = openFileAsRdkit("EZ_ref2.mol2", removeHs=False)
        mol = openFileAsRdkit("EZ_mol2.mol2", removeHs=False)
        results = getMCSMap(ref, mol)
        assert all([len(x) == 28 for x in results])
        results_rev = getMCSMap(mol, ref)
        assert all([len(x) == 28 for x in results_rev])

        # test mapping of a double bond to a single bond with unfavourable
        # conformation
        ref = openFileAsRdkit("EZ_ref3.mol2", removeHs=False)
        mol = openFileAsRdkit("EZ_mol3.mol2", removeHs=False)
        results = getMCSMap(ref, mol, timeout=1, two_way_matching=True)
        assert all([len(x) == 38 for x in results])
        results = getMCSMap(ref, mol, timeout=1, two_way_matching=False)
        assert all([len(x) == 27 for x in results])

        # test mapping of a double bond to a single bond with favourable
        # conformation
        ref = openFileAsRdkit("EZ_ref4.mol2", removeHs=False)
        mol = openFileAsRdkit("EZ_mol4.mol2", removeHs=False)
        results = getMCSMap(ref, mol, timeout=1, two_way_matching=True)
        assert all([len(x) == 41 for x in results])
        results = getMCSMap(ref, mol, timeout=1, two_way_matching=False)
        assert all([len(x) == 27 for x in results])

        # test the recursive algorithm for > 1 mismatching bonds
        ref = openFileAsRdkit("EZ_ref5.mol2", removeHs=False)
        mol = openFileAsRdkit("EZ_mol5.mol2", removeHs=False)
        results = getMCSMap(ref, mol, timeout=1)
        assert all([len(x) == 17 for x in results])
        results_rev = getMCSMap(ref, mol, timeout=1)
        assert all([len(x) == 17 for x in results_rev])

        # test the recursive algorithm for > 1 mismatching bonds and matching
        # a double bond onto a single bond
        ref = openFileAsRdkit("EZ_ref6.mol2", removeHs=False)
        mol = openFileAsRdkit("EZ_mol6.mol2", removeHs=False)
        results = getMCSMap(ref, mol, timeout=1, two_way_matching=True)
        assert all([len(x) == 17 for x in results])

        # test mapping of mismatching amide onto an ester
        ref = openFileAsRdkit("EZ_ref7.mol2", removeHs=False)
        mol = openFileAsRdkit("EZ_mol7.mol2", removeHs=False)
        results = getMCSMap(ref, mol, timeout=1)
        assert all([len(x) == 30 for x in results])


def test_MCS_RS():
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


def test_align():
    from ProtoCaller.Wrappers.rdkitwrapper import _nonMCSDihedrals
    with Dir(PC.TESTDIR + "/shared"):
        # test conservation of dihedrals
        ref = openFileAsRdkit("Align_ref1.mol2", removeHs=False)
        mol = openFileAsRdkit("Align_mol1.mol2", removeHs=False)
        mcs = getMCSMap(ref, mol)[0]
        dihedrals = _nonMCSDihedrals(mol, mcs)
        assert len(dihedrals) == 3
        mol_new, mcs_new = alignTwoMolecules(ref, mol, minimise_score=False)
        assert mcs == mcs_new
        dihedrals_new = _nonMCSDihedrals(mol_new, mcs_new)
        assert dihedrals.keys() == dihedrals_new.keys()
        assert approx(dihedrals.values(), dihedrals.values())