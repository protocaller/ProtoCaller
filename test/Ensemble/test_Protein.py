import ProtoCaller as PC
from ProtoCaller.Ensemble import Protein, Ligand
from ProtoCaller.Parametrise import Params
from ProtoCaller.Utils.fileio import Dir

if PC.BIOSIMSPACE:
    import BioSimSpace as BSS

import pytest


def test_prepare_1bji():
    with Dir(PC.TESTDIR + "/Ensemble"):
        with Dir("temp", temp=True):
            protein = Protein("1bji", ligand_ref="479G")
            assert isinstance(protein.ligand_ref, Ligand)
            assert len(protein.ligands) == 9
            assert protein._pdb_obj.numberOfAtoms == 3259

            protein.filter(waters=None, include_mols=["480"], ligands=None)
            assert isinstance(protein.ligand_ref, Ligand)
            assert len(protein.ligands) == 0
            assert protein._pdb_obj.numberOfAtoms == 3069

            # this time keep all waters and parametrise with some end state checks
            protein = Protein("1bji", ligand_ref="479G")
            protein.filter(ligands=None, waters="all")
            protein.prepare()
            protein.parametrise(Params())

            if PC.BIOSIMSPACE and isinstance(protein.complex_template, BSS._SireWrappers._system.System):
                assert protein.complex_template.nAtoms() == 6543
                assert pytest.approx(protein.complex_template.charge().magnitude()) == 2
            else:
                # TODO: add support for ParmEd topologies when the Morph class is fully functional
                pass
