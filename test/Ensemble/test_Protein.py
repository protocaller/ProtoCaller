import ProtoCaller as PC
from ProtoCaller.Ensemble.Protein import Protein
from ProtoCaller.Parametrise import Params
from ProtoCaller.Utils.fileio import Dir

if PC.BIOSIMSPACE:
    import BioSimSpace as BSS

from rdkit.Chem.rdchem import Mol
import pytest


def test_prepare_1bji():
    with Dir(PC.TESTDIR + "/Ensemble"):
        with Dir("temp", temp=True):
            protein = Protein("1bji", ligand_ref="479G")
            assert isinstance(protein.ligand_ref, Mol)
            assert len(protein.ligands) == 9
            assert protein._protein_obj.numberOfAtoms == 3259

            protein.filter(waters=None, include_mols=["resSeq==480"], exclude_mols=["type=='simple_cation'"],
                           ligands=None)
            assert isinstance(protein.ligand_ref, Mol)
            assert len(protein.ligands) == 0
            assert protein._protein_obj.numberOfAtoms == 3068

            # this time keep only site waters and parametrise with some end state checks
            protein = Protein("1bji", ligand_ref="479G")
            protein.filter(ligands=None, waters="site")
            protein.prepare()
            protein.parametrise(Params())

            if PC.BIOSIMSPACE and isinstance(protein.complex_template, BSS._SireWrappers._system.System):
                assert protein.complex_template.nAtoms() == 6003
                assert pytest.approx(protein.complex_template.charge().magnitude()) == -1
            else:
                # TODO: add support for ParmEd topologies when the Morph class is fully functional
                pass
