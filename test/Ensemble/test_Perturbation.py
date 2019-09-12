import ProtoCaller as PC
from ProtoCaller.Ensemble.Ligand import Ligand
from ProtoCaller.Ensemble.Perturbation import Perturbation
from ProtoCaller.Utils.fileio import Dir


def test_align():
    with Dir(PC.TESTDIR + "/shared"):
        with Dir("temp", temp=True):
            lig_ref = Ligand("../toluene.sdf", protonated=True, minimise=False)
            lig1 = Ligand("c1ccccc1CC", minimise=False)
            lig2 = Ligand("c1ccccc1CCC", minimise=False)

            morph = Perturbation(lig1, lig2)
            mol, mcs = morph.alignAndCreateMorph(lig_ref)
            assert len(mcs) == 18

            positions = {
                (-32.355, 7.263, 2.207),
                (-33.181, 6.606, 3.301),
                (-33.238, 7.160, 4.565),
                (-33.992, 6.571, 5.570),
                (-34.711, 5.417, 5.332),
                (-34.638, 4.878, 4.060),
                (-33.894, 5.447, 3.048),
                (-31.882, 8.141, 2.595),
                (-32.993, 7.533, 1.392),
                (-31.608, 6.578, 1.864),
                (-32.712, 8.025, 4.765),
                (-34.016, 7.001, 6.508),
                (-35.279, 4.975, 6.072),
                (-35.163, 4.013, 3.860),
                (-33.869, 5.013, 2.112),
            }

            vec1 = mol._sire_object.property("coordinates0").toVector()
            positions_A = {(x[0], x[1], x[2]) for x in vec1}
            vec2 = mol._sire_object.property("coordinates1").toVector()
            positions_B = {(x[0], x[1], x[2]) for x in vec2}

            assert positions.issubset(positions_A)
            assert positions.issubset(positions_B)
            assert len(positions_A) == len(positions_B) == 21
            assert positions_A == positions_B
