import ProtoCaller as PC
from ProtoCaller.Ensemble.Ligand import Ligand
from ProtoCaller.Ensemble.Perturbation import Perturbation
from ProtoCaller.Utils.fileio import Dir


def test_align():
    with Dir(PC.TESTDIR + "/shared"):
        with Dir("temp", temp=True):
            lig_ref = Ligand("../toluene.sdf", protonated=False)
            lig1 = Ligand("c1ccccc1CC")
            lig2 = Ligand("c1ccccc1CCC")

            morph = Perturbation(lig1, lig2)
            mol, mcs = morph.alignAndCreateMorph(lig_ref)

            positions = [(-33.238, 7.16, 4.565), (-33.992, 6.571, 5.57), (-34.711, 5.417, 5.332),
                         (-34.638, 4.878, 4.06), (-33.894, 5.447, 3.048), (-33.181, 6.606, 3.301),
                         (-32.355, 7.263, 2.207), (-30.927, 6.741, 2.192), (-32.686, 8.07, 4.789),
                         (-34.019, 7.027, 6.557), (-35.305, 4.95, 6.11), (-35.195, 3.966, 3.849),
                         (-33.872, 4.983, 2.066), (-32.828, 7.084, 1.234), (-32.354, 8.35, 2.351),
                         (-30.422, 6.939, 3.144), (-30.901, 5.66, 2.017), (-30.354, 7.228, 1.397),
                         (-30.715, 6.347, 4.026), (-30.416, 7.977, 3.515), (-29.339, 6.718, 3.124)]

            vec1 = mol._sire_molecule.property("coordinates0").toVector()
            vec2 = mol._sire_molecule.property("coordinates1").toVector()

            # need element by element check because Vector and Tuple are different types but they both support []
            for item1, item2, item3 in zip(positions, vec1, vec2):
                for i in range(3):
                    assert item1[i] == item2[i] == item3[i]

            assert len(vec1) == len(vec2) == len(positions)
