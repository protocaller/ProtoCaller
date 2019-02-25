import ProtoCaller as PC
from ProtoCaller.Ensemble.Ligand import Ligand
from ProtoCaller.Ensemble.Perturbation import Perturbation
from ProtoCaller.Utils.fileio import Dir


def test_align():
    with Dir(PC.TESTDIR + "/shared"):
        with Dir("temp", temp=True):
            lig_ref = Ligand("../toluene.sdf", protonated=False, minimise=False)
            lig1 = Ligand("c1ccccc1CC", minimise=False)
            lig2 = Ligand("c1ccccc1CCC", minimise=False)

            morph = Perturbation(lig1, lig2)
            mol, mcs = morph.alignAndCreateMorph(lig_ref)

            positions = [(-33.238, 7.16, 4.565), (-33.992, 6.571, 5.57), (-34.711, 5.417, 5.332),
                         (-34.638, 4.878, 4.06), (-33.894, 5.447, 3.048), (-33.181, 6.606, 3.301),
                         (-32.355, 7.263, 2.207), (-33.192, 8.279, 1.437), (-32.691, 8.069, 4.784),
                         (-34.021, 7.024, 6.553), (-35.304, 4.952, 6.109), (-35.192, 3.972, 3.85),
                         (-33.876, 4.987, 2.067), (-31.974, 6.493, 1.502), (-31.472, 7.774, 2.648),
                         (-33.555, 9.075, 2.122), (-34.063, 7.779, 0.962), (-32.573, 8.746, 0.643),
                         (-32.758, 9.628, 2.683), (-34.361, 8.772, 2.834), (-34.071, 9.887, 1.543)]

            vec1 = mol._sire_molecule.property("coordinates0").toVector()
            vec2 = mol._sire_molecule.property("coordinates1").toVector()

            # need element by element check because Vector and Tuple are different types but they both support []
            for item1, item2, item3 in zip(positions, vec1, vec2):
                for i in range(3):
                    assert item1[i] == item2[i] == item3[i]

            assert len(vec1) == len(vec2) == len(positions)
