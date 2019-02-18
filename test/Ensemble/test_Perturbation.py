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
                         (-32.355, 7.263, 2.207), (-30.927, 6.728, 2.218), (-33.477, 8.244, 4.515), 
                         (-33.434, 6.55, 6.532), (-34.505, 4.631, 6.091), (-35.606, 4.442, 3.728), 
                         (-34.521, 5.582, 2.139), (-32.81, 7.067, 1.212), (-32.333, 8.365, 2.351), 
                         (-30.34, 7.218, 1.413), (-30.443, 6.944, 3.194), (-30.926, 5.631, 2.042), 
                         (-30.711, 7.014, 0.374), (-29.284, 6.843, 1.397), (-30.227, 8.328, 1.526)]

            vec1 = mol._sire_molecule.property("coordinates0").toVector()
            vec2 = mol._sire_molecule.property("coordinates1").toVector()

            # need element by element check because Vector and Tuple are different types but they both support []
            for item1, item2, item3 in zip(positions, vec1, vec2):
                for i in range(3):
                    assert item1[i] == item2[i] == item3[i]

            assert len(vec1) == len(vec2) == len(positions)
