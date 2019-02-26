from .Ligand import Ligand as _Ligand
from .Perturbation import Perturbation as _Perturbation
from ProtoCaller.Utils.ConditionalList import ConditionalList as _ConditionalList


class PerturbationList(_ConditionalList):
    """
    A ProtoCaller.Utils.ConditionalList.ConditionalList of ProtoCaller.Ensemble.Perturbation.Perturbation
    """
    def __init__(self, perturbations):
        _ConditionalList.__init__(self, perturbations, transformfunc=self._transformMorph)

    @staticmethod
    def _transformMorph(morph):
        if isinstance(morph, _Perturbation):
            return morph
        elif isinstance(morph, list):
            if len(morph) != 2:
                raise ValueError("Need a list of 2 Ligands for a perturbation")
            return _Perturbation(*[x if isinstance(x, _Ligand) else _Ligand(x) for x in morph])
