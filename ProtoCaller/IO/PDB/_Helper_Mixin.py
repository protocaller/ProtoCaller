class HelperMixin:
    """A helper class containing some common helper functions."""
    def sameChain(self, obj):
        """bool: Returns whether the object is in the same chain as obj."""
        from . import Chain as _Chain
        try:
            return all([getattr(self, prop) == getattr(obj, prop) for prop in _Chain._common_properties])
        except AttributeError:
            print("Need to pass a valid object with {} attributes".format(_Chain._common_properties))

    def sameResidue(self, obj):
        """bool: Returns whether the object is in the same residue as obj."""
        from . import Residue as _Residue
        try:
            return all([getattr(self, prop) == getattr(obj, prop) for prop in _Residue._common_properties])
        except AttributeError:
            print("Need to pass a valid object with {} attributes".format(_Residue._common_properties))
