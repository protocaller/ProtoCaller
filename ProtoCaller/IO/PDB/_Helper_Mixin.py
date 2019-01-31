class HelperMixin:
    def sameChain(self, obj):
        from . import Chain as _Chain
        try:
            return all([getattr(self, prop) == getattr(obj, prop) for prop in _Chain._common_properties])
        except AttributeError:
            print("Need to pass a valid object with {} attributes".format(_Chain._common_properties))

    def sameResidue(self, obj):
        from . import Residue as _Residue
        try:
            return all([getattr(self, prop) == getattr(obj, prop) for prop in _Residue._common_properties])
        except AttributeError:
            print("Need to pass a valid object with {} attributes".format(_Residue._common_properties))
