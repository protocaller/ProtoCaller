ProtoCaller.Parametrise package
===============================

This package is responsible for parametrising molecules to a particular
force field. So far, ProtoCaller only supports automatic parametrisation
via AmberTools (antechamber and tleap). Proteins, ligands, water molecules
and simple hard anions and cations can be readily parametrised using Amber\
tools. However, cofactors and transition metal complexes are not trivial to
parametrise automatically and as such support for these systems is currently
limited.

The most important elements in this module are the :class:`~ProtoCaller\
.Parametrise.Params` class, which is a simple, yet needed by ProtoCaller
class, which describes the desired force fields and the :func:`~ProtoCaller.\
Parametrise.parametriseFile` function, which handles the actual parametrisaton.

Submodules
----------

.. toctree::

   ProtoCaller.Parametrise.amber

Module contents
---------------

.. automodule:: ProtoCaller.Parametrise
    :members:
    :undoc-members:
    :show-inheritance:
