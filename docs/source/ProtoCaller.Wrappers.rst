ProtoCaller.Wrappers package
============================

This package is responsible for most of the external linkage between
ProtoCaller and other pieces of software. While these packages have been
extensively incorporated into the more user-friendly :mod:`Ensemble` package,
it is still possible that the user might find these useful on their own.

There are two types of wrappers here. The first type is the black box wrapper,
which calls an external program to modify an object and then returns a modified
version thereof. Such submodules are :mod:`~ProtoCaller.Wrappers.\
PDB2PQRwrapper`, :mod:`~ProtoCaller.Wrappers.charmmguiwrapper`,
:mod:`~ProtoCaller.Wrappers.modellerwrapper`,
:mod:`~ProtoCaller.Wrappers.pdbfixerwrapper`,
:mod:`~ProtoCaller.Wrappers.babelwrapper` and the not-yet-functional
:mod:`~ProtoCaller.Wrappers.protosswrapper`. The second type is a wrapper
that provides additional functionality. E.g. in :mod:`~ProtoCaller.Wrappers.\
BioSimSpacewrapper`, :mod:`~ProtoCaller.Wrappers.parmedwrapper` and
:mod:`~ProtoCaller.Wrappers.rdkitwrapper` one can find functions which extend
the functionalities of the corresponding packages. A more notable such function
is the :func:`~ProtoCaller.Wrappers.rdkitwrapper.getMCSMap` function, which
expands on the Maximum Common Substructure (MCS) capabilities in RDKit.

Submodules
----------

.. toctree::

   ProtoCaller.Wrappers.BioSimSpacewrapper
   ProtoCaller.Wrappers.PDB2PQRwrapper
   ProtoCaller.Wrappers.babelwrapper
   ProtoCaller.Wrappers.charmmguiwrapper
   ProtoCaller.Wrappers.modellerwrapper
   ProtoCaller.Wrappers.parmedwrapper
   ProtoCaller.Wrappers.pdbfixerwrapper
   ProtoCaller.Wrappers.protosswrapper
   ProtoCaller.Wrappers.rdkitwrapper

Module contents
---------------

.. automodule:: ProtoCaller.Wrappers
    :members:
    :undoc-members:
    :show-inheritance:
