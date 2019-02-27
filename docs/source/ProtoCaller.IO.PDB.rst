ProtoCaller.IO.PDB package
==========================

This module includes a PDB parser which takes care of many details that are
ignored by other PDB parsers. Here we take care of missing residues, missing
atoms, disulfide bonds, modified residues and site residues in addition to
the regular atom information. By doing so, ProtoCaller makes sure that the use
of external tools which modify PDB files will not result in any data loss
when e.g. feeding in a protein into Modeller and then PDB2PQR to prepare it
while still keeping track of disulfide bonds.

The :class:`~ProtoCaller.IO.PDB.PDB` class is highly structured into
:class:`~ProtoCaller.IO.PDB.Chain.Chain`, :class:`~ProtoCaller.IO.PDB.Chain.\
Chain` into :class:`~ProtoCaller.IO.PDB.Residue.Residue` and
:class:`~ProtoCaller.IO.PDB.Residue.Residue` into :class:`~ProtoCaller.IO.PDB.\
Atom.Atom` classes. In addition, there are the :class:`~ProtoCaller.IO.PDB.\
Missing.MissingResidue` and :class:`~ProtoCaller.IO.PDB.Missing.MissingAtoms`
classes which act as placeholders for missing elements in the PDB file.

This structure enables sophisticated manipulations of the PDB file inside
ProtoCaller when the input PDB file is not readily usable and manually editing
it would be non-trivial. An example would be the preparation of :doc:`Examples\
.Sialyltransferase`, where the PDB file contains a modified selenomethionine
residues which the user might want to convert into methionine residues due to
lack of parameters. This example also shows that in the case where the same
atom might have two alternative positions, the user could choose one of them
from inside ProtoCaller.

Submodules
----------

.. toctree::

   ProtoCaller.IO.PDB.Atom
   ProtoCaller.IO.PDB.Chain
   ProtoCaller.IO.PDB.Missing
   ProtoCaller.IO.PDB.Residue

Module contents
---------------

.. automodule:: ProtoCaller.IO.PDB
    :members:
    :undoc-members:
    :show-inheritance:
