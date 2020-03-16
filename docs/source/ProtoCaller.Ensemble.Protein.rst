ProtoCaller.Ensemble.Protein module
===================================

This module handles the :class:`~ProtoCaller.Ensemble.Protein.Protein` class which
keeps track of protein instantiation, missing residue/atom addition, PDB file
modification, protonation and parametrisation. More broadly, it contains the
PDB file of the protein, the corresponding SDF files of any relevant extra
ligands/cofactors and a file containing the FASTA sequence, all downloaded from
the `Protein Data Bank <https://www.rcsb.org/>`_. If the user supplies these
manually, they need to make sure that they provide ProtoCaller with all relevant
information.

This class handles calls to various wrappers in the :mod:`~ProtoCaller.Wrappers`
module and provides a filtering function, which is a convenient way to remove
specific parts of the PDB file if needed. Alternatively, the user can change
the :const:`~ProtoCaller.Ensemble.Protein.Protein.pdb_obj` manually using
the :mod:`~ProtoCaller.IO` module's capabilities. In this case the PDB file
has to be overwritten after any changes to the PDB object.

.. automodule:: ProtoCaller.Ensemble.Protein
    :members:
    :undoc-members:
    :show-inheritance:
