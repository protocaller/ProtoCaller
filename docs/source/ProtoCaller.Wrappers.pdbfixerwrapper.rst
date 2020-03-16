ProtoCaller.Wrappers.pdbfixerwrapper module
===========================================

This is a light wrapper around pdbfixer. Here the user can add missing residues
and/or atoms to their structure and ProtoCaller keeps track of any PDB information
lost in the process. It is also possible to replace nonstandard residues with
their standard counterparts. This is highly recommended, since in most cases there
are no reliable parameters for modified amino acid residues. If the user runs
an external modelling code, they need to make sure that no missing atom/disulfide
bond information is lost from the PDB file before passing it to ProtoCaller.

.. automodule:: ProtoCaller.Wrappers.pdbfixerwrapper
    :members:
    :undoc-members:
    :show-inheritance:
