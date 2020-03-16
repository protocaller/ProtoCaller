ProtoCaller.Wrappers.modellerwrapper module
===========================================

This is the wrapper responsible for handling Modeller. Here the user can add missing
residues and optionally missing atoms. Note that it is not yet possible to add missing
atoms if there are no missing residues. Additionally, the user can perform an optional
refinement from within Modeller. More information can be found
`here <https://salilab.org/modeller/9v8/manual/node44.html>`_.

Whenever PDB information is lost during the transformation, ProtoCaller automatically
restores it. If the user is running Modeller outside ProtoCaller, they need to make
sure that appropriate information about missing atoms, disulfide bonds and chain IDs
is retained.

Note that in all cases you will need to `install <https://anaconda.org/salilab/modeller>`_
Modeller and have a valid `license <https://salilab.org/modeller/registration.html>`_.

.. automodule:: ProtoCaller.Wrappers.modellerwrapper
    :members:
    :undoc-members:
    :show-inheritance:
