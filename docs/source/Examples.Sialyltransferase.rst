Sialyltransferase
=================

This is a more advanced example, where the PDB file is not as straightforward
to use out of the box and the ligands are tougher to handle. Beware that this
script may take a long time to execute, due to the size of the ligands.

Here we input the ligands as stereospecific SMILES strings. It is important to
make sure the stereochemistry is right, since one wrong (or unspecified)
stereocentre can result in a very wrong structure. Always check your minimised
input structures before you start a simulation, no matter how confident you are
in your input files or ProtoCaller! Alternatively, for molecules with many
stereocentres it might be easier to input the ligand as a file in a supported
input format (e.g. SDF). For more information on how to do that, please refer
to the documentation.

NOTE: You will need a Modeller license to run this example, since there are
missing residues in the PDB file.

.. literalinclude:: ../../examples/Sialyltransferase_1.py
    :language: python

Alternatively, we can use a custom PDB file and ligand files to avoid some
troubles with the preparation. In this case we don't need Modeller and we
get better alignment. Note that in this case the original files may be
overwritten by ProtoCaller during the workflow.

.. literalinclude:: ../../examples/Sialyltransferase_2.py
    :language: python