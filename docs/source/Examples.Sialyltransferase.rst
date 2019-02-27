Sialyltransferase
=================

This is a more advanced example, where the PDB file is not as straightforward
to use out of the box and the ligands are tougher to handle. Beware that this
script may take a long time to execute, due to the size of the ligands.

NOTE: You will need a Modeller license to run this example, since there are
missing residues in the PDB file.

.. literalinclude:: ../../examples/Sialyltransferase.py
    :language: python

We can see that ProtoCaller provides us with an opportunity to edit our system
inside Python in a scenario where an ordinary black box approach would fail.