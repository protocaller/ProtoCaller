Dihydrofolate Reductase
=======================

This example is similar to the factor Xa example in that we set up
a number of perturbations from custom SDF files whilst keeping
the parametrised files.

In this example we show how to prepare the same ligands with different
crystal structures. In the first case we use the original ligand
orientation, while in the second one we use a custom one provided
by an external SDF file that we generated using PyMOL.

We also demonstrate that it is also possible to edit the PDB file
inside ProtoCaller by creating two batches of structures - one from
the altLoc A coordinates and one from the altLoc B coordinates in
the original PDB file. . This functionality can be very useful when
the PDB files are not perfect and one needs to edit them in a
controlled manner.

Finally, DHFR has a cofactor (NADPH) close to the active site and as
such needs modelling. In this case we use one of the parameter sets
included in ProtoCaller from the `AMBER parameter database website
<http://research.bmh.manchester.ac.uk/bryce/amber/>`_. This means
that ProtoCaller can currently handle systems with some of the
most common cofactors. The full list of cofactors supported by
ProtoCaller can be accessed from the global variable
:const:`ProtoCaller.COFACTORNAMES`.

.. literalinclude:: ../../examples/DHFR.py
    :language: python
