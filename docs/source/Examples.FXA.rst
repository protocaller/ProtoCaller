Factor Xa
==========

This is another example where we demonstrate how to load already
parametrised ligands into memory. Initialising a ligand from already
generated GROMACS coordinate and topology files is straightforward,
but in order to do this, we need bond connectivity information in the
shape of a SMILES string or a ligand file, such as an SDF file.

After that we set up a number of perturbations and only keep the
crystal waters and the calcium cation. Modelling the protein can be
challenging and you will need a Modeller license to complete this
example.

Finally we create perturbed structures where we scale the dummy atom
bond lengths by a half in order to improve phase space overlap while
running the simulation.


.. literalinclude:: ../../examples/FXA.py
    :language: python