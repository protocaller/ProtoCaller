ProtoCaller
===========

About
-----

ProtoCaller is a Python library which enables controlled automation of relative protein-ligand binding free energy
calculations in GROMACS. ProtoCaller uses a variety of tools to automate the free energy calculation process,
such as: Biopython, BioSimSpace, CHARMM-GUI, (optionally) Modeller, Open Babel, ParmEd, PDB2PQR, pdbfixer, RDKit.

ProtoCaller can be run on both Linux and macOS. Installation is easy and performed through Conda. Please check the
other sections for further information.


Installation
------------

This package is distributed via Conda. To install it, run the following command:

.. code-block:: bash

    conda install -c conda-forge -c omnia -c michellab -c essexlab protocaller

The development version can be installed with this command (use with caution):

.. code-block:: bash

    conda install -c conda-forge -c omnia -c michellab -c essexlab/label/dev protocaller

The newest version (1.1.1) is highly recommended due to some compatibility issues with the most recent versions of
BioSimSpace and Modeller.

IMPORTANT: Please note that some recent (10.04.2020) changes to the Protein Data Bank broke compatibility with
ProtoCaller. As there seem to be some unresolved issues with the database, there is a temporary hotfix for
ProtoCaller that will be subject to change as the matter progresses further. Until then, please use the current
development version:

.. code-block:: bash

    conda install -c conda-forge -c omnia -c michellab -c essexlab/label/dev protocaller=1.1.2


Getting Started
---------------

Full docstring documentation can be found `here <https://protocaller.readthedocs.io/en/latest/ProtoCaller.html>`_.
There are also a few `examples <https://protocaller.readthedocs.io/en/latest/Examples.html>`_ which you can run to
see how ProtoCaller works.