ProtoCaller
===========

About
-----


Fully automates high-throughput relative protein-ligand binding free energy calculations

This tool uses a variety of tools to automate the free energy calculation process.


Installation
------------

This package is distributed via Conda. To install it, run the following commands:

.. code-block:: bash

    conda install -c michellab/label/dev -c rdkit -c conda-forge -c omnia -c bioconda -c essexlab protocaller

ProtoCaller depends on BioSimSpace. If it is run on MacOS, there is a known
issue with BioSimSpace which means that instead of python one has to run
the sire_python interpreter instead.


Getting Started
---------------

There are a few examples which you can find either in the examples/ subfolder or in the "Examples" section of the documentation which you can run to see how ProtoCaller works.