ProtoCaller
===========

About
-----


Fully automates high-throughput relative protein-ligand binding free energy calculations

This tool uses a variety of tools to automate the free energy calculation process.


Installation
------------

This code is dependent on the BioSimSpace package. If you don't yet have it, you can install it from the git repository:
https://github.com/michellab/BioSimSpace.

Alternatively, you can run:

.. code-block:: bash

    source ./install.sh

for an automatic download of BioSimSpace from the official server and subsequent installation of both BioSimSpace and ProtoCaller.

If you already have installed BioSimSpace, you need to first activate the conda environment and then run setup.py. E.g. if BioSimSpace is installed in ~/biosimspace.app:

.. code-block:: bash

    source activate ~/biosimspace.app
    python setup.py install

Be patient during installation! ProtoCaller will tell you if it has been successfully installed.


Getting Started
---------------

There are a few examples in the examples/ subfolder which you can run to see how ProtoCaller works. More detailed documentation for the software will be added in the future.