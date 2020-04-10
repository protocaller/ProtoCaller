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

System preparation constitutes the larger part of ProtoCaller. Central to system preparation is the
:class:`~ProtoCaller.Ensemble.Ensemble` class, which is partitioned into a :class:`~ProtoCaller.Ensemble.Protein`
and a :class:`~ProtoCaller.Ensemble.PerturbationList`. A :class:`~ProtoCaller.Ensemble.PerturbationList` consists of
multiple :class:`~ProtoCaller.Ensemble.Perturbation` s and each :class:`~ProtoCaller.Ensemble.Perturbation` contains two
:class:`~ProtoCaller.Ensemble.Ligand` s. Each :class:`~ProtoCaller.Ensemble.Protein` and
:class:`~ProtoCaller.Ensemble.Ligand` depends on the force field :class:`~ProtoCaller.Parametrise.Params`. Finally,
a :class:`~ProtoCaller.Ensemble.Protein` also contains a :class:`~ProtoCaller.IO.PDB.PDB` object, which is divided into
:class:`~ProtoCaller.IO.PDB.Chain.Chain` s. :class:`~ProtoCaller.IO.PDB.Chain.Chain` s are split into
:class:`~ProtoCaller.IO.PDB.Residue.Residue` s and :class:`~ProtoCaller.IO.PDB.Residue.Residue` s into
:class:`~ProtoCaller.IO.PDB.Atom.Atom` s, containing all the important information from a PDB file. Sometimes a
:class:`~ProtoCaller.IO.PDB.PDB` also contains one or more :class:`~ProtoCaller.IO.PDB.Missing.MissingResidue` s and/or
:class:`~ProtoCaller.IO.PDB.Missing.MissingAtoms`. The whole dependency tree can be visualised here:

.. graphviz::
    :align: center

    digraph A {
        Ens [label = "Ensemble", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.Ensemble.html#ProtoCaller.Ensemble.Ensemble"];
        Lig [label = "Ligand", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.Ensemble.Ligand.html#ProtoCaller.Ensemble.Ligand.Ligand"]
        PertList [label = "PerturbationList", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.Ensemble.PerturbationList.html#ProtoCaller.Ensemble.PerturbationList.PerturbationList"]
        Pert [label = "Perturbation", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.Ensemble.Perturbation.html#ProtoCaller.Ensemble.Perturbation.Perturbation"]
        Par [label = "Params", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.Parametrise.html#ProtoCaller.Parametrise.Params"]
        Pro [label = "Protein", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.Ensemble.Protein.html#ProtoCaller.Ensemble.Protein.Protein"]
        PDB [label = "PDB", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.IO.PDB.html#ProtoCaller.IO.PDB.PDB"]
        Chain [label = "Chain", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.IO.PDB.Chain.html#ProtoCaller.IO.PDB.Chain.Chain"]
        Residue [label = "Residue", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.IO.PDB.Residue.html#ProtoCaller.IO.PDB.Residue.Residue"]
        Atom [label = "Atom", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.IO.PDB.Atom.html#ProtoCaller.IO.PDB.Atom.Atom"]
        MisRes [label = "MissingResidue", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.IO.PDB.Missing.html#ProtoCaller.IO.PDB.Missing.MissingResidue"]
        MisAtom [label = "MissingAtoms", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.IO.PDB.Missing.html#ProtoCaller.IO.PDB.Missing.MissingAtoms"]

        Ens -> PertList;
        Ens -> Pro;
        PertList -> Pert;
        Pert -> Lig;
        Pert -> Lig;
        Pro -> Par;
        Pro -> PDB;
        Lig -> Par;
        PDB -> Chain;
        PDB -> MisRes;
        PDB -> MisAtom;
        Chain -> Residue;
        Residue -> Atom;
    }

|
Afterwards, it is possible to run the generated files through some simple wrappers. Simulations in GROMACS are done
through the :class:`~ProtoCaller.Simulation.RunGMX` class, which consists of several
:class:`~ProtoCaller.Protocol.Protocol` s which describe the different stages of the free energy simulation,
e.g. minimisation, equilibration and production.

.. graphviz::
    :align: center

    digraph B {
        rankdir = "LR";

        RunGMX [label = "RunGMX", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.Simulation.html#ProtoCaller.Simulation.RunGMX"]
        Protocol [label = "Protocol", href = "https://protocaller.readthedocs.io/en/latest/ProtoCaller.Protocol.html#ProtoCaller.Protocol.Protocol"]

        RunGMX -> Protocol
    }

|
For more information you can look at the `full docstring documentation <https://protocaller.readthedocs.io/en/latest/ProtoCaller.html>`_.
There are also a few `examples <https://protocaller.readthedocs.io/en/latest/Examples.html>`_ which you can run to
see how ProtoCaller works.