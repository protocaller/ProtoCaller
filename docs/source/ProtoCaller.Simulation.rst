ProtoCaller.Simulation package
==============================

This package automates free energy simulations based on the input files
created from :mod:`~ProtoCaller.Ensemble`. Here the user can either run
sequential simulations in GROMACS using :class:`~ProtoCaller.Simulation.\
GMXSerialRuns` or REST2 simulations using :class:`~ProtoCaller.Simulation.\
GMX_REST_FEP_Runs`. The latter is of particular interest, since GROMACS does
not provide native support for it. Hence, one needs to use GROMACS patched
with PLUMED in order to use this functionality.

It is likely that this module will be subject to changes. Support for free
energy calculations in other MD engines will be added in the future.

Module contents
---------------

.. automodule:: ProtoCaller.Simulation
    :members:
    :undoc-members:
    :show-inheritance:
