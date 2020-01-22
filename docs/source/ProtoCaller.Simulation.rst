ProtoCaller.Simulation package
==============================

This package automates free energy simulations based on the input files
created from :mod:`~ProtoCaller.Ensemble`. The module currently consists of
a single class: :class:`~ProtoCaller.Simulation.RunGMX`. Through this class,
it is possible to streamline free energy calculations in an accessible way.
One can either instantiate this class with a single GRO and TOP file or
with several topologies and/or positions. In any case the lambda values
have to be passed to the initialiser.

The main function of the method is :func:`~ProtoCaller.Simulation.RunGMX.\
runSimulation`. It enables running simulations based on protocols generated
through :mod:`~ProtoCaller.Protocol`. The different lambda windows can
be run either sequentially or in parallel where the user can optionally
specify OpenMP and MPI threads. If the detected GROMACS version is patched
with PLUMED, the user can also perform replica exchange between adjacent
lambda windows.

The GROMACS binary is detected automatically from the bash aliases "gmx"
and "gmx_mpi". If the used binaries are not satisfactory, the user can
set the environment variable :const:`GROMACSHOME` to the path where the
bin/ directory is located. Alternatively, one can set the
:const:`ProtoCaller.GROMACSEXE` and/or the :const:`ProtoCaller.GROMACSMPIEXE`
variables to the absolute path of the desired binary.

A useful example of this module in action can be found :doc:`here <Examples.Run>`.

Module contents
---------------

.. automodule:: ProtoCaller.Simulation
    :members:
    :undoc-members:
    :show-inheritance:
