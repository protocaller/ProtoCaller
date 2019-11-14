Performing Free Energy Calculations in GROMACS
==============================================

So far we have shown a number of examples which can hopefully
get you started with ProtoCaller. However, all of the examples
covered system preparation. While system preparation is arguably
the most sophisticated part of the free energy workflow, it
is also useful to streamline the subsequent simulations.

In this tutorial we show an example run script which accepts
the working directory and the stage of the simulation (bound or
solvated leg) as arguments. This means that this script can
be universally used for your structures prepared with
ProtoCaller either on your desktop or on a computing cluster.

The below script defines a number of van der Waals, Coulomb
and bonded lambda values from 0 to 1 (we assume that we are
shrinking the ligand). After that it runs minimisation,
equilibration in the NVT ensemble, equilibration in the NPT
ensemble and final production for 1 ns per lambda window.
All of the parameters and presets can be changed by the user,
including running the lambda windows in parallel. These
extra options can be found in the
:class:`~ProtoCaller.Protocol.Protocol` and
:class:`~ProtoCaller.Simulation.RunGMX` classes.

.. literalinclude:: ../../examples/run.py
    :language: python
