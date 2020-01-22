ProtoCaller.Ensemble.Ligand module
==================================

The :class:`~ProtoCaller.Ensemble.Ligand.Ligand` class keeps track of ligand
instantiation, minimisation, protonation and parametrisation. It consists of
an RDKit object and, when parametrised, a link to the parameter files.

This means that this class has to be provided with compelte atom and bond
information. This can be done using: SMILES strings, InChI strings, files
which keep track of bond connectivity (such as SDF, MOL, MOL2) or an RDKit
object. It is also highly advisable to not use SMILES/InChI strings when
possible, since these provide no 3D information which is crucial if there
are multiple conformational minima. Different conformations will lead to
different force field charges and therefore non-reproducible results.
Moreover, 1D identifiers are usually not very good at providing protonation
state/tautomer information which can be important, depending on the system.

The user can also instantiate the ligand with parameter coordinate and topology
files but still needs to provide some bond connectivity information. In this
case SMILES/InChI strings are not an issue.

.. automodule:: ProtoCaller.Ensemble.Ligand
    :members:
    :undoc-members:
    :show-inheritance:
