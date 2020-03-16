ProtoCaller.Wrappers.charmmguiwrapper module
============================================

This module accesses CHARMM-GUI (http://charmm-gui.org) using an external headless
driver which simulates a web browser and automatically performs clicks on the
webserver. This means that the user needs a common web browser (firefox or chrome)
and in some cases (depending on the operating system) an extra chromedriver which
can be installed from your favourite package manager.

This wrapper can upload a PDB file to the server and add missing residues using
CHARMM-GUI's defaults. Any PDB information lost by the CHARMM-GUI wrapper is
automatically restored by ProtoCaller. If this procedure doesn't work, the user
should either go through the procedure manually, or resort to using other modelling
codes. In this case, the user should make sure that no PDB file information is lost.

There are also wrappers for CHARMM-GUI's PDB Reader and Ligand Reader which are
gateways for parametrising proteins and ligands using the CHARMM force field.
However, ProtoCaller can not yet handle this force field and as such these
routines are not yet very useful.

.. automodule:: ProtoCaller.Wrappers.charmmguiwrapper
    :members:
    :undoc-members:
    :show-inheritance:
