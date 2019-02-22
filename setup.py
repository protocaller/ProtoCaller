from setuptools import setup, find_packages

setup(
    name="ProtoCaller",
    version="0.1dev",
    description="ProtoCaller: Full automation of relative protein-ligand binding free energy simulations.",
    author="Miroslav Suruzhon",
    email="M.Suruzhon@soton.ac.uk",
    packages=find_packages(),
    include_package_data=True,
)

try:
    import BioSimSpace
    import Sire.Base
    bin_dir = Sire.Base.getBinDir()
except:
    raise ImportError("Module BioSimSpace not found. Please install the latest development version of BioSimSpace "
                      "from https://github.com/michellab/BioSimSpace")

import os, subprocess, sys
if "install" in sys.argv:
    stdout = open("setup.out", "w")
    stderr = open("setup.err", "w")

    command = "%s/conda clean --yes --all" % bin_dir
    subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

    # fixes an issue with the openbabel used in BioSimSpace
    command = "%s/conda remove -q -y openbabel" % bin_dir
    subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

    command = "%s/conda env update -q -f requirements.yml --prune" % bin_dir
    subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

    command = "%s/conda clean --yes --all" % bin_dir
    subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

    stdout.close()
    stderr.close()

    try:
        import ProtoCaller
        os.remove("setup.out")
        os.remove("setup.err")
        print("ProtoCaller has been successfully installed")
    except:
        print("Installation failed. Please check output in 'setup.out' and 'setup.err'")



