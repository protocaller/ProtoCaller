# GROMACS is detected automatically.
# Set the GROMACSHOME variable to a custom location otherwise
# import os
# os.environ["GROMACSHOME"] = "CUSTOM_LOCATION/gromacs-2018.4"

import logging
logging.basicConfig(level=logging.INFO)

import numpy as np
import ProtoCaller.Simulation as Simulation
from ProtoCaller.Utils.fileio import Dir
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mode", type=str,
                    help="Whether to run bound or solvated leg")
parser.add_argument("-w", "--workdir", type=str, help="Path to where "
                    "complex_final.* and/or morph.* are located")
args = parser.parse_args()

workdir = Dir(args.workdir) if args.workdir else Dir(".")

if args.mode == "bound":
    dirname = "Run_Bound"
    gro = "complex_final.gro"
    top = "complex_final.top"
elif args.mode == "solvated":
    dirname = "Run_Solvated"
    gro = "morph.gro"
    top = "morph.top"
else:
    raise ValueError("Please choose either 'bound' or 'solvated'")

# determine the lambda values
vdw_lambdas = [0] * 9 + [round(x, 2) for x in np.linspace(0.00, 0.95, num=26)] \
              + [0.97, 0.98, 0.99, 0.999, 1.00]
bonded_lambdas = vdw_lambdas[:]
coulomb_lambdas = [round(x, 2) for x in np.linspace(0.00, 1.00, num=10)] + \
                  [1] * 30

with workdir:
    # here we show how to pass extra parameters to GROMACS.
    # In this case these are not needed, however.
    extra_params = {"fep-lambdas": [0] * 40, "calc-lambda-neighbors": "-1"}
    # initialise a Simulation object
    run = Simulation.RunGMX(dirname, gro, top, coulomb_lambdas=coulomb_lambdas,
                            vdw_lambdas=vdw_lambdas, bonded_lambdas=bonded_lambdas)
    # run minimisation
    run.runSimulation(name="Minimisation", use_preset="minimisation",
                      n_steps=25000, softcore_coulomb=False, constraint="no",
                      extra_params=extra_params)
    # run NVT equilibration
    run.runSimulation(name="Equilibration_NVT", use_preset="equilibration_nvt",
                      n_steps=50000, timestep=0.001, softcore_coulomb=False,
                      extra_params=extra_params)
    # run NPT equilibration
    run.runSimulation(name="Equilibration_NPT", use_preset="equilibration_npt",
                      n_steps=50000, timestep=0.001, softcore_coulomb=False,
                      extra_params=extra_params)
    # run production in NPT
    run.runSimulation(name="Production", use_preset="production",
                      integrator="stochastic", n_steps=500000, timestep=0.002,
                      softcore_coulomb=False, extra_params=extra_params,
                      skip_positions=2500, skip_energies=2500)

