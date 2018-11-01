from collections.abc import Iterable as _Iterable

import BioSimSpace as _BSS

class Protocol:
    def __init__(self, use_preset=None, **kwargs):
        #integrators
        self.integrator = _BSS.Gateway.String(help="Which integrator to use",
                                              allowed=["leapfrog", "velocity_verlet", "steep", "l-bfgs", "stochastic"])
        self.timestep = _BSS.Gateway.Float(help="What timestep to use in ps")
        self.n_steps = _BSS.Gateway.Integer(help="How many timesteps to run the simulation for", minimum=0)

        #output options
        self.skip_positions = _BSS.Gateway.Integer(help="How many timesteps to skip when writing the positions",
                                                   minimum=0)
        self.skip_velocities = _BSS.Gateway.Integer(help="How many timesteps to skip when writing the velocities",
                                                    minimum=0)
        self.skip_forces = _BSS.Gateway.Integer(help="How many timesteps to skip when writing the forces",
                                                minimum=0)
        self.skip_energies = _BSS.Gateway.Integer(help="How many timesteps to skip when writing the energies",
                                                  minimum=0)

        #neigbours
        self.periodic_boundary_conditions = _BSS.Gateway.String(
            help="Whether and how to implement periodic boundary conditions",
            allowed=["3d", "2d", "no"])
        self.neighbour_cutoff = _BSS.Gateway.Float(help="Cutoff in nm for neighbour list",
                                                   minimum=-1)
        self.neighbour_frequency = _BSS.Gateway.Integer(help="How often to update the neighbour list in timesteps",
                                                        minimum=-1)

        #electrostatics options
        self.coulomb_type = _BSS.Gateway.String(help="What type of Coulomb interactions to use",
                                                allowed=["cutoff", "ewald", "pme"])
        self.pme_order = _BSS.Gateway.Integer(help="What order of PME to use. Ignored if coulomb_type is not PME.",
                                              allowed=[4, 6, 8, 10])
        self.coulomb_cutoff = _BSS.Gateway.Float(help="Cutoff in nm for electrostatics",
                                                 minimum=0)

        #van der Waals options
        self.vdw_type = _BSS.Gateway.String(help="What type of van der Waals interactions to use",
                                            allowed=["cutoff", "pme"])
        self.vdw_corr = _BSS.Gateway.String(help="What type of long-distance correction to apply",
                                            allowed=["no", "energy", "energy_pressure"])
        self.vdw_cutoff = _BSS.Gateway.Float(help="Cutoff in nm for van der Waals interactions",
                                             minimum=0)

        #temperature options
        self.thermostat = _BSS.Gateway.String(help="What thermostat to use",
                                              allowed=["no", "berendsen", "nose-hoover", "andersen"])
        self.temp_frequency = _BSS.Gateway.Integer(help="Thermostat friction coefficient / collision frequency in THz.",
                                                   minimum=-1)
        self.temp_time_const = _BSS.Gateway.Float(
            help="Time constant for thermostat coupling in ps. -1 means no coupling",
            minimum=-1)
        self.temperature = _BSS.Gateway.Integer(help="Simulation temperature",
                                                minimum=0)

        #pressure options
        self.barostat = _BSS.Gateway.String(help="What barostat to use",
                                            allowed=["no", "berendsen", "parrinello-rahman"])
        self.pres_frequency = _BSS.Gateway.Integer(help="Barostat friction coefficient / collision frequency in THz.",
                                                   minimum=-1)
        self.pres_time_const = _BSS.Gateway.Float(
            help="Time constant for barostat coupling in ps. -1 means no coupling",
            minimum=-1)
        self.pressure = _BSS.Gateway.Float(help="Simulation pressure",
                                           minimum=0.0)

        #initial velocity options
        self.random_velocities = _BSS.Gateway.Boolean(help="Whether to generate random velocities")
        self.random_velocities_temperature = _BSS.Gateway.Integer(help="Temperature to sample velocities from",
                                                                  minimum=0)
        self.random_velocities_seed = _BSS.Gateway.Integer(help="Seed for random velocity sampling. -1 is random seed",
                                                           minimum=-1)

        #constraint options
        self.constraint = _BSS.Gateway.String(help="Which constraints to apply",
                                              allowed=["no", "h_bonds", "h_angles", "all_bonds", "all_angles"])
        self.constraint_type = _BSS.Gateway.String(help="Which constraint algorithm to use",
                                                   allowed=["lincs", "shake"])

        #free energy options
        self.free_energy = _BSS.Gateway.Boolean(help="Whether this is a free energy calculation")
        self.current_lambda = _BSS.Gateway.Integer(help="The current lambda. Indices start from 0",
                                                   minimum=0)
        self.coulomb_lambdas = []
        self.vdw_lambdas = []
        self.bonded_lambdas = []
        self.restraint_lambdas = []
        self.mass_lambdas = []
        self.temperature_lambdas = []
        self.write_derivatives = _BSS.Gateway.Boolean(help="Whether to write dH/dλ")

        use_preset = use_preset.strip().lower() if use_preset is not None else ""
        if use_preset == "default":
            self._generateGenericParams()
        elif use_preset == "minimisation":
            self._generateMinimisationParams()
        elif use_preset == "equilibration_nvt":
            self._generateNVTEquilibrationParams()
        elif use_preset == "equilibration_npt":
            self._generateNPTEquilibrationParams()
        elif use_preset == "production":
            self._generateProductionParams()

        for name, value in kwargs.items():
            self.__setattr__(name, value)

    def __getattr__(self, name):
        try:
            val = self.__dict__[name]
            if "Gateway" in str(type(val)):
                return val.getValue()
            else:
                return val
        except:
            raise AttributeError("No attribute {} found. You need to set it first.".format(name))

    def __setattr__(self, name, value):
        if isinstance(value, str):
            value = value.strip().lower()
        if name not in self.__dict__.keys():
            self.__dict__[name] = value
        else:
            val_old = self.__dict__[name]
            if "Gateway" not in str(type(val_old)):
                if not isinstance(value, type(val_old)):
                    raise TypeError("Type of argument ({}) does not agree with requested type ({}).".format(
                        str(type(value), str(type(val_old)))))
                else:
                    self.__dict__[name] = value
            else:
                self.__dict__[name].setValue(value)

    def write(self, engine, filebase="protocol"):
        engine = engine.strip().upper()
        if engine == "GROMACS":
            return self._writeToGROMACS(filebase)

    def _generateGenericParams(self):
        self.timestep = 0.002
        self.skip_positions = 500
        self.skip_velocities = 500
        self.skip_forces = 0
        self.skip_energies = 500
        self.periodic_boundary_conditions = "3d"
        self.neighbour_cutoff = 1.2
        self.neighbour_frequency = 20
        self.coulomb_type = "pme"
        self.coulomb_cutoff = 1.2
        self.vdw_type = "cutoff"
        self.vdw_cutoff = 1.2
        self.vdw_corr = "energy_pressure"
        self.constraint = "h_bonds"
        self.constraint_type = "lincs"

    def _generateMinimisationParams(self):
        self._generateGenericParams()
        self.integrator = "steep"
        self.n_steps = 5000
        self.thermostat = "no"
        self.barostat = "no"
        self.random_velocities = False

    def _generateNVTEquilibrationParams(self):
        self._generateGenericParams()
        self.integrator = "stochastic"
        self.n_steps = 25000
        self.thermostat = "berendsen"
        self.temp_time_const = 1
        self.temperature = 298
        self.barostat = "no"
        self.random_velocities = True
        self.random_velocities_temperature = 298

    def _generateNPTEquilibrationParams(self):
        self._generateGenericParams()
        self.integrator = "stochastic"
        self.n_steps = 25000
        self.thermostat = "berendsen"
        self.temp_time_const = 1
        self.temperature = 298
        self.barostat = "berendsen"
        self.pres_time_const = 1
        self.random_velocities = False

    def _generateProductionParams(self):
        self._generateNPTEquilibrationParams()
        self.integrator = "velocity_verlet"
        self.n_steps = 2500000
        self.thermostat = "nose-hoover"
        self.barostat = "parrinello-rahman"

    def _writeToGROMACS(self, filebase):
        name_dict = {
            "integrator": "integrator",
            "timestep": "dt",
            "n_steps": "nsteps",

            "skip_positions": "nstxout",
            "skip_velocities": "nstvout",
            "skip_forces": "nstfout",
            "skip_energies": "nstenergy",

            "periodic_boundary_conditions": "pbc",
            "neighbour_cutoff": "rlist",
            "neighbour_frequency": "nstlist",

            "coulomb_type": "coulombtype",
            "pme_order": "pme-order",
            "coulomb_cutoff": "rcoulomb",

            "vdw_type": "vdwtype",
            "vdw_corr": "DispCorr",
            "vdw_cutoff": "rvdw",

            "thermostat": "tcoupl",
            "temp_frequency": "nsttcouple",
            "temp_time_const": "tau-t",
            "temperature": "ref-t",

            "barostat": "pcoupl",
            "pres_frequency": "nstpcouple",
            "pres_time_const": "tau-p",
            "pressure": "ref-p",

            "random_velocities": "gen-vel",
            "random_velocities_temperature": "gen-temp",
            "random_velocities_seed": "gen-seed",

            "constraint": "constraints",
            "constraint_type": "constraint-algorithm",

            "free_energy": "free-energy",
            "current_lambda": "init-lambda",
            "coulomb_lambdas": "coul-lambdas",
            "vdw_lambdas": "vdw-lambdas",
            "bonded_lambdas": "bonded-lambdas",
            "restraint_lambdas": "restraint-lambdas",
            "mass_lambdas": "mass-lambdas",
            "temperature_lambdas": "temperature-lambdas",
            "write_derivatives": "dhdl-derivatives",
        }

        value_dict = {}

        value_dict["integrator"] = {
            "leapfrog": "md",
            "velocity_verlet": "md-vv",
            "steep": "steep",
            "l-bfgs": "l-bfgs",
            "stochastic": "sd",
        }

        value_dict["periodic_boundary_conditions"] = {
            "3d": "xyz",
            "2d": "xy",
            "no": "no"
        }

        value_dict["coulomb_type"] = {
            "cutoff": "Cut-off",
            "ewald": "Ewald",
            "pme": "PME",
        }

        value_dict["vdw_type"] = {
            "cutoff": "Cut-off",
            "pme": "PME",
        }

        value_dict["vdw_corr"] = {
            "no": "no",
            "energy": "Ener",
            "energy_pressure": "EnerPres",
        }

        value_dict["thermostat"] = {
            "no": "no",
            "berendsen": "berendsen",
            "nose-hoover": "nose-hoover",
            "andersen": "andersen",
        }

        value_dict["barostat"] = {
            "no": "no",
            "berendsen": "berendsen",
            "parrinello-rahman": "Parrinello-Rahman",
        }

        value_dict["constraint"] = {
            "no": "none",
            "h_bonds": "h-bonds",
            "h_angles": "h-angles",
            "all_bonds": "all-bonds",
            "all_angles": "all-angles",
        }

        value_dict["constraint_type"] = {
            "lincs": "LINCS",
            "shake": "SHAKE",
        }

        filename = filebase + ".mdp"
        with open(filename, "w") as file:
            for name, value in sorted(self.__dict__.items()):
                if name in name_dict.keys():
                    name_str = name_dict[name]
                elif name[0] == "_":
                    name_str = name[1:]
                else:
                    name_str = name

                if isinstance(value, _Iterable) and not isinstance(value, str):
                    value_str = "\t".join(value)
                elif isinstance(value, _BSS.Gateway.Boolean):
                    if value.getValue is True:
                        value_str = "yes"
                    elif value.getValue is False:
                        value_str = "no"
                    else:
                        value_str = ""
                elif isinstance(value, _BSS.Gateway.String):
                    if name in value_dict.keys() and value.getValue() in value_dict[name].keys():
                        value_str = value_dict[name][value.getValue()]
                    else:
                        value_str = value.getValue()
                elif "Gateway" in str(type(value)):
                    value_str = str(value.getValue())
                elif value not in [[], None, ""]:
                    value_str = str(value)
                else:
                    continue

                if name_str and value_str:
                    file.write("{:<30} = {}\n".format(name_str, value_str))
        return filename













































'''
def MDPgen(Lambda, filename="mdp.mdp", temperature = True, pressure = False, *Lambdalists, **kwargs):
    RUN_DEFAULT = {
        "integrator": "sd",
        "tinit": "0",
        "dt": "0.002",
        "nsteps": "2500000",
        "nstcomm": "100",
    }

    OUT_DEFAULT = {
        "nstxout": "500",
        "nstvout": "500",
        "nstfout": "0",
        "nstlog": "500",
        "nstenergy": "500",
        "nstxout-compressed": "0",
    }

    NEIGHBOUR_DEFAULT = {
        "cutoff-scheme": "verlet",
        "nstlist": "20",
        "ns-type": "grid",
        "pbc": "xyz",
        "rlist": "1.2",
    }

    COULOMB_DEFAULT = {
        "coulombtype": "PME",
        "rcoulomb": "1.2",
    }

    VDW_DEFAULT = {
        "vdwtype": "cutoff",
        "vdw-modifier": "potential-switch",
        "rvdw-switch": "1.0",
        "rvdw": "1.2",
    }

    COULOMB_LR_DEFAULT = {
        "fourierspacing": "0.12",
        "pme-order": "6",
        "ewald-rtol": "1e-06",
        "epsilon-surface": "0",
    }

    VDW_LR_DEFAULT = {
        "DispCorr": "EnerPres",
    }

    TEMPERATURE_DEFAULT = {
        "tcoupl": "berendsen",
        "tc-grps": "system",
        "tau-t": "1.0",
        "ref-t": "298",
    }

    PRESSURE_DEFAULT = {
        "pcoupl": "Parrinello-Rahman",
        "tau-p": "1.0",
        "compressibility": "4.5e-05",
        "ref-p": "1.0",
    }

    FE_DEFAULT = {
        "free-energy": "yes",
        "init-lambda-state": "0",
        "delta-lambda": "0",
        "calc-lambda-neighbors": "1",
    }

    SC_DEFAULT = {
        "sc-alpha": "0.5",
        "sc-coul": "no",
        "sc-power": "1",
        "sc-sigma": "0.3",
        "nstdhdl": "10",
    }

    CONSTRAINTS_DEFAULT = {
        "constraints": "h-bonds",
        "constraint-algorithm": "lincs",
        "lincs-order": "12",
    }

    MISC_DEFAULT = {
        "continuation": "no",
    }

    USERDICT = {}

    DICTLIST = [USERDICT, RUN_DEFAULT, OUT_DEFAULT, NEIGHBOUR_DEFAULT, COULOMB_DEFAULT, VDW_DEFAULT, COULOMB_LR_DEFAULT, VDW_LR_DEFAULT, TEMPERATURE_DEFAULT,
        PRESSURE_DEFAULT, FE_DEFAULT, SC_DEFAULT, CONSTRAINTS_DEFAULT, MISC_DEFAULT]

    FE_DEFAULT["init-lambda-state"] = "%d" % Lambda

    for key, value in kwargs.iteritems():
        found = False
        for dict in DICTLIST:
            if(found): break
            if(key in dict):
                found = True
                dict[key] = value
        if(found == False):
            USERDICT[key] = value

    with open(filename, 'w') as file:
        if(temperature == False):
            TEMPERATURE_DEFAULT = {"tcoupl": "no"}
        if(pressure == False):
            PRESSURE_DEFAULT = {"pcoupl": "no"}

        DICTLIST = [USERDICT, RUN_DEFAULT, OUT_DEFAULT, NEIGHBOUR_DEFAULT, COULOMB_DEFAULT, VDW_DEFAULT, COULOMB_LR_DEFAULT, VDW_LR_DEFAULT, TEMPERATURE_DEFAULT,
            PRESSURE_DEFAULT, FE_DEFAULT, SC_DEFAULT, CONSTRAINTS_DEFAULT, MISC_DEFAULT]

        for i, dict in enumerate(DICTLIST):
            for key, value in dict.iteritems():
                if(value):
                    file.write("{0: <45} = {1}\n".format(key, value))
            if(i == 10):
                for Lambdalist in Lambdalists:
                    if(len(Lambdalist) > 1):
                        file.write("{0: <45} = ".format(Lambdalist[0]))
                        for value in Lambdalist[1:]:
                            file.write("%s  " % value)
                        file.write("\n")
            file.write("\n")
class MDPpreset:
    params = {}
    Lambda = 0;
    vdw_list = ["vdw-lambdas"]
    coulomb_list = ["coul-lambdas"]
    bonded_list = ["bonded-lambdas"]
    restraint_list = ["restraint-lambdas"]
    filename = "mdp.mdp"
    temperature = True
    pressure = False

    def WriteMDP(self):
        MDPgen(self.Lambda, self.filename, self.temperature, self.pressure, self.vdw_list, self.coulomb_list, self.bonded_list, self.restraint_list, **self.params)

steep = MDPpreset()
steep.params = {
    "integrator": "steep",
    "tinit": "",
    "dt": "",
    "nsteps": "5000",
    "nstcomm": "",
    "emtol": "100",
    "emstep": "0.01",
    "niter": "20",
    "nstlog": "1",
    "nstenergy": "1",
    "nstlist": "1",
}
steep.temperature = False

steep_vac = MDPpreset()
steep_vac.params = {
    "integrator": "steep",
    "tinit": "",
    "dt": "",
    "nsteps": "5000",
    "nstcomm": "",
    "emtol": "100",
    "emstep": "0.01",
    "niter": "20",
    "nstlog": "1",
    "nstenergy": "1",
    "nstlist": "1",
    "constraints": "none",
    "cutoff-scheme": "group",
    "nstlist": "0",
    "ns-type": "simple",
    "pbc": "no",
    "rlist": "0",
    "coulombtype": "Cut-off",
    "rcoulomb": "0",
    "vdw-modifier": "None",
    "rvdw": "0",
    "DispCorr": "no",
}
steep_vac.temperature = False

l_bfgs = MDPpreset()
l_bfgs.params = {
    "integrator": "l-bfgs",
    "tinit": "",
    "dt": "",
    "nsteps": "5000",
    "nstcomm": "",
    "emtol": "100",
    "emstep": "0.01",
    "niter": "20",
    "nbfgscorr": "10",
    "define": "-DFLEXIBLE",
    "nstlog": "1",
    "nstenergy": "1",
    "nstlist": "1",
}
l_bfgs.temperature = False

nvt = MDPpreset()
nvt.params = {
    "gen-vel": "yes",
    "gen-temp": "300",
    "gen-seed": "-1",
}

npt = MDPpreset()
npt.pressure = True

npt_vac = MDPpreset()
npt_vac.params = {
    "cutoff-scheme": "group",
    "nstlist": "0",
    "ns-type": "simple",
    "pbc": "no",
    "rlist": "0",
    "coulombtype": "Cut-off",
    "rcoulomb": "0",
    "vdw-modifier": "None",
    "rvdw": "0",
    "DispCorr": "no",
}
npt_vac.pressure = True

class Subdir():
    def __init__(self, dirname, copydirname=None, overwrite=False):
        self.dirname = dirname
        self.copydirname = copydirname
        self.overwrite = overwrite
    def __enter__(self):
        self.workdirname = os.getcwd()
        if(self.overwrite and os.path.exists(self.dirname)):
            shutil.rmtree("%s/%s" % (self.workdirname, self.dirname))
        if(self.overwrite and self.copydirname and os.path.exists(self.copydirname)):
            shutil.copytree("%s/%s" % (self.workdirname, self.copydirname), "%s/%s" % (self.workdirname, self.dirname))
        elif(not os.path.exists(self.dirname)):
            os.mkdir(self.dirname)
        os.chdir(self.dirname)
    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.workdirname)

def WriteMDPs():
    global preset_list, nsteps_list, checkpoint_list, dt_list, filename_list, coulomb_list, vdw_list, bonded_list
    if(not os.path.exists("MDP")):
        with Subdir("MDP"):
            for Lambda in range(len(coulomb_list)):
                for i, preset in enumerate(preset_list):
                    preset.Lambda = Lambda
                    preset.params["dt"] = dt_list[i]
                    preset.params["nsteps"] = nsteps_list[i]
                    preset.coulomb_list = ["coul-lambdas"] + coulomb_list
                    preset.vdw_list = ["vdw-lambdas"] + vdw_list
                    preset.bonded_list = ["bonded-lambdas"] + bonded_list
                    preset.filename = "%s_%d.mdp" % (filename_list[i], Lambda)
                    preset.params["continuation"] = "yes" if checkpoint_list[i] else "no"
                    preset.WriteMDP()

def WriteLambdas(Lambdalist):
    global preset_list, nsteps_list, checkpoint_list, dt_list, filename_list, coulomb_list, vdw_list, bonded_list, WORKDIR, FILENAME
    with Subdir("Lambdas"):
        for Lambda in Lambdalist:
            with Subdir("Lambda_%s" % Lambda):
                for i, name in enumerate(filename_list):
                    with Subdir(name.upper()):
                        if(i == 0):
                            subprocess.check_call(["gmx grompp -f {0}/MDP/{1}_{2}.mdp -c {0}/{5}.gro -p {0}/{5}.top -o {1}_{2}.tpr -maxwarn 10"
                            .format(WORKDIR, name, Lambda, filename_list[i-1].upper(), filename_list[i-1], FILENAME)], shell="True")
                        elif(checkpoint_list[i]):
                            subprocess.check_call(["gmx grompp -f {0}/MDP/{1}_{2}.mdp -c ../{3}/{4}_{2}.tpr -p {0}/{5}.top -t ../{3}/{4}_{2}.cpt -o {1}_{2}.tpr -maxwarn 10"
                            .format(WORKDIR, name, Lambda, filename_list[i-1].upper(), filename_list[i-1], FILENAME)], shell="True")
                        else:
                            subprocess.check_call(["gmx grompp -f {0}/MDP/{1}_{2}.mdp -c ../{3}/{4}_{2}.tpr -p {0}/{5}.top -o {1}_{2}.tpr -maxwarn 10"
                            .format(WORKDIR, name, Lambda, filename_list[i-1].upper(), filename_list[i-1], FILENAME)], shell="True")
                        subprocess.check_call(["gmx mdrun -s {0}_{1}.tpr -deffnm {0}_{1}".format(name, Lambda)], shell=True)

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", default="complex_final", help="name of input .gro and .top files (must be the same name). Default name: \"complex_final\"")
parser.add_argument("-l", "--lambda", type = int, nargs="+", help="lambda value to calculate. Default: calculate for every lambda sequentially")
args = parser.parse_args()

coulomb_list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.6, 0.8, 1.0]
vdw_list =     [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0]
bonded_list =  [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
preset_list =  [steep, nvt, nvt, npt, npt]
nsteps_list =  [5000, 10000, 500000, 500000, 2500000]
checkpoint_list = [False, False, True, True, True] #whether or not to pass checkpoint files to gmx grompp. Coincides with continue parameter
dt_list = ["", 0.0002, 0.002, 0.002, 0.002]
filename_list = ["steep", "nvt_acc", "nvt", "npt", "production_md"]

#subprocess.check_call(["module load gromacs/5.1.4/intel"], shell=True)
WORKDIR = os.path.dirname(os.path.abspath(__file__))
FILENAME = args.output
LAMBDALIST = [i for i in range(len(coulomb_list))] if vars(args)["lambda"] is None else vars(args)["lambda"]

WriteMDPs()
WriteLambdas(LAMBDALIST)'''