from collections.abc import Iterable as _Iterable
from collections import OrderedDict
import warnings as _warnings

import BioSimSpace as _BSS

class Protocol:
    """
    A wrapper for different simulation parameters.

    Parameters
    ----------
    use_preset : str, None
        Which default preset to use. One of: "minimisation", "equilibration_nvt", "equilibration_npt", "production"
        "vacuum".
    extra_params : dict
        Used to set custom attributes with illegal Python characters in them.
    kwargs:
        Used to set any of the attributes, custom or not.

    Attributes
    ----------

    integrator : string
        Which integrator to use.
    timestep : float
        What timestep to use in ps.
    n_steps : int
        How many timesteps to run the simulation for.
    skip_positions : int
        How many timesteps to skip when writing the positions.
    skip_velocities : int
        How many timesteps to skip when writing the velocities.
    skip_forces : int
        How many timesteps to skip when writing the forces.
    skip_energies : int
        How many timesteps to skip when writing the energies.
    periodic_boundary_conditions : string
        Whether and how to implement periodic boundary conditions.
    neighbour_cutoff : float
        Cutoff in nm for neighbour list.
    neighbour_frequency : int
        How often to update the neighbour list in timesteps.
    coulomb_type : string
        What type of Coulomb interactions to use.
    pme_order : int
        What order of PME to use. Ignored if coulomb_type is not PME.
    coulomb_cutoff : float
        Cutoff in nm for electrostatics.
    vdw_type : string
        What type of van der Waals interactions to use.
    vdw_corr : string
        What type of long-distance correction to apply.
    vdw_cutoff : float
        Cutoff in nm for van der Waals interactions.
    thermostat : string
        What thermostat to use.
    temp_frequency : int
        Thermostat friction coefficient / collision frequency in THz.
    temp_time_coupling : float
        Time constant for thermostat coupling in ps. -1 means no coupling.
    temperature : int
        Simulation temperature.
    temp_groups : string
        Which parts of the system to heat up.
    barostat : string
        What barostat to use.
    pres_frequency : int
        Barostat friction coefficient / collision frequency in THz.
    pres_time_coupling : float
        Time constant for barostat coupling in ps. -1 means no coupling.
    pressure : float
        Simulation pressure.
    compressibility : float
        System compressibility.
    random_velocities : bool
        Whether to generate random velocities.
    random_velocities_temperature : int
        Temperature to sample velocities from.
    random_velocities_seed : int
        Seed for random velocity sampling. -1 is random seed.
    constraint : string
        Which constraints to apply.
    constraint_type : string
        Which constraint algorithm to use.
    free_energy : bool
        Whether this is a free energy calculation.
    current_lambda : int
        The current lambda. Indices start from 0.
    coulomb_lambdas : [float]
        Initialises coulomb_lambdas.
    vdw_lambdas : [float]
        Initialises vdw_lambdas.
    bonded_lambdas : [float]
        Initialises bonded_lambdas.
    restraint_lambdas : [float]
        Initialises restraint_lambdas.
    mass_lambdas : [float]
        Initialises mass_lambdas.
    temperature_lambdas : [float]
        Initialises temperature_lambdas.
    write_derivatives : bool
        Whether to write dH/dλ.
    softcore_vdw : bool
        Apply soft-core potential to van der Waals interactions.
    softcore_coulomb : bool
        Apply soft-core potential to Coulomb interactions.
    softcore_alpha : float
        Soft-core alpha parameter, default is 0.5.
    softcore_lambda_power : int
        Power for lambda in the soft-core function.
    """
    def __init__(self, use_preset=None, extra_params=None, **kwargs):
        self._attrs = OrderedDict()
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
        self.pme_order = _BSS.Gateway.Integer(help="What order of PME to use. Ignored if coulomb_type is not PME",
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
        self.temp_frequency = _BSS.Gateway.Integer(help="Thermostat friction coefficient / collision frequency in THz",
                                                   minimum=-1)
        self.temp_time_coupling = _BSS.Gateway.Float(
            help="Time constant for thermostat coupling in ps. -1 means no coupling",
            minimum=-1)
        self.temperature = _BSS.Gateway.Integer(help="Simulation temperature",
                                                minimum=0)
        self.temp_groups = _BSS.Gateway.String(help="Which parts of the system to heat up",
                                               allowed=["all"])

        #pressure options
        self.barostat = _BSS.Gateway.String(help="What barostat to use",
                                            allowed=["no", "berendsen", "parrinello-rahman", "mttk"])
        self.pres_frequency = _BSS.Gateway.Integer(help="Barostat friction coefficient / collision frequency in THz",
                                                   minimum=-1)
        self.pres_time_coupling = _BSS.Gateway.Float(
            help="Time constant for barostat coupling in ps. -1 means no coupling",
            minimum=-1)
        self.pressure = _BSS.Gateway.Float(help="Simulation pressure",
                                           minimum=0)
        self.compressibility = _BSS.Gateway.Float(help="System compressibility",
                                                  minimum=0)

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

        #softcore options
        self.softcore_vdw = _BSS.Gateway.Boolean(help="Apply soft-core potential to van der Waals interactions")
        self.softcore_coulomb = _BSS.Gateway.Boolean(help="Apply soft-core potential to Coulomb interactions")
        self.softcore_alpha = _BSS.Gateway.Float(help="Soft-core alpha parameter, default is 0.5")
        self.softcore_lambda_power = _BSS.Gateway.Integer(help="Power for lambda in the soft-core function",
                                                          allowed=[1, 2])

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
        elif use_preset == "vacuum":
            self._generateVacuumParams()
        else:
            _warnings.warn("Invalid preset name: %s" % use_preset)

        extra_params = extra_params if extra_params is not None else {}
        all_kwargs = {**extra_params, **kwargs}
        for name, value in all_kwargs.items():
            self.__setattr__(name, value)

    def __copy__(self):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def __getattr__(self, name):
        if name == "_Protocol__attrs":
            return self._attrs
        elif "__" in name:
            return self.name
        try:
            val = self.__attrs[name]
            if "Gateway" in str(type(val)):
                return val.getValue()
            else:
                return val
        except:
            raise AttributeError("No attribute {} found. You need to set it first.".format(name))

    def __setattr__(self, name, value):
        if isinstance(value, str):
            value = value.strip()

        if "lambdas" in name and value:
            self.free_energy = True

        if name == "_attrs":
            super(Protocol, self).__setattr__(name, value)
            return

        if name not in self.__attrs.keys():
            self.__attrs[name] = value
        else:
            val_old = self.__attrs[name]
            if "Gateway" not in str(type(val_old)):
                if not isinstance(value, type(val_old)):
                    raise TypeError("Type of argument ({}) does not agree with requested type ({}).".format(
                        str(type(value)), str(type(val_old))))
                else:
                    self.__attrs[name] = value
            else:
                self.__attrs[name].setValue(value, name)

    def write(self, engine, filebase="protocol"):
        """
        Writes the protocol to a file for a specific MD engine.

        Parameters
        ----------
        engine : str
            Determines the format of the file. One of: "GROMACS".
        filebase : str, optional
            The base of the output file.

        Returns
        -------
        filename : str
            The absolute path to the output parameter file.
        """
        engine = engine.strip().upper()
        if engine == "GROMACS":
            return self._writeToGROMACS(filebase)

    def _generateGenericParams(self):
        self.timestep = 0.002
        self.skip_positions = 500
        self.skip_velocities = 0
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
        self.softcore_vdw = True
        self.softcore_coulomb = True
        self.softcore_alpha = 0.5
        self.softcore_lambda_power = 1
        self.constraint = "h_bonds"
        self.constraint_type = "lincs"

    def _generateMinimisationParams(self):
        self._generateGenericParams()
        self.skip_positions = 0
        self.integrator = "steep"
        self.n_steps = 5000
        self.thermostat = "no"
        self.barostat = "no"
        self.random_velocities = False

    def _generateNVTEquilibrationParams(self):
        self._generateGenericParams()
        self.skip_positions = 5000
        self.skip_energies = 5000
        self.integrator = "stochastic"
        self.n_steps = 25000
        self.thermostat = "berendsen"
        self.temp_time_coupling = 1
        self.temperature = 298
        self.temp_groups = "all"
        self.barostat = "no"
        self.random_velocities = True
        self.random_velocities_temperature = 298

    def _generateNPTEquilibrationParams(self):
        self._generateGenericParams()
        self.skip_positions = 5000
        self.skip_energies = 5000
        self.integrator = "stochastic"
        self.n_steps = 25000
        self.thermostat = "berendsen"
        self.temp_time_coupling = 1
        self.temperature = 298
        self.temp_groups = "all"
        self.barostat = "berendsen"
        self.pres_time_coupling = 1
        self.pressure = 1
        self.compressibility = 4.5 * 10**-5
        self.random_velocities = False

    def _generateProductionParams(self):
        self._generateNPTEquilibrationParams()
        self.skip_positions = 10000
        self.skip_energies = 10000
        self.integrator = "stochastic"
        self.n_steps = 2500000
        self.barostat = "parrinello-rahman"

    def _generateVacuumParams(self):
        self.timestep = 0.002
        self.skip_positions = 10000
        self.skip_velocities = 0
        self.skip_forces = 0
        self.skip_energies = 10000
        self.periodic_boundary_conditions = "no"
        self.neighbour_cutoff = 0
        self.neighbour_frequency = 0
        self.coulomb_type = "cutoff"
        self.coulomb_cutoff = 0
        self.vdw_type = "cutoff"
        self.vdw_cutoff = 0
        self.vdw_corr = "no"
        self.constraint = "h_bonds"
        self.constraint_type = "lincs"
        self.integrator = "stochastic"
        self.n_steps = 250000
        self.temp_time_coupling = 1
        self.temperature = 298
        self.temp_groups = "all"
        self.barostat = "no"
        self.random_velocities = True
        self.random_velocities_temperature = 298
        self.__setattr__("cutoff-scheme", "group")
        self.__setattr__("comm_mode", "ANGULAR")

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
            "temp_time_coupling": "tau-t",
            "temperature": "ref-t",

            "barostat": "pcoupl",
            "pres_frequency": "nstpcouple",
            "pres_time_coupling": "tau-p",
            "pressure": "ref-p",
            "compressibility": "compressibility",
            "temp_groups" : "tc-grps",

            "random_velocities": "gen-vel",
            "random_velocities_temperature": "gen-temp",
            "random_velocities_seed": "gen-seed",

            "constraint": "constraints",
            "constraint_type": "constraint-algorithm",

            "free_energy": "free-energy",
            "current_lambda": "init-lambda-state",
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

        value_dict["temp_groups"] = {
            "all" : "system",
        }

        value_dict["barostat"] = {
            "no": "no",
            "berendsen": "berendsen",
            "parrinello-rahman": "Parrinello-Rahman",
            "mttk" : "MTTK"
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
            #bool flags keeping track of special cases
            write_softcore = True

            for name, value in self.__attrs.items():
                #screen for special cases:

                # SPECIAL CASE: constraints - don't write constraint type of there is no constraint
                if name == "constraint_type" and self.constraint == "no":
                    continue

                #SPECIAL CASE: softcore parameters
                if name in ["softcore_vdw", "softcore_coulomb", "softcore_lambda_power", "softcore_alpha"]:
                    if write_softcore:
                        if self.free_energy:
                            if self.softcore_vdw:
                                file.write("{:<30} = {}\n".format("sc-alpha", self.softcore_alpha))
                                file.write("{:<30} = {}\n".format("sc-power", self.softcore_lambda_power))
                            if self.softcore_coulomb:
                                coulomb_str = "yes" if self.softcore_coulomb else "no"
                                file.write("{:<30} = {}\n".format("sc-coul", coulomb_str))
                        write_softcore = False
                    continue

                #translate the parameter name in GROMACS language
                if name in name_dict.keys():
                    name_str = name_dict[name]
                elif name[0] == "_":
                    name_str = name[1:]
                else:
                    name_str = name

                #translate the parameter value in GROMACS language
                if isinstance(value, _BSS.Gateway.Boolean):
                    if value.getValue() is True:
                        value_str = "yes"
                    elif value.getValue() is False:
                        value_str = "no"
                    else:
                        value_str = ""
                elif isinstance(value, _BSS.Gateway.String):
                    if name in value_dict.keys() and value.getValue() in value_dict[name].keys():
                        value_str = value_dict[name][value.getValue()]
                    else:
                        value_str = value.getValue()
                elif "Gateway" in str(type(value)):
                    value_str = str(value.getValue()) if value.getValue() not in [None, ""] else ""
                elif isinstance(value, _Iterable) and not isinstance(value, str):
                    value_str = "\t".join(str(x) for x in value)
                elif value not in [[], None, ""]:
                    value_str = str(value)
                else:
                    continue

                #write a parameter name and value pair
                if name_str and value_str:
                    file.write("{:<30} = {}\n".format(name_str, value_str))
        return filename
