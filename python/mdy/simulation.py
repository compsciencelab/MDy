from simtk.openmm.app import *
import simtk.openmm.app.simulation
from simtk.openmm import *
from simtk.unit import *
from simtk.unit import nanometer
from simtk.unit import angstroms
from simtk.unit import picoseconds
from simtk.unit import picosecond
from simtk.unit import femtoseconds
from sys import stdout
from mdy.bincoor import *
from mdy.xsc import *
from mdy.reporters import *
import mdy.openmmfix.charmmparameterset
import mdy.openmmfix.amberprmtopfile
from natsort import natsorted
from htmd.molecule.molecule import Molecule


class Simulation:
    def __init__(self, config, device=None, platform=None):

        # Now initialize the OpenMM platform
        self.setup_platform(device=device, platform=platform)
        self._stdoutlog = None
        self.restarted = False
        self.setup_rigid_bonds(config)
        self.setup_nonbonded(config)

        if (config.forcefield and not config.structure and not config.parameters and not config.parmfile):
            self.setup_openmm(config)
        elif (not config.forcefield and config.parmfile and not config.structure):
            self.setup_amber(config)
        elif (not config.forcefield and config.parameters and config.structure):
            self.setup_charmm(config)
        else:
            raise NameError("Incompatible parameters found. See 'forcefield' for configuration rules")
        self.setup_integrator(config)
        self.setup_barostat(config)
        self.setup_constraints(config)
        self.setup_flatbottom(config)
        self.setup_simulation(config)

        self.setup_cell_dimensions(config)
        self.setup_coordinates(config)

        self.setup_try_restart(config)
        # NB -- setup_v has to be after the coordinates are set
        # or the energies go nan
        self.setup_velocities(config)

        self.setup_reporters(config)


        # do any requested minimization if this isn't a restart
        if (config.minimize and not self.restarted):
            print("\nMinimizing...\n");
            self.simulation.minimizeEnergy()



        # run the simulation

        if (config.run):
            try:
                print("\nRunning for " + str(config.run) + " steps\n")
                if (self._stdoutlog):
                    self._stdoutlog.headers()
                # self.simulation.step(10)

                # self.simulation.saveCheckpoint( 'foo.chk')
                # self.simulation.loadCheckpoint( 'foo.chk')
                self.simulation.step(config.run)

            except KeyboardInterrupt as e:
                print("\n Terminating at user request during simulation\n")
                return

        # done - write out terminal output files
        if (config.outputname):
            state = self.simulation.context.getState(getPositions=True, getVelocities=True)
            pos = state.getPositions()
            vel = state.getVelocities()
            box = state.getPeriodicBoxVectors(asNumpy=True)
            bincoor_write(pos, config.outputname + ".coor", angstrom)
            bincoor_write(vel, config.outputname + ".vel", angstrom / picosecond)
            xsc_write(box, config.outputname + ".xsc")
        return

    def setup_flatbottom(self, config):
        if not config.flatbottom:
            # Flatbottom constraints are not enabeld
            return

        cfile = config.flatbottom_file
        if not cfile and config.coordinates.endswith(".pdb"):
            cfile = config.coordinates
        if not cfile:
            raise NameError("Flatbottom requires a PDB file - either 'flatbottom-file' or 'coordinates'")

        self._flatbottom_force = None

        m = Molecule(cfile)

        s = config.flatbottom_reference
        if not s:
            s = config.flatbottom_target

        for potentialset in [s]:
            sel = m.atomselect(potentialset, indexes=True)
            coords = m.get("coords", sel=s)

            com = [0, 0, 0] * angstroms
            mass = 0.
            for i in sel:
                xx = m.coords[i, :, 0] * angstroms
                com = com + xx
                mass = mass + 1
            com = com / mass
            r = np.array(config.flatbottom_range) * angstroms
            print(" Applying flatbottom constraints to " + str(len(sel)) + " atoms using: ")
            print("  Target       : '" + s + "'")
            if config.flatbottom_reference:
                print("  Reference    : '" + config.flatbottom_reference + "'")

            d = r.value_in_unit(angstrom)
            d[0] = d[0] + config.flatbottom_offset[0]
            d[1] = d[1] + config.flatbottom_offset[1]
            d[2] = d[2] + config.flatbottom_offset[2]

            target = "target"
            if config.flatbottom_reference:
                target = "reference"
            print("  Bounding box : [%.2f %.2f %.2f] about %s" % (
                (d[0]),
                (d[1]),
                (d[2]),
                target
            )
                  )
            # Create the custom force restraints

            force = CustomExternalForce("""
                 k* ( ( step( periodicdistance(x, y0, z0, x0, y0, z0) - bx0 ) * ( periodicdistance(x, y0, z0, x0, y0, z0) - bx0 )  )^2 +
                    ( step( periodicdistance(x0, y, z0, x0, y0, z0) - by0 ) * ( periodicdistance(x0, y, z0, x0, y0, z0) - by0 )  )^2 +
                    ( step( periodicdistance(x0, y0, z, x0, y0, z0) - bz0 ) * ( periodicdistance(x0, y0, z, x0, y0, z0) - bz0 )  )^2
		 )
		""")

            #       force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)" )
            # NB don't change the order of these - reports.FlatbottomUpdate relies on this ordering
            force.addGlobalParameter("k", config.flatbottom_strength * kilocalories_per_mole / angstroms ** 2)
            force.addGlobalParameter("x0", d[0] * .1)
            force.addGlobalParameter("y0", d[1] * .1)
            force.addGlobalParameter("z0", d[2] * .1)

            force.addGlobalParameter("bx0", config.flatbottom_range[0] * 0.1 / 2.)
            force.addGlobalParameter("by0", config.flatbottom_range[1] * 0.1 / 2.)
            force.addGlobalParameter("bz0", config.flatbottom_range[2] * 0.1 / 2.)

            for i in sel:
                xx = m.coords[i, :, 0] * angstroms
                force.addParticle(int(i))

            self.system.addForce(force)
            self._flatbottom_force = force

        pass

    def setup_constraints(self, config):
        if not config.constraints:
            # Constraints are not enabeld
            return

        cfile = config.constraints_file
        if not cfile and config.coordinates.endswith(".pdb"):
            cfile = config.coordinates
        if not cfile:
            raise NameError("Constraints require a PDB file - either constraints-file or coordinates")

        m = Molecule(cfile)
        sel = m.atomselect(config.constraints_selection, indexes=True)
        coords = m.get("coords", sel=config.constraints_selection)

        print(" Applying constraints to " + str(
            len(sel)) + " atoms using criterion '" + config.constraints_selection + "'")
        # Create the custom force restraints

        #        Temporarily use non-periodic version pending fix for #1192 propagating into release
        force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
        #       force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)" )
        force.addGlobalParameter("k", config.constraints_strength * kilocalories_per_mole / angstroms ** 2)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")

        for i in sel:
            xx = m.coords[i, :, 0] * angstroms
            force.addParticle(int(i), xx.value_in_unit(nanometer))

        self.system.addForce(force)
        self._constraints_force = force

    def setup_try_restart(self, config):
        #        print("\n Warning: Restarting is currently unimplemented. This simulation will start from the beginning.\n" )
        #        return
        if config.restart and os.path.isfile(config.restartfile) and os.path.getsize(config.restartfile) > 0:
            # Try to load a checkpoint file
            try:
                print(" Restarting from checkpoint file '" + config.restartfile + "'")
                # with open( config.restartfile, 'rb') as f:
                #    l = bytearray(f.read()).decode("latin-1")
                #    self.simulation.context.loadCheckpoint( l )
                self.simulation.loadCheckpoint(str(config.restartfile))
                self.restarted = True
            #    print(" Resuming from step " + str(self.simulation.currentStep))
            #    print(self.simulation)
            except:
                # raise
                raise NameError("Restarting from checkpoint file '" + config.restartfile + "' failed")
        else:
            # no checkpoint file exists -- do nothing
            if (not config.coordinates):
                raise NameError("'coordinates' must be set")
            else:
                self.load_coordinates(config.coordinates)
            pass

    def setup_velocities(self, config):
        if not self.restarted:
            if (config.velocities):
                self.load_velocities(config.velocities)
            elif config.temperature is not None:
                self.simulation.context.setVelocitiesToTemperature(config.temperature * kelvin)
            else:
                raise NameError("No initial temperature / velocity field set")
        else:
            # Restarted so vels come from checkpoint file
            pass

    def load_coordinates(self, filename):
        if filename.endswith(".pdb"):
            pos = PDBFile(filename).positions
            pass
        elif filename.endswith(".coor") or filename.endswith(".bin"):
            pos = bincoor_read(filename, angstrom)
            pass
        elif filename.endswith(".inpcrd") or filename.endswith(".crd"):
            pos = AmberInpcrdFile(filename).positions
        else:
            raise NameError("Cannot recognise the file format for '" + filename + "'")
        return pos

    def setup_coordinates(self, config):
        # print(pos)
        self.simulation.context.setPositions(self.load_coordinates(config.coordinates))

    def load_velocities(self, filename):
        if filename.endswith(".vel") or filename.endswith(".bin"):
            vel = bincoor_read(filename, angstrom / picosecond)
            pass
        else:
            raise NameError("Cannot recognise the file format for '" + filename + "'")
        self.simulation.context.setVelocities(vel)

    def setup_barostat(self, config):
        if (config.barostat):
            if (not config.thermostat):
                raise NameError("The barostat requires the thermostat to be enabled")
            if config.barostat_mode == "iso":
                self.system.addForce(
                    MonteCarloBarostat(config.barostat_pressure * bar, config.thermostat_temperature * kelvin))
            elif config.barostat_mode == "aniso":
                self.system.addForce(MonteCarloAnisotropicBarostat(config.barostat_pressure * bar,
                                                                   config.thermostat_temperature * kelvin))
            elif config.barostat_mode == "membrane":
                self.system.addForce( MonteCarloMembraneBarostat(
                 config.barostat_pressure * bar ,
                 config.barostat_surface_tension * bar * nanometer,
                 config.thermostat_temperature * kelvin ,
                 MonteCarloMembraneBarostat.XYIsotropic ,
                 MonteCarloMembraneBarostat.ZFree ) )
            else:
                raise NameError("Unknown barostat mode '" + config.barostat_mode + "'");

    def setup_reporters(self, config):
        if not self.restarted:
            try:
                os.unlink(config.trajfile)
            except:
                pass

        # Set up trajectory writing

        if (config.trajfreq and config.trajfile):
            if (config.trajfile.endswith(".xtc")):

                self.simulation.reporters.append(XTCReporter(config.trajfile, config.trajfreq))
            elif (config.trajfile.endswith(".dcd")):
                if self.restarted:
                    raise NameError("Can't append to DCDs when restarting. Use XTC instead")
                self.simulation.reporters.append(DCDReporter(config.trajfile, config.trajfreq))
            else:
                raise NameError("trajfile format not supported. Filename must end with .xtc or .dcd");

        # Set up logging to stdout
        if (config.energyfreq):
            self._stdoutlog = StdoutLogReporter(
                config.energyfreq, config.run
            )
            self.simulation.reporters.append(self._stdoutlog)

        # Set up checkpointing
        if (config.restartfile):
            freq = 25000
            if (config.trajfile):
                freq = config.trajfreq
            self.simulation.reporters.append(CheckpointReporter(config.restartfile, freq));

        if config.constraints and config.constraints_decay > 0.:
            self.simulation.reporters.append(
                ConstraintsDecay(config.constraints_strength, config.constraints_decay, self._constraints_force,
                                 config.timestep))

        if config.flatbottom and (config.flatbottom_decay or config.flatbottom_reference):
            self.simulation.reporters.append(FlatbottomUpdate( self._flatbottom_force, config))

    def setup_simulation(self, config):
        self.simulation = simtk.openmm.app.Simulation(
            self.topology,
            self.system,
            self.integrator,
            self.platform,
            self.properties
        )

    def setup_integrator(self, config):
        if (config.thermostat):
            self.integrator = LangevinIntegrator(
                config.thermostat_temperature * kelvin,
                config.thermostat_damping / picoseconds,
                config.timestep * femtoseconds
            )
            self.integrator.setConstraintTolerance(0.00001)

        else:
            self.integrator = VerletIntegrator(
                config.timestep * femtoseconds
            )

        self.integrator.setConstraintTolerance(1e-5)

    def setup_nonbonded(self, config):
        if (config.pme):
            self.nbmethod = PME
        else:
            self.nbmethod = CutoffPeriodic

        self.cutoff = config.cutoff * angstroms
        self.switchdist = config.switchdist * angstroms
        if (self.switchdist >= self.cutoff):
            # Turn it off
            self.switchdist = 0.

    def setup_openmm(self, config):
        ff = ForceField(*config.forcefield)
        if not config.coordinates.endswith(".pdb"):
            raise NameError("'coordinates' must be a PDB when starting with 'forcefield'")
        self.topology = PDBFile(config.coordinates).topology
        self.system = ff.createSystem(
            self.topology,
            nonbondedMethod=self.nbmethod,
            nonbondedCutoff=self.cutoff,
            constraints=self.constraints,
            rigidWater=self.rigidwater,
            removeCMMotion=True,
            hydrogenMass=self.hmass,
            switchDistance=self.switchdist
        )

    def setup_amber(self, config):
        self.parameters = mdy.openmmfix.amberprmtopfile.AmberPrmtopFile(config.parmfile)

        self.system = self.parameters.createSystem(
            nonbondedMethod=self.nbmethod,
            nonbondedCutoff=self.cutoff,
            constraints=self.constraints,
            hydrogenMass=self.hmass,
            rigidWater=self.rigidwater,
            removeCMMotion=True,
            switchDistance=self.switchdist
        )
        self.topology = self.parameters.topology

    def setup_charmm(self, config):

        psf = CharmmPsfFile(config.structure)

        dim = self.get_cell_dimensions(config)
        psf.setBox(
            dim[0][0].value_in_unit(nanometer),
            dim[1][1].value_in_unit(nanometer),
            dim[2][2].value_in_unit(nanometer)
        )

        self.parameters = mdy.openmmfix.charmmparameterset.CharmmParameterSet(*config.parameters, permissive=True)
        #        for pfile in  config.parameters:
        #           self.parameters.readParameterFile(pfile ) #, permissive=True)
        #           self.parameters.readStreamFile(pfile ) #, permissive=True)

        self.system = psf.createSystem(
            self.parameters,
            nonbondedMethod=self.nbmethod,
            nonbondedCutoff=self.cutoff,
            constraints=self.constraints,
            hydrogenMass=self.hmass,
            rigidWater=self.rigidwater,
            removeCMMotion=True,
            switchDistance=self.switchdist
        )
        self.topology = psf.topology

    def setup_rigid_bonds(self, config):

        # First configure the constraints based on the timestep
        # The logic is:
        # dt<=1 - no constraints
        # dt<=2.5 - hbonds, rigidwater
        # dt>2.5  - hbonds, rigidwater, hmass scaling

        self.timestep = config.timestep
        self.rigidwater = False
        self.hmass = None
        self.constraints = None
        if (self.timestep > 1):
            self.constraints = HBonds
            self.rigidwater = True
        if (self.timestep > 2.5):
            self.hmass = 4 * amu

    def setup_cell_dimensions(self, config):
        dim = self.get_cell_dimensions(config)
        self.simulation.context.setPeriodicBoxVectors(dim[0], dim[1], dim[2])

    def get_cell_dimensions(self, config):
        if (config.extendedsystem):
            return xsc_read(config.extendedsystem)
            # read xsc
            pass
        elif (config.celldimension):
            return [
                       Vec3(float(config.celldimension[0]), 0., 0.),
                       Vec3(0., float(config.celldimension[1]), 0.),
                       Vec3(0., 0., float(config.celldimension[2]))
                   ] * angstroms
            pass
        elif config.coordinates and config.coordinates.endswith(".pdb"):
            pdb = PDBFile(config.coordinates)
            vec = pdb.topology.getPeriodicBoxVectors()
            if not vec or len(vec) != 3:
                raise NameError(
                    "PDB has invalid unit cell dimensions . Either add CRYST line to PDB or set 'celldimension' or 'extendedsystem'")
            return vec
        elif config.coordinates and (config.coordinates.endswith(".inpcrd") or config.coordinates.endswith(".crd")):
            coor = AmberInpcrdFile(config.coordinates)
            vec = coor.getBoxVectors()
            if len(vec) != 3:
                raise NameError("Inpcrd has invalid unit cell dimensions")
            return vec

        else:
            raise NameError(
                "Can't determine periodic cell size - either 'celldimension' or 'extendedsystem' must be set")

    def setup_platform(self, device=None, platform=None):
        self.platform = None

        version = Platform.getOpenMMVersion()
        plugindir = Platform.getDefaultPluginsDirectory()
        print(" OpenMM details:")
        print("  Version     : OpenMM " + str(version))
        print("  Plugin dir  : " + str(plugindir))

        # Try loading the plugins and checking for errors
        # Platform.loadPluginsFromDirectory( plugindir )
        errs = Platform.getPluginLoadFailures()
        if (len(errs)):
            print("\n  Some errors were found loading plugins. Some platforms may not be available: \n");
            for x in errs:
                print(x)
            print("")
        #           sys.exit(1)

        num = Platform.getNumPlatforms();
        speed = []
        name_by_speed = dict()
        for i in range(num):
            pp = Platform.getPlatform(i)
            s = pp.getSpeed()
            speed.append(s)
            name_by_speed[s] = pp.getName()

        speed = natsorted(speed)
        speed.reverse()
        print("  Platforms   : ", end="")
        for i in speed:
            print(name_by_speed[i] + " ", end="")
        #          for k in ( pp.getPropertyNames() ):
        #             print( "    " + k  )
        print("")

        self.properties = {}

        print("")
        if platform:
            try:
                self.platform = Platform.getPlatformByName(platform)
                #              print(self.platform.getPropertyNames())
                if "CudaPrecision" in self.platform.getPropertyNames():
                    self.properties["CudaPrecision"] = "mixed"
                if "OpenCLPrecision" in self.platform.getPropertyNames():
                    self.properties["OpenCLPrecision"] = "mixed"
                if device and "CudaDeviceIndex" in self.platform.getPropertyNames():
                    self.properties["CudaDeviceIndex"] = str(device)
                if device and "OpenCLDeviceIndex" in self.platform.getPropertyNames():
                    self.properties["OpenCLDeviceIndex"] = str(device)
                if "CpuThreads" in self.platform.getPropertyNames():
                    self.properties["CpuThreads"] = str(self.num_cpus())

                print("  Using platform " + platform)
                return
            except:
                raise NameError("Could not initialize OpenMM platform '" + platform + "'")

        for i in speed:
            name = name_by_speed[i];
            if not self.platform:
                try:
                    self.platform = Platform.getPlatformByName(name)
                    print("  Using platform " + name)
                except:
                    print("Could not initialise platform")
                    pass
        print("")
        if not self.platform:
            raise NameError("Could not initialize OpenMM platform")

    def num_cpus(self):
        # See if ncpus envvar is set (eg a cluster job)
        ncpus = os.environ['NCPUS']
        if ncpus:
            return ncpus

        import multiprocessing
        return multiprocessing.cpu_count()
