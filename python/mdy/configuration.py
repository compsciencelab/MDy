from htmd.protocols.protocolinterface import ProtocolInterface, TYPE_FLOAT, TYPE_INT, RANGE_0POS, RANGE_POS
import re
import os


class Configuration(ProtocolInterface):
    def __init__(self, config=None):
        super().__init__()
        self._cmdDeprecated('amber')
        self._cmdDeprecated('bincoordinates', 'coordinates')
        self._cmdDeprecated('binvelocities', 'velocities')
        self._cmdDeprecated('hydrogenscale')
        self._cmdDeprecated('rigidbonds')
        self._cmdDeprecated('1-4scaling')
        self._cmdDeprecated('scaling1-4')
        self._cmdDeprecated('switching')
        self._cmdDeprecated('exclude')
        self._cmdDeprecated('pmefreq')
        self._cmdDeprecated('fullelectfrequency')
        self._cmdDeprecated('dcdfreq', 'trajfreq')
        self._cmdDeprecated('dcdfile', 'trajfile')
        self._cmdDeprecated('xtcfreq', 'trajfreq')
        self._cmdDeprecated('xtcfile', 'trajfile')
        self._cmdDeprecated('langevin', 'thermostat')
        self._cmdDeprecated('langevintemp', 'thermostat-temperature')
        self._cmdDeprecated('langevindamping', 'thermostat-damping')
        self._cmdDeprecated('berendsenpressure', 'barostat')
        self._cmdDeprecated('berendsenpressuretarget', 'barostat-pressure')
        self._cmdDeprecated('pmegridsizex')
        self._cmdDeprecated('pmegridsizey')
        self._cmdDeprecated('pmegridsizez')
        self._cmdDeprecated('pmegridspacing')
        self._cmdDeprecated('consref', 'constraints-file')
        self._cmdDeprecated('constraintscaling', 'constraints-strength')
        self._cmdDeprecated('restartfreq')
        self._cmdDeprecated('celldimensions', 'celldimension')

        self._cmdFile('coordinates', 'str', os.path.join('help', 'en', 'coordinates'), None, exist=True)
        self._cmdFile('velocities', 'str', os.path.join('help', 'en', 'velocities'), None, exist=True)
        self._cmdValue('temperature', 'float', os.path.join('help', 'en', 'temperature'), None, TYPE_FLOAT, RANGE_0POS)
        self._cmdFile('parameters', 'str', os.path.join('help', 'en', 'parameters'), None, exist=True, multiple=True)
        self._cmdFile('outputname', 'str', os.path.join('help', 'en', 'outputname'), "output", writable=True)
        self._cmdFile('structure', 'str', os.path.join('help', 'en', 'structure'), None, exist=True)
        self._cmdFile('parmfile', 'str', os.path.join('help', 'en', 'parmfile'), None, exist=True)
        self._cmdValue('timestep', 'float', os.path.join('help', 'en', 'timestep'), 4.0, TYPE_FLOAT, RANGE_0POS, units='fs')
        self._cmdBinary('pme', 'bool', os.path.join('help', 'en', 'binary'), True)
        self._cmdValue('cutoff', 'float', os.path.join('help', 'en', 'cutoff'), 9.0, TYPE_FLOAT, RANGE_POS, units='angstrom')
        self._cmdValue('switchdist', 'float', os.path.join('help', 'en', 'switchdist'), 7.5, TYPE_FLOAT, RANGE_0POS, units='angstrom')
        self._cmdValue('trajfreq', 'int', os.path.join('help', 'en', 'trajfreq'), 25000, TYPE_INT, RANGE_0POS, units='steps')
        self._cmdFile('trajfile', 'str', os.path.join('help', 'en', 'trajfile'), "output.xtc", writable=True)
        self._cmdBinary('minimize', 'bool', os.path.join('help', 'en', 'minimize'), False)
        self._cmdTimestep('run', 'float', os.path.join('help', 'en', 'run'), 1000)
        self._cmdBinary('thermostat', 'bool', os.path.join('help', 'en', 'thermostat'), True)
        self._cmdValue('thermostat-damping', 'float', os.path.join('help', 'en', 'thermostat-damping'), 0.1, TYPE_FLOAT, RANGE_POS, units='/ps')
        self._cmdValue('hermostat-temperature', 'float', os.path.join('help', 'en', 'thermostat-temperature'), 300., TYPE_FLOAT, RANGE_0POS, units='K')
        self._cmdBinary('barostat', 'bool', os.path.join('help', 'en', 'barostat'), False)
        self._cmdValue('barostat-pressure', 'float', os.path.join('help', 'en', 'barostat-pressure'), 1.01325, TYPE_FLOAT, RANGE_0POS, units='bar')
        # Command.commands[ 'barostat-damping' ]  = Command.Value( 0.1, Command.TYPE_FLOAT, Command.RANGE_POS, 1 )
        self._cmdList('barostat-mode', 'str', os.path.join('help', 'en', 'barostat-mode'), 'iso', ["iso", "aniso", "membrane"])
        self._cmdValue('barostat-surface-tension', 'float', os.path.join('help', 'en', 'barostat-surface-tension'), 0, TYPE_FLOAT, RANGE_0POS, units='bar.nanometre')
        self._cmdBinary('restart', 'bool', os.path.join('help', 'en', 'restart'), True)
        self._cmdFile('restartfile', 'str', os.path.join('help', 'en', 'restartfile'), 'restart.chk', writable=True)
        # Command.commands[ 'restartfreq' ]    = Command.Value( 25000, Command.TYPE_INT , Command.RANGE_0POS, 1 )
        self._cmdValue('energyfreq', 'float', os.path.join('help', 'en', 'energyfreq'), 5000, TYPE_INT, RANGE_POS, units='steps')
        self._cmdValue('celldimension', 'float', os.path.join('help', 'en', 'celldimension'), 5000, TYPE_FLOAT, RANGE_POS, 3, units='angstrom')
        self._cmdFile('extendedsystem', 'str', os.path.join('help', 'en', 'extendedsystem'), None, exist=True)
        self._cmdBinary('constraints', 'bool', os.path.join('help', 'en', 'constraints'), False)
        self._cmdFile('constraints-file', 'str', os.path.join('help', 'en', 'constraints-file'), None, exist=True)
        self._cmdValue('constraints-strength', 'float', os.path.join('help', 'en', 'constraints-strength'), 1.0, TYPE_FLOAT, RANGE_0POS, 1, units='kcal/mol/angstrom**2')
        self._cmdString('constraints-selection', 'str', os.path.join('help', 'en', 'constraints-selection'), "beta > 0")
        self._cmdValue('constraints-decay', 'float', os.path.join('help', 'en', 'constraints-decay'), 0, TYPE_FLOAT, RANGE_0POS, 1, units='kcal/mol/ps"')
        self._cmdBinary('flatbottom', 'bool', os.path.join('help', 'en', 'flatbottom'), False)
        self._cmdFile('flatbottom-file', 'str', os.path.join('help', 'en', 'flatbottom-file'), None, exist=True)
        self._cmdValue('flatbottom-strength', 'float', os.path.join('help', 'en', 'flatbottom-strength'), 1.0, TYPE_FLOAT, RANGE_0POS, 1, units='kcal/mol/angstrom**2')
        self._cmdString('flatbottom-target', 'str', os.path.join('help', 'en', 'flatbottom-target'), "resname MOL and noh")
        self._cmdString('flatbottom-reference', 'str', os.path.join('help', 'en', 'flatbottom-reference'), None)
        self._cmdValue('flatbottom-range', 'float', os.path.join('help', 'en', 'flatbottom-range'), [20, 20, 20], TYPE_FLOAT, RANGE_0POS, 3, units='angstrom')
        self._cmdValue('flatbottom-offset', 'float', os.path.join('help', 'en', 'flatbottom-offset'), [0, 0, 0], TYPE_FLOAT, RANGE_0POS, 3, units='angstrom')
        self._cmdValue('flatbottom-decay', 'float', os.path.join('help', 'en', 'flatbottom-decay'), 0, TYPE_FLOAT, RANGE_0POS, 1, units='kcal/mol/ps')
        self._cmdFile('forcefield', 'str', os.path.join('help', 'en', 'forcefield'), None, exist=True, multiple=True) # , check=check_ff)

        inputfile = None

        if (config):
            if os.path.isfile(config):
                self._basedir = os.path.dirname((config))
                inputfile = config
            elif os.path.isdir(config):
                self._basedir = config  # (os.path.abspath( config ))
                inputfile = os.path.join(self._basedir, "input")
                #              print(inputfile)
                if not os.path.isfile(inputfile):
                    raise NameError("Input file not found in config directory")
            else:
                raise NameError("Invalid input. 'config' must be an  input file or a directory containing 'input'")

        for i in c.keys():
            self.__dict__[i] = c[i]

        if not inputfile:
            return

        # Read in key-value pairs, excluding any comments beginning #
        f = open(inputfile, "r")
        ll = 1
        for linex in f:
            line = linex
            line = re.sub(r'\t', ' ', line)
            line = re.sub(r'\n', '', line)
            line = re.sub(r'\r', '', line)
            line = re.sub(r'#.*$', '', line, flags=re.M | re.S)
            # Remove any leading or trailing whitespace
            line = line.strip()
            # Remove comment
            if (len(line) > 0):
                key = re.sub('\s.*$', '', line, flags=re.M | re.S)
                val = re.sub('^[^\s]+\s', '', line, flags=re.M | re.S)
                key = key.strip()
                val = val.strip()
                val = re.split(r'\s+', val)
                try:
                    self.__setattr__(key, val)
                except NameError as e:
                    raise NameError(" Line " + str(ll) + "\t: " + linex + "        \t: " + str(e))
            ll = ll + 1
        f.close()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        s = ""
        for cmd in sorted(self.__dict__):
            if not cmd.startswith("_"):
                s += "%25s %s\n" % (cmd, str(self.__dict__[cmd]))
        return s

    def set(self, key, value):
        self.__setattr__(key, value)

    def save(self, filename):
        """Write out configuration to a file"""
        with open(filename, "w") as fh:
            for cmd in sorted(self.__dict__):
                if self.__dict__[cmd] is not None:
                    print("%25s %s" % (cmd, str(self.__dict__[cmd])), file=fh)

    def stage(self, directory):
        directory = os.path.abspath(directory)
        os.makedirs(directory)
        self.save(os.path.join(directory, "input"))


if __name__ == "__main__":
    c = Configuration()
    c.set('barostat', 'on')
    c.barostat = True
    print(c)
