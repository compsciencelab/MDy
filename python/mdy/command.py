from difflib import get_close_matches
import os
import sys
import re
import shutil
import inspect


class Command:
    setup = False
    commands = {}

    RANGE_ANY = 0
    RANGE_POS = 1
    RANGE_NEG = 2
    RANGE_0POS = 3
    RANGE_0NEG = 4

    TYPE_INT = 0
    TYPE_FLOAT = 1

    @staticmethod
    def test_deprecation(key, value):
        if (key in Command.deprecated):
            newkey = Command.deprecated[key]
            if (not newkey):
                print(" Command '" + key + "'\t is deprecated and is no longer required")
            else:
                print(" Command '" + key + "'\t is deprecated and replaced by '" + newkey + "'")
                return newkey
        else:
            return key

    @staticmethod
    def add_commands():
        if Command.setup:
            return
        Command.deprecated = {
            "amber": None,
            "bincoordinates": "coordinates",
            "binvelocities": "velocities",
            "hydrogenscale": None,
            "rigidbonds": None,
            "1-4scaling": None,
            "scaling1-4": None,
            "switching": None,
            "exclude": None,
            "pmefreq": None,
            "fullelectfrequency": None,
            "dcdfreq": "trajfreq",
            "dcdfile": "trajfile",
            "xtcfreq": "trajfreq",
            "xtcfile": "trajfile",
            "langevin": "thermostat",
            "langevin": "thermostat",
            "langevintemp": "thermostat-temperature",
            "langevindamping": "thermostat-damping",
            "berendsenpressure": "barostat",
            "berendsenpressuretarget": "barostat-pressure",
            "pmegridsizex": None,
            "pmegridsizey": None,
            "pmegridsizez": None,
            "pmegridspacing": None,
            "consref": "constraints-file",
            "constraintscaling": "constraints-strength",
            "restartfreq": None,
            "celldimensions": "celldimension"
        }

        Command.commands['coordinates'] = Command.File(None, exist=True, )
        Command.commands['velocities'] = Command.File(None, exist=True)
        Command.commands['temperature'] = Command.Value(None, Command.TYPE_FLOAT, Command.RANGE_0POS, 1)
        Command.commands['parameters'] = Command.File(None, exist=True, multiple=True)
        Command.commands['outputname'] = Command.File("output", writable=True)
        Command.commands['structure'] = Command.File(None, exist=True)
        Command.commands['parmfile'] = Command.File(None, exist=True)
        Command.commands['timestep'] = Command.Value(4.0, Command.TYPE_FLOAT, Command.RANGE_0POS, 1, units="fs")
        Command.commands['pme'] = Command.Binary(True)
        Command.commands['cutoff'] = Command.Value(9.0, Command.TYPE_FLOAT, Command.RANGE_POS, 1, units="angstrom")
        Command.commands['switchdist'] = Command.Value(7.5, Command.TYPE_FLOAT, Command.RANGE_0POS, 1, units="angstrom")
        Command.commands['trajfreq'] = Command.Value(25000, Command.TYPE_INT, Command.RANGE_0POS, 1, units="steps")
        Command.commands['trajfile'] = Command.File("output.xtc", writable=True)
        Command.commands['minimize'] = Command.Binary(False)
        Command.commands['run'] = Command.Timestep(1000)
        Command.commands['thermostat'] = Command.Binary(True)
        Command.commands['thermostat-damping'] = Command.Value(0.1, Command.TYPE_FLOAT, Command.RANGE_POS, 1,
                                                               units="/ps")
        Command.commands['thermostat-temperature'] = Command.Value(300., Command.TYPE_FLOAT, Command.RANGE_0POS, 1,
                                                                   units="K")

        Command.commands['barostat'] = Command.Binary(False)
        #		Command.commands[ 'barostat-damping' ]  = Command.Value( 0.1, Command.TYPE_FLOAT, Command.RANGE_POS, 1 )
        Command.commands['barostat-pressure'] = Command.Value(1.01325, Command.TYPE_FLOAT, Command.RANGE_0POS, 1,
                                                              units="bar")
        Command.commands['barostat-mode'] = Command.List("iso", ["iso", "aniso", "membrane"])
        Command.commands['barostat-surface-tension'] = Command.Value(0, Command.TYPE_FLOAT, Command.RANGE_0POS, 1,
                                                                     units="bar.nanometre")

        Command.commands['restart'] = Command.Binary(True)
        Command.commands['restartfile'] = Command.File("restart.chk", writable=True)
        #		Command.commands[ 'restartfreq' ]    = Command.Value( 25000, Command.TYPE_INT , Command.RANGE_0POS, 1 )
        Command.commands['energyfreq'] = Command.Value(5000, Command.TYPE_INT, Command.RANGE_POS, 1, units="steps")
        Command.commands['celldimension'] = Command.Value(None, Command.TYPE_FLOAT, Command.RANGE_POS, 3,
                                                          units="angstrom")
        Command.commands['extendedsystem'] = Command.File(None, exist=True)
        Command.commands['constraints'] = Command.Binary(False)
        Command.commands['constraints-file'] = Command.File(None, exist=True)
        Command.commands['constraints-strength'] = Command.Value(1.0, Command.TYPE_FLOAT, Command.RANGE_0POS, 1,
                                                                 units="kcal/mol/angstrom**2")
        Command.commands['constraints-selection'] = Command.String("beta > 0", "Atom selection")
        Command.commands['constraints-decay'] = Command.Value(0, Command.TYPE_FLOAT, Command.RANGE_0POS, 1,
                                                              units="kcal/mol/ps")
        Command.commands['flatbottom'] = Command.Binary(False)
        Command.commands['flatbottom-file'] = Command.File(None, exist=True)
        Command.commands['flatbottom-strength'] = Command.Value(1.0, Command.TYPE_FLOAT, Command.RANGE_0POS, 1,
                                                                units="kcal/mol/angstrom**2")
        Command.commands['flatbottom-target'] = Command.String("resname MOL and noh", "Atom selection")
        Command.commands['flatbottom-reference'] = Command.String(None, "Atom selection")
        Command.commands['flatbottom-range'] = Command.Value([20, 20, 20], Command.TYPE_FLOAT, Command.RANGE_0POS, 3,
                                                             units="angstrom")
        Command.commands['flatbottom-offset'] = Command.Value([0, 0, 0], Command.TYPE_FLOAT, Command.RANGE_ANY, 3,
                                                              units="angstrom")
        Command.commands['flatbottom-decay'] = Command.Value(0, Command.TYPE_FLOAT, Command.RANGE_0POS, 1,
                                                             units="kcal/mol/ps")
        Command.commands['forcefield'] = Command.File(None, exist=False, multiple=True, check=check_ff)
        Command.setup = True

    @staticmethod
    def get_help_string(cmd):
        libdir = (os.path.dirname(inspect.getfile(Command)))

        libdir = (os.path.join(libdir, "help"))

        libdir = (os.path.join(libdir, "en"))

        helpfile = (os.path.join(libdir, os.path.basename(cmd)))

        helpfile = helpfile + ".txt"

        if not os.path.isfile(helpfile):
            return "No help found"
        else:
            with open(helpfile) as ff:
                return ff.read().strip()

    @staticmethod
    def pretty_print(cmd):
        print(cmd)
        return

        width = shutil.get_terminal_size(fallback=(80, 25))
        width = width.columns
        if width < 40:
            width = 80
        tok = re.split(r'\s+', cmd);
        l = 4
        print("   ", end="")
        for t in tok:
            if (l + len(t) + 1 <= width):
                print(" " + t, end="");
                l = l + len(t) + 1
            else:
                l = 4;
                print("\n    " + t, end="");
                l = l + len(t) + 1

    @staticmethod
    def help(cmd):
        Command.add_commands()
        if (not cmd):
            print("\n  Valid configuration file commands:\n");
            for a in sorted(Command.commands.keys()):
                print("    " + a)
            print("\n  acemd --command [command] for detailed help on a specific command\n\n");
            # print sections
            pass
        else:
            match = get_close_matches(cmd, Command.commands)
            if not match:
                print("\nNo matching command for '" + cmd + "' found\n")
                return
            else:
                helpstr = Command.get_help_string(match[0])
                print("\n   " + match[0] + " " + Command.commands[match[0]].args() + "\n")

                Command.pretty_print(helpstr);

                units = Command.commands[match[0]].units
                if not units:
                    units = ""
                print("\n   Default: " + str(Command.commands[match[0]].default) + " " + units + "\n");



                # find

    @staticmethod
    def get_default_configuration():
        Command.add_commands()
        ret = {}
        for i in Command.commands.keys():
            ret[i] = Command.commands[i].default
        return ret

    @staticmethod
    def validate(key, value, basedir=None):
        try:
            cmd = Command.commands[key]
        except:
            strerror = "Command '" + key + "' not found.";
            match = get_close_matches(key, Command.commands)
            if match:
                strerror = strerror + " Try '" + match[0] + "'";
            raise NameError(strerror)

        return cmd.validate(value, basedir=basedir)

    class String:
        def __init__(self, default, args=None):
            self.value = default
            self.default = default
            self.args = args
            self.units = None

        def args(self):
            if self.args:
                return args
            else:
                return "free string"

        def validate(self, value_list, basedir=None):
            import re
            if type(value_list) is not list and type(value_list) is not tuple:
                value_list = [value_list]
            value = ""
            for v in value_list:
                value = value + " " + v
            self.value = value.strip()
            self.value = re.sub('^"', '', self.value)
            self.value = re.sub('"$', '', self.value)
            self.value = re.sub("^'", '', self.value)
            self.value = re.sub("'$", '', self.value)
            return self.value

    class Timestep:
        def __init__(self, default):
            self.value = default
            self.default = default
            self.units = None

        def args(self):
            return "[ number (ps|ns|us) ]"

        def validate(self, value_list, basedir=None):
            if type(value_list) is not list and type(value_list) is not tuple:
                value_list = [value_list]
            value = value_list[0]
            try:
                if len(value_list) == 1:
                    self.value = int(value_list[0])
                elif len(value_list) == 2:
                    self.value = float(value_list[0])
                    if value_list[1].lower() == "ps":
                        self.value = self.value * 1000 / 4  # TODO FIXME -- note the hard-coded assumption about the timestep
                    if value_list[1].lower() == "ns":
                        self.value = self.value * 1000000 / 4
                    if value_list[1].lower() == "us":
                        self.value = self.value * 1000000000 / 4
                    self.value = int(self.value)
                else:
                    raise Error()
            except:
                raise NameError("Value must be >= 0 and may have a unit suffix [ps,ns,us]")
            return self.value

    class List:
        def __init__(self, default, sel):
            self.value = default
            self.default = default
            self.sel = sel
            self.units = None

        def args(self):
            return str(self.sel)

        def validate(self, value_list, basedir=None):
            if type(value_list) is not list and type(value_list) is not tuple:
                value_list = [value_list]
            value = value_list[0]
            try:
                if value in self.sel:
                    self.value = value
                return value
            except:
                raise NameError("Value must be one of " + str(self.sel))

    class Binary:
        def __init__(self, default):
            self.value = default
            self.default = default
            self.units = None

        def args(self):
            return "[ on | off ]"

        def validate(self, value_list, basedir=None):

            if type(value_list) is not list and type(value_list) is not tuple:
                value_list = [value_list]
            value = value_list[0]

            try:
                if (value == True):
                    value = True
                    return value
                if (value == False or value == None):
                    value = False
                    return value

                value = value.lower()

                if (value == "yes" or value == "on" or value == "1" or value == "all" or value == "true"):
                    value = True
                elif (value == "no" or value == "off" or value == "0" or value == "none" or value == "false"):
                    value = False
                else:
                    raise Error("")
            except:
                raise NameError("Value must be binary ( on|off )")
            return value

    class File:

        def __init__(self, default, exist=False, writable=False, multiple=False, check=None):
            self.multiple = multiple
            self.value = default
            self.must_exist = exist
            self.writable = writable
            self.default = default
            self.check = check

        def args(self):
            dd = "file"
            if self.must_exist:
                dd = "input file"
            if self.writable:
                dd = "output file"
            if self.multiple:
                dd = "list of " + dd + "s"
            return "[ " + dd + " ]"

        def validate(self, value_list, basedir=None):
            if type(value_list) is not list and type(value_list) is not tuple:
                value_list = [value_list]

            l = []
            for value in value_list:
                if (self.must_exist):
                    found = False
                    if basedir and os.path.isfile(os.path.join(basedir, value)):
                        value = os.path.join(basedir, value)
                        found = True
                    if not found and os.path.isfile(value):
                        found = True
                    if not found:
                        raise NameError("File '" + value + "' does not exist")
                    if not os.access(value, os.R_OK):
                        raise NameError("File '" + value + "' cannot be read")
                if (self.writable):
                    if basedir:
                        value = os.path.join(basedir, value)
                    try:
                        f = open(value, "a+");
                        f.close();
                    except:
                        raise NameError("File '" + value + "' is not writable")
                if (self.check):
                    value = self.check(value, basedir=basedir)

                l.append(value)

            # Test to see if the file exists
            if len(l) != 1 and not self.multiple:
                raise NameError("Value must be a single filename")
            if not self.multiple:
                self.value = l[0]
            else:
                self.value = l
            return self.value

    class Value:
        def __init__(self, default, datatype, valid_range, list_len, units=None):
            self.default = default
            self.valid_range = valid_range
            self.list_len = list_len
            self.data_type = datatype
            self.units = units

        def validate(self, value_list, basedir=None):
            valid_range = self.valid_range

            if type(value_list) is not list and type(value_list) is not tuple:
                value_list = [value_list]

            l = []
            if self.list_len != len(value_list):
                if (self.list_len == 0):
                    raise NameError("Value must be a scalar")
                else:
                    raise NameError("Value must be a vector with " + str(self.list_len) + " elements")

            for value in value_list:
                if (self.data_type == Command.TYPE_FLOAT):
                    strfudge = "a number"
                elif (self.data_type == Command.TYPE_INT):
                    strfudge = "an integer"

                try:
                    try:
                        if (self.data_type == Command.TYPE_FLOAT):
                            value = float(value)
                        if (self.data_type == Command.TYPE_INT):
                            value = int(value)
                    except:
                        raise NameError("Value is not " + strfudge);

                    if (valid_range == Command.RANGE_POS and value <= 0):
                        raise NameError("Value must be " + strfudge + " > 0")
                    if (valid_range == Command.RANGE_0POS and value < 0):
                        raise NameError("Value must be " + strfudge + " >=0")
                    if (valid_range == Command.RANGE_NEG and value >= 0):
                        raise NameError("Value must be " + strfudge + " < 0")
                    if (valid_range == Command.RANGE_0NEG and value > 0):
                        raise NameError("Value must be " + strfudge + " <=0")

                except NameError as e:
                    raise NameError(e)
                if (self.list_len == 1):
                    l = value
                else:
                    l.append(value)

            self.value = l
            return self.value

        def args(self):
            if (self.data_type == Command.TYPE_FLOAT):
                strfudge = "number "
            elif (self.data_type == Command.TYPE_INT):
                strfudge = "integer "
            if (self.valid_range == Command.RANGE_POS):
                strfudge = strfudge + "> 0"
            if (self.valid_range == Command.RANGE_0POS):
                strfudge = strfudge + ">= 0"
            if (self.valid_range == Command.RANGE_NEG):
                strfudge = strfudge + "< 0"
            if (self.valid_range == Command.RANGE_0NEG):
                strfudge = strfudge + "<= 0"
            if self.list_len > 1:
                strfudge = str(self.list_len) + " x " + strfudge;

            xx = "[ " + strfudge + " ]"

            return xx


def check_ff(value, basedir=None):
    import inspect
    import simtk.openmm.app.forcefield
    import os
    import natsort
    import re

    if (basedir):
        vv = os.path.join(basedir, value)
        if (os.path.isfile(vv)):
            return vv
        if (os.path.isfile(vv + ".xml")):
            return vv + ".xml"

    if (os.path.isfile(value)):
        return value
    if (os.path.isfile(value + ".xml")):
        return value + ".xml"

    pp = (os.path.dirname(inspect.getfile(simtk.openmm.app.forcefield.ForceField)))
    pp = os.path.join(pp, "data")
    i = 0
    files = os.listdir(pp)
    if (value in files):
        return value
    if (value + ".xml" in files):
        return value + ".xml"

    raise NameError("Forcefield parameter file is not found")
