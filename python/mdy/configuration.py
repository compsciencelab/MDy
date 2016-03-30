from mdy.command import Command
from difflib import get_close_matches
import re
import os


class Configuration:
    def __init__(self, config=None):
        c = Command.get_default_configuration()

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

    def __setattr__(self, key, value):
        if key.startswith("_"):
            self.__dict__[key] = value
            return

        key = key.replace('_', '-')
        key = key.lower()
        key = Command.test_deprecation(key, value)
        if key:
            value = Command.validate(key, value, basedir=self._basedir)
            self.__dict__[key] = value;

    def __getattr__(self, key):
        if key.startswith("_"):
            self.__dict__[key] = value
            return
        try:
            key = key.replace('_', '-');
            key = key.lower()
            return self.__dict__[key]
        except:
            errstr = "Command '" + key + "' not found.";
            match = get_close_matches(key, Command.commands)
            if match:
                errstr = errstr + " Try '" + match[0] + "'";
            raise NameError(errstr)

    def __repr__(self):
        return (self.__str__())

    def __str__(self):
        s = ""
        for cmd in sorted(self.__dict__):
            if not cmd.startswith("_"):
                s = s + "%25s %s\n" % (cmd, str(self.__dict__[cmd]))
        return s

    def set(self, key, value):
        self.__setattr__(key, value)

    @staticmethod
    def help(command):
        return Command.help(command)

    def save(self, filename):
        """Write out configuration to a file"""
        with open(filename, "w") as fh:
            for cmd in sorted(self.__dict__):
                if (self.__dict__[cmd] is not None):
                    print("%25s %s" % (cmd, str(self.__dict__[cmd])), file=fh)

    def stage(self, directory):
        directory = os.path.abspath(directory)
        os.mkdirs(directory)
        self.save(os.path.join(directory, "input"))


if __name__ == "__main__":
    c = Configuration()
    c.set('barostat', 'on')
    c.barostat = True
    print(c)
