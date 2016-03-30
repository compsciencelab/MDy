#!/usr/bin/env python

import os
import argparse
import re
import sys
from mdy.configuration import Configuration
from mdy.command import Command
from mdy.simulation import Simulation
from simtk.openmm import Platform
from htmdx.cli import *


def syntax():
    print("")
    print(" HTMD Molecular Dynamics 2016 (c) Acellera\n")
    print("  Syntax: cli [args]:\n")
    print("    --input [inputfile]      : Specify the input file. Default value 'input'")
    print("    --device [device-number] : Set the GPU to run on ")
    print("    --platform [name]        : OpenMM platform. Defaults to fastest available")
    print("    --command ([command])    : Show help for an input-file command")
    print("    --verbose                : Show full simulation configuration before run")
    print("    --ff                     : Show supported force-field parameter files")
    print("    --license                : Show the software license")
    print("    --help")
    print("")
    print(" HTMD Molecular Dynamics incorporates OpenMM. See http://wiki.simtk.org/openmm/License")
    print(" No re-distribution in whole or part")
    print("")


def main_cli():
    input_file = None
    device = None
    verbose = False
    platform = None
    debug = False
    check_registration(product="md")
    a = 1
    try:
        while (a < len(sys.argv)):
            if (sys.argv[a] == "--input"):
                if a < len(sys.argv) - 1:
                    input_file = sys.argv[a + 1]
                    a = a + 1
                else:
                    raise NameError()
            elif (sys.argv[a] == "--ff" or sys.argv[a] == "--forcefield"):
                show_ffs()
                sys.exit(0)
            elif (sys.argv[a] == "--license"):
                show_license(product="md")
                sys.exit(0)
            elif (sys.argv[a] == "--verbose"):
                verbose = True;
                a = a + 1
            elif (sys.argv[a] == "--device"):
                if a < len(sys.argv) - 1:
                    device = int(sys.argv[a + 1])
                    a = a + 1
                else:
                    raise NameError("--device  requires an argument")
            elif (sys.argv[a] == "--platform"):
                if (a < len(sys.argv) - 1):
                    platform = sys.argv[a + 1]
                    a = a + 1
                else:
                    show_platforms()
                    sys.exit(0)
            elif (sys.argv[a] == "--command"):
                if a < len(sys.argv) - 1:
                    Command.help(sys.argv[a + 1])
                    sys.exit(1)
                    a = a + 1
                else:
                    Command.help(None)
                    sys.exit(1)

            elif (sys.argv[a] == "--help" or sys.argv[a] == "-h"):
                syntax()
                sys.exit(0)
            elif (sys.argv[a] == "--debug"):
                a = a + 1
                debug = True
            else:
                if not input_file:
                    input_file = sys.argv[a]
                else:
                    raise NameError()
            a = a + 1
    except NameError as e:
        raise (e)
        syntax()
        sys.exit(1)

    if (not input_file):
        input_file = "input"

    if not input_file or not os.path.isfile(input_file):
        syntax();
        sys.exit(1)

    print("\n === HTMD Molecular Dynamics 2016 ===")
    print("      (c) Acellera\n")
    try:
        config = Configuration(config=input_file)
    except NameError as e:
        print("\n Failed to parse input file : \n\n  " + str(e) + "\n\n")
        if debug: raise e
        sys.exit(2)

    if verbose:
        print("\n Simulation configuration:")
        print(config)

    try:
        Simulation(config, device=device, platform=platform)
    except NameError as e:
        print("\n Failed to run simulation : " + str(e) + "\n")
        if debug: raise e
        sys.exit(3)
    except ValueError as e:
        print("\n Error setting up the simulation : " + str(e) + "\n")
        if debug: raise e
        sys.exit(4)
    except KeyboardInterrupt as e:
        print("\n Terminating at user request before simulation started\n");
        if debug: raise e
        sys.exit(5)

    sys.exit(0)


def show_platforms():
    for i in range(Platform.getNumPlatforms()):
        print(Platform.getPlatform(i).getName())


def show_ffs():
    import inspect
    import simtk.openmm.app.forcefield
    import os
    import natsort
    import re
    pp = (os.path.dirname(inspect.getfile(simtk.openmm.app.forcefield.ForceField)))
    pp = os.path.join(pp, "data")
    i = 0
    print("")
    for a in (natsort.natsorted(os.listdir(pp))):
        a = re.sub(".xml$", "", a)
        print("%20s" % (a), end="")
        i = i + 1
        if (not i % 4): print("")

    print("\n")


if __name__ == "__main__":
    main_cli()
