__author__ = 'mjharvey'
from simtk.unit import *
import re
from simtk.openmm.vec3 import Vec3


def xsc_write(box, filename):
    f = open(filename, "w")
    print("%f %f %f" % (
    box[0][0].value_in_unit(angstrom), box[1][1].value_in_unit(angstrom), box[2][2].value_in_unit(angstrom)), file=f)
    f.close()


def xsc_read(filename):
    try:
        f = open(filename, "r")
        r = f.readline()
        b = re.split(r'\s+', r)
        box = [
            Vec3(float(b[0]), 0., 0.),
            Vec3(0., float(b[1]), 0.),
            Vec3(0., 0., float(b[2])),
        ]
        f.close()
        return box * angstrom
    except:
        raise NameError("Invalid extendedsystem file:")
