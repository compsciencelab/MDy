from simtk.unit import *
from simtk.openmm.vec3 import *
import numpy
import os

import struct




def bincoor_read(filename, unit):
    f = open(filename, 'rb')
    data = bytearray(f.read())
    f.close()


    n = struct.unpack( 'i' * 1     , data[0:4] )
    natoms=n[0]

    coords = struct.unpack( 'd' * natoms * 3, data[4:len(data)] )

    vec3coords=[]
    for i in range(natoms):
        vec3coords.append( Vec3( coords[i*3+0] , coords[i*3+1] , coords[i*3+2] ) )

    vec3coords = vec3coords * unit
    #print (vec3coords)

    return vec3coords





def bincoor_write(vec3coords, filename, unit):
    natoms = [ 0 ]
    natoms[0]  = len(vec3coords)
    coords     = numpy.empty(( natoms[0], 3), dtype=numpy.float64)

    for i in range(natoms[0]):


        coords[i,0] = vec3coords[i][0].value_in_unit( unit )
        coords[i,1] = vec3coords[i][1].value_in_unit( unit )
        coords[i,2] = vec3coords[i][2].value_in_unit( unit )

    f = open(filename, 'wb')

    dat = coords[:, :]
    dat = dat.reshape(dat.shape[0] * 3).astype(numpy.float64)

    fmt1 = 'i' * 1  #natoms[0]
    bin1 = struct.pack(fmt1, *natoms)
    fmt2 = 'd' * dat.shape[0]
    bin2 = struct.pack(fmt2, *dat)
    f.write(bin1)
    f.write(bin2)
    f.close()


if __name__ == "__main__":
    x = bincoor_read( 'test.coor')
    bincoor_write( x, 'out.coor' )
