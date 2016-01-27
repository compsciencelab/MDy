"""
This module contains classes for parsing and processing CHARMM parameter,
topology, and stream files. It only extracts atom properties from the
topology files and extracts all parameters from the parameter files

This file is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of Biological
Structures at Stanford, funded under the NIH Roadmap for Medical Research,
grant U54 GM072970. See https://simtk.org.  This code was originally part of
the ParmEd program and was ported for use with OpenMM.

Copyright (c) 2014 the Authors

Author: Jason M. Swails
Contributors:
Date: Sep. 17, 2014

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import absolute_import
import os
from simtk.openmm.app.internal.charmm._charmmfile import (
            CharmmFile, CharmmStreamFile)
from simtk.openmm.app.internal.charmm.topologyobjects import (
            AtomType, BondType, AngleType, DihedralType, ImproperType, CmapType,
            UreyBradleyType, NoUreyBradley)
from simtk.openmm.app.internal.charmm.exceptions import CharmmFileError
from mdy.openmmfix.element import Element, get_by_symbol
import simtk.unit as u
import warnings

class CharmmParameterSet(object):
    """
    Stores a parameter set defined by CHARMM files. It stores the equivalent of
    the information found in the MASS section of the CHARMM topology file
    (TOP/RTF) and all of the information in the parameter files (PAR)

    Parameters:
        - filenames : List of topology, parameter, and stream files to load into
          the parameter set. The following file type suffixes are recognized.
          Unrecognized file types raise a TypeError
            .rtf, .top -- Residue topology file
            .par, .prm -- Parameter file
            .str -- Stream file
            .inp -- If "par" is in the file name, it is a parameter file. If
                    "top" is in the file name, it is a topology file. Otherwise,
                    raise TypeError

    Attributes:
        All type lists are dictionaries whose keys are tuples (with however
        many elements are needed to define that type of parameter). The types
        that can be in any order are SORTED.
        
        - atom_types_str
        - atom_types_int
        - atom_types_tuple
        - bond_types
        - angle_types
        - urey_bradley_types
        - dihedral_types
        - improper_types
        - cmap_types
        - nbfix_types

        The dihedral types can be multiterm, so the values for each dict key is
        actually a list of DihedralType instances. The atom_types are dicts that
        match the name (str), number (int), or (name, number) tuple (tuple) to
        the atom type. The tuple is guaranteed to be the most robust, although
        when only the integer or string is available the other dictionaries are
        helpful

    Example:
    >>> params = CharmmParameterSet('charmm22.top', 'charmm22.par', 'file.str')
    """

    @staticmethod
    def _convert(data, type, msg=''):
        """
        Converts a data type to a desired type, raising CharmmFileError if it
        fails
        """
        try:
            return type(data)
        except ValueError:
            raise CharmmFileError('Could not convert %s to %s' % (msg, type))

    def __init__(self, *args, **kwargs):
        # Instantiate the list types
        self.atom_types_str = dict()
        self.atom_types_int = dict()
        self.atom_types_tuple = dict()
        self.bond_types = dict()
        self.angle_types = dict()
        self.urey_bradley_types = dict()
        self.dihedral_types = dict()
        self.improper_types = dict()
        self.cmap_types = dict()
        self.nbfix_types = dict()
        self.parametersets = []
        self.setup_charmm_types()
        
        # Load all of the files
        tops, pars, strs = [], [], []
        for arg in args:
            if arg.endswith('.rtf') or arg.endswith('.top'):
                tops.append(arg)
            elif arg.endswith('.par') or arg.endswith('.prm'):
                pars.append(arg)
            elif arg.endswith('.str'):
                strs.append(arg)
            elif arg.endswith('.inp'):
                # Only consider the file name (since the directory is likely
                # "toppar" and will screw up file type detection)
                fname = os.path.split(arg)[1]
                if 'par' in fname:
                    pars.append(arg)
                elif 'top' in fname:
                    tops.append(arg)
                else:
                    raise TypeError('Unrecognized file type: %s' % arg)
            else:
                raise TypeError('Unrecognized file type: %s' % arg)

        permissive=kwargs.pop("permissive", False)
        if len(kwargs):
            raise TypeError('Unrecognised named argument')

        for top in tops: self.readTopologyFile(top)
        for par in pars: self.readParameterFile(par, permissive=permissive)
        for strf in strs: self.readStreamFile(strf)

    @classmethod
    def loadSet(cls, tfile=None, pfile=None, sfiles=[], permissive=False):
        """
        Instantiates a CharmmParameterSet from a Topology file and a Parameter
        file (or just a Parameter file if it has all information)

        Parameters:
            - tfile (str) : Name of the Topology (RTF/TOP) file
            - pfile (str) : Name of the Parameter (PAR) file
            - sfiles (list of str) : List or tuple of stream (STR) file names.
            - permissive (bool) : Accept non-bonbded parameters for undefined 
                                  atom types (default False)

        Returns:
            New CharmmParameterSet populated with the parameters found in the
            provided files.
            
        Notes:
            The RTF file is read first (if provided), followed by the PAR file,
            followed by the list of stream files (in the order they are
            provided). Parameters in each stream file will overwrite those that
            came before (or simply append to the existing set if they are
            different)
        """
        inst = cls()
        if tfile is not None:
            inst.readTopologyFile(tfile)
        if pfile is not None:
            inst.readParameterFile(pfile, permissive=permissive)
        if isinstance(sfiles, str):
            # The API docstring requests a list, but allow for users to pass a
            # string with a single filename instead
            inst.readStreamFile(sfiles)
        elif sfiles is not None:
            for sfile in sfiles:
                inst.readStreamFile(sfile)
        return inst

    def readParameterFile(self, pfile, permissive=False):
        """
        Reads all of the parameters from a parameter file. Versions 36 and
        later of the CHARMM force field files have an ATOMS section defining
        all of the atom types.  Older versions need to load this information
        from the RTF/TOP files.

        Parameters:
            - pfile (str) : Name of the CHARMM PARameter file to read
            - permissive (bool) : Accept non-bonbded parameters for undefined 
                                  atom types (default False)

        Notes: The atom types must all be loaded by the end of this routine.
        Either supply a PAR file with atom definitions in them or read in a
        RTF/TOP file first. Failure to do so will result in a raised
        RuntimeError.
        """
        conv = CharmmParameterSet._convert
        if isinstance(pfile, str):
            own_handle = True
            f = CharmmFile(pfile)
        else:
            own_handle = False
            f = pfile
        # What section are we parsing?
        section = None
        # The current cmap we are building (these span multiple lines)
        current_cmap = None
        current_cmap_data = []
        current_cmap_res = 0
        nonbonded_types = dict() # Holder
        parameterset = None
        read_first_nonbonded = False
        for line in f:
            line = line.strip()
            if not line:
                # This is a blank line
                continue
            if parameterset is None and line.strip().startswith('*>>'):
                parameterset = line.strip()[1:78]
                continue
            # Set section if this is a section header
            if line.startswith('ATOMS'):
                section = 'ATOMS'
                continue
            if line.startswith('BONDS'):
                section = 'BONDS'
                continue
            if line.startswith('ANGLES'):
                section = 'ANGLES'
                continue
            if line.startswith('DIHEDRALS'):
                section = 'DIHEDRALS'
                continue
            if line.startswith('IMPROPER'):
                section = 'IMPROPER'
                continue
            if line.startswith('CMAP'):
                section = 'CMAP'
                continue
            if line.startswith('NONBONDED'):
                read_first_nonbonded = False
                section = 'NONBONDED'
                continue
            if line.startswith('NBFIX'):
                section = 'NBFIX'
                continue
            if line.startswith('HBOND'):
                section = None
                continue
            # It seems like files? sections? can be terminated with 'END'
            if line.startswith('END'): # should this be case-insensitive?
                section = None
                continue
            # If we have no section, skip
            if section is None: continue
            # Now handle each section specifically
            if section == 'ATOMS':
                if not line.startswith('MASS'): continue # Should this happen?
                words = line.split()
                try:
                    idx = conv(words[1], int, 'atom type')
                    name = words[2]
                    mass = conv(words[3], float, 'atom mass')
                except IndexError:
                    raise CharmmFileError('Could not parse MASS section.')
                # The parameter file might or might not have an element name
                try:
                    elem = words[4]
                    atomic_number = get_by_symbol(elem).atomic_number
                except (IndexError, KeyError):
                    # Figure it out from the mass
                    masselem = Element.getByMass(mass)
                    if masselem is None:
                        atomic_number = 0 # Extra point or something
                    else:
                        atomic_number = masselem.atomic_number
                atype = AtomType(name=name, number=idx, mass=mass,
                                 atomic_number=atomic_number)
                self.atom_types_str[atype.name] = atype
                self.atom_types_int[atype.number] = atype
                self.atom_types_tuple[(atype.name, atype.number)] = atype
                continue
            if section == 'BONDS':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    k = conv(words[2], float, 'bond force constant')
                    req = conv(words[3], float, 'bond equilibrium dist')
                except IndexError:
                    raise CharmmFileError('Could not parse bonds.')
                key = (min(type1, type2), max(type1, type2))
                self.bond_types[key] = BondType(k, req)
                continue
            if section == 'ANGLES':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    k = conv(words[3], float, 'angle force constant')
                    theteq = conv(words[4], float, 'angle equilibrium value')
                except IndexError:
                    raise CharmmFileError('Could not parse angles.')
                key = (min(type1, type3), type2, max(type1, type3))
                self.angle_types[key] = AngleType(k, theteq)
                # See if we have a urey-bradley
                try:
                    ubk = conv(words[5], float, 'Urey-Bradley force constant')
                    ubeq = conv(words[6], float, 'Urey-Bradley equil. value')
                    ubtype = UreyBradleyType(ubk, ubeq)
                except IndexError:
                    ubtype = NoUreyBradley
                self.urey_bradley_types[key] = ubtype
                continue
            if section == 'DIHEDRALS':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    type4 = words[3]
                    k = conv(words[4], float, 'dihedral force constant')
                    n = conv(words[5], float, 'dihedral periodicity')
                    phase = conv(words[6], float, 'dihedral phase')
                except IndexError:
                    raise CharmmFileError('Could not parse dihedrals.')
                # Torsion can be in either direction. Sort by end groups first,
                # then sort by middle 2
                if type1 < type4:
                    key = (type1, type2, type3, type4)
                elif type1 > type4:
                    key = (type4, type3, type2, type1)
                else:
                    # OK, we need to sort by the middle atoms now
                    if type2 < type3:
                        key = (type1, type2, type3, type4)
                    else:
                        key = (type4, type3, type2, type1)
                # See if this is a second (or more) term of the dihedral group
                # that's already present.
                dihedral = DihedralType(k, n, phase)
                if key in self.dihedral_types:
                    # See if the existing dihedral type list has a term with
                    # the same periodicity -- If so, replace it
                    replaced = False
                    for i, dtype in enumerate(self.dihedral_types[key]):
                        if dtype.per == dihedral.per:
                            # Replace. Warn if they are different
                            if dtype != dihedral:
                                warnings.warn('Replacing dihedral %r with %r' % 
                                              (dtype, dihedral))
                            self.dihedral_types[key]
                            replaced = True
                            break
                    if not replaced:
                        self.dihedral_types[key].append(dihedral)
                else: # key not present
                    self.dihedral_types[key] = [dihedral]
                continue
            if section == 'IMPROPER':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    type4 = words[3]
                    k = conv(words[4], float, 'improper force constant')
                    theteq = conv(words[5], float, 'improper equil. value')
                except IndexError:
                    raise CharmmFileError('Could not parse dihedrals.')
                # If we have a 7th column, that is the real psi0 (and the 6th
                # is just a dummy 0)
                try:
                    tmp = conv(words[6], float, 'improper equil. value')
                    theteq = tmp
                except IndexError:
                    pass # Do nothing
                # Improper types seem not to have the central atom defined in
                # the first place, so just have the key a fully sorted list. We
                # still depend on the PSF having properly ordered improper atoms
                key = tuple(sorted([type1, type2, type3, type4]))
                self.improper_types[key] = ImproperType(k, theteq)
                continue
            if section == 'CMAP':
                # This is the most complicated part, since cmap parameters span
                # many lines. We won't do much error catching here.
                words = line.split()
                try:
                    holder = [float(w) for w in words]
                    current_cmap_data.extend(holder)
                except ValueError:
                    # We assume this is a definition of a new CMAP, so
                    # terminate the last CMAP if applicable
                    if current_cmap is not None:
                        # We have a map to terminate
                        ty = CmapType(current_cmap_res, current_cmap_data)
                        self.cmap_types[current_cmap] = ty
                    try:
                        type1 = words[0]
                        type2 = words[1]
                        type3 = words[2]
                        type4 = words[3]
                        type5 = words[4]
                        type6 = words[5]
                        type7 = words[6]
                        type8 = words[7]
                        res = conv(words[8], int, 'CMAP resolution')
                    except IndexError:
                        raise CharmmFileError('Could not parse CMAP data.')
                    # order the torsions independently
                    k1 = [type1, type2, type3, type4]
                    k2 = [type4, type3, type2, type1]
                    key1 = min(k1, k2)
                    k1 = [type5, type6, type7, type8]
                    k2 = [type8, type7, type6, type5]
                    key2 = min(k1, k2)
                    current_cmap = tuple(key1 + key2)
                    current_cmap_res = res
                    current_cmap_data = []
                continue
            if section == 'NONBONDED':
                # Now get the nonbonded values
                words = line.split()
                try:
                    atype = words[0]
                    # 1st column is ignored
                    epsilon = conv(words[2], float, 'vdW epsilon term')
                    rmin = conv(words[3], float, 'vdW Rmin/2 term')
                except IndexError:
                    # If we haven't read our first nonbonded term yet, we may
                    # just be parsing the settings that should be used. So
                    # soldier on
                    if not read_first_nonbonded: continue
                    raise CharmmFileError('Could not parse nonbonded terms.')
                except CharmmFileError as e:
                    if not read_first_nonbonded: continue
                    raise CharmmFileError(str(e))
                else:
                    # OK, we've read our first nonbonded section for sure now
                    read_first_nonbonded = True
                # See if we have 1-4 parameters
                try:
                    # 4th column is ignored
                    eps14 = conv(words[5], float, '1-4 vdW epsilon term')
                    rmin14 = conv(words[6], float, '1-4 vdW Rmin/2 term')
                except IndexError:
                    eps14 = rmin14 = None
                nonbonded_types[atype] = [epsilon, rmin, eps14, rmin14]
                continue
            if section == 'NBFIX':
                words = line.split()
                try:
                    at1 = words[0]
                    at2 = words[1]
                    emin = abs(conv(words[2], float, 'NBFIX Emin'))
                    rmin = conv(words[3], float, 'NBFIX Rmin')
                    try:
                        emin14 = abs(conv(words[4], float, 'NBFIX Emin 1-4'))
                        rmin14 = conv(words[5], float, 'NBFIX Rmin 1-4')
                    except IndexError:
                        emin14 = rmin14 = None
                    try:
                        self.atom_types_str[at1].add_nbfix(at2, rmin, emin,
                                                           rmin14, emin14)
                        self.atom_types_str[at2].add_nbfix(at1, rmin, emin,
                                                           rmin14, emin14)
                    except KeyError:
                        # Some stream files define NBFIX terms with an atom that
                        # is defined in another toppar file that does not
                        # necessarily have to be loaded. As a result, not every
                        # NBFIX found here will necessarily need to be applied.
                        # If we can't find a particular atom type, don't bother
                        # adding that nbfix and press on
                        pass
                except IndexError:
                    raise CharmmFileError('Could not parse NBFIX terms.')
                self.nbfix_types[(min(at1, at2), max(at1, at2))] = (emin, rmin)
        # If there were any CMAP terms stored in the parameter set, the last one
        # defined will not have been added to the set. Add it now.
        if current_cmap is not None:
            ty = CmapType(current_cmap_res, current_cmap_data)
            self.cmap_types[current_cmap] = ty

        # If in permissive mode create an atomtype for every type used in  
        # the nonbonded parameters. This is a work-around for when all that's
        # available is a CHARMM22 inp file, which has no ATOM/MASS fields

        if permissive:
            try:
               idx = max(self.atom_types_int.keys())+1000
            except ValueError:
               idx = 10000
            for key in nonbonded_types:
                if not key in self.atom_types_str:
                  if key in self.charmm_types:
                      tt = self.charmm_types[key]
                      self.atom_types_str[key] = tt
                      self.atom_types_int[tt.number] = tt
                  else:
                    atype =AtomType(name=key, number=idx, mass= float('NaN'), atomic_number= 0 )
                    self.atom_types_str[key] = atype 
                    self.atom_types_int[idx] = atype
                    idx=idx+1

        # Now we're done. Load the nonbonded types into the relevant AtomType
        # instances. In order for this to work, all keys in nonbonded_types
        # must be in the self.atom_types_str dict. Raise a RuntimeError if this
        # is not satisfied
        try:
            for key in nonbonded_types:
                self.atom_types_str[key].set_lj_params(*nonbonded_types[key])
        except KeyError:
            raise RuntimeError('Atom type %s not present in AtomType list' %
                               key)

        if parameterset is not None: self.parametersets.append(parameterset)
        if own_handle: f.close()

    def readTopologyFile(self, tfile):
        """
        Reads _only_ the atom type definitions from a topology file. This is
        unnecessary for versions 36 and later of the CHARMM force field.

        Parameters:
            - tfile (str) : Name of the CHARMM TOPology file to read

        Note: The CHARMM TOPology file is also called a Residue Topology File
        """
        conv = CharmmParameterSet._convert
        if isinstance(tfile, str):
            own_handle = True
            f = CharmmFile(tfile)
        else:
            own_handle = False
            f = tfile
        for line in f:
            line = line.strip()
            if line[:4] != 'MASS': continue
            words = line.split()
            try:
                idx = conv(words[1], int, 'atom type')
                name = words[2]
                mass = conv(words[3], float, 'atom mass')
            except IndexError:
                raise CharmmFileError('Could not parse MASS section of %s' %
                                      tfile)
            # The parameter file might or might not have an element name
            try:
                elem = words[4]
                atomic_number = get_by_symbol(elem).atomic_number
            except (IndexError, KeyError):
                # Figure it out from the mass
#                print("MJH --- mass -- " + str(mass) )
                masselem = Element.getByMass(mass)
                if masselem is None:
                    atomic_number = 0 # Extra point or something
                else:
                    atomic_number = masselem.atomic_number
            atype = AtomType(name=name, number=idx, mass=mass,
                             atomic_number=atomic_number)
            self.atom_types_str[atype.name] = atype
            self.atom_types_int[atype.number] = atype
            self.atom_types_tuple[(atype.name, atype.number)] = atype
        if own_handle: f.close()

    def readStreamFile(self, sfile):
        """
        Reads RTF and PAR sections from a stream file and dispatches the
        sections to readTopologyFile or readParameterFile

        Parameters:
            - sfile (str or CharmmStreamFile) : Stream file to parse
        """
        if isinstance(sfile, CharmmStreamFile):
            f = sfile
        else:
            f = CharmmStreamFile(sfile)

        title, section = f.next_section()
        while title is not None and section is not None:
            words = title.lower().split()
            if words[1] == 'rtf':
                # This is a Residue Topology File section.
                self.readTopologyFile(section)
            elif words[1].startswith('para'):
                # This is a Parameter file section
                self.readParameterFile(section)
            title, section = f.next_section()

    def condense(self):
        """
        This function goes through each of the parameter type dicts and
        eliminates duplicate types. After calling this function, every unique
        bond, angle, dihedral, improper, or cmap type will pair with EVERY key
        in the type mapping dictionaries that points to the equivalent type

        Returns:
            - Returns the instance that is being condensed.

        Notes:
            The return value allows you to condense the types at construction
            time.

        Example:
        >>> params = CharmmParameterSet('charmm.prm').condense()
        """
        # First scan through all of the bond types
        self._condense_types(self.bond_types)
        self._condense_types(self.angle_types)
        self._condense_types(self.urey_bradley_types)
        self._condense_types(self.improper_types)
        self._condense_types(self.cmap_types)
        # Dihedrals have to be handled separately, since each key corresponds to
        # a list of (potentially multiterm) dihedral terms. Since all terms in a
        # multiterm dihedral have to have a DIFFERENT periodicity, we don't have
        # to condense _within_ a single list of torsions assigned to the same
        # key (they're guaranteed to be different)
        keylist = list(self.dihedral_types.keys())
        for i in range(len(keylist) - 1):
            key1 = keylist[i]
            for dihedral in self.dihedral_types[key1]:
                for j in range(i+1, len(keylist)):
                    key2 = keylist[j]
                    for jj, dihedral2 in enumerate(self.dihedral_types[key2]):
                        if dihedral2 == dihedral:
                            self.dihedral_types[key2][jj] = dihedral
        return self

    @staticmethod
    def _condense_types(typedict):
        """
        Loops through the given dict and condenses all types.

        Parameter:
            - typedict : Type dictionary to condense
        """
        keylist = list(typedict.keys())
        for i in range(len(keylist) - 1):
            key1 = keylist[i]
            for j in range(i+1, len(keylist)):
                key2 = keylist[j]
                if typedict[key1] == typedict[key2]:
                    typedict[key2] = typedict[key1]

    def setup_charmm_types( self ):
        self.charmm_types = dict()
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        self.charmm_types["H"]  = AtomType(name="H", number=1, mass=1.00800, atomic_number=1 )
        self.charmm_types["HC"]  = AtomType(name="HC", number=2, mass=1.00800, atomic_number=1 )
        self.charmm_types["HA"]  = AtomType(name="HA", number=3, mass=1.00800, atomic_number=1 )
        self.charmm_types["HT"]  = AtomType(name="HT", number=4, mass=1.00800, atomic_number=1 )
        self.charmm_types["HP"]  = AtomType(name="HP", number=5, mass=1.00800, atomic_number=1 )
        self.charmm_types["HB"]  = AtomType(name="HB", number=6, mass=1.00800, atomic_number=1 )
        self.charmm_types["HR1"]  = AtomType(name="HR1", number=7, mass=1.00800, atomic_number=1 )
        self.charmm_types["HR2"]  = AtomType(name="HR2", number=8, mass=1.00800, atomic_number=1 )
        self.charmm_types["HR3"]  = AtomType(name="HR3", number=9, mass=1.00800, atomic_number=1 )
        self.charmm_types["HS"]  = AtomType(name="HS", number=10, mass=1.00800, atomic_number=1 )
        self.charmm_types["HE1"]  = AtomType(name="HE1", number=11, mass=1.00800, atomic_number=1 )
        self.charmm_types["HE2"]  = AtomType(name="HE2", number=12, mass=1.00800, atomic_number=1 )
        self.charmm_types["HA1"]  = AtomType(name="HA1", number=13, mass=1.00800, atomic_number=1 )
        self.charmm_types["HA2"]  = AtomType(name="HA2", number=14, mass=1.00800, atomic_number=1 )
        self.charmm_types["HA3"]  = AtomType(name="HA3", number=15, mass=1.00800, atomic_number=1 )
        self.charmm_types["HF1"]  = AtomType(name="HF1", number=16, mass=1.00800, atomic_number=1 )
        self.charmm_types["HF2"]  = AtomType(name="HF2", number=17, mass=1.00800, atomic_number=1 )
        self.charmm_types["C"]  = AtomType(name="C", number=20, mass=12.01100, atomic_number=6 )
        self.charmm_types["CA"]  = AtomType(name="CA", number=21, mass=12.01100, atomic_number=6 )
        self.charmm_types["CT1"]  = AtomType(name="CT1", number=22, mass=12.01100, atomic_number=6 )
        self.charmm_types["CT2"]  = AtomType(name="CT2", number=23, mass=12.01100, atomic_number=6 )
        self.charmm_types["CT3"]  = AtomType(name="CT3", number=24, mass=12.01100, atomic_number=6 )
        self.charmm_types["CPH1"]  = AtomType(name="CPH1", number=25, mass=12.01100, atomic_number=6 )
        self.charmm_types["CPH2"]  = AtomType(name="CPH2", number=26, mass=12.01100, atomic_number=6 )
        self.charmm_types["CPT"]  = AtomType(name="CPT", number=27, mass=12.01100, atomic_number=6 )
        self.charmm_types["CY"]  = AtomType(name="CY", number=28, mass=12.01100, atomic_number=6 )
        self.charmm_types["CP1"]  = AtomType(name="CP1", number=29, mass=12.01100, atomic_number=6 )
        self.charmm_types["CP2"]  = AtomType(name="CP2", number=30, mass=12.01100, atomic_number=6 )
        self.charmm_types["CP3"]  = AtomType(name="CP3", number=31, mass=12.01100, atomic_number=6 )
        self.charmm_types["CC"]  = AtomType(name="CC", number=32, mass=12.01100, atomic_number=6 )
        self.charmm_types["CD"]  = AtomType(name="CD", number=33, mass=12.01100, atomic_number=6 )
        self.charmm_types["CPA"]  = AtomType(name="CPA", number=34, mass=12.01100, atomic_number=6 )
        self.charmm_types["CPB"]  = AtomType(name="CPB", number=35, mass=12.01100, atomic_number=6 )
        self.charmm_types["CPM"]  = AtomType(name="CPM", number=36, mass=12.01100, atomic_number=6 )
        self.charmm_types["CM"]  = AtomType(name="CM", number=37, mass=12.01100, atomic_number=6 )
        self.charmm_types["CS"]  = AtomType(name="CS", number=38, mass=12.01100, atomic_number=6 )
        self.charmm_types["CE1"]  = AtomType(name="CE1", number=39, mass=12.01100, atomic_number=6 )
        self.charmm_types["CE2"]  = AtomType(name="CE2", number=40, mass=12.01100, atomic_number=6 )
        self.charmm_types["CST"]  = AtomType(name="CST", number=41, mass=12.01100, atomic_number=6 )
        self.charmm_types["CT"]  = AtomType(name="CT", number=42, mass=12.01100, atomic_number=6 )
        self.charmm_types["CT1x"]  = AtomType(name="CT1x", number=43, mass=12.01100, atomic_number=6 )
        self.charmm_types["CT2x"]  = AtomType(name="CT2x", number=44, mass=12.01100, atomic_number=6 )
        self.charmm_types["CT3x"]  = AtomType(name="CT3x", number=45, mass=12.01100, atomic_number=6 )
        self.charmm_types["CN"]  = AtomType(name="CN", number=46, mass=12.01100, atomic_number=6 )
        self.charmm_types["CAP"]  = AtomType(name="CAP", number=47, mass=12.01100, atomic_number=6 )
        self.charmm_types["COA"]  = AtomType(name="COA", number=48, mass=12.01100, atomic_number=6 )
        self.charmm_types["C3"]  = AtomType(name="C3", number=49, mass=12.01100, atomic_number=6 )
        self.charmm_types["N"]  = AtomType(name="N", number=50, mass=14.00700, atomic_number=7 )
        self.charmm_types["NR1"]  = AtomType(name="NR1", number=51, mass=14.00700, atomic_number=7 )
        self.charmm_types["NR2"]  = AtomType(name="NR2", number=52, mass=14.00700, atomic_number=7 )
        self.charmm_types["NR3"]  = AtomType(name="NR3", number=53, mass=14.00700, atomic_number=7 )
        self.charmm_types["NH1"]  = AtomType(name="NH1", number=54, mass=14.00700, atomic_number=7 )
        self.charmm_types["NH2"]  = AtomType(name="NH2", number=55, mass=14.00700, atomic_number=7 )
        self.charmm_types["NH3"]  = AtomType(name="NH3", number=56, mass=14.00700, atomic_number=7 )
        self.charmm_types["NC2"]  = AtomType(name="NC2", number=57, mass=14.00700, atomic_number=7 )
        self.charmm_types["NY"]  = AtomType(name="NY", number=58, mass=14.00700, atomic_number=7 )
        self.charmm_types["NP"]  = AtomType(name="NP", number=59, mass=14.00700, atomic_number=7 )
        self.charmm_types["NPH"]  = AtomType(name="NPH", number=60, mass=14.00700, atomic_number=7 )
        self.charmm_types["NC"]  = AtomType(name="NC", number=61, mass=14.00700, atomic_number=7 )
        self.charmm_types["O"]  = AtomType(name="O", number=70, mass=15.99900, atomic_number=8 )
        self.charmm_types["OB"]  = AtomType(name="OB", number=71, mass=15.99900, atomic_number=8 )
        self.charmm_types["OC"]  = AtomType(name="OC", number=72, mass=15.99900, atomic_number=8 )
        self.charmm_types["OH1"]  = AtomType(name="OH1", number=73, mass=15.99900, atomic_number=8 )
        self.charmm_types["OS"]  = AtomType(name="OS", number=74, mass=8, atomic_number=8 )
        self.charmm_types["OT"]  = AtomType(name="OT", number=75, mass=8, atomic_number=8 )
        self.charmm_types["OM"]  = AtomType(name="OM", number=76, mass=15.99900, atomic_number=8 )
        self.charmm_types["OST"]  = AtomType(name="OST", number=77, mass=15.99900, atomic_number=8 )
        self.charmm_types["OCA"]  = AtomType(name="OCA", number=78, mass=15.99900, atomic_number=8 )
        self.charmm_types["S"]  = AtomType(name="S", number=81, mass=16, atomic_number=16 )
        self.charmm_types["SM"]  = AtomType(name="SM", number=82, mass=16, atomic_number=16 )
        self.charmm_types["SS"]  = AtomType(name="SS", number=83, mass=16, atomic_number=16 )
        self.charmm_types["HE"]  = AtomType(name="HE", number=85, mass=2, atomic_number=2 )
        self.charmm_types["NE"]  = AtomType(name="NE", number=86, mass=10, atomic_number=10 )
        self.charmm_types["CF1"]  = AtomType(name="CF1", number=87, mass=12.01100, atomic_number=6 )
        self.charmm_types["CF2"]  = AtomType(name="CF2", number=88, mass=12.01100, atomic_number=6 )
        self.charmm_types["CF3"]  = AtomType(name="CF3", number=89, mass=12.01100, atomic_number=6 )
        self.charmm_types["FE"]  = AtomType(name="FE", number=90, mass=26, atomic_number=26 )
        self.charmm_types["CLAL"]  = AtomType(name="CLAL", number=91, mass=17, atomic_number=17 )
        self.charmm_types["FA"]  = AtomType(name="FA", number=92, mass=9, atomic_number=9 )
        self.charmm_types["F1"]  = AtomType(name="F1", number=93, mass=9, atomic_number=9 )
        self.charmm_types["F2"]  = AtomType(name="F2", number=94, mass=9, atomic_number=9 )
        self.charmm_types["F3"]  = AtomType(name="F3", number=95, mass=9, atomic_number=9 )
        self.charmm_types["DUM"]  = AtomType(name="DUM", number=99, mass=0.00000, atomic_number=0.00000 )
        self.charmm_types["SOD"]  = AtomType(name="SOD", number=100, mass=11, atomic_number=11 )
        self.charmm_types["MG"]  = AtomType(name="MG", number=101, mass=12, atomic_number=12 )
        self.charmm_types["POT"]  = AtomType(name="POT", number=102, mass=19, atomic_number=19 )
        self.charmm_types["CES"]  = AtomType(name="CES", number=103, mass=55, atomic_number=55 )
        self.charmm_types["CAL"]  = AtomType(name="CAL", number=104, mass=20, atomic_number=20 )
        self.charmm_types["CLA"]  = AtomType(name="CLA", number=105, mass=17, atomic_number=17 )
        self.charmm_types["ZN"]  = AtomType(name="ZN", number=106, mass=30, atomic_number=30 )
        self.charmm_types["CC1A"]  = AtomType(name="CC1A", number=112, mass=12.01100, atomic_number=6 )
        self.charmm_types["CC1B"]  = AtomType(name="CC1B", number=113, mass=12.01100, atomic_number=6 )
        self.charmm_types["CC2"]  = AtomType(name="CC2", number=114, mass=12.01100, atomic_number=6 )
        self.charmm_types["NS1"]  = AtomType(name="NS1", number=120, mass=14.00700, atomic_number=7 )
        self.charmm_types["NS2"]  = AtomType(name="NS2", number=121, mass=14.00700, atomic_number=7 )
