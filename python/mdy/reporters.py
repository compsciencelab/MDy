from simtk.unit import *
from simtk.openmm import *
from simtk.unit import nanometer
from simtk.unit import angstrom
from simtk.unit import picosecond
from simtk.unit import femtosecond

from htmd.molecule.xtc import XTCwrite
import numpy as np
import time
import math
import os
import sys
from htmd.molecule.molecule import Molecule


class StdoutLogReporter:
    def __init__(self, reportInterval, totalSteps):
        self._reportInterval = reportInterval
        self._inited = False
        self._totalSteps = totalSteps
        self._lastvol = None

    def headers(self):
        print("  %10s %11s %11s %11s %8s %8s %8s %8s %11s" % (
        "Step", "PE", "KE", "Total E", "Temp", "Volume", "Fluct.", "Speed", "Completion"))
        print("  %10s %11s %11s %11s %8s %8s %8s %8s %11s" % (
        "", "kJ/mol", "kJ/mol", "kJ/mol", "K", "nm^3", "%", "ns/day", "dd:hh:mm:ss"))

    def _init(self, simulation, system, state):
        # Compute the number of degrees of freedom.
        dof = 0
        for i in range(system.getNumParticles()):
            if system.getParticleMass(i) > 0 * unit.dalton:
                dof += 3
        dof -= system.getNumConstraints()
        if any(type(system.getForce(i)) == CMMotionRemover for i in range(system.getNumForces())):
            dof -= 3
        self._dof = dof
        self._initialSteps = simulation.currentStep
        self._initialClockTime = time.time()
        self._initialSimulationTime = state.getTime()
        self._inited = True

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, True)

    def report(self, simulation, state):
        if not self._inited:
            self._init(simulation, simulation.system, state)
        clockTime = time.time()
        step = simulation.currentStep;
        #timex = state.getTime().value_in_unit(picosecond)
        pe = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        ke = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
        te = pe + ke
        temp = (2 * state.getKineticEnergy() / (self._dof * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)

        elapsedSeconds = clockTime - self._initialClockTime
        elapsedDays = 0
        if (elapsedSeconds):
            elapsedDays = (elapsedSeconds) / 86400.0
        elapsedSteps = simulation.currentStep - self._initialSteps
        elapsedNs = (state.getTime() - self._initialSimulationTime).value_in_unit(unit.nanosecond)
        if (elapsedDays):
            speed = elapsedNs / elapsedDays

        if elapsedSteps:
            estimatedTotalSeconds = (self._totalSteps - self._initialSteps) * elapsedSeconds / elapsedSteps
        else:
            estimatedTotalSeconds = 0
        remainingSeconds = int(estimatedTotalSeconds - elapsedSeconds)
        remainingDays = remainingSeconds // 86400
        remainingSeconds -= remainingDays * 86400
        remainingHours = remainingSeconds // 3600
        remainingSeconds -= remainingHours * 3600
        remainingMinutes = remainingSeconds // 60
        remainingSeconds -= remainingMinutes * 60
        remianingString = "--"
        if remainingDays > 0:
            remainingString = "%d:%d:%02d:%02d" % (remainingDays, remainingHours, remainingMinutes, remainingSeconds)
        elif remainingHours > 0:
            remainingString = "%d:%02d:%02d" % (remainingHours, remainingMinutes, remainingSeconds)
        elif remainingMinutes > 0:
            remainingString = "%d:%02d" % (remainingMinutes, remainingSeconds)
        else:
            remainingString = "0:%02d" % remainingSeconds

        box = state.getPeriodicBoxVectors()
        v = box[0][0] * box[1][1] * box[2][2]
        volume = v.value_in_unit(nanometer ** 3)
        #    print(pe)
        #    print(ke)
        #    print(te)
        #    print(temp)
        if self._lastvol:
            fluctuation = 100. * (volume / self._lastvol - 1.)
        else:
            fluctuation = 0.
        self._lastvol = volume

        print("# %10ld %11.2f %11.2f %11.2f %8.2f %8.2f %8.2f %8.2f %11s" % (
        step, pe, ke, te, temp, volume, fluctuation, speed, remainingString))

        if (math.isnan(pe) or math.isnan(ke) or math.isnan(temp)):
            raise ValueError("Simulation has become unstable. Aborted")


class ConstraintsDecay:
    def __init__(self, k, rate, customforce, dt):
        self._k = k
        self._rate = rate
        self._customforce = customforce
        self._dt = dt
        self._reportInterval = int(1000 / dt)  # aim up update every ps
        self._done = False

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        if self._done: return

        kk = self._k - self._rate * (simulation.currentStep / self._reportInterval)
        if (kk < 0.):
            kk = 0
            self._done = True
        self._customforce.setGlobalParameterDefaultValue(0, kk * kilocalories_per_mole / angstroms ** 2)
        self._customforce.updateParametersInContext(simulation.context)


class XTCReporter:
    def __init__(self, file, reportInterval):
        self._reportInterval = reportInterval
        self._out = file

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters:
         - simulation (Simulation) The Simulation to generate a report for
        Returns: A five element tuple.  The first element is the number of steps until the
        next report.  The remaining elements specify whether that report will require
        positions, velocities, forces, and energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False)

    def report(self, simulation, state):
        step = simulation.currentStep
        time = step * simulation.integrator.getStepSize().value_in_unit(femtosecond)
        a, b, c = state.getPeriodicBoxVectors()
        coords = state.getPositions(asNumpy=True).value_in_unit(angstrom)

        box = np.empty((3, 1))
        box[0, 0] = a[0].value_in_unit(nanometer)
        box[1, 0] = b[1].value_in_unit(nanometer)
        box[2, 0] = c[2].value_in_unit(nanometer)
        XTCwrite(coords.reshape(coords.shape[0], 3, 1), box, self._out, time=[time], step=[step])


class FlatbottomUpdate:
    def __init__(self, forceobj, config):
        self._flatbottom_force = forceobj
        self._config = config
        self._reportInterval = int(1000 / config.timestep)  # aim up update every ps

        f = config.coordinates
        if config.flatbottom_file:
            f = config.flatbottom_file

        self._molecule = Molecule(f)

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False)

    def report(self, simulation, state):

        #if self._done: return
        update = False
        if self._config.flatbottom_decay:
            kk = self._config.flatbottom_strength - self._config.flatbottom_decay * (
            simulation.currentStep / self._reportInterval)
            if (kk < 0.):
                kk = 0
                self._done = True
            self._flatbottom_force.setGlobalParameterDefaultValue(0, kk * kilocalories_per_mole / angstroms ** 2)
            update = True

        # IF flatbottom_reference is set, then the geometric centre of the
        # potential well is time-variant. Update date it

        if self._config.flatbottom_reference:
            coords = state.getPositions(asNumpy=True).value_in_unit(angstrom)
            self._molecule.coords = coords
            sel = self._molecule.atomselect(self._config, indexes=True)

            com = [0, 0, 0] * angstroms
            mass = 0.
            for i in sel:
                xx = coords[i, :] * angstroms
                com = com + xx
                mass = mass + 1
            com = com / mass

            r = np.array(self._config.flatbottom_range) * angstroms
            #print(" Applying flatbottom constraints to " + str(len(sel)) + " atoms using: ")
            #print("  Reference    : '" + self._config.flatbottom_reference + "'")
            #print("  Bounding box : [%.2f %.2f %.2f]" % (
            #    (com[0].value_in_unit(angstrom)),
            #    (com[1].value_in_unit(angstrom)),
            #    (com[2].value_in_unit(angstrom)))
            #      )
        d = r.value_in_unit(angstrom)
        d[0] = d[0] + self._config.flatbottom_offset[0]
        d[1] = d[1] + self._config.flatbottom_offset[1]
        d[2] = d[2] + self._config.flatbottom_offset[2]

        self._flatbottom_force.setGlobalParameterDefaultValue(1, d[0] * .1)
        self._flatbottom_force.setGlobalParameterDefaultValue(2, d[1] * .1)
        self._flatbottom_force.setGlobalParameterDefaultValue(3, d[2] * .1)
        update = True

        if update:
            self._flatbottom_force.updateParametersInContext(simulation.context)
