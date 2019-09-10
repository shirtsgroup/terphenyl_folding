import numpy as np
from simtk import unit
from simtk.openmm.openmm import LangevinIntegrator
from simtk.openmm.app.simulation import Simulation
from simtk.openmm.app.topology import Topology
from simtk.openmm.app.pdbfile import PDBFile
from simtk.openmm.app.pdbreporter import PDBReporter
from simtk.openmm.app.statedatareporter import StateDataReporter
from simtk.openmm.app.forcefield import ForceField


# Test of OpenMM simulation of terphenyl oligomers.

pdb_mm_obj = PDBFile('new.pdb')
positions = pdb_mm_obj.getPositions()
topology = pdb_mm_obj.getTopology()
temperature=300.0 * unit.kelvin
simulation_time_step=0.2 * unit.femtosecond
forcefield = ForceField('amber14-all.xml','/home/gmeek/Foldamers/terphenyl_folding/setup_files/tcm.xml')
system = forcefield.createSystem(topology, nonbondedMethod=2, nonbondedCutoff=1.0 * unit.nanometer)
integrator = LangevinIntegrator(temperature._value,friction,simulation_time_step.in_units_of(unit.picosecond)._value)
simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)
simulation.reporters.append(PDBReporter("test_openmm.pdb",1))
simulation.reporters.append(StateDataReporter("test_openmm.dat",1, \
        step=True, totalEnergy=True, potentialEnergy=True, kineticEnergy=True, temperature=True))
simulation.minimizeEnergy()
simulation.step(100)

exit()
