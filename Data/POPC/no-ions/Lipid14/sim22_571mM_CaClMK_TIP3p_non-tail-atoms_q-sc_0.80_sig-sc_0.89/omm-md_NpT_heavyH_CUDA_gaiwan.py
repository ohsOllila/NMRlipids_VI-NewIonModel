from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os
from datetime import datetime
#import openmmplumed
import IPython
#import omm_vfswitch
#import omm_readparams

seconds_per_day = 86400
NSTEPS = 250000 # 1ns with 4fs time step


dt = 0.004*picoseconds
temp = 313.15

print "Read coordinates and parameters"
gro = GromacsGroFile('conf_init.gro')
top = GromacsTopFile('topol.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds, hydrogenMass=4*atomic_mass_unit)

print "Created system"
print "Set-up barostat and integrator"

# barostat adopted from charmm-gui scripts
barostat = MonteCarloMembraneBarostat(1*bar, 0*bar*nanometer, temp*kelvin, MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree, 100)
system.addForce(barostat)

integrator = LangevinIntegrator(temp*kelvin, 1/picosecond, dt)

#%%
########## SIMULATION ##########

print "Set Simulation class"
#simulation_run = Simulation(top.topology, system, integrator, platform=openmm.Platform_getPlatformByName('CUDA'))   #, platformProperties={"CudaCompiler":"nvcc -isystem /usr/local/gcc/5.4.0/include/c++/5.4.0/ -ccbin /usr/local/gcc/5.4.0/bin/gcc"}) ##'OpenCL')) ##'CPU')) #
simulation_run = Simulation(top.topology, system, integrator, platform=openmm.Platform_getPlatformByName('CUDA'), platformProperties={"CudaCompiler":"nvcc -isystem /usr/local/gcc/5.4.0/include/c++/5.4.0/ -ccbin /usr/local/gcc/5.4.0/bin/gcc"}) ##'OpenCL')) ##'CPU')) #
simulation_run.context.setPositions(gro.positions)
simulation_run.loadState('state_checkpoint.openmm')

print "Set Reporter"
simulation_run.reporters.append(StateDataReporter("ener_openMM.log",  2500, step=True, potentialEnergy=True, temperature=True, speed=True))

# Minimize energy before MD?
#print "Minimize energy before MD"
#simulation_run.minimizeEnergy(maxIterations=2500)
print "Equilibrate before MD"
simulation_run.step(1.0*NSTEPS)

#%%

print "Set DCD Reporter" # on SCRATCH"
simulation_run.reporters.append(DCDReporter("traj_openMM.dcd", 5000))  # every 20ps
#simulation_run.reporters.append(PDBReporter('traj_openMM_every1ns.pdb', 250000))  # every 1ns

#IPython.embed()

#start = datetime.now()
for i in range(200):
	print('running part %d ...' % (i) )
	simulation_run.step(NSTEPS)
	#end = datetime.now()
	#time_elapsed = end -start
	# write checkpoint for restarting simulation
	simulation_run.saveState('state_checkpoint.openmm')

