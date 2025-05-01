from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout,argv
import numpy as np
import argparse
import json
import numpy as np
import tqdm
import gc
#{'inpcrd': '', 'prmtop': '', 'output_root': 'out', 'temperatures': [], 'start': None, 'end': None, 'step': None, 'step_time_fs': 2, 'eq_steps': 100000, 'eq_temperature':300,'steps_per_temp': 10000,'dcd_interval':1000}

if len(argv)<2:
    print('error: missing config file (json). exiting')
    exit()
fname=argv[1]
f=open(fname,'r')
config=json.load(f)
f.close()
print(config)
output_base=config['output_root']
steps_per_temp=config['steps_per_temp']
step_size=config['step_time_fs']
eq_steps=config['eq_steps']
step_size=config['step_time_fs']
dcd_interval=config['dcd_interval']

print(f'equilibration for {eq_steps}, or {eq_steps*step_size*1e-6}ns')
print(f'production run steps per temperature: {steps_per_temp} ({steps_per_temp*step_size*1e-6} ns per temp)')
if config['start'] is not None:
    temps=np.arange(config['start'],config['end'],config['step'])
else:
    temps=config['temperatures']
print(f'running temperatures: {temps}')
print(f'capturing dcd every {dcd_interval} steps, or {dcd_interval*step_size*1e-6}ns')
inpcrd = AmberInpcrdFile(config['inpcrd'])
prmtop = AmberPrmtopFile(config['prmtop'], periodicBoxVectors=inpcrd.boxVectors)
#platform=openmm.Platform.getPlatformByName('OpenCL')
#platform=openmm.Platform.getPlatformByName('CPU')
method=PME if inpcrd.boxVectors is not None else NoCutoff
system = prmtop.createSystem(implicitSolvent=OBC2,nonbondedMethod=method, constraints=HBonds, implicitSolventSaltConc=0.15*moles/liter)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 1e-3*config['step_time_fs']*picoseconds)
# DO NOT SET PLATFORM ON M2 MACBOOK
# OPENCL IS SLOW AND SETTING AS CPU MAKES IT HANG
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
with open(config['output_root']+"_initial.pdb", 'w') as outfile:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=outfile, keepIds=True)

print('starting minimization')
simulation.minimizeEnergy()
#simulation.reporters.append(PDBReporter(argv[1].split('.')[0]+'.csv', 1000))
#simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
#        potentialEnergy=True, temperature=True))
simulation.context.setVelocitiesToTemperature(config['eq_temperature'])
print('starting equilibration for',eq_steps,'steps')
null_out = open(os.devnull, 'w')

# Set up StateDataReporter to write to /dev/null
simulation.reporters.append(StateDataReporter(null_out, 10000, step=True, temperature=True))
for i in tqdm.tqdm(range(0,eq_steps//1000)):
    simulation.step(1000)
    #gc.collect()
    #if i%10000==0:
    #    state=simulation.context.getState()
    #    del state
#simulation.reporters.pop()
with open(config['output_root']+'_eqed.pdb', 'w') as outfile:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=outfile, keepIds=True)

print('starting production run')
dcdf=open(output_base+'_traj.dcd','wb')
dcd=DCDFile(dcdf,simulation.topology,dcd_interval)
try:
    for temp in temps:
        print('setting temp at',temp)
        integrator.setTemperature(temp*kelvin)
        print(f"now at {temp}K")
        for i in range(0,steps_per_temp//dcd_interval):
            integrator.step(dcd_interval)
            dcd.writeModel(simulation.context.getState(getPositions=True).getPositions())
except Exception as e:
    dcdf.close()
    print('exception happened')
    raise e
