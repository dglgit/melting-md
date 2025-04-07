from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout,argv
import numpy as np
import argparse
import json

parser=argparse.ArgumentParser()
parser.add_argument('-i','--inpcrd')
parser.add_argument('-p','--prmtop')
parser.add_argument('-t','--perTemp',default=50000)
parser.add_argument('-o')
parser.add_argument('-d','--dcdinterval',default=1000)
parser.add_argument('-b','--begin',default=300)
parser.add_argument('-e','--end')
parser.add_argument('-s','--step',default=1)
args=parser.parse_args()
inpcrd = AmberInpcrdFile(args.inpcrd)
prmtop = AmberPrmtopFile(args.prmtop, periodicBoxVectors=inpcrd.boxVectors)
method=PME if inpcrd.boxVectors is not None else NoCutoff
system = prmtop.createSystem(implicitSolvent=OBC2,nonbondedMethod=method, constraints=HBonds, implicitSolventSaltConc=0.15*moles/liter)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
with open(argv[3]+"_initial.pdb", 'w') as outfile:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=outfile, keepIds=True)



print('starting minimization')
simulation.minimizeEnergy()
#simulation.reporters.append(PDBReporter(argv[1].split('.')[0]+'.csv', 1000))
#simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
#        potentialEnergy=True, temperature=True))
simulation.context.setVelocitiesToTemperature(300)
print('starting equilibration')
simulation.step(10000)
with open(argv[3]+'_eqed.pdb', 'w') as outfile:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=outfile, keepIds=True)
output_base=argv[3].split('.')[0]

print('starting production run')

start=args.begin
stop=args.end
step=args.step
temps=np.arange(start,stop,step)

#perTemp=total_steps//len(temps)
perTemp=args.pertemp
print(perTemp)
#simulation.reporters.append(DCDReporter(output_base+"_traj.dcd", perTemp))
dcdInterval=args.dcdinterval
dcdf=open(output_base+'_traj.dcd','wb')
dcd=DCDFile(dcdf,simulation.topology,dcdInterval)
with open('outputBase'+"_meta.json",'w') as f:
    json.dump(args,f)
for temp in temps:
    print('setting temp at',temp)
    integrator.setTemperature(temp*kelvin)
    print(f"now at {temp}K")
    for i in range(0,perTemp//dcdInterval):
        integrator.step(dcdInterval)
        dcd.writeModel(simulation.context.getState(getPositions=True).getPositions())


