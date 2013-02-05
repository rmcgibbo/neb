from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np

numParticles = 3
numImages = 5
mass = 1.0
numSteps = 7500
stepSize = 0.001
springConstant = 100

def linear_interpolate(initial, final, n_images):
    """Linear interpolation between initial and final
    
    Example
    -------
    >>> linear_interpolate([0,0], [1,1], 5)
    [[ 0.    0.  ]
     [ 0.25  0.25]
     [ 0.5   0.5 ]
     [ 0.75  0.75]
     [ 1.    1.  ]]
    """
    delta = np.array(final) - np.array(initial)
    images = np.ones((n_images,) + delta.shape, dtype=np.float)
    # goes from 0 to 1
    lambdas = np.arange(n_images, dtype=np.float) / (n_images-1)
    for i in range(n_images):
        images[i] = initial + lambdas[i]*delta
        
    return images


system = System()
for i in range(numParticles):
    system.addParticle(mass)

length = 1.0 * nanometers
bond_k = 1.0 * kilojoules_per_mole / nanometer

forcefield = HarmonicBondForce()
forcefield.addBond(0, 1, length, bond_k)
forcefield.addBond(0, 2, length, bond_k)
forcefield.addBond(1, 2, length, bond_k)

system.addForce(forcefield)
        
integrator = NEBIntegrator(numImages, stepSize, springConstant)
context = Context(system, integrator)
initial = Quantity([[0,-0.5, 0], [0, 0, 0.5],  [0, 0.5, 0]], unit=nanometer)
final   = Quantity([[0,-0.5, 0], [0, 0, -0.5], [0, 0.5, 0]], unit=nanometer)


def minimize(pos):
    context.setPositions(pos)    
    LocalEnergyMinimizer.minimize(context, 1e-5)
    pos = context.getState(getPositions=True).getPositions(asNumpy=True)
    print context.getState(getForces=True).getForces()
    return pos

initial = minimize(initial)
final = minimize(final)

print 'initial', initial
print 'final', final

path = linear_interpolate(initial.value_in_unit(nanometers),
                         final.value_in_unit(nanometers), numImages)

import matplotlib.pyplot as pp
pp.title('Initial Path')
pp.xlim(-1,1)
pp.ylim(-1,1)
colors = ['r', 'g', 'b', 'y', 'm', 'c']
for i, image in enumerate(path):
    pp.scatter(image[:,1], image[:,2], label=str(i), color=colors[i])
    pp.plot(image[:,1], image[:,2], color=colors[i])
pp.legend()
#pp.savefig('initial.png')

for i, image in enumerate(path):
    integrator.setPositions(i, Quantity(image, unit=nanometer))

# set all of the positions to the same
for i in range(numImages):
    print i
    integrator.setPositions(i, Quantity([[0,0,0], [0,0.5,0], [0,1,0]], unit=nanometer))

# run the integrator    
integrator.step(numSteps)

print 'Halfway through'
for i in range(numImages):
    print 'Image', i
    print integrator.getState(i, getPositions=True).getPositions(asNumpy=True)

integrator.step(numSteps)


pp.figure()
pp.title('Final Path')
pp.xlim(-1,1)
pp.ylim(-1,1)
for i in range(numImages):
    print 'Image', i
    state = integrator.getState(i, getPositions=True, getForces=True)
    image = state.getPositions(asNumpy=True)
    print image
    #print state.getForces(asNumpy=True)
    
    pp.scatter(image[:,1], image[:,2], label=str(i), color=colors[i])
    pp.plot(image[:,1], image[:,2], color=colors[i])
pp.legend()
#pp.savefig('final.png')

# at the end
pp.show()
