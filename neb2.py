from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import matplotlib.pyplot as pp
np.random.seed(42)

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

    # square the lambdas!
    lambdas **= 2

    for i in range(n_images):
        images[i] = initial + lambdas[i]*delta
        
    return images

numParticles = 3
numImages = 4
mass = 1.0
numSteps = 500
stepSize = 0.001
springConstant = 10

system = System()
for i in range(numParticles):
    system.addParticle(mass)

# create the forcefield
length = 1.0 * nanometer
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

for i in range(numImages):
    integrator.setPositions(i, Quantity(path[i], unit=nanometers))

colors = ['r', 'b', 'g', 'y']
for r in range(2):
    pp.figure()
    pp.title('Round %s' % r)
    for i in range(numImages):
        state = integrator.getState(i, getPositions=True)
        image = state.getPositions(asNumpy=True)
        pp.plot(image[:,1], image[:,2], '-o', color=colors[i], label=str(i), markersize=30)

    pp.legend()
    integrator.step(numSteps)


    print 'On iteration', i
    print 'Positions of image 1'
    print integrator.getState(1, getPositions=True).getPositions(asNumpy=True)


pp.show()
