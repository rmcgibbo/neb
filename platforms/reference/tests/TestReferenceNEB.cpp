#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/LocalEnergyMinimizer.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/NEBIntegrator.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include "openmm/NonbondedForce.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

int main() {
    Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
    std::cout << "Start of the test (TestReferenceNEB)\n";
    
    const int numParticles = 3;
    const int numImages = 5;
    const double temperature = 300.0;
    const double mass = 1.0;
    const int numSteps = 100;
    const double stepSize = 0.0001;
    const double springConstant = 1;
    
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(mass);
    
    // they'll attract with coulomb
    NonbondedForce* forceField = new NonbondedForce();

    // epsilon 1, sigma 1
    forceField->addParticle(0, 1, 1);
    forceField->addParticle(0, 1, 1);
    forceField->addParticle(0, 1, 1);

    system.addForce(forceField);
        
    NEBIntegrator integ(numImages, stepSize, springConstant);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
   

    // random number context/
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    vector<Vec3> velocities(numParticles);

    for (int i = 0; i < numImages; i++) {
      // two particles: first is at the origin
      // second goes from z=1 to z=2
      positions[0] = Vec3(0, 0, 0);
      positions[1] = Vec3(1, 0, 0);
      positions[2] = Vec3(0.5, 0.5 - 1.0 * i / (numImages - 1), 0);

      integ.setPositions(i, positions);
      integ.setVelocities(i, velocities);
    }

    //LocalEnergyMinimizer::minimize(context);
    //positions = context.getState(State::Positions).getPositions();

    std::cout << "Positions: " << endl;
    for (int j = 0; j < numParticles; j++) {
      std::cout << positions[j] << std::endl;
    }

    std::cout << "Stepping...\n";
    integ.step(numSteps);
    
    for (int i = 0; i < numImages; i++) {
      std::cout << "Image: " << i << endl;
      State state = integ.getState(i, State::Positions | State::Velocities | State::Forces);
      positions = state.getPositions();
      velocities = state.getVelocities();

      std::cout << "Positions: " << endl;
      for (int j = 0; j < numParticles; j++) {
	std::cout << positions[j] << std::endl;
      }
    }
}
