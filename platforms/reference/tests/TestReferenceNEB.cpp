#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/HarmonicBondForce.h"
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
    
    const int numParticles = 2;
    const int numImages = 5;
    const double temperature = 300.0;
    const double mass = 1.0;
    const int numSteps = 10;
    
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(mass);
    
    // they'll attract with coulomb
    NonbondedForce* forceField = new NonbondedForce();
    forceField->addParticle(0.5, 1, 0);
    forceField->addParticle(-1.5, 1, 0);
    system.addForce(forceField);
        
    NEBIntegrator integ(numImages, temperature, 10.0, 0.001, 1.0);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);

    // random number context/
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numImages; i++) {
        for (int j = 0; j < numParticles; j++)
            positions[j] = Vec3(0.02*genrand_real2(sfmt), 0.02*genrand_real2(sfmt), 0.02*genrand_real2(sfmt));
        integ.setPositions(i, positions);
    }
    
    std::cout << "Stepping...\n";
    integ.step(numSteps);
    std::cout << "Done!\n";
}
