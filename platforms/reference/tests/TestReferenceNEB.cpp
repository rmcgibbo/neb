#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/NEBIntegrator.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

int main() {
    Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
    std::cout << "Hello World! (TestReferenceNEB)\n";
    
    const int numParticles = 3;
    const int numImages = 3;
    const double temperature = 300.0;
    const double mass = 1.0;
    
    System system;
    for (int i = 0; i < numParticles; i++)
        system.addParticle(mass);
    
    NEBIntegrator integ(numImages, temperature, 10.0, 0.001);
}