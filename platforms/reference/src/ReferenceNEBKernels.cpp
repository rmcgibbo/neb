#include "ReferenceNEBKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

void ReferenceIntegrateNEBStepKernel::initialize(const System& system, const NEBIntegrator& integrator) {
    std::cout << "Initialize the ReferenceIntegrateNEBStepKernel. Release the hounds!";
}

void ReferenceIntegrateNEBStepKernel::execute(ContextImpl& context, const NEBIntegrator& integrator, bool forcesAreValid) {
        std::cout << "Execute the ReferenceIntegrateNEBStepKernel. Release the hounds!";
}

void ReferenceIntegrateNEBStepKernel::setPositions(int copy, const std::vector<Vec3>& positions) {
}

void ReferenceIntegrateNEBStepKernel::setVelocities(int copy, const std::vector<Vec3>& velocities) {
    
}
void ReferenceIntegrateNEBStepKernel::copyToContext(int copy, ContextImpl& context) {
    
}

