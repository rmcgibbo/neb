#include "openmm/NEBIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/NEBKernels.h"
#include <ctime>
#include <string>
#include <iostream>

using namespace OpenMM;
using namespace std;


NEBIntegrator::NEBIntegrator(int numImages, double temperature, double frictionCoeff, double stepSize) {
    std::cout << "Constructor!";
}

void NEBIntegrator::initialize(ContextImpl& contextRef) {
    std::cout << "Initialize!";
}

void NEBIntegrator::stateChanged(State::DataType changed) {
    forcesAreValid = false;
}

vector<string> NEBIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateNEBStepKernel::Name());
    return names;
}

void NEBIntegrator::setPositions(int image, const vector<Vec3>& positions) {
    dynamic_cast<IntegrateNEBStepKernel&>(kernel.getImpl()).setPositions(image, positions);
    hasSetPosition = true;
}

void NEBIntegrator::setVelocities(int image, const vector<Vec3>& velocities) {
    dynamic_cast<IntegrateNEBStepKernel&>(kernel.getImpl()).setVelocities(image, velocities);
    hasSetVelocity = true;
}

State NEBIntegrator::getState(int image, int types, bool enforcePeriodicBox, int groups) {
    dynamic_cast<IntegrateNEBStepKernel&>(kernel.getImpl()).copyToContext(image, *context);
    return context->getOwner().getState(types, enforcePeriodicBox, groups);
}

void NEBIntegrator::step(int steps) {
    // if (!hasSetPosition) {
    //     // Initialize the positions from the context.
    //     
    //     State s = context->getOwner().getState(State::Positions);
    //     for (int i = 0; i < numImages; i++) {
    //         // i don't tho
    //         setPositions(i, s.getPositions());
    //     }
    // }
    // if (!hasSetVelocity) {
    //     // Initialize the velocities from the context.
    //     
    //     State s = context->getOwner().getState(State::Velocities);
    //     for (int i = 0; i < numCopies; i++)
    //         setVelocities(i, s.getVelocities());
    // }
    for (int i = 0; i < steps; ++i) {
        dynamic_cast<IntegrateNEBStepKernel&>(kernel.getImpl()).execute(*context, *this, forcesAreValid);
        forcesAreValid = true;
    }
}