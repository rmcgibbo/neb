#include "openmm/NEBIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/NEBKernels.h"
#include <ctime>
#include <string>
#include <iostream>

using namespace OpenMM;
using std::string;
using std::vector;


NEBIntegrator::NEBIntegrator(int numImages, double temperature, double frictionCoeff, double stepSize, double springConstant) :
    owner(NULL), numImages(numImages), springConstant(springConstant), forcesAreValid(false), hasSetPosition(false), hasSetVelocity(false){
    std::cout << "Constructor! (NEBIntegrator)\n";
    setTemperature(temperature);
    setFriction(frictionCoeff);
    setStepSize(stepSize);
    setConstraintTolerance(1e-4);
    setRandomNumberSeed((int) time(NULL));
}

void NEBIntegrator::initialize(ContextImpl& contextRef) {
    std::cout << "Initialize! (NEBIntegrator)\n";

    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    if (contextRef.getSystem().getNumConstraints() > 0)
        throw OpenMMException("NEBIntegrator cannot be used with Systems that include constraints");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateNEBStepKernel::Name(), contextRef);
    dynamic_cast<IntegrateNEBStepKernel&>(kernel.getImpl()).initialize(contextRef.getSystem(), *this);
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
    std::cout << "NEBIntegrator::setPositions...\n";
    dynamic_cast<IntegrateNEBStepKernel&>(kernel.getImpl()).setPositions(image, positions);
    hasSetPosition = true;
    std::cout << "NEBIntegrator::setPositions... finished\n";
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
    if (!hasSetPosition) {
        // Initialize the positions from the context.
        throw "Error! No positions!";
    }
    if (!hasSetVelocity) {
        // Initialize the velocities from the context.
        
        State s = context->getOwner().getState(State::Velocities);
        for (int i = 0; i < numImages; i++)
            setVelocities(i, s.getVelocities());
    }
    std::cout << "NEBIntegrator::step\n";
    
    for (int i = 0; i < steps; ++i) {
        dynamic_cast<IntegrateNEBStepKernel&>(kernel.getImpl()).execute(*context, *this, forcesAreValid);
        forcesAreValid = true;
    }
}