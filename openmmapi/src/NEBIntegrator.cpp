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