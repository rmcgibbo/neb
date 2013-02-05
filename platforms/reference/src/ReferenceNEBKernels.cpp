#include "ReferenceNEBKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

// 
static void OptGradientDescent(std::vector<RealVec>& positions, std::vector<RealVec>& forces, const RealOpenMM stepSize) {
    int numParticles = positions.size();
    for (int i = 0; i < numParticles; i++)
      positions[i] += forces[i] * stepSize;
}

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->velocities);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

void ReferenceIntegrateNEBStepKernel::initialize(const System& system, const NEBIntegrator& integrator) {
    std::cout << "Initialize the ReferenceIntegrateNEBStepKernel. Release the hounds!\n";
    int numImages = integrator.getNumImages();
    int numParticles = system.getNumParticles();
    positions.resize(numImages);
    velocities.resize(numImages);
    forces.resize(numImages);

    for (int i = 0; i < numImages; i++) {
	positions[i].resize(numParticles);
	velocities[i].resize(numParticles);
	forces[i].resize(numParticles);
    }	 
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateNEBStepKernel::execute(ContextImpl& context, const NEBIntegrator& integrator, bool forcesAreValid) {
    const int numImages = positions.size();
    const RealOpenMM stepSize = integrator.getStepSize();
    const int numParticles = positions[0].size();
    const RealOpenMM k = integrator.getSpringConstant();
    const System& system = context.getSystem();
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& vel = extractVelocities(context);
    vector<RealVec>& f = extractForces(context);
    
    //compute the forces on each of the images
    if (!forcesAreValid) {
	for (int i = 0; i < numImages; i++) {
	    pos = positions[i];
	    vel = velocities[i];
	    context.calcForcesAndEnergy(true, false);
	    forces[i] = f;
	}
    }
    
    //std::cout << "ReferenceIntegrateNEBStepKernel execute" << endl;
    // compute the spring forces on each of the 
    // the force on the endpoints is going to be zeroed out.
    
    // declare a 3d tensor of size numImages x numPartices x 3 that
    // will hold the force due to the spring between adjacent images
    vector< vector< Vec3 > > springForces(numImages);
    for (int i = 1; i < numImages - 1; i++) {
	// Compute the tangent vector from image[i-1] to image[i+1]
	vector<Vec3> unit_tangent(numParticles);
	RealOpenMM sq_norm_of_unit_tangent = 0; // square of norm
	for (int j = 0; j < numParticles; j++) {
	    for (int k = 0; k < 3; k++) {
		unit_tangent[j][k] = positions[i+1][j][k] - positions[i-1][j][k];
		sq_norm_of_unit_tangent += unit_tangent[j][k] * unit_tangent[j][k];
	    }
	}
	// normalize the tangent vector to length 1
	RealOpenMM norm_of_unit_tangent = sqrt(sq_norm_of_unit_tangent);

	for (int j = 0; j < numParticles; j++)
	    for (int k = 0; k < 3; k++)
	    unit_tangent[j][k] /= norm_of_unit_tangent;

	/*std::cout << "Unit Tangent" << endl;
	for (int j = 0; j < numParticles; j++)
	std::cout << unit_tangent[j] << endl;*/
	
	vector<Vec3> spring_force(numParticles);
	for (int j = 0; j < numParticles; j++) {
	    RealOpenMM f_x = -2*positions[i][j][0] + positions[i-1][j][0] + positions[i+1][j][0];
	    RealOpenMM f_y = -2*positions[i][j][1] + positions[i-1][j][1] + positions[i+1][j][1];
	    RealOpenMM f_z = -2*positions[i][j][2] + positions[i-1][j][2] + positions[i+1][j][2];
	    spring_force[j] = Vec3(numImages * k * f_x,
			      numImages * k * f_y,
			      numImages * k * f_z);
	}
	// compute the dot product between the forces on this image and the unit
	// tangent vector -- this is the length of the projection
	RealOpenMM dot_product_with_spring_force = 0;
	for (int j = 0; j < numParticles; j++)
	    for (int k = 0; k < 3; k++)
		dot_product_with_spring_force += spring_force[j][k]*unit_tangent[j][k];
	
	// set the spring_force to be just the parallel projection along the
	// tangent
	//std::cout << "dot product with spring force" << dot_product_with_spring_force;
	for (int j = 0; j < numParticles; j++)
	    for (int k = 0; k < 3; k++)
		spring_force[j][k] = unit_tangent[j][k] * dot_product_with_spring_force;
	
	/*std::cout << "Spring Force" << endl;
	for (int j = 0; j < numParticles; j++)
	std::cout << spring_force[j] << endl;*/

	// take only the perpindicular component of the force due to the potential
	// energy surface
	RealOpenMM dot_product_with_potential_force = 0;
	for (int j = 0; j < numParticles; j++)
	    for (int k = 0; k < 3; k++)
		dot_product_with_potential_force += forces[i][j][k]*unit_tangent[j][k];

	/*std::cout << "Raw forces" << endl;
        for (int j = 0; j < numParticles; j++)
	std::cout << forces[i][j] << endl;*/

	//for this force, we're subtracting off the parallel component
	for (int j = 0; j < numParticles; j++) {
	    for (int k = 0; k < 3; k++) {
	      forces[i][j][k] = forces[i][j][k] - unit_tangent[j][k] * dot_product_with_potential_force;
	      
	      // add the two forces together
	      //forces[i][j][k] += spring_force[j][k];
	    }
	}

	/*std::cout << "All forces" << endl;
	for (int j = 0; j < numParticles; j++)
	std::cout << forces[i][j] << endl;*/
    }
    
    // null out forces on first and last images
    for (int j = 0; j < numParticles; j++) {
	for (int k = 0; k < 3; k++) {
	    forces[0][j][k] = 0;
	    forces[numImages-1][j][k] = 0;
	}
    }

    /*for (int i = 0; i < numImages; i++)
      for (int j = 0; j < numParticles; j++)
      std::cout << "final_forces[" << j << "] " << forces[i][j] << endl;*/



    // run the actual integrator with the modified forces
    // for each system (skip the end images)
    // std::cout << "Force Norm" << endl;
    for (int i = 1; i < numImages - 1; i++) {
      double force_norm = 0.0;
      for (int j = 0; j < numParticles; j++) {
	for (int k = 0; k < 3; k++)
	  force_norm += forces[i][j][k] * forces[i][j][k];
	
      }
      //std::cout << force_norm << ", ";

      // by modifying pos, vel, and f, we change which coordinates
      // the context object is pointing to
      OptGradientDescent(positions[i], forces[i], stepSize);
    }
    //std::cout << endl;
    

}

void ReferenceIntegrateNEBStepKernel::setPositions(int image, const std::vector<Vec3>& pos) {
    int numParticles = positions[image].size();
    for (int i = 0; i < numParticles; i++)
	positions[image][i] = pos[i];
}

void ReferenceIntegrateNEBStepKernel::setVelocities(int image, const std::vector<Vec3>& vel) {
    int numParticles = velocities[image].size();
    for (int i = 0; i < numParticles; i++)
	velocities[image][i] = vel[i];
}

void ReferenceIntegrateNEBStepKernel::copyToContext(int image, ContextImpl& context) {
    //std::cout << "copyToContext ReferenceIntegrateNEBStepKernel\n";
    extractPositions(context) = positions[image];
    extractVelocities(context) = velocities[image];    

    /*std::cout << "Forces:" << endl;    
    for (int j = 0; j < forces[image].size(); j++) {
      std::cout << forces[image][j] << std::endl;
      }*/
}

