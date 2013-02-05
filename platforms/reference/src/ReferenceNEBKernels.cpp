#include "ReferenceNEBKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

/* 
 * Steepest descent optimization. Advance the positions along the forces by stepSize
 */ 
static void OptGradientDescent(std::vector<RealVec>& positions, std::vector<RealVec>& forces, const RealOpenMM stepSize) {
    int numParticles = positions.size();
    for (int i = 0; i < numParticles; i++)
      positions[i] += forces[i] * stepSize;
}

static std::vector<RealVec> RealSpaceTangent(int image, std::vector< std::vector<RealVec> >& positions, std::vector<RealOpenMM>& potentials) {
  // Tangent vector to the path at this image, translated from QE's fortran code
  int numImages = positions.size();
  int numParticles = positions[0].size();
  std::vector<RealVec> tangent(numParticles);
  
  //first image: take the forward difference
  if (image == 0) {
    for (int i = 0; i < numParticles; i++)
      tangent[i] = positions[image+1][i] - positions[image][i];
  
  // last image: take the backward difference
  } else if (image == numImages - 1) {
    for (int i = 0; i < numParticles; i++)
      tangent[i] = positions[image][i] - positions[image-1][i];
  } else {

    RealOpenMM V_previous = potentials[image - 1];
    RealOpenMM V_actual = potentials[image];
    RealOpenMM V_next = potentials[image + 1];

    if ((V_next > V_actual) && (V_actual > V_previous)) {
      //cout << "forward" << endl;
      //forward difference if energy is increasing
      for (int i = 0; i < numParticles; i++)
	tangent[i] = positions[image+1][i] - positions[image][i];
    } else if ((V_next < V_actual) && (V_actual < V_previous)) {
      //cout << "backward" << endl;
      //backwards difference if energy is decreasing
      for (int i = 0; i < numParticles; i++)
	tangent[i] = positions[image][i] - positions[image-1][i];
    } else {
      // otherwise it's more complicated...
      RealOpenMM abs_next = abs(V_next - V_actual);
      RealOpenMM abs_previous = abs(V_actual - V_previous);
      RealOpenMM delta_V_max = max(abs_next, abs_previous);
      RealOpenMM delta_V_min = min(abs_next, abs_previous);
      /*cout << "Mixed" << endl;
      cout << V_previous << " " << V_next << " " <<  V_actual << endl;
      cout << abs_next << " " << abs_previous << endl;
      cout << delta_V_max << " " << delta_V_min << endl;*/
      
      if (V_next > V_previous) {
	// sort of a weighted sum of the forward and backwards differences
	for (int i = 0; i < numParticles; i++)
	  tangent[i] = (positions[image+1][i] - positions[image][i]) * delta_V_max + \
	    (positions[image][i] - positions[image-1][i]) * delta_V_min;
      } else if (V_next < V_previous) {
	// reverse the weights
	for (int i = 0; i < numParticles; i++)
	  tangent[i] = (positions[image+1][i] - positions[image][i]) * delta_V_min + \
	    (positions[image][i] - positions[image-1][i]) * delta_V_max;
      } else {
	// totally symetric
	for (int i = 0; i < numParticles; i++)
	  tangent[i] = positions[image+1][i] - positions[image-1][i];
      }
    }
  }

  // now we need to normalize the tangent to unit length
  RealOpenMM sq_norm_of_tangent = 0; // square of norm                                                                                      
  for (int i = 0; i < numParticles; i++)
    for (int j = 0; j < 3; j++)
      sq_norm_of_tangent += tangent[i][j] * tangent[i][j];

  // normalize the tangent vector to length 1
  RealOpenMM norm_of_tangent = sqrt(sq_norm_of_tangent);
  
  for (int i = 0; i < numParticles; i++)
    for (int j = 0; j < 3; j++)
      tangent[i][j] /= norm_of_tangent;

  /*cout << "RealSpaceTangent(" << image << "):" << endl;
  for (int i = 0; i < numParticles; i++)
  cout << tangent[i] << endl;*/
  
  return tangent;
}

static void NEBForce(std::vector< std::vector <RealVec> >& positions,
		   std::vector< std::vector < RealVec> >& forces,
		   std::vector<RealOpenMM>& springConstants,
		   std::vector<RealOpenMM>& potentials) {
  // Calculte the NEB Force, modifying the forces array in place
  // translated from Q-E's neb_gradient()

  int numImages = positions.size();
  int numParticles = positions[0].size();
  RealOpenMM forwards_norm, backwards_norm, forwards_k_sum, backwards_k_sum;
  RealOpenMM forces_tangent_dot;
  std::vector< std::vector<RealVec> > elastic_grad(numImages);

  for (int i = 1; i < numImages - 1; i++) {
    backwards_norm = forwards_norm = 0;
    for (int j = 0; j < numParticles; j++) {
      for (int k = 0; k < 3; k++) {
	backwards_norm += (positions[i][j][k] - positions[i-1][j][k]) * (positions[i][j][k] - positions[i-1][j][k]);
	forwards_norm  += (positions[i+1][j][k] - positions[i][j][k]) * (positions[i+1][j][k] - positions[i][j][k]);
      }
    }

    //cout << "got norms" << endl;

    backwards_k_sum = springConstants[i] + springConstants[i-1];
    forwards_k_sum = springConstants[i] + springConstants[i+1];
    std::vector< RealVec> tangent = RealSpaceTangent(i, positions, potentials);

    //gradient of the spring potential acting along the tangent
    std::vector< RealVec > elastic_grad_i(numParticles);
    for (int j = 0; j < numParticles; j++)
      for (int k = 0; k < 3; k++)
	elastic_grad_i[j][k] = 0.5 * tangent[j][k] * (backwards_k_sum * backwards_norm - forwards_k_sum * forwards_norm);
    
    //cout << "got elastic grad" << endl;
    
    forces_tangent_dot = 0;
    for (int j = 0; j < numParticles; j++)
      for (int k = 0; k < 3; k++)
	forces_tangent_dot += forces[i][j][k] * tangent[j][k];

    for (int j = 0; j < numParticles; j++)
      for (int k = 0; k < 3; k++)
	forces[i][j][k] = (forces[i][j][k] - tangent[j][k] * forces_tangent_dot) - elastic_grad_i[j][k];

    //cout << "mod forces" << endl;
  }

  /*for (int i = 0; i < numImages; i++) {
    std::cout << "Forces " << i << ": " << endl;
    for (int j = 0; j < numParticles; j++) {
      std::cout << forces[i][j] << endl;
    }
    std::cout << endl;
    }*/

  return;
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
    //std::cout << "Initialize the ReferenceIntegrateNEBStepKernel. Release the hounds!!!!\n";
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

void ReferenceIntegrateNEBStepKernel::execute(ContextImpl& context, const NEBIntegrator& integrator) {
    const int numImages = positions.size();
    const RealOpenMM stepSize = integrator.getStepSize();
    const int numParticles = positions[0].size();
    const System& system = context.getSystem();
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& vel = extractVelocities(context);
    vector<RealVec>& f = extractForces(context);
    vector<RealOpenMM> potentials(numImages);

    //cout << "execute" << endl;
    
    // this is a temp. thing. really we should be using adapive spring forces
    // or just setting it in initialize() as a const
    vector<RealOpenMM> springConstants(numImages);
    for (int i = 0; i < numImages; i++)
      springConstants[i] = integrator.getSpringConstant();
      
    //compute the forces on each of the images
    for (int i = 0; i < numImages; i++) {
      pos = positions[i];
      vel = velocities[i];
      potentials[i] = context.calcForcesAndEnergy(true, true);
      forces[i] = f;
    }

    //modifies the forces in place
    NEBForce(positions, forces, springConstants, potentials);

    //cout << "execute2" << endl;

    //run the opt step
    for (int i = 0; i < numImages; i++) 
      OptGradientDescent(positions[i], forces[i], stepSize);
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

