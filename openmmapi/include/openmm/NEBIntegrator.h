#ifndef OPENMM_NEBINTEGRATOR_H_
#define OPENMM_NEBINTEGRATOR_H_


#include "openmm/Integrator.h"
#include "openmm/Kernel.h"
#include "openmm/State.h"
#include "openmm/Vec3.h"
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

class OPENMM_EXPORT NEBIntegrator : public Integrator {
public:
    /**
     * Create a NEBIntegrator.
     *
     * @param numImages      the number of copies of the system that should be simulated
     * @param stepSize       the step size with which to integrator the system (in picoseconds)
     * @param springConstant the spring constant between the images
     */
     NEBIntegrator(int numImages, double stepSize, double springConstant);
     
     int getNumImages() const {
         return numImages;
     }
     
     double getSpringConstant() const {
         return springConstant;
     }

     double getStepSize() const {
          return stepSize;
     }

     void setSpringConstant(double k) {
         springConstant = k;
     }
          
     int getRandomNumberSeed() const {
         return randomNumberSeed;
     }
     
     void setRandomNumberSeed(int seed) {
         randomNumberSeed = seed;
     }
     /**
      * Set the positions of all particles in one copy of the system.
      * 
      * @param copy      the index of the copy for which to set positions
      * @param positions the positions of all particles in the system
      */
     void setPositions(int copy, const std::vector<Vec3>& positions);

     /**
      * Get the velocities of all particles in one copy of the system.
      * 
      * @param copy       the index of the copy for which to set velocities
      * @param velocities the velocities of all particles in the system
      */
     void setVelocities(int copy, const std::vector<Vec3>& velocities);

     /**
      * Get a State object recording the current state information about one copy of the system.
      * 
      * @param copy  the index of the copy for which to retrieve state information
      * @param types the set of data types which should be stored in the State object.  This
      * should be a union of DataType values, e.g. (State::Positions | State::Velocities).
      * @param enforcePeriodicBox if false, the position of each particle will be whatever position
      * is stored by the integrator, regardless of periodic boundary conditions.  If true, particle
      * positions will be translated so the center of every molecule lies in the same periodic box.
      * @param groups a set of bit flags for which force groups to include when computing forces
      * and energies.  Group i will be included if (groups&(1<<i)) != 0.  The default value includes all groups.
      */
     State getState(int copy, int types, bool enforcePeriodicBox=false, int groups=0xFFFFFFFF);

     /**
      * Advance a simulation through time by taking a series of time steps.
      *
      * @param steps   the number of time steps to take
      */
     void step(int steps);

protected:
     /**
      * This will be called by the Context when it is created.  It informs the Integrator
      * of what context it will be integrating, and gives it a chance to do any necessary initialization.
      * It will also get called again if the application calls reinitialize() on the Context.
      */
     void initialize(ContextImpl& context);
     /**
      * When the user modifies the state, we need to mark that the forces need to be recalculated.
      */
     void stateChanged(State::DataType changed);
     /**
      * Get the names of all Kernels used by this Integrator.
      */
     std::vector<std::string> getKernelNames();

private:
    double temperature, friction, springConstant, stepSize;
    int numImages, randomNumberSeed;
    bool forcesAreValid, hasSetPosition, hasSetVelocity;
    ContextImpl* context;
    Context* owner;
    Kernel kernel;
};

} // namespace OpenMM


#endif  /*OPENMM_NEBINTEGRATOR_H_*/
