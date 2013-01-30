#ifndef NEB_KERNELS_H_
#define NEB_KERNELS_H_


#include "openmm/NEBIntegrator.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/Vec3.h"
#include <string>
#include <vector>

namespace OpenMM {

class IntegrateNEBStepKernel : public KernelImpl {
public:
    static std::string Name() {
        return "IntegrateNEBStep";
    }
    IntegrateNEBStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the NEBIntegrator this kernel will be used for
     */
    virtual void initialize(const System& system, const NEBIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the NEBIntegrator this kernel is being used for
     * @param forcesAreValid if the context has been modified since the last time step, this will be
     *                       false to show that cached forces are invalid and must be recalculated
     */
    virtual void execute(ContextImpl& context, const NEBIntegrator& integrator, bool forcesAreValid) = 0;
    /**
     * Get the positions of all particles in one copy of the system.
     */
    virtual void setPositions(int copy, const std::vector<Vec3>& positions) = 0;
    /**
     * Get the velocities of all particles in one copy of the system.
     */
    virtual void setVelocities(int copy, const std::vector<Vec3>& velocities) = 0;
    /**
     * Copy positions and velocities for one copy into the context.
     */
    virtual void copyToContext(int copy, ContextImpl& context) = 0;
};

} // namespace OpenMM

#endif /*NEB_KERNELS_H_*/
