#ifndef REFERENCE_NEB_KERNELS_H_
#define REFERENCE_NEB_KERNELS_H_

#include "ReferencePlatform.h"
#include "ReferenceKernels.h"
#include "openmm/kernels.h"
#include "openmm/NEBKernels.h"
#include "SimTKUtilities/RealVec.h"

namespace OpenMM {

/**
 * This kernel is invoked by NEBIntegrator to take one time step, and to get and
 * set the state of system copies.
 */
class ReferenceIntegrateNEBStepKernel : public IntegrateNEBStepKernel {
public:
    ReferenceIntegrateNEBStepKernel(std::string name, const Platform& platform, ReferencePlatform::PlatformData& data) :
            IntegrateNEBStepKernel(name, platform), data(data) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the NEBIntegrator this kernel will be used for
     */
    void initialize(const System& system, const NEBIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the NEBIntegrator this kernel is being used for
     * @param forcesAreValid if the context has been modified since the last time step, this will be
     *                       false to show that cached forces are invalid and must be recalculated
     */
    void execute(ContextImpl& context, const NEBIntegrator& integrator, bool forcesAreValid);
    /**
     * Get the positions of all particles in one copy of the system.
     */
    void setPositions(int copy, const std::vector<Vec3>& positions);
    /**
     * Get the velocities of all particles in one copy of the system.
     */
    void setVelocities(int copy, const std::vector<Vec3>& velocities);
    /**
     * Copy positions and velocities for one copy into the context.
     */
    void copyToContext(int copy, ContextImpl& context);

private:
    ReferencePlatform::PlatformData& data;
    ReferenceIntegrateVerletStepKernel* dynamicsStepKernel; 
    VerletIntegrator* stepIntegrator;
    std::vector<std::vector<RealVec> > positions;
    std::vector<std::vector<RealVec> > velocities;
    std::vector<std::vector<RealVec> > forces;
};

} // namespace OpenMM

#endif