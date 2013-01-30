#include "ReferenceNEBKernelFactory.h"
#include "ReferenceNEBKernels.h"
#include "ReferencePlatform.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

extern "C" void __attribute__((constructor)) initNEBReferenceKernels();

extern "C" void initNEBReferenceKernels() {
    Platform& platform = Platform::getPlatformByName("Reference");
    ReferenceRpmdKernelFactory* factory = new ReferenceNEBKernelFactory();
    platform.registerKernelFactory(IntegrateNEBStepKernel::Name(), factory);
}

KernelImpl* ReferenceNEBKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    if (name == IntegrateNEBStepKernel::Name())
        return new ReferenceIntegrateNEBStepKernel(name, platform);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
