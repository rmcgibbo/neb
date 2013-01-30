#include "ReferenceNEBKernelFactory.h"
#include "ReferenceNEBKernels.h"
#include "ReferencePlatform.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

#if defined(WIN32)
    #include <windows.h>
    extern "C" void initNEBReferenceKernels();
    BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
        if (ul_reason_for_call == DLL_PROCESS_ATTACH)
            initNEBReferenceKernels();
        return TRUE;
    }
#else
    extern "C" void __attribute__((constructor)) initNEBReferenceKernels();
#endif

extern "C" void initNEBReferenceKernels() {
    Platform& platform = Platform::getPlatformByName("Reference");
    ReferenceNEBKernelFactory* factory = new ReferenceNEBKernelFactory();
}

KernelImpl* ReferenceNEBKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    if (name == IntegrateNEBStepKernel::Name())
        return new ReferenceIntegrateNEBStepKernel(name, platform);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
