#ifndef OPENMM_REFERENCERNEBKERNELFACTORY_H_
#define OPENMM_REFERENCENEBKERNELFACTORY_H_

#include "openmm/KernelFactory.h"

namespace OpenMM {

class ReferenceNEBKernelFactory : public KernelFactory {
public:
    KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const;
};

} // namespace OpenMM

#endif
