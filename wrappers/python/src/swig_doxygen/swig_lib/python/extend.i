%extend OpenMM::Context {
  PyObject *_getStateAsLists(int getPositions,
                            int getVelocities,
                            int getForces,
                            int getEnergy,
                            int getParameters,
                            int enforcePeriodic,
                            int groups) {
    State state;
    Py_BEGIN_ALLOW_THREADS
    int types = 0;
    if (getPositions) types |= State::Positions;
    if (getVelocities) types |= State::Velocities;
    if (getForces) types |= State::Forces;
    if (getEnergy) types |= State::Energy;
    if (getParameters) types |= State::Parameters;
    state = self->getState(types, enforcePeriodic, groups);
    Py_END_ALLOW_THREADS
    return _convertStateToLists(state);
  }


  %pythoncode {
    def getState(self,
                 getPositions=False,
                 getVelocities=False,
                 getForces=False,
                 getEnergy=False,
                 getParameters=False,
                 enforcePeriodicBox=False,
                 groups=-1):
        """
        getState(self,
                 getPositions = False,
                 getVelocities = False,
                 getForces = False,
                 getEnergy = False,
                 getParameters = False,
                 enforcePeriodicBox = False,
                 groups = -1)
              -> State
        
        Get a State object recording the current state information stored in this context.
        
        Parameters:
         - getPositions (bool=False) whether to store particle positions in the State
         - getVelocities (bool=False) whether to store particle velocities in the State
         - getForces (bool=False) whether to store the forces acting on particles in the State
         - getEnergy (bool=False) whether to store potential and kinetic energy in the State
         - getParameter (bool=False) whether to store context parameters in the State
         - enforcePeriodicBox (bool=False) if false, the position of each particle will be whatever position is stored in the Context, regardless of periodic boundary conditions.  If true, particle positions will be translated so the center of every molecule lies in the same periodic box.
         - groups (int=-1) a set of bit flags for which force groups to include when computing forces and energies.  Group i will be included if (groups&(1<<i)) != 0.  The default value includes all groups.
        """
        
        if getPositions: getP=1
        else: getP=0
        if getVelocities: getV=1
        else: getV=0
        if getForces: getF=1
        else: getF=0
        if getEnergy: getE=1
        else: getE=0
        if getParameters: getPa=1
        else: getPa=0
        if enforcePeriodicBox: enforcePeriodic=1
        else: enforcePeriodic=0

        (simTime, periodicBoxVectorsList, energy, coordList, velList,
         forceList, paramMap) = \
            self._getStateAsLists(getP, getV, getF, getE, getPa, enforcePeriodic, groups)
        
        state = State(simTime=simTime,
                      energy=energy,
                      coordList=coordList,
                      velList=velList,
                      forceList=forceList,
                      periodicBoxVectorsList=periodicBoxVectorsList,
                      paramMap=paramMap)
        return state
  }

  %feature("docstring") createCheckpoint1 "Create a checkpoint recording the current state of the Context.
This should be treated as an opaque block of binary data.  See loadCheckpoint() for more details.

Returns: a string containing the checkpoint data
"
  std::string createCheckpoint1() {
    std::stringstream stream(std::ios_base::in | std::ios_base::out | std::ios_base::binary);
    self->createCheckpoint(stream);
    return stream.str();
  }

  %feature ("docstring") loadCheckpoint "Load a checkpoint that was written by createCheckpoint().

A checkpoint contains not only publicly visible data such as the particle positions and
velocities, but also internal data such as the states of random number generators.  Ideally,
loading a checkpoint should restore the Context to an identical state to when it was written,
such that continuing the simulation will produce an identical trajectory.  This is not strictly
guaranteed to be true, however, and should not be relied on.  For most purposes, however, the
internal state should be close enough to be reasonably considered equivalent.

A checkpoint contains data that is highly specific to the Context from which it was created.
It depends on the details of the System, the Platform being used, and the hardware and software
of the computer it was created on.  If you try to load it on a computer with different hardware,
or for a System that is different in any way, loading is likely to fail.  Checkpoints created
with different versions of OpenMM are also often incompatible.  If a checkpoint cannot be loaded,
that is signaled by throwing an exception.

Parameters:
 - checkpoint (string) the checkpoint data to load
"
  void loadCheckpoint(std::string checkpoint) {
    std::stringstream stream(std::ios_base::in | std::ios_base::out | std::ios_base::binary);
    stream << checkpoint;
    self->loadCheckpoint(stream);
  }
 }

%extend OpenMM::RPMDIntegrator {
  PyObject *_getStateAsLists(int copy,
                            int getPositions,
                            int getVelocities,
                            int getForces,
                            int getEnergy,
                            int getParameters,
                            int enforcePeriodic,
                            int groups) {
    State state;
    Py_BEGIN_ALLOW_THREADS
    int types = 0;
    if (getPositions) types |= State::Positions;
    if (getVelocities) types |= State::Velocities;
    if (getForces) types |= State::Forces;
    if (getEnergy) types |= State::Energy;
    if (getParameters) types |= State::Parameters;
    state = self->getState(copy, types, enforcePeriodic, groups);
    Py_END_ALLOW_THREADS
    return _convertStateToLists(state);
  }


  %pythoncode {
    def getState(self,
                 copy,
                 getPositions=False,
                 getVelocities=False,
                 getForces=False,
                 getEnergy=False,
                 getParameters=False,
                 enforcePeriodicBox=False,
                 groups=-1):
        """
        getState(self,
                 copy,
                 getPositions = False,
                 getVelocities = False,
                 getForces = False,
                 getEnergy = False,
                 getParameters = False,
                 enforcePeriodicBox = False,
                 groups = -1)
              -> State
        
        Get a State object recording the current state information about one copy of the system.
        
        Parameters:
         - copy (int) the index of the copy for which to retrieve state information
         - getPositions (bool=False) whether to store particle positions in the State
         - getVelocities (bool=False) whether to store particle velocities in the State
         - getForces (bool=False) whether to store the forces acting on particles in the State
         - getEnergy (bool=False) whether to store potential and kinetic energy in the State
         - getParameter (bool=False) whether to store context parameters in the State
         - enforcePeriodicBox (bool=False) if false, the position of each particle will be whatever position is stored in the Context, regardless of periodic boundary conditions.  If true, particle positions will be translated so the center of every molecule lies in the same periodic box.
         - groups (int=-1) a set of bit flags for which force groups to include when computing forces and energies.  Group i will be included if (groups&(1<<i)) != 0.  The default value includes all groups.
        """
        
        if getPositions: getP=1
        else: getP=0
        if getVelocities: getV=1
        else: getV=0
        if getForces: getF=1
        else: getF=0
        if getEnergy: getE=1
        else: getE=0
        if getParameters: getPa=1
        else: getPa=0
        if enforcePeriodicBox: enforcePeriodic=1
        else: enforcePeriodic=0

        (simTime, periodicBoxVectorsList, energy, coordList, velList,
         forceList, paramMap) = \
            self._getStateAsLists(copy, getP, getV, getF, getE, getPa, enforcePeriodic, groups)
        
        state = State(simTime=simTime,
                      energy=energy,
                      coordList=coordList,
                      velList=velList,
                      forceList=forceList,
                      periodicBoxVectorsList=periodicBoxVectorsList,
                      paramMap=paramMap)
        return state
  }
}


%extend OpenMM::NonbondedForce {
  %pythoncode {
    def addParticle_usingRVdw(self, charge, rVDW, epsilon):
        """Add particle using elemetrary charge.  Rvdw and epsilon,
           which is consistent with AMBER parameter file usage.
           Note that the sum of the radii of the two interacting atoms is
           the minimum energy point in the Lennard Jones potential and
           is often called rMin.  The conversion from sigma follows:
           rVDW = 2^1/6 * sigma/2
        """
        return self.addParticle(charge, rVDW/RVDW_PER_SIGMA, epsilon)

    def addException_usingRMin(self, particle1, particle2,
                               chargeProd, rMin, epsilon):
        """Add interaction exception using the product of the two atoms'
           elementary charges, rMin and epsilon, which is standard for AMBER
           force fields.  Note that rMin is the minimum energy point in the
           Lennard Jones potential.  The conversion from sigma is:
           rMin = 2^1/6 * sigma.
        """
        return self.addException(particle1, particle2,
                                 chargeProd, rMin/RMIN_PER_SIGMA, epsilon)
  }
}


%extend OpenMM::System {
  %pythoncode {
    def __getstate__(self):
        serializationString = XmlSerializer.serializeSystem(self)
        return serializationString

    def __setstate__(self, serializationString):
        system = XmlSerializer.deserializeSystem(serializationString)
        self.this = system.this
    def getForces(self):
        """Get the list of Forces in this System"""
        return [self.getForce(i) for i in range(self.getNumForces())]
  }
}

%extend OpenMM::CustomIntegrator {
    PyObject* getPerDofVariable(int index) const {
        std::vector<Vec3> values;
        self->getPerDofVariable(index, values);
        return copyVVec3ToList(values);
    }
}

%extend OpenMM::Force {
  %pythoncode {
    def __copy__(self):
        copy = self.__class__.__new__(self.__class__)
        copy.__init__(self)
        return copy

    def __deepcopy__(self, memo):
        return self.__copy__()
  }
}

%extend OpenMM::NEBIntegrator {
  PyObject *_getStateAsLists(int copy,
                            int getPositions,
                            int getVelocities,
                            int getForces,
                            int getEnergy,
                            int getParameters,
                            int enforcePeriodic,
                            int groups) {
    State state;
    Py_BEGIN_ALLOW_THREADS
    int types = 0;
    if (getPositions) types |= State::Positions;
    if (getVelocities) types |= State::Velocities;
    if (getForces) types |= State::Forces;
    if (getEnergy) types |= State::Energy;
    if (getParameters) types |= State::Parameters;
    state = self->getState(copy, types, enforcePeriodic, groups);
    Py_END_ALLOW_THREADS
    return _convertStateToLists(state);
  }

  %pythoncode {
    def getState(self,
                 copy,
                 getPositions=False,
                 getVelocities=False,
                 getForces=False,
                 getEnergy=False,
                 getParameters=False,
                 enforcePeriodicBox=False,
                 groups=-1):
        """
        """
        if getPositions: getP=1
        else: getP=0
        if getVelocities: getV=1
        else: getV=0
        if getForces: getF=1
        else: getF=0
        if getEnergy: getE=1
        else: getE=0
        if getParameters: getPa=1
        else: getPa=0
        if enforcePeriodicBox: enforcePeriodic=1
        else: enforcePeriodic=0

        (simTime, periodicBoxVectorsList, energy, coordList, velList,
         forceList, paramMap) = \
            self._getStateAsLists(copy, getP, getV, getF, getE, getPa, enforcePeriodic, groups)
        
        state = State(simTime=simTime,
                      energy=energy,
                      coordList=coordList,
                      velList=velList,
                      forceList=forceList,
                      periodicBoxVectorsList=periodicBoxVectorsList,
                      paramMap=paramMap)
        return state
  }
}
