# neb
_OpenMM Nudged Elastic Band Plugin_

## installation

0) Clone this project into the plugins directory, `<openmm_source_dir>/plugins/`.

1) Add the following lines to the root `CMakeLists.txt` in `<OpenMMSrc>/CMakeLists.txt`.
This instructs the OpenMM CMake build system to build the NEB plugin.

```
# Nudged Elastic Band (NEB) Plugin
SET(OPENMM_BUILD_NEB_PLUGIN ON CACHE BOOL "Build NEB plugin")
SET(OPENMM_BUILD_NEB_PATH)
IF(OPENMM_BUILD_NEB_PLUGIN)
   SET(OPENMM_BUILD_NEB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/plugins/neb)
   ADD_SUBDIRECTORY(plugins/neb)
ENDIF(OPENMM_BUILD_NEB_PLUGIN)
```

2) Copy over the wrappers directory from this git repository into `<openmm_source_dir>`, replacing
the existing wrappers directory. This 

3) Run cmake, and then `make` and `make install`.

4) If the python layer doesn't build automatically, `cd` into `<open_build_directory>/python`
and run `python setup.py install`. You will probably have to set the `OPENMM_LIB_PATH` and
`OPENMM_INCLUDE_PATH`.


