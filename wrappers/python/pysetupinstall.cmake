set(ENV{OPENMM_INCLUDE_PATH} "/home/rmcgibbo/OpenMM5.0-Source/./include;/home/rmcgibbo/OpenMM5.0-Source/./include/openmm;/home/rmcgibbo/OpenMM5.0-Source/./include/openmm/internal;/home/rmcgibbo/OpenMM5.0-Source/openmmapi/include;/home/rmcgibbo/OpenMM5.0-Source/openmmapi/include/openmm;/home/rmcgibbo/OpenMM5.0-Source/openmmapi/include/openmm/internal;/home/rmcgibbo/OpenMM5.0-Source/olla/include;/home/rmcgibbo/OpenMM5.0-Source/olla/include/openmm;/home/rmcgibbo/OpenMM5.0-Source/olla/include/openmm/internal;/home/rmcgibbo/OpenMM5.0-Source/serialization/include;/home/rmcgibbo/OpenMM5.0-Source/serialization/include/openmm;/home/rmcgibbo/OpenMM5.0-Source/serialization/include/openmm/internal;/home/rmcgibbo/OpenMM5.0-Source/plugins/amoeba/openmmapi/include;/home/rmcgibbo/OpenMM5.0-Source/plugins/amoeba/openmmapi/include/openmm;/home/rmcgibbo/OpenMM5.0-Source/plugins/amoeba/openmmapi/include/openmm/internal;/home/rmcgibbo/OpenMM5.0-Source/plugins/rpmd/openmmapi,/include;/home/rmcgibbo/OpenMM5.0-Source/plugins/rpmd/openmmapi,/include/openmm;/home/rmcgibbo/OpenMM5.0-Source/plugins/rpmd/openmmapi,/include/openmm/internal;/home/rmcgibbo/OpenMM5.0-Source/plugins/neb/openmmapi/include;/home/rmcgibbo/OpenMM5.0-Source/plugins/neb/openmmapi/include/openmm;/home/rmcgibbo/OpenMM5.0-Source/plugins/neb/openmmapi/include/openmm/internal")
file(TO_NATIVE_PATH "/home/rmcgibbo/opt/openmm/lib" OPENMM_LIB_PATH)
set(ENV{OPENMM_LIB_PATH} "${OPENMM_LIB_PATH}")
message("OPENMM_LIB_PATH = " $ENV{OPENMM_LIB_PATH})
message("OPENMM_INCLUDE_PATH = " $ENV{OPENMM_INCLUDE_PATH})
execute_process(
    COMMAND "/usr/bin/python" setup.py install
    WORKING_DIRECTORY ""
)
