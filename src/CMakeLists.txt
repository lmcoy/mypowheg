set(LIB_NAME mypowheg)

include_directories(SYSTEM ${LHAPDF_INCLUDE_DIR})
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/src)

set(HEADERS libconfig.h)

add_subdirectory(util)
add_subdirectory(config)
add_subdirectory(phasespace)
add_subdirectory(math)
add_subdirectory(physics)
add_subdirectory(fks)
add_subdirectory(powheg)
add_subdirectory(process)
add_subdirectory(random)
add_subdirectory(lhe)
add_subdirectory(pdf)

set(SOURCES
    ${CONFIG_SOURCES}
    ${UTIL_SOURCES}
    ${FKS_SOURCES}
    ${POWHEG_SOURCES}
    ${MATH_SOURCES}
    ${PHASESPACE_SOURCES}
    ${PHYSICS_SOURCES}
    ${LHE_SOURCES}
    ${PROCESS_SOURCES}
    ${RANDOM_SOURCES}
    ${PDF_SOURCES}
)

set(HEADERS
    ${CONFIG_HEADERS}
    ${UTIL_HEADERS}
    ${FKS_HEADERS}
    ${POWHEG_HEADERS}
    ${MATH_HEADERS}
    ${PHASESPACE_HEADERS}
    ${PHYSICS_HEADERS}
    ${LHE_HEADERS}
    ${PROCESS_HEADERS}
    ${RANDOM_HEADERS}
    ${PDF_HEADERS}
)

add_library(${LIB_NAME} ${SOURCES} ${HEADERS})

# group files for Xcode
source_group(config\\Sources FILES ${CONFIG_SOURCES})
source_group(config\\Headers FILES ${CONFIG_HEADERS})

source_group(util\\Sources FILES ${UTIL_SOURCES})
source_group(util\\Headers FILES ${UTIL_HEADERS})

source_group(fks\\Sources FILES ${FKS_SOURCES})
source_group(fks\\Headers FILES ${FKS_HEADERS})

source_group(powheg\\Sources FILES ${POWHEG_SOURCES})
source_group(powheg\\Headers FILES ${POWHEG_HEADERS})

source_group(math\\Sources FILES ${MATH_SOURCES})
source_group(math\\Headers FILES ${MATH_HEADERS})

source_group(phasespace\\Sources FILES ${PHASESPACE_SOURCES})
source_group(phasespace\\Headers FILES ${PHASESPACE_HEADERS})

source_group(physics\\Sources FILES ${PHYSICS_SOURCES})
source_group(physics\\Headers FILES ${PHYSICS_HEADERS})

source_group(lhe\\Sources FILES ${LHE_SOURCES})
source_group(lhe\\Headers FILES ${LHE_HEADERS})

source_group(process\\Sources FILES ${PROCESS_SOURCES})
source_group(process\\Headers FILES ${PROCESS_HEADERS})

source_group(random\\Sources FILES ${RANDOM_SOURCES})
source_group(random\\Headers FILES ${RANDOM_HEADERS})

source_group(pdf\\Sources FILES ${PDF_SOURCES})
source_group(pdf\\Headers FILES ${PDF_HEADERS})