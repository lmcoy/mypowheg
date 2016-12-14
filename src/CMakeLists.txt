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

add_library(${LIB_NAME} ${SOURCEFILES})

