add_subdirectory(me)
# --- Config ---

include_directories(SYSTEM ${LHAPDF_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/vegas ${PROJECT_SOURCE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})

set(HEADERS scales.h cuts.h matrixelement.h zdata.h bornpsgenerator.h breitwignermapping.h)
set(SRC main.cpp cuts.cpp matrixelement.cpp zdata.cpp)

# copy config file
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/config.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_executable(z ${SRC} ${HEADERS})
set_target_properties(z PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set_target_properties(z PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
target_link_libraries(z mypowheg vegas z_me powheggen ${LHAPDF_LIBRARIES} ${MPI_CXX_LIBRARIES_FLAGS} )
set_property(TARGET z PROPERTY FOLDER "processes/z")


