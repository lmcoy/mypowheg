add_subdirectory(me)

include_directories(SYSTEM ${LHAPDF_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/vegas ${PROJECT_SOURCE_DIR})

set(HEADERS  matrixelement.h cuts.h scales.h breitwignermapping.h bornpsgenerator.h myhistograms.h wdata.h)
set(SRC main.cpp matrixelement.cpp cuts.cpp wdata.cpp)

execute_process(COMMAND mpic++ -showme:link OUTPUT_VARIABLE MPI_CXX_LINK OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND mpic++ -showme:compile OUTPUT_VARIABLE MPI_CXX_COMPILE OUTPUT_STRIP_TRAILING_WHITESPACE)

add_executable(w ${SRC} ${HEADERS})
set_target_properties(w PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE}")
target_link_libraries(w mew mypowheg vegas powheggen ${LHAPDF_LIBRARIES} ${MPI_CXX_LINK} ) 
set_property(TARGET w PROPERTY FOLDER "processes/w")
