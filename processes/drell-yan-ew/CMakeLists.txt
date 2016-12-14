add_subdirectory(me)
# --- Config ---

include_directories(SYSTEM ${LHAPDF_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/vegas)

set(HEADERS cuts.h processes.h processlist.h  genbornps.h scales.h)
set(SRC main.cpp data.cpp processes.cpp  cuts.cpp genbornps.cpp)

set(SRC_POWHEG powheg.cpp data.cpp processes.cpp cuts.cpp genbornps.cpp)

execute_process(COMMAND mpic++ -showme:link OUTPUT_VARIABLE MPI_CXX_LINK OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND mpic++ -showme:compile OUTPUT_VARIABLE MPI_CXX_COMPILE OUTPUT_STRIP_TRAILING_WHITESPACE)

add_executable(dy_ew ${SRC} ${HEADERS})
set_target_properties(dy_ew PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE}")
target_link_libraries(dy_ew mypowheg vegas me_dy_ew ${LHAPDF_LIBRARIES} ${MPI_CXX_LINK} ) 

add_executable(powheg_ex ${SRC_POWHEG} ${HEADERS})
set_target_properties(powheg_ex PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE}")
target_link_libraries(powheg_ex mypowheg vegas me_dy_ew ${LHAPDF_LIBRARIES} ${MPI_CXX_LINK}) 

