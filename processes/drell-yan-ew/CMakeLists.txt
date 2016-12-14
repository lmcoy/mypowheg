add_subdirectory(me)
# --- Config ---

include_directories(SYSTEM ${LHAPDF_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/vegas)

set(HEADERS cuts.h processlist.h  genbornps.h scales.h matrixelement.h zdata.h)
set(SRC main.cpp  cuts.cpp genbornps.cpp matrixelement.cpp zdata.cpp)

set(SRC_POWHEG powheg.cpp cuts.cpp genbornps.cpp matrixelement.cpp zdata.cpp)

execute_process(COMMAND mpic++ -showme:link OUTPUT_VARIABLE MPI_CXX_LINK OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND mpic++ -showme:compile OUTPUT_VARIABLE MPI_CXX_COMPILE OUTPUT_STRIP_TRAILING_WHITESPACE)

add_executable(dy_ew ${SRC} ${HEADERS})
set_target_properties(dy_ew PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE}")
target_link_libraries(dy_ew mypowheg vegas me_dy_ew ${LHAPDF_LIBRARIES} ${MPI_CXX_LINK} )
  set_property(TARGET dy_ew PROPERTY FOLDER "processes/drell-yan-ew")

add_executable(powheg_ex ${SRC_POWHEG} ${HEADERS})
set_target_properties(powheg_ex PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE}")
target_link_libraries(powheg_ex mypowheg vegas me_dy_ew ${LHAPDF_LIBRARIES} ${MPI_CXX_LINK})
  set_property(TARGET powheg_ex PROPERTY FOLDER "processes/drell-yan-ew")

