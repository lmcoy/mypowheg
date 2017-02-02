add_subdirectory(me)
include_directories(SYSTEM ${LHAPDF_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/vegas ${PROJECT_SOURCE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})


set(LLIBRARIES aloha_cm mypowheg vegas mewj powheggen)


set(HEADERS  phasespace.h matrixelement.h wjdata.h cuts.h myhistograms.h scales.h breitwignermapping.h)
set(SRC main.cpp matrixelement.cpp phasespace.cpp cuts.cpp wjdata.cpp)

add_executable(wj ${SRC} ${HEADERS})
set_target_properties(wj PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
target_link_libraries(wj ${LLIBRARIES} ${LHAPDF_LIBRARIES} ${MPI_CXX_LIBRARIES})

set(SRC_DEBUG debug_integrand.cpp matrixelement.cpp phasespace.cpp cuts.cpp wjdata.cpp)
add_executable(wj_debug ${SRC_DEBUG} ${HEADERS})
set_target_properties(wj_debug PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
target_link_libraries(wj_debug ${LLIBRARIES} ${LHAPDF_LIBRARIES} ${MPI_CXX_LIBRARIES})

add_executable(wj_genps genps_for_tests.cpp phasespace.cpp ${HEADERS})
set_target_properties(wj_genps PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
target_link_libraries(wj_genps ${LLIBRARIES} ${LHAPDF_LIBRARIES} ${MPI_CXX_LIBRARIES})

# find and build tests
include_directories(${gtest_SOURCE_DIR}/include)

include(ProcessorCount)
ProcessorCount(N)
if(N EQUAL 0)
  set(N 1)
endif()

set(test phasespace_test.cpp)
set(testname wj_phasespace_test)
  add_executable(${testname} ${test} phasespace.cpp)
  set_target_properties(${testname} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  target_link_libraries(${testname} ${LLIBRARIES} ${LHAPDF_LIBRARIES} ${MPI_CXX_LIBRARIES}  gtest gtest_main) 
  add_test(NAME ${testname} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${N} $<TARGET_FILE:${testname}> )
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")