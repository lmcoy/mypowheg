# --- Config ---

# name of the library
set(LIB_NAME "vegas")

# files
set(SOURCEFILES vegas_mpi.c)
set(HEADERS vegas.h)

include_directories(${MPI_C_INCLUDE_PATH})

add_library(${LIB_NAME} ${SOURCEFILES} ${HEADERS})
set_target_properties(${testname} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
target_link_libraries(${LIB_NAME} ${MPI_C_LIBRARIES})

include_directories(${gtest_SOURCE_DIR}/include)

# add test vegas
set(testname vegas_test)
set(test_srcs vegas_test.cpp)

add_executable(${testname} ${test_srcs})
set_target_properties(${testname} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
target_link_libraries(${testname} ${LIB_NAME} ${MPI_CXX_LIBRARIES} gtest gtest_main) 
add_test(NAME "${testname} 1 procs" COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${testname} ${MPIEXEC_POSTFLAGS})
add_test(NAME "${testname} 4 procs" COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${testname} ${MPIEXEC_POSTFLAGS})

