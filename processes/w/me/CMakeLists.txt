add_subdirectory(virtual)
# --- Config ---

include_directories(${PROJECT_SOURCE_DIR}/src ${LHAPDF_INCLUDE_DIR})

# name of the library
set(LIB_NAME "mew")
set(LLIBRARIES aloha_cm mypowheg wvirtual)

# --- Build ---

set(SOURCEFILES virtual.cpp wg.cpp wa.cpp wqx.cpp wq.cpp parameters_sm.cpp)

set(HEADERS virtual.h wg.h wa.h wqx.h wq.h parameters_sm.h)


add_library(${LIB_NAME} ${SOURCEFILES} ${HEADERS})
target_link_libraries(${LIB_NAME} ${LLIBRARIES})
set_property(TARGET ${LIB_NAME} PROPERTY FOLDER "processes/w/me")

# find and build tests
include_directories(${gtest_SOURCE_DIR}/include)

set(test wg_test.cpp)
set(testname wg_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/w/tests")

set(test wa_test.cpp)
set(testname wa_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/w/tests")

set(test wqx_test.cpp)
set(testname wqx_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/w/tests")

set(test wq_test.cpp)
set(testname wq_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/w/tests")

set(test virtual_test.cpp)
set(testname w_virtual_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/w/tests")





