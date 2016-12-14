add_subdirectory(VirtualQED)

# --- Config ---

include_directories(${PROJECT_SOURCE_DIR}/src ${LHAPDF_INCLUDE_DIR})

# name of the library
set(LIB_NAME "z_me")
set(LLIBRARIES aloha_cm mypowheg z_me_virt)



# --- Build ---

set(SOURCEFILES virtual_qed.cpp ddx_mupmuma.cpp uux_mupmuma.cpp 
parameters_sm.cpp uux_mupmumg.cpp ddx_mupmumg.cpp gd_mupmumd.cpp 
gu_mupmumu.cpp gux_mupmumux.cpp gdx_mupmumdx.cpp)

set(HEADERS virtual_qed.h ddx_mupmuma.h uux_mupmuma.h 
parameters_sm.h uux_mupmumg.h ddx_mupmumg.h gd_mupmumd.h 
gu_mupmumu.h gux_mupmumux.h gdx_mupmumdx.h)


add_library(${LIB_NAME} ${SOURCEFILES} ${HEADERS})
target_link_libraries(${LIB_NAME} ${LLIBRARIES})
set_property(TARGET ${LIB_NAME} PROPERTY FOLDER "processes/z/me")

# find and build tests
include_directories(${gtest_SOURCE_DIR}/include)

set(test uux_mupmuma_test.cpp)
set(testname z_uux_mupmuma_test)
  add_executable(ew_${testname} ${test})
  target_link_libraries(ew_${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(ew_${testname} ew_${testname})
  set_property(TARGET ew_${testname} PROPERTY FOLDER "processes/z/tests")

set(test uux_mupmumg_test.cpp)
set(testname z_uux_mupmumg_test)
  add_executable(ew_${testname} ${test})
  target_link_libraries(ew_${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(ew_${testname} ew_${testname})
  set_property(TARGET ew_${testname} PROPERTY FOLDER "processes/z/tests")

set(test ddx_mupmuma_test.cpp)
set(testname z_ddx_mupmuma_test)
  add_executable(ew_${testname} ${test})
  target_link_libraries(ew_${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(ew_${testname} ew_${testname})
  set_property(TARGET ew_${testname} PROPERTY FOLDER "processes/z/tests")

set(test ddx_mupmumg_test.cpp)
set(testname z_ddx_mupmumg_test)
  add_executable(ew_${testname} ${test})
  target_link_libraries(ew_${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(ew_${testname} ew_${testname})
  set_property(TARGET ew_${testname} PROPERTY FOLDER "processes/z/tests")

set(test gd_mupmumd_test.cpp)
set(testname z_gd_mupmumd_test)
  add_executable(ew_${testname} ${test})
  target_link_libraries(ew_${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(ew_${testname} ew_${testname})
  set_property(TARGET ew_${testname} PROPERTY FOLDER "processes/z/tests")

set(test gu_mupmumu_test.cpp)
set(testname z_gu_mupmumu_test)
  add_executable(ew_${testname} ${test})
  target_link_libraries(ew_${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(ew_${testname} ew_${testname})
  set_property(TARGET ew_${testname} PROPERTY FOLDER "processes/z/tests")

set(test gux_mupmumux_test.cpp)
set(testname z_gux_mupmumux_test)
  add_executable(ew_${testname} ${test})
  target_link_libraries(ew_${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(ew_${testname} ew_${testname})
  set_property(TARGET ew_${testname} PROPERTY FOLDER "processes/z/tests")

set(test gdx_mupmumdx_test.cpp)
set(testname z_gdx_mupmumdx_test)
  add_executable(ew_${testname} ${test})
  target_link_libraries(ew_${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(ew_${testname} ew_${testname})
  set_property(TARGET ew_${testname} PROPERTY FOLDER "processes/z/tests")

set(test virtual_qed_test.cpp)
set(testname z_virtual_qed_test)
  add_executable(ew_${testname} ${test})
  target_link_libraries(ew_${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
  add_test(ew_${testname} ew_${testname})
  set_property(TARGET ew_${testname} PROPERTY FOLDER "processes/z/tests")




