add_subdirectory(virtual)
# --- Config ---

include_directories(${PROJECT_SOURCE_DIR}/src ${LHAPDF_INCLUDE_DIR})

# name of the library
set(LIB_NAME "mewj")
set(LLIBRARIES aloha_cm mypowheg wjvirtual)

# --- Build ---

set(SOURCEFILES
  ddx_dux.cpp # 7
  dxdx_uxdx.cpp # 9
  gdx_gux.cpp # 19
  gg_dux.cpp # 20
  gu_gd.cpp # 18
  uc_us.cpp # 10
  ud_dd.cpp # 3
  udx_gg.cpp # 1
  udx_uux.cpp # 5
  uu_ud.cpp # 2
  uux_scx.cpp # 12
  udx_ccx.cpp # 14
  udx_ddx.cpp # 6
  usx_ucx.cpp # 15
  uux_dux.cpp # 4
  uxdx_uxux.cpp # 8
  ucx_dcx.cpp # 13
  uxsx_uxcx.cpp # 16
  udx_ga.cpp
  gu_da.cpp
  gdx_uxa.cpp
  wg.cpp
  wq.cpp
  wqx.cpp
  virtual.cpp
  parameters_sm.cpp
)

set(HEADERS
  ddx_dux.h
  dxdx_uxdx.h
  gdx_gux.h
  gg_dux.h
  gu_gd.h
  uc_us.h
  ud_dd.h
  udx_gg.h
  udx_uux.h
  uu_ud.h
  uux_scx.h
  udx_ccx.h
  udx_ddx.h
  usx_ucx.h
  uux_dux.h
  uxdx_uxux.h
  ucx_dcx.h
  uxsx_uxcx.h
  udx_ga.h
  gu_da.h
  gdx_uxa.h
  vvv1p0_1.h
  ffv1p0_3.h
  wg.h
  wq.h
  wqx.h
  virtual.h
  parameters_sm.h
)



add_library(${LIB_NAME} ${SOURCEFILES} ${HEADERS})
target_link_libraries(${LIB_NAME} ${LLIBRARIES})
set_property(TARGET ${LIB_NAME} PROPERTY FOLDER "processes/wj/me")

# find and build tests
include_directories(${gtest_SOURCE_DIR}/include)

set(test ddx_dux_test.cpp)
set(testname ddx_dux_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test dxdx_uxdx_test.cpp)
set(testname dxdx_uxdx_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test gdx_gux_test.cpp)
set(testname gdx_gux_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test gg_dux_test.cpp)
set(testname gg_dux_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test gu_gd_test.cpp)
set(testname gu_gd_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test uc_us_test.cpp)
set(testname uc_us_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test ud_dd_test.cpp)
set(testname ud_dd_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test udx_gg_test.cpp)
set(testname udx_gg_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test udx_uux_test.cpp)
set(testname udx_uux_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test uu_ud_test.cpp)
set(testname uu_ud_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test uux_scx_test.cpp)
set(testname uux_scx_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test udx_ccx_test.cpp)
set(testname udx_ccx_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test udx_ddx_test.cpp)
set(testname udx_ddx_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test usx_ucx_test.cpp)
set(testname usx_ucx_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test uux_dux_test.cpp)
set(testname uux_dux_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test uxdx_uxux_test.cpp)
set(testname uxdx_uxux_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test ucx_dcx_test.cpp)
set(testname ucx_dcx_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test uxsx_uxcx_test.cpp)
set(testname uxsx_uxcx_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test wg_test.cpp)
set(testname wj_wg_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test wq_test.cpp)
set(testname wj_wq_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test wqx_test.cpp)
set(testname wj_wqx_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test udx_ga_test.cpp)
set(testname wj_udx_ga_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test gu_da_test.cpp)
set(testname wj_gu_da_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test gdx_uxa_test.cpp)
set(testname wj_gdx_uxa_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")

set(test virtual_test.cpp)
set(testname wj_virtual_test)
  add_executable(${testname} ${test})
  target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main)
  add_test(${testname} ${testname})
  set_property(TARGET ${testname} PROPERTY FOLDER "processes/wj/tests")
