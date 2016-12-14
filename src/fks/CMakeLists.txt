# files
set(SOURCEFILES
    ${SOURCEFILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/limits.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/luminosity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/remnants.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sfunctions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/splitting.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/virtual.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ximax.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/subtraction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fks_g.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xsec.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/radiationregion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/phasespaces.cpp
    PARENT_SCOPE
)

set(HEADERS 
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/limits.h
    ${CMAKE_CURRENT_SOURCE_DIR}/luminosity.h
    ${CMAKE_CURRENT_SOURCE_DIR}/color.h 
    ${CMAKE_CURRENT_SOURCE_DIR}/pdf.h.h
    ${CMAKE_CURRENT_SOURCE_DIR}/regions.h
    ${CMAKE_CURRENT_SOURCE_DIR}/process.h
    ${CMAKE_CURRENT_SOURCE_DIR}/remnants.h
    ${CMAKE_CURRENT_SOURCE_DIR}/sfunctions.h
    ${CMAKE_CURRENT_SOURCE_DIR}/splitting.h
    ${CMAKE_CURRENT_SOURCE_DIR}/virtual.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ximax.h
    ${CMAKE_CURRENT_SOURCE_DIR}/subtraction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/fks_g.h
    ${CMAKE_CURRENT_SOURCE_DIR}/xsec.h
    ${CMAKE_CURRENT_SOURCE_DIR}/radiationregion.h
    ${CMAKE_CURRENT_SOURCE_DIR}/phasespaces.h
    PARENT_SCOPE
)

#------ add tests -----------------------------------------------------------
include_directories(${gtest_SOURCE_DIR}/include ${LHAPDF_INCLUDE_DIR})

set(testname limits_test)
set(test limits_test.cpp)
add_executable(${testname} ${test})
target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
add_test(${testname} ${testname})

set(testname sfunctions_test)
set(test sfunctions_test.cpp)
add_executable(${testname} ${test})
target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
add_test(${testname} ${testname})

set(testname virtual_test)
set(test virtual_test.cpp)
add_executable(${testname} ${test})
target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
add_test(${testname} ${testname})

set(testname ximax_test)
set(test ximax_test.cpp)
add_executable(${testname} ${test})
target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
add_test(${testname} ${testname})

