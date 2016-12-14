# files
set(PHASESPACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/realphasespace.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/twoparticlegenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/recombination.cpp
    PARENT_SCOPE
)

set(PHASESPACE_HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/realphasespace.h
    ${CMAKE_CURRENT_SOURCE_DIR}/twoparticlegenerator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/internal.h
    ${CMAKE_CURRENT_SOURCE_DIR}/phasespace.h 
    ${CMAKE_CURRENT_SOURCE_DIR}/recombination.h
    PARENT_SCOPE
)

# find and build tests
include_directories(${gtest_SOURCE_DIR}/include)

set(testname twoparticlegenerator_test)
set(test twoparticlegenerator_test.cpp)
add_executable(${testname} ${test})
target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
add_test(${testname} ${testname})
set_property(TARGET ${testname} PROPERTY FOLDER "mypowheg/tests/phasespace")


set(testname realphasespace_test)
set(test realphasespace_test.cpp)
add_executable(${testname} ${test})
target_link_libraries(${testname} ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
add_test(${testname} ${testname})
set_property(TARGET ${testname} PROPERTY FOLDER "mypowheg/tests/phasespace")
