# files
set(SOURCEFILES
    ${SOURCEFILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/stringutil.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/histogram.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/databuffer.cpp
    PARENT_SCOPE
)

set(HEADERS 
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/stringutil.h
    ${CMAKE_CURRENT_SOURCE_DIR}/matrix.h
    ${CMAKE_CURRENT_SOURCE_DIR}/histogram.h
    ${CMAKE_CURRENT_SOURCE_DIR}/byte.h databuffer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/staticmatrix.h
    PARENT_SCOPE
)

include_directories(${gtest_SOURCE_DIR}/include)

# add test histogram_test.cpp
set(testname histogram_test)
set(test_srcs histogram_test.cpp)
add_executable(${testname} ${test_srcs})
target_link_libraries(${testname} ${LIB_NAME} gtest gtest_main) 
add_test(${testname} ${testname})

# add test stringutil_test.cpp
set(testname stringutil_test)
set(test_srcs stringutil_test.cpp)
add_executable(${testname} ${test_srcs})
target_link_libraries(${testname} ${LIB_NAME} gtest gtest_main) 
add_test(${testname} ${testname})

# add test matrix_test.cpp
set(testname matrix_test)
set(test_srcs matrix_test.cpp)
add_executable(${testname} ${test_srcs})
target_link_libraries(${testname} ${LIB_NAME} gtest gtest_main) 
add_test(${testname} ${testname})