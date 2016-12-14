# files
set(SOURCEFILES
    ${SOURCEFILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/fourmomentum.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/util.cpp
    PARENT_SCOPE
)

set(HEADERS 
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/fourmomentum.h
    ${CMAKE_CURRENT_SOURCE_DIR}/util.h
    ${CMAKE_CURRENT_SOURCE_DIR}/dilog.h
    ${CMAKE_CURRENT_SOURCE_DIR}/math.h 
    PARENT_SCOPE
)


include_directories(${gtest_SOURCE_DIR}/include)

# add test dilog
set(testname dilog_test)
set(test_srcs dilog_test.cpp)
add_executable(${testname} ${test_srcs})
target_link_libraries(${testname} ${LIB_NAME} ${LLIB} gtest gtest_main) 
add_test(${testname} ${testname})

# add test util
set(testname util_test)
set(test_srcs util_test.cpp)
add_executable(${testname} ${test_srcs})
target_link_libraries(${testname} ${LIB_NAME} ${LLIB} gtest gtest_main) 
add_test(${testname} ${testname})

# add test fourmomentum
set(testname fourmomentum_test)
set(test_srcs fourmomentum_test.cpp)
add_executable(${testname} ${test_srcs})
target_link_libraries(${testname} ${LIB_NAME} ${LLIB} gtest gtest_main) 
add_test(${testname} ${testname})