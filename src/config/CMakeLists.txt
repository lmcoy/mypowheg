# files
set(CONFIG_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/file.cpp
    PARENT_SCOPE
)

set(CONFIG_HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/file.h
    PARENT_SCOPE
)


include_directories(${gtest_SOURCE_DIR}/include)

# add test file_test.cpp
set(testname file_test)
set(test_srcs file_test.cpp)
add_executable(${testname} ${test_srcs})
target_link_libraries(${testname} ${LIB_NAME} ${LLIB} gtest gtest_main) 
add_test(${testname} ${testname})
set_property(TARGET ${testname} PROPERTY FOLDER "mypowheg/test/config")
