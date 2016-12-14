# files
set(PHYSICS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/pdgcode.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/alphas.cpp
    PARENT_SCOPE
)

set(PHYSICS_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/pdgcode.h
    ${CMAKE_CURRENT_SOURCE_DIR}/alphas.h
    PARENT_SCOPE
)

# build tests
include_directories(${gtest_SOURCE_DIR}/include)

add_executable(pdgcode_test pdgcode_test.cpp)
target_link_libraries(pdgcode_test ${LIB_NAME} ${LLIBRARIES} gtest gtest_main) 
add_test(pdgcode_test pdgcode_test)
set_property(TARGET pdgcode_test PROPERTY FOLDER "mypowheg/tests/physics")

