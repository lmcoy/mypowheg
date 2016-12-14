cmake_minimum_required (VERSION 2.6)

if(USE_CLANG)
  set(CMAKE_C_COMPILER clang)
  set(CMAKE_CXX_COMPILER clang++)
endif(USE_CLANG)

project(drell_yan_fks C CXX Fortran)

set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

option(USE_GTEST "enable unit tests" ON)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# get info from git repo
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --tags --always
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "DY_BUILD_VERSION"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} diff 
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "DY_GIT_DIFF"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    message(STATUS "Git version: ${DY_BUILD_VERSION}")
    if(DY_GIT_DIFF)
      string(REPLACE "\\"  "\\\\" DY_GIT_DIFF ${DY_GIT_DIFF})
      string(REPLACE "\""  "\\\"" DY_GIT_DIFF ${DY_GIT_DIFF})
      string(REPLACE "\n"  "\\n\" \\\n\"" DY_GIT_DIFF ${DY_GIT_DIFF})
      message(STATUS "Git modified: yes")
    else(DY_GIT_DIFF)
      message(STATUS "Git modified: no")
    endif(DY_GIT_DIFF)

  else(GIT_FOUND)
    set(DY_BUILD_VERSION 0)
    set(DY_GIT_DIFF "")
  endif(GIT_FOUND)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/libconfig.h.in ${CMAKE_CURRENT_BINARY_DIR}/libconfig.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

if(USE_GTEST)
  enable_testing()
  # Google test support
  add_subdirectory(external/gtest-1.7.0)
endif(USE_GTEST)

find_package(MPI REQUIRED)

find_package(Boost REQUIRED)

#--- LHAPDF -------------------------------------------------------------------------------------------------
find_program(LHAPDF_CONFIG_EXEC NAMES lhapdf-config)
if (NOT LHAPDF_CONFIG_EXEC)
    set(LHAPDF_FOUND FALSE)
    message(FATAL_ERROR "lhapdf-config not found")
else()
    set(LHAPDF_FOUND TRUE)
endif()
execute_process(COMMAND ${LHAPDF_CONFIG_EXEC} --ldflags OUTPUT_VARIABLE LHAPDF_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LHAPDF_CONFIG_EXEC} --cppflags OUTPUT_VARIABLE LHAPDF_CPPFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LHAPDF_CONFIG_EXEC} --incdir OUTPUT_VARIABLE LHAPDF_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LHAPDF_CONFIG_EXEC} --libdir OUTPUT_VARIABLE LHAPDF_LIBRARY_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
#------------------------------------------------------------------------------------------------------------

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall")
set(CMAKE_C_FLAGS "-Wall")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ggdb")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -ggdb")
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG "-ggdb -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -march=native")

set(CMAKE_C_FLAGS_DEBUG "-ggdb -O0")
set(CMAKE_C_FLAGS_RELEASE "-O3 -ffast-math -march=native")

add_subdirectory(external/aloha_cm)
add_subdirectory(external/LoopLib)
add_subdirectory(vegas)
add_subdirectory(src)
add_subdirectory(powheggen)
add_subdirectory(processes)




