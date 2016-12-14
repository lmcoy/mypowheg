# --- Config ---

# name of the library
set(LIB_NAME "powheggen")

include_directories(SYSTEM ${LHAPDF_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/vegas)

# files
set(SOURCEFILES eventgenerator.cpp)
set(HEADERS phasespacegenerator.h integraltransformation.h eventgenerator.h run.h)

include_directories(${MPI_CXX_INCLUDE_PATH})

add_library(${LIB_NAME} ${SOURCEFILES} ${HEADERS})
set_target_properties(${LIB_NAME} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
target_link_libraries(${LIB_NAME} vegas mypowheg ${MPI_CXX_LIBRARIES})




