# files
set(LHE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/eventbuffer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lhefile.cpp
    PARENT_SCOPE
)

set(LHE_HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/event.h
    ${CMAKE_CURRENT_SOURCE_DIR}/eventbuffer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/lhefile.h
    PARENT_SCOPE
)
