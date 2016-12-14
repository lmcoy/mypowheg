# files
set(POWHEG_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/roverb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/upperbounding.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pickelement.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/btilde.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/generateevents.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/generateradiation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/reshuffle_momenta.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/unweighting.cpp
    PARENT_SCOPE
)

set(POWHEG_HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/roverb.h
    ${CMAKE_CURRENT_SOURCE_DIR}/upperbounding.h
    ${CMAKE_CURRENT_SOURCE_DIR}/pickelement.h
    ${CMAKE_CURRENT_SOURCE_DIR}/btilde.h 
    ${CMAKE_CURRENT_SOURCE_DIR}/generateevents.h
    ${CMAKE_CURRENT_SOURCE_DIR}/generateradiation.h
    ${CMAKE_CURRENT_SOURCE_DIR}/qcdfsr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/qcdisr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/qedfsr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/qedisr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/fsr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/isr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/solverqcdisr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/solverqedisr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ptgenerator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/radiationgenerator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/resonance.h
    ${CMAKE_CURRENT_SOURCE_DIR}/reshuffle_momenta.h
    ${CMAKE_CURRENT_SOURCE_DIR}/unweighting.h
    PARENT_SCOPE
)



