set(SOURCE coli_aux.F  coli_b0.F  coli_c0.F  coli_d0.F  coli_d0reg.F  coli_oint.F  ctoliserg.f  ctolis.f  intstd.f) 

set(HEADER checkparams_coli.h  common_coli.h  global_coli.h  params_coli.h)

set (CMAKE_Fortran_FLAGS_DEBUG ${CMAKE_Fortran_FLAGS} "-g")
set (CMAKE_Fortran_FLAGS "-O2")

add_library(loop ${SOURCE} ${HEADER})
