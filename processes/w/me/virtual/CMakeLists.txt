
set(SRC
abbrampsfacctEW1.F
abbrampsfacctQCD1.F
abbrborn1.F
abbrconst.F
abbrrealQED1.F
abbrrealQED2.F
abbrrealQED3.F
abbrrealQED4.F
abbrrealQED5.F
abbrrealQED6.F
abbrvirtEW1.F
abbrvirtQCD1.F
ampsfacctEW1.F
ampsfacctQCD1.F
born1.F
kinematics.F
loop_interface.F
realQED1.F
realQED2.F
realQED3.F
realQED4.F
realQED5.F
realQED6.F
squaredME.F
virtEW1.F
virtQCD1.F
)

set (CMAKE_Fortran_FLAGS_DEBUG ${CMAKE_Fortran_FLAGS} "-g")
set (CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} "-O2 -ffixed-line-length-200 -xf77-cpp-input -fno-align-commons")

add_library(wvirtual ${SRC})
target_link_libraries(wvirtual loop)
set_property(TARGET wvirtual PROPERTY FOLDER "processes/w/me/virtual")