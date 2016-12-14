

set(EW_VirtSrc  abbrampsfacctEWnonphotonic1.F  abbrampsfacctQCD1.F       abbrvirtEWnonphotonic2.F  abbrvirtEWphotonicISR2.F   ampsfacctEWphotonicFSR2.F  ampsfacctQCD2.F   virtEWnonphotonic1.F  virtEWphotonicISR1.F
abbrampsfacctEWnonphotonic2.F  abbrampsfacctQCD2.F       abbrvirtEWphotonicFSR1.F  abbrvirtQCD1.F             ampsfacctEWphotonicINT1.F  born1.F           virtEWnonphotonic2.F  virtEWphotonicISR2.F
abbrampsfacctEWphotonicFSR1.F  abbrborn1.F               abbrvirtEWphotonicFSR2.F  abbrvirtQCD2.F             ampsfacctEWphotonicINT2.F  born2.F           virtEWphotonicFSR1.F  virtQCD1.F
abbrampsfacctEWphotonicFSR2.F  abbrborn2.F               abbrvirtEWphotonicINT1.F  ampsfacctEWnonphotonic1.F  ampsfacctEWphotonicISR1.F  kinematics.F      virtEWphotonicFSR2.F  virtQCD2.F
abbrampsfacctEWphotonicISR1.F  abbrconst.F               abbrvirtEWphotonicINT2.F  ampsfacctEWnonphotonic2.F  ampsfacctEWphotonicISR2.F  loop_interface.F  virtEWphotonicINT1.F
abbrampsfacctEWphotonicISR2.F  abbrvirtEWnonphotonic1.F  abbrvirtEWphotonicISR1.F  ampsfacctEWphotonicFSR1.F  ampsfacctQCD1.F            squaredME.F       virtEWphotonicINT2.F)

set (CMAKE_Fortran_FLAGS_DEBUG ${CMAKE_Fortran_FLAGS} "-g")
set (CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} "-O2 -ffixed-line-length-200 -xf77-cpp-input -fno-align-commons")

add_library(ew_virtqed ${EW_VirtSrc})
target_link_libraries(ew_virtqed loop)
set_property(TARGET ew_virtqed PROPERTY FOLDER "processes/drell-yan-ew/me/VirtualQED")
