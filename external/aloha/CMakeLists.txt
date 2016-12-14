set(SOURCEFILES aloha_aux_functions.cc  FFS1_3_2.cc  FFS2_3.cc  FFS4_1.cc  FFV1_3.cc     FFV2_3_2.cc  FFV2_4_3.cc  FFV3_1.cc  FFV4_3.cc   oxxxxx.cc  SSSS1_1.cc  VSS1_0.cc  VVS1_2.cc   VVSS1_4.cc  VVVV1_1.cc  VVVV2_2.cc  VVVV3_3.cc  VVVV4_4.cc  vxxxxx.cc
FFS1_0.cc                FFS1_3_3.cc  FFS3_0.cc  FFS4_2.cc  FFV2_0.cc     FFV2_3_3.cc  FFV2_5_0.cc  FFV3_2.cc  FFV5_0.cc   SSS1_0.cc   SSSS1_2.cc  VSS1_1.cc  VVS1_3.cc   VVV1_0.cc   VVVV1_2.cc  VVVV2_3.cc  VVVV3_4.cc  VVVV5_0.cc
FFS1_1.cc                FFS1_3.cc     FFS3_1.cc  FFS4_3.cc  FFV2_1.cc     FFV2_3.cc     FFV2_5_1.cc  FFV3_3.cc  FFV5_1.cc   SSS1_1.cc   SSSS1_3.cc  VSS1_2.cc  VVSS1_0.cc  VVV1_1.cc   VVVV1_3.cc  VVVV2_4.cc  VVVV4_0.cc  VVVV5_1.cc
FFS1_2.cc                FFS2_0.cc     FFS3_2.cc  FFV1_0.cc  FFV2_2.cc     FFV2_4_0.cc  FFV2_5_2.cc  FFV4_0.cc  FFV5_2.cc   SSS1_2.cc   SSSS1_4.cc  VSS1_3.cc  VVSS1_1.cc  VVV1_2.cc   VVVV1_4.cc  VVVV3_0.cc  VVVV4_1.cc  VVVV5_2.cc
FFS1_3_0.cc             FFS2_1.cc     FFS3_3.cc  FFV1_1.cc  FFV2_3_0.cc  FFV2_4_1.cc  FFV2_5_3.cc  FFV4_1.cc  FFV5_3.cc   SSS1_3.cc   sxxxxx.cc  VVS1_0.cc  VVSS1_2.cc  VVV1_3.cc   VVVV2_0.cc  VVVV3_1.cc  VVVV4_2.cc  VVVV5_3.cc
FFS1_3_1.cc             FFS2_2.cc     FFS4_0.cc  FFV1_2.cc  FFV2_3_1.cc  FFV2_4_2.cc  FFV3_0.cc     FFV4_2.cc  ixxxxx.cc  SSSS1_0.cc  txxxxx.cc  VVS1_1.cc  VVSS1_3.cc  VVVV1_0.cc  VVVV2_1.cc  VVVV3_2.cc  VVVV4_3.cc  VVVV5_4.cc)
set(HEADERS aloha_aux_functions.h  FFS1_3_1.h  FFS2_2.h  FFS4_0.h  FFV1_2.h    FFV2_3_1.h  FFV2_4_2.h  FFV3_0.h  FFV4_2.h  ixxxxx.h  SSSS1_0.h  txxxxx.h  VVS1_1.h   VVSS1_3.h  VVVV1_0.h  VVVV2_1.h  VVVV3_2.h  VVVV4_3.h  VVVV5_4.h
aloha.h                FFS1_3_2.h  FFS2_3.h  FFS4_1.h  FFV1_3.h    FFV2_3_2.h  FFV2_4_3.h  FFV3_1.h  FFV4_3.h  oxxxxx.h  SSSS1_1.h  VSS1_0.h  VVS1_2.h   VVSS1_4.h  VVVV1_1.h  VVVV2_2.h  VVVV3_3.h  VVVV4_4.h  vxxxxx.h
FFS1_0.h               FFS1_3_3.h  FFS3_0.h  FFS4_2.h  FFV2_0.h    FFV2_3_3.h  FFV2_5_0.h  FFV3_2.h  FFV5_0.h  SSS1_0.h  SSSS1_2.h  VSS1_1.h  VVS1_3.h   VVV1_0.h   VVVV1_2.h  VVVV2_3.h  VVVV3_4.h  VVVV5_0.h
FFS1_1.h               FFS1_3.h    FFS3_1.h  FFS4_3.h  FFV2_1.h    FFV2_3.h    FFV2_5_1.h  FFV3_3.h  FFV5_1.h  SSS1_1.h  SSSS1_3.h  VSS1_2.h  VVSS1_0.h  VVV1_1.h   VVVV1_3.h  VVVV2_4.h  VVVV4_0.h  VVVV5_1.h
FFS1_2.h               FFS2_0.h    FFS3_2.h  FFV1_0.h  FFV2_2.h    FFV2_4_0.h  FFV2_5_2.h  FFV4_0.h  FFV5_2.h  SSS1_2.h  SSSS1_4.h  VSS1_3.h  VVSS1_1.h  VVV1_2.h   VVVV1_4.h  VVVV3_0.h  VVVV4_1.h  VVVV5_2.h
FFS1_3_0.h             FFS2_1.h    FFS3_3.h  FFV1_1.h  FFV2_3_0.h  FFV2_4_1.h  FFV2_5_3.h  FFV4_1.h  FFV5_3.h  SSS1_3.h  sxxxxx.h   VVS1_0.h  VVSS1_2.h  VVV1_3.h   VVVV2_0.h  VVVV3_1.h  VVVV4_2.h  VVVV5_3.h)

add_library(aloha ${SOURCEFILES} ${HEADERS})
