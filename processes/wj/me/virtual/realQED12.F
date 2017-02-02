      subroutine realQED12(amp,resnr)
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "params.h" 

      integer helnr,colnr,resnr,procnr

      double complex amp(nproc,ncol,nhel)

      double precision Alfas,Alfas2,GS

      common /my_alfas_common/Alfas,Alfas2,GS
#include "amps/ampvars.h" 
#include "amps/ampcom.h" 

*The Amplitude is evaluated as amp(procnr,colnr,helnr), where
*
*   Resonances:
*     resnr=0:   non-resonant
*
*   Partonic Subprocesses:
*     procnr=1:   {~up  ,  up  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=2:   {~up  , ~do  ,  ve  , ~el  , ~up  , ~up  }
*     procnr=3:   { up  ,  up  ,  ve  , ~el  ,  up  ,  do  }
*     procnr=4:   { up  , ~ch  ,  ve  , ~el  , ~ch  ,  do  }
*     procnr=5:   { up  ,  ch  ,  ve  , ~el  ,  ch  ,  do  }
*     procnr=6:   { up  , ~do  ,  ve  , ~el  , ~up  ,  up  }
*     procnr=7:   { up  , ~do  ,  ve  , ~el  , ~ch  ,  ch  }
*     procnr=8:   { up  , ~do  ,  ve  , ~el  , ~do  ,  do  }
*     procnr=9:   { up  , ~do  ,  ve  , ~el  , ~st  ,  st  }
*     procnr=10:   { up  , ~do  ,  ve  , ~el  , ~bo  ,  bo  }
*     procnr=11:   { up  , ~do  ,  ve  , ~el  ,  gl  ,  ga  }
*     procnr=12:   { up  , ~do  ,  ve  , ~el  ,  gl  ,  gl  }
*     procnr=13:   { up  ,  do  ,  ve  , ~el  ,  do  ,  do  }
*     procnr=14:   { up  , ~st  ,  ve  , ~el  ,  do  , ~st  }
*     procnr=15:   { up  ,  st  ,  ve  , ~el  ,  do  ,  st  }
*     procnr=16:   { up  ,  ga  ,  ve  , ~el  ,  do  ,  gl  }
*     procnr=17:   { up  ,  gl  ,  ve  , ~el  ,  do  ,  ga  }
*     procnr=18:   { up  ,  gl  ,  ve  , ~el  ,  do  ,  gl  }
*     procnr=19:   {~ch  ,  ch  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=20:   {~ch  , ~do  ,  ve  , ~el  , ~up  , ~ch  }
*     procnr=21:   { ch  , ~do  ,  ve  , ~el  , ~up  ,  ch  }
*     procnr=22:   {~do  , ~do  ,  ve  , ~el  , ~up  , ~do  }
*     procnr=23:   {~do  ,  do  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=24:   {~do  , ~st  ,  ve  , ~el  , ~up  , ~st  }
*     procnr=25:   {~do  ,  st  ,  ve  , ~el  , ~up  ,  st  }
*     procnr=26:   {~do  ,  ga  ,  ve  , ~el  , ~up  ,  gl  }
*     procnr=27:   {~do  ,  gl  ,  ve  , ~el  , ~up  ,  ga  }
*     procnr=28:   {~do  ,  gl  ,  ve  , ~el  , ~up  ,  gl  }
*     procnr=29:   {~st  ,  st  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=30:   { ga  ,  gl  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=31:   { gl  ,  gl  ,  ve  , ~el  , ~up  ,  do  }
*
*   Helicities:
*     helnr=1:   {-1, -1, -1, 1, -1, -1}
*     helnr=2:   {-1, 1, -1, 1, 1, -1}
*     helnr=3:   {1, -1, -1, 1, 1, -1}
*     helnr=4:   {-1, 1, -1, 1, -1, 1}
*     helnr=5:   {1, 1, -1, 1, 1, 1}
*     helnr=6:   {-1, 1, -1, 1, -1, -1}
*     helnr=7:   {-1, 1, -1, 1, 1, 1}
*     helnr=8:   {-1, -1, -1, 1, -1, 1}
*     helnr=9:   {1, -1, -1, 1, 1, 1}
*     helnr=10:   {1, 1, -1, 1, 1, -1}
*     helnr=11:   {-1, -1, -1, 1, 1, -1}
*
*   Color Structures:
*     colnr=1:   SUNT[Col1, Col3] SUNT[Col2, Col4]
*     colnr=2:   SUNT[Col1, Col2] SUNT[Col3, Col4]
*     colnr=3:   SUNT[Col1, Col4] SUNT[Col2, Col3]
*     colnr=4:   SUNT[Glu1, Col2, Col1]
*     colnr=5:   SUNT[Glu1, Glu2, Col2, Col1]
*     colnr=6:   SUNT[Glu2, Glu1, Col2, Col1]
*     colnr=7:   SUNT[Glu1, Col1, Col2]
*     colnr=8:   SUNT[Glu1, Glu2, Col1, Col2]
*     colnr=9:   SUNT[Glu2, Glu1, Col1, Col2]


      if(resnr.eq.0)then

        tmp1134 = (32*Alfas*Cpl122*Prp136*Prp137*
     -      (WProd1076*WProd1077 - WProd1077*WProd1078 + 
     -        WProd1077*WProd1079 - WProd1080*WProd1081 + 
     -        WProd1080*WProd1082 + WProd1076*WProd1083 - 
     -        WProd1078*WProd1083 + WProd1079*WProd1083 - 
     -        WProd1076*WProd1085 + WProd1078*WProd1085 - 
     -        WProd1079*WProd1085 - WProd1084*WProd1085 + 
     -        WProd1080*WProd1085*WProd1086 - 
     -        WProd1076*WProd1087 + WProd1078*WProd1087 - 
     -        WProd1079*WProd1087 - WProd1084*WProd1087 + 
     -        WProd1080*WProd1086*WProd1087 - 
     -        WProd1080*WProd1088 - WProd1080*WProd1089 - 
     -        WProd1080*WProd1090 - 2*WProd1075*WProd1092 - 
     -        2*WProd1091*WProd1092 + 2*WProd1092*WProd1093 + 
     -        2*WProd1092*WProd1094 + 2*WProd1092*WProd1095 + 
     -        2*WProd1092*WProd1096 + 2*WProd1092*WProd1097 + 
     -        WProd1092*WProd1098 - WProd1092*WProd1099 - 
     -        WProd1092*WProd1100 - WProd1101 - WProd1102 + 
     -        2*WProd1075*WProd1103 + 2*WProd1091*WProd1103 - 
     -        2*WProd1093*WProd1103 - 2*WProd1094*WProd1103 - 
     -        2*WProd1095*WProd1103 - 2*WProd1096*WProd1103 - 
     -        2*WProd1097*WProd1103 - WProd1098*WProd1103 + 
     -        WProd1099*WProd1103 + WProd1100*WProd1103 - 
     -        WProd1104 + WProd1105 + WProd1106 + WProd1107))/
     -    WProd1108 + Prp138*
     -    ((-32*Alfas*Cpl122*Prp136*
     -         (-(WProd1077*WProd1079) - WProd1079*WProd1083 + 
     -           WProd1079*WProd1085 + WProd1079*WProd1087 + 
     -           WProd1080*WProd1089 - 2*WProd1075*WProd1103 - 
     -           2*WProd1091*WProd1103 + 2*WProd1093*WProd1103 + 
     -           2*WProd1094*WProd1103 + 2*WProd1095*WProd1103 + 
     -           2*WProd1096*WProd1103 + 2*WProd1097*WProd1103 + 
     -           2*WProd1098*WProd1103 - WProd1099*WProd1103 - 
     -           WProd1100*WProd1103 - WProd1105 - WProd1106 + 
     -           WProd1077*WProd1110 + WProd1083*WProd1110 - 
     -           WProd1085*WProd1110 - WProd1087*WProd1110 - 
     -           WProd1080*WProd1111 - WProd1112))/WProd1108 - 
     -      (16*Alfas*Cpl122*Prp139*
     -         (-2*WProd1077*WProd1079 - 2*WProd1079*WProd1083 + 
     -           2*WProd1079*WProd1085 + 2*WProd1079*WProd1087 + 
     -           2*WProd1080*WProd1089 + 2*WProd1077*WProd1110 + 
     -           2*WProd1083*WProd1110 - 2*WProd1085*WProd1110 - 
     -           2*WProd1087*WProd1110 - 2*WProd1080*WProd1111 + 
     -           2*WProd1089*WProd1113 - 2*WProd1111*WProd1113 - 
     -           2*WProd1089*WProd1114 + 2*WProd1111*WProd1114 - 
     -           2*WProd1089*WProd1115 + 2*WProd1111*WProd1115 - 
     -           4*WProd1075*WProd1086*WProd1116 - 
     -           5*WProd1086*WProd1091*WProd1116 + 
     -           4*WProd1086*WProd1093*WProd1116 + 
     -           5*WProd1086*WProd1094*WProd1116 + 
     -           4*WProd1086*WProd1095*WProd1116 + 
     -           5*WProd1086*WProd1096*WProd1116 + 
     -           5*WProd1086*WProd1097*WProd1116 + 
     -           5*WProd1086*WProd1098*WProd1116 - 
     -           3*WProd1086*WProd1099*WProd1116 - 
     -           3*WProd1086*WProd1100*WProd1116 - 
     -           2*WProd1116*WProd1117 - 2*WProd1116*WProd1118 - 
     -           2*WProd1116*WProd1119))/WProd1108)
        tmp1134 = tmp1134 + 
     -   (32*Alfas*Cpl122*Prp137*Prp140*
     -      (WProd1076*WProd1077 + WProd1077*WProd1079 - 
     -        WProd1080*WProd1081 + WProd1076*WProd1083 + 
     -        WProd1079*WProd1083 - WProd1076*WProd1085 + 
     -        WProd1078*WProd1085 - WProd1079*WProd1085 - 
     -        WProd1084*WProd1085 + 
     -        WProd1080*WProd1085*WProd1086 - 
     -        WProd1076*WProd1087 + WProd1078*WProd1087 - 
     -        WProd1079*WProd1087 - WProd1084*WProd1087 + 
     -        WProd1080*WProd1086*WProd1087 - 
     -        WProd1080*WProd1088 - WProd1080*WProd1089 - 
     -        WProd1080*WProd1090 - 3*WProd1075*WProd1092 - 
     -        2*WProd1091*WProd1092 + 2*WProd1092*WProd1093 + 
     -        3*WProd1092*WProd1094 + 2*WProd1092*WProd1095 + 
     -        3*WProd1092*WProd1096 + 2*WProd1092*WProd1097 + 
     -        2*WProd1092*WProd1098 - 2*WProd1092*WProd1099 - 
     -        2*WProd1092*WProd1100 - WProd1101 - WProd1102 + 
     -        3*WProd1075*WProd1103 + 2*WProd1091*WProd1103 - 
     -        2*WProd1093*WProd1103 - 3*WProd1094*WProd1103 - 
     -        2*WProd1095*WProd1103 - 3*WProd1096*WProd1103 - 
     -        2*WProd1097*WProd1103 - 2*WProd1098*WProd1103 + 
     -        2*WProd1099*WProd1103 + 2*WProd1100*WProd1103 - 
     -        WProd1104 + WProd1105 + WProd1081*WProd1114 - 
     -        WProd1082*WProd1114 - 
     -        WProd1077*WProd1086*WProd1114 - 
     -        WProd1083*WProd1086*WProd1114 + 
     -        WProd1088*WProd1114 + WProd1089*WProd1114 + 
     -        WProd1090*WProd1114 + WProd1081*WProd1115 - 
     -        WProd1082*WProd1115 - 
     -        WProd1077*WProd1086*WProd1115 - 
     -        WProd1083*WProd1086*WProd1115 + 
     -        WProd1088*WProd1115 + WProd1089*WProd1115 + 
     -        WProd1090*WProd1115 + WProd1077*WProd1121 + 
     -        WProd1083*WProd1121 - WProd1085*WProd1121 - 
     -        WProd1087*WProd1121 + WProd1077*WProd1122 + 
     -        WProd1083*WProd1122 - WProd1085*WProd1122 - 
     -        WProd1087*WProd1122 + 3*WProd1075*WProd1123 + 
     -        2*WProd1091*WProd1123 - 2*WProd1093*WProd1123 - 
     -        3*WProd1094*WProd1123 - 2*WProd1095*WProd1123 - 
     -        3*WProd1096*WProd1123 - 2*WProd1097*WProd1123 - 
     -        2*WProd1098*WProd1123 + 2*WProd1099*WProd1123 + 
     -        2*WProd1100*WProd1123 - WProd1085*WProd1124 - 
     -        WProd1087*WProd1124 + WProd1125 + 
     -        WProd1114*WProd1126 + WProd1115*WProd1126 + 
     -        WProd1127 + WProd1114*WProd1128 + 
     -        WProd1115*WProd1128 - WProd1129))/WProd1108 + 
     -   (16*Alfas*Cpl122*Prp139*Prp140*
     -      (2*WProd1076*WProd1077 + 2*WProd1077*WProd1079 - 
     -        2*WProd1080*WProd1081 + 2*WProd1076*WProd1083 + 
     -        2*WProd1079*WProd1083 - 2*WProd1076*WProd1085 - 
     -        2*WProd1079*WProd1085 + 
     -        2*WProd1080*WProd1085*WProd1086 - 
     -        2*WProd1076*WProd1087 - 2*WProd1079*WProd1087 + 
     -        2*WProd1080*WProd1086*WProd1087 - 
     -        2*WProd1080*WProd1088 - 2*WProd1080*WProd1089 - 
     -        2*WProd1080*WProd1090 - 2*WProd1081*WProd1113 + 
     -        2*WProd1085*WProd1086*WProd1113 + 
     -        2*WProd1086*WProd1087*WProd1113 - 
     -        2*WProd1088*WProd1113 - 2*WProd1089*WProd1113 - 
     -        2*WProd1090*WProd1113 + 2*WProd1081*WProd1114 - 
     -        2*WProd1077*WProd1086*WProd1114 - 
     -        2*WProd1083*WProd1086*WProd1114 + 
     -        2*WProd1088*WProd1114 + 2*WProd1089*WProd1114 + 
     -        2*WProd1090*WProd1114 + 2*WProd1081*WProd1115 - 
     -        2*WProd1077*WProd1086*WProd1115 - 
     -        2*WProd1083*WProd1086*WProd1115 + 
     -        2*WProd1088*WProd1115 + 2*WProd1089*WProd1115 + 
     -        2*WProd1090*WProd1115 + 
     -        4*WProd1075*WProd1086*WProd1116 + 
     -        3*WProd1086*WProd1091*WProd1116 - 
     -        3*WProd1086*WProd1093*WProd1116 - 
     -        4*WProd1086*WProd1094*WProd1116 - 
     -        3*WProd1086*WProd1095*WProd1116 - 
     -        4*WProd1086*WProd1096*WProd1116 - 
     -        3*WProd1086*WProd1097*WProd1116 - 
     -        3*WProd1086*WProd1098*WProd1116 + 
     -        3*WProd1086*WProd1099*WProd1116 + 
     -        3*WProd1086*WProd1100*WProd1116 + 
     -        2*WProd1116*WProd1118 + 2*WProd1077*WProd1121 + 
     -        2*WProd1083*WProd1121 - 2*WProd1085*WProd1121 - 
     -        2*WProd1087*WProd1121 + 2*WProd1077*WProd1122 + 
     -        2*WProd1083*WProd1122 - 2*WProd1085*WProd1122 - 
     -        2*WProd1087*WProd1122 - 2*WProd1116*WProd1130 - 
     -        2*WProd1116*WProd1131))/WProd1108
	amp(12,5,6) = Prp135*tmp1134

#ifdef DEBUG
	print *, 'amp(12,5,6) =', amp(12,5,6)
#endif

        tmp1135 = (32*Alfas*Cpl122*Prp136*Prp137*
     -      (WProd1076*WProd1133 - WProd1078*WProd1133 + 
     -        WProd1079*WProd1133 - WProd1080*WProd1134 + 
     -        WProd1080*WProd1135 + WProd1076*WProd1136 - 
     -        WProd1078*WProd1136 + WProd1079*WProd1136 - 
     -        WProd1076*WProd1137 + WProd1078*WProd1137 - 
     -        WProd1079*WProd1137 - WProd1084*WProd1137 + 
     -        WProd1080*WProd1086*WProd1137 - 
     -        WProd1076*WProd1138 + WProd1078*WProd1138 - 
     -        WProd1079*WProd1138 - WProd1084*WProd1138 + 
     -        WProd1080*WProd1086*WProd1138 - 
     -        WProd1080*WProd1139 - WProd1080*WProd1140 - 
     -        WProd1080*WProd1141 - 2*WProd1075*WProd1142 - 
     -        2*WProd1091*WProd1142 + 2*WProd1093*WProd1142 + 
     -        2*WProd1094*WProd1142 + 2*WProd1095*WProd1142 + 
     -        2*WProd1096*WProd1142 + 2*WProd1097*WProd1142 + 
     -        WProd1098*WProd1142 - WProd1099*WProd1142 - 
     -        WProd1100*WProd1142 - WProd1143 - WProd1144 + 
     -        2*WProd1075*WProd1145 + 2*WProd1091*WProd1145 - 
     -        2*WProd1093*WProd1145 - 2*WProd1094*WProd1145 - 
     -        2*WProd1095*WProd1145 - 2*WProd1096*WProd1145 - 
     -        2*WProd1097*WProd1145 - WProd1098*WProd1145 + 
     -        WProd1099*WProd1145 + WProd1100*WProd1145 - 
     -        WProd1146 + WProd1147 + WProd1148 + WProd1149))/
     -    WProd1132 + (16*Alfas*Cpl122*Prp139*Prp140*
     -      (2*WProd1076*WProd1133 + 2*WProd1079*WProd1133 - 
     -        2*WProd1086*WProd1114*WProd1133 - 
     -        2*WProd1086*WProd1115*WProd1133 + 
     -        2*WProd1121*WProd1133 + 2*WProd1122*WProd1133 - 
     -        2*WProd1080*WProd1134 - 2*WProd1113*WProd1134 + 
     -        2*WProd1114*WProd1134 + 2*WProd1115*WProd1134 + 
     -        2*WProd1076*WProd1136 + 2*WProd1079*WProd1136 - 
     -        2*WProd1086*WProd1114*WProd1136 - 
     -        2*WProd1086*WProd1115*WProd1136 + 
     -        2*WProd1121*WProd1136 + 2*WProd1122*WProd1136 - 
     -        2*WProd1076*WProd1137 - 2*WProd1079*WProd1137 + 
     -        2*WProd1080*WProd1086*WProd1137 + 
     -        2*WProd1086*WProd1113*WProd1137 - 
     -        2*WProd1121*WProd1137 - 2*WProd1122*WProd1137 - 
     -        2*WProd1076*WProd1138 - 2*WProd1079*WProd1138 + 
     -        2*WProd1080*WProd1086*WProd1138 + 
     -        2*WProd1086*WProd1113*WProd1138 - 
     -        2*WProd1121*WProd1138 - 2*WProd1122*WProd1138 - 
     -        2*WProd1080*WProd1139 - 2*WProd1113*WProd1139 + 
     -        2*WProd1114*WProd1139 + 2*WProd1115*WProd1139 - 
     -        2*WProd1080*WProd1140 - 2*WProd1113*WProd1140 + 
     -        2*WProd1114*WProd1140 + 2*WProd1115*WProd1140 - 
     -        2*WProd1080*WProd1141 - 2*WProd1113*WProd1141 + 
     -        2*WProd1114*WProd1141 + 2*WProd1115*WProd1141 + 
     -        4*WProd1075*WProd1086*WProd1152 + 
     -        3*WProd1086*WProd1091*WProd1152 - 
     -        3*WProd1086*WProd1093*WProd1152 - 
     -        4*WProd1086*WProd1094*WProd1152 - 
     -        3*WProd1086*WProd1095*WProd1152 - 
     -        4*WProd1086*WProd1096*WProd1152 - 
     -        3*WProd1086*WProd1097*WProd1152 - 
     -        3*WProd1086*WProd1098*WProd1152 + 
     -        3*WProd1086*WProd1099*WProd1152 + 
     -        3*WProd1086*WProd1100*WProd1152 + 
     -        2*WProd1118*WProd1152 - 2*WProd1130*WProd1152 - 
     -        2*WProd1131*WProd1152))/WProd1132 + 
     -   Prp138*((-32*Alfas*Cpl122*Prp136*
     -         (-(WProd1079*WProd1133) + WProd1110*WProd1133 - 
     -           WProd1079*WProd1136 + WProd1110*WProd1136 + 
     -           WProd1079*WProd1137 - WProd1110*WProd1137 + 
     -           WProd1079*WProd1138 - WProd1110*WProd1138 + 
     -           WProd1080*WProd1140 - 2*WProd1075*WProd1145 - 
     -           2*WProd1091*WProd1145 + 2*WProd1093*WProd1145 + 
     -           2*WProd1094*WProd1145 + 2*WProd1095*WProd1145 + 
     -           2*WProd1096*WProd1145 + 2*WProd1097*WProd1145 + 
     -           2*WProd1098*WProd1145 - WProd1099*WProd1145 - 
     -           WProd1100*WProd1145 - WProd1147 - WProd1148 - 
     -           WProd1080*WProd1150 - WProd1151))/WProd1132 - 
     -      (16*Alfas*Cpl122*Prp139*
     -         (-2*WProd1079*WProd1133 + 2*WProd1110*WProd1133 - 
     -           2*WProd1079*WProd1136 + 2*WProd1110*WProd1136 + 
     -           2*WProd1079*WProd1137 - 2*WProd1110*WProd1137 + 
     -           2*WProd1079*WProd1138 - 2*WProd1110*WProd1138 + 
     -           2*WProd1080*WProd1140 + 2*WProd1113*WProd1140 - 
     -           2*WProd1114*WProd1140 - 2*WProd1115*WProd1140 - 
     -           2*WProd1080*WProd1150 - 2*WProd1113*WProd1150 + 
     -           2*WProd1114*WProd1150 + 2*WProd1115*WProd1150 - 
     -           4*WProd1075*WProd1086*WProd1152 - 
     -           5*WProd1086*WProd1091*WProd1152 + 
     -           4*WProd1086*WProd1093*WProd1152 + 
     -           5*WProd1086*WProd1094*WProd1152 + 
     -           4*WProd1086*WProd1095*WProd1152 + 
     -           5*WProd1086*WProd1096*WProd1152 + 
     -           5*WProd1086*WProd1097*WProd1152 + 
     -           5*WProd1086*WProd1098*WProd1152 - 
     -           3*WProd1086*WProd1099*WProd1152 - 
     -           3*WProd1086*WProd1100*WProd1152 - 
     -           2*WProd1117*WProd1152 - 2*WProd1118*WProd1152 - 
     -           2*WProd1119*WProd1152))/WProd1132)
        tmp1135 = tmp1135 + 
     -   (32*Alfas*Cpl122*Prp137*Prp140*
     -      (WProd1076*WProd1133 + WProd1079*WProd1133 - 
     -        WProd1086*WProd1114*WProd1133 - 
     -        WProd1086*WProd1115*WProd1133 + 
     -        WProd1121*WProd1133 + WProd1122*WProd1133 - 
     -        WProd1080*WProd1134 + WProd1114*WProd1134 + 
     -        WProd1115*WProd1134 - WProd1114*WProd1135 - 
     -        WProd1115*WProd1135 + WProd1076*WProd1136 + 
     -        WProd1079*WProd1136 - 
     -        WProd1086*WProd1114*WProd1136 - 
     -        WProd1086*WProd1115*WProd1136 + 
     -        WProd1121*WProd1136 + WProd1122*WProd1136 - 
     -        WProd1076*WProd1137 + WProd1078*WProd1137 - 
     -        WProd1079*WProd1137 - WProd1084*WProd1137 + 
     -        WProd1080*WProd1086*WProd1137 - 
     -        WProd1121*WProd1137 - WProd1122*WProd1137 - 
     -        WProd1124*WProd1137 - WProd1076*WProd1138 + 
     -        WProd1078*WProd1138 - WProd1079*WProd1138 - 
     -        WProd1084*WProd1138 + 
     -        WProd1080*WProd1086*WProd1138 - 
     -        WProd1121*WProd1138 - WProd1122*WProd1138 - 
     -        WProd1124*WProd1138 - WProd1080*WProd1139 + 
     -        WProd1114*WProd1139 + WProd1115*WProd1139 - 
     -        WProd1080*WProd1140 + WProd1114*WProd1140 + 
     -        WProd1115*WProd1140 - WProd1080*WProd1141 + 
     -        WProd1114*WProd1141 + WProd1115*WProd1141 - 
     -        3*WProd1075*WProd1142 - 2*WProd1091*WProd1142 + 
     -        2*WProd1093*WProd1142 + 3*WProd1094*WProd1142 + 
     -        2*WProd1095*WProd1142 + 3*WProd1096*WProd1142 + 
     -        2*WProd1097*WProd1142 + 2*WProd1098*WProd1142 - 
     -        2*WProd1099*WProd1142 - 2*WProd1100*WProd1142 - 
     -        WProd1143 - WProd1144 + 3*WProd1075*WProd1145 + 
     -        2*WProd1091*WProd1145 - 2*WProd1093*WProd1145 - 
     -        3*WProd1094*WProd1145 - 2*WProd1095*WProd1145 - 
     -        3*WProd1096*WProd1145 - 2*WProd1097*WProd1145 - 
     -        2*WProd1098*WProd1145 + 2*WProd1099*WProd1145 + 
     -        2*WProd1100*WProd1145 - WProd1146 + WProd1147 + 
     -        3*WProd1075*WProd1153 + 2*WProd1091*WProd1153 - 
     -        2*WProd1093*WProd1153 - 3*WProd1094*WProd1153 - 
     -        2*WProd1095*WProd1153 - 3*WProd1096*WProd1153 - 
     -        2*WProd1097*WProd1153 - 2*WProd1098*WProd1153 + 
     -        2*WProd1099*WProd1153 + 2*WProd1100*WProd1153 + 
     -        WProd1154 + WProd1114*WProd1155 + 
     -        WProd1115*WProd1155 + WProd1156 + 
     -        WProd1114*WProd1157 + WProd1115*WProd1157 - WProd1158
     -        ))/WProd1132
	amp(12,5,4) = Prp135*tmp1135

#ifdef DEBUG
	print *, 'amp(12,5,4) =', amp(12,5,4)
#endif

        tmp1136 = (32*Alfas*Cpl122*Prp136*Prp137*
     -      (WProd1077*WProd1159 + WProd1083*WProd1159 - 
     -        WProd1085*WProd1159 - WProd1087*WProd1159 + 
     -        WProd1077*WProd1160 + WProd1083*WProd1160 - 
     -        WProd1085*WProd1160 - WProd1087*WProd1160 - 
     -        WProd1077*WProd1161 - WProd1083*WProd1161 + 
     -        WProd1085*WProd1161 + WProd1087*WProd1161 - 
     -        WProd1081*WProd1162 + WProd1082*WProd1162 + 
     -        WProd1085*WProd1086*WProd1162 + 
     -        WProd1086*WProd1087*WProd1162 - 
     -        WProd1088*WProd1162 - WProd1089*WProd1162 - 
     -        WProd1090*WProd1162 - 2*WProd1075*WProd1163 - 
     -        2*WProd1091*WProd1163 + 2*WProd1093*WProd1163 + 
     -        2*WProd1094*WProd1163 + 2*WProd1095*WProd1163 + 
     -        2*WProd1096*WProd1163 + 2*WProd1097*WProd1163 + 
     -        WProd1098*WProd1163 - WProd1099*WProd1163 - 
     -        WProd1100*WProd1163 - WProd1164 - WProd1165 + 
     -        2*WProd1075*WProd1166 + 2*WProd1091*WProd1166 - 
     -        2*WProd1093*WProd1166 - 2*WProd1094*WProd1166 - 
     -        2*WProd1095*WProd1166 - 2*WProd1096*WProd1166 - 
     -        2*WProd1097*WProd1166 - WProd1098*WProd1166 + 
     -        WProd1099*WProd1166 + WProd1100*WProd1166 - 
     -        WProd1167 + WProd1168 + WProd1169 + WProd1170))/
     -    WProd1171 + Prp138*
     -    ((-32*Alfas*Cpl122*Prp136*
     -         (-(WProd1077*WProd1160) - WProd1083*WProd1160 + 
     -           WProd1085*WProd1160 + WProd1087*WProd1160 + 
     -           WProd1089*WProd1162 - WProd1111*WProd1162 - 
     -           2*WProd1075*WProd1166 - 2*WProd1091*WProd1166 + 
     -           2*WProd1093*WProd1166 + 2*WProd1094*WProd1166 + 
     -           2*WProd1095*WProd1166 + 2*WProd1096*WProd1166 + 
     -           2*WProd1097*WProd1166 + 2*WProd1098*WProd1166 - 
     -           WProd1099*WProd1166 - WProd1100*WProd1166 - 
     -           WProd1168 - WProd1169 + WProd1077*WProd1172 + 
     -           WProd1083*WProd1172 - WProd1085*WProd1172 - 
     -           WProd1087*WProd1172 - WProd1173))/WProd1171 - 
     -      (16*Alfas*Cpl122*Prp139*
     -         (-2*WProd1077*WProd1160 - 2*WProd1083*WProd1160 + 
     -           2*WProd1085*WProd1160 + 2*WProd1087*WProd1160 + 
     -           2*WProd1089*WProd1162 - 2*WProd1111*WProd1162 + 
     -           2*WProd1077*WProd1172 + 2*WProd1083*WProd1172 - 
     -           2*WProd1085*WProd1172 - 2*WProd1087*WProd1172 + 
     -           2*WProd1089*WProd1174 - 2*WProd1111*WProd1174 - 
     -           2*WProd1089*WProd1175 + 2*WProd1111*WProd1175 - 
     -           2*WProd1089*WProd1176 + 2*WProd1111*WProd1176 - 
     -           4*WProd1075*WProd1086*WProd1177 - 
     -           5*WProd1086*WProd1091*WProd1177 + 
     -           4*WProd1086*WProd1093*WProd1177 + 
     -           5*WProd1086*WProd1094*WProd1177 + 
     -           4*WProd1086*WProd1095*WProd1177 + 
     -           5*WProd1086*WProd1096*WProd1177 + 
     -           5*WProd1086*WProd1097*WProd1177 + 
     -           5*WProd1086*WProd1098*WProd1177 - 
     -           3*WProd1086*WProd1099*WProd1177 - 
     -           3*WProd1086*WProd1100*WProd1177 - 
     -           2*WProd1117*WProd1177 - 2*WProd1118*WProd1177 - 
     -           2*WProd1119*WProd1177))/WProd1171) + 
     -   (16*Alfas*Cpl122*Prp139*Prp140*
     -      (2*WProd1077*WProd1159 + 2*WProd1083*WProd1159 - 
     -        2*WProd1085*WProd1159 - 2*WProd1087*WProd1159 + 
     -        2*WProd1077*WProd1160 + 2*WProd1083*WProd1160 - 
     -        2*WProd1085*WProd1160 - 2*WProd1087*WProd1160 - 
     -        2*WProd1081*WProd1162 + 
     -        2*WProd1085*WProd1086*WProd1162 + 
     -        2*WProd1086*WProd1087*WProd1162 - 
     -        2*WProd1088*WProd1162 - 2*WProd1089*WProd1162 - 
     -        2*WProd1090*WProd1162 - 2*WProd1081*WProd1174 + 
     -        2*WProd1085*WProd1086*WProd1174 + 
     -        2*WProd1086*WProd1087*WProd1174 - 
     -        2*WProd1088*WProd1174 - 2*WProd1089*WProd1174 - 
     -        2*WProd1090*WProd1174 + 2*WProd1081*WProd1175 - 
     -        2*WProd1077*WProd1086*WProd1175 - 
     -        2*WProd1083*WProd1086*WProd1175 + 
     -        2*WProd1088*WProd1175 + 2*WProd1089*WProd1175 + 
     -        2*WProd1090*WProd1175 + 2*WProd1081*WProd1176 - 
     -        2*WProd1077*WProd1086*WProd1176 - 
     -        2*WProd1083*WProd1086*WProd1176 + 
     -        2*WProd1088*WProd1176 + 2*WProd1089*WProd1176 + 
     -        2*WProd1090*WProd1176 + 
     -        4*WProd1075*WProd1086*WProd1177 + 
     -        3*WProd1086*WProd1091*WProd1177 - 
     -        3*WProd1086*WProd1093*WProd1177 - 
     -        4*WProd1086*WProd1094*WProd1177 - 
     -        3*WProd1086*WProd1095*WProd1177 - 
     -        4*WProd1086*WProd1096*WProd1177 - 
     -        3*WProd1086*WProd1097*WProd1177 - 
     -        3*WProd1086*WProd1098*WProd1177 + 
     -        3*WProd1086*WProd1099*WProd1177 + 
     -        3*WProd1086*WProd1100*WProd1177 + 
     -        2*WProd1118*WProd1177 - 2*WProd1130*WProd1177 - 
     -        2*WProd1131*WProd1177 + 2*WProd1077*WProd1178 + 
     -        2*WProd1083*WProd1178 - 2*WProd1085*WProd1178 - 
     -        2*WProd1087*WProd1178 + 2*WProd1077*WProd1179 + 
     -        2*WProd1083*WProd1179 - 2*WProd1085*WProd1179 - 
     -        2*WProd1087*WProd1179))/WProd1171
        tmp1136 = tmp1136 + 
     -   (32*Alfas*Cpl122*Prp137*Prp140*
     -      (WProd1077*WProd1159 + WProd1083*WProd1159 - 
     -        WProd1085*WProd1159 - WProd1087*WProd1159 + 
     -        WProd1077*WProd1160 + WProd1083*WProd1160 - 
     -        WProd1085*WProd1160 - WProd1087*WProd1160 + 
     -        WProd1085*WProd1161 + WProd1087*WProd1161 - 
     -        WProd1081*WProd1162 + 
     -        WProd1085*WProd1086*WProd1162 + 
     -        WProd1086*WProd1087*WProd1162 - 
     -        WProd1088*WProd1162 - WProd1089*WProd1162 - 
     -        WProd1090*WProd1162 - 3*WProd1075*WProd1163 - 
     -        2*WProd1091*WProd1163 + 2*WProd1093*WProd1163 + 
     -        3*WProd1094*WProd1163 + 2*WProd1095*WProd1163 + 
     -        3*WProd1096*WProd1163 + 2*WProd1097*WProd1163 + 
     -        2*WProd1098*WProd1163 - 2*WProd1099*WProd1163 - 
     -        2*WProd1100*WProd1163 - WProd1164 - WProd1165 + 
     -        3*WProd1075*WProd1166 + 2*WProd1091*WProd1166 - 
     -        2*WProd1093*WProd1166 - 3*WProd1094*WProd1166 - 
     -        2*WProd1095*WProd1166 - 3*WProd1096*WProd1166 - 
     -        2*WProd1097*WProd1166 - 2*WProd1098*WProd1166 + 
     -        2*WProd1099*WProd1166 + 2*WProd1100*WProd1166 - 
     -        WProd1167 + WProd1168 + WProd1081*WProd1175 - 
     -        WProd1082*WProd1175 - 
     -        WProd1077*WProd1086*WProd1175 - 
     -        WProd1083*WProd1086*WProd1175 + 
     -        WProd1088*WProd1175 + WProd1089*WProd1175 + 
     -        WProd1090*WProd1175 + WProd1126*WProd1175 + 
     -        WProd1128*WProd1175 + WProd1081*WProd1176 - 
     -        WProd1082*WProd1176 - 
     -        WProd1077*WProd1086*WProd1176 - 
     -        WProd1083*WProd1086*WProd1176 + 
     -        WProd1088*WProd1176 + WProd1089*WProd1176 + 
     -        WProd1090*WProd1176 + WProd1126*WProd1176 + 
     -        WProd1128*WProd1176 + WProd1077*WProd1178 + 
     -        WProd1083*WProd1178 - WProd1085*WProd1178 - 
     -        WProd1087*WProd1178 + WProd1077*WProd1179 + 
     -        WProd1083*WProd1179 - WProd1085*WProd1179 - 
     -        WProd1087*WProd1179 + 3*WProd1075*WProd1180 + 
     -        2*WProd1091*WProd1180 - 2*WProd1093*WProd1180 - 
     -        3*WProd1094*WProd1180 - 2*WProd1095*WProd1180 - 
     -        3*WProd1096*WProd1180 - 2*WProd1097*WProd1180 - 
     -        2*WProd1098*WProd1180 + 2*WProd1099*WProd1180 + 
     -        2*WProd1100*WProd1180 - WProd1085*WProd1181 - 
     -        WProd1087*WProd1181 + WProd1182 + WProd1183 - 
     -        WProd1184))/WProd1171
	amp(12,5,2) = Prp135*tmp1136

#ifdef DEBUG
	print *, 'amp(12,5,2) =', amp(12,5,2)
#endif

        tmp1137 = (32*Alfas*Cpl122*Prp136*Prp137*
     -      (WProd1133*WProd1159 + WProd1136*WProd1159 - 
     -        WProd1137*WProd1159 - WProd1138*WProd1159 + 
     -        WProd1133*WProd1160 + WProd1136*WProd1160 - 
     -        WProd1137*WProd1160 - WProd1138*WProd1160 - 
     -        WProd1133*WProd1161 - WProd1136*WProd1161 + 
     -        WProd1137*WProd1161 + WProd1138*WProd1161 - 
     -        WProd1134*WProd1162 + WProd1135*WProd1162 + 
     -        WProd1086*WProd1137*WProd1162 + 
     -        WProd1086*WProd1138*WProd1162 - 
     -        WProd1139*WProd1162 - WProd1140*WProd1162 - 
     -        WProd1141*WProd1162 - 2*WProd1075*WProd1186 - 
     -        2*WProd1091*WProd1186 + 2*WProd1093*WProd1186 + 
     -        2*WProd1094*WProd1186 + 2*WProd1095*WProd1186 + 
     -        2*WProd1096*WProd1186 + 2*WProd1097*WProd1186 + 
     -        WProd1098*WProd1186 - WProd1099*WProd1186 - 
     -        WProd1100*WProd1186 - WProd1187 - WProd1188 + 
     -        2*WProd1075*WProd1189 + 2*WProd1091*WProd1189 - 
     -        2*WProd1093*WProd1189 - 2*WProd1094*WProd1189 - 
     -        2*WProd1095*WProd1189 - 2*WProd1096*WProd1189 - 
     -        2*WProd1097*WProd1189 - WProd1098*WProd1189 + 
     -        WProd1099*WProd1189 + WProd1100*WProd1189 - 
     -        WProd1190 + WProd1191 + WProd1192 + WProd1193))/
     -    WProd1185 + (16*Alfas*Cpl122*Prp139*Prp140*
     -      (2*WProd1133*WProd1159 + 2*WProd1136*WProd1159 - 
     -        2*WProd1137*WProd1159 - 2*WProd1138*WProd1159 + 
     -        2*WProd1133*WProd1160 + 2*WProd1136*WProd1160 - 
     -        2*WProd1137*WProd1160 - 2*WProd1138*WProd1160 - 
     -        2*WProd1134*WProd1162 + 
     -        2*WProd1086*WProd1137*WProd1162 + 
     -        2*WProd1086*WProd1138*WProd1162 - 
     -        2*WProd1139*WProd1162 - 2*WProd1140*WProd1162 - 
     -        2*WProd1141*WProd1162 - 2*WProd1134*WProd1174 + 
     -        2*WProd1086*WProd1137*WProd1174 + 
     -        2*WProd1086*WProd1138*WProd1174 - 
     -        2*WProd1139*WProd1174 - 2*WProd1140*WProd1174 - 
     -        2*WProd1141*WProd1174 - 
     -        2*WProd1086*WProd1133*WProd1175 + 
     -        2*WProd1134*WProd1175 - 
     -        2*WProd1086*WProd1136*WProd1175 + 
     -        2*WProd1139*WProd1175 + 2*WProd1140*WProd1175 + 
     -        2*WProd1141*WProd1175 - 
     -        2*WProd1086*WProd1133*WProd1176 + 
     -        2*WProd1134*WProd1176 - 
     -        2*WProd1086*WProd1136*WProd1176 + 
     -        2*WProd1139*WProd1176 + 2*WProd1140*WProd1176 + 
     -        2*WProd1141*WProd1176 + 2*WProd1133*WProd1178 + 
     -        2*WProd1136*WProd1178 - 2*WProd1137*WProd1178 - 
     -        2*WProd1138*WProd1178 + 2*WProd1133*WProd1179 + 
     -        2*WProd1136*WProd1179 - 2*WProd1137*WProd1179 - 
     -        2*WProd1138*WProd1179 + 
     -        4*WProd1075*WProd1086*WProd1195 + 
     -        3*WProd1086*WProd1091*WProd1195 - 
     -        3*WProd1086*WProd1093*WProd1195 - 
     -        4*WProd1086*WProd1094*WProd1195 - 
     -        3*WProd1086*WProd1095*WProd1195 - 
     -        4*WProd1086*WProd1096*WProd1195 - 
     -        3*WProd1086*WProd1097*WProd1195 - 
     -        3*WProd1086*WProd1098*WProd1195 + 
     -        3*WProd1086*WProd1099*WProd1195 + 
     -        3*WProd1086*WProd1100*WProd1195 + 
     -        2*WProd1118*WProd1195 - 2*WProd1130*WProd1195 - 
     -        2*WProd1131*WProd1195))/WProd1185 + 
     -   Prp138*((-32*Alfas*Cpl122*Prp136*
     -         (-(WProd1133*WProd1160) - WProd1136*WProd1160 + 
     -           WProd1137*WProd1160 + WProd1138*WProd1160 + 
     -           WProd1140*WProd1162 - WProd1150*WProd1162 + 
     -           WProd1133*WProd1172 + WProd1136*WProd1172 - 
     -           WProd1137*WProd1172 - WProd1138*WProd1172 - 
     -           2*WProd1075*WProd1189 - 2*WProd1091*WProd1189 + 
     -           2*WProd1093*WProd1189 + 2*WProd1094*WProd1189 + 
     -           2*WProd1095*WProd1189 + 2*WProd1096*WProd1189 + 
     -           2*WProd1097*WProd1189 + 2*WProd1098*WProd1189 - 
     -           WProd1099*WProd1189 - WProd1100*WProd1189 - 
     -           WProd1191 - WProd1192 - WProd1194))/WProd1185 - 
     -      (16*Alfas*Cpl122*Prp139*
     -         (-2*WProd1133*WProd1160 - 2*WProd1136*WProd1160 + 
     -           2*WProd1137*WProd1160 + 2*WProd1138*WProd1160 + 
     -           2*WProd1140*WProd1162 - 2*WProd1150*WProd1162 + 
     -           2*WProd1133*WProd1172 + 2*WProd1136*WProd1172 - 
     -           2*WProd1137*WProd1172 - 2*WProd1138*WProd1172 + 
     -           2*WProd1140*WProd1174 - 2*WProd1150*WProd1174 - 
     -           2*WProd1140*WProd1175 + 2*WProd1150*WProd1175 - 
     -           2*WProd1140*WProd1176 + 2*WProd1150*WProd1176 - 
     -           4*WProd1075*WProd1086*WProd1195 - 
     -           5*WProd1086*WProd1091*WProd1195 + 
     -           4*WProd1086*WProd1093*WProd1195 + 
     -           5*WProd1086*WProd1094*WProd1195 + 
     -           4*WProd1086*WProd1095*WProd1195 + 
     -           5*WProd1086*WProd1096*WProd1195 + 
     -           5*WProd1086*WProd1097*WProd1195 + 
     -           5*WProd1086*WProd1098*WProd1195 - 
     -           3*WProd1086*WProd1099*WProd1195 - 
     -           3*WProd1086*WProd1100*WProd1195 - 
     -           2*WProd1117*WProd1195 - 2*WProd1118*WProd1195 - 
     -           2*WProd1119*WProd1195))/WProd1185)
        tmp1137 = tmp1137 + 
     -   (32*Alfas*Cpl122*Prp137*Prp140*
     -      (WProd1133*WProd1159 + WProd1136*WProd1159 - 
     -        WProd1137*WProd1159 - WProd1138*WProd1159 + 
     -        WProd1133*WProd1160 + WProd1136*WProd1160 - 
     -        WProd1137*WProd1160 - WProd1138*WProd1160 + 
     -        WProd1137*WProd1161 + WProd1138*WProd1161 - 
     -        WProd1134*WProd1162 + 
     -        WProd1086*WProd1137*WProd1162 + 
     -        WProd1086*WProd1138*WProd1162 - 
     -        WProd1139*WProd1162 - WProd1140*WProd1162 - 
     -        WProd1141*WProd1162 - 
     -        WProd1086*WProd1133*WProd1175 + 
     -        WProd1134*WProd1175 - WProd1135*WProd1175 - 
     -        WProd1086*WProd1136*WProd1175 + 
     -        WProd1139*WProd1175 + WProd1140*WProd1175 + 
     -        WProd1141*WProd1175 + WProd1155*WProd1175 + 
     -        WProd1157*WProd1175 - 
     -        WProd1086*WProd1133*WProd1176 + 
     -        WProd1134*WProd1176 - WProd1135*WProd1176 - 
     -        WProd1086*WProd1136*WProd1176 + 
     -        WProd1139*WProd1176 + WProd1140*WProd1176 + 
     -        WProd1141*WProd1176 + WProd1155*WProd1176 + 
     -        WProd1157*WProd1176 + WProd1133*WProd1178 + 
     -        WProd1136*WProd1178 - WProd1137*WProd1178 - 
     -        WProd1138*WProd1178 + WProd1133*WProd1179 + 
     -        WProd1136*WProd1179 - WProd1137*WProd1179 - 
     -        WProd1138*WProd1179 - WProd1137*WProd1181 - 
     -        WProd1138*WProd1181 - 3*WProd1075*WProd1186 - 
     -        2*WProd1091*WProd1186 + 2*WProd1093*WProd1186 + 
     -        3*WProd1094*WProd1186 + 2*WProd1095*WProd1186 + 
     -        3*WProd1096*WProd1186 + 2*WProd1097*WProd1186 + 
     -        2*WProd1098*WProd1186 - 2*WProd1099*WProd1186 - 
     -        2*WProd1100*WProd1186 - WProd1187 - WProd1188 + 
     -        3*WProd1075*WProd1189 + 2*WProd1091*WProd1189 - 
     -        2*WProd1093*WProd1189 - 3*WProd1094*WProd1189 - 
     -        2*WProd1095*WProd1189 - 3*WProd1096*WProd1189 - 
     -        2*WProd1097*WProd1189 - 2*WProd1098*WProd1189 + 
     -        2*WProd1099*WProd1189 + 2*WProd1100*WProd1189 - 
     -        WProd1190 + WProd1191 + 3*WProd1075*WProd1196 + 
     -        2*WProd1091*WProd1196 - 2*WProd1093*WProd1196 - 
     -        3*WProd1094*WProd1196 - 2*WProd1095*WProd1196 - 
     -        3*WProd1096*WProd1196 - 2*WProd1097*WProd1196 - 
     -        2*WProd1098*WProd1196 + 2*WProd1099*WProd1196 + 
     -        2*WProd1100*WProd1196 + WProd1197 + WProd1198 - 
     -        WProd1199))/WProd1185
	amp(12,5,7) = Prp135*tmp1137

#ifdef DEBUG
	print *, 'amp(12,5,7) =', amp(12,5,7)
#endif

        tmp1138 = Prp138*
     -    ((16*Alfas*Cpl122*Prp139*
     -         (-2*WProd1077*WProd1079 - 2*WProd1079*WProd1083 + 
     -           2*WProd1079*WProd1085 + 2*WProd1079*WProd1087 + 
     -           2*WProd1080*WProd1089 + 2*WProd1077*WProd1110 + 
     -           2*WProd1083*WProd1110 - 2*WProd1085*WProd1110 - 
     -           2*WProd1087*WProd1110 - 2*WProd1080*WProd1111 + 
     -           2*WProd1089*WProd1113 - 2*WProd1111*WProd1113 - 
     -           2*WProd1089*WProd1114 + 2*WProd1111*WProd1114 - 
     -           2*WProd1089*WProd1115 + 2*WProd1111*WProd1115 - 
     -           4*WProd1075*WProd1086*WProd1116 - 
     -           5*WProd1086*WProd1091*WProd1116 + 
     -           4*WProd1086*WProd1093*WProd1116 + 
     -           5*WProd1086*WProd1094*WProd1116 + 
     -           4*WProd1086*WProd1095*WProd1116 + 
     -           5*WProd1086*WProd1096*WProd1116 + 
     -           5*WProd1086*WProd1097*WProd1116 + 
     -           5*WProd1086*WProd1098*WProd1116 - 
     -           3*WProd1086*WProd1099*WProd1116 - 
     -           3*WProd1086*WProd1100*WProd1116 - 
     -           2*WProd1116*WProd1117 - 2*WProd1116*WProd1118 - 
     -           2*WProd1116*WProd1119))/WProd1108 - 
     -      (32*Alfas*Cpl122*Prp142*
     -         (-2*WProd1075*WProd1103 - 3*WProd1091*WProd1103 + 
     -           2*WProd1093*WProd1103 + 3*WProd1094*WProd1103 + 
     -           2*WProd1095*WProd1103 + 3*WProd1096*WProd1103 + 
     -           3*WProd1097*WProd1103 + 3*WProd1098*WProd1103 - 
     -           2*WProd1099*WProd1103 - 2*WProd1100*WProd1103 - 
     -           WProd1105 - WProd1106 - WProd1112 - 
     -           WProd1089*WProd1113 + WProd1111*WProd1113 + 
     -           WProd1089*WProd1114 - WProd1111*WProd1114 + 
     -           WProd1089*WProd1115 - WProd1111*WProd1115 + 
     -           2*WProd1075*WProd1086*WProd1116 + 
     -           3*WProd1086*WProd1091*WProd1116 - 
     -           2*WProd1086*WProd1093*WProd1116 - 
     -           3*WProd1086*WProd1094*WProd1116 - 
     -           2*WProd1086*WProd1095*WProd1116 - 
     -           3*WProd1086*WProd1096*WProd1116 - 
     -           3*WProd1086*WProd1097*WProd1116 - 
     -           3*WProd1086*WProd1098*WProd1116 + 
     -           2*WProd1086*WProd1099*WProd1116 + 
     -           2*WProd1086*WProd1100*WProd1116 + 
     -           WProd1116*WProd1117 + WProd1116*WProd1118 + 
     -           WProd1116*WProd1119))/WProd1108) - 
     -   (16*Alfas*Cpl122*Prp139*Prp140*
     -      (2*WProd1076*WProd1077 + 2*WProd1077*WProd1079 - 
     -        2*WProd1080*WProd1081 + 2*WProd1076*WProd1083 + 
     -        2*WProd1079*WProd1083 - 2*WProd1076*WProd1085 - 
     -        2*WProd1079*WProd1085 + 
     -        2*WProd1080*WProd1085*WProd1086 - 
     -        2*WProd1076*WProd1087 - 2*WProd1079*WProd1087 + 
     -        2*WProd1080*WProd1086*WProd1087 - 
     -        2*WProd1080*WProd1088 - 2*WProd1080*WProd1089 - 
     -        2*WProd1080*WProd1090 - 2*WProd1081*WProd1113 + 
     -        2*WProd1085*WProd1086*WProd1113 + 
     -        2*WProd1086*WProd1087*WProd1113 - 
     -        2*WProd1088*WProd1113 - 2*WProd1089*WProd1113 - 
     -        2*WProd1090*WProd1113 + 2*WProd1081*WProd1114 - 
     -        2*WProd1077*WProd1086*WProd1114 - 
     -        2*WProd1083*WProd1086*WProd1114 + 
     -        2*WProd1088*WProd1114 + 2*WProd1089*WProd1114 + 
     -        2*WProd1090*WProd1114 + 2*WProd1081*WProd1115 - 
     -        2*WProd1077*WProd1086*WProd1115 - 
     -        2*WProd1083*WProd1086*WProd1115 + 
     -        2*WProd1088*WProd1115 + 2*WProd1089*WProd1115 + 
     -        2*WProd1090*WProd1115 + 
     -        4*WProd1075*WProd1086*WProd1116 + 
     -        3*WProd1086*WProd1091*WProd1116 - 
     -        3*WProd1086*WProd1093*WProd1116 - 
     -        4*WProd1086*WProd1094*WProd1116 - 
     -        3*WProd1086*WProd1095*WProd1116 - 
     -        4*WProd1086*WProd1096*WProd1116 - 
     -        3*WProd1086*WProd1097*WProd1116 - 
     -        3*WProd1086*WProd1098*WProd1116 + 
     -        3*WProd1086*WProd1099*WProd1116 + 
     -        3*WProd1086*WProd1100*WProd1116 + 
     -        2*WProd1116*WProd1118 + 2*WProd1077*WProd1121 + 
     -        2*WProd1083*WProd1121 - 2*WProd1085*WProd1121 - 
     -        2*WProd1087*WProd1121 + 2*WProd1077*WProd1122 + 
     -        2*WProd1083*WProd1122 - 2*WProd1085*WProd1122 - 
     -        2*WProd1087*WProd1122 - 2*WProd1116*WProd1130 - 
     -        2*WProd1116*WProd1131))/WProd1108
        tmp1138 = tmp1138 + 
     -   (32*Alfas*Cpl122*Prp140*Prp141*
     -      (WProd1078*WProd1085 - WProd1084*WProd1085 + 
     -        WProd1078*WProd1087 - WProd1084*WProd1087 - 
     -        WProd1075*WProd1092 - WProd1091*WProd1092 + 
     -        WProd1092*WProd1093 + WProd1092*WProd1094 + 
     -        WProd1092*WProd1095 + WProd1092*WProd1096 + 
     -        WProd1092*WProd1097 + WProd1092*WProd1098 - 
     -        WProd1092*WProd1099 - WProd1092*WProd1100 - 
     -        WProd1101 - WProd1102 + WProd1075*WProd1103 + 
     -        WProd1091*WProd1103 - WProd1093*WProd1103 - 
     -        WProd1094*WProd1103 - WProd1095*WProd1103 - 
     -        WProd1096*WProd1103 - WProd1097*WProd1103 - 
     -        WProd1098*WProd1103 + WProd1099*WProd1103 + 
     -        WProd1100*WProd1103 - WProd1104 + WProd1105 + 
     -        WProd1081*WProd1113 - 
     -        WProd1085*WProd1086*WProd1113 - 
     -        WProd1086*WProd1087*WProd1113 + 
     -        WProd1088*WProd1113 + WProd1089*WProd1113 + 
     -        WProd1090*WProd1113 - WProd1082*WProd1114 - 
     -        WProd1082*WProd1115 - 
     -        WProd1075*WProd1086*WProd1116 - 
     -        WProd1086*WProd1091*WProd1116 + 
     -        WProd1086*WProd1093*WProd1116 + 
     -        WProd1086*WProd1094*WProd1116 + 
     -        WProd1086*WProd1095*WProd1116 + 
     -        WProd1086*WProd1096*WProd1116 + 
     -        WProd1086*WProd1097*WProd1116 + 
     -        WProd1086*WProd1098*WProd1116 - 
     -        WProd1086*WProd1099*WProd1116 - 
     -        WProd1086*WProd1100*WProd1116 - 
     -        WProd1116*WProd1118 + WProd1075*WProd1123 + 
     -        WProd1091*WProd1123 - WProd1093*WProd1123 - 
     -        WProd1094*WProd1123 - WProd1095*WProd1123 - 
     -        WProd1096*WProd1123 - WProd1097*WProd1123 - 
     -        WProd1098*WProd1123 + WProd1099*WProd1123 + 
     -        WProd1100*WProd1123 - WProd1085*WProd1124 - 
     -        WProd1087*WProd1124 + WProd1125 + 
     -        WProd1114*WProd1126 + WProd1115*WProd1126 + 
     -        WProd1127 + WProd1114*WProd1128 + 
     -        WProd1115*WProd1128 - WProd1129 + 
     -        WProd1116*WProd1130 + WProd1116*WProd1131))/WProd1108
     -     - (32*Alfas*Cpl122*Prp141*Prp142*
     -      (WProd1097*WProd1103 - WProd1105 + WProd1107 - 
     -        WProd1112 - WProd1082*WProd1113 - 
     -        WProd1089*WProd1113 + WProd1111*WProd1113 + 
     -        WProd1082*WProd1114 + WProd1082*WProd1115 - 
     -        WProd1086*WProd1097*WProd1116 + 
     -        WProd1116*WProd1117 + WProd1116*WProd1118 + 
     -        WProd1097*WProd1123 - WProd1125 - WProd1127 - 
     -        WProd1114*WProd1128 - WProd1115*WProd1128 + 
     -        WProd1129 - WProd1116*WProd1131 - WProd1116*WProd1200
     -        ))/WProd1108
	amp(12,6,6) = Prp135*tmp1138

#ifdef DEBUG
	print *, 'amp(12,6,6) =', amp(12,6,6)
#endif

        tmp1139 = (-16*Alfas*Cpl122*Prp139*Prp140*
     -      (2*WProd1076*WProd1133 + 2*WProd1079*WProd1133 - 
     -        2*WProd1086*WProd1114*WProd1133 - 
     -        2*WProd1086*WProd1115*WProd1133 + 
     -        2*WProd1121*WProd1133 + 2*WProd1122*WProd1133 - 
     -        2*WProd1080*WProd1134 - 2*WProd1113*WProd1134 + 
     -        2*WProd1114*WProd1134 + 2*WProd1115*WProd1134 + 
     -        2*WProd1076*WProd1136 + 2*WProd1079*WProd1136 - 
     -        2*WProd1086*WProd1114*WProd1136 - 
     -        2*WProd1086*WProd1115*WProd1136 + 
     -        2*WProd1121*WProd1136 + 2*WProd1122*WProd1136 - 
     -        2*WProd1076*WProd1137 - 2*WProd1079*WProd1137 + 
     -        2*WProd1080*WProd1086*WProd1137 + 
     -        2*WProd1086*WProd1113*WProd1137 - 
     -        2*WProd1121*WProd1137 - 2*WProd1122*WProd1137 - 
     -        2*WProd1076*WProd1138 - 2*WProd1079*WProd1138 + 
     -        2*WProd1080*WProd1086*WProd1138 + 
     -        2*WProd1086*WProd1113*WProd1138 - 
     -        2*WProd1121*WProd1138 - 2*WProd1122*WProd1138 - 
     -        2*WProd1080*WProd1139 - 2*WProd1113*WProd1139 + 
     -        2*WProd1114*WProd1139 + 2*WProd1115*WProd1139 - 
     -        2*WProd1080*WProd1140 - 2*WProd1113*WProd1140 + 
     -        2*WProd1114*WProd1140 + 2*WProd1115*WProd1140 - 
     -        2*WProd1080*WProd1141 - 2*WProd1113*WProd1141 + 
     -        2*WProd1114*WProd1141 + 2*WProd1115*WProd1141 + 
     -        4*WProd1075*WProd1086*WProd1152 + 
     -        3*WProd1086*WProd1091*WProd1152 - 
     -        3*WProd1086*WProd1093*WProd1152 - 
     -        4*WProd1086*WProd1094*WProd1152 - 
     -        3*WProd1086*WProd1095*WProd1152 - 
     -        4*WProd1086*WProd1096*WProd1152 - 
     -        3*WProd1086*WProd1097*WProd1152 - 
     -        3*WProd1086*WProd1098*WProd1152 + 
     -        3*WProd1086*WProd1099*WProd1152 + 
     -        3*WProd1086*WProd1100*WProd1152 + 
     -        2*WProd1118*WProd1152 - 2*WProd1130*WProd1152 - 
     -        2*WProd1131*WProd1152))/WProd1132 + 
     -   Prp138*((16*Alfas*Cpl122*Prp139*
     -         (-2*WProd1079*WProd1133 + 2*WProd1110*WProd1133 - 
     -           2*WProd1079*WProd1136 + 2*WProd1110*WProd1136 + 
     -           2*WProd1079*WProd1137 - 2*WProd1110*WProd1137 + 
     -           2*WProd1079*WProd1138 - 2*WProd1110*WProd1138 + 
     -           2*WProd1080*WProd1140 + 2*WProd1113*WProd1140 - 
     -           2*WProd1114*WProd1140 - 2*WProd1115*WProd1140 - 
     -           2*WProd1080*WProd1150 - 2*WProd1113*WProd1150 + 
     -           2*WProd1114*WProd1150 + 2*WProd1115*WProd1150 - 
     -           4*WProd1075*WProd1086*WProd1152 - 
     -           5*WProd1086*WProd1091*WProd1152 + 
     -           4*WProd1086*WProd1093*WProd1152 + 
     -           5*WProd1086*WProd1094*WProd1152 + 
     -           4*WProd1086*WProd1095*WProd1152 + 
     -           5*WProd1086*WProd1096*WProd1152 + 
     -           5*WProd1086*WProd1097*WProd1152 + 
     -           5*WProd1086*WProd1098*WProd1152 - 
     -           3*WProd1086*WProd1099*WProd1152 - 
     -           3*WProd1086*WProd1100*WProd1152 - 
     -           2*WProd1117*WProd1152 - 2*WProd1118*WProd1152 - 
     -           2*WProd1119*WProd1152))/WProd1132 - 
     -      (32*Alfas*Cpl122*Prp142*
     -         (-(WProd1113*WProd1140) + WProd1114*WProd1140 + 
     -           WProd1115*WProd1140 - 2*WProd1075*WProd1145 - 
     -           3*WProd1091*WProd1145 + 2*WProd1093*WProd1145 + 
     -           3*WProd1094*WProd1145 + 2*WProd1095*WProd1145 + 
     -           3*WProd1096*WProd1145 + 3*WProd1097*WProd1145 + 
     -           3*WProd1098*WProd1145 - 2*WProd1099*WProd1145 - 
     -           2*WProd1100*WProd1145 - WProd1147 - WProd1148 + 
     -           WProd1113*WProd1150 - WProd1114*WProd1150 - 
     -           WProd1115*WProd1150 - WProd1151 + 
     -           2*WProd1075*WProd1086*WProd1152 + 
     -           3*WProd1086*WProd1091*WProd1152 - 
     -           2*WProd1086*WProd1093*WProd1152 - 
     -           3*WProd1086*WProd1094*WProd1152 - 
     -           2*WProd1086*WProd1095*WProd1152 - 
     -           3*WProd1086*WProd1096*WProd1152 - 
     -           3*WProd1086*WProd1097*WProd1152 - 
     -           3*WProd1086*WProd1098*WProd1152 + 
     -           2*WProd1086*WProd1099*WProd1152 + 
     -           2*WProd1086*WProd1100*WProd1152 + 
     -           WProd1117*WProd1152 + WProd1118*WProd1152 + 
     -           WProd1119*WProd1152))/WProd1132)
        tmp1139 = tmp1139 + 
     -   (32*Alfas*Cpl122*Prp140*Prp141*
     -      (WProd1113*WProd1134 - WProd1114*WProd1135 - 
     -        WProd1115*WProd1135 + WProd1078*WProd1137 - 
     -        WProd1084*WProd1137 - 
     -        WProd1086*WProd1113*WProd1137 - 
     -        WProd1124*WProd1137 + WProd1078*WProd1138 - 
     -        WProd1084*WProd1138 - 
     -        WProd1086*WProd1113*WProd1138 - 
     -        WProd1124*WProd1138 + WProd1113*WProd1139 + 
     -        WProd1113*WProd1140 + WProd1113*WProd1141 - 
     -        WProd1075*WProd1142 - WProd1091*WProd1142 + 
     -        WProd1093*WProd1142 + WProd1094*WProd1142 + 
     -        WProd1095*WProd1142 + WProd1096*WProd1142 + 
     -        WProd1097*WProd1142 + WProd1098*WProd1142 - 
     -        WProd1099*WProd1142 - WProd1100*WProd1142 - 
     -        WProd1143 - WProd1144 + WProd1075*WProd1145 + 
     -        WProd1091*WProd1145 - WProd1093*WProd1145 - 
     -        WProd1094*WProd1145 - WProd1095*WProd1145 - 
     -        WProd1096*WProd1145 - WProd1097*WProd1145 - 
     -        WProd1098*WProd1145 + WProd1099*WProd1145 + 
     -        WProd1100*WProd1145 - WProd1146 + WProd1147 - 
     -        WProd1075*WProd1086*WProd1152 - 
     -        WProd1086*WProd1091*WProd1152 + 
     -        WProd1086*WProd1093*WProd1152 + 
     -        WProd1086*WProd1094*WProd1152 + 
     -        WProd1086*WProd1095*WProd1152 + 
     -        WProd1086*WProd1096*WProd1152 + 
     -        WProd1086*WProd1097*WProd1152 + 
     -        WProd1086*WProd1098*WProd1152 - 
     -        WProd1086*WProd1099*WProd1152 - 
     -        WProd1086*WProd1100*WProd1152 - 
     -        WProd1118*WProd1152 + WProd1130*WProd1152 + 
     -        WProd1131*WProd1152 + WProd1075*WProd1153 + 
     -        WProd1091*WProd1153 - WProd1093*WProd1153 - 
     -        WProd1094*WProd1153 - WProd1095*WProd1153 - 
     -        WProd1096*WProd1153 - WProd1097*WProd1153 - 
     -        WProd1098*WProd1153 + WProd1099*WProd1153 + 
     -        WProd1100*WProd1153 + WProd1154 + 
     -        WProd1114*WProd1155 + WProd1115*WProd1155 + 
     -        WProd1156 + WProd1114*WProd1157 + 
     -        WProd1115*WProd1157 - WProd1158))/WProd1132 - 
     -   (32*Alfas*Cpl122*Prp141*Prp142*
     -      (-(WProd1113*WProd1135) + WProd1114*WProd1135 + 
     -        WProd1115*WProd1135 - WProd1113*WProd1140 + 
     -        WProd1097*WProd1145 - WProd1147 + WProd1149 + 
     -        WProd1113*WProd1150 - WProd1151 - 
     -        WProd1086*WProd1097*WProd1152 + 
     -        WProd1117*WProd1152 + WProd1118*WProd1152 - 
     -        WProd1131*WProd1152 + WProd1097*WProd1153 - 
     -        WProd1154 - WProd1156 - WProd1114*WProd1157 - 
     -        WProd1115*WProd1157 + WProd1158 - WProd1152*WProd1200
     -        ))/WProd1132
	amp(12,6,4) = Prp135*tmp1139

#ifdef DEBUG
	print *, 'amp(12,6,4) =', amp(12,6,4)
#endif

        tmp1140 = Prp138*
     -    ((16*Alfas*Cpl122*Prp139*
     -         (-2*WProd1077*WProd1160 - 2*WProd1083*WProd1160 + 
     -           2*WProd1085*WProd1160 + 2*WProd1087*WProd1160 + 
     -           2*WProd1089*WProd1162 - 2*WProd1111*WProd1162 + 
     -           2*WProd1077*WProd1172 + 2*WProd1083*WProd1172 - 
     -           2*WProd1085*WProd1172 - 2*WProd1087*WProd1172 + 
     -           2*WProd1089*WProd1174 - 2*WProd1111*WProd1174 - 
     -           2*WProd1089*WProd1175 + 2*WProd1111*WProd1175 - 
     -           2*WProd1089*WProd1176 + 2*WProd1111*WProd1176 - 
     -           4*WProd1075*WProd1086*WProd1177 - 
     -           5*WProd1086*WProd1091*WProd1177 + 
     -           4*WProd1086*WProd1093*WProd1177 + 
     -           5*WProd1086*WProd1094*WProd1177 + 
     -           4*WProd1086*WProd1095*WProd1177 + 
     -           5*WProd1086*WProd1096*WProd1177 + 
     -           5*WProd1086*WProd1097*WProd1177 + 
     -           5*WProd1086*WProd1098*WProd1177 - 
     -           3*WProd1086*WProd1099*WProd1177 - 
     -           3*WProd1086*WProd1100*WProd1177 - 
     -           2*WProd1117*WProd1177 - 2*WProd1118*WProd1177 - 
     -           2*WProd1119*WProd1177))/WProd1171 - 
     -      (32*Alfas*Cpl122*Prp142*
     -         (-2*WProd1075*WProd1166 - 3*WProd1091*WProd1166 + 
     -           2*WProd1093*WProd1166 + 3*WProd1094*WProd1166 + 
     -           2*WProd1095*WProd1166 + 3*WProd1096*WProd1166 + 
     -           3*WProd1097*WProd1166 + 3*WProd1098*WProd1166 - 
     -           2*WProd1099*WProd1166 - 2*WProd1100*WProd1166 - 
     -           WProd1168 - WProd1169 - WProd1173 - 
     -           WProd1089*WProd1174 + WProd1111*WProd1174 + 
     -           WProd1089*WProd1175 - WProd1111*WProd1175 + 
     -           WProd1089*WProd1176 - WProd1111*WProd1176 + 
     -           2*WProd1075*WProd1086*WProd1177 + 
     -           3*WProd1086*WProd1091*WProd1177 - 
     -           2*WProd1086*WProd1093*WProd1177 - 
     -           3*WProd1086*WProd1094*WProd1177 - 
     -           2*WProd1086*WProd1095*WProd1177 - 
     -           3*WProd1086*WProd1096*WProd1177 - 
     -           3*WProd1086*WProd1097*WProd1177 - 
     -           3*WProd1086*WProd1098*WProd1177 + 
     -           2*WProd1086*WProd1099*WProd1177 + 
     -           2*WProd1086*WProd1100*WProd1177 + 
     -           WProd1117*WProd1177 + WProd1118*WProd1177 + 
     -           WProd1119*WProd1177))/WProd1171) - 
     -   (16*Alfas*Cpl122*Prp139*Prp140*
     -      (2*WProd1077*WProd1159 + 2*WProd1083*WProd1159 - 
     -        2*WProd1085*WProd1159 - 2*WProd1087*WProd1159 + 
     -        2*WProd1077*WProd1160 + 2*WProd1083*WProd1160 - 
     -        2*WProd1085*WProd1160 - 2*WProd1087*WProd1160 - 
     -        2*WProd1081*WProd1162 + 
     -        2*WProd1085*WProd1086*WProd1162 + 
     -        2*WProd1086*WProd1087*WProd1162 - 
     -        2*WProd1088*WProd1162 - 2*WProd1089*WProd1162 - 
     -        2*WProd1090*WProd1162 - 2*WProd1081*WProd1174 + 
     -        2*WProd1085*WProd1086*WProd1174 + 
     -        2*WProd1086*WProd1087*WProd1174 - 
     -        2*WProd1088*WProd1174 - 2*WProd1089*WProd1174 - 
     -        2*WProd1090*WProd1174 + 2*WProd1081*WProd1175 - 
     -        2*WProd1077*WProd1086*WProd1175 - 
     -        2*WProd1083*WProd1086*WProd1175 + 
     -        2*WProd1088*WProd1175 + 2*WProd1089*WProd1175 + 
     -        2*WProd1090*WProd1175 + 2*WProd1081*WProd1176 - 
     -        2*WProd1077*WProd1086*WProd1176 - 
     -        2*WProd1083*WProd1086*WProd1176 + 
     -        2*WProd1088*WProd1176 + 2*WProd1089*WProd1176 + 
     -        2*WProd1090*WProd1176 + 
     -        4*WProd1075*WProd1086*WProd1177 + 
     -        3*WProd1086*WProd1091*WProd1177 - 
     -        3*WProd1086*WProd1093*WProd1177 - 
     -        4*WProd1086*WProd1094*WProd1177 - 
     -        3*WProd1086*WProd1095*WProd1177 - 
     -        4*WProd1086*WProd1096*WProd1177 - 
     -        3*WProd1086*WProd1097*WProd1177 - 
     -        3*WProd1086*WProd1098*WProd1177 + 
     -        3*WProd1086*WProd1099*WProd1177 + 
     -        3*WProd1086*WProd1100*WProd1177 + 
     -        2*WProd1118*WProd1177 - 2*WProd1130*WProd1177 - 
     -        2*WProd1131*WProd1177 + 2*WProd1077*WProd1178 + 
     -        2*WProd1083*WProd1178 - 2*WProd1085*WProd1178 - 
     -        2*WProd1087*WProd1178 + 2*WProd1077*WProd1179 + 
     -        2*WProd1083*WProd1179 - 2*WProd1085*WProd1179 - 
     -        2*WProd1087*WProd1179))/WProd1171
        tmp1140 = tmp1140 + 
     -   (32*Alfas*Cpl122*Prp140*Prp141*
     -      (WProd1085*WProd1161 + WProd1087*WProd1161 - 
     -        WProd1075*WProd1163 - WProd1091*WProd1163 + 
     -        WProd1093*WProd1163 + WProd1094*WProd1163 + 
     -        WProd1095*WProd1163 + WProd1096*WProd1163 + 
     -        WProd1097*WProd1163 + WProd1098*WProd1163 - 
     -        WProd1099*WProd1163 - WProd1100*WProd1163 - 
     -        WProd1164 - WProd1165 + WProd1075*WProd1166 + 
     -        WProd1091*WProd1166 - WProd1093*WProd1166 - 
     -        WProd1094*WProd1166 - WProd1095*WProd1166 - 
     -        WProd1096*WProd1166 - WProd1097*WProd1166 - 
     -        WProd1098*WProd1166 + WProd1099*WProd1166 + 
     -        WProd1100*WProd1166 - WProd1167 + WProd1168 + 
     -        WProd1081*WProd1174 - 
     -        WProd1085*WProd1086*WProd1174 - 
     -        WProd1086*WProd1087*WProd1174 + 
     -        WProd1088*WProd1174 + WProd1089*WProd1174 + 
     -        WProd1090*WProd1174 - WProd1082*WProd1175 + 
     -        WProd1126*WProd1175 + WProd1128*WProd1175 - 
     -        WProd1082*WProd1176 + WProd1126*WProd1176 + 
     -        WProd1128*WProd1176 - 
     -        WProd1075*WProd1086*WProd1177 - 
     -        WProd1086*WProd1091*WProd1177 + 
     -        WProd1086*WProd1093*WProd1177 + 
     -        WProd1086*WProd1094*WProd1177 + 
     -        WProd1086*WProd1095*WProd1177 + 
     -        WProd1086*WProd1096*WProd1177 + 
     -        WProd1086*WProd1097*WProd1177 + 
     -        WProd1086*WProd1098*WProd1177 - 
     -        WProd1086*WProd1099*WProd1177 - 
     -        WProd1086*WProd1100*WProd1177 - 
     -        WProd1118*WProd1177 + WProd1130*WProd1177 + 
     -        WProd1131*WProd1177 + WProd1075*WProd1180 + 
     -        WProd1091*WProd1180 - WProd1093*WProd1180 - 
     -        WProd1094*WProd1180 - WProd1095*WProd1180 - 
     -        WProd1096*WProd1180 - WProd1097*WProd1180 - 
     -        WProd1098*WProd1180 + WProd1099*WProd1180 + 
     -        WProd1100*WProd1180 - WProd1085*WProd1181 - 
     -        WProd1087*WProd1181 + WProd1182 + WProd1183 - 
     -        WProd1184))/WProd1171 - 
     -   (32*Alfas*Cpl122*Prp141*Prp142*
     -      (WProd1097*WProd1166 - WProd1168 + WProd1170 - 
     -        WProd1173 - WProd1082*WProd1174 - 
     -        WProd1089*WProd1174 + WProd1111*WProd1174 + 
     -        WProd1082*WProd1175 - WProd1128*WProd1175 + 
     -        WProd1082*WProd1176 - WProd1128*WProd1176 - 
     -        WProd1086*WProd1097*WProd1177 + 
     -        WProd1117*WProd1177 + WProd1118*WProd1177 - 
     -        WProd1131*WProd1177 + WProd1097*WProd1180 - 
     -        WProd1182 - WProd1183 + WProd1184 - 
     -        WProd1177*WProd1200))/WProd1171
	amp(12,6,2) = Prp135*tmp1140

#ifdef DEBUG
	print *, 'amp(12,6,2) =', amp(12,6,2)
#endif

        tmp1141 = (-16*Alfas*Cpl122*Prp139*Prp140*
     -      (2*WProd1133*WProd1159 + 2*WProd1136*WProd1159 - 
     -        2*WProd1137*WProd1159 - 2*WProd1138*WProd1159 + 
     -        2*WProd1133*WProd1160 + 2*WProd1136*WProd1160 - 
     -        2*WProd1137*WProd1160 - 2*WProd1138*WProd1160 - 
     -        2*WProd1134*WProd1162 + 
     -        2*WProd1086*WProd1137*WProd1162 + 
     -        2*WProd1086*WProd1138*WProd1162 - 
     -        2*WProd1139*WProd1162 - 2*WProd1140*WProd1162 - 
     -        2*WProd1141*WProd1162 - 2*WProd1134*WProd1174 + 
     -        2*WProd1086*WProd1137*WProd1174 + 
     -        2*WProd1086*WProd1138*WProd1174 - 
     -        2*WProd1139*WProd1174 - 2*WProd1140*WProd1174 - 
     -        2*WProd1141*WProd1174 - 
     -        2*WProd1086*WProd1133*WProd1175 + 
     -        2*WProd1134*WProd1175 - 
     -        2*WProd1086*WProd1136*WProd1175 + 
     -        2*WProd1139*WProd1175 + 2*WProd1140*WProd1175 + 
     -        2*WProd1141*WProd1175 - 
     -        2*WProd1086*WProd1133*WProd1176 + 
     -        2*WProd1134*WProd1176 - 
     -        2*WProd1086*WProd1136*WProd1176 + 
     -        2*WProd1139*WProd1176 + 2*WProd1140*WProd1176 + 
     -        2*WProd1141*WProd1176 + 2*WProd1133*WProd1178 + 
     -        2*WProd1136*WProd1178 - 2*WProd1137*WProd1178 - 
     -        2*WProd1138*WProd1178 + 2*WProd1133*WProd1179 + 
     -        2*WProd1136*WProd1179 - 2*WProd1137*WProd1179 - 
     -        2*WProd1138*WProd1179 + 
     -        4*WProd1075*WProd1086*WProd1195 + 
     -        3*WProd1086*WProd1091*WProd1195 - 
     -        3*WProd1086*WProd1093*WProd1195 - 
     -        4*WProd1086*WProd1094*WProd1195 - 
     -        3*WProd1086*WProd1095*WProd1195 - 
     -        4*WProd1086*WProd1096*WProd1195 - 
     -        3*WProd1086*WProd1097*WProd1195 - 
     -        3*WProd1086*WProd1098*WProd1195 + 
     -        3*WProd1086*WProd1099*WProd1195 + 
     -        3*WProd1086*WProd1100*WProd1195 + 
     -        2*WProd1118*WProd1195 - 2*WProd1130*WProd1195 - 
     -        2*WProd1131*WProd1195))/WProd1185 + 
     -   Prp138*((16*Alfas*Cpl122*Prp139*
     -         (-2*WProd1133*WProd1160 - 2*WProd1136*WProd1160 + 
     -           2*WProd1137*WProd1160 + 2*WProd1138*WProd1160 + 
     -           2*WProd1140*WProd1162 - 2*WProd1150*WProd1162 + 
     -           2*WProd1133*WProd1172 + 2*WProd1136*WProd1172 - 
     -           2*WProd1137*WProd1172 - 2*WProd1138*WProd1172 + 
     -           2*WProd1140*WProd1174 - 2*WProd1150*WProd1174 - 
     -           2*WProd1140*WProd1175 + 2*WProd1150*WProd1175 - 
     -           2*WProd1140*WProd1176 + 2*WProd1150*WProd1176 - 
     -           4*WProd1075*WProd1086*WProd1195 - 
     -           5*WProd1086*WProd1091*WProd1195 + 
     -           4*WProd1086*WProd1093*WProd1195 + 
     -           5*WProd1086*WProd1094*WProd1195 + 
     -           4*WProd1086*WProd1095*WProd1195 + 
     -           5*WProd1086*WProd1096*WProd1195 + 
     -           5*WProd1086*WProd1097*WProd1195 + 
     -           5*WProd1086*WProd1098*WProd1195 - 
     -           3*WProd1086*WProd1099*WProd1195 - 
     -           3*WProd1086*WProd1100*WProd1195 - 
     -           2*WProd1117*WProd1195 - 2*WProd1118*WProd1195 - 
     -           2*WProd1119*WProd1195))/WProd1185 - 
     -      (32*Alfas*Cpl122*Prp142*
     -         (-(WProd1140*WProd1174) + WProd1150*WProd1174 + 
     -           WProd1140*WProd1175 - WProd1150*WProd1175 + 
     -           WProd1140*WProd1176 - WProd1150*WProd1176 - 
     -           2*WProd1075*WProd1189 - 3*WProd1091*WProd1189 + 
     -           2*WProd1093*WProd1189 + 3*WProd1094*WProd1189 + 
     -           2*WProd1095*WProd1189 + 3*WProd1096*WProd1189 + 
     -           3*WProd1097*WProd1189 + 3*WProd1098*WProd1189 - 
     -           2*WProd1099*WProd1189 - 2*WProd1100*WProd1189 - 
     -           WProd1191 - WProd1192 - WProd1194 + 
     -           2*WProd1075*WProd1086*WProd1195 + 
     -           3*WProd1086*WProd1091*WProd1195 - 
     -           2*WProd1086*WProd1093*WProd1195 - 
     -           3*WProd1086*WProd1094*WProd1195 - 
     -           2*WProd1086*WProd1095*WProd1195 - 
     -           3*WProd1086*WProd1096*WProd1195 - 
     -           3*WProd1086*WProd1097*WProd1195 - 
     -           3*WProd1086*WProd1098*WProd1195 + 
     -           2*WProd1086*WProd1099*WProd1195 + 
     -           2*WProd1086*WProd1100*WProd1195 + 
     -           WProd1117*WProd1195 + WProd1118*WProd1195 + 
     -           WProd1119*WProd1195))/WProd1185)
        tmp1141 = tmp1141 + 
     -   (32*Alfas*Cpl122*Prp140*Prp141*
     -      (WProd1137*WProd1161 + WProd1138*WProd1161 + 
     -        WProd1134*WProd1174 - 
     -        WProd1086*WProd1137*WProd1174 - 
     -        WProd1086*WProd1138*WProd1174 + 
     -        WProd1139*WProd1174 + WProd1140*WProd1174 + 
     -        WProd1141*WProd1174 - WProd1135*WProd1175 + 
     -        WProd1155*WProd1175 + WProd1157*WProd1175 - 
     -        WProd1135*WProd1176 + WProd1155*WProd1176 + 
     -        WProd1157*WProd1176 - WProd1137*WProd1181 - 
     -        WProd1138*WProd1181 - WProd1075*WProd1186 - 
     -        WProd1091*WProd1186 + WProd1093*WProd1186 + 
     -        WProd1094*WProd1186 + WProd1095*WProd1186 + 
     -        WProd1096*WProd1186 + WProd1097*WProd1186 + 
     -        WProd1098*WProd1186 - WProd1099*WProd1186 - 
     -        WProd1100*WProd1186 - WProd1187 - WProd1188 + 
     -        WProd1075*WProd1189 + WProd1091*WProd1189 - 
     -        WProd1093*WProd1189 - WProd1094*WProd1189 - 
     -        WProd1095*WProd1189 - WProd1096*WProd1189 - 
     -        WProd1097*WProd1189 - WProd1098*WProd1189 + 
     -        WProd1099*WProd1189 + WProd1100*WProd1189 - 
     -        WProd1190 + WProd1191 - 
     -        WProd1075*WProd1086*WProd1195 - 
     -        WProd1086*WProd1091*WProd1195 + 
     -        WProd1086*WProd1093*WProd1195 + 
     -        WProd1086*WProd1094*WProd1195 + 
     -        WProd1086*WProd1095*WProd1195 + 
     -        WProd1086*WProd1096*WProd1195 + 
     -        WProd1086*WProd1097*WProd1195 + 
     -        WProd1086*WProd1098*WProd1195 - 
     -        WProd1086*WProd1099*WProd1195 - 
     -        WProd1086*WProd1100*WProd1195 - 
     -        WProd1118*WProd1195 + WProd1130*WProd1195 + 
     -        WProd1131*WProd1195 + WProd1075*WProd1196 + 
     -        WProd1091*WProd1196 - WProd1093*WProd1196 - 
     -        WProd1094*WProd1196 - WProd1095*WProd1196 - 
     -        WProd1096*WProd1196 - WProd1097*WProd1196 - 
     -        WProd1098*WProd1196 + WProd1099*WProd1196 + 
     -        WProd1100*WProd1196 + WProd1197 + WProd1198 - 
     -        WProd1199))/WProd1185 - 
     -   (32*Alfas*Cpl122*Prp141*Prp142*
     -      (-(WProd1135*WProd1174) - WProd1140*WProd1174 + 
     -        WProd1150*WProd1174 + WProd1135*WProd1175 - 
     -        WProd1157*WProd1175 + WProd1135*WProd1176 - 
     -        WProd1157*WProd1176 + WProd1097*WProd1189 - 
     -        WProd1191 + WProd1193 - WProd1194 - 
     -        WProd1086*WProd1097*WProd1195 + 
     -        WProd1117*WProd1195 + WProd1118*WProd1195 - 
     -        WProd1131*WProd1195 + WProd1097*WProd1196 - 
     -        WProd1197 - WProd1198 + WProd1199 - 
     -        WProd1195*WProd1200))/WProd1185
	amp(12,6,7) = Prp135*tmp1141

#ifdef DEBUG
	print *, 'amp(12,6,7) =', amp(12,6,7)
#endif

      endif

      end


