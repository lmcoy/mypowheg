      subroutine abbrrealQED12(resnr)
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "params.h" 

      integer resnr,delta
      double complex myA0,myB0,myB1,myB00,myB11,myB001,myB111,myDB0
      double complex myDB1,myDB00,myDB11,myBval,myCval,myDval,myEval
      double complex Prp,myDeltaIR,myDeltaUV

      integer myBget,myCget,myDget,myEget

#include "amps/looptools.h" 
#include "amps/ampvars.h" 
#include "amps/ampcommon.h" 
#include "amps/ampcom.h" 
      if(resnr.eq.0)then

      
	WProd1185 = Prod(5,2,5,1)*Prod(6,2,6,1)
	WProd1091 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1093 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1094 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd1086 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd1095 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd1096 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd1075 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd1097 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd1174 = Prod(5,2,1,1)*ProdC(5,1,1,1)
        WProd1117 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd1172 = Prod(3,1,1,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd1130 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd1159 = Prod(3,1,2,1)*Prod(5,2,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd1098 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd1162 = Prod(5,2,2,1)*ProdC(5,1,2,1)
	WProd1099 = Prod(5,1,3,1)*ProdC(5,1,3,1)
	WProd1175 = Prod(5,2,3,1)*ProdC(5,1,3,1)
        WProd1118 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
	WProd1100 = Prod(5,1,4,1)*ProdC(5,1,4,1)
	WProd1176 = Prod(5,2,4,1)*ProdC(5,1,4,1)
        WProd1131 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1178 = Prod(3,1,2,1)*Prod(5,2,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1119 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1160 = -(Prod(4,1,3,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1179 = Prod(4,1,2,1)*Prod(5,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1200 = Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1161 = Prod(5,1,3,1)*Prod(5,2,2,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1181 = Prod(5,1,2,1)*Prod(5,2,3,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1084 = Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,5,2)
	WProd1113 = Prod(5,1,1,1)*ProdC(5,2,1,1)
        WProd1110 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,1,1)
        WProd1076 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,2,1,1)
        WProd1078 = Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -   ProdC(5,2,1,1)
	WProd1080 = Prod(5,1,2,1)*ProdC(5,2,2,1)
	WProd1114 = Prod(5,1,3,1)*ProdC(5,2,3,1)
	WProd1115 = Prod(5,1,4,1)*ProdC(5,2,4,1)
        WProd1121 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,2,4,1)
        WProd1079 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,2,4,1))
        WProd1122 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,4,1)
        WProd1124 = Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(5,2,4,1)
	WProd1132 = Prod(6,2,6,1)*ProdC(5,2,5,1)
	WProd1109 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd1133 = Prod(6,2,1,1)*ProdC(6,1,1,1)
        WProd1150 = Prod(3,1,1,1)*Prod(6,2,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd1134 = Prod(3,1,2,1)*Prod(6,2,2,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd1135 = Prod(5,1,3,1)*Prod(6,2,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,1,1)
        WProd1196 = Prod(5,2,3,1)*Prod(6,2,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,1,1)
        WProd1153 = Prod(5,1,3,1)*Prod(6,2,2,1)*ProdC(5,2,4,1)*
     -   ProdC(6,1,1,1)
	WProd1201 = Prod(6,1,2,1)*ProdC(6,1,2,1)
	WProd1136 = Prod(6,2,2,1)*ProdC(6,1,2,1)
	WProd1137 = Prod(6,2,3,1)*ProdC(6,1,3,1)
        WProd1197 = Prod(5,1,3,1)*Prod(5,2,3,1)*Prod(6,2,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,3,1)
        WProd1154 = Prod(5,1,3,1)**2*Prod(6,2,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,1,3,1)
	WProd1138 = Prod(6,2,4,1)*ProdC(6,1,4,1)
        WProd1139 = Prod(3,1,2,1)*Prod(6,2,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1140 = -(Prod(4,1,3,1)*Prod(6,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd1141 = Prod(4,1,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1155 = Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(5,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1186 = Prod(5,2,2,1)*Prod(6,2,3,1)*ProdC(5,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1187 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,3,1)*ProdC(6,1,4,1)
        WProd1198 = Prod(5,1,4,1)*Prod(5,2,3,1)*Prod(6,2,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,4,1)
        WProd1188 = Prod(5,1,4,1)*Prod(5,2,2,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,4,1)
        WProd1142 = Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(5,2,1,1)*
     -   ProdC(6,1,4,1)
        WProd1143 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,3,1)*ProdC(6,1,4,1)
        WProd1156 = Prod(5,1,3,1)*Prod(5,1,4,1)*Prod(6,2,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,1,4,1)
        WProd1144 = Prod(5,1,2,1)*Prod(5,1,4,1)*Prod(6,2,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,1,4,1)
	WProd1120 = Prod(6,1,5,1)*ProdC(6,1,5,1)
	WProd1195 = Prod(6,2,5,2)*ProdC(6,1,5,1)
        WProd1157 = -(Prod(5,1,3,1)*Prod(6,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,5,1))
        WProd1189 = Prod(5,2,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd1194 = Prod(3,1,1,1)*Prod(5,2,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,1,5,1)
        WProd1190 = Prod(3,1,2,1)*Prod(5,2,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,1,5,1)
        WProd1191 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,2,3,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,3,1)*ProdC(6,1,5,1)
        WProd1199 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,2,3,1)*
     -   ProdC(3,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,1)
        WProd1192 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,2,4,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,1)
        WProd1193 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,2,5,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,1)
	WProd1152 = Prod(6,2,5,1)*ProdC(6,1,5,2)
        WProd1145 = Prod(5,1,2,1)*Prod(6,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,2)
        WProd1151 = Prod(3,1,1,1)*Prod(5,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,1,5,2)
        WProd1146 = Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,2,5,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,1,5,2)
        WProd1147 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,3,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,3,1)*ProdC(6,1,5,2)
        WProd1158 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,3,1)*
     -   ProdC(3,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,2)
        WProd1148 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,4,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,2)
        WProd1149 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,2,5,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,1,5,2)
	WProd1077 = Prod(6,1,1,1)*ProdC(6,2,1,1)
        WProd1111 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,1,1)
        WProd1081 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(6,2,1,1)
        WProd1082 = Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,1,1)
        WProd1180 = Prod(5,2,3,1)*Prod(6,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,2,1,1)
        WProd1123 = Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,2,4,1)*
     -   ProdC(6,2,1,1)
	WProd1083 = Prod(6,1,2,1)*ProdC(6,2,2,1)
	WProd1085 = Prod(6,1,3,1)*ProdC(6,2,3,1)
        WProd1182 = Prod(5,1,3,1)*Prod(5,2,3,1)*Prod(6,1,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,3,1)
        WProd1125 = Prod(5,1,3,1)**2*Prod(6,1,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,2,3,1)
	WProd1087 = Prod(6,1,4,1)*ProdC(6,2,4,1)
        WProd1088 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,2,4,1)
        WProd1089 = -(Prod(4,1,3,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,4,1))
        WProd1090 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,4,1)
        WProd1126 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(6,2,4,1)
        WProd1163 = Prod(5,2,2,1)*Prod(6,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(6,2,4,1)
        WProd1164 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,3,1)*ProdC(6,2,4,1)
        WProd1183 = Prod(5,1,4,1)*Prod(5,2,3,1)*Prod(6,1,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,4,1)
        WProd1165 = Prod(5,1,4,1)*Prod(5,2,2,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,4,1)
        WProd1092 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(5,2,1,1)*
     -   ProdC(6,2,4,1)
        WProd1101 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,3,1)*ProdC(6,2,4,1)
        WProd1127 = Prod(5,1,3,1)*Prod(5,1,4,1)*Prod(6,1,2,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,2,4,1)
        WProd1102 = Prod(5,1,2,1)*Prod(5,1,4,1)*Prod(6,1,3,1)*
     -   ProdC(5,1,1,1)*ProdC(5,2,4,1)*ProdC(6,2,4,1)
	WProd1177 = Prod(6,1,5,2)*ProdC(6,2,5,1)
        WProd1128 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,2,5,1))
        WProd1166 = Prod(5,2,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,5,1)
        WProd1173 = Prod(3,1,1,1)*Prod(5,2,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,2,5,1)
        WProd1167 = Prod(3,1,2,1)*Prod(5,2,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,2,5,1)
        WProd1168 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,1,3,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,3,1)*ProdC(6,2,5,1)
        WProd1184 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,1,3,1)*
     -   ProdC(3,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,1)
        WProd1169 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,1,4,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,1)
        WProd1170 = Prod(5,1,3,1)*Prod(5,2,2,1)*Prod(6,1,5,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,1)
	WProd1116 = Prod(6,1,5,1)*ProdC(6,2,5,2)
        WProd1103 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,2,5,2)
        WProd1112 = Prod(3,1,1,1)*Prod(5,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,1,1)*ProdC(6,2,5,2)
        WProd1104 = Prod(3,1,2,1)*Prod(5,1,2,1)*Prod(6,1,5,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,1,1)*ProdC(6,2,5,2)
        WProd1105 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,3,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,3,1)*ProdC(6,2,5,2)
        WProd1129 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,3,1)*
     -   ProdC(3,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,2)
        WProd1106 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,4,1)*
     -   ProdC(4,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,2)
        WProd1107 = Prod(5,1,2,1)*Prod(5,1,3,1)*Prod(6,1,5,1)*
     -   ProdC(5,1,1,1)*ProdC(5,1,4,1)*ProdC(6,2,5,2)
	WProd1171 = Prod(5,2,5,1)*ProdC(6,2,6,1)
	WProd1108 = ProdC(5,2,5,1)*ProdC(6,2,6,1)
      
      
	Prp135 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1075)
	Prp138 = 1.D0/(WProd1075 - WProd1093 - WProd1095)
	Prp141 = -(1.D0/WProd1097)
	Prp136 = -(1.D0/WProd1098)
	Prp137 = -(1.D0/WProd1109)
	Prp139 = 1.D0/WProd1120
	Prp140 = 1.D0/(-WProd1097 - WProd1109 + WProd1120)
	Prp142 = -(1.D0/WProd1201)
      
      
      
      endif

      end


