      subroutine abbrrealQED13(resnr)
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

      
	WProd1208 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1214 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1215 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd1203 = -(Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd1219 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,1,1)
        WProd1204 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,1,1)
        WProd1205 = -(Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1206 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd1216 = Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd1217 = Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
	WProd1202 = Prod(4,1,3,1)*ProdC(4,1,3,1)
        WProd1218 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,3,1)
	WProd1209 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd1227 = -(Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,1,1))
        WProd1233 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,1,1))
        WProd1228 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(5,1,1,1))
	WProd1210 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd1229 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd1207 = -(Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1211 = Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1230 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1231 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1232 = Prod(3,1,2,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1234 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -     ProdC(5,1,4,1))
	WProd1212 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd1220 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,1,1))
        WProd1221 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(6,1,1,1))
        WProd1226 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,1,1))
	WProd1213 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd1222 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,3,1)
        WProd1223 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1224 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1225 = Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
      
      
	Prp143 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1202)
	Prp144 = 1.D0/(WProd1208 - WProd1209 - WProd1210)
	Prp145 = -(1.D0/WProd1210)
	Prp146 = 1.D0/(WProd1208 - WProd1212 - WProd1213)
	Prp147 = -(1.D0/WProd1213)
	Prp148 = 1.D0/(WProd1202 - WProd1214 - WProd1215)
      
      
      
      endif

      end


