      subroutine abbrrealQED15(resnr)
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

      
	WProd1268 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1264 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1265 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd1258 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,1,1)
        WProd1259 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,1,1)
        WProd1260 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd1261 = Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd1262 = Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
	WProd1257 = Prod(4,1,3,1)*ProdC(4,1,3,1)
        WProd1263 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,3,1)
        WProd1267 = Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd1269 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd1270 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,1,1))
        WProd1271 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(6,1,1,1))
        WProd1276 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,1,1))
	WProd1266 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd1272 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,3,1)
        WProd1273 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1274 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1275 = Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
      
      
	Prp153 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1257)
	Prp154 = 1.D0/(WProd1257 - WProd1264 - WProd1265)
	Prp155 = -(1.D0/WProd1266)
	Prp156 = 1.D0/(-WProd1266 + WProd1268 - WProd1269)
      
      
      
      endif

      end


