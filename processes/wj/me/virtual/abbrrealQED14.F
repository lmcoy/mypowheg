      subroutine abbrrealQED14(resnr)
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

      
	WProd1248 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1242 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1243 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd1236 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,1,1)
        WProd1245 = -(Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -     ProdC(4,1,1,1))
        WProd1237 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,1,1)
        WProd1238 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd1246 = -(Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd1239 = Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd1240 = Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
	WProd1235 = Prod(4,1,3,1)*ProdC(4,1,3,1)
        WProd1241 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,3,1)
        WProd1247 = Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd1249 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd1250 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,1,1))
        WProd1251 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(6,1,1,1))
        WProd1256 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,1,1))
	WProd1244 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd1252 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,3,1)
        WProd1253 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1254 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd1255 = Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
      
      
	Prp149 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1235)
	Prp150 = 1.D0/(WProd1235 - WProd1242 - WProd1243)
	Prp151 = -(1.D0/WProd1244)
	Prp152 = 1.D0/(-WProd1244 + WProd1248 - WProd1249)
      
      
      
      endif

      end


