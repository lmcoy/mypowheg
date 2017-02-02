      subroutine abbrrealQED29(resnr)
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

      
	WProd2305 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd2306 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd2313 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd2307 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd2308 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd2319 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd2309 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd2314 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd2320 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd2310 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd2321 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd2322 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd2323 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd2311 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd2315 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd2325 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd2316 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd2326 = Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd2312 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd2324 = Prod(5,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,1)*
     -   ProdC(5,1,4,1)
	WProd2317 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd2318 = Prod(6,1,2,1)*ProdC(6,1,2,1)
      
      
	Prp241 = 1.D0/WProd2305
	Prp242 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd2306)
	Prp243 = 1.D0/(WProd2305 - WProd2313 - WProd2314)
	Prp244 = 1.D0/(WProd2305 - WProd2317 - WProd2318)
      
      
      
      endif

      end


