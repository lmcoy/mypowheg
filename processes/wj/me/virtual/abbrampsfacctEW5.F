      subroutine abbrampsfacctEW5(resnr)
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

      
	WProd477 = Prod(2,2,2,1)
	WProd490 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd482 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd494 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd491 = ProdC(2,2,2,1)
	WProd486 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd487 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd484 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd495 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd488 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd489 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd485 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd496 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd478 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd479 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd480 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd481 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd492 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
        WProd493 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
	WProd483 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
        WProd497 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
      
      
	Prp54 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd478)
	Prp55 = -(1.D0/WProd479)
	Prp56 = 1.D0/WProd490
      
      
      
      endif

      end


