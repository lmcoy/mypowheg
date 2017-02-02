      subroutine abbrampsfacctEW3(resnr)
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

      
	WProd456 = Prod(2,2,2,1)
	WProd469 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd460 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd471 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd470 = ProdC(2,2,2,1)
	WProd465 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd466 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd461 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd472 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd467 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd459 = -(Prod(5,1,3,1)*ProdC(4,1,1,1))
        WProd473 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(4,1,1,1)
	WProd468 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd462 = Prod(4,1,2,2)*ProdC(4,1,2,1)
        WProd463 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd464 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd475 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,2,1))
	WProd474 = Prod(4,1,2,1)*ProdC(4,1,2,2)
        WProd476 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,2))
	WProd457 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd458 = Prod(5,1,2,1)*ProdC(5,1,2,1)
      
      
	Prp51 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd457)
	Prp52 = -(1.D0/WProd458)
	Prp53 = 1.D0/WProd469
      
      
      
      endif

      end


