      subroutine abbrampsfacctQCD5(resnr)
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

      
	WProd596 = Prod(2,2,2,1)
	WProd604 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd597 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd611 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd610 = ProdC(2,2,2,1)
	WProd598 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd599 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd607 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd615 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd600 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd601 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd608 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd616 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd602 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd609 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd605 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd606 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd613 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
        WProd614 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
	WProd603 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
        WProd612 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
      
      
	Prp73 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd602)
	Prp74 = 1.D0/WProd604
	Prp75 = -(1.D0/WProd609)
      
      
      
      endif

      end


