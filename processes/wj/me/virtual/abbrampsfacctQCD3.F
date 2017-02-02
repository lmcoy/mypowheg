      subroutine abbrampsfacctQCD3(resnr)
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

      
	WProd548 = Prod(2,2,2,1)
	WProd558 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd549 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd563 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd562 = ProdC(2,2,2,1)
	WProd553 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd554 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd559 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd567 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd555 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd550 = -(Prod(5,1,3,1)*ProdC(4,1,1,1))
        WProd564 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(4,1,1,1)
	WProd556 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd560 = Prod(4,1,2,2)*ProdC(4,1,2,1)
        WProd551 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd552 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd565 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,2,1))
	WProd568 = Prod(4,1,2,1)*ProdC(4,1,2,2)
        WProd566 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,2))
	WProd557 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd561 = Prod(5,1,2,1)*ProdC(5,1,2,1)
      
      
	Prp65 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd557)
	Prp66 = 1.D0/WProd558
	Prp67 = -(1.D0/WProd561)
      
      
      
      endif

      end


