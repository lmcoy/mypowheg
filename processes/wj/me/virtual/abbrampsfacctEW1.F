      subroutine abbrampsfacctEW1(resnr)
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

      
	WProd448 = Prod(5,2,5,1)
	WProd447 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd445 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd439 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd446 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd433 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd438 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd451 = Prod(5,2,1,1)*ProdC(5,1,1,1)
        WProd449 = Prod(3,1,1,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd452 = Prod(3,1,2,1)*Prod(5,2,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd435 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd453 = Prod(5,2,2,1)*ProdC(5,1,2,1)
        WProd454 = Prod(3,1,2,1)*Prod(5,2,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd450 = -(Prod(4,1,3,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd455 = Prod(4,1,2,1)*Prod(5,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd440 = Prod(5,1,1,1)*ProdC(5,2,1,1)
        WProd436 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,1,1)
        WProd441 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,2,1,1)
	WProd442 = Prod(5,1,2,1)*ProdC(5,2,2,1)
        WProd443 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,2,4,1)
        WProd437 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,2,4,1))
        WProd444 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,4,1)
	WProd434 = ProdC(5,2,5,1)
      
      
	Prp48 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd433)
	Prp49 = -(1.D0/WProd435)
	Prp50 = -(1.D0/WProd438)
      
      
      
      endif

      end


