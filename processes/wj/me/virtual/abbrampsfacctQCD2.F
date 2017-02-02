      subroutine abbrampsfacctQCD2(resnr)
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

      
	WProd521 = Prod(2,2,2,1)
	WProd531 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd522 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd539 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd538 = ProdC(2,2,2,1)
	WProd526 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd527 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd536 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd547 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd528 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd523 = -(Prod(5,1,3,1)*ProdC(4,1,1,1))
        WProd543 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -     ProdC(4,1,1,1))
        WProd540 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(4,1,1,1)
	WProd529 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd532 = Prod(4,1,2,2)*ProdC(4,1,2,1)
        WProd524 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd525 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd541 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,2,1))
        WProd533 = -(Prod(3,1,2,2)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd534 = -(Prod(3,1,2,2)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
	WProd544 = Prod(4,1,2,1)*ProdC(4,1,2,2)
        WProd542 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd545 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,2))
        WProd546 = -(Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,2))
	WProd530 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd535 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd537 = Prod(5,1,2,1)*ProdC(5,1,2,1)
      
      
	Prp62 = -(1.D0/WProd529)
	Prp60 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd530)
	Prp61 = 1.D0/WProd531
	Prp63 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r - WProd535)
	Prp64 = -(1.D0/WProd537)
      
      
      
      endif

      end


