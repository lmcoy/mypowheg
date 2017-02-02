      subroutine abbrampsfacctQCD4(resnr)
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

      
	WProd569 = Prod(2,2,2,1)
	WProd577 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd570 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd587 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd586 = ProdC(2,2,2,1)
	WProd571 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd572 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd584 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd595 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd573 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd574 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd580 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd591 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd575 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd582 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd585 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd578 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd583 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd594 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
        WProd589 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
        WProd579 = -(Prod(3,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,3,1))
        WProd590 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,3,1))
	WProd576 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
        WProd592 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd588 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
        WProd581 = -(Prod(3,1,2,2)*Prod(4,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd593 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,4,1))
      
      
	Prp70 = -(1.D0/WProd574)
	Prp68 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd575)
	Prp69 = 1.D0/WProd577
	Prp71 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r - WProd582)
	Prp72 = -(1.D0/WProd585)
      
      
      
      endif

      end


