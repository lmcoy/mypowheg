      subroutine abbrrealQED2(resnr)
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

      
	WProd663 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd654 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd655 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd675 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,1,1))
	WProd664 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd676 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd669 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd670 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd677 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd671 = -(Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd678 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd679 = Prod(3,1,2,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd672 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd673 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd667 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd656 = -(Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,1,1))
        WProd657 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(6,1,1,1))
        WProd658 = -(Prod(5,1,2,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,1,1))
	WProd674 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd687 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,2,1)
        WProd665 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,3,1)
        WProd685 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,3,1)
        WProd659 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd660 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd661 = -(Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd666 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd662 = Prod(3,1,2,1)*Prod(5,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd686 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,3,1)*
     -   ProdC(6,1,4,1)
        WProd680 = Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd681 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
	WProd668 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd682 = Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd690 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(6,1,5,1))
        WProd683 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,5,1)
        WProd688 = Prod(3,1,2,1)*Prod(5,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
        WProd689 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,5,1))
        WProd684 = Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
      
      
	Prp82 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd654)
	Prp83 = -(1.D0/WProd655)
	Prp84 = 1.D0/(-WProd655 + WProd663 - WProd664)
	Prp86 = -(1.D0/WProd667)
	Prp85 = 1.D0/(-WProd655 - WProd667 + WProd668)
	Prp87 = 1.D0/(WProd663 - WProd667 - WProd674)
      
      
      
      endif

      end


