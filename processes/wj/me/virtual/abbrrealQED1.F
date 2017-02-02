      subroutine abbrrealQED1(resnr)
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

      
	WProd626 = Prod(2,1,1,1)*ProdC(2,1,1,1)
        WProd619 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,1,1)
        WProd620 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,1,1)
        WProd621 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd622 = -(Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd628 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd629 = Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd623 = -(Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd624 = Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
	WProd617 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd618 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd632 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd633 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd641 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd634 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd627 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd649 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd642 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd635 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd650 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd643 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd644 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd651 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,3,1)
        WProd645 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd636 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd625 = Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd638 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd647 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd639 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd648 = Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd653 = Prod(3,1,1,1)*Prod(6,1,4,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd637 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd652 = -(Prod(5,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd646 = Prod(5,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,1)*
     -   ProdC(5,1,4,1)
	WProd630 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd640 = Prod(6,1,2,1)*ProdC(6,1,2,1)
	WProd631 = Prod(6,1,5,1)*ProdC(6,1,5,1)
      
      
	Prp76 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd617)
	Prp77 = -(1.D0/WProd618)
	Prp80 = 1.D0/WProd626
	Prp78 = 1.D0/(-WProd618 + WProd626 - WProd627)
	Prp79 = 1.D0/(-WProd618 - WProd630 + WProd631)
	Prp81 = 1.D0/(WProd626 - WProd630 - WProd640)
      
      
      
      endif

      end


