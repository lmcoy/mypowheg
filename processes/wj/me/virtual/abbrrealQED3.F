      subroutine abbrrealQED3(resnr)
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

      
	WProd703 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd698 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd699 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd692 = -(Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd705 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,1,1)
        WProd693 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,1,1)
        WProd701 = -(Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd706 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd694 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd695 = Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd696 = Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
	WProd691 = Prod(4,1,3,1)*ProdC(4,1,3,1)
        WProd697 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,3,1)
	WProd704 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd710 = -(Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,1,1))
        WProd716 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,1,1))
        WProd711 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(5,1,1,1))
	WProd700 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd718 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd719 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd712 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd720 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,3,1)
        WProd702 = -(Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd707 = Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd713 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd721 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd714 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd715 = Prod(3,1,2,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd723 = Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd724 = Prod(3,1,1,1)*Prod(6,1,4,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd717 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd722 = -(Prod(5,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
	WProd708 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd709 = Prod(6,1,5,1)*ProdC(6,1,5,1)
      
      
	Prp88 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd691)
	Prp89 = 1.D0/(WProd691 - WProd698 - WProd699)
	Prp90 = -(1.D0/WProd700)
	Prp91 = 1.D0/(-WProd700 + WProd703 - WProd704)
	Prp92 = -(1.D0/WProd704)
	Prp93 = 1.D0/(-WProd704 - WProd708 + WProd709)
      
      
      
      endif

      end


