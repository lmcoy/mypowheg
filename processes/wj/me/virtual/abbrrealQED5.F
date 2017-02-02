      subroutine abbrrealQED5(resnr)
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

      
	WProd761 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd756 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd757 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd750 = -(Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd751 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,1,1)
        WProd759 = -(Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd752 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd753 = Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd754 = Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
	WProd749 = Prod(4,1,3,1)*ProdC(4,1,3,1)
        WProd755 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,3,1)
	WProd762 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd763 = -(Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,1,1))
        WProd769 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,1,1))
        WProd764 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(5,1,1,1))
	WProd758 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd765 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd760 = -(Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd766 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd767 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd768 = Prod(3,1,2,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd770 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -     ProdC(5,1,4,1))
      
      
	Prp98 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd749)
	Prp99 = 1.D0/(WProd749 - WProd756 - WProd757)
	Prp100 = -(1.D0/WProd758)
	Prp101 = 1.D0/(-WProd758 + WProd761 - WProd762)
      
      
      
      endif

      end


