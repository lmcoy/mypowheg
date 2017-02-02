      subroutine abbrrealQED21(resnr)
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

      
	WProd1747 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1740 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd1752 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd1749 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,1,1))
        WProd1750 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd1742 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1743 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1744 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1751 = Prod(3,1,2,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1745 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1746 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd1741 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd1748 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd1754 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,2,1)
        WProd1755 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,3,1))
        WProd1759 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,3,1)
        WProd1760 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,3,1)*
     -   ProdC(6,1,4,1)
        WProd1761 = Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd1756 = -(Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,4,1))
        WProd1762 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
	WProd1753 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd1763 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(6,1,5,1))
        WProd1757 = Prod(3,1,2,1)*Prod(5,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
        WProd1758 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,5,1))
      
      
	Prp191 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1740)
	Prp192 = -(1.D0/WProd1741)
	Prp193 = 1.D0/(-WProd1741 + WProd1747 - WProd1748)
	Prp194 = 1.D0/(-WProd1741 - WProd1752 + WProd1753)
      
      
      
      endif

      end


