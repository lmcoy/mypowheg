      subroutine abbrrealQED22(resnr)
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

      
	WProd1771 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1786 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1787 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd1764 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd1778 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd1773 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,1,1))
        WProd1780 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd1781 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -   ProdC(5,1,2,1)
        WProd1774 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd1782 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1766 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1767 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1775 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1788 = -(Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1783 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1768 = -(Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1776 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1777 = Prod(3,1,2,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1789 = -(Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1784 = Prod(3,1,1,1)*Prod(6,1,4,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1769 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1770 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1790 = -(Prod(5,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1785 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(5,1,2,1)*
     -   ProdC(5,1,4,1)
	WProd1765 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd1791 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,1,1))
	WProd1772 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd1799 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,2,1)
        WProd1792 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,3,1)
        WProd1802 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,3,1)*
     -   ProdC(6,1,4,1)
        WProd1797 = Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd1793 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
	WProd1779 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd1794 = Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd1803 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(6,1,5,1))
        WProd1795 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,3,1)*
     -   ProdC(6,1,5,1)
        WProd1800 = Prod(3,1,2,1)*Prod(5,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
        WProd1796 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
        WProd1801 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,5,1))
        WProd1798 = Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
      
      
	Prp195 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1764)
	Prp196 = -(1.D0/WProd1765)
	Prp197 = 1.D0/(-WProd1765 + WProd1771 - WProd1772)
	Prp199 = -(1.D0/WProd1772)
	Prp198 = 1.D0/(-WProd1765 - WProd1778 + WProd1779)
	Prp200 = 1.D0/(WProd1764 - WProd1786 - WProd1787)
      
      
      
      endif

      end


