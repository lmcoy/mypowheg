      subroutine abbrrealQED23(resnr)
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

      
	WProd1804 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1830 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1831 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd1805 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd1812 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd1806 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd1807 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd1818 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd1808 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd1813 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd1826 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd1819 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd1809 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd1820 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd1821 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd1827 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -   ProdC(5,1,2,1)
        WProd1822 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1810 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1814 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1824 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1815 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1825 = Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1828 = Prod(3,1,1,1)*Prod(6,1,4,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1811 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1832 = -(Prod(5,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1823 = Prod(5,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1829 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(5,1,2,1)*
     -   ProdC(5,1,4,1)
	WProd1816 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd1833 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,1,1))
	WProd1817 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd1839 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,3,1))
        WProd1834 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,3,1)
        WProd1840 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,4,1))
        WProd1835 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd1836 = Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd1837 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,3,1)*
     -   ProdC(6,1,5,1)
        WProd1838 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
        WProd1841 = Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
      
      
	Prp201 = 1.D0/WProd1804
	Prp202 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1805)
	Prp203 = 1.D0/(WProd1804 - WProd1812 - WProd1813)
	Prp204 = 1.D0/(WProd1804 - WProd1816 - WProd1817)
	Prp205 = -(1.D0/WProd1817)
	Prp206 = 1.D0/(WProd1805 - WProd1830 - WProd1831)
      
      
      
      endif

      end


