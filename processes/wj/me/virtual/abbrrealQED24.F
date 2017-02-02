      subroutine abbrrealQED24(resnr)
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

      
	WProd1855 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1849 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1850 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd1842 = Prod(4,1,3,1)*ProdC(4,1,3,1)
        WProd1843 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd1844 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -   ProdC(5,1,2,1)
        WProd1845 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1852 = -(Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1846 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1853 = -(Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1847 = Prod(3,1,1,1)*Prod(6,1,4,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1854 = -(Prod(5,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1848 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(5,1,2,1)*
     -   ProdC(5,1,4,1)
	WProd1856 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd1857 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,1,1))
	WProd1851 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd1858 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,3,1)
        WProd1863 = Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd1859 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd1860 = Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd1861 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,3,1)*
     -   ProdC(6,1,5,1)
        WProd1862 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
        WProd1864 = Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
      
      
	Prp207 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1842)
	Prp208 = 1.D0/(WProd1842 - WProd1849 - WProd1850)
	Prp209 = -(1.D0/WProd1851)
	Prp210 = 1.D0/(-WProd1851 + WProd1855 - WProd1856)
      
      
      
      endif

      end


