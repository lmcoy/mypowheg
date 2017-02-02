      subroutine abbrrealQED25(resnr)
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

      
	WProd1877 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1872 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd1873 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd1865 = Prod(4,1,3,1)*ProdC(4,1,3,1)
        WProd1866 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd1867 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -   ProdC(5,1,2,1)
        WProd1868 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1869 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1875 = Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1870 = Prod(3,1,1,1)*Prod(6,1,4,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1876 = -(Prod(5,1,3,1)*Prod(6,1,1,1)*ProdC(5,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd1871 = Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(5,1,2,1)*
     -   ProdC(5,1,4,1)
	WProd1878 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd1879 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,1,1))
	WProd1874 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd1885 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,3,1))
        WProd1880 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,3,1)
        WProd1886 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,4,1))
        WProd1881 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd1882 = Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,5,1)
        WProd1883 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,3,1)*
     -   ProdC(6,1,5,1)
        WProd1884 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
        WProd1887 = Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
      
      
	Prp211 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1865)
	Prp212 = 1.D0/(WProd1865 - WProd1872 - WProd1873)
	Prp213 = -(1.D0/WProd1874)
	Prp214 = 1.D0/(-WProd1874 + WProd1877 - WProd1878)
      
      
      
      endif

      end


