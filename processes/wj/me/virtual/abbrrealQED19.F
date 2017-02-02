      subroutine abbrrealQED19(resnr)
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

      
	WProd1694 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1695 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd1702 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd1696 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd1697 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd1708 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd1698 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd1703 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd1709 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd1699 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd1710 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd1711 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd1712 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1700 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1704 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1714 = Prod(3,1,1,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1705 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1715 = Prod(4,1,1,1)*Prod(6,1,3,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,4,1)
        WProd1701 = Prod(5,1,2,1)*Prod(6,1,3,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1713 = Prod(5,1,1,1)*Prod(6,1,3,1)*ProdC(5,1,2,1)*
     -   ProdC(5,1,4,1)
	WProd1706 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd1707 = Prod(6,1,2,1)*ProdC(6,1,2,1)
      
      
	Prp183 = 1.D0/WProd1694
	Prp184 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1695)
	Prp185 = 1.D0/(WProd1694 - WProd1702 - WProd1703)
	Prp186 = 1.D0/(WProd1694 - WProd1706 - WProd1707)
      
      
      
      endif

      end


