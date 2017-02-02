      subroutine abbrrealQED20(resnr)
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

      
	WProd1723 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd1716 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd1730 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd1725 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,1,1))
        WProd1726 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd1718 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1719 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1727 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1720 = -(Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1728 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1729 = Prod(3,1,2,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd1721 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd1722 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd1717 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd1724 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd1732 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,2,1)
        WProd1733 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,3,1)
        WProd1737 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,3,1)*
     -   ProdC(6,1,4,1)
        WProd1738 = Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd1734 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
	WProd1731 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd1739 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(6,1,5,1))
        WProd1735 = Prod(3,1,2,1)*Prod(5,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
        WProd1736 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,5,1))
      
      
	Prp187 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1716)
	Prp188 = -(1.D0/WProd1717)
	Prp189 = 1.D0/(-WProd1717 + WProd1723 - WProd1724)
	Prp190 = 1.D0/(-WProd1717 - WProd1730 + WProd1731)
      
      
      
      endif

      end


