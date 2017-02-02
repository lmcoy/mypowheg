      subroutine abbrrealQED7(resnr)
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

      
	WProd818 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd819 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd811 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd824 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd826 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd827 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd832 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd828 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(5,1,1,1))
        WProd829 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd833 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd830 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd831 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd825 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd812 = Prod(3,1,2,1)*Prod(5,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd813 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd821 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd814 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(6,1,1,1))
        WProd815 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,2,1)
        WProd822 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,3,1))
        WProd816 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd823 = -(Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd817 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
	WProd820 = Prod(6,1,5,1)*ProdC(6,1,5,1)
      
      
	Prp108 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd811)
	Prp109 = 1.D0/(WProd811 - WProd818 - WProd819)
	Prp110 = 1.D0/WProd820
	Prp111 = 1.D0/(WProd820 - WProd824 - WProd825)
      
      
      
      endif

      end


