      subroutine abbrrealQED4(resnr)
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

      
	WProd739 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd732 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd733 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd726 = -(Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd735 = -(Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -     ProdC(4,1,1,1))
        WProd727 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,1,1)
        WProd736 = -(Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd728 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd729 = Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd737 = -(Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd730 = Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
	WProd725 = Prod(4,1,3,1)*ProdC(4,1,3,1)
        WProd731 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,3,1)
	WProd740 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd741 = -(Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,1,1))
        WProd747 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,1,1))
        WProd742 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(5,1,1,1))
	WProd734 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd743 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd738 = -(Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd744 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd745 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd746 = Prod(3,1,2,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd748 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -     ProdC(5,1,4,1))
      
      
	Prp94 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd725)
	Prp95 = 1.D0/(WProd725 - WProd732 - WProd733)
	Prp96 = -(1.D0/WProd734)
	Prp97 = 1.D0/(-WProd734 + WProd739 - WProd740)
      
      
      
      endif

      end


