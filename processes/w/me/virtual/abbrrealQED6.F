      subroutine abbrrealQED6(resnr)
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "params.h" 

      integer resnr,delta
      double complex myA0,myB0,myB1,myB00,myB11,myB001,myB111,myDB0
      double complex myDB1,myDB00,myDB11,myBval,myCval,myDval,myEval
      double complex Prp,myDeltaIR,myDeltaUV

      integer myBget,myCget,myDget,myEget

#include "looptools.h" 
#include "amps/ampcommon.h" 
#include "amps/ampvars.h" 
#include "amps/ampcom.h" 
      if(resnr.eq.0)then

      
	WProd126 = Prod(2,2,2,1)
	WProd129 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd127 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd137 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd136 = ProdC(2,2,2,1)
	WProd133 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd141 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd134 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd142 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd130 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd135 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd131 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd132 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd139 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
        WProd140 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
	WProd128 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
        WProd138 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
      
      
	Prp27 = 1.D0/WProd129
	Prp28 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd130)
	Prp29 = -(1.D0/WProd135)
      
      
      
      endif

      end


