      subroutine abbrborn5(resnr)
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

      
	WProd86 = Prod(2,2,2,1)
	WProd89 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd87 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd97 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd96 = ProdC(2,2,2,1)
	WProd93 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd101 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd94 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd102 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd90 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd95 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd91 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd92 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd99 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
        WProd100 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
	WProd88 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
        WProd98 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
      
      
	Prp17 = 1.D0/WProd89
	Prp18 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd90)
	Prp19 = -(1.D0/WProd95)
      
      
      
      endif

      end


