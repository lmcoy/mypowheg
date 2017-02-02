      subroutine abbrborn1(resnr)
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

      
	WProd13 = Prod(5,2,5,1)
	WProd2 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd1 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd10 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd14 = Prod(5,2,1,1)*ProdC(5,1,1,1)
        WProd20 = Prod(3,1,1,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd15 = Prod(3,1,2,1)*Prod(5,2,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd12 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd16 = Prod(5,2,2,1)*ProdC(5,1,2,1)
        WProd17 = Prod(3,1,2,1)*Prod(5,2,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd18 = -(Prod(4,1,3,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd19 = Prod(4,1,2,1)*Prod(5,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd3 = Prod(5,1,1,1)*ProdC(5,2,1,1)
        WProd11 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,1,1)
        WProd4 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,2,1,1)
	WProd5 = Prod(5,1,2,1)*ProdC(5,2,2,1)
        WProd6 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,2,4,1)
        WProd7 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,2,4,1))
        WProd8 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,4,1)
	WProd9 = ProdC(5,2,5,1)
      
      
	Prp1 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd1)
	Prp2 = -(1.D0/WProd10)
	Prp3 = -(1.D0/WProd12)
      
      
      
      endif

      end


