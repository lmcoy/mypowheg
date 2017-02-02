      subroutine abbrborn3(resnr)
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

      
	WProd45 = Prod(2,2,2,1)
	WProd50 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd46 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd56 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd55 = ProdC(2,2,2,1)
	WProd52 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd60 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd47 = -(Prod(5,1,3,1)*ProdC(4,1,1,1))
        WProd57 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(4,1,1,1)
	WProd53 = Prod(4,1,2,2)*ProdC(4,1,2,1)
        WProd48 = -(Prod(3,1,2,2)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd49 = -(Prod(3,1,2,1)*Prod(5,1,2,2)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd58 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,2,1,1)*
     -     ProdC(4,1,2,1))
	WProd61 = Prod(4,1,2,1)*ProdC(4,1,2,2)
        WProd59 = -(Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,2))
	WProd51 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd54 = Prod(5,1,2,1)*ProdC(5,1,2,1)
      
      
	Prp9 = 1.D0/WProd50
	Prp10 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd51)
	Prp11 = -(1.D0/WProd54)
      
      
      
      endif

      end


