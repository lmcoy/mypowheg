      subroutine abbrrealQED2(resnr)
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

      
	WProd53 = Prod(5,2,5,1)
	WProd42 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd41 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd50 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd54 = Prod(5,2,1,1)*ProdC(5,1,1,1)
        WProd60 = Prod(3,1,1,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd55 = Prod(3,1,2,1)*Prod(5,2,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd52 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd56 = Prod(5,2,2,1)*ProdC(5,1,2,1)
        WProd57 = Prod(3,1,2,1)*Prod(5,2,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd58 = -(Prod(4,1,3,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd59 = Prod(4,1,2,1)*Prod(5,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd43 = Prod(5,1,1,1)*ProdC(5,2,1,1)
        WProd51 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,1,1)
        WProd44 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,2,1,1)
	WProd45 = Prod(5,1,2,1)*ProdC(5,2,2,1)
        WProd46 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,2,4,1)
        WProd47 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,2,4,1))
        WProd48 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,4,1)
	WProd49 = ProdC(5,2,5,1)
      
      
	Prp11 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd41)
	Prp12 = -(1.D0/WProd50)
	Prp13 = -(1.D0/WProd52)
      
      
      
      endif

      end


