      subroutine abbrrealQED1(resnr)
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

      
	WProd32 = Prod(5,2,5,1)
	WProd29 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd18 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd17 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd26 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd33 = Prod(5,2,1,1)*ProdC(5,1,1,1)
        WProd39 = Prod(3,1,1,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd40 = -(Prod(2,1,1,1)*Prod(5,2,3,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,1,1))
        WProd34 = Prod(3,1,2,1)*Prod(5,2,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd28 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd35 = Prod(5,2,2,1)*ProdC(5,1,2,1)
	WProd31 = Prod(5,1,4,1)*ProdC(5,1,4,1)
        WProd36 = Prod(3,1,2,1)*Prod(5,2,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd37 = -(Prod(4,1,3,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd38 = Prod(4,1,2,1)*Prod(5,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd19 = Prod(5,1,1,1)*ProdC(5,2,1,1)
        WProd27 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,1,1)
        WProd30 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(5,2,1,1))
        WProd20 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,2,1,1)
	WProd21 = Prod(5,1,2,1)*ProdC(5,2,2,1)
        WProd22 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,2,4,1)
        WProd23 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,2,4,1))
        WProd24 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,4,1)
	WProd25 = ProdC(5,2,5,1)
      
      
	Prp6 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd17)
	Prp7 = -(1.D0/WProd26)
	Prp8 = -(1.D0/WProd28)
	Prp9 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd29)
	Prp10 = 1.D0/WProd31
      
      
      
      endif

      end


