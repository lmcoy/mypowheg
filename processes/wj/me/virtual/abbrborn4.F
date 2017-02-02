      subroutine abbrborn4(resnr)
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

      
	WProd62 = Prod(2,2,2,1)
	WProd65 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd63 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd77 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd76 = ProdC(2,2,2,1)
	WProd74 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd85 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd71 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd69 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd81 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd66 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd72 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd75 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd67 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd73 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd84 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
        WProd79 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
        WProd68 = -(Prod(3,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,3,1))
        WProd80 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,3,1))
	WProd64 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
        WProd82 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd78 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
        WProd70 = -(Prod(3,1,2,2)*Prod(4,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd83 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,4,1))
      
      
	Prp12 = 1.D0/WProd65
	Prp13 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd66)
	Prp14 = -(1.D0/WProd71)
	Prp15 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r - WProd72)
	Prp16 = -(1.D0/WProd75)
      
      
      
      endif

      end


