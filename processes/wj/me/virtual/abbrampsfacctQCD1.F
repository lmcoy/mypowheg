      subroutine abbrampsfacctQCD1(resnr)
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

      
	WProd513 = Prod(5,2,5,1)
	WProd510 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd499 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd501 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd500 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd498 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd509 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd514 = Prod(5,2,1,1)*ProdC(5,1,1,1)
        WProd520 = Prod(3,1,1,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd515 = Prod(3,1,2,1)*Prod(5,2,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
	WProd512 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd516 = Prod(5,2,2,1)*ProdC(5,1,2,1)
        WProd517 = Prod(3,1,2,1)*Prod(5,2,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd518 = -(Prod(4,1,3,1)*Prod(5,2,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd519 = Prod(4,1,2,1)*Prod(5,2,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd502 = Prod(5,1,1,1)*ProdC(5,2,1,1)
        WProd511 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,1,1)
        WProd503 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,2,1,1)
	WProd504 = Prod(5,1,2,1)*ProdC(5,2,2,1)
        WProd505 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,2,4,1)
        WProd506 = -(Prod(4,1,3,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,2,4,1))
        WProd507 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,2,4,1)
	WProd508 = ProdC(5,2,5,1)
      
      
	Prp57 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd498)
	Prp58 = -(1.D0/WProd509)
	Prp59 = -(1.D0/WProd512)
      
      
      
      endif

      end


