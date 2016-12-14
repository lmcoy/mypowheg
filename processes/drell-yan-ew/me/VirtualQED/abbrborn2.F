      subroutine abbrborn2(resnr)
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

      
	WProd7 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd8 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd10 = Prod(4,1,1,1)*ProdC(3,1,2,1)
	WProd6 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd9 = Prod(3,1,1,1)*ProdC(4,1,2,1)
      
      
	Prp3 = 1.D0/WProd7
	Prp4 = 1.D0/((0.D0,1.D0)*GZ*Sqrt(MZ2r) - MZ2r + WProd7)
      
      
      
      endif

      end


