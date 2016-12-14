      subroutine abbrvirtQCD1(resnr)
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

      
	WProd4 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd3 = Prod(3,1,2,1)*ProdC(4,1,1,1)
      
      
	Prp2 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd4)
      
	IInt53 = myBget(WProd4,0.D0,0.D0,0)
	LInt1 = myBval(bb0,IInt53)
	IInt54 = myCget(0.D0,WProd4,0.D0,0.D0,0.D0,0.D0,2)
	LInt2 = myCval(cc0,IInt54)
	LInt3 = myCval(cc00,IInt54)
	LInt4 = myCval(cc1,IInt54)
	LInt5 = myCval(cc2,IInt54)
      
      
      endif

      end


