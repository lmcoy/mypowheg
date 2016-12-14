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

#include "amps/looptools.h" 
#include "amps/ampcommon.h"
#include "amps/ampvars.h" 
#include "amps/ampcom.h" 
      if(resnr.eq.0)then

      
	WProd12 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd13 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd15 = Prod(4,1,1,1)*ProdC(3,1,2,1)
	WProd11 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd14 = Prod(3,1,1,1)*ProdC(4,1,2,1)
      
      
	Prp5 = 1.D0/WProd12
	Prp6 = 1.D0/((0.D0,1.D0)*GZ*Sqrt(MZ2r) - MZ2r + WProd12)
      
	IInt57 = myBget(WProd12,0.D0,0.D0,0)
	LInt1 = myBval(bb0,IInt57)
	IInt58 = myCget(0.D0,WProd12,0.D0,0.D0,0.D0,0.D0,2)
	LInt2 = myCval(cc0,IInt58)
	LInt3 = myCval(cc00,IInt58)
	LInt4 = myCval(cc1,IInt58)
	LInt5 = myCval(cc2,IInt58)   
      endif

      end


