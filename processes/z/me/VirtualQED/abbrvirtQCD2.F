      subroutine abbrvirtQCD2(resnr)
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

      
	WProd17 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd18 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd20 = Prod(4,1,1,1)*ProdC(3,1,2,1)
	WProd16 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd19 = Prod(3,1,1,1)*ProdC(4,1,2,1)
      
      
	Prp7 = 1.D0/WProd17
	Prp8 = 1.D0/((0.D0,1.D0)*GZ*Sqrt(MZ2r) - MZ2r + WProd17)
      
	IInt93 = myBget(WProd17,0.D0,0.D0,0)
	LInt6 = myBval(bb0,IInt93)
	IInt92 = myCget(0.D0,WProd17,0.D0,0.D0,0.D0,0.D0,2)
	LInt10 = myCval(cc2,IInt92)
	LInt7 = myCval(cc0,IInt92)
	LInt8 = myCval(cc00,IInt92)
	LInt9 = myCval(cc1,IInt92)
      
      endif

      end


