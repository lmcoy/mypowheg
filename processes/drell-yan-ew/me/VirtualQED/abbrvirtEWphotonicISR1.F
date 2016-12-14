      subroutine abbrvirtEWphotonicISR1(resnr)
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)


#include "params.h" 

      integer resnr,delta
      double complex myA0,myB0,myB1,myB00,myB11,myB001,myB111,myDB0
      double complex myDB1,myDB00,myDB11,myBval,myCval,myDval,myEval
      double complex Prp,myDeltaIR,myDeltaUV
      double precision dr
      
      integer myBget,myCget,myDget,myEget

#include "amps/looptools.h" 
#include "amps/ampcommon.h"
#include "amps/ampvars.h" 
#include "amps/ampcom.h" 
      if(resnr.eq.0)then

      
	WProd22 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd23 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd25 = Prod(4,1,1,1)*ProdC(3,1,2,1)
	WProd21 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd24 = Prod(3,1,1,1)*ProdC(4,1,2,1)
      
      
	Prp9 = 1.D0/WProd22
	Prp10 = 1.D0/((0.D0,1.D0)*GZ*Sqrt(MZ2r) - MZ2r + WProd22)
      
	dr=1D0
	If(dr.eq.0D0) then
	IInt55 = myBget(WProd22,MU2,MU2,0)
	LInt81 = myBval(bb0,IInt55)
	IInt56 = myCget(MU2,WProd22,MU2,0.D0,MU2,MU2,2)
	LInt82 = myCval(cc0,IInt56)
	LInt83 = myCval(cc00,IInt56)
	LInt84 = myCval(cc1,IInt56)
	LInt85 = myCval(cc2,IInt56)
	Else
	IInt55 = myBget(WProd22,0D0,0D0,0)
	LInt81 = myBval(bb0,IInt55)
	IInt56 = myCget(0D0,WProd22,0D0,0.D0,0D0,0D0,2)
	LInt82 = myCval(cc0,IInt56)
	LInt83 = myCval(cc00,IInt56)
	LInt84 = myCval(cc1,IInt56)
	LInt85 = myCval(cc2,IInt56)
	Endif
      
      
      endif

      end


