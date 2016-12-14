      subroutine abbrvirtEWphotonicFSR1(resnr)
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

      
	WProd32 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd33 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd35 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd34 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd37 = Prod(4,1,1,1)*ProdC(3,1,2,1)
	WProd31 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd36 = Prod(3,1,1,1)*ProdC(4,1,2,1)
      
      
	Prp13 = 1.D0/WProd32
	Prp14 = 1.D0/((0.D0,1.D0)*GZ*Sqrt(MZ2r) - MZ2r + WProd32)
      
	dr=1D0
	If(dr.eq.0D0) then
	IInt44 = myBget(WProd32,ME2,ME2,0)
	LInt91 = myBval(bb0,IInt44)
	IInt45 = myCget(ME2,WProd32,ME2,0.D0,ME2,ME2,2)
	LInt92 = myCval(cc0,IInt45)
	LInt93 = myCval(cc00,IInt45)
	LInt94 = myCval(cc1,IInt45)
	LInt95 = myCval(cc2,IInt45)
	Else
	IInt44 = myBget(WProd32,0D0,0D0,0)
	LInt91 = myBval(bb0,IInt44)
	IInt45 = myCget(0D0,WProd32,0D0,0.D0,0D0,0D0,2)
	LInt92 = myCval(cc0,IInt45)
	LInt93 = myCval(cc00,IInt45)
	LInt94 = myCval(cc1,IInt45)
	LInt95 = myCval(cc2,IInt45)
	Endif
      
      
      endif

      end


