      subroutine abbrvirtEWphotonicFSR2(resnr)
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

      
	WProd39 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd40 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd42 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd41 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd44 = Prod(4,1,1,1)*ProdC(3,1,2,1)
	WProd38 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd43 = Prod(3,1,1,1)*ProdC(4,1,2,1)
      
      
	Prp15 = 1.D0/WProd39
	Prp16 = 1.D0/((0.D0,1.D0)*GZ*Sqrt(MZ2r) - MZ2r + WProd39)
      
	dr=1D0
	If(dr.eq.0D0) then
	IInt80 = myBget(WProd39,ME2,ME2,0)
	LInt96 = myBval(bb0,IInt80)
	IInt79 = myCget(ME2,WProd39,ME2,0.D0,ME2,ME2,2)
	LInt100 = myCval(cc2,IInt79)
	LInt97 = myCval(cc0,IInt79)
	LInt98 = myCval(cc00,IInt79)
	LInt99 = myCval(cc1,IInt79)
	Else
	IInt80 = myBget(WProd39,0D0,0D0,0)
	LInt96 = myBval(bb0,IInt80)
	IInt79 = myCget(0D0,WProd39,0D0,0.D0,0D0,0D0,2)
	LInt100 = myCval(cc2,IInt79)
	LInt97 = myCval(cc0,IInt79)
	LInt98 = myCval(cc00,IInt79)
	LInt99 = myCval(cc1,IInt79)
	Endif      
      
      endif

      end


