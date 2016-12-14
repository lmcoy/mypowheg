      subroutine abbrvirtEWphotonicISR2(resnr)
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

      
	WProd27 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd28 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd30 = Prod(4,1,1,1)*ProdC(3,1,2,1)
	WProd26 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd29 = Prod(3,1,1,1)*ProdC(4,1,2,1)
      
      
	Prp11 = 1.D0/WProd27
	Prp12 = 1.D0/((0.D0,1.D0)*GZ*Sqrt(MZ2r) - MZ2r + WProd27)
      
	dr=1D0
	If(dr.eq.0D0) then
	IInt90 = myBget(WProd27,MD2,MD2,0)
	LInt86 = myBval(bb0,IInt90)
	IInt91 = myCget(MD2,WProd27,MD2,0.D0,MD2,MD2,2)
	LInt87 = myCval(cc0,IInt91)
	LInt88 = myCval(cc00,IInt91)
	LInt89 = myCval(cc1,IInt91)
	LInt90 = myCval(cc2,IInt91)
	Else
	IInt90 = myBget(WProd27,0D0,0D0,0)
	LInt86 = myBval(bb0,IInt90)
	IInt91 = myCget(0D0,WProd27,0D0,0.D0,0D0,0D0,2)
	LInt87 = myCval(cc0,IInt91)
	LInt88 = myCval(cc00,IInt91)
	LInt89 = myCval(cc1,IInt91)
	LInt90 = myCval(cc2,IInt91)
	Endif
      
      
      endif

      end


