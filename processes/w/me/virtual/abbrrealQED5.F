      subroutine abbrrealQED5(resnr)
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

      
	WProd102 = Prod(2,2,2,1)
	WProd105 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd103 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd117 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd116 = ProdC(2,2,2,1)
	WProd114 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd125 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd111 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd109 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd121 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd106 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd112 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd115 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd107 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd113 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd124 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
        WProd119 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
        WProd108 = -(Prod(3,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,3,1))
        WProd120 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,3,1))
	WProd104 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
        WProd122 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd118 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
        WProd110 = -(Prod(3,1,2,2)*Prod(4,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd123 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,4,1))
      
      
	Prp22 = 1.D0/WProd105
	Prp23 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd106)
	Prp24 = -(1.D0/WProd111)
	Prp25 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r - WProd112)
	Prp26 = -(1.D0/WProd115)
      
      
      
      endif

      end


