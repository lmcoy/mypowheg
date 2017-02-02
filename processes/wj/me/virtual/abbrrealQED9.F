      subroutine abbrrealQED9(resnr)
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

      
	WProd880 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd881 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd873 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd886 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd888 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd889 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd894 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd890 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(5,1,1,1))
        WProd891 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd895 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd892 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd893 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd887 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd874 = Prod(3,1,2,1)*Prod(5,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd875 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd883 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd876 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(6,1,1,1))
        WProd877 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,2,1)
        WProd884 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,3,1))
        WProd878 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd885 = -(Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd879 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
	WProd882 = Prod(6,1,5,1)*ProdC(6,1,5,1)
      
      
	Prp118 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd873)
	Prp119 = 1.D0/(WProd873 - WProd880 - WProd881)
	Prp120 = 1.D0/WProd882
	Prp121 = 1.D0/(WProd882 - WProd886 - WProd887)
      
      
      
      endif

      end


