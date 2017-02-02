      subroutine abbrrealQED10(resnr)
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

      
	WProd903 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd904 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd896 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd909 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd911 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd912 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd917 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd913 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(5,1,1,1))
        WProd914 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd918 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd915 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd916 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd910 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd897 = Prod(3,1,2,1)*Prod(5,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd898 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd906 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd899 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(6,1,1,1))
        WProd900 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,2,1)
        WProd907 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,3,1))
        WProd901 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd908 = -(Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd902 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
	WProd905 = Prod(6,1,5,1)*ProdC(6,1,5,1)
      
      
	Prp122 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd896)
	Prp123 = 1.D0/(WProd896 - WProd903 - WProd904)
	Prp124 = 1.D0/WProd905
	Prp125 = 1.D0/(WProd905 - WProd909 - WProd910)
      
      
      
      endif

      end


