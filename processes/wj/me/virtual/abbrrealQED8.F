      subroutine abbrrealQED8(resnr)
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

      
	WProd848 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd841 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd842 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd835 = -(Prod(3,1,1,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,1,1))
        WProd844 = -(Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -     ProdC(4,1,1,1))
        WProd836 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,1,1)
        WProd845 = -(Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd837 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd838 = Prod(5,1,4,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
        WProd846 = -(Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -     ProdC(4,1,2,1))
        WProd839 = Prod(5,1,3,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(4,1,2,1)
	WProd834 = Prod(4,1,3,1)*ProdC(4,1,3,1)
        WProd840 = Prod(5,1,3,1)*Prod(6,1,3,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,3,1)
	WProd849 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd869 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd861 = -(Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,1,1))
        WProd870 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd862 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,1,1))
        WProd872 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd866 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(5,1,1,1))
	WProd843 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd871 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd863 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd847 = -(Prod(5,1,3,1)*Prod(6,1,5,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd864 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd867 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd868 = Prod(3,1,2,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd865 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -     ProdC(5,1,4,1))
	WProd860 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd850 = Prod(3,1,2,1)*Prod(5,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd851 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd857 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd852 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(6,1,1,1))
        WProd853 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,2,1)
        WProd858 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,3,1))
        WProd854 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd859 = -(Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd855 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
	WProd856 = Prod(6,1,5,1)*ProdC(6,1,5,1)
      
      
	Prp112 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd834)
	Prp113 = 1.D0/(WProd834 - WProd841 - WProd842)
	Prp114 = -(1.D0/WProd843)
	Prp115 = 1.D0/(-WProd843 + WProd848 - WProd849)
	Prp116 = 1.D0/WProd856
	Prp117 = 1.D0/(-WProd849 + WProd856 - WProd860)
      
      
      
      endif

      end


