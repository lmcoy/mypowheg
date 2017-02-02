      subroutine abbrrealQED6(resnr)
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

      
	WProd797 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd778 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd779 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd771 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd784 = Prod(5,1,1,1)*ProdC(5,1,1,1)
        WProd786 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd787 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd799 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,1,1))
        WProd792 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd788 = -(Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(5,1,1,1))
        WProd789 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd793 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,3,1)
        WProd794 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd790 = Prod(3,1,2,1)*Prod(6,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd791 = Prod(4,1,2,1)*Prod(6,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd800 = Prod(3,1,2,1)*Prod(6,1,4,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd795 = -(Prod(5,1,3,1)*Prod(6,1,2,1)*ProdC(5,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd796 = Prod(3,1,2,1)*Prod(6,1,5,1)*ProdC(5,1,1,1)*
     -   ProdC(5,1,4,1)
	WProd785 = Prod(6,1,1,1)*ProdC(6,1,1,1)
        WProd772 = Prod(3,1,2,1)*Prod(5,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd773 = Prod(3,1,1,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,1,1)
        WProd781 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(6,1,1,1)
        WProd774 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,3,1)*
     -     ProdC(6,1,1,1))
	WProd798 = Prod(6,1,2,1)*ProdC(6,1,2,1)
        WProd775 = Prod(3,1,2,1)*Prod(5,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,2,1)
        WProd801 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,2,1)
        WProd782 = -(Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,3,1))
        WProd802 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,3,1))
        WProd806 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,3,1)
        WProd776 = Prod(3,1,2,1)*Prod(5,1,3,1)*ProdC(3,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd783 = -(Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -     ProdC(6,1,4,1))
        WProd777 = Prod(4,1,2,1)*Prod(5,1,3,1)*ProdC(4,1,1,1)*
     -   ProdC(6,1,4,1)
        WProd807 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(5,1,3,1)*
     -   ProdC(6,1,4,1)
        WProd808 = Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
        WProd803 = -(Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,4,1))
        WProd809 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,4,1)
	WProd780 = Prod(6,1,5,1)*ProdC(6,1,5,1)
        WProd810 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(6,1,5,1))
        WProd804 = Prod(3,1,2,1)*Prod(5,1,1,1)*ProdC(5,1,4,1)*
     -   ProdC(6,1,5,1)
        WProd805 = -(Prod(2,1,1,1)*Prod(5,1,3,1)*ProdC(5,1,4,1)*
     -     ProdC(6,1,5,1))
      
      
	Prp102 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd771)
	Prp103 = 1.D0/(WProd771 - WProd778 - WProd779)
	Prp104 = 1.D0/WProd780
	Prp105 = 1.D0/(WProd780 - WProd784 - WProd785)
	Prp106 = -(1.D0/WProd785)
	Prp107 = 1.D0/(-WProd785 + WProd797 - WProd798)
      
      
      
      endif

      end


