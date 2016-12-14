      subroutine abbrvirtEWnonphotonic2(resnr)
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

      
	WProd100 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd101 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd107 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd108 = Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(3,1,1,1)**2
	WProd103 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd113 = Prod(4,1,1,1)*ProdC(3,1,2,1)
        WProd114 = Prod(2,1,1,1)*Prod(4,1,1,1)*ProdC(2,1,1,1)*
     -   ProdC(3,1,2,1)
        WProd115 = Prod(3,1,1,1)*Prod(4,1,1,1)*ProdC(3,1,1,1)*
     -   ProdC(3,1,2,1)
        WProd109 = Prod(3,1,2,1)*Prod(4,1,2,1)*ProdC(3,1,1,1)*
     -   ProdC(3,1,2,1)
	WProd116 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(3,1,2,1)**2
	WProd99 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd102 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd104 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,1,1)
        WProd105 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,1,1)
	WProd111 = Prod(3,1,1,1)**2*ProdC(3,1,2,1)*ProdC(4,1,1,1)
	WProd110 = Prod(3,1,1,1)*ProdC(4,1,2,1)
	WProd106 = Prod(3,1,2,1)**2*ProdC(3,1,1,1)*ProdC(4,1,2,1)
        WProd112 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,2,1)
      
      
	Prp19 = 1.D0/WProd100
	Prp20 = 1.D0/((0.D0,1.D0)*GZ*Sqrt(MZ2r) - MZ2r + WProd100)
      
	IInt59 = myCget(WProd100,ME2,ME2,MZ2r,MZ2r,ME2,0)
	LInt283 = myCval(cc0,IInt59)
	IInt63 = myBget(WProd100,MD2,MD2,0)
	LInt306 = myBval(bb0,IInt63)
	IInt64 = myBget(WProd100,ME2,ME2,0)
	LInt307 = myBval(bb0,IInt64)
	IInt71 = myBget(WProd100,MU2,MU2,0)
	LInt334 = myBval(bb0,IInt71)
	IInt74 = myCget(WProd100,ME2,ME2,MW2r,MW2r,0.D0,0)
	LInt342 = myCval(cc0,IInt74)
      

	do Gen1=1,3
	IInt77(Gen1) = myBget(WProd100,MLE2(Gen1),MLE2(Gen1),1)
	LInt323(Gen1) = myBval(bb1,IInt77(Gen1))
	IInt78(Gen1) = myBget(WProd100,MQD2(Gen1),MQD2(Gen1),1)
	LInt324(Gen1) = myBval(bb1,IInt78(Gen1))
        enddo

	IInt67 = myBget(WProd100,MH2,MZ2r,2)
	LInt318 = myBval(bb0,IInt67)
	LInt319 = myBval(bb00,IInt67)
	IInt68 = myBget(WProd100,0.D0,0.D0,1)
	LInt320 = myBval(bb1,IInt68)
	LInt326 = myBval(bb0,IInt68)
      

	do Gen1=1,3
	IInt76(Gen1) = myBget(WProd100,MQU2(Gen1),MQU2(Gen1),1)
	LInt322(Gen1) = myBval(bb0,IInt76(Gen1))
	LInt325(Gen1) = myBval(bb1,IInt76(Gen1))
        enddo

	IInt62 = myBget(WProd100,MW2r,MW2r,2)
	LInt303 = myBval(bb0,IInt62)
	LInt304 = myBval(bb00,IInt62)
	LInt305 = myBval(bb1,IInt62)
	IInt70 = myCget(ME2,WProd100,ME2,0.D0,MW2r,MW2r,2)
	LInt328 = myCval(cc00,IInt70)
	LInt330 = myCval(cc1,IInt70)
	LInt332 = myCval(cc2,IInt70)
	IInt73 = myCget(MD2,WProd100,MD2,MU2,MW2r,MW2r,2)
	LInt336 = myCval(cc00,IInt73)
	LInt338 = myCval(cc1,IInt73)
	LInt340 = myCval(cc2,IInt73)
	IInt65 = myCget(MD2,WProd100,MD2,MZ2r,MD2,MD2,2)
	LInt308 = myCval(cc0,IInt65)
	LInt310 = myCval(cc00,IInt65)
	LInt312 = myCval(cc1,IInt65)
	LInt314 = myCval(cc2,IInt65)
	IInt66 = myCget(ME2,WProd100,ME2,MZ2r,ME2,ME2,2)
	LInt309 = myCval(cc0,IInt66)
	LInt311 = myCval(cc00,IInt66)
	LInt313 = myCval(cc1,IInt66)
	LInt315 = myCval(cc2,IInt66)
	IInt69 = myCget(ME2,WProd100,ME2,MW2r,0.D0,0.D0,2)
	LInt327 = myCval(cc0,IInt69)
	LInt329 = myCval(cc00,IInt69)
	LInt331 = myCval(cc1,IInt69)
	LInt333 = myCval(cc2,IInt69)
	IInt72 = myCget(MD2,WProd100,MD2,MW2r,MU2,MU2,2)
	LInt335 = myCval(cc0,IInt72)
	LInt337 = myCval(cc00,IInt72)
	LInt339 = myCval(cc1,IInt72)
	LInt341 = myCval(cc2,IInt72)
        IInt60 = myDget(MD2,WProd100,ME2,-WProd101,MD2,ME2,MD2,
     -   MZ2r,MZ2r,ME2,2)
	LInt284 = myDval(dd00,IInt60)
	LInt286 = myDval(dd1,IInt60)
	LInt288 = myDval(dd11,IInt60)
	LInt290 = myDval(dd12,IInt60)
	LInt292 = myDval(dd13,IInt60)
	LInt294 = myDval(dd2,IInt60)
	LInt296 = myDval(dd23,IInt60)
	LInt298 = myDval(dd3,IInt60)
	LInt300 = myDval(dd33,IInt60)
        IInt61 = myDget(MD2,WProd100,ME2,-WProd102,MD2,ME2,MD2,
     -   MZ2r,MZ2r,ME2,2)
	LInt285 = myDval(dd00,IInt61)
	LInt287 = myDval(dd1,IInt61)
	LInt289 = myDval(dd11,IInt61)
	LInt291 = myDval(dd12,IInt61)
	LInt293 = myDval(dd13,IInt61)
	LInt295 = myDval(dd2,IInt61)
	LInt297 = myDval(dd23,IInt61)
	LInt299 = myDval(dd3,IInt61)
	LInt301 = myDval(dd33,IInt61)
        IInt75 = myDget(MD2,WProd100,ME2,-WProd101,MD2,ME2,MU2,
     -   MW2r,MW2r,0.D0,2)
	LInt343 = myDval(dd00,IInt75)
	LInt344 = myDval(dd1,IInt75)
	LInt345 = myDval(dd11,IInt75)
	LInt346 = myDval(dd12,IInt75)
	LInt347 = myDval(dd13,IInt75)
	LInt348 = myDval(dd2,IInt75)
	LInt349 = myDval(dd23,IInt75)
	LInt350 = myDval(dd3,IInt75)
	LInt351 = myDval(dd33,IInt75)
      
      
      endif

      end


