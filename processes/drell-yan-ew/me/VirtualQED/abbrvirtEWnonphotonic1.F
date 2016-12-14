      subroutine abbrvirtEWnonphotonic1(resnr)
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

      
	WProd82 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd83 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd89 = Prod(4,1,2,1)*ProdC(3,1,1,1)
	WProd90 = Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(3,1,1,1)**2
	WProd85 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd95 = Prod(4,1,1,1)*ProdC(3,1,2,1)
        WProd96 = Prod(2,1,1,1)*Prod(4,1,1,1)*ProdC(2,1,1,1)*
     -   ProdC(3,1,2,1)
        WProd97 = Prod(3,1,1,1)*Prod(4,1,1,1)*ProdC(3,1,1,1)*
     -   ProdC(3,1,2,1)
        WProd91 = Prod(3,1,2,1)*Prod(4,1,2,1)*ProdC(3,1,1,1)*
     -   ProdC(3,1,2,1)
	WProd98 = Prod(3,1,1,1)*Prod(4,1,2,1)*ProdC(3,1,2,1)**2
	WProd81 = Prod(3,1,2,1)*ProdC(4,1,1,1)
	WProd84 = Prod(4,1,1,1)*ProdC(4,1,1,1)
        WProd86 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,1,1)*
     -   ProdC(4,1,1,1)
        WProd87 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(3,1,1,1)*
     -   ProdC(4,1,1,1)
	WProd93 = Prod(3,1,1,1)**2*ProdC(3,1,2,1)*ProdC(4,1,1,1)
	WProd92 = Prod(3,1,1,1)*ProdC(4,1,2,1)
	WProd88 = Prod(3,1,2,1)**2*ProdC(3,1,1,1)*ProdC(4,1,2,1)
        WProd94 = Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(3,1,2,1)*
     -   ProdC(4,1,2,1)
      
      
	Prp17 = 1.D0/WProd82
	Prp18 = 1.D0/((0.D0,1.D0)*GZ*Sqrt(MZ2r) - MZ2r + WProd82)
      
	IInt25 = myCget(WProd82,ME2,ME2,MZ2r,MZ2r,ME2,0)
	LInt215 = myCval(cc0,IInt25)
	IInt29 = myBget(WProd82,ME2,ME2,0)
	LInt238 = myBval(bb0,IInt29)
	IInt30 = myBget(WProd82,MU2,MU2,0)
	LInt239 = myBval(bb0,IInt30)
	IInt37 = myBget(WProd82,MD2,MD2,0)
	LInt266 = myBval(bb0,IInt37)
      

	do Gen1=1,3
	IInt42(Gen1) = myBget(WProd82,MLE2(Gen1),MLE2(Gen1),1)
	LInt255(Gen1) = myBval(bb1,IInt42(Gen1))
	IInt43(Gen1) = myBget(WProd82,MQD2(Gen1),MQD2(Gen1),1)
	LInt256(Gen1) = myBval(bb1,IInt43(Gen1))
        enddo

	IInt33 = myBget(WProd82,MH2,MZ2r,2)
	LInt250 = myBval(bb0,IInt33)
	LInt251 = myBval(bb00,IInt33)
	IInt34 = myBget(WProd82,0.D0,0.D0,1)
	LInt252 = myBval(bb1,IInt34)
	LInt258 = myBval(bb0,IInt34)
      

	do Gen1=1,3
	IInt41(Gen1) = myBget(WProd82,MQU2(Gen1),MQU2(Gen1),1)
	LInt254(Gen1) = myBval(bb0,IInt41(Gen1))
	LInt257(Gen1) = myBval(bb1,IInt41(Gen1))
        enddo

	IInt28 = myBget(WProd82,MW2r,MW2r,2)
	LInt235 = myBval(bb0,IInt28)
	LInt236 = myBval(bb00,IInt28)
	LInt237 = myBval(bb1,IInt28)
	IInt36 = myCget(ME2,WProd82,ME2,0.D0,MW2r,MW2r,2)
	LInt260 = myCval(cc00,IInt36)
	LInt262 = myCval(cc1,IInt36)
	LInt264 = myCval(cc2,IInt36)
	IInt39 = myCget(MU2,WProd82,MU2,MD2,MW2r,MW2r,2)
	LInt268 = myCval(cc00,IInt39)
	LInt270 = myCval(cc1,IInt39)
	LInt272 = myCval(cc2,IInt39)
	IInt31 = myCget(ME2,WProd82,ME2,MZ2r,ME2,ME2,2)
	LInt240 = myCval(cc0,IInt31)
	LInt242 = myCval(cc00,IInt31)
	LInt244 = myCval(cc1,IInt31)
	LInt246 = myCval(cc2,IInt31)
	IInt32 = myCget(MU2,WProd82,MU2,MZ2r,MU2,MU2,2)
	LInt241 = myCval(cc0,IInt32)
	LInt243 = myCval(cc00,IInt32)
	LInt245 = myCval(cc1,IInt32)
	LInt247 = myCval(cc2,IInt32)
	IInt35 = myCget(ME2,WProd82,ME2,MW2r,0.D0,0.D0,2)
	LInt259 = myCval(cc0,IInt35)
	LInt261 = myCval(cc00,IInt35)
	LInt263 = myCval(cc1,IInt35)
	LInt265 = myCval(cc2,IInt35)
	IInt38 = myCget(MU2,WProd82,MU2,MW2r,MD2,MD2,2)
	LInt267 = myCval(cc0,IInt38)
	LInt269 = myCval(cc00,IInt38)
	LInt271 = myCval(cc1,IInt38)
	LInt273 = myCval(cc2,IInt38)
        IInt26 = myDget(MU2,WProd82,ME2,-WProd83,MU2,ME2,MU2,
     -   MZ2r,MZ2r,ME2,2)
	LInt216 = myDval(dd00,IInt26)
	LInt218 = myDval(dd1,IInt26)
	LInt220 = myDval(dd11,IInt26)
	LInt222 = myDval(dd12,IInt26)
	LInt224 = myDval(dd13,IInt26)
	LInt226 = myDval(dd2,IInt26)
	LInt228 = myDval(dd23,IInt26)
	LInt230 = myDval(dd3,IInt26)
	LInt232 = myDval(dd33,IInt26)
        IInt27 = myDget(MU2,WProd82,ME2,-WProd84,MU2,ME2,MU2,
     -   MZ2r,MZ2r,ME2,2)
	LInt217 = myDval(dd00,IInt27)
	LInt219 = myDval(dd1,IInt27)
	LInt221 = myDval(dd11,IInt27)
	LInt223 = myDval(dd12,IInt27)
	LInt225 = myDval(dd13,IInt27)
	LInt227 = myDval(dd2,IInt27)
	LInt229 = myDval(dd23,IInt27)
	LInt231 = myDval(dd3,IInt27)
	LInt233 = myDval(dd33,IInt27)
        IInt40 = myDget(MU2,WProd82,ME2,-WProd84,MU2,ME2,MD2,
     -   MW2r,MW2r,0.D0,2)
	LInt274 = myDval(dd00,IInt40)
	LInt275 = myDval(dd1,IInt40)
	LInt276 = myDval(dd11,IInt40)
	LInt277 = myDval(dd12,IInt40)
	LInt278 = myDval(dd13,IInt40)
	LInt279 = myDval(dd2,IInt40)
	LInt280 = myDval(dd23,IInt40)
	LInt281 = myDval(dd3,IInt40)
	LInt282 = myDval(dd33,IInt40)
      
      
      endif

      end


