      subroutine abbrvirtQCD4(resnr)
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

      
	WProd189 = Prod(2,2,2,1)
	WProd196 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd205 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd219 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd209 = ProdC(2,2,2,1)
	WProd197 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd199 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd204 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd217 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd198 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd200 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd192 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd212 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd201 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd195 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd206 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd207 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd202 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd208 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd190 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd203 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd215 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
        WProd216 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,2))
        WProd210 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
        WProd191 = -(Prod(3,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,3,1))
        WProd211 = -(Prod(3,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,3,1))
	WProd193 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
        WProd213 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd218 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
        WProd194 = -(Prod(3,1,2,2)*Prod(4,1,1,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,4,1))
        WProd214 = -(Prod(3,1,2,1)*Prod(4,1,1,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,4,1))
      
      
	Prp31 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r - WProd195)
	Prp34 = 1.D0/WProd196
	Prp32 = -(1.D0/WProd200)
	Prp33 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd201)
	Prp35 = -(1.D0/WProd206)
      
	IInt239 = myBget(-WProd195,0.D0,0.D0,0)
	LInt292 = myBval(bb0,IInt239)
	IInt245 = myBget(WProd196,0.D0,0.D0,0)
	LInt340 = myBval(bb0,IInt245)
	IInt246 = myBget(WProd201,0.D0,0.D0,0)
	LInt341 = myBval(bb0,IInt246)
	IInt249 = myBget(WProd196,0.D0,0.D0,1)
	LInt351 = myBval(bb1,IInt249)
	IInt250 = myBget(-WProd206,0.D0,0.D0,0)
	LInt352 = myBval(bb0,IInt250)
	IInt253 = myBget(-WProd206,0.D0,0.D0,1)
	LInt364 = myBval(bb1,IInt253)
	IInt241 = myCget(-WProd195,0.D0,WProd201,0.D0,0.D0,0.D0,1)
	LInt297 = myCval(cc0,IInt241)
	LInt299 = myCval(cc1,IInt241)
	LInt301 = myCval(cc2,IInt241)
	IInt242 = myCget(WProd196,WProd201,0.D0,0.D0,0.D0,0.D0,1)
	LInt298 = myCval(cc0,IInt242)
	LInt300 = myCval(cc1,IInt242)
	LInt302 = myCval(cc2,IInt242)
	IInt240 = myCget(0.D0,-WProd195,0.D0,0.D0,0.D0,0.D0,2)
	LInt293 = myCval(cc0,IInt240)
	LInt294 = myCval(cc00,IInt240)
	LInt295 = myCval(cc1,IInt240)
	LInt296 = myCval(cc2,IInt240)
	IInt252 = myCget(0.D0,-WProd206,0.D0,0.D0,0.D0,0.D0,2)
	LInt354 = myCval(cc00,IInt252)
	LInt358 = myCval(cc12,IInt252)
	LInt360 = myCval(cc2,IInt252)
	LInt362 = myCval(cc22,IInt252)
	IInt248 = myCget(0.D0,WProd196,0.D0,0.D0,0.D0,0.D0,2)
	LInt343 = myCval(cc00,IInt248)
	LInt346 = myCval(cc11,IInt248)
	LInt347 = myCval(cc12,IInt248)
	LInt365 = myCval(cc1,IInt248)
	LInt366 = myCval(cc2,IInt248)
	IInt247 = myCget(0.D0,WProd201,WProd196,0.D0,0.D0,0.D0,2)
	LInt342 = myCval(cc0,IInt247)
	LInt344 = myCval(cc00,IInt247)
	LInt345 = myCval(cc1,IInt247)
	LInt348 = myCval(cc12,IInt247)
	LInt349 = myCval(cc2,IInt247)
	LInt350 = myCval(cc22,IInt247)
	IInt251 = myCget(0.D0,WProd201,-WProd206,0.D0,0.D0,0.D0,2)
	LInt353 = myCval(cc0,IInt251)
	LInt355 = myCval(cc00,IInt251)
	LInt356 = myCval(cc1,IInt251)
	LInt357 = myCval(cc11,IInt251)
	LInt359 = myCval(cc12,IInt251)
	LInt361 = myCval(cc2,IInt251)
	LInt363 = myCval(cc22,IInt251)
        IInt244 = myDget(0.D0,WProd196,WProd201,-WProd195,0.D0,0.D0,0.D0,
     -   0.D0,0.D0,0.D0,3)
	LInt305 = myDval(dd00,IInt244)
	LInt307 = myDval(dd001,IInt244)
	LInt309 = myDval(dd002,IInt244)
	LInt311 = myDval(dd003,IInt244)
	LInt313 = myDval(dd1,IInt244)
	LInt315 = myDval(dd11,IInt244)
	LInt316 = myDval(dd112,IInt244)
	LInt319 = myDval(dd12,IInt244)
	LInt320 = myDval(dd122,IInt244)
	LInt322 = myDval(dd123,IInt244)
	LInt324 = myDval(dd13,IInt244)
	LInt327 = myDval(dd2,IInt244)
	LInt330 = myDval(dd223,IInt244)
	LInt332 = myDval(dd23,IInt244)
	LInt334 = myDval(dd233,IInt244)
	LInt336 = myDval(dd3,IInt244)
	LInt338 = myDval(dd33,IInt244)
        IInt243 = myDget(0.D0,-WProd195,0.D0,-WProd206,0.D0,WProd201,
     -   0.D0,0.D0,0.D0,0.D0,3)
	LInt303 = myDval(dd0,IInt243)
	LInt304 = myDval(dd00,IInt243)
	LInt306 = myDval(dd001,IInt243)
	LInt308 = myDval(dd002,IInt243)
	LInt310 = myDval(dd003,IInt243)
	LInt312 = myDval(dd1,IInt243)
	LInt314 = myDval(dd11,IInt243)
	LInt317 = myDval(dd113,IInt243)
	LInt318 = myDval(dd12,IInt243)
	LInt321 = myDval(dd123,IInt243)
	LInt323 = myDval(dd13,IInt243)
	LInt325 = myDval(dd133,IInt243)
	LInt326 = myDval(dd2,IInt243)
	LInt328 = myDval(dd22,IInt243)
	LInt329 = myDval(dd223,IInt243)
	LInt331 = myDval(dd23,IInt243)
	LInt333 = myDval(dd233,IInt243)
	LInt335 = myDval(dd3,IInt243)
	LInt337 = myDval(dd33,IInt243)
	LInt339 = myDval(dd333,IInt243)
      
      
      endif

      end


