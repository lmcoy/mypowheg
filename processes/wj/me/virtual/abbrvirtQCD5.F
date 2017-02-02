      subroutine abbrvirtQCD5(resnr)
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

      
	WProd221 = Prod(2,2,2,1)
	WProd229 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd225 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd242 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd238 = ProdC(2,2,2,1)
	WProd234 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd231 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd227 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd243 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd235 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd232 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd228 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd244 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd220 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd230 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd233 = Prod(5,1,2,1)*ProdC(5,1,2,1)
        WProd236 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd222 = -(Prod(2,2,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,1))
        WProd237 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd223 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd224 = -(Prod(2,1,1,1)*Prod(3,1,2,2)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,1))
        WProd239 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,2)*
     -     ProdC(5,1,2,1))
        WProd240 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*ProdC(4,1,1,1)*
     -     ProdC(5,1,2,2))
        WProd241 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(4,1,2,1)*
     -     ProdC(5,1,2,2))
	WProd226 = -(Prod(3,1,1,1)*ProdC(5,1,4,1))
        WProd245 = Prod(2,1,1,1)*Prod(3,1,2,1)*ProdC(2,1,2,2)*
     -   ProdC(5,1,4,1)
      
      
	Prp36 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd220)
	Prp37 = 1.D0/WProd229
	Prp38 = -(1.D0/WProd233)
      
	IInt343 = myBget(WProd220,0.D0,0.D0,0)
	LInt435 = myBval(bb0,IInt343)
	IInt344 = myBget(WProd229,0.D0,0.D0,0)
	LInt436 = myBval(bb0,IInt344)
	IInt348 = myBget(WProd229,0.D0,0.D0,1)
	LInt452 = myBval(bb1,IInt348)
	IInt349 = myBget(-WProd233,0.D0,0.D0,0)
	LInt453 = myBval(bb0,IInt349)
	IInt353 = myBget(-WProd233,0.D0,0.D0,1)
	LInt471 = myBval(bb1,IInt353)
	IInt337 = myCget(WProd229,0.D0,WProd220,0.D0,0.D0,0.D0,1)
	LInt367 = myCval(cc0,IInt337)
	LInt370 = myCval(cc1,IInt337)
	LInt373 = myCval(cc2,IInt337)
	IInt338 = myCget(WProd229,WProd220,0.D0,0.D0,0.D0,0.D0,1)
	LInt368 = myCval(cc0,IInt338)
	LInt371 = myCval(cc1,IInt338)
	LInt374 = myCval(cc2,IInt338)
	IInt339 = myCget(-WProd230,0.D0,WProd220,0.D0,0.D0,0.D0,1)
	LInt369 = myCval(cc0,IInt339)
	LInt372 = myCval(cc1,IInt339)
	LInt375 = myCval(cc2,IInt339)
	IInt352 = myCget(0.D0,-WProd233,0.D0,0.D0,0.D0,0.D0,2)
	LInt457 = myCval(cc00,IInt352)
	LInt463 = myCval(cc12,IInt352)
	LInt466 = myCval(cc2,IInt352)
	LInt469 = myCval(cc22,IInt352)
	IInt347 = myCget(0.D0,WProd229,0.D0,0.D0,0.D0,0.D0,2)
	LInt440 = myCval(cc00,IInt347)
	LInt445 = myCval(cc11,IInt347)
	LInt447 = myCval(cc12,IInt347)
	LInt472 = myCval(cc1,IInt347)
	LInt473 = myCval(cc2,IInt347)
	IInt345 = myCget(0.D0,WProd229,0.D0,0.D0,0.D0,0.D0,2)
	LInt437 = myCval(cc0,IInt345)
	LInt439 = myCval(cc00,IInt345)
	LInt442 = myCval(cc1,IInt345)
	LInt444 = myCval(cc11,IInt345)
	LInt446 = myCval(cc12,IInt345)
	LInt449 = myCval(cc2,IInt345)
	IInt346 = myCget(0.D0,WProd220,WProd229,0.D0,0.D0,0.D0,2)
	LInt438 = myCval(cc0,IInt346)
	LInt441 = myCval(cc00,IInt346)
	LInt443 = myCval(cc1,IInt346)
	LInt448 = myCval(cc12,IInt346)
	LInt450 = myCval(cc2,IInt346)
	LInt451 = myCval(cc22,IInt346)
	IInt350 = myCget(0.D0,-WProd233,0.D0,0.D0,0.D0,0.D0,2)
	LInt454 = myCval(cc0,IInt350)
	LInt456 = myCval(cc00,IInt350)
	LInt459 = myCval(cc1,IInt350)
	LInt462 = myCval(cc12,IInt350)
	LInt465 = myCval(cc2,IInt350)
	LInt468 = myCval(cc22,IInt350)
	IInt351 = myCget(0.D0,WProd220,-WProd233,0.D0,0.D0,0.D0,2)
	LInt455 = myCval(cc0,IInt351)
	LInt458 = myCval(cc00,IInt351)
	LInt460 = myCval(cc1,IInt351)
	LInt461 = myCval(cc11,IInt351)
	LInt464 = myCval(cc12,IInt351)
	LInt467 = myCval(cc2,IInt351)
	LInt470 = myCval(cc22,IInt351)
        IInt342 = myDget(0.D0,WProd229,WProd220,-WProd230,0.D0,0.D0,0.D0,
     -   0.D0,0.D0,0.D0,3)
	LInt379 = myDval(dd00,IInt342)
	LInt382 = myDval(dd001,IInt342)
	LInt385 = myDval(dd002,IInt342)
	LInt388 = myDval(dd003,IInt342)
	LInt391 = myDval(dd1,IInt342)
	LInt394 = myDval(dd11,IInt342)
	LInt397 = myDval(dd112,IInt342)
	LInt401 = myDval(dd12,IInt342)
	LInt404 = myDval(dd122,IInt342)
	LInt406 = myDval(dd123,IInt342)
	LInt409 = myDval(dd13,IInt342)
	LInt414 = myDval(dd2,IInt342)
	LInt419 = myDval(dd223,IInt342)
	LInt422 = myDval(dd23,IInt342)
	LInt425 = myDval(dd233,IInt342)
	LInt428 = myDval(dd3,IInt342)
	LInt431 = myDval(dd33,IInt342)
        IInt341 = myDget(0.D0,-WProd230,0.D0,-WProd233,0.D0,WProd220,
     -   0.D0,0.D0,0.D0,0.D0,3)
	LInt377 = myDval(dd0,IInt341)
	LInt380 = myDval(dd00,IInt341)
	LInt383 = myDval(dd001,IInt341)
	LInt386 = myDval(dd002,IInt341)
	LInt389 = myDval(dd003,IInt341)
	LInt392 = myDval(dd1,IInt341)
	LInt395 = myDval(dd11,IInt341)
	LInt399 = myDval(dd113,IInt341)
	LInt402 = myDval(dd12,IInt341)
	LInt407 = myDval(dd123,IInt341)
	LInt410 = myDval(dd13,IInt341)
	LInt412 = myDval(dd133,IInt341)
	LInt415 = myDval(dd2,IInt341)
	LInt417 = myDval(dd22,IInt341)
	LInt420 = myDval(dd223,IInt341)
	LInt423 = myDval(dd23,IInt341)
	LInt426 = myDval(dd233,IInt341)
	LInt429 = myDval(dd3,IInt341)
	LInt432 = myDval(dd33,IInt341)
	LInt434 = myDval(dd333,IInt341)
        IInt340 = myDget(0.D0,WProd229,0.D0,-WProd233,0.D0,WProd220,0.D0,
     -   0.D0,0.D0,0.D0,3)
	LInt376 = myDval(dd0,IInt340)
	LInt378 = myDval(dd00,IInt340)
	LInt381 = myDval(dd001,IInt340)
	LInt384 = myDval(dd002,IInt340)
	LInt387 = myDval(dd003,IInt340)
	LInt390 = myDval(dd1,IInt340)
	LInt393 = myDval(dd11,IInt340)
	LInt396 = myDval(dd112,IInt340)
	LInt398 = myDval(dd113,IInt340)
	LInt400 = myDval(dd12,IInt340)
	LInt403 = myDval(dd122,IInt340)
	LInt405 = myDval(dd123,IInt340)
	LInt408 = myDval(dd13,IInt340)
	LInt411 = myDval(dd133,IInt340)
	LInt413 = myDval(dd2,IInt340)
	LInt416 = myDval(dd22,IInt340)
	LInt418 = myDval(dd223,IInt340)
	LInt421 = myDval(dd23,IInt340)
	LInt424 = myDval(dd233,IInt340)
	LInt427 = myDval(dd3,IInt340)
	LInt430 = myDval(dd33,IInt340)
	LInt433 = myDval(dd333,IInt340)
      
      
      endif

      end


