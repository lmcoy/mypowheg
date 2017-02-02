      subroutine abbrrealQED31(resnr)
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

      
	WProd2499 = Prod(2,2,2,1)
	WProd2428 = Prod(1,2,1,1)*Prod(2,2,2,1)
	WProd2487 = Prod(2,2,2,1)*ProdC(1,2,1,1)
	WProd2426 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd2434 = Prod(2,1,1,2)*ProdC(2,1,1,1)
	WProd2429 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd2432 = Prod(2,2,1,2)*ProdC(2,1,1,1)
	WProd2490 = Prod(2,1,1,1)*ProdC(2,1,1,2)
	WProd2489 = Prod(2,2,1,1)*ProdC(2,1,1,2)
	WProd2470 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd2471 = Prod(2,1,1,2)*ProdC(2,2,1,1)
	WProd2509 = Prod(2,1,1,1)*ProdC(2,2,1,2)
	WProd2510 = ProdC(2,2,2,1)
	WProd2469 = Prod(1,2,1,1)*ProdC(2,2,2,1)
	WProd2508 = ProdC(1,2,1,1)*ProdC(2,2,2,1)
	WProd2443 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd2456 = Prod(3,1,1,2)*ProdC(3,1,1,1)
	WProd2497 = Prod(3,1,1,1)*ProdC(3,1,1,2)
	WProd2444 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd2457 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd2479 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd2445 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd2458 = Prod(4,1,1,2)*ProdC(4,1,1,1)
	WProd2498 = Prod(4,1,1,1)*ProdC(4,1,1,2)
	WProd2446 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd2459 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd2480 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd2427 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd2447 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd2448 = Prod(5,1,1,2)*ProdC(5,1,1,1)
        WProd2430 = Prod(3,1,1,2)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd2431 = Prod(3,1,1,1)*Prod(6,1,1,2)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd2488 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,2)*
     -   ProdC(5,1,1,1)
        WProd2433 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd2435 = Prod(3,1,2,2)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd2436 = Prod(3,1,2,1)*Prod(6,1,1,2)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd2462 = Prod(3,1,2,2)*Prod(6,1,1,2)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd2463 = Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,1,2)*
     -   ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(5,1,1,1)
        WProd2501 = Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -   ProdC(2,1,1,2)*ProdC(4,1,2,1)*ProdC(5,1,1,1)
        WProd2481 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,2)*
     -   ProdC(2,2,1,1)*ProdC(4,1,2,1)*ProdC(5,1,1,1)
        WProd2512 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -   ProdC(2,2,1,2)*ProdC(4,1,2,1)*ProdC(5,1,1,1)
        WProd2472 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,2)*
     -   ProdC(5,1,1,1)
        WProd2482 = Prod(3,1,2,1)*Prod(6,1,1,2)*ProdC(4,1,2,2)*
     -   ProdC(5,1,1,1)
	WProd2494 = Prod(5,1,1,1)*ProdC(5,1,1,2)
        WProd2491 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,2)
        WProd2492 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,2)
        WProd2502 = Prod(3,1,2,2)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,2)
        WProd2513 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,2)*
     -   ProdC(5,1,1,2)
	WProd2449 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd2450 = Prod(5,1,2,2)*ProdC(5,1,2,1)
        WProd2437 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd2438 = Prod(3,1,1,2)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd2439 = Prod(3,1,1,1)*Prod(6,1,2,2)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd2464 = Prod(3,1,1,2)*Prod(6,1,2,2)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd2465 = -(Prod(2,1,2,2)*Prod(3,1,1,1)*Prod(6,1,1,2)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(5,1,2,1))
        WProd2466 = Prod(2,1,1,1)*Prod(3,1,1,2)*Prod(6,1,2,2)*
     -   ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(5,1,2,1)
        WProd2503 = -(Prod(2,1,2,2)*Prod(3,1,1,1)*Prod(6,1,1,1)*
     -     ProdC(2,1,1,2)*ProdC(4,1,1,1)*ProdC(5,1,2,1))
        WProd2483 = Prod(2,1,1,1)*Prod(3,1,1,2)*Prod(6,1,2,1)*
     -   ProdC(2,2,1,1)*ProdC(4,1,1,1)*ProdC(5,1,2,1)
        WProd2493 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,2)*
     -   ProdC(5,1,2,1)
        WProd2504 = Prod(3,1,1,1)*Prod(6,1,2,2)*ProdC(4,1,1,2)*
     -   ProdC(5,1,2,1)
        WProd2505 = Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,2)*
     -   ProdC(2,1,1,1)*ProdC(4,1,1,2)*ProdC(5,1,2,1)
        WProd2514 = Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -   ProdC(2,2,1,1)*ProdC(4,1,1,2)*ProdC(5,1,2,1)
        WProd2500 = Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,1,1,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,2,1)
        WProd2440 = Prod(3,1,2,2)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd2441 = Prod(3,1,2,1)*Prod(6,1,2,2)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd2511 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -   ProdC(4,1,2,2)*ProdC(5,1,2,1)
        WProd2473 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,2)*
     -   ProdC(5,1,2,1)
	WProd2476 = Prod(5,1,2,1)*ProdC(5,1,2,2)
        WProd2474 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,2)
        WProd2484 = Prod(3,1,1,2)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,2)
        WProd2515 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,2)*
     -   ProdC(5,1,2,2)
        WProd2475 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,2)
	WProd2451 = Prod(5,1,3,1)*ProdC(5,1,3,1)
	WProd2452 = Prod(5,1,4,1)*ProdC(5,1,4,1)
	WProd2442 = Prod(6,1,3,1)*ProdC(5,1,4,1)
        WProd2495 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(1,1,1,2)*
     -     ProdC(5,1,4,1))
        WProd2453 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd2454 = -(Prod(3,1,2,1)*Prod(6,1,1,2)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd2467 = -(Prod(3,1,2,2)*Prod(6,1,1,2)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd2455 = Prod(3,1,1,1)*Prod(6,1,2,2)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd2468 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,1,2)*
     -     ProdC(2,1,1,1)**2*ProdC(5,1,4,1))
        WProd2496 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(2,1,1,2)*
     -     ProdC(5,1,4,1))
        WProd2506 = -(Prod(3,1,2,2)*Prod(6,1,1,1)*ProdC(2,1,1,2)*
     -     ProdC(5,1,4,1))
        WProd2507 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(2,1,1,1)*ProdC(2,1,1,2)*ProdC(5,1,4,1))
        WProd2477 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd2485 = -(Prod(3,1,2,1)*Prod(6,1,1,2)*ProdC(2,2,1,1)*
     -     ProdC(5,1,4,1))
        WProd2478 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(2,2,1,1)*
     -   ProdC(5,1,4,1)
        WProd2486 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,2)*
     -     ProdC(2,1,1,1)*ProdC(2,2,1,1)*ProdC(5,1,4,1))
        WProd2516 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(2,2,1,2)*
     -     ProdC(5,1,4,1))
        WProd2517 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(2,1,1,1)*ProdC(2,2,1,2)*ProdC(5,1,4,1))
	WProd2460 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd2461 = Prod(6,1,2,1)*ProdC(6,1,2,1)
      
      
	Prp254 = 1.D0/WProd2426
	Prp255 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd2427)
	Prp260 = -(1.D0/WProd2447)
	Prp256 = 1.D0/(WProd2426 - WProd2447 - WProd2449)
	Prp258 = -(1.D0/WProd2449)
	Prp259 = -(1.D0/WProd2460)
	Prp257 = 1.D0/(WProd2426 - WProd2460 - WProd2461)
	Prp261 = -(1.D0/WProd2461)
      
      
      
      endif

      end


