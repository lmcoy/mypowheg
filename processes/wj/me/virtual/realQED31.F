      subroutine realQED31(amp,resnr)
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "params.h" 

      integer helnr,colnr,resnr,procnr

      double complex amp(nproc,ncol,nhel)

      double precision Alfas,Alfas2,GS

      common /my_alfas_common/Alfas,Alfas2,GS
#include "amps/ampvars.h" 
#include "amps/ampcom.h" 

*The Amplitude is evaluated as amp(procnr,colnr,helnr), where
*
*   Resonances:
*     resnr=0:   non-resonant
*
*   Partonic Subprocesses:
*     procnr=1:   {~up  ,  up  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=2:   {~up  , ~do  ,  ve  , ~el  , ~up  , ~up  }
*     procnr=3:   { up  ,  up  ,  ve  , ~el  ,  up  ,  do  }
*     procnr=4:   { up  , ~ch  ,  ve  , ~el  , ~ch  ,  do  }
*     procnr=5:   { up  ,  ch  ,  ve  , ~el  ,  ch  ,  do  }
*     procnr=6:   { up  , ~do  ,  ve  , ~el  , ~up  ,  up  }
*     procnr=7:   { up  , ~do  ,  ve  , ~el  , ~ch  ,  ch  }
*     procnr=8:   { up  , ~do  ,  ve  , ~el  , ~do  ,  do  }
*     procnr=9:   { up  , ~do  ,  ve  , ~el  , ~st  ,  st  }
*     procnr=10:   { up  , ~do  ,  ve  , ~el  , ~bo  ,  bo  }
*     procnr=11:   { up  , ~do  ,  ve  , ~el  ,  gl  ,  ga  }
*     procnr=12:   { up  , ~do  ,  ve  , ~el  ,  gl  ,  gl  }
*     procnr=13:   { up  ,  do  ,  ve  , ~el  ,  do  ,  do  }
*     procnr=14:   { up  , ~st  ,  ve  , ~el  ,  do  , ~st  }
*     procnr=15:   { up  ,  st  ,  ve  , ~el  ,  do  ,  st  }
*     procnr=16:   { up  ,  ga  ,  ve  , ~el  ,  do  ,  gl  }
*     procnr=17:   { up  ,  gl  ,  ve  , ~el  ,  do  ,  ga  }
*     procnr=18:   { up  ,  gl  ,  ve  , ~el  ,  do  ,  gl  }
*     procnr=19:   {~ch  ,  ch  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=20:   {~ch  , ~do  ,  ve  , ~el  , ~up  , ~ch  }
*     procnr=21:   { ch  , ~do  ,  ve  , ~el  , ~up  ,  ch  }
*     procnr=22:   {~do  , ~do  ,  ve  , ~el  , ~up  , ~do  }
*     procnr=23:   {~do  ,  do  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=24:   {~do  , ~st  ,  ve  , ~el  , ~up  , ~st  }
*     procnr=25:   {~do  ,  st  ,  ve  , ~el  , ~up  ,  st  }
*     procnr=26:   {~do  ,  ga  ,  ve  , ~el  , ~up  ,  gl  }
*     procnr=27:   {~do  ,  gl  ,  ve  , ~el  , ~up  ,  ga  }
*     procnr=28:   {~do  ,  gl  ,  ve  , ~el  , ~up  ,  gl  }
*     procnr=29:   {~st  ,  st  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=30:   { ga  ,  gl  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=31:   { gl  ,  gl  ,  ve  , ~el  , ~up  ,  do  }
*
*   Helicities:
*     helnr=1:   {-1, -1, -1, 1, -1, -1}
*     helnr=2:   {-1, 1, -1, 1, 1, -1}
*     helnr=3:   {1, -1, -1, 1, 1, -1}
*     helnr=4:   {-1, 1, -1, 1, -1, 1}
*     helnr=5:   {1, 1, -1, 1, 1, 1}
*     helnr=6:   {-1, 1, -1, 1, -1, -1}
*     helnr=7:   {-1, 1, -1, 1, 1, 1}
*     helnr=8:   {-1, -1, -1, 1, -1, 1}
*     helnr=9:   {1, -1, -1, 1, 1, 1}
*     helnr=10:   {1, 1, -1, 1, 1, -1}
*     helnr=11:   {-1, -1, -1, 1, 1, -1}
*
*   Color Structures:
*     colnr=1:   SUNT[Col1, Col3] SUNT[Col2, Col4]
*     colnr=2:   SUNT[Col1, Col2] SUNT[Col3, Col4]
*     colnr=3:   SUNT[Col1, Col4] SUNT[Col2, Col3]
*     colnr=4:   SUNT[Glu1, Col2, Col1]
*     colnr=5:   SUNT[Glu1, Glu2, Col2, Col1]
*     colnr=6:   SUNT[Glu2, Glu1, Col2, Col1]
*     colnr=7:   SUNT[Glu1, Col1, Col2]
*     colnr=8:   SUNT[Glu1, Glu2, Col1, Col2]
*     colnr=9:   SUNT[Glu2, Glu1, Col1, Col2]


      if(resnr.eq.0)then

        amp(31,5,11) = Prp254*Prp255*
     -    ((-16*Alfas*Cpl122*Prp256*
     -         (2*WProd2429*WProd2430 - 2*WProd2429*WProd2431 + 
     -           2*WProd2432*WProd2433 - 2*WProd2434*WProd2435 - 
     -           2*WProd2429*WProd2436 - 2*WProd2432*WProd2437 + 
     -           2*WProd2429*WProd2438 + 2*WProd2434*WProd2439 - 
     -           2*WProd2434*WProd2440 + 2*WProd2434*WProd2441 + 
     -           4*WProd2426*WProd2432*WProd2442 + 
     -           3*WProd2427*WProd2432*WProd2442 - 
     -           3*WProd2432*WProd2442*WProd2443 - 
     -           3*WProd2432*WProd2442*WProd2444 - 
     -           3*WProd2432*WProd2442*WProd2445 - 
     -           3*WProd2432*WProd2442*WProd2446 - 
     -           4*WProd2432*WProd2442*WProd2447 + 
     -           2*WProd2429*WProd2442*WProd2448 - 
     -           2*WProd2432*WProd2442*WProd2449 - 
     -           2*WProd2434*WProd2442*WProd2450 + 
     -           3*WProd2432*WProd2442*WProd2451 + 
     -           3*WProd2432*WProd2442*WProd2452 - 
     -           2*WProd2432*WProd2453 + 2*WProd2429*WProd2454 - 
     -           2*WProd2434*WProd2455))/WProd2428 - 
     -      (16*Alfas*Cpl122*Prp257*
     -         (4*WProd2426*WProd2432*WProd2442 + 
     -           3*WProd2427*WProd2432*WProd2442 - 
     -           4*WProd2432*WProd2442*WProd2443 - 
     -           2*WProd2432*WProd2442*WProd2444 - 
     -           4*WProd2432*WProd2442*WProd2445 - 
     -           2*WProd2432*WProd2442*WProd2446 - 
     -           4*WProd2432*WProd2442*WProd2447 + 
     -           2*WProd2429*WProd2442*WProd2448 - 
     -           2*WProd2432*WProd2442*WProd2449 - 
     -           2*WProd2434*WProd2442*WProd2450 + 
     -           3*WProd2432*WProd2442*WProd2451 + 
     -           3*WProd2432*WProd2442*WProd2452 - 
     -           2*WProd2432*WProd2453 + 2*WProd2429*WProd2454 - 
     -           2*WProd2434*WProd2455 + 
     -           2*WProd2429*WProd2442*WProd2456 - 
     -           2*WProd2434*WProd2442*WProd2457 + 
     -           2*WProd2429*WProd2442*WProd2458 - 
     -           2*WProd2434*WProd2442*WProd2459))/WProd2428) + 
     -   Prp255*((32*Alfas*Cpl122*Prp256*Prp258*
     -         (WProd2429*WProd2436 - WProd2429*WProd2438 + 
     -           WProd2434*WProd2440 - WProd2434*WProd2441 - 
     -           WProd2440*WProd2448 + WProd2441*WProd2448 - 
     -           WProd2430*WProd2450 + WProd2431*WProd2450 + 
     -           WProd2434*WProd2442*WProd2450 - 
     -           WProd2442*WProd2448*WProd2450 - 
     -           WProd2429*WProd2454 + 2*WProd2426*WProd2462 + 
     -           2*WProd2427*WProd2462 - 2*WProd2443*WProd2462 - 
     -           2*WProd2444*WProd2462 - 2*WProd2445*WProd2462 - 
     -           2*WProd2446*WProd2462 - 2*WProd2447*WProd2462 - 
     -           WProd2449*WProd2462 + 2*WProd2451*WProd2462 + 
     -           2*WProd2452*WProd2462 - WProd2463 - 
     -           2*WProd2426*WProd2464 - 2*WProd2427*WProd2464 + 
     -           2*WProd2443*WProd2464 + 2*WProd2444*WProd2464 + 
     -           2*WProd2445*WProd2464 + 2*WProd2446*WProd2464 + 
     -           2*WProd2447*WProd2464 + WProd2449*WProd2464 - 
     -           2*WProd2451*WProd2464 - 2*WProd2452*WProd2464 + 
     -           WProd2465 + WProd2466 - 2*WProd2426*WProd2467 - 
     -           2*WProd2427*WProd2467 + 2*WProd2443*WProd2467 + 
     -           2*WProd2444*WProd2467 + 2*WProd2445*WProd2467 + 
     -           2*WProd2446*WProd2467 + 2*WProd2447*WProd2467 + 
     -           WProd2449*WProd2467 - 2*WProd2451*WProd2467 - 
     -           2*WProd2452*WProd2467 + WProd2468))/WProd2428 + 
     -      (32*Alfas*Cpl122*Prp258*Prp259*
     -         (WProd2429*WProd2436 + WProd2434*WProd2440 - 
     -           WProd2434*WProd2441 - WProd2440*WProd2448 + 
     -           WProd2441*WProd2448 + 
     -           WProd2434*WProd2442*WProd2450 - 
     -           WProd2442*WProd2448*WProd2450 - 
     -           WProd2429*WProd2454 - WProd2440*WProd2456 + 
     -           WProd2441*WProd2456 - 
     -           WProd2442*WProd2450*WProd2456 - 
     -           WProd2440*WProd2458 + WProd2441*WProd2458 - 
     -           WProd2442*WProd2450*WProd2458 + 
     -           2*WProd2426*WProd2462 + 2*WProd2427*WProd2462 - 
     -           2*WProd2443*WProd2462 - 2*WProd2444*WProd2462 - 
     -           2*WProd2445*WProd2462 - 2*WProd2446*WProd2462 - 
     -           2*WProd2447*WProd2462 - 2*WProd2449*WProd2462 + 
     -           2*WProd2451*WProd2462 + 2*WProd2452*WProd2462 - 
     -           WProd2463 - 2*WProd2426*WProd2467 - 
     -           2*WProd2427*WProd2467 + 2*WProd2443*WProd2467 + 
     -           2*WProd2444*WProd2467 + 2*WProd2445*WProd2467 + 
     -           2*WProd2446*WProd2467 + 2*WProd2447*WProd2467 + 
     -           2*WProd2449*WProd2467 - 2*WProd2451*WProd2467 - 
     -           2*WProd2452*WProd2467 + WProd2468))/WProd2428 + 
     -      (32*Alfas*Cpl122*Prp257*Prp259*
     -         (WProd2434*WProd2442*WProd2450 - 
     -           WProd2442*WProd2448*WProd2450 - 
     -           WProd2429*WProd2454 - 
     -           WProd2442*WProd2450*WProd2456 + 
     -           WProd2434*WProd2442*WProd2457 - 
     -           WProd2442*WProd2448*WProd2457 - 
     -           WProd2442*WProd2456*WProd2457 - 
     -           WProd2442*WProd2450*WProd2458 - 
     -           WProd2442*WProd2457*WProd2458 + 
     -           WProd2434*WProd2442*WProd2459 - 
     -           WProd2442*WProd2448*WProd2459 - 
     -           WProd2442*WProd2456*WProd2459 - 
     -           WProd2442*WProd2458*WProd2459 - 
     -           2*WProd2426*WProd2467 - WProd2427*WProd2467 + 
     -           2*WProd2443*WProd2467 + WProd2444*WProd2467 + 
     -           2*WProd2445*WProd2467 + WProd2446*WProd2467 + 
     -           2*WProd2447*WProd2467 + WProd2449*WProd2467 - 
     -           WProd2451*WProd2467 - WProd2452*WProd2467 + 
     -           WProd2468))/WProd2428)

#ifdef DEBUG
	print *, 'amp(31,5,11) =', amp(31,5,11)
#endif

        amp(31,5,2) = Prp254*Prp255*
     -   ((-16*Alfas*Cpl122*Prp256*
     -        (2*WProd2430*WProd2470 - 2*WProd2431*WProd2470 - 
     -          2*WProd2436*WProd2470 + 2*WProd2438*WProd2470 + 
     -          2*WProd2442*WProd2448*WProd2470 + 
     -          2*WProd2454*WProd2470 + 2*WProd2433*WProd2471 - 
     -          2*WProd2437*WProd2471 + 
     -          4*WProd2426*WProd2442*WProd2471 + 
     -          3*WProd2427*WProd2442*WProd2471 - 
     -          3*WProd2442*WProd2443*WProd2471 - 
     -          3*WProd2442*WProd2444*WProd2471 - 
     -          3*WProd2442*WProd2445*WProd2471 - 
     -          3*WProd2442*WProd2446*WProd2471 - 
     -          4*WProd2442*WProd2447*WProd2471 - 
     -          2*WProd2442*WProd2449*WProd2471 + 
     -          3*WProd2442*WProd2451*WProd2471 + 
     -          3*WProd2442*WProd2452*WProd2471 - 
     -          2*WProd2453*WProd2471 - 2*WProd2434*WProd2472 - 
     -          2*WProd2434*WProd2473 + 2*WProd2434*WProd2474 + 
     -          2*WProd2434*WProd2475 - 
     -          2*WProd2434*WProd2442*WProd2476 - 
     -          2*WProd2434*WProd2477 - 2*WProd2434*WProd2478))/
     -      WProd2469 - (16*Alfas*Cpl122*Prp257*
     -        (2*WProd2442*WProd2448*WProd2470 + 
     -          2*WProd2454*WProd2470 + 
     -          2*WProd2442*WProd2456*WProd2470 + 
     -          2*WProd2442*WProd2458*WProd2470 + 
     -          4*WProd2426*WProd2442*WProd2471 + 
     -          3*WProd2427*WProd2442*WProd2471 - 
     -          4*WProd2442*WProd2443*WProd2471 - 
     -          2*WProd2442*WProd2444*WProd2471 - 
     -          4*WProd2442*WProd2445*WProd2471 - 
     -          2*WProd2442*WProd2446*WProd2471 - 
     -          4*WProd2442*WProd2447*WProd2471 - 
     -          2*WProd2442*WProd2449*WProd2471 + 
     -          3*WProd2442*WProd2451*WProd2471 + 
     -          3*WProd2442*WProd2452*WProd2471 - 
     -          2*WProd2453*WProd2471 - 
     -          2*WProd2434*WProd2442*WProd2476 - 
     -          2*WProd2434*WProd2477 - 2*WProd2434*WProd2478 - 
     -          2*WProd2434*WProd2442*WProd2479 - 
     -          2*WProd2434*WProd2442*WProd2480))/WProd2469)
        amp(31,5,2) = Prp255*
     -    ((32*Alfas*Cpl122*Prp256*Prp258*
     -         (WProd2436*WProd2470 - WProd2438*WProd2470 - 
     -           WProd2454*WProd2470 + WProd2434*WProd2473 - 
     -           WProd2448*WProd2473 - WProd2434*WProd2475 + 
     -           WProd2448*WProd2475 - WProd2430*WProd2476 + 
     -           WProd2431*WProd2476 + 
     -           WProd2434*WProd2442*WProd2476 - 
     -           WProd2442*WProd2448*WProd2476 + 
     -           WProd2434*WProd2477 - WProd2448*WProd2477 - 
     -           WProd2481 + 2*WProd2426*WProd2482 + 
     -           2*WProd2427*WProd2482 - 2*WProd2443*WProd2482 - 
     -           2*WProd2444*WProd2482 - 2*WProd2445*WProd2482 - 
     -           2*WProd2446*WProd2482 - 2*WProd2447*WProd2482 - 
     -           WProd2449*WProd2482 + 2*WProd2451*WProd2482 + 
     -           2*WProd2452*WProd2482 + WProd2483 - 
     -           2*WProd2426*WProd2484 - 2*WProd2427*WProd2484 + 
     -           2*WProd2443*WProd2484 + 2*WProd2444*WProd2484 + 
     -           2*WProd2445*WProd2484 + 2*WProd2446*WProd2484 + 
     -           2*WProd2447*WProd2484 + WProd2449*WProd2484 - 
     -           2*WProd2451*WProd2484 - 2*WProd2452*WProd2484 - 
     -           2*WProd2426*WProd2485 - 2*WProd2427*WProd2485 + 
     -           2*WProd2443*WProd2485 + 2*WProd2444*WProd2485 + 
     -           2*WProd2445*WProd2485 + 2*WProd2446*WProd2485 + 
     -           2*WProd2447*WProd2485 + WProd2449*WProd2485 - 
     -           2*WProd2451*WProd2485 - 2*WProd2452*WProd2485 + 
     -           WProd2486))/WProd2469 + 
     -      (32*Alfas*Cpl122*Prp258*Prp259*
     -         (WProd2436*WProd2470 - WProd2454*WProd2470 + 
     -           WProd2434*WProd2473 - WProd2448*WProd2473 - 
     -           WProd2456*WProd2473 - WProd2458*WProd2473 - 
     -           WProd2434*WProd2475 + WProd2448*WProd2475 + 
     -           WProd2456*WProd2475 + WProd2458*WProd2475 + 
     -           WProd2434*WProd2442*WProd2476 - 
     -           WProd2442*WProd2448*WProd2476 - 
     -           WProd2442*WProd2456*WProd2476 - 
     -           WProd2442*WProd2458*WProd2476 + 
     -           WProd2434*WProd2477 - WProd2448*WProd2477 - 
     -           WProd2456*WProd2477 - WProd2458*WProd2477 - 
     -           WProd2481 + 2*WProd2426*WProd2482 + 
     -           2*WProd2427*WProd2482 - 2*WProd2443*WProd2482 - 
     -           2*WProd2444*WProd2482 - 2*WProd2445*WProd2482 - 
     -           2*WProd2446*WProd2482 - 2*WProd2447*WProd2482 - 
     -           2*WProd2449*WProd2482 + 2*WProd2451*WProd2482 + 
     -           2*WProd2452*WProd2482 - 2*WProd2426*WProd2485 - 
     -           2*WProd2427*WProd2485 + 2*WProd2443*WProd2485 + 
     -           2*WProd2444*WProd2485 + 2*WProd2445*WProd2485 + 
     -           2*WProd2446*WProd2485 + 2*WProd2447*WProd2485 + 
     -           2*WProd2449*WProd2485 - 2*WProd2451*WProd2485 - 
     -           2*WProd2452*WProd2485 + WProd2486))/WProd2469 + 
     -      (32*Alfas*Cpl122*Prp257*Prp259*
     -         (-(WProd2454*WProd2470) + 
     -           WProd2434*WProd2442*WProd2476 - 
     -           WProd2442*WProd2448*WProd2476 - 
     -           WProd2442*WProd2456*WProd2476 - 
     -           WProd2442*WProd2458*WProd2476 + 
     -           WProd2434*WProd2477 - WProd2448*WProd2477 - 
     -           WProd2456*WProd2477 - WProd2458*WProd2477 + 
     -           WProd2434*WProd2442*WProd2479 - 
     -           WProd2442*WProd2448*WProd2479 - 
     -           WProd2442*WProd2456*WProd2479 - 
     -           WProd2442*WProd2458*WProd2479 + 
     -           WProd2434*WProd2442*WProd2480 - 
     -           WProd2442*WProd2448*WProd2480 - 
     -           WProd2442*WProd2456*WProd2480 - 
     -           WProd2442*WProd2458*WProd2480 - 
     -           2*WProd2426*WProd2485 - WProd2427*WProd2485 + 
     -           2*WProd2443*WProd2485 + WProd2444*WProd2485 + 
     -           2*WProd2445*WProd2485 + WProd2446*WProd2485 + 
     -           2*WProd2447*WProd2485 + WProd2449*WProd2485 - 
     -           WProd2451*WProd2485 - WProd2452*WProd2485 + 
     -           WProd2486))/WProd2469) + amp(31,5,2)

#ifdef DEBUG
	print *, 'amp(31,5,2) =', amp(31,5,2)
#endif

        amp(31,5,3) = Prp254*Prp255*
     -   ((-16*Alfas*Cpl122*Prp256*
     -        (2*WProd2429*WProd2488 + 2*WProd2433*WProd2489 - 
     -          2*WProd2437*WProd2489 + 
     -          4*WProd2426*WProd2442*WProd2489 + 
     -          3*WProd2427*WProd2442*WProd2489 - 
     -          3*WProd2442*WProd2443*WProd2489 - 
     -          3*WProd2442*WProd2444*WProd2489 - 
     -          3*WProd2442*WProd2445*WProd2489 - 
     -          3*WProd2442*WProd2446*WProd2489 - 
     -          4*WProd2442*WProd2447*WProd2489 - 
     -          2*WProd2442*WProd2449*WProd2489 + 
     -          3*WProd2442*WProd2451*WProd2489 + 
     -          3*WProd2442*WProd2452*WProd2489 - 
     -          2*WProd2453*WProd2489 - 2*WProd2435*WProd2490 + 
     -          2*WProd2439*WProd2490 - 2*WProd2440*WProd2490 + 
     -          2*WProd2441*WProd2490 - 
     -          2*WProd2442*WProd2450*WProd2490 - 
     -          2*WProd2455*WProd2490 - 2*WProd2429*WProd2491 - 
     -          2*WProd2429*WProd2492 + 2*WProd2429*WProd2493 + 
     -          2*WProd2429*WProd2442*WProd2494 + 
     -          2*WProd2429*WProd2495 + 2*WProd2429*WProd2496))/
     -      WProd2487 - (16*Alfas*Cpl122*Prp257*
     -        (4*WProd2426*WProd2442*WProd2489 + 
     -          3*WProd2427*WProd2442*WProd2489 - 
     -          4*WProd2442*WProd2443*WProd2489 - 
     -          2*WProd2442*WProd2444*WProd2489 - 
     -          4*WProd2442*WProd2445*WProd2489 - 
     -          2*WProd2442*WProd2446*WProd2489 - 
     -          4*WProd2442*WProd2447*WProd2489 - 
     -          2*WProd2442*WProd2449*WProd2489 + 
     -          3*WProd2442*WProd2451*WProd2489 + 
     -          3*WProd2442*WProd2452*WProd2489 - 
     -          2*WProd2453*WProd2489 - 
     -          2*WProd2442*WProd2450*WProd2490 - 
     -          2*WProd2455*WProd2490 - 
     -          2*WProd2442*WProd2457*WProd2490 - 
     -          2*WProd2442*WProd2459*WProd2490 + 
     -          2*WProd2429*WProd2442*WProd2494 + 
     -          2*WProd2429*WProd2495 + 2*WProd2429*WProd2496 + 
     -          2*WProd2429*WProd2442*WProd2497 + 
     -          2*WProd2429*WProd2442*WProd2498))/WProd2487)
        amp(31,5,3) = Prp255*
     -    (Prp258*Prp259*((-32*Alfas*Cpl122*WProd2500)/WProd2499 - 
     -         (32*Alfas*Cpl122*
     -            (-(WProd2440*WProd2490) + WProd2441*WProd2490 - 
     -              WProd2442*WProd2450*WProd2490 - 
     -              WProd2429*WProd2492 + WProd2440*WProd2494 - 
     -              WProd2441*WProd2494 + 
     -              WProd2442*WProd2450*WProd2494 + 
     -              WProd2450*WProd2495 + WProd2429*WProd2496 + 
     -              WProd2440*WProd2497 - WProd2441*WProd2497 + 
     -              WProd2442*WProd2450*WProd2497 + 
     -              WProd2440*WProd2498 - WProd2441*WProd2498 + 
     -              WProd2442*WProd2450*WProd2498 + WProd2501 - 
     -              2*WProd2426*WProd2502 - 
     -              2*WProd2427*WProd2502 + 
     -              2*WProd2443*WProd2502 + 
     -              2*WProd2444*WProd2502 + 
     -              2*WProd2445*WProd2502 + 
     -              2*WProd2446*WProd2502 + 
     -              2*WProd2447*WProd2502 + 
     -              2*WProd2449*WProd2502 - 
     -              2*WProd2451*WProd2502 - 
     -              2*WProd2452*WProd2502 + 
     -              2*WProd2426*WProd2506 + 
     -              2*WProd2427*WProd2506 - 
     -              2*WProd2443*WProd2506 - 
     -              2*WProd2444*WProd2506 - 
     -              2*WProd2445*WProd2506 - 
     -              2*WProd2446*WProd2506 - 
     -              2*WProd2447*WProd2506 - 
     -              2*WProd2449*WProd2506 + 
     -              2*WProd2451*WProd2506 + 
     -              2*WProd2452*WProd2506 - WProd2507))/WProd2487)+
     -        Prp256*Prp258*
     -       ((-32*Alfas*Cpl122*WProd2500)/WProd2499 - 
     -         (32*Alfas*Cpl122*
     -            (WProd2450*WProd2488 - WProd2440*WProd2490 + 
     -              WProd2441*WProd2490 - 
     -              WProd2442*WProd2450*WProd2490 - 
     -              WProd2450*WProd2491 - WProd2429*WProd2492 + 
     -              WProd2429*WProd2493 + WProd2440*WProd2494 - 
     -              WProd2441*WProd2494 + 
     -              WProd2442*WProd2450*WProd2494 + 
     -              WProd2450*WProd2495 + WProd2429*WProd2496 + 
     -              WProd2501 - 2*WProd2426*WProd2502 - 
     -              2*WProd2427*WProd2502 + 
     -              2*WProd2443*WProd2502 + 
     -              2*WProd2444*WProd2502 + 
     -              2*WProd2445*WProd2502 + 
     -              2*WProd2446*WProd2502 + 
     -              2*WProd2447*WProd2502 + WProd2449*WProd2502 - 
     -              2*WProd2451*WProd2502 - 
     -              2*WProd2452*WProd2502 - WProd2503 + 
     -              2*WProd2426*WProd2504 + 
     -              2*WProd2427*WProd2504 - 
     -              2*WProd2443*WProd2504 - 
     -              2*WProd2444*WProd2504 - 
     -              2*WProd2445*WProd2504 - 
     -              2*WProd2446*WProd2504 - 
     -              2*WProd2447*WProd2504 - WProd2449*WProd2504 + 
     -              2*WProd2451*WProd2504 + 
     -              2*WProd2452*WProd2504 - WProd2505 + 
     -              2*WProd2426*WProd2506 + 
     -              2*WProd2427*WProd2506 - 
     -              2*WProd2443*WProd2506 - 
     -              2*WProd2444*WProd2506 - 
     -              2*WProd2445*WProd2506 - 
     -              2*WProd2446*WProd2506 - 
     -              2*WProd2447*WProd2506 - WProd2449*WProd2506 + 
     -              2*WProd2451*WProd2506 + 
     -              2*WProd2452*WProd2506 - WProd2507))/WProd2487)+
     -        (32*Alfas*Cpl122*Prp257*Prp259*
     -         (WProd2442*WProd2450*WProd2490 + 
     -           WProd2442*WProd2457*WProd2490 + 
     -           WProd2442*WProd2459*WProd2490 - 
     -           WProd2442*WProd2450*WProd2494 - 
     -           WProd2442*WProd2457*WProd2494 - 
     -           WProd2442*WProd2459*WProd2494 - 
     -           WProd2450*WProd2495 - WProd2457*WProd2495 - 
     -           WProd2459*WProd2495 - WProd2429*WProd2496 - 
     -           WProd2442*WProd2450*WProd2497 - 
     -           WProd2442*WProd2457*WProd2497 - 
     -           WProd2442*WProd2459*WProd2497 - 
     -           WProd2442*WProd2450*WProd2498 - 
     -           WProd2442*WProd2457*WProd2498 - 
     -           WProd2442*WProd2459*WProd2498 - 
     -           2*WProd2426*WProd2506 - WProd2427*WProd2506 + 
     -           2*WProd2443*WProd2506 + WProd2444*WProd2506 + 
     -           2*WProd2445*WProd2506 + WProd2446*WProd2506 + 
     -           2*WProd2447*WProd2506 + WProd2449*WProd2506 - 
     -           WProd2451*WProd2506 - WProd2452*WProd2506 + 
     -           WProd2507))/WProd2487) + amp(31,5,3)

#ifdef DEBUG
	print *, 'amp(31,5,3) =', amp(31,5,3)
#endif

        amp(31,5,10) = Prp254*Prp255*
     -   ((-16*Alfas*Cpl122*Prp256*
     -        (2*WProd2470*WProd2488 - 2*WProd2472*WProd2490 - 
     -          2*WProd2473*WProd2490 + 2*WProd2474*WProd2490 + 
     -          2*WProd2475*WProd2490 - 
     -          2*WProd2442*WProd2476*WProd2490 - 
     -          2*WProd2477*WProd2490 - 2*WProd2478*WProd2490 - 
     -          2*WProd2470*WProd2491 - 2*WProd2470*WProd2492 + 
     -          2*WProd2470*WProd2493 + 
     -          2*WProd2442*WProd2470*WProd2494 + 
     -          2*WProd2470*WProd2495 + 2*WProd2470*WProd2496 + 
     -          2*WProd2433*WProd2509 - 2*WProd2437*WProd2509 + 
     -          4*WProd2426*WProd2442*WProd2509 + 
     -          3*WProd2427*WProd2442*WProd2509 - 
     -          3*WProd2442*WProd2443*WProd2509 - 
     -          3*WProd2442*WProd2444*WProd2509 - 
     -          3*WProd2442*WProd2445*WProd2509 - 
     -          3*WProd2442*WProd2446*WProd2509 - 
     -          4*WProd2442*WProd2447*WProd2509 - 
     -          2*WProd2442*WProd2449*WProd2509 + 
     -          3*WProd2442*WProd2451*WProd2509 + 
     -          3*WProd2442*WProd2452*WProd2509 - 
     -          2*WProd2453*WProd2509))/WProd2508 - 
     -     (16*Alfas*Cpl122*Prp257*
     -        (-2*WProd2442*WProd2476*WProd2490 - 
     -          2*WProd2477*WProd2490 - 2*WProd2478*WProd2490 - 
     -          2*WProd2442*WProd2479*WProd2490 - 
     -          2*WProd2442*WProd2480*WProd2490 + 
     -          2*WProd2442*WProd2470*WProd2494 + 
     -          2*WProd2470*WProd2495 + 2*WProd2470*WProd2496 + 
     -          2*WProd2442*WProd2470*WProd2497 + 
     -          2*WProd2442*WProd2470*WProd2498 + 
     -          4*WProd2426*WProd2442*WProd2509 + 
     -          3*WProd2427*WProd2442*WProd2509 - 
     -          4*WProd2442*WProd2443*WProd2509 - 
     -          2*WProd2442*WProd2444*WProd2509 - 
     -          4*WProd2442*WProd2445*WProd2509 - 
     -          2*WProd2442*WProd2446*WProd2509 - 
     -          4*WProd2442*WProd2447*WProd2509 - 
     -          2*WProd2442*WProd2449*WProd2509 + 
     -          3*WProd2442*WProd2451*WProd2509 + 
     -          3*WProd2442*WProd2452*WProd2509 - 
     -          2*WProd2453*WProd2509))/WProd2508)
        amp(31,5,10) = Prp255*
     -    (Prp258*Prp259*((-32*Alfas*Cpl122*WProd2511)/WProd2510 - 
     -         (32*Alfas*Cpl122*
     -            (-(WProd2473*WProd2490) + WProd2475*WProd2490 - 
     -              WProd2442*WProd2476*WProd2490 - 
     -              WProd2477*WProd2490 - WProd2470*WProd2492 + 
     -              WProd2473*WProd2494 - WProd2475*WProd2494 + 
     -              WProd2442*WProd2476*WProd2494 + 
     -              WProd2477*WProd2494 + WProd2476*WProd2495 + 
     -              WProd2470*WProd2496 + WProd2473*WProd2497 - 
     -              WProd2475*WProd2497 + 
     -              WProd2442*WProd2476*WProd2497 + 
     -              WProd2477*WProd2497 + WProd2473*WProd2498 - 
     -              WProd2475*WProd2498 + 
     -              WProd2442*WProd2476*WProd2498 + 
     -              WProd2477*WProd2498 + WProd2512 - 
     -              2*WProd2426*WProd2513 - 
     -              2*WProd2427*WProd2513 + 
     -              2*WProd2443*WProd2513 + 
     -              2*WProd2444*WProd2513 + 
     -              2*WProd2445*WProd2513 + 
     -              2*WProd2446*WProd2513 + 
     -              2*WProd2447*WProd2513 + 
     -              2*WProd2449*WProd2513 - 
     -              2*WProd2451*WProd2513 - 
     -              2*WProd2452*WProd2513 + 
     -              2*WProd2426*WProd2516 + 
     -              2*WProd2427*WProd2516 - 
     -              2*WProd2443*WProd2516 - 
     -              2*WProd2444*WProd2516 - 
     -              2*WProd2445*WProd2516 - 
     -              2*WProd2446*WProd2516 - 
     -              2*WProd2447*WProd2516 - 
     -              2*WProd2449*WProd2516 + 
     -              2*WProd2451*WProd2516 + 
     -              2*WProd2452*WProd2516 - WProd2517))/WProd2508)+
     -        Prp256*Prp258*
     -       ((-32*Alfas*Cpl122*WProd2511)/WProd2510 - 
     -         (32*Alfas*Cpl122*
     -            (WProd2476*WProd2488 - WProd2473*WProd2490 + 
     -              WProd2475*WProd2490 - 
     -              WProd2442*WProd2476*WProd2490 - 
     -              WProd2477*WProd2490 - WProd2476*WProd2491 - 
     -              WProd2470*WProd2492 + WProd2470*WProd2493 + 
     -              WProd2473*WProd2494 - WProd2475*WProd2494 + 
     -              WProd2442*WProd2476*WProd2494 + 
     -              WProd2477*WProd2494 + WProd2476*WProd2495 + 
     -              WProd2470*WProd2496 + WProd2512 - 
     -              2*WProd2426*WProd2513 - 
     -              2*WProd2427*WProd2513 + 
     -              2*WProd2443*WProd2513 + 
     -              2*WProd2444*WProd2513 + 
     -              2*WProd2445*WProd2513 + 
     -              2*WProd2446*WProd2513 + 
     -              2*WProd2447*WProd2513 + WProd2449*WProd2513 - 
     -              2*WProd2451*WProd2513 - 
     -              2*WProd2452*WProd2513 - WProd2514 + 
     -              2*WProd2426*WProd2515 + 
     -              2*WProd2427*WProd2515 - 
     -              2*WProd2443*WProd2515 - 
     -              2*WProd2444*WProd2515 - 
     -              2*WProd2445*WProd2515 - 
     -              2*WProd2446*WProd2515 - 
     -              2*WProd2447*WProd2515 - WProd2449*WProd2515 + 
     -              2*WProd2451*WProd2515 + 
     -              2*WProd2452*WProd2515 + 
     -              2*WProd2426*WProd2516 + 
     -              2*WProd2427*WProd2516 - 
     -              2*WProd2443*WProd2516 - 
     -              2*WProd2444*WProd2516 - 
     -              2*WProd2445*WProd2516 - 
     -              2*WProd2446*WProd2516 - 
     -              2*WProd2447*WProd2516 - WProd2449*WProd2516 + 
     -              2*WProd2451*WProd2516 + 
     -              2*WProd2452*WProd2516 - WProd2517))/WProd2508)+
     -        (32*Alfas*Cpl122*Prp257*Prp259*
     -         (WProd2442*WProd2476*WProd2490 + 
     -           WProd2477*WProd2490 + 
     -           WProd2442*WProd2479*WProd2490 + 
     -           WProd2442*WProd2480*WProd2490 - 
     -           WProd2442*WProd2476*WProd2494 - 
     -           WProd2477*WProd2494 - 
     -           WProd2442*WProd2479*WProd2494 - 
     -           WProd2442*WProd2480*WProd2494 - 
     -           WProd2476*WProd2495 - WProd2479*WProd2495 - 
     -           WProd2480*WProd2495 - WProd2470*WProd2496 - 
     -           WProd2442*WProd2476*WProd2497 - 
     -           WProd2477*WProd2497 - 
     -           WProd2442*WProd2479*WProd2497 - 
     -           WProd2442*WProd2480*WProd2497 - 
     -           WProd2442*WProd2476*WProd2498 - 
     -           WProd2477*WProd2498 - 
     -           WProd2442*WProd2479*WProd2498 - 
     -           WProd2442*WProd2480*WProd2498 - 
     -           2*WProd2426*WProd2516 - WProd2427*WProd2516 + 
     -           2*WProd2443*WProd2516 + WProd2444*WProd2516 + 
     -           2*WProd2445*WProd2516 + WProd2446*WProd2516 + 
     -           2*WProd2447*WProd2516 + WProd2449*WProd2516 - 
     -           WProd2451*WProd2516 - WProd2452*WProd2516 + 
     -           WProd2517))/WProd2508) + amp(31,5,10)

#ifdef DEBUG
	print *, 'amp(31,5,10) =', amp(31,5,10)
#endif

        amp(31,6,11) = Prp254*Prp255*
     -   ((16*Alfas*Cpl122*Prp256*
     -        (2*WProd2429*WProd2430 - 2*WProd2429*WProd2431 + 
     -          2*WProd2432*WProd2433 - 2*WProd2434*WProd2435 - 
     -          2*WProd2429*WProd2436 - 2*WProd2432*WProd2437 + 
     -          2*WProd2429*WProd2438 + 2*WProd2434*WProd2439 - 
     -          2*WProd2434*WProd2440 + 2*WProd2434*WProd2441 + 
     -          4*WProd2426*WProd2432*WProd2442 + 
     -          3*WProd2427*WProd2432*WProd2442 - 
     -          3*WProd2432*WProd2442*WProd2443 - 
     -          3*WProd2432*WProd2442*WProd2444 - 
     -          3*WProd2432*WProd2442*WProd2445 - 
     -          3*WProd2432*WProd2442*WProd2446 - 
     -          4*WProd2432*WProd2442*WProd2447 + 
     -          2*WProd2429*WProd2442*WProd2448 - 
     -          2*WProd2432*WProd2442*WProd2449 - 
     -          2*WProd2434*WProd2442*WProd2450 + 
     -          3*WProd2432*WProd2442*WProd2451 + 
     -          3*WProd2432*WProd2442*WProd2452 - 
     -          2*WProd2432*WProd2453 + 2*WProd2429*WProd2454 - 
     -          2*WProd2434*WProd2455))/WProd2428 + 
     -     (16*Alfas*Cpl122*Prp257*
     -        (4*WProd2426*WProd2432*WProd2442 + 
     -          3*WProd2427*WProd2432*WProd2442 - 
     -          4*WProd2432*WProd2442*WProd2443 - 
     -          2*WProd2432*WProd2442*WProd2444 - 
     -          4*WProd2432*WProd2442*WProd2445 - 
     -          2*WProd2432*WProd2442*WProd2446 - 
     -          4*WProd2432*WProd2442*WProd2447 + 
     -          2*WProd2429*WProd2442*WProd2448 - 
     -          2*WProd2432*WProd2442*WProd2449 - 
     -          2*WProd2434*WProd2442*WProd2450 + 
     -          3*WProd2432*WProd2442*WProd2451 + 
     -          3*WProd2432*WProd2442*WProd2452 - 
     -          2*WProd2432*WProd2453 + 2*WProd2429*WProd2454 - 
     -          2*WProd2434*WProd2455 + 
     -          2*WProd2429*WProd2442*WProd2456 - 
     -          2*WProd2434*WProd2442*WProd2457 + 
     -          2*WProd2429*WProd2442*WProd2458 - 
     -          2*WProd2434*WProd2442*WProd2459))/WProd2428)
        amp(31,6,11) = Prp255*
     -    ((32*Alfas*Cpl122*Prp257*Prp261*
     -         (2*WProd2426*WProd2432*WProd2442 + 
     -           2*WProd2427*WProd2432*WProd2442 - 
     -           2*WProd2432*WProd2442*WProd2443 - 
     -           WProd2432*WProd2442*WProd2444 - 
     -           2*WProd2432*WProd2442*WProd2445 - 
     -           WProd2432*WProd2442*WProd2446 - 
     -           2*WProd2432*WProd2442*WProd2447 + 
     -           WProd2429*WProd2442*WProd2448 - 
     -           WProd2432*WProd2442*WProd2449 - 
     -           WProd2442*WProd2448*WProd2450 + 
     -           2*WProd2432*WProd2442*WProd2451 + 
     -           2*WProd2432*WProd2442*WProd2452 - 
     -           WProd2432*WProd2453 - WProd2434*WProd2455 + 
     -           WProd2429*WProd2442*WProd2456 - 
     -           WProd2442*WProd2450*WProd2456 - 
     -           WProd2442*WProd2448*WProd2457 - 
     -           WProd2442*WProd2456*WProd2457 + 
     -           WProd2429*WProd2442*WProd2458 - 
     -           WProd2442*WProd2450*WProd2458 - 
     -           WProd2442*WProd2457*WProd2458 - 
     -           WProd2442*WProd2448*WProd2459 - 
     -           WProd2442*WProd2456*WProd2459 - 
     -           WProd2442*WProd2458*WProd2459 - 
     -           2*WProd2426*WProd2467 - 2*WProd2427*WProd2467 + 
     -           2*WProd2443*WProd2467 + WProd2444*WProd2467 + 
     -           2*WProd2445*WProd2467 + WProd2446*WProd2467 + 
     -           2*WProd2447*WProd2467 + WProd2449*WProd2467 - 
     -           2*WProd2451*WProd2467 - 2*WProd2452*WProd2467 + 
     -           WProd2468))/WProd2428 + 
     -      (32*Alfas*Cpl122*Prp260*Prp261*
     -         (WProd2429*WProd2430 - WProd2429*WProd2431 - 
     -           WProd2432*WProd2437 + WProd2434*WProd2439 + 
     -           3*WProd2426*WProd2432*WProd2442 + 
     -           2*WProd2427*WProd2432*WProd2442 - 
     -           2*WProd2432*WProd2442*WProd2443 - 
     -           2*WProd2432*WProd2442*WProd2444 - 
     -           2*WProd2432*WProd2442*WProd2445 - 
     -           2*WProd2432*WProd2442*WProd2446 - 
     -           2*WProd2432*WProd2442*WProd2447 + 
     -           WProd2429*WProd2442*WProd2448 - 
     -           2*WProd2432*WProd2442*WProd2449 - 
     -           WProd2430*WProd2450 + WProd2431*WProd2450 - 
     -           WProd2442*WProd2448*WProd2450 + 
     -           2*WProd2432*WProd2442*WProd2451 + 
     -           2*WProd2432*WProd2442*WProd2452 - 
     -           WProd2432*WProd2453 - WProd2434*WProd2455 - 
     -           WProd2430*WProd2457 + WProd2431*WProd2457 - 
     -           WProd2442*WProd2448*WProd2457 - 
     -           WProd2430*WProd2459 + WProd2431*WProd2459 - 
     -           WProd2442*WProd2448*WProd2459 - 
     -           3*WProd2426*WProd2464 - 2*WProd2427*WProd2464 + 
     -           2*WProd2443*WProd2464 + 2*WProd2444*WProd2464 + 
     -           2*WProd2445*WProd2464 + 2*WProd2446*WProd2464 + 
     -           2*WProd2447*WProd2464 + 2*WProd2449*WProd2464 - 
     -           2*WProd2451*WProd2464 - 2*WProd2452*WProd2464 + 
     -           WProd2465 + WProd2466 - 3*WProd2426*WProd2467 - 
     -           2*WProd2427*WProd2467 + 2*WProd2443*WProd2467 + 
     -           2*WProd2444*WProd2467 + 2*WProd2445*WProd2467 + 
     -           2*WProd2446*WProd2467 + 2*WProd2447*WProd2467 + 
     -           2*WProd2449*WProd2467 - 2*WProd2451*WProd2467 - 
     -           2*WProd2452*WProd2467 + WProd2468))/WProd2428 + 
     -      (32*Alfas*Cpl122*Prp256*Prp260*
     -         (WProd2429*WProd2430 - WProd2429*WProd2431 + 
     -           WProd2432*WProd2433 - WProd2434*WProd2435 - 
     -           WProd2432*WProd2437 + WProd2434*WProd2439 + 
     -           2*WProd2426*WProd2432*WProd2442 + 
     -           WProd2427*WProd2432*WProd2442 - 
     -           WProd2432*WProd2442*WProd2443 - 
     -           WProd2432*WProd2442*WProd2444 - 
     -           WProd2432*WProd2442*WProd2445 - 
     -           WProd2432*WProd2442*WProd2446 - 
     -           2*WProd2432*WProd2442*WProd2447 - 
     -           WProd2440*WProd2448 + WProd2441*WProd2448 + 
     -           WProd2429*WProd2442*WProd2448 - 
     -           WProd2432*WProd2442*WProd2449 - 
     -           WProd2430*WProd2450 + WProd2431*WProd2450 - 
     -           WProd2442*WProd2448*WProd2450 + 
     -           WProd2432*WProd2442*WProd2451 + 
     -           WProd2432*WProd2442*WProd2452 - 
     -           WProd2432*WProd2453 - WProd2434*WProd2455 + 
     -           2*WProd2426*WProd2462 + WProd2427*WProd2462 - 
     -           WProd2443*WProd2462 - WProd2444*WProd2462 - 
     -           WProd2445*WProd2462 - WProd2446*WProd2462 - 
     -           2*WProd2447*WProd2462 - WProd2449*WProd2462 + 
     -           WProd2451*WProd2462 + WProd2452*WProd2462 - 
     -           WProd2463 - 2*WProd2426*WProd2464 - 
     -           WProd2427*WProd2464 + WProd2443*WProd2464 + 
     -           WProd2444*WProd2464 + WProd2445*WProd2464 + 
     -           WProd2446*WProd2464 + 2*WProd2447*WProd2464 + 
     -           WProd2449*WProd2464 - WProd2451*WProd2464 - 
     -           WProd2452*WProd2464 + WProd2465 + WProd2466 - 
     -           2*WProd2426*WProd2467 - WProd2427*WProd2467 + 
     -           WProd2443*WProd2467 + WProd2444*WProd2467 + 
     -           WProd2445*WProd2467 + WProd2446*WProd2467 + 
     -           2*WProd2447*WProd2467 + WProd2449*WProd2467 - 
     -           WProd2451*WProd2467 - WProd2452*WProd2467 + 
     -           WProd2468))/WProd2428) + amp(31,6,11)

#ifdef DEBUG
	print *, 'amp(31,6,11) =', amp(31,6,11)
#endif

        amp(31,6,2) = Prp254*Prp255*
     -   ((16*Alfas*Cpl122*Prp256*
     -        (2*WProd2430*WProd2470 - 2*WProd2431*WProd2470 - 
     -          2*WProd2436*WProd2470 + 2*WProd2438*WProd2470 + 
     -          2*WProd2442*WProd2448*WProd2470 + 
     -          2*WProd2454*WProd2470 + 2*WProd2433*WProd2471 - 
     -          2*WProd2437*WProd2471 + 
     -          4*WProd2426*WProd2442*WProd2471 + 
     -          3*WProd2427*WProd2442*WProd2471 - 
     -          3*WProd2442*WProd2443*WProd2471 - 
     -          3*WProd2442*WProd2444*WProd2471 - 
     -          3*WProd2442*WProd2445*WProd2471 - 
     -          3*WProd2442*WProd2446*WProd2471 - 
     -          4*WProd2442*WProd2447*WProd2471 - 
     -          2*WProd2442*WProd2449*WProd2471 + 
     -          3*WProd2442*WProd2451*WProd2471 + 
     -          3*WProd2442*WProd2452*WProd2471 - 
     -          2*WProd2453*WProd2471 - 2*WProd2434*WProd2472 - 
     -          2*WProd2434*WProd2473 + 2*WProd2434*WProd2474 + 
     -          2*WProd2434*WProd2475 - 
     -          2*WProd2434*WProd2442*WProd2476 - 
     -          2*WProd2434*WProd2477 - 2*WProd2434*WProd2478))/
     -      WProd2469 + (16*Alfas*Cpl122*Prp257*
     -        (2*WProd2442*WProd2448*WProd2470 + 
     -          2*WProd2454*WProd2470 + 
     -          2*WProd2442*WProd2456*WProd2470 + 
     -          2*WProd2442*WProd2458*WProd2470 + 
     -          4*WProd2426*WProd2442*WProd2471 + 
     -          3*WProd2427*WProd2442*WProd2471 - 
     -          4*WProd2442*WProd2443*WProd2471 - 
     -          2*WProd2442*WProd2444*WProd2471 - 
     -          4*WProd2442*WProd2445*WProd2471 - 
     -          2*WProd2442*WProd2446*WProd2471 - 
     -          4*WProd2442*WProd2447*WProd2471 - 
     -          2*WProd2442*WProd2449*WProd2471 + 
     -          3*WProd2442*WProd2451*WProd2471 + 
     -          3*WProd2442*WProd2452*WProd2471 - 
     -          2*WProd2453*WProd2471 - 
     -          2*WProd2434*WProd2442*WProd2476 - 
     -          2*WProd2434*WProd2477 - 2*WProd2434*WProd2478 - 
     -          2*WProd2434*WProd2442*WProd2479 - 
     -          2*WProd2434*WProd2442*WProd2480))/WProd2469)
        amp(31,6,2) = Prp255*
     -    ((32*Alfas*Cpl122*Prp257*Prp261*
     -         (WProd2442*WProd2448*WProd2470 + 
     -           WProd2442*WProd2456*WProd2470 + 
     -           WProd2442*WProd2458*WProd2470 + 
     -           2*WProd2426*WProd2442*WProd2471 + 
     -           2*WProd2427*WProd2442*WProd2471 - 
     -           2*WProd2442*WProd2443*WProd2471 - 
     -           WProd2442*WProd2444*WProd2471 - 
     -           2*WProd2442*WProd2445*WProd2471 - 
     -           WProd2442*WProd2446*WProd2471 - 
     -           2*WProd2442*WProd2447*WProd2471 - 
     -           WProd2442*WProd2449*WProd2471 + 
     -           2*WProd2442*WProd2451*WProd2471 + 
     -           2*WProd2442*WProd2452*WProd2471 - 
     -           WProd2453*WProd2471 - 
     -           WProd2442*WProd2448*WProd2476 - 
     -           WProd2442*WProd2456*WProd2476 - 
     -           WProd2442*WProd2458*WProd2476 - 
     -           WProd2448*WProd2477 - WProd2456*WProd2477 - 
     -           WProd2458*WProd2477 - WProd2434*WProd2478 - 
     -           WProd2442*WProd2448*WProd2479 - 
     -           WProd2442*WProd2456*WProd2479 - 
     -           WProd2442*WProd2458*WProd2479 - 
     -           WProd2442*WProd2448*WProd2480 - 
     -           WProd2442*WProd2456*WProd2480 - 
     -           WProd2442*WProd2458*WProd2480 - 
     -           2*WProd2426*WProd2485 - 2*WProd2427*WProd2485 + 
     -           2*WProd2443*WProd2485 + WProd2444*WProd2485 + 
     -           2*WProd2445*WProd2485 + WProd2446*WProd2485 + 
     -           2*WProd2447*WProd2485 + WProd2449*WProd2485 - 
     -           2*WProd2451*WProd2485 - 2*WProd2452*WProd2485 + 
     -           WProd2486))/WProd2469 + 
     -      (32*Alfas*Cpl122*Prp260*Prp261*
     -         (WProd2430*WProd2470 - WProd2431*WProd2470 + 
     -           WProd2442*WProd2448*WProd2470 - 
     -           WProd2437*WProd2471 + 
     -           3*WProd2426*WProd2442*WProd2471 + 
     -           2*WProd2427*WProd2442*WProd2471 - 
     -           2*WProd2442*WProd2443*WProd2471 - 
     -           2*WProd2442*WProd2444*WProd2471 - 
     -           2*WProd2442*WProd2445*WProd2471 - 
     -           2*WProd2442*WProd2446*WProd2471 - 
     -           2*WProd2442*WProd2447*WProd2471 - 
     -           2*WProd2442*WProd2449*WProd2471 + 
     -           2*WProd2442*WProd2451*WProd2471 + 
     -           2*WProd2442*WProd2452*WProd2471 - 
     -           WProd2453*WProd2471 + WProd2434*WProd2474 - 
     -           WProd2430*WProd2476 + WProd2431*WProd2476 - 
     -           WProd2442*WProd2448*WProd2476 - 
     -           WProd2448*WProd2477 - WProd2434*WProd2478 - 
     -           WProd2430*WProd2479 + WProd2431*WProd2479 - 
     -           WProd2442*WProd2448*WProd2479 - 
     -           WProd2430*WProd2480 + WProd2431*WProd2480 - 
     -           WProd2442*WProd2448*WProd2480 + WProd2483 - 
     -           3*WProd2426*WProd2484 - 2*WProd2427*WProd2484 + 
     -           2*WProd2443*WProd2484 + 2*WProd2444*WProd2484 + 
     -           2*WProd2445*WProd2484 + 2*WProd2446*WProd2484 + 
     -           2*WProd2447*WProd2484 + 2*WProd2449*WProd2484 - 
     -           2*WProd2451*WProd2484 - 2*WProd2452*WProd2484 - 
     -           3*WProd2426*WProd2485 - 2*WProd2427*WProd2485 + 
     -           2*WProd2443*WProd2485 + 2*WProd2444*WProd2485 + 
     -           2*WProd2445*WProd2485 + 2*WProd2446*WProd2485 + 
     -           2*WProd2447*WProd2485 + 2*WProd2449*WProd2485 - 
     -           2*WProd2451*WProd2485 - 2*WProd2452*WProd2485 + 
     -           WProd2486))/WProd2469 + 
     -      (32*Alfas*Cpl122*Prp256*Prp260*
     -         (WProd2430*WProd2470 - WProd2431*WProd2470 + 
     -           WProd2442*WProd2448*WProd2470 + 
     -           WProd2433*WProd2471 - WProd2437*WProd2471 + 
     -           2*WProd2426*WProd2442*WProd2471 + 
     -           WProd2427*WProd2442*WProd2471 - 
     -           WProd2442*WProd2443*WProd2471 - 
     -           WProd2442*WProd2444*WProd2471 - 
     -           WProd2442*WProd2445*WProd2471 - 
     -           WProd2442*WProd2446*WProd2471 - 
     -           2*WProd2442*WProd2447*WProd2471 - 
     -           WProd2442*WProd2449*WProd2471 + 
     -           WProd2442*WProd2451*WProd2471 + 
     -           WProd2442*WProd2452*WProd2471 - 
     -           WProd2453*WProd2471 - WProd2434*WProd2472 - 
     -           WProd2448*WProd2473 + WProd2434*WProd2474 + 
     -           WProd2448*WProd2475 - WProd2430*WProd2476 + 
     -           WProd2431*WProd2476 - 
     -           WProd2442*WProd2448*WProd2476 - 
     -           WProd2448*WProd2477 - WProd2434*WProd2478 - 
     -           WProd2481 + 2*WProd2426*WProd2482 + 
     -           WProd2427*WProd2482 - WProd2443*WProd2482 - 
     -           WProd2444*WProd2482 - WProd2445*WProd2482 - 
     -           WProd2446*WProd2482 - 2*WProd2447*WProd2482 - 
     -           WProd2449*WProd2482 + WProd2451*WProd2482 + 
     -           WProd2452*WProd2482 + WProd2483 - 
     -           2*WProd2426*WProd2484 - WProd2427*WProd2484 + 
     -           WProd2443*WProd2484 + WProd2444*WProd2484 + 
     -           WProd2445*WProd2484 + WProd2446*WProd2484 + 
     -           2*WProd2447*WProd2484 + WProd2449*WProd2484 - 
     -           WProd2451*WProd2484 - WProd2452*WProd2484 - 
     -           2*WProd2426*WProd2485 - WProd2427*WProd2485 + 
     -           WProd2443*WProd2485 + WProd2444*WProd2485 + 
     -           WProd2445*WProd2485 + WProd2446*WProd2485 + 
     -           2*WProd2447*WProd2485 + WProd2449*WProd2485 - 
     -           WProd2451*WProd2485 - WProd2452*WProd2485 + 
     -           WProd2486))/WProd2469) + amp(31,6,2)

#ifdef DEBUG
	print *, 'amp(31,6,2) =', amp(31,6,2)
#endif

        amp(31,6,3) = Prp254*Prp255*
     -   ((16*Alfas*Cpl122*Prp256*
     -        (2*WProd2429*WProd2488 + 2*WProd2433*WProd2489 - 
     -          2*WProd2437*WProd2489 + 
     -          4*WProd2426*WProd2442*WProd2489 + 
     -          3*WProd2427*WProd2442*WProd2489 - 
     -          3*WProd2442*WProd2443*WProd2489 - 
     -          3*WProd2442*WProd2444*WProd2489 - 
     -          3*WProd2442*WProd2445*WProd2489 - 
     -          3*WProd2442*WProd2446*WProd2489 - 
     -          4*WProd2442*WProd2447*WProd2489 - 
     -          2*WProd2442*WProd2449*WProd2489 + 
     -          3*WProd2442*WProd2451*WProd2489 + 
     -          3*WProd2442*WProd2452*WProd2489 - 
     -          2*WProd2453*WProd2489 - 2*WProd2435*WProd2490 + 
     -          2*WProd2439*WProd2490 - 2*WProd2440*WProd2490 + 
     -          2*WProd2441*WProd2490 - 
     -          2*WProd2442*WProd2450*WProd2490 - 
     -          2*WProd2455*WProd2490 - 2*WProd2429*WProd2491 - 
     -          2*WProd2429*WProd2492 + 2*WProd2429*WProd2493 + 
     -          2*WProd2429*WProd2442*WProd2494 + 
     -          2*WProd2429*WProd2495 + 2*WProd2429*WProd2496))/
     -      WProd2487 + (16*Alfas*Cpl122*Prp257*
     -        (4*WProd2426*WProd2442*WProd2489 + 
     -          3*WProd2427*WProd2442*WProd2489 - 
     -          4*WProd2442*WProd2443*WProd2489 - 
     -          2*WProd2442*WProd2444*WProd2489 - 
     -          4*WProd2442*WProd2445*WProd2489 - 
     -          2*WProd2442*WProd2446*WProd2489 - 
     -          4*WProd2442*WProd2447*WProd2489 - 
     -          2*WProd2442*WProd2449*WProd2489 + 
     -          3*WProd2442*WProd2451*WProd2489 + 
     -          3*WProd2442*WProd2452*WProd2489 - 
     -          2*WProd2453*WProd2489 - 
     -          2*WProd2442*WProd2450*WProd2490 - 
     -          2*WProd2455*WProd2490 - 
     -          2*WProd2442*WProd2457*WProd2490 - 
     -          2*WProd2442*WProd2459*WProd2490 + 
     -          2*WProd2429*WProd2442*WProd2494 + 
     -          2*WProd2429*WProd2495 + 2*WProd2429*WProd2496 + 
     -          2*WProd2429*WProd2442*WProd2497 + 
     -          2*WProd2429*WProd2442*WProd2498))/WProd2487)
        amp(31,6,3) = Prp255*
     -    ((32*Alfas*Cpl122*Prp257*Prp261*
     -         (2*WProd2426*WProd2442*WProd2489 + 
     -           2*WProd2427*WProd2442*WProd2489 - 
     -           2*WProd2442*WProd2443*WProd2489 - 
     -           WProd2442*WProd2444*WProd2489 - 
     -           2*WProd2442*WProd2445*WProd2489 - 
     -           WProd2442*WProd2446*WProd2489 - 
     -           2*WProd2442*WProd2447*WProd2489 - 
     -           WProd2442*WProd2449*WProd2489 + 
     -           2*WProd2442*WProd2451*WProd2489 + 
     -           2*WProd2442*WProd2452*WProd2489 - 
     -           WProd2453*WProd2489 - WProd2455*WProd2490 + 
     -           WProd2429*WProd2442*WProd2494 - 
     -           WProd2442*WProd2450*WProd2494 - 
     -           WProd2442*WProd2457*WProd2494 - 
     -           WProd2442*WProd2459*WProd2494 + 
     -           WProd2429*WProd2495 - WProd2450*WProd2495 - 
     -           WProd2457*WProd2495 - WProd2459*WProd2495 + 
     -           WProd2429*WProd2442*WProd2497 - 
     -           WProd2442*WProd2450*WProd2497 - 
     -           WProd2442*WProd2457*WProd2497 - 
     -           WProd2442*WProd2459*WProd2497 + 
     -           WProd2429*WProd2442*WProd2498 - 
     -           WProd2442*WProd2450*WProd2498 - 
     -           WProd2442*WProd2457*WProd2498 - 
     -           WProd2442*WProd2459*WProd2498 - 
     -           2*WProd2426*WProd2506 - 2*WProd2427*WProd2506 + 
     -           2*WProd2443*WProd2506 + WProd2444*WProd2506 + 
     -           2*WProd2445*WProd2506 + WProd2446*WProd2506 + 
     -           2*WProd2447*WProd2506 + WProd2449*WProd2506 - 
     -           2*WProd2451*WProd2506 - 2*WProd2452*WProd2506 + 
     -           WProd2507))/WProd2487 + 
     -      (32*Alfas*Cpl122*Prp260*Prp261*
     -         (WProd2429*WProd2488 - WProd2450*WProd2488 - 
     -           WProd2457*WProd2488 - WProd2459*WProd2488 - 
     -           WProd2437*WProd2489 + 
     -           3*WProd2426*WProd2442*WProd2489 + 
     -           2*WProd2427*WProd2442*WProd2489 - 
     -           2*WProd2442*WProd2443*WProd2489 - 
     -           2*WProd2442*WProd2444*WProd2489 - 
     -           2*WProd2442*WProd2445*WProd2489 - 
     -           2*WProd2442*WProd2446*WProd2489 - 
     -           2*WProd2442*WProd2447*WProd2489 - 
     -           2*WProd2442*WProd2449*WProd2489 + 
     -           2*WProd2442*WProd2451*WProd2489 + 
     -           2*WProd2442*WProd2452*WProd2489 - 
     -           WProd2453*WProd2489 + WProd2439*WProd2490 - 
     -           WProd2455*WProd2490 - WProd2429*WProd2491 + 
     -           WProd2450*WProd2491 + WProd2457*WProd2491 + 
     -           WProd2459*WProd2491 + 
     -           WProd2429*WProd2442*WProd2494 - 
     -           WProd2442*WProd2450*WProd2494 - 
     -           WProd2442*WProd2457*WProd2494 - 
     -           WProd2442*WProd2459*WProd2494 + 
     -           WProd2429*WProd2495 - WProd2450*WProd2495 - 
     -           WProd2457*WProd2495 - WProd2459*WProd2495 + 
     -           WProd2503 - 3*WProd2426*WProd2504 - 
     -           2*WProd2427*WProd2504 + 2*WProd2443*WProd2504 + 
     -           2*WProd2444*WProd2504 + 2*WProd2445*WProd2504 + 
     -           2*WProd2446*WProd2504 + 2*WProd2447*WProd2504 + 
     -           2*WProd2449*WProd2504 - 2*WProd2451*WProd2504 - 
     -           2*WProd2452*WProd2504 + WProd2505 - 
     -           3*WProd2426*WProd2506 - 2*WProd2427*WProd2506 + 
     -           2*WProd2443*WProd2506 + 2*WProd2444*WProd2506 + 
     -           2*WProd2445*WProd2506 + 2*WProd2446*WProd2506 + 
     -           2*WProd2447*WProd2506 + 2*WProd2449*WProd2506 - 
     -           2*WProd2451*WProd2506 - 2*WProd2452*WProd2506 + 
     -           WProd2507))/WProd2487 + 
     -      Prp256*Prp260*((-32*Alfas*Cpl122*WProd2500)/
     -          WProd2499 + 
     -         (32*Alfas*Cpl122*
     -            (WProd2429*WProd2488 - WProd2450*WProd2488 + 
     -              WProd2433*WProd2489 - WProd2437*WProd2489 + 
     -              2*WProd2426*WProd2442*WProd2489 + 
     -              WProd2427*WProd2442*WProd2489 - 
     -              WProd2442*WProd2443*WProd2489 - 
     -              WProd2442*WProd2444*WProd2489 - 
     -              WProd2442*WProd2445*WProd2489 - 
     -              WProd2442*WProd2446*WProd2489 - 
     -              2*WProd2442*WProd2447*WProd2489 - 
     -              WProd2442*WProd2449*WProd2489 + 
     -              WProd2442*WProd2451*WProd2489 + 
     -              WProd2442*WProd2452*WProd2489 - 
     -              WProd2453*WProd2489 - WProd2435*WProd2490 + 
     -              WProd2439*WProd2490 - WProd2455*WProd2490 - 
     -              WProd2429*WProd2491 + WProd2450*WProd2491 - 
     -              WProd2440*WProd2494 + WProd2441*WProd2494 + 
     -              WProd2429*WProd2442*WProd2494 - 
     -              WProd2442*WProd2450*WProd2494 + 
     -              WProd2429*WProd2495 - WProd2450*WProd2495 - 
     -              WProd2501 + 2*WProd2426*WProd2502 + 
     -              WProd2427*WProd2502 - WProd2443*WProd2502 - 
     -              WProd2444*WProd2502 - WProd2445*WProd2502 - 
     -              WProd2446*WProd2502 - 2*WProd2447*WProd2502 - 
     -              WProd2449*WProd2502 + WProd2451*WProd2502 + 
     -              WProd2452*WProd2502 + WProd2503 - 
     -              2*WProd2426*WProd2504 - WProd2427*WProd2504 + 
     -              WProd2443*WProd2504 + WProd2444*WProd2504 + 
     -              WProd2445*WProd2504 + WProd2446*WProd2504 + 
     -              2*WProd2447*WProd2504 + WProd2449*WProd2504 - 
     -              WProd2451*WProd2504 - WProd2452*WProd2504 + 
     -              WProd2505 - 2*WProd2426*WProd2506 - 
     -              WProd2427*WProd2506 + WProd2443*WProd2506 + 
     -              WProd2444*WProd2506 + WProd2445*WProd2506 + 
     -              WProd2446*WProd2506 + 2*WProd2447*WProd2506 + 
     -              WProd2449*WProd2506 - WProd2451*WProd2506 - 
     -              WProd2452*WProd2506 + WProd2507))/WProd2487))+
     -     amp(31,6,3)

#ifdef DEBUG
	print *, 'amp(31,6,3) =', amp(31,6,3)
#endif

        amp(31,6,10) = Prp254*Prp255*
     -   ((16*Alfas*Cpl122*Prp256*
     -        (2*WProd2470*WProd2488 - 2*WProd2472*WProd2490 - 
     -          2*WProd2473*WProd2490 + 2*WProd2474*WProd2490 + 
     -          2*WProd2475*WProd2490 - 
     -          2*WProd2442*WProd2476*WProd2490 - 
     -          2*WProd2477*WProd2490 - 2*WProd2478*WProd2490 - 
     -          2*WProd2470*WProd2491 - 2*WProd2470*WProd2492 + 
     -          2*WProd2470*WProd2493 + 
     -          2*WProd2442*WProd2470*WProd2494 + 
     -          2*WProd2470*WProd2495 + 2*WProd2470*WProd2496 + 
     -          2*WProd2433*WProd2509 - 2*WProd2437*WProd2509 + 
     -          4*WProd2426*WProd2442*WProd2509 + 
     -          3*WProd2427*WProd2442*WProd2509 - 
     -          3*WProd2442*WProd2443*WProd2509 - 
     -          3*WProd2442*WProd2444*WProd2509 - 
     -          3*WProd2442*WProd2445*WProd2509 - 
     -          3*WProd2442*WProd2446*WProd2509 - 
     -          4*WProd2442*WProd2447*WProd2509 - 
     -          2*WProd2442*WProd2449*WProd2509 + 
     -          3*WProd2442*WProd2451*WProd2509 + 
     -          3*WProd2442*WProd2452*WProd2509 - 
     -          2*WProd2453*WProd2509))/WProd2508 + 
     -     (16*Alfas*Cpl122*Prp257*
     -        (-2*WProd2442*WProd2476*WProd2490 - 
     -          2*WProd2477*WProd2490 - 2*WProd2478*WProd2490 - 
     -          2*WProd2442*WProd2479*WProd2490 - 
     -          2*WProd2442*WProd2480*WProd2490 + 
     -          2*WProd2442*WProd2470*WProd2494 + 
     -          2*WProd2470*WProd2495 + 2*WProd2470*WProd2496 + 
     -          2*WProd2442*WProd2470*WProd2497 + 
     -          2*WProd2442*WProd2470*WProd2498 + 
     -          4*WProd2426*WProd2442*WProd2509 + 
     -          3*WProd2427*WProd2442*WProd2509 - 
     -          4*WProd2442*WProd2443*WProd2509 - 
     -          2*WProd2442*WProd2444*WProd2509 - 
     -          4*WProd2442*WProd2445*WProd2509 - 
     -          2*WProd2442*WProd2446*WProd2509 - 
     -          4*WProd2442*WProd2447*WProd2509 - 
     -          2*WProd2442*WProd2449*WProd2509 + 
     -          3*WProd2442*WProd2451*WProd2509 + 
     -          3*WProd2442*WProd2452*WProd2509 - 
     -          2*WProd2453*WProd2509))/WProd2508)
        amp(31,6,10) = Prp255*
     -    ((32*Alfas*Cpl122*Prp257*Prp261*
     -         (-(WProd2478*WProd2490) + 
     -           WProd2442*WProd2470*WProd2494 - 
     -           WProd2442*WProd2476*WProd2494 - 
     -           WProd2477*WProd2494 - 
     -           WProd2442*WProd2479*WProd2494 - 
     -           WProd2442*WProd2480*WProd2494 + 
     -           WProd2470*WProd2495 - WProd2476*WProd2495 - 
     -           WProd2479*WProd2495 - WProd2480*WProd2495 + 
     -           WProd2442*WProd2470*WProd2497 - 
     -           WProd2442*WProd2476*WProd2497 - 
     -           WProd2477*WProd2497 - 
     -           WProd2442*WProd2479*WProd2497 - 
     -           WProd2442*WProd2480*WProd2497 + 
     -           WProd2442*WProd2470*WProd2498 - 
     -           WProd2442*WProd2476*WProd2498 - 
     -           WProd2477*WProd2498 - 
     -           WProd2442*WProd2479*WProd2498 - 
     -           WProd2442*WProd2480*WProd2498 + 
     -           2*WProd2426*WProd2442*WProd2509 + 
     -           2*WProd2427*WProd2442*WProd2509 - 
     -           2*WProd2442*WProd2443*WProd2509 - 
     -           WProd2442*WProd2444*WProd2509 - 
     -           2*WProd2442*WProd2445*WProd2509 - 
     -           WProd2442*WProd2446*WProd2509 - 
     -           2*WProd2442*WProd2447*WProd2509 - 
     -           WProd2442*WProd2449*WProd2509 + 
     -           2*WProd2442*WProd2451*WProd2509 + 
     -           2*WProd2442*WProd2452*WProd2509 - 
     -           WProd2453*WProd2509 - 2*WProd2426*WProd2516 - 
     -           2*WProd2427*WProd2516 + 2*WProd2443*WProd2516 + 
     -           WProd2444*WProd2516 + 2*WProd2445*WProd2516 + 
     -           WProd2446*WProd2516 + 2*WProd2447*WProd2516 + 
     -           WProd2449*WProd2516 - 2*WProd2451*WProd2516 - 
     -           2*WProd2452*WProd2516 + WProd2517))/WProd2508 + 
     -      (32*Alfas*Cpl122*Prp260*Prp261*
     -         (WProd2470*WProd2488 - WProd2476*WProd2488 - 
     -           WProd2479*WProd2488 - WProd2480*WProd2488 + 
     -           WProd2474*WProd2490 - WProd2478*WProd2490 - 
     -           WProd2470*WProd2491 + WProd2476*WProd2491 + 
     -           WProd2479*WProd2491 + WProd2480*WProd2491 + 
     -           WProd2442*WProd2470*WProd2494 - 
     -           WProd2442*WProd2476*WProd2494 - 
     -           WProd2477*WProd2494 - 
     -           WProd2442*WProd2479*WProd2494 - 
     -           WProd2442*WProd2480*WProd2494 + 
     -           WProd2470*WProd2495 - WProd2476*WProd2495 - 
     -           WProd2479*WProd2495 - WProd2480*WProd2495 - 
     -           WProd2437*WProd2509 + 
     -           3*WProd2426*WProd2442*WProd2509 + 
     -           2*WProd2427*WProd2442*WProd2509 - 
     -           2*WProd2442*WProd2443*WProd2509 - 
     -           2*WProd2442*WProd2444*WProd2509 - 
     -           2*WProd2442*WProd2445*WProd2509 - 
     -           2*WProd2442*WProd2446*WProd2509 - 
     -           2*WProd2442*WProd2447*WProd2509 - 
     -           2*WProd2442*WProd2449*WProd2509 + 
     -           2*WProd2442*WProd2451*WProd2509 + 
     -           2*WProd2442*WProd2452*WProd2509 - 
     -           WProd2453*WProd2509 + WProd2514 - 
     -           3*WProd2426*WProd2515 - 2*WProd2427*WProd2515 + 
     -           2*WProd2443*WProd2515 + 2*WProd2444*WProd2515 + 
     -           2*WProd2445*WProd2515 + 2*WProd2446*WProd2515 + 
     -           2*WProd2447*WProd2515 + 2*WProd2449*WProd2515 - 
     -           2*WProd2451*WProd2515 - 2*WProd2452*WProd2515 - 
     -           3*WProd2426*WProd2516 - 2*WProd2427*WProd2516 + 
     -           2*WProd2443*WProd2516 + 2*WProd2444*WProd2516 + 
     -           2*WProd2445*WProd2516 + 2*WProd2446*WProd2516 + 
     -           2*WProd2447*WProd2516 + 2*WProd2449*WProd2516 - 
     -           2*WProd2451*WProd2516 - 2*WProd2452*WProd2516 + 
     -           WProd2517))/WProd2508 + 
     -      Prp256*Prp260*((-32*Alfas*Cpl122*WProd2511)/
     -          WProd2510 + 
     -         (32*Alfas*Cpl122*
     -            (WProd2470*WProd2488 - WProd2476*WProd2488 - 
     -              WProd2472*WProd2490 + WProd2474*WProd2490 - 
     -              WProd2478*WProd2490 - WProd2470*WProd2491 + 
     -              WProd2476*WProd2491 + 
     -              WProd2442*WProd2470*WProd2494 - 
     -              WProd2473*WProd2494 + WProd2475*WProd2494 - 
     -              WProd2442*WProd2476*WProd2494 - 
     -              WProd2477*WProd2494 + WProd2470*WProd2495 - 
     -              WProd2476*WProd2495 + WProd2433*WProd2509 - 
     -              WProd2437*WProd2509 + 
     -              2*WProd2426*WProd2442*WProd2509 + 
     -              WProd2427*WProd2442*WProd2509 - 
     -              WProd2442*WProd2443*WProd2509 - 
     -              WProd2442*WProd2444*WProd2509 - 
     -              WProd2442*WProd2445*WProd2509 - 
     -              WProd2442*WProd2446*WProd2509 - 
     -              2*WProd2442*WProd2447*WProd2509 - 
     -              WProd2442*WProd2449*WProd2509 + 
     -              WProd2442*WProd2451*WProd2509 + 
     -              WProd2442*WProd2452*WProd2509 - 
     -              WProd2453*WProd2509 - WProd2512 + 
     -              2*WProd2426*WProd2513 + WProd2427*WProd2513 - 
     -              WProd2443*WProd2513 - WProd2444*WProd2513 - 
     -              WProd2445*WProd2513 - WProd2446*WProd2513 - 
     -              2*WProd2447*WProd2513 - WProd2449*WProd2513 + 
     -              WProd2451*WProd2513 + WProd2452*WProd2513 + 
     -              WProd2514 - 2*WProd2426*WProd2515 - 
     -              WProd2427*WProd2515 + WProd2443*WProd2515 + 
     -              WProd2444*WProd2515 + WProd2445*WProd2515 + 
     -              WProd2446*WProd2515 + 2*WProd2447*WProd2515 + 
     -              WProd2449*WProd2515 - WProd2451*WProd2515 - 
     -              WProd2452*WProd2515 - 2*WProd2426*WProd2516 - 
     -              WProd2427*WProd2516 + WProd2443*WProd2516 + 
     -              WProd2444*WProd2516 + WProd2445*WProd2516 + 
     -              WProd2446*WProd2516 + 2*WProd2447*WProd2516 + 
     -              WProd2449*WProd2516 - WProd2451*WProd2516 - 
     -              WProd2452*WProd2516 + WProd2517))/WProd2508))+
     -     amp(31,6,10)

#ifdef DEBUG
	print *, 'amp(31,6,10) =', amp(31,6,10)
#endif

      endif

      end


