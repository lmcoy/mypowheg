      subroutine realQED30(amp,resnr)
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

        tmp1164 = Prp246*
     -    ((-4*Cpl2*GS*Prp247*
     -         (2*WProd2331*WProd2332 - 2*WProd2332*WProd2333 - 
     -           4*WProd2327*WProd2335 - 4*WProd2328*WProd2335 + 
     -           4*WProd2329*WProd2335 + 4*WProd2334*WProd2335 - 
     -           4*WProd2335*WProd2336 - 4*WProd2335*WProd2337 - 
     -           4*WProd2335*WProd2338 - 3*WProd2335*WProd2339 - 
     -           2*WProd2340 - 2*WProd2342 - 
     -           2*WProd2343*WProd2344 + 2*WProd2343*WProd2345 + 
     -           4*WProd2335*WProd2346 + 4*WProd2335*WProd2347 - 
     -           2*WProd2332*WProd2343*WProd2348 + 
     -           2*WProd2349*WProd2360 - 2*WProd2344*WProd2365 + 
     -           2*WProd2345*WProd2365 - 
     -           2*WProd2332*WProd2348*WProd2365))/WProd2330 + 
     -      (4*Cpl2*GS*Prp248*
     -         (-2*WProd2331*WProd2332 + 2*WProd2332*WProd2333 + 
     -           4*WProd2327*WProd2335 + 4*WProd2328*WProd2335 - 
     -           4*WProd2329*WProd2335 - 5*WProd2334*WProd2335 + 
     -           3*WProd2335*WProd2336 + 3*WProd2335*WProd2337 + 
     -           4*WProd2335*WProd2338 + 3*WProd2335*WProd2339 + 
     -           2*WProd2340 + 2*WProd2342 - 
     -           4*WProd2335*WProd2346 - 4*WProd2335*WProd2347 + 
     -           2*WProd2332*WProd2343*WProd2348 + 
     -           2*WProd2331*WProd2349 - 2*WProd2333*WProd2349 - 
     -           2*WProd2343*WProd2348*WProd2349 - 
     -           2*WProd2331*WProd2350 + 2*WProd2333*WProd2350 + 
     -           2*WProd2343*WProd2348*WProd2350 - 
     -           2*WProd2331*WProd2351 + 2*WProd2333*WProd2351 + 
     -           2*WProd2343*WProd2348*WProd2351 - 
     -           2*WProd2352*WProd2353 + 2*WProd2354*WProd2355 + 
     -           2*WProd2332*WProd2348*WProd2365 - 
     -           2*WProd2348*WProd2349*WProd2365 + 
     -           2*WProd2348*WProd2350*WProd2365 + 
     -           2*WProd2348*WProd2351*WProd2365))/WProd2330) + 
     -   5.333333333333333D0*
     -    (Cpl2*GS*Prp247*Prp250*
     -       (-(WProd2331*WProd2332) + WProd2332*WProd2333 + 
     -         2*WProd2327*WProd2335 + 2*WProd2328*WProd2335 - 
     -         2*WProd2329*WProd2335 - 2*WProd2334*WProd2335 + 
     -         2*WProd2335*WProd2336 + 2*WProd2335*WProd2337 + 
     -         2*WProd2335*WProd2338 + WProd2335*WProd2339 + 
     -         WProd2340 + WProd2342 - 2*WProd2335*WProd2346 - 
     -         2*WProd2335*WProd2347 + WProd2344*WProd2354 - 
     -         WProd2345*WProd2354 + 
     -         WProd2332*WProd2348*WProd2354 + 
     -         WProd2349*WProd2357 - 2*WProd2327*WProd2358 - 
     -         2*WProd2328*WProd2358 + 2*WProd2329*WProd2358 + 
     -         2*WProd2334*WProd2358 - 2*WProd2336*WProd2358 - 
     -         2*WProd2337*WProd2358 - 2*WProd2338*WProd2358 - 
     -         WProd2339*WProd2358 + 2*WProd2346*WProd2358 + 
     -         2*WProd2347*WProd2358 - WProd2359 - 
     -         WProd2349*WProd2360 - WProd2344*WProd2361 + 
     -         WProd2345*WProd2361 - 
     -         WProd2332*WProd2348*WProd2361 - 
     -         WProd2349*WProd2362 + 2*WProd2327*WProd2363 + 
     -         2*WProd2328*WProd2363 - 2*WProd2329*WProd2363 - 
     -         2*WProd2334*WProd2363 + 2*WProd2336*WProd2363 + 
     -         2*WProd2337*WProd2363 + 2*WProd2338*WProd2363 + 
     -         WProd2339*WProd2363 - 2*WProd2346*WProd2363 - 
     -         2*WProd2347*WProd2363 + WProd2364))/WProd2330 - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp247*Prp253*
     -       (-(WProd2343*WProd2344) + WProd2343*WProd2345 - 
     -         WProd2332*WProd2343*WProd2348 + 
     -         WProd2344*WProd2354 - WProd2345*WProd2354 + 
     -         WProd2332*WProd2348*WProd2354 + 
     -         WProd2349*WProd2357 - 2*WProd2327*WProd2358 - 
     -         2*WProd2328*WProd2358 + 2*WProd2329*WProd2358 + 
     -         2*WProd2334*WProd2358 - 2*WProd2336*WProd2358 - 
     -         2*WProd2337*WProd2358 - 2*WProd2338*WProd2358 - 
     -         2*WProd2339*WProd2358 + 2*WProd2346*WProd2358 + 
     -         2*WProd2347*WProd2358 - WProd2359 - 
     -         WProd2344*WProd2361 + WProd2345*WProd2361 - 
     -         WProd2332*WProd2348*WProd2361 - 
     -         WProd2349*WProd2362 + 2*WProd2327*WProd2363 + 
     -         2*WProd2328*WProd2363 - 2*WProd2329*WProd2363 - 
     -         2*WProd2334*WProd2363 + 2*WProd2336*WProd2363 + 
     -         2*WProd2337*WProd2363 + 2*WProd2338*WProd2363 + 
     -         2*WProd2339*WProd2363 - 2*WProd2346*WProd2363 - 
     -         2*WProd2347*WProd2363 + WProd2364 - 
     -         WProd2344*WProd2365 + WProd2345*WProd2365 - 
     -         WProd2332*WProd2348*WProd2365))/WProd2330 - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp252*Prp253*
     -       (-(WProd2332*WProd2343*WProd2348) - 
     -         WProd2343*WProd2348*WProd2350 - 
     -         WProd2343*WProd2348*WProd2351 + 
     -         WProd2332*WProd2348*WProd2354 + 
     -         WProd2348*WProd2350*WProd2354 + 
     -         WProd2348*WProd2351*WProd2354 - 
     -         WProd2332*WProd2348*WProd2361 - 
     -         WProd2348*WProd2350*WProd2361 - 
     -         WProd2348*WProd2351*WProd2361 - 
     -         WProd2349*WProd2362 + 2*WProd2327*WProd2363 + 
     -         2*WProd2328*WProd2363 - WProd2329*WProd2363 - 
     -         2*WProd2334*WProd2363 + WProd2336*WProd2363 + 
     -         WProd2337*WProd2363 + 2*WProd2338*WProd2363 + 
     -         WProd2339*WProd2363 - WProd2346*WProd2363 - 
     -         WProd2347*WProd2363 + WProd2364 - 
     -         WProd2332*WProd2348*WProd2365 - 
     -         WProd2348*WProd2350*WProd2365 - 
     -         WProd2348*WProd2351*WProd2365))/WProd2330
        tmp1164 = tmp1164 - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp248*Prp252*
     -       (-(WProd2332*WProd2343*WProd2348) + 
     -         WProd2343*WProd2348*WProd2349 - 
     -         WProd2343*WProd2348*WProd2350 - 
     -         WProd2343*WProd2348*WProd2351 - 
     -         2*WProd2327*WProd2348*WProd2352 - 
     -         2*WProd2328*WProd2348*WProd2352 + 
     -         2*WProd2329*WProd2348*WProd2352 + 
     -         2*WProd2334*WProd2348*WProd2352 - 
     -         WProd2336*WProd2348*WProd2352 - 
     -         WProd2337*WProd2348*WProd2352 - 
     -         2*WProd2338*WProd2348*WProd2352 - 
     -         WProd2339*WProd2348*WProd2352 + 
     -         2*WProd2346*WProd2348*WProd2352 + 
     -         2*WProd2347*WProd2348*WProd2352 - 
     -         WProd2332*WProd2348*WProd2361 + 
     -         WProd2348*WProd2349*WProd2361 - 
     -         WProd2348*WProd2350*WProd2361 - 
     -         WProd2348*WProd2351*WProd2361 + 
     -         2*WProd2327*WProd2363 + 2*WProd2328*WProd2363 - 
     -         2*WProd2329*WProd2363 - 2*WProd2334*WProd2363 + 
     -         WProd2336*WProd2363 + WProd2337*WProd2363 + 
     -         2*WProd2338*WProd2363 + WProd2339*WProd2363 - 
     -         2*WProd2346*WProd2363 - 2*WProd2347*WProd2363 + 
     -         WProd2364 - WProd2332*WProd2348*WProd2365 + 
     -         WProd2348*WProd2349*WProd2365 - 
     -         WProd2348*WProd2350*WProd2365 - 
     -         WProd2348*WProd2351*WProd2365 - 
     -         WProd2352*WProd2367 - WProd2354*WProd2368))/
     -     WProd2330 + Prp251*
     -    (5.333333333333333D0*
     -       (Cpl2*GS*Prp248*
     -          (-(WProd2331*WProd2332) + WProd2332*WProd2333 + 
     -            2*WProd2327*WProd2335 + 2*WProd2328*WProd2335 - 
     -            2*WProd2329*WProd2335 - 3*WProd2334*WProd2335 + 
     -            2*WProd2335*WProd2336 + 2*WProd2335*WProd2337 + 
     -            2*WProd2335*WProd2338 + 2*WProd2335*WProd2339 + 
     -            WProd2340 + WProd2342 - 2*WProd2335*WProd2346 - 
     -            2*WProd2335*WProd2347 + WProd2331*WProd2349 - 
     -            WProd2333*WProd2349 - WProd2331*WProd2350 + 
     -            WProd2333*WProd2350 - WProd2331*WProd2351 + 
     -            WProd2333*WProd2351 - 
     -            2*WProd2327*WProd2348*WProd2352 - 
     -            2*WProd2328*WProd2348*WProd2352 + 
     -            2*WProd2329*WProd2348*WProd2352 + 
     -            3*WProd2334*WProd2348*WProd2352 - 
     -            2*WProd2336*WProd2348*WProd2352 - 
     -            2*WProd2337*WProd2348*WProd2352 - 
     -            2*WProd2338*WProd2348*WProd2352 - 
     -            2*WProd2339*WProd2348*WProd2352 + 
     -            2*WProd2346*WProd2348*WProd2352 + 
     -            2*WProd2347*WProd2348*WProd2352 - 
     -            WProd2352*WProd2353 + WProd2354*WProd2355 - 
     -            WProd2332*WProd2348*WProd2361 + 
     -            WProd2348*WProd2349*WProd2361 - 
     -            WProd2348*WProd2350*WProd2361 - 
     -            WProd2348*WProd2351*WProd2361 + 
     -            2*WProd2327*WProd2363 + 2*WProd2328*WProd2363 - 
     -            2*WProd2329*WProd2363 - 3*WProd2334*WProd2363 + 
     -            2*WProd2336*WProd2363 + 2*WProd2337*WProd2363 + 
     -            2*WProd2338*WProd2363 + 2*WProd2339*WProd2363 - 
     -            2*WProd2346*WProd2363 - 2*WProd2347*WProd2363 + 
     -            WProd2364 - WProd2352*WProd2367 - 
     -            WProd2354*WProd2368))/WProd2330 + 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp250*
     -          (-(WProd2331*WProd2332) + WProd2332*WProd2333 + 
     -            WProd2327*WProd2335 + WProd2328*WProd2335 - 
     -            WProd2329*WProd2335 - 2*WProd2334*WProd2335 + 
     -            WProd2335*WProd2336 + WProd2335*WProd2337 + 
     -            2*WProd2335*WProd2338 + WProd2335*WProd2339 + 
     -            WProd2340 + WProd2342 - WProd2335*WProd2346 - 
     -            WProd2335*WProd2347 + WProd2331*WProd2349 - 
     -            WProd2333*WProd2349 - 
     -            WProd2327*WProd2348*WProd2352 - 
     -            WProd2328*WProd2348*WProd2352 + 
     -            WProd2329*WProd2348*WProd2352 + 
     -            2*WProd2334*WProd2348*WProd2352 - 
     -            WProd2336*WProd2348*WProd2352 - 
     -            WProd2337*WProd2348*WProd2352 - 
     -            2*WProd2338*WProd2348*WProd2352 - 
     -            WProd2339*WProd2348*WProd2352 + 
     -            WProd2346*WProd2348*WProd2352 + 
     -            WProd2347*WProd2348*WProd2352 - 
     -            WProd2352*WProd2353 + WProd2354*WProd2355 - 
     -            WProd2327*WProd2358 - WProd2328*WProd2358 + 
     -            WProd2329*WProd2358 + 2*WProd2334*WProd2358 - 
     -            WProd2336*WProd2358 - WProd2337*WProd2358 - 
     -            2*WProd2338*WProd2358 - WProd2339*WProd2358 + 
     -            WProd2346*WProd2358 + WProd2347*WProd2358 - 
     -            WProd2359 - WProd2344*WProd2361 + 
     -            WProd2345*WProd2361 - 
     -            WProd2332*WProd2348*WProd2361 + 
     -            WProd2348*WProd2349*WProd2361 + 
     -            WProd2327*WProd2363 + WProd2328*WProd2363 - 
     -            WProd2329*WProd2363 - 2*WProd2334*WProd2363 + 
     -            WProd2336*WProd2363 + WProd2337*WProd2363 + 
     -            2*WProd2338*WProd2363 + WProd2339*WProd2363 - 
     -            WProd2346*WProd2363 - WProd2347*WProd2363 + 
     -            WProd2364 - WProd2352*WProd2367 - 
     -            WProd2354*WProd2368 + WProd2352*WProd2369 - 
     -            WProd2354*WProd2370))/WProd2330)
        amp(30,4,11) = Prp249*tmp1164 + 
     -   Prp245*Prp246*((-8*Cpl2*GS*Prp247*
     -         (WProd2331*WProd2332 - WProd2332*WProd2333 - 
     -           4*WProd2327*WProd2335 - 4*WProd2328*WProd2335 + 
     -           4*WProd2329*WProd2335 + 4*WProd2334*WProd2335 - 
     -           4*WProd2335*WProd2336 - 4*WProd2335*WProd2337 - 
     -           4*WProd2335*WProd2338 - 4*WProd2335*WProd2339 - 
     -           WProd2340 + WProd2341 - WProd2342 - 
     -           WProd2343*WProd2344 + WProd2343*WProd2345 + 
     -           4*WProd2335*WProd2346 + 4*WProd2335*WProd2347 - 
     -           WProd2332*WProd2343*WProd2348))/WProd2330 + 
     -      (8*Cpl2*GS*Prp248*
     -         (-(WProd2331*WProd2332) + WProd2332*WProd2333 + 
     -           2*WProd2327*WProd2335 + 2*WProd2328*WProd2335 - 
     -           2*WProd2329*WProd2335 - 3*WProd2334*WProd2335 + 
     -           2*WProd2335*WProd2336 + 2*WProd2335*WProd2337 + 
     -           2*WProd2335*WProd2338 + 2*WProd2335*WProd2339 + 
     -           WProd2340 + WProd2342 - 2*WProd2335*WProd2346 - 
     -           2*WProd2335*WProd2347 + 
     -           WProd2332*WProd2343*WProd2348 + 
     -           WProd2331*WProd2349 - WProd2333*WProd2349 - 
     -           WProd2343*WProd2348*WProd2349 - 
     -           WProd2331*WProd2350 + WProd2333*WProd2350 + 
     -           WProd2343*WProd2348*WProd2350 - 
     -           WProd2331*WProd2351 + WProd2333*WProd2351 + 
     -           WProd2343*WProd2348*WProd2351 - 
     -           WProd2352*WProd2353 + WProd2354*WProd2355))/
     -       WProd2330)

#ifdef DEBUG
	print *, 'amp(30,4,11) =', amp(30,4,11)
#endif

        tmp1165 = Prp246*
     -    ((-4*Cpl2*GS*Prp247*
     -         (-2*WProd2343*WProd2372 - 2*WProd2365*WProd2372 + 
     -           2*WProd2331*WProd2373 - 2*WProd2333*WProd2373 - 
     -           2*WProd2343*WProd2348*WProd2373 - 
     -           2*WProd2348*WProd2365*WProd2373 - 
     -           4*WProd2327*WProd2374 - 4*WProd2328*WProd2374 + 
     -           4*WProd2329*WProd2374 + 4*WProd2334*WProd2374 - 
     -           4*WProd2336*WProd2374 - 4*WProd2337*WProd2374 - 
     -           4*WProd2338*WProd2374 - 3*WProd2339*WProd2374 + 
     -           4*WProd2346*WProd2374 + 4*WProd2347*WProd2374 + 
     -           2*WProd2343*WProd2375 + 2*WProd2365*WProd2375 - 
     -           2*WProd2343*WProd2376 - 2*WProd2365*WProd2376 + 
     -           2*WProd2360*WProd2377 - 2*WProd2381))/WProd2371 + 
     -      (4*Cpl2*GS*Prp248*
     -         (-2*WProd2331*WProd2373 + 2*WProd2333*WProd2373 + 
     -           2*WProd2343*WProd2348*WProd2373 + 
     -           2*WProd2348*WProd2365*WProd2373 + 
     -           4*WProd2327*WProd2374 + 4*WProd2328*WProd2374 - 
     -           4*WProd2329*WProd2374 - 5*WProd2334*WProd2374 + 
     -           3*WProd2336*WProd2374 + 3*WProd2337*WProd2374 + 
     -           4*WProd2338*WProd2374 + 3*WProd2339*WProd2374 - 
     -           4*WProd2346*WProd2374 - 4*WProd2347*WProd2374 + 
     -           2*WProd2343*WProd2376 + 2*WProd2365*WProd2376 + 
     -           2*WProd2331*WProd2377 - 2*WProd2333*WProd2377 - 
     -           2*WProd2343*WProd2348*WProd2377 - 
     -           2*WProd2348*WProd2365*WProd2377 - 
     -           2*WProd2331*WProd2378 + 2*WProd2333*WProd2378 + 
     -           2*WProd2343*WProd2348*WProd2378 + 
     -           2*WProd2348*WProd2365*WProd2378 - 
     -           2*WProd2331*WProd2379 + 2*WProd2333*WProd2379 + 
     -           2*WProd2343*WProd2348*WProd2379 + 
     -           2*WProd2348*WProd2365*WProd2379 - 
     -           2*WProd2353*WProd2380 + 2*WProd2381 + 
     -           2*WProd2354*WProd2382))/WProd2371) + 
     -   5.333333333333333D0*
     -    (Cpl2*GS*Prp247*Prp250*
     -       (WProd2354*WProd2372 - WProd2361*WProd2372 - 
     -         WProd2331*WProd2373 + WProd2333*WProd2373 + 
     -         WProd2348*WProd2354*WProd2373 - 
     -         WProd2348*WProd2361*WProd2373 + 
     -         2*WProd2327*WProd2374 + 2*WProd2328*WProd2374 - 
     -         2*WProd2329*WProd2374 - 2*WProd2334*WProd2374 + 
     -         2*WProd2336*WProd2374 + 2*WProd2337*WProd2374 + 
     -         2*WProd2338*WProd2374 + WProd2339*WProd2374 - 
     -         2*WProd2346*WProd2374 - 2*WProd2347*WProd2374 - 
     -         WProd2354*WProd2375 + WProd2361*WProd2375 + 
     -         WProd2354*WProd2376 - WProd2361*WProd2376 + 
     -         WProd2357*WProd2377 - WProd2360*WProd2377 - 
     -         WProd2362*WProd2377 + WProd2381 - WProd2383 - 
     -         2*WProd2327*WProd2384 - 2*WProd2328*WProd2384 + 
     -         2*WProd2329*WProd2384 + 2*WProd2334*WProd2384 - 
     -         2*WProd2336*WProd2384 - 2*WProd2337*WProd2384 - 
     -         2*WProd2338*WProd2384 - WProd2339*WProd2384 + 
     -         2*WProd2346*WProd2384 + 2*WProd2347*WProd2384 + 
     -         2*WProd2327*WProd2385 + 2*WProd2328*WProd2385 - 
     -         2*WProd2329*WProd2385 - 2*WProd2334*WProd2385 + 
     -         2*WProd2336*WProd2385 + 2*WProd2337*WProd2385 + 
     -         2*WProd2338*WProd2385 + WProd2339*WProd2385 - 
     -         2*WProd2346*WProd2385 - 2*WProd2347*WProd2385 + 
     -         WProd2386))/WProd2371 - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp247*Prp253*
     -       (-(WProd2343*WProd2372) + WProd2354*WProd2372 - 
     -         WProd2361*WProd2372 - WProd2365*WProd2372 - 
     -         WProd2343*WProd2348*WProd2373 + 
     -         WProd2348*WProd2354*WProd2373 - 
     -         WProd2348*WProd2361*WProd2373 - 
     -         WProd2348*WProd2365*WProd2373 + 
     -         WProd2343*WProd2375 - WProd2354*WProd2375 + 
     -         WProd2361*WProd2375 + WProd2365*WProd2375 - 
     -         WProd2343*WProd2376 + WProd2354*WProd2376 - 
     -         WProd2361*WProd2376 - WProd2365*WProd2376 + 
     -         WProd2357*WProd2377 - WProd2362*WProd2377 - 
     -         WProd2383 - 2*WProd2327*WProd2384 - 
     -         2*WProd2328*WProd2384 + 2*WProd2329*WProd2384 + 
     -         2*WProd2334*WProd2384 - 2*WProd2336*WProd2384 - 
     -         2*WProd2337*WProd2384 - 2*WProd2338*WProd2384 - 
     -         2*WProd2339*WProd2384 + 2*WProd2346*WProd2384 + 
     -         2*WProd2347*WProd2384 + 2*WProd2327*WProd2385 + 
     -         2*WProd2328*WProd2385 - 2*WProd2329*WProd2385 - 
     -         2*WProd2334*WProd2385 + 2*WProd2336*WProd2385 + 
     -         2*WProd2337*WProd2385 + 2*WProd2338*WProd2385 + 
     -         2*WProd2339*WProd2385 - 2*WProd2346*WProd2385 - 
     -         2*WProd2347*WProd2385 + WProd2386))/WProd2371
        tmp1165 = tmp1165 - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp252*Prp253*
     -       (-(WProd2343*WProd2348*WProd2373) + 
     -         WProd2348*WProd2354*WProd2373 - 
     -         WProd2348*WProd2361*WProd2373 - 
     -         WProd2348*WProd2365*WProd2373 - 
     -         WProd2343*WProd2376 + WProd2354*WProd2376 - 
     -         WProd2361*WProd2376 - WProd2365*WProd2376 - 
     -         WProd2362*WProd2377 - 
     -         WProd2343*WProd2348*WProd2378 + 
     -         WProd2348*WProd2354*WProd2378 - 
     -         WProd2348*WProd2361*WProd2378 - 
     -         WProd2348*WProd2365*WProd2378 - 
     -         WProd2343*WProd2348*WProd2379 + 
     -         WProd2348*WProd2354*WProd2379 - 
     -         WProd2348*WProd2361*WProd2379 - 
     -         WProd2348*WProd2365*WProd2379 + 
     -         2*WProd2327*WProd2385 + 2*WProd2328*WProd2385 - 
     -         WProd2329*WProd2385 - 2*WProd2334*WProd2385 + 
     -         WProd2336*WProd2385 + WProd2337*WProd2385 + 
     -         2*WProd2338*WProd2385 + WProd2339*WProd2385 - 
     -         WProd2346*WProd2385 - WProd2347*WProd2385 + 
     -         WProd2386))/WProd2371 - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp248*Prp252*
     -       (-(WProd2343*WProd2348*WProd2373) - 
     -         WProd2348*WProd2361*WProd2373 - 
     -         WProd2348*WProd2365*WProd2373 - 
     -         WProd2343*WProd2376 - WProd2361*WProd2376 - 
     -         WProd2365*WProd2376 + 
     -         WProd2343*WProd2348*WProd2377 + 
     -         WProd2348*WProd2361*WProd2377 + 
     -         WProd2348*WProd2365*WProd2377 - 
     -         WProd2343*WProd2348*WProd2378 - 
     -         WProd2348*WProd2361*WProd2378 - 
     -         WProd2348*WProd2365*WProd2378 - 
     -         WProd2343*WProd2348*WProd2379 - 
     -         WProd2348*WProd2361*WProd2379 - 
     -         WProd2348*WProd2365*WProd2379 - 
     -         2*WProd2327*WProd2348*WProd2380 - 
     -         2*WProd2328*WProd2348*WProd2380 + 
     -         2*WProd2329*WProd2348*WProd2380 + 
     -         2*WProd2334*WProd2348*WProd2380 - 
     -         WProd2336*WProd2348*WProd2380 - 
     -         WProd2337*WProd2348*WProd2380 - 
     -         2*WProd2338*WProd2348*WProd2380 - 
     -         WProd2339*WProd2348*WProd2380 + 
     -         2*WProd2346*WProd2348*WProd2380 + 
     -         2*WProd2347*WProd2348*WProd2380 - 
     -         WProd2367*WProd2380 + 2*WProd2327*WProd2385 + 
     -         2*WProd2328*WProd2385 - 2*WProd2329*WProd2385 - 
     -         2*WProd2334*WProd2385 + WProd2336*WProd2385 + 
     -         WProd2337*WProd2385 + 2*WProd2338*WProd2385 + 
     -         WProd2339*WProd2385 - 2*WProd2346*WProd2385 - 
     -         2*WProd2347*WProd2385 + WProd2386 - 
     -         WProd2354*WProd2387))/WProd2371
        tmp1165 = tmp1165 + 
     -   Prp251*(5.333333333333333D0*
     -       (Cpl2*GS*Prp248*
     -          (-(WProd2331*WProd2373) + WProd2333*WProd2373 - 
     -            WProd2348*WProd2361*WProd2373 + 
     -            2*WProd2327*WProd2374 + 2*WProd2328*WProd2374 - 
     -            2*WProd2329*WProd2374 - 3*WProd2334*WProd2374 + 
     -            2*WProd2336*WProd2374 + 2*WProd2337*WProd2374 + 
     -            2*WProd2338*WProd2374 + 2*WProd2339*WProd2374 - 
     -            2*WProd2346*WProd2374 - 2*WProd2347*WProd2374 - 
     -            WProd2361*WProd2376 + WProd2331*WProd2377 - 
     -            WProd2333*WProd2377 + 
     -            WProd2348*WProd2361*WProd2377 - 
     -            WProd2331*WProd2378 + WProd2333*WProd2378 - 
     -            WProd2348*WProd2361*WProd2378 - 
     -            WProd2331*WProd2379 + WProd2333*WProd2379 - 
     -            WProd2348*WProd2361*WProd2379 - 
     -            2*WProd2327*WProd2348*WProd2380 - 
     -            2*WProd2328*WProd2348*WProd2380 + 
     -            2*WProd2329*WProd2348*WProd2380 + 
     -            3*WProd2334*WProd2348*WProd2380 - 
     -            2*WProd2336*WProd2348*WProd2380 - 
     -            2*WProd2337*WProd2348*WProd2380 - 
     -            2*WProd2338*WProd2348*WProd2380 - 
     -            2*WProd2339*WProd2348*WProd2380 + 
     -            2*WProd2346*WProd2348*WProd2380 + 
     -            2*WProd2347*WProd2348*WProd2380 - 
     -            WProd2353*WProd2380 - WProd2367*WProd2380 + 
     -            WProd2381 + WProd2354*WProd2382 + 
     -            2*WProd2327*WProd2385 + 2*WProd2328*WProd2385 - 
     -            2*WProd2329*WProd2385 - 3*WProd2334*WProd2385 + 
     -            2*WProd2336*WProd2385 + 2*WProd2337*WProd2385 + 
     -            2*WProd2338*WProd2385 + 2*WProd2339*WProd2385 - 
     -            2*WProd2346*WProd2385 - 2*WProd2347*WProd2385 + 
     -            WProd2386 - WProd2354*WProd2387))/WProd2371 + 
     -      5.333333333333333D0*
     -       (Cpl2*GS*Prp250*
     -          (-(WProd2361*WProd2372) - WProd2331*WProd2373 + 
     -            WProd2333*WProd2373 - 
     -            WProd2348*WProd2361*WProd2373 + 
     -            WProd2327*WProd2374 + WProd2328*WProd2374 - 
     -            WProd2329*WProd2374 - 2*WProd2334*WProd2374 + 
     -            WProd2336*WProd2374 + WProd2337*WProd2374 + 
     -            2*WProd2338*WProd2374 + WProd2339*WProd2374 - 
     -            WProd2346*WProd2374 - WProd2347*WProd2374 + 
     -            WProd2361*WProd2375 - WProd2361*WProd2376 + 
     -            WProd2331*WProd2377 - WProd2333*WProd2377 + 
     -            WProd2348*WProd2361*WProd2377 - 
     -            WProd2327*WProd2348*WProd2380 - 
     -            WProd2328*WProd2348*WProd2380 + 
     -            WProd2329*WProd2348*WProd2380 + 
     -            2*WProd2334*WProd2348*WProd2380 - 
     -            WProd2336*WProd2348*WProd2380 - 
     -            WProd2337*WProd2348*WProd2380 - 
     -            2*WProd2338*WProd2348*WProd2380 - 
     -            WProd2339*WProd2348*WProd2380 + 
     -            WProd2346*WProd2348*WProd2380 + 
     -            WProd2347*WProd2348*WProd2380 - 
     -            WProd2353*WProd2380 - WProd2367*WProd2380 + 
     -            WProd2369*WProd2380 + WProd2381 + 
     -            WProd2354*WProd2382 - WProd2383 - 
     -            WProd2327*WProd2384 - WProd2328*WProd2384 + 
     -            WProd2329*WProd2384 + 2*WProd2334*WProd2384 - 
     -            WProd2336*WProd2384 - WProd2337*WProd2384 - 
     -            2*WProd2338*WProd2384 - WProd2339*WProd2384 + 
     -            WProd2346*WProd2384 + WProd2347*WProd2384 + 
     -            WProd2327*WProd2385 + WProd2328*WProd2385 - 
     -            WProd2329*WProd2385 - 2*WProd2334*WProd2385 + 
     -            WProd2336*WProd2385 + WProd2337*WProd2385 + 
     -            2*WProd2338*WProd2385 + WProd2339*WProd2385 - 
     -            WProd2346*WProd2385 - WProd2347*WProd2385 + 
     -            WProd2386 - WProd2354*WProd2387 - 
     -            WProd2354*WProd2388))/WProd2371)
        amp(30,4,2) = Prp249*tmp1165 + 
     -   Prp245*Prp246*((8*Cpl2*GS*Prp247*
     -         (WProd2343*WProd2372 - WProd2331*WProd2373 + 
     -           WProd2333*WProd2373 + 
     -           WProd2343*WProd2348*WProd2373 + 
     -           4*WProd2327*WProd2374 + 4*WProd2328*WProd2374 - 
     -           4*WProd2329*WProd2374 - 4*WProd2334*WProd2374 + 
     -           4*WProd2336*WProd2374 + 4*WProd2337*WProd2374 + 
     -           4*WProd2338*WProd2374 + 4*WProd2339*WProd2374 - 
     -           4*WProd2346*WProd2374 - 4*WProd2347*WProd2374 - 
     -           WProd2343*WProd2375 + WProd2343*WProd2376))/
     -       WProd2371 + (8*Cpl2*GS*Prp248*
     -         (-(WProd2331*WProd2373) + WProd2333*WProd2373 + 
     -           WProd2343*WProd2348*WProd2373 + 
     -           2*WProd2327*WProd2374 + 2*WProd2328*WProd2374 - 
     -           2*WProd2329*WProd2374 - 3*WProd2334*WProd2374 + 
     -           2*WProd2336*WProd2374 + 2*WProd2337*WProd2374 + 
     -           2*WProd2338*WProd2374 + 2*WProd2339*WProd2374 - 
     -           2*WProd2346*WProd2374 - 2*WProd2347*WProd2374 + 
     -           WProd2343*WProd2376 + WProd2331*WProd2377 - 
     -           WProd2333*WProd2377 - 
     -           WProd2343*WProd2348*WProd2377 - 
     -           WProd2331*WProd2378 + WProd2333*WProd2378 + 
     -           WProd2343*WProd2348*WProd2378 - 
     -           WProd2331*WProd2379 + WProd2333*WProd2379 + 
     -           WProd2343*WProd2348*WProd2379 - 
     -           WProd2353*WProd2380 + WProd2381 + 
     -           WProd2354*WProd2382))/WProd2371)

#ifdef DEBUG
	print *, 'amp(30,4,2) =', amp(30,4,2)
#endif

        tmp1166 = Prp247*Prp250*
     -    (-(5.333333333333333D0*(Cpl2*GS*WProd2402)/WProd2389) - 
     -      5.333333333333333D0*
     -       (Cpl2*GS*(WProd2332*WProd2392 - WProd2332*WProd2393 - 
     -            WProd2394 - 2*WProd2327*WProd2396 - 
     -            2*WProd2328*WProd2396 + 2*WProd2329*WProd2396 + 
     -            2*WProd2334*WProd2396 - 2*WProd2336*WProd2396 - 
     -            2*WProd2337*WProd2396 - 2*WProd2338*WProd2396 - 
     -            WProd2339*WProd2396 + 2*WProd2346*WProd2396 + 
     -            2*WProd2347*WProd2396 - WProd2397 - 
     -            WProd2344*WProd2401 + WProd2345*WProd2401 - 
     -            WProd2332*WProd2348*WProd2401 + WProd2403 - 
     -            WProd2349*WProd2404 + 2*WProd2327*WProd2405 + 
     -            2*WProd2328*WProd2405 - 2*WProd2329*WProd2405 - 
     -            2*WProd2334*WProd2405 + 2*WProd2336*WProd2405 + 
     -            2*WProd2337*WProd2405 + 2*WProd2338*WProd2405 + 
     -            WProd2339*WProd2405 - 2*WProd2346*WProd2405 - 
     -            2*WProd2347*WProd2405 + WProd2349*WProd2406 + 
     -            WProd2344*WProd2407 - WProd2345*WProd2407 + 
     -            WProd2332*WProd2348*WProd2407 + 
     -            WProd2332*WProd2408 + WProd2349*WProd2409 - 
     -            2*WProd2327*WProd2410 - 2*WProd2328*WProd2410 + 
     -            2*WProd2329*WProd2410 + 2*WProd2334*WProd2410 - 
     -            2*WProd2336*WProd2410 - 2*WProd2337*WProd2410 - 
     -            2*WProd2338*WProd2410 - WProd2339*WProd2410 + 
     -            2*WProd2346*WProd2410 + 2*WProd2347*WProd2410 - 
     -            WProd2411))/WProd2391) + 
     -   Prp251*(5.333333333333333D0*
     -       (Cpl2*GS*Prp248*
     -          (-(WProd2332*WProd2392) + WProd2349*WProd2392 - 
     -            WProd2350*WProd2392 - WProd2351*WProd2392 + 
     -            WProd2332*WProd2393 - WProd2349*WProd2393 + 
     -            WProd2350*WProd2393 + WProd2351*WProd2393 + 
     -            WProd2394 + 2*WProd2327*WProd2396 + 
     -            2*WProd2328*WProd2396 - 2*WProd2329*WProd2396 - 
     -            3*WProd2334*WProd2396 + 2*WProd2336*WProd2396 + 
     -            2*WProd2337*WProd2396 + 2*WProd2338*WProd2396 + 
     -            2*WProd2339*WProd2396 - 2*WProd2346*WProd2396 - 
     -            2*WProd2347*WProd2396 + WProd2397 - 
     -            2*WProd2327*WProd2348*WProd2400 - 
     -            2*WProd2328*WProd2348*WProd2400 + 
     -            2*WProd2329*WProd2348*WProd2400 + 
     -            3*WProd2334*WProd2348*WProd2400 - 
     -            2*WProd2336*WProd2348*WProd2400 - 
     -            2*WProd2337*WProd2348*WProd2400 - 
     -            2*WProd2338*WProd2348*WProd2400 - 
     -            2*WProd2339*WProd2348*WProd2400 + 
     -            2*WProd2346*WProd2348*WProd2400 + 
     -            2*WProd2347*WProd2348*WProd2400 - 
     -            WProd2353*WProd2400 - WProd2367*WProd2400 + 
     -            WProd2355*WProd2401 - WProd2368*WProd2401 - 
     -            WProd2332*WProd2348*WProd2407 + 
     -            WProd2348*WProd2349*WProd2407 - 
     -            WProd2348*WProd2350*WProd2407 - 
     -            WProd2348*WProd2351*WProd2407 - 
     -            WProd2332*WProd2408 + WProd2349*WProd2408 - 
     -            WProd2350*WProd2408 - WProd2351*WProd2408 + 
     -            2*WProd2327*WProd2410 + 2*WProd2328*WProd2410 - 
     -            2*WProd2329*WProd2410 - 3*WProd2334*WProd2410 + 
     -            2*WProd2336*WProd2410 + 2*WProd2337*WProd2410 + 
     -            2*WProd2338*WProd2410 + 2*WProd2339*WProd2410 - 
     -            2*WProd2346*WProd2410 - 2*WProd2347*WProd2410 + 
     -            WProd2411))/WProd2391 + 
     -      Prp250*(-(5.333333333333333D0*
     -            (Cpl2*GS*WProd2402)/WProd2389) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(-(WProd2332*WProd2392) + 
     -               WProd2349*WProd2392 + WProd2332*WProd2393 - 
     -               WProd2349*WProd2393 + WProd2394 + 
     -               WProd2327*WProd2396 + WProd2328*WProd2396 - 
     -               WProd2329*WProd2396 - 2*WProd2334*WProd2396 + 
     -               WProd2336*WProd2396 + WProd2337*WProd2396 + 
     -               2*WProd2338*WProd2396 + WProd2339*WProd2396 - 
     -               WProd2346*WProd2396 - WProd2347*WProd2396 + 
     -               WProd2397 - WProd2327*WProd2348*WProd2400 - 
     -               WProd2328*WProd2348*WProd2400 + 
     -               WProd2329*WProd2348*WProd2400 + 
     -               2*WProd2334*WProd2348*WProd2400 - 
     -               WProd2336*WProd2348*WProd2400 - 
     -               WProd2337*WProd2348*WProd2400 - 
     -               2*WProd2338*WProd2348*WProd2400 - 
     -               WProd2339*WProd2348*WProd2400 + 
     -               WProd2346*WProd2348*WProd2400 + 
     -               WProd2347*WProd2348*WProd2400 - 
     -               WProd2353*WProd2400 - WProd2367*WProd2400 + 
     -               WProd2369*WProd2400 + WProd2355*WProd2401 - 
     -               WProd2368*WProd2401 - WProd2370*WProd2401 - 
     -               WProd2403 - WProd2327*WProd2405 - 
     -               WProd2328*WProd2405 + WProd2329*WProd2405 + 
     -               2*WProd2334*WProd2405 - WProd2336*WProd2405 - 
     -               WProd2337*WProd2405 - 2*WProd2338*WProd2405 - 
     -               WProd2339*WProd2405 + WProd2346*WProd2405 + 
     -               WProd2347*WProd2405 - WProd2344*WProd2407 + 
     -               WProd2345*WProd2407 - 
     -               WProd2332*WProd2348*WProd2407 + 
     -               WProd2348*WProd2349*WProd2407 - 
     -               WProd2332*WProd2408 + WProd2349*WProd2408 + 
     -               WProd2327*WProd2410 + WProd2328*WProd2410 - 
     -               WProd2329*WProd2410 - 2*WProd2334*WProd2410 + 
     -               WProd2336*WProd2410 + WProd2337*WProd2410 + 
     -               2*WProd2338*WProd2410 + WProd2339*WProd2410 - 
     -               WProd2346*WProd2410 - WProd2347*WProd2410 + 
     -               WProd2411))/WProd2391))
        tmp1166 = tmp1166 + 
     -   Prp247*Prp253*(2.6666666666666665D0*
     -       (Cpl2*GS*WProd2402)/WProd2389 + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*(WProd2344*WProd2398 - WProd2345*WProd2398 + 
     -            WProd2332*WProd2348*WProd2398 - 
     -            WProd2344*WProd2401 + WProd2345*WProd2401 - 
     -            WProd2332*WProd2348*WProd2401 + WProd2403 - 
     -            WProd2349*WProd2404 + 2*WProd2327*WProd2405 + 
     -            2*WProd2328*WProd2405 - 2*WProd2329*WProd2405 - 
     -            2*WProd2334*WProd2405 + 2*WProd2336*WProd2405 + 
     -            2*WProd2337*WProd2405 + 2*WProd2338*WProd2405 + 
     -            2*WProd2339*WProd2405 - 2*WProd2346*WProd2405 - 
     -            2*WProd2347*WProd2405 + WProd2344*WProd2407 - 
     -            WProd2345*WProd2407 + 
     -            WProd2332*WProd2348*WProd2407 + 
     -            WProd2332*WProd2408 + WProd2349*WProd2409 - 
     -            2*WProd2327*WProd2410 - 2*WProd2328*WProd2410 + 
     -            2*WProd2329*WProd2410 + 2*WProd2334*WProd2410 - 
     -            2*WProd2336*WProd2410 - 2*WProd2337*WProd2410 - 
     -            2*WProd2338*WProd2410 - 2*WProd2339*WProd2410 + 
     -            2*WProd2346*WProd2410 + 2*WProd2347*WProd2410 - 
     -            WProd2411 + WProd2344*WProd2412 - 
     -            WProd2345*WProd2412 + 
     -            WProd2332*WProd2348*WProd2412))/WProd2391) - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp252*Prp253*
     -       (-(WProd2332*WProd2348*WProd2398) - 
     -         WProd2348*WProd2350*WProd2398 - 
     -         WProd2348*WProd2351*WProd2398 + 
     -         WProd2332*WProd2348*WProd2401 + 
     -         WProd2348*WProd2350*WProd2401 + 
     -         WProd2348*WProd2351*WProd2401 - 
     -         WProd2332*WProd2348*WProd2407 - 
     -         WProd2348*WProd2350*WProd2407 - 
     -         WProd2348*WProd2351*WProd2407 - 
     -         WProd2332*WProd2408 - WProd2350*WProd2408 - 
     -         WProd2351*WProd2408 - WProd2349*WProd2409 + 
     -         2*WProd2327*WProd2410 + 2*WProd2328*WProd2410 - 
     -         WProd2329*WProd2410 - 2*WProd2334*WProd2410 + 
     -         WProd2336*WProd2410 + WProd2337*WProd2410 + 
     -         2*WProd2338*WProd2410 + WProd2339*WProd2410 - 
     -         WProd2346*WProd2410 - WProd2347*WProd2410 + 
     -         WProd2411 - WProd2332*WProd2348*WProd2412 - 
     -         WProd2348*WProd2350*WProd2412 - 
     -         WProd2348*WProd2351*WProd2412))/WProd2391 - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp248*Prp252*
     -       (-(WProd2332*WProd2348*WProd2398) + 
     -         WProd2348*WProd2349*WProd2398 - 
     -         WProd2348*WProd2350*WProd2398 - 
     -         WProd2348*WProd2351*WProd2398 - 
     -         2*WProd2327*WProd2348*WProd2400 - 
     -         2*WProd2328*WProd2348*WProd2400 + 
     -         2*WProd2329*WProd2348*WProd2400 + 
     -         2*WProd2334*WProd2348*WProd2400 - 
     -         WProd2336*WProd2348*WProd2400 - 
     -         WProd2337*WProd2348*WProd2400 - 
     -         2*WProd2338*WProd2348*WProd2400 - 
     -         WProd2339*WProd2348*WProd2400 + 
     -         2*WProd2346*WProd2348*WProd2400 + 
     -         2*WProd2347*WProd2348*WProd2400 - 
     -         WProd2367*WProd2400 - WProd2368*WProd2401 - 
     -         WProd2332*WProd2348*WProd2407 + 
     -         WProd2348*WProd2349*WProd2407 - 
     -         WProd2348*WProd2350*WProd2407 - 
     -         WProd2348*WProd2351*WProd2407 - 
     -         WProd2332*WProd2408 + WProd2349*WProd2408 - 
     -         WProd2350*WProd2408 - WProd2351*WProd2408 + 
     -         2*WProd2327*WProd2410 + 2*WProd2328*WProd2410 - 
     -         2*WProd2329*WProd2410 - 2*WProd2334*WProd2410 + 
     -         WProd2336*WProd2410 + WProd2337*WProd2410 + 
     -         2*WProd2338*WProd2410 + WProd2339*WProd2410 - 
     -         2*WProd2346*WProd2410 - 2*WProd2347*WProd2410 + 
     -         WProd2411 - WProd2332*WProd2348*WProd2412 + 
     -         WProd2348*WProd2349*WProd2412 - 
     -         WProd2348*WProd2350*WProd2412 - 
     -         WProd2348*WProd2351*WProd2412))/WProd2391
        tmp1166 = tmp1166 + 
     -   Prp246*((-4*Cpl2*GS*Prp247*
     -         (2*WProd2332*WProd2392 - 2*WProd2332*WProd2393 - 
     -           2*WProd2394 - 4*WProd2327*WProd2396 - 
     -           4*WProd2328*WProd2396 + 4*WProd2329*WProd2396 + 
     -           4*WProd2334*WProd2396 - 4*WProd2336*WProd2396 - 
     -           4*WProd2337*WProd2396 - 4*WProd2338*WProd2396 - 
     -           3*WProd2339*WProd2396 + 4*WProd2346*WProd2396 + 
     -           4*WProd2347*WProd2396 - 2*WProd2397 - 
     -           2*WProd2344*WProd2398 + 2*WProd2345*WProd2398 - 
     -           2*WProd2332*WProd2348*WProd2398 + 
     -           2*WProd2349*WProd2406 - 2*WProd2344*WProd2412 + 
     -           2*WProd2345*WProd2412 - 
     -           2*WProd2332*WProd2348*WProd2412))/WProd2391 + 
     -      (4*Cpl2*GS*Prp248*
     -         (-2*WProd2332*WProd2392 + 2*WProd2349*WProd2392 - 
     -           2*WProd2350*WProd2392 - 2*WProd2351*WProd2392 + 
     -           2*WProd2332*WProd2393 - 2*WProd2349*WProd2393 + 
     -           2*WProd2350*WProd2393 + 2*WProd2351*WProd2393 + 
     -           2*WProd2394 + 4*WProd2327*WProd2396 + 
     -           4*WProd2328*WProd2396 - 4*WProd2329*WProd2396 - 
     -           5*WProd2334*WProd2396 + 3*WProd2336*WProd2396 + 
     -           3*WProd2337*WProd2396 + 4*WProd2338*WProd2396 + 
     -           3*WProd2339*WProd2396 - 4*WProd2346*WProd2396 - 
     -           4*WProd2347*WProd2396 + 2*WProd2397 + 
     -           2*WProd2332*WProd2348*WProd2398 - 
     -           2*WProd2348*WProd2349*WProd2398 + 
     -           2*WProd2348*WProd2350*WProd2398 + 
     -           2*WProd2348*WProd2351*WProd2398 - 
     -           2*WProd2353*WProd2400 + 2*WProd2355*WProd2401 + 
     -           2*WProd2332*WProd2348*WProd2412 - 
     -           2*WProd2348*WProd2349*WProd2412 + 
     -           2*WProd2348*WProd2350*WProd2412 + 
     -           2*WProd2348*WProd2351*WProd2412))/WProd2391)
        amp(30,4,3) = Prp249*tmp1166 + 
     -   Prp245*Prp246*(Prp247*
     -       ((-8*Cpl2*GS*WProd2390)/WProd2389 - 
     -         (8*Cpl2*GS*(WProd2332*WProd2392 - 
     -              WProd2332*WProd2393 - WProd2394 + WProd2395 - 
     -              4*WProd2327*WProd2396 - 
     -              4*WProd2328*WProd2396 + 
     -              4*WProd2329*WProd2396 + 
     -              4*WProd2334*WProd2396 - 
     -              4*WProd2336*WProd2396 - 
     -              4*WProd2337*WProd2396 - 
     -              4*WProd2338*WProd2396 - 
     -              4*WProd2339*WProd2396 + 
     -              4*WProd2346*WProd2396 + 
     -              4*WProd2347*WProd2396 - WProd2397 - 
     -              WProd2344*WProd2398 + WProd2345*WProd2398 - 
     -              WProd2332*WProd2348*WProd2398 - 
     -              WProd2332*WProd2399))/WProd2391) + 
     -      (8*Cpl2*GS*Prp248*
     -         (-(WProd2332*WProd2392) + WProd2349*WProd2392 - 
     -           WProd2350*WProd2392 - WProd2351*WProd2392 + 
     -           WProd2332*WProd2393 - WProd2349*WProd2393 + 
     -           WProd2350*WProd2393 + WProd2351*WProd2393 + 
     -           WProd2394 + 2*WProd2327*WProd2396 + 
     -           2*WProd2328*WProd2396 - 2*WProd2329*WProd2396 - 
     -           3*WProd2334*WProd2396 + 2*WProd2336*WProd2396 + 
     -           2*WProd2337*WProd2396 + 2*WProd2338*WProd2396 + 
     -           2*WProd2339*WProd2396 - 2*WProd2346*WProd2396 - 
     -           2*WProd2347*WProd2396 + WProd2397 + 
     -           WProd2332*WProd2348*WProd2398 - 
     -           WProd2348*WProd2349*WProd2398 + 
     -           WProd2348*WProd2350*WProd2398 + 
     -           WProd2348*WProd2351*WProd2398 + 
     -           WProd2332*WProd2399 - WProd2349*WProd2399 + 
     -           WProd2350*WProd2399 + WProd2351*WProd2399 - 
     -           WProd2353*WProd2400 + WProd2355*WProd2401))/
     -       WProd2391)

#ifdef DEBUG
	print *, 'amp(30,4,3) =', amp(30,4,3)
#endif

        tmp1167 = Prp246*
     -    ((-4*Cpl2*GS*Prp247*
     -         (2*WProd2373*WProd2392 - 2*WProd2373*WProd2393 - 
     -           2*WProd2372*WProd2398 - 
     -           2*WProd2348*WProd2373*WProd2398 + 
     -           2*WProd2375*WProd2398 - 2*WProd2376*WProd2398 + 
     -           2*WProd2377*WProd2406 - 2*WProd2372*WProd2412 - 
     -           2*WProd2348*WProd2373*WProd2412 + 
     -           2*WProd2375*WProd2412 - 2*WProd2376*WProd2412 - 
     -           2*WProd2418 - 4*WProd2327*WProd2419 - 
     -           4*WProd2328*WProd2419 + 4*WProd2329*WProd2419 + 
     -           4*WProd2334*WProd2419 - 4*WProd2336*WProd2419 - 
     -           4*WProd2337*WProd2419 - 4*WProd2338*WProd2419 - 
     -           3*WProd2339*WProd2419 + 4*WProd2346*WProd2419 + 
     -           4*WProd2347*WProd2419))/WProd2416 + 
     -      (4*Cpl2*GS*Prp248*
     -         (-2*WProd2373*WProd2392 + 2*WProd2377*WProd2392 - 
     -           2*WProd2378*WProd2392 - 2*WProd2379*WProd2392 + 
     -           2*WProd2373*WProd2393 - 2*WProd2377*WProd2393 + 
     -           2*WProd2378*WProd2393 + 2*WProd2379*WProd2393 + 
     -           2*WProd2348*WProd2373*WProd2398 + 
     -           2*WProd2376*WProd2398 - 
     -           2*WProd2348*WProd2377*WProd2398 + 
     -           2*WProd2348*WProd2378*WProd2398 + 
     -           2*WProd2348*WProd2379*WProd2398 + 
     -           2*WProd2382*WProd2401 + 
     -           2*WProd2348*WProd2373*WProd2412 + 
     -           2*WProd2376*WProd2412 - 
     -           2*WProd2348*WProd2377*WProd2412 + 
     -           2*WProd2348*WProd2378*WProd2412 + 
     -           2*WProd2348*WProd2379*WProd2412 + 2*WProd2418 + 
     -           4*WProd2327*WProd2419 + 4*WProd2328*WProd2419 - 
     -           4*WProd2329*WProd2419 - 5*WProd2334*WProd2419 + 
     -           3*WProd2336*WProd2419 + 3*WProd2337*WProd2419 + 
     -           4*WProd2338*WProd2419 + 3*WProd2339*WProd2419 - 
     -           4*WProd2346*WProd2419 - 4*WProd2347*WProd2419 - 
     -           2*WProd2353*WProd2420))/WProd2416) + 
     -   Prp247*Prp253*(2.6666666666666665D0*
     -       (Cpl2*GS*WProd2421)/WProd2413 + 
     -      2.6666666666666665D0*
     -       (Cpl2*GS*(WProd2372*WProd2398 + 
     -            WProd2348*WProd2373*WProd2398 - 
     -            WProd2375*WProd2398 + WProd2376*WProd2398 - 
     -            WProd2372*WProd2401 - 
     -            WProd2348*WProd2373*WProd2401 + 
     -            WProd2375*WProd2401 - WProd2376*WProd2401 - 
     -            WProd2377*WProd2404 + WProd2372*WProd2407 + 
     -            WProd2348*WProd2373*WProd2407 - 
     -            WProd2375*WProd2407 + WProd2376*WProd2407 + 
     -            WProd2373*WProd2408 + WProd2377*WProd2409 + 
     -            WProd2372*WProd2412 + 
     -            WProd2348*WProd2373*WProd2412 - 
     -            WProd2375*WProd2412 + WProd2376*WProd2412 + 
     -            WProd2422 + 2*WProd2327*WProd2423 + 
     -            2*WProd2328*WProd2423 - 2*WProd2329*WProd2423 - 
     -            2*WProd2334*WProd2423 + 2*WProd2336*WProd2423 + 
     -            2*WProd2337*WProd2423 + 2*WProd2338*WProd2423 + 
     -            2*WProd2339*WProd2423 - 2*WProd2346*WProd2423 - 
     -            2*WProd2347*WProd2423 - 2*WProd2327*WProd2424 - 
     -            2*WProd2328*WProd2424 + 2*WProd2329*WProd2424 + 
     -            2*WProd2334*WProd2424 - 2*WProd2336*WProd2424 - 
     -            2*WProd2337*WProd2424 - 2*WProd2338*WProd2424 - 
     -            2*WProd2339*WProd2424 + 2*WProd2346*WProd2424 + 
     -            2*WProd2347*WProd2424 - WProd2425))/WProd2416) + 
     -   Prp247*Prp250*(-(5.333333333333333D0*
     -         (Cpl2*GS*WProd2421)/WProd2413) - 
     -      5.333333333333333D0*
     -       (Cpl2*GS*(WProd2373*WProd2392 - WProd2373*WProd2393 - 
     -            WProd2372*WProd2401 - 
     -            WProd2348*WProd2373*WProd2401 + 
     -            WProd2375*WProd2401 - WProd2376*WProd2401 - 
     -            WProd2377*WProd2404 + WProd2377*WProd2406 + 
     -            WProd2372*WProd2407 + 
     -            WProd2348*WProd2373*WProd2407 - 
     -            WProd2375*WProd2407 + WProd2376*WProd2407 + 
     -            WProd2373*WProd2408 + WProd2377*WProd2409 - 
     -            WProd2418 - 2*WProd2327*WProd2419 - 
     -            2*WProd2328*WProd2419 + 2*WProd2329*WProd2419 + 
     -            2*WProd2334*WProd2419 - 2*WProd2336*WProd2419 - 
     -            2*WProd2337*WProd2419 - 2*WProd2338*WProd2419 - 
     -            WProd2339*WProd2419 + 2*WProd2346*WProd2419 + 
     -            2*WProd2347*WProd2419 + WProd2422 + 
     -            2*WProd2327*WProd2423 + 2*WProd2328*WProd2423 - 
     -            2*WProd2329*WProd2423 - 2*WProd2334*WProd2423 + 
     -            2*WProd2336*WProd2423 + 2*WProd2337*WProd2423 + 
     -            2*WProd2338*WProd2423 + WProd2339*WProd2423 - 
     -            2*WProd2346*WProd2423 - 2*WProd2347*WProd2423 - 
     -            2*WProd2327*WProd2424 - 2*WProd2328*WProd2424 + 
     -            2*WProd2329*WProd2424 + 2*WProd2334*WProd2424 - 
     -            2*WProd2336*WProd2424 - 2*WProd2337*WProd2424 - 
     -            2*WProd2338*WProd2424 - WProd2339*WProd2424 + 
     -            2*WProd2346*WProd2424 + 2*WProd2347*WProd2424 - 
     -            WProd2425))/WProd2416)
        tmp1167 = tmp1167 - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp248*Prp252*
     -       (-(WProd2348*WProd2373*WProd2398) - 
     -         WProd2376*WProd2398 + 
     -         WProd2348*WProd2377*WProd2398 - 
     -         WProd2348*WProd2378*WProd2398 - 
     -         WProd2348*WProd2379*WProd2398 - 
     -         WProd2387*WProd2401 - 
     -         WProd2348*WProd2373*WProd2407 - 
     -         WProd2376*WProd2407 + 
     -         WProd2348*WProd2377*WProd2407 - 
     -         WProd2348*WProd2378*WProd2407 - 
     -         WProd2348*WProd2379*WProd2407 - 
     -         WProd2373*WProd2408 + WProd2377*WProd2408 - 
     -         WProd2378*WProd2408 - WProd2379*WProd2408 - 
     -         WProd2348*WProd2373*WProd2412 - 
     -         WProd2376*WProd2412 + 
     -         WProd2348*WProd2377*WProd2412 - 
     -         WProd2348*WProd2378*WProd2412 - 
     -         WProd2348*WProd2379*WProd2412 - 
     -         2*WProd2327*WProd2348*WProd2420 - 
     -         2*WProd2328*WProd2348*WProd2420 + 
     -         2*WProd2329*WProd2348*WProd2420 + 
     -         2*WProd2334*WProd2348*WProd2420 - 
     -         WProd2336*WProd2348*WProd2420 - 
     -         WProd2337*WProd2348*WProd2420 - 
     -         2*WProd2338*WProd2348*WProd2420 - 
     -         WProd2339*WProd2348*WProd2420 + 
     -         2*WProd2346*WProd2348*WProd2420 + 
     -         2*WProd2347*WProd2348*WProd2420 - 
     -         WProd2367*WProd2420 + 2*WProd2327*WProd2424 + 
     -         2*WProd2328*WProd2424 - 2*WProd2329*WProd2424 - 
     -         2*WProd2334*WProd2424 + WProd2336*WProd2424 + 
     -         WProd2337*WProd2424 + 2*WProd2338*WProd2424 + 
     -         WProd2339*WProd2424 - 2*WProd2346*WProd2424 - 
     -         2*WProd2347*WProd2424 + WProd2425))/WProd2416 - 
     -   2.6666666666666665D0*
     -    (Cpl2*GS*Prp252*Prp253*
     -       (-(WProd2348*WProd2373*WProd2398) - 
     -         WProd2376*WProd2398 - 
     -         WProd2348*WProd2378*WProd2398 - 
     -         WProd2348*WProd2379*WProd2398 + 
     -         WProd2348*WProd2373*WProd2401 + 
     -         WProd2376*WProd2401 + 
     -         WProd2348*WProd2378*WProd2401 + 
     -         WProd2348*WProd2379*WProd2401 - 
     -         WProd2348*WProd2373*WProd2407 - 
     -         WProd2376*WProd2407 - 
     -         WProd2348*WProd2378*WProd2407 - 
     -         WProd2348*WProd2379*WProd2407 - 
     -         WProd2373*WProd2408 - WProd2378*WProd2408 - 
     -         WProd2379*WProd2408 - WProd2377*WProd2409 - 
     -         WProd2348*WProd2373*WProd2412 - 
     -         WProd2376*WProd2412 - 
     -         WProd2348*WProd2378*WProd2412 - 
     -         WProd2348*WProd2379*WProd2412 + 
     -         2*WProd2327*WProd2424 + 2*WProd2328*WProd2424 - 
     -         WProd2329*WProd2424 - 2*WProd2334*WProd2424 + 
     -         WProd2336*WProd2424 + WProd2337*WProd2424 + 
     -         2*WProd2338*WProd2424 + WProd2339*WProd2424 - 
     -         WProd2346*WProd2424 - WProd2347*WProd2424 + 
     -         WProd2425))/WProd2416
        tmp1167 = tmp1167 + 
     -   Prp251*(5.333333333333333D0*
     -       (Cpl2*GS*Prp248*
     -          (-(WProd2373*WProd2392) + WProd2377*WProd2392 - 
     -            WProd2378*WProd2392 - WProd2379*WProd2392 + 
     -            WProd2373*WProd2393 - WProd2377*WProd2393 + 
     -            WProd2378*WProd2393 + WProd2379*WProd2393 + 
     -            WProd2382*WProd2401 - WProd2387*WProd2401 - 
     -            WProd2348*WProd2373*WProd2407 - 
     -            WProd2376*WProd2407 + 
     -            WProd2348*WProd2377*WProd2407 - 
     -            WProd2348*WProd2378*WProd2407 - 
     -            WProd2348*WProd2379*WProd2407 - 
     -            WProd2373*WProd2408 + WProd2377*WProd2408 - 
     -            WProd2378*WProd2408 - WProd2379*WProd2408 + 
     -            WProd2418 + 2*WProd2327*WProd2419 + 
     -            2*WProd2328*WProd2419 - 2*WProd2329*WProd2419 - 
     -            3*WProd2334*WProd2419 + 2*WProd2336*WProd2419 + 
     -            2*WProd2337*WProd2419 + 2*WProd2338*WProd2419 + 
     -            2*WProd2339*WProd2419 - 2*WProd2346*WProd2419 - 
     -            2*WProd2347*WProd2419 - 
     -            2*WProd2327*WProd2348*WProd2420 - 
     -            2*WProd2328*WProd2348*WProd2420 + 
     -            2*WProd2329*WProd2348*WProd2420 + 
     -            3*WProd2334*WProd2348*WProd2420 - 
     -            2*WProd2336*WProd2348*WProd2420 - 
     -            2*WProd2337*WProd2348*WProd2420 - 
     -            2*WProd2338*WProd2348*WProd2420 - 
     -            2*WProd2339*WProd2348*WProd2420 + 
     -            2*WProd2346*WProd2348*WProd2420 + 
     -            2*WProd2347*WProd2348*WProd2420 - 
     -            WProd2353*WProd2420 - WProd2367*WProd2420 + 
     -            2*WProd2327*WProd2424 + 2*WProd2328*WProd2424 - 
     -            2*WProd2329*WProd2424 - 3*WProd2334*WProd2424 + 
     -            2*WProd2336*WProd2424 + 2*WProd2337*WProd2424 + 
     -            2*WProd2338*WProd2424 + 2*WProd2339*WProd2424 - 
     -            2*WProd2346*WProd2424 - 2*WProd2347*WProd2424 + 
     -            WProd2425))/WProd2416 + 
     -      Prp250*(-(5.333333333333333D0*
     -            (Cpl2*GS*WProd2421)/WProd2413) + 
     -         5.333333333333333D0*
     -          (Cpl2*GS*(-(WProd2373*WProd2392) + 
     -               WProd2377*WProd2392 + WProd2373*WProd2393 - 
     -               WProd2377*WProd2393 + WProd2382*WProd2401 - 
     -               WProd2387*WProd2401 - WProd2388*WProd2401 - 
     -               WProd2372*WProd2407 - 
     -               WProd2348*WProd2373*WProd2407 + 
     -               WProd2375*WProd2407 - WProd2376*WProd2407 + 
     -               WProd2348*WProd2377*WProd2407 - 
     -               WProd2373*WProd2408 + WProd2377*WProd2408 + 
     -               WProd2418 + WProd2327*WProd2419 + 
     -               WProd2328*WProd2419 - WProd2329*WProd2419 - 
     -               2*WProd2334*WProd2419 + WProd2336*WProd2419 + 
     -               WProd2337*WProd2419 + 2*WProd2338*WProd2419 + 
     -               WProd2339*WProd2419 - WProd2346*WProd2419 - 
     -               WProd2347*WProd2419 - 
     -               WProd2327*WProd2348*WProd2420 - 
     -               WProd2328*WProd2348*WProd2420 + 
     -               WProd2329*WProd2348*WProd2420 + 
     -               2*WProd2334*WProd2348*WProd2420 - 
     -               WProd2336*WProd2348*WProd2420 - 
     -               WProd2337*WProd2348*WProd2420 - 
     -               2*WProd2338*WProd2348*WProd2420 - 
     -               WProd2339*WProd2348*WProd2420 + 
     -               WProd2346*WProd2348*WProd2420 + 
     -               WProd2347*WProd2348*WProd2420 - 
     -               WProd2353*WProd2420 - WProd2367*WProd2420 + 
     -               WProd2369*WProd2420 - WProd2422 - 
     -               WProd2327*WProd2423 - WProd2328*WProd2423 + 
     -               WProd2329*WProd2423 + 2*WProd2334*WProd2423 - 
     -               WProd2336*WProd2423 - WProd2337*WProd2423 - 
     -               2*WProd2338*WProd2423 - WProd2339*WProd2423 + 
     -               WProd2346*WProd2423 + WProd2347*WProd2423 + 
     -               WProd2327*WProd2424 + WProd2328*WProd2424 - 
     -               WProd2329*WProd2424 - 2*WProd2334*WProd2424 + 
     -               WProd2336*WProd2424 + WProd2337*WProd2424 + 
     -               2*WProd2338*WProd2424 + WProd2339*WProd2424 - 
     -               WProd2346*WProd2424 - WProd2347*WProd2424 + 
     -               WProd2425))/WProd2416))
        amp(30,4,10) = Prp249*tmp1167 + 
     -   Prp245*Prp246*(Prp247*
     -       ((-8*Cpl2*GS*(WProd2414 - WProd2415))/WProd2413 - 
     -         (8*Cpl2*GS*(WProd2373*WProd2392 - 
     -              WProd2373*WProd2393 - WProd2372*WProd2398 - 
     -              WProd2348*WProd2373*WProd2398 + 
     -              WProd2375*WProd2398 - WProd2376*WProd2398 - 
     -              WProd2373*WProd2399 + WProd2417 - WProd2418 - 
     -              4*WProd2327*WProd2419 - 
     -              4*WProd2328*WProd2419 + 
     -              4*WProd2329*WProd2419 + 
     -              4*WProd2334*WProd2419 - 
     -              4*WProd2336*WProd2419 - 
     -              4*WProd2337*WProd2419 - 
     -              4*WProd2338*WProd2419 - 
     -              4*WProd2339*WProd2419 + 
     -              4*WProd2346*WProd2419 + 4*WProd2347*WProd2419))
     -           /WProd2416) + 
     -      Prp248*((8*Cpl2*GS*WProd2415)/WProd2413 + 
     -         (8*Cpl2*GS*(-(WProd2373*WProd2392) + 
     -              WProd2377*WProd2392 - WProd2378*WProd2392 - 
     -              WProd2379*WProd2392 + WProd2373*WProd2393 - 
     -              WProd2377*WProd2393 + WProd2378*WProd2393 + 
     -              WProd2379*WProd2393 + 
     -              WProd2348*WProd2373*WProd2398 + 
     -              WProd2376*WProd2398 - 
     -              WProd2348*WProd2377*WProd2398 + 
     -              WProd2348*WProd2378*WProd2398 + 
     -              WProd2348*WProd2379*WProd2398 + 
     -              WProd2373*WProd2399 - WProd2377*WProd2399 + 
     -              WProd2378*WProd2399 + WProd2379*WProd2399 + 
     -              WProd2382*WProd2401 + WProd2418 + 
     -              2*WProd2327*WProd2419 + 
     -              2*WProd2328*WProd2419 - 
     -              2*WProd2329*WProd2419 - 
     -              3*WProd2334*WProd2419 + 
     -              2*WProd2336*WProd2419 + 
     -              2*WProd2337*WProd2419 + 
     -              2*WProd2338*WProd2419 + 
     -              2*WProd2339*WProd2419 - 
     -              2*WProd2346*WProd2419 - 
     -              2*WProd2347*WProd2419 - WProd2353*WProd2420))/
     -          WProd2416))

#ifdef DEBUG
	print *, 'amp(30,4,10) =', amp(30,4,10)
#endif

      endif

      end


