      subroutine abbrrealQED30(resnr)
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

      
	WProd2389 = Prod(2,2,2,1)
	WProd2330 = Prod(1,2,1,1)*Prod(2,2,2,1)
	WProd2391 = Prod(2,2,2,1)*ProdC(1,2,1,1)
	WProd2334 = Prod(2,1,1,1)*ProdC(2,1,1,1)
	WProd2354 = Prod(2,1,1,2)*ProdC(2,1,1,1)
	WProd2349 = Prod(2,2,1,1)*ProdC(2,1,1,1)
	WProd2352 = Prod(2,2,1,2)*ProdC(2,1,1,1)
	WProd2401 = Prod(2,1,1,1)*ProdC(2,1,1,2)
	WProd2400 = Prod(2,2,1,1)*ProdC(2,1,1,2)
	WProd2377 = Prod(2,1,1,1)*ProdC(2,2,1,1)
	WProd2380 = Prod(2,1,1,2)*ProdC(2,2,1,1)
	WProd2420 = Prod(2,1,1,1)*ProdC(2,2,1,2)
	WProd2413 = ProdC(2,2,2,1)
	WProd2371 = Prod(1,2,1,1)*ProdC(2,2,2,1)
	WProd2416 = ProdC(1,2,1,1)*ProdC(2,2,2,1)
	WProd2328 = Prod(3,1,1,1)*ProdC(3,1,1,1)
	WProd2365 = Prod(3,1,1,2)*ProdC(3,1,1,1)
	WProd2412 = Prod(3,1,1,1)*ProdC(3,1,1,2)
	WProd2336 = Prod(3,1,2,1)*ProdC(3,1,2,1)
	WProd2350 = Prod(3,1,2,2)*ProdC(3,1,2,1)
	WProd2378 = Prod(3,1,2,1)*ProdC(3,1,2,2)
	WProd2327 = Prod(4,1,1,1)*ProdC(4,1,1,1)
	WProd2343 = Prod(4,1,1,2)*ProdC(4,1,1,1)
	WProd2398 = Prod(4,1,1,1)*ProdC(4,1,1,2)
	WProd2337 = Prod(4,1,2,1)*ProdC(4,1,2,1)
	WProd2351 = Prod(4,1,2,2)*ProdC(4,1,2,1)
	WProd2379 = Prod(4,1,2,1)*ProdC(4,1,2,2)
	WProd2329 = Prod(4,1,3,1)*ProdC(4,1,3,1)
	WProd2338 = Prod(5,1,1,1)*ProdC(5,1,1,1)
	WProd2361 = Prod(5,1,1,2)*ProdC(5,1,1,1)
        WProd2331 = Prod(3,1,1,2)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd2333 = Prod(3,1,1,1)*Prod(6,1,1,2)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,1)
        WProd2392 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,2)*
     -   ProdC(5,1,1,1)
        WProd2369 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd2370 = Prod(3,1,2,2)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd2357 = Prod(3,1,2,1)*Prod(6,1,1,2)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd2358 = Prod(3,1,2,2)*Prod(6,1,1,2)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,1)
        WProd2359 = Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,1,2)*
     -   ProdC(2,1,1,1)*ProdC(4,1,2,1)*ProdC(5,1,1,1)
        WProd2403 = Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -   ProdC(2,1,1,2)*ProdC(4,1,2,1)*ProdC(5,1,1,1)
        WProd2383 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,2)*
     -   ProdC(2,2,1,1)*ProdC(4,1,2,1)*ProdC(5,1,1,1)
        WProd2422 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -   ProdC(2,2,1,2)*ProdC(4,1,2,1)*ProdC(5,1,1,1)
        WProd2388 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,2)*
     -   ProdC(5,1,1,1)
        WProd2384 = Prod(3,1,2,1)*Prod(6,1,1,2)*ProdC(4,1,2,2)*
     -   ProdC(5,1,1,1)
	WProd2407 = Prod(5,1,1,1)*ProdC(5,1,1,2)
        WProd2393 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,1,2)
        WProd2404 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,2)
        WProd2405 = Prod(3,1,2,2)*Prod(6,1,1,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,1,2)
        WProd2423 = Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(4,1,2,2)*
     -   ProdC(5,1,1,2)
	WProd2339 = Prod(5,1,2,1)*ProdC(5,1,2,1)
	WProd2332 = Prod(5,1,2,2)*ProdC(5,1,2,1)
        WProd2353 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd2360 = Prod(3,1,1,2)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd2355 = Prod(3,1,1,1)*Prod(6,1,2,2)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd2335 = Prod(3,1,1,2)*Prod(6,1,2,2)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,1)
        WProd2340 = -(Prod(2,1,2,2)*Prod(3,1,1,1)*Prod(6,1,1,2)*
     -     ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(5,1,2,1))
        WProd2341 = Prod(2,2,1,1)*Prod(3,1,1,2)*Prod(6,1,2,1)*
     -   ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(5,1,2,1)
        WProd2342 = Prod(2,1,1,1)*Prod(3,1,1,2)*Prod(6,1,2,2)*
     -   ProdC(2,1,1,1)*ProdC(4,1,1,1)*ProdC(5,1,2,1)
        WProd2394 = -(Prod(2,1,2,2)*Prod(3,1,1,1)*Prod(6,1,1,1)*
     -     ProdC(2,1,1,2)*ProdC(4,1,1,1)*ProdC(5,1,2,1))
        WProd2395 = Prod(2,2,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -   ProdC(2,1,1,2)*ProdC(4,1,1,1)*ProdC(5,1,2,1)
        WProd2381 = Prod(2,1,1,1)*Prod(3,1,1,2)*Prod(6,1,2,1)*
     -   ProdC(2,2,1,1)*ProdC(4,1,1,1)*ProdC(5,1,2,1)
        WProd2417 = Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -   ProdC(2,2,1,2)*ProdC(4,1,1,1)*ProdC(5,1,2,1)
        WProd2406 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,2)*
     -   ProdC(5,1,2,1)
        WProd2396 = Prod(3,1,1,1)*Prod(6,1,2,2)*ProdC(4,1,1,2)*
     -   ProdC(5,1,2,1)
        WProd2397 = Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,2)*
     -   ProdC(2,1,1,1)*ProdC(4,1,1,2)*ProdC(5,1,2,1)
        WProd2418 = Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -   ProdC(2,2,1,1)*ProdC(4,1,1,2)*ProdC(5,1,2,1)
        WProd2402 = Prod(2,1,1,1)*Prod(3,1,2,2)*Prod(6,1,1,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,2,1)
        WProd2344 = Prod(3,1,2,2)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd2390 = Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,2)*
     -   ProdC(4,1,2,1)*ProdC(5,1,2,1)
        WProd2345 = Prod(3,1,2,1)*Prod(6,1,2,2)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,1)
        WProd2421 = Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -   ProdC(4,1,2,2)*ProdC(5,1,2,1)
        WProd2372 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,2)*
     -   ProdC(5,1,2,1)
	WProd2373 = Prod(5,1,2,1)*ProdC(5,1,2,2)
        WProd2382 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,2)
        WProd2374 = Prod(3,1,1,2)*Prod(6,1,2,1)*ProdC(4,1,1,1)*
     -   ProdC(5,1,2,2)
        WProd2419 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(4,1,1,2)*
     -   ProdC(5,1,2,2)
        WProd2414 = Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -   ProdC(4,1,2,1)*ProdC(5,1,2,2)
        WProd2375 = Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(4,1,2,1)*
     -   ProdC(5,1,2,2)
	WProd2346 = Prod(5,1,3,1)*ProdC(5,1,3,1)
	WProd2347 = Prod(5,1,4,1)*ProdC(5,1,4,1)
	WProd2348 = Prod(6,1,3,1)*ProdC(5,1,4,1)
        WProd2408 = -(Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(1,1,1,2)*
     -     ProdC(5,1,4,1))
        WProd2399 = Prod(3,1,1,1)*Prod(6,1,1,1)*ProdC(1,1,1,2)*
     -   ProdC(5,1,4,1)
        WProd2367 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd2362 = -(Prod(3,1,2,1)*Prod(6,1,1,2)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd2363 = -(Prod(3,1,2,2)*Prod(6,1,1,2)*ProdC(2,1,1,1)*
     -     ProdC(5,1,4,1))
        WProd2368 = Prod(3,1,1,1)*Prod(6,1,2,2)*ProdC(2,1,1,1)*
     -   ProdC(5,1,4,1)
        WProd2364 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,1,2)*
     -     ProdC(2,1,1,1)**2*ProdC(5,1,4,1))
        WProd2409 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(2,1,1,2)*
     -     ProdC(5,1,4,1))
        WProd2410 = -(Prod(3,1,2,2)*Prod(6,1,1,1)*ProdC(2,1,1,2)*
     -     ProdC(5,1,4,1))
        WProd2411 = -(Prod(2,2,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(2,1,1,1)*ProdC(2,1,1,2)*ProdC(5,1,4,1))
        WProd2415 = -(Prod(2,1,1,1)*Prod(3,1,1,1)*Prod(6,1,2,1)*
     -     ProdC(2,1,2,2)*ProdC(5,1,4,1))
        WProd2376 = -(Prod(3,1,2,1)*Prod(6,1,2,1)*ProdC(2,1,2,2)*
     -     ProdC(5,1,4,1))
        WProd2385 = -(Prod(3,1,2,1)*Prod(6,1,1,2)*ProdC(2,2,1,1)*
     -     ProdC(5,1,4,1))
        WProd2387 = Prod(3,1,1,1)*Prod(6,1,2,1)*ProdC(2,2,1,1)*
     -   ProdC(5,1,4,1)
        WProd2386 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,2)*
     -     ProdC(2,1,1,1)*ProdC(2,2,1,1)*ProdC(5,1,4,1))
        WProd2424 = -(Prod(3,1,2,1)*Prod(6,1,1,1)*ProdC(2,2,1,2)*
     -     ProdC(5,1,4,1))
        WProd2425 = -(Prod(2,1,1,1)*Prod(3,1,2,1)*Prod(6,1,1,1)*
     -     ProdC(2,1,1,1)*ProdC(2,2,1,2)*ProdC(5,1,4,1))
	WProd2366 = Prod(6,1,1,1)*ProdC(6,1,1,1)
	WProd2356 = Prod(6,1,2,1)*ProdC(6,1,2,1)
      
      
	Prp245 = -(1.D0/WProd2327)
	Prp249 = 1.D0/((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r + WProd2329)
        Prp246 = 1.D0/
     -   ((0.D0,1.D0)*GW*Sqrt(MW2r) - MW2r - WProd2327 - WProd2328 + 
     -     WProd2329)
	Prp251 = -(1.D0/WProd2338)
	Prp250 = 1.D0/(WProd2334 - WProd2338 - WProd2339)
	Prp247 = -(1.D0/WProd2339)
	Prp248 = -(1.D0/WProd2356)
	Prp252 = 1.D0/(WProd2334 - WProd2356 - WProd2366)
	Prp253 = -(1.D0/WProd2366)
      
      
      
      endif

      end


