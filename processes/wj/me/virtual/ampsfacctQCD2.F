      subroutine ampsfacctQCD2(amp,resnr)
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
*     procnr=1:   { up  , ~do  ,  ve  , ~el  ,  gl  }
*     procnr=2:   { up  ,  ga  ,  ve  , ~el  ,  do  }
*     procnr=3:   { up  ,  gl  ,  ve  , ~el  ,  do  }
*     procnr=4:   {~do  ,  ga  ,  ve  , ~el  , ~up  }
*     procnr=5:   {~do  ,  gl  ,  ve  , ~el  , ~up  }
*
*   Helicities:
*     helnr=1:   {-1, 1, -1, 1, -1}
*     helnr=2:   {-1, 1, -1, 1, 1}
*     helnr=3:   {-1, -1, -1, 1, -1}
*     helnr=4:   {1, -1, -1, 1, 1}
*     helnr=5:   {1, 1, -1, 1, 1}
*
*   Color Structures:
*     colnr=1:   SUNT[Glu1, Col2, Col1]
*     colnr=2:   SUNT[Col1, Col2]
*     colnr=3:   SUNT[Glu1, Col1, Col2]


      if(resnr.eq.0)then

	tmp9 = Alfas*dZfL1dimreg(3,1,1) + Alfas*dZfL1dimreg(4,1,1)
        amp(2,2,3) = -(2.8284271247461903D0*
     -      (Cpl2*Prp62*Prp63*tmp9*
     -         (WProd525 + WProd523*WProd532 - WProd533 - WProd534)
     -         )/WProd521) + 
     -   3.7712361663282534D0*
     -    (Alfas*Cpl2*Prp60*Prp61**2*
     -       (WProd522*WProd523 - WProd524 + WProd525)*
     -       (WProd526 + WProd527 + WProd528 + WProd529 - 
     -         WProd530)*dZfL1dimreg(3,1,1))/WProd521 - 
     -   1.8856180831641267D0*
     -    (Cpl2*Prp60*Prp61*
     -       (WProd522*WProd523 - WProd524 + WProd525)*
     -       (3*Alfas*dZfL1dimreg(3,1,1) + 
     -         Alfas*dZfL1dimreg(4,1,1)))/WProd521 + 
     -   Prp60*(-(2.8284271247461903D0*
     -         (Cpl2*Prp63*tmp9*
     -            (-WProd524 + WProd525 + WProd523*WProd532 + 
     -              WProd523*WProd536))/WProd521) + 
     -      1.8856180831641267D0*
     -       (Alfas*Cpl2*Prp64**2*WProd523*
     -          (-WProd527 - WProd529 + WProd531)*
     -          (WProd522 - WProd532 - WProd536)*dZfL1dimreg(4,1,1)
     -          )/WProd521 + 
     -      0.9428090415820634D0*
     -       (Cpl2*Prp64*WProd523*(WProd522 - WProd532 - WProd536)*
     -          (Alfas*dZfL1dimreg(3,1,1) + 
     -            3*Alfas*dZfL1dimreg(4,1,1)))/WProd521)

#ifdef DEBUG
	print *, 'amp(2,2,3) =', amp(2,2,3)
#endif

	tmp10 = WProd523*WProd539 - WProd540 + WProd541 - WProd542
        tmp11 = WProd523*WProd539 - WProd540 - 
     -   WProd523*WProd544 - WProd523*WProd547
	tmp12 = Alfas*dZfL1dimreg(3,1,1) + Alfas*dZfL1dimreg(4,1,1)
        amp(2,2,1) = -(2.8284271247461903D0*
     -      (Cpl2*Prp62*Prp63*tmp12*
     -         (WProd541 + WProd543 + WProd523*WProd544 - 
     -           WProd545 - WProd546))/WProd538) + 
     -   3.7712361663282534D0*
     -    (Alfas*Cpl2*Prp60*Prp61**2*tmp10*
     -       (WProd526 + WProd527 + WProd528 + WProd529 - 
     -         WProd530)*dZfL1dimreg(3,1,1))/WProd538 - 
     -   1.8856180831641267D0*
     -    (Cpl2*Prp60*Prp61*tmp10*
     -       (3*Alfas*dZfL1dimreg(3,1,1) + 
     -         Alfas*dZfL1dimreg(4,1,1)))/WProd538 + 
     -   Prp60*(-(2.8284271247461903D0*
     -         (Cpl2*Prp63*tmp12*
     -            (WProd541 - WProd542 + WProd523*WProd544 + 
     -              WProd523*WProd547))/WProd538) + 
     -      1.8856180831641267D0*
     -       (Alfas*Cpl2*Prp64**2*tmp11*
     -          (-WProd527 - WProd529 + WProd531)*
     -          dZfL1dimreg(4,1,1))/WProd538 + 
     -      0.9428090415820634D0*
     -       (Cpl2*Prp64*tmp11*
     -          (Alfas*dZfL1dimreg(3,1,1) + 
     -            3*Alfas*dZfL1dimreg(4,1,1)))/WProd538)

#ifdef DEBUG
	print *, 'amp(2,2,1) =', amp(2,2,1)
#endif

      endif

      end


