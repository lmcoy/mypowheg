      subroutine ampsfacctEW5(amp,resnr)
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

        tmp5 = WProd480 - WProd481 - WProd482*WProd483 + 
     -   WProd483*WProd484 + WProd483*WProd485
        amp(5,3,4) = (Prp54**2*Prp55*
     -      (5.656854249492381D0*(Cpl105*dZW1*GS*tmp5) + 
     -        5.656854249492381D0*
     -         (Cpl1*GS*tmp5*(dMWsq1 - dZW1*WProd478))))/WProd477+
     -     5.656854249492381D0*
     -    (Cpl1*GS*Prp54*Prp56**2*WProd482*WProd483*
     -       (-WProd478 + WProd486 + WProd487 + WProd488 + 
     -         WProd489)*dZfL1(4,1,1))/WProd477 + 
     -   Prp54*(5.656854249492381D0*
     -       (Cpl1*GS*Prp55**2*tmp5*
     -          (-WProd487 - WProd489 + WProd490)*dZfL1(3,1,1))/
     -        WProd477 + (Prp55*
     -         (-(11.313708498984761D0*(Cpl106*dSW1*GS*tmp5)) + 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*tmp5*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + 3*dZfL1(3,1,1) + dZfL1(4,1,1))
     -              )))/WProd477) + 
     -   Prp56*((Prp54**2*(-(5.656854249492381D0*
     -              (Cpl105*dZW1*GS*WProd482*WProd483)) - 
     -           5.656854249492381D0*
     -            (Cpl1*GS*(dMWsq1 - dZW1*WProd478)*WProd482*
     -              WProd483)))/WProd477 + 
     -      (Prp54*(11.313708498984761D0*
     -            (Cpl106*dSW1*GS*WProd482*WProd483) - 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*WProd482*WProd483*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + dZfL1(3,1,1) + 3*dZfL1(4,1,1))
     -              )))/WProd477)

#ifdef DEBUG
	print *, 'amp(5,3,4) =', amp(5,3,4)
#endif

        tmp6 = WProd492 - WProd493 + WProd483*WProd494 - 
     -   WProd483*WProd495 - WProd483*WProd496 + WProd497
        amp(5,3,5) = (Prp54**2*Prp55*
     -      (-(5.656854249492381D0*(Cpl105*dZW1*GS*tmp6)) - 
     -        5.656854249492381D0*
     -         (Cpl1*GS*tmp6*(dMWsq1 - dZW1*WProd478))))/WProd491+
     -     5.656854249492381D0*
     -    (Cpl1*GS*Prp54*Prp56**2*
     -       (-WProd478 + WProd486 + WProd487 + WProd488 + 
     -         WProd489)*(WProd483*WProd494 + WProd497)*
     -       dZfL1(4,1,1))/WProd491 + 
     -   Prp54*(-(5.656854249492381D0*
     -         (Cpl1*GS*Prp55**2*tmp6*
     -            (-WProd487 - WProd489 + WProd490)*dZfL1(3,1,1))/
     -          WProd491) + 
     -      (Prp55*(11.313708498984761D0*(Cpl106*dSW1*GS*tmp6) - 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*tmp6*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + 3*dZfL1(3,1,1) + dZfL1(4,1,1))
     -              )))/WProd491) + 
     -   Prp56*((Prp54**2*(-(5.656854249492381D0*
     -              (Cpl105*dZW1*GS*(WProd483*WProd494 + WProd497))
     -              ) - 5.656854249492381D0*
     -            (Cpl1*GS*(dMWsq1 - dZW1*WProd478)*
     -              (WProd483*WProd494 + WProd497))))/WProd491 + 
     -      (Prp54*(11.313708498984761D0*
     -            (Cpl106*dSW1*GS*(WProd483*WProd494 + WProd497))-
     -             2.8284271247461903D0*
     -            (Cpl1*GS*(WProd483*WProd494 + WProd497)*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + dZfL1(3,1,1) + 3*dZfL1(4,1,1))
     -              )))/WProd491)

#ifdef DEBUG
	print *, 'amp(5,3,5) =', amp(5,3,5)
#endif

      endif

      end


