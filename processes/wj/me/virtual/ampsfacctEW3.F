      subroutine ampsfacctEW3(amp,resnr)
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

        amp(3,1,3) = (Prp51**2*Prp52*
     -      (-(5.656854249492381D0*
     -           (Cpl105*dZW1*GS*WProd459*
     -             (WProd460 - WProd461 - WProd462))) - 
     -        5.656854249492381D0*
     -         (Cpl1*GS*(dMWsq1 - dZW1*WProd457)*WProd459*
     -           (WProd460 - WProd461 - WProd462))))/WProd456 + 
     -   5.656854249492381D0*
     -    (Cpl1*GS*Prp51*Prp53**2*
     -       (WProd459*WProd460 - WProd463 + WProd464)*
     -       (-WProd457 + WProd465 + WProd466 + WProd467 + 
     -         WProd468)*dZfL1(3,1,1))/WProd456 + 
     -   Prp53*((Prp51**2*(-(5.656854249492381D0*
     -              (Cpl105*dZW1*GS*
     -                (WProd459*WProd460 - WProd463 + WProd464)))-
     -             5.656854249492381D0*
     -            (Cpl1*GS*(dMWsq1 - dZW1*WProd457)*
     -              (WProd459*WProd460 - WProd463 + WProd464))))/
     -       WProd456 + (Prp51*
     -         (11.313708498984761D0*
     -            (Cpl106*dSW1*GS*
     -              (WProd459*WProd460 - WProd463 + WProd464)) - 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*
     -              (WProd459*WProd460 - WProd463 + WProd464)*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + 3*dZfL1(3,1,1) + dZfL1(4,1,1))
     -              )))/WProd456) + 
     -   Prp51*(-(5.656854249492381D0*
     -         (Cpl1*GS*Prp52**2*WProd459*
     -            (WProd460 - WProd461 - WProd462)*
     -            (-WProd466 - WProd468 + WProd469)*dZfL1(4,1,1))/
     -          WProd456) + 
     -      (Prp52*(11.313708498984761D0*
     -            (Cpl106*dSW1*GS*WProd459*
     -              (WProd460 - WProd461 - WProd462)) - 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*WProd459*
     -              (WProd460 - WProd461 - WProd462)*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + dZfL1(3,1,1) + 3*dZfL1(4,1,1))
     -              )))/WProd456)

#ifdef DEBUG
	print *, 'amp(3,1,3) =', amp(3,1,3)
#endif

        tmp3 = WProd459*WProd471 - WProd459*WProd472 - 
     -   WProd473 - WProd459*WProd474
	tmp4 = WProd459*WProd471 - WProd473 + WProd475 - WProd476
        amp(3,1,1) = (Prp51**2*Prp52*
     -      (-(5.656854249492381D0*(Cpl105*dZW1*GS*tmp3)) - 
     -        5.656854249492381D0*
     -         (Cpl1*GS*tmp3*(dMWsq1 - dZW1*WProd457))))/WProd470+
     -     5.656854249492381D0*
     -    (Cpl1*GS*Prp51*Prp53**2*tmp4*
     -       (-WProd457 + WProd465 + WProd466 + WProd467 + 
     -         WProd468)*dZfL1(3,1,1))/WProd470 + 
     -   Prp53*((Prp51**2*(-(5.656854249492381D0*
     -              (Cpl105*dZW1*GS*tmp4)) - 
     -           5.656854249492381D0*
     -            (Cpl1*GS*tmp4*(dMWsq1 - dZW1*WProd457))))/
     -       WProd470 + (Prp51*
     -         (11.313708498984761D0*(Cpl106*dSW1*GS*tmp4) - 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*tmp4*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + 3*dZfL1(3,1,1) + dZfL1(4,1,1))
     -              )))/WProd470) + 
     -   Prp51*(-(5.656854249492381D0*
     -         (Cpl1*GS*Prp52**2*tmp3*
     -            (-WProd466 - WProd468 + WProd469)*dZfL1(4,1,1))/
     -          WProd470) + 
     -      (Prp52*(11.313708498984761D0*(Cpl106*dSW1*GS*tmp3) - 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*tmp3*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + dZfL1(3,1,1) + 3*dZfL1(4,1,1))
     -              )))/WProd470)

#ifdef DEBUG
	print *, 'amp(3,1,1) =', amp(3,1,1)
#endif

      endif

      end


