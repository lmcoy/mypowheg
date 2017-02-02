      subroutine ampsfacctEW1(amp,resnr)
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

        tmp1 = -WProd437 + WProd439*WProd440 - WProd441 + 
     -   WProd439*WProd442 - WProd443 - WProd444
        amp(1,1,1) = Prp48**2*
     -    ((Prp50*(-(5.656854249492381D0*(Cpl105*dZW1*GS*tmp1)) - 
     -           5.656854249492381D0*
     -            (Cpl1*GS*tmp1*(dMWsq1 - dZW1*WProd433))))/
     -       WProd434 + (Prp49*
     -         (-(5.656854249492381D0*
     -              (Cpl105*dZW1*GS*(WProd436 - WProd437))) - 
     -           5.656854249492381D0*
     -            (Cpl1*GS*(dMWsq1 - dZW1*WProd433)*
     -              (WProd436 - WProd437))))/WProd434) + 
     -   Prp48*(-(5.656854249492381D0*
     -         (Cpl1*GS*Prp50**2*tmp1*
     -            (-WProd433 + WProd445 + WProd446)*dZfL1(3,1,1))/
     -          WProd434) - 
     -      5.656854249492381D0*
     -       (Cpl1*GS*Prp49**2*(WProd436 - WProd437)*
     -          (-WProd445 - WProd446 + WProd447)*dZfL1(4,1,1))/
     -        WProd434 + (Prp50*
     -         (11.313708498984761D0*(Cpl106*dSW1*GS*tmp1) - 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*tmp1*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + 3*dZfL1(3,1,1) + dZfL1(4,1,1))
     -              )))/WProd434 + 
     -      (Prp49*(11.313708498984761D0*
     -            (Cpl106*dSW1*GS*(WProd436 - WProd437)) - 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*(WProd436 - WProd437)*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + dZfL1(3,1,1) + 3*dZfL1(4,1,1))
     -              )))/WProd434)

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        tmp2 = -WProd450 + WProd439*WProd451 - WProd452 + 
     -   WProd439*WProd453 - WProd454 - WProd455
        amp(1,1,2) = Prp48**2*
     -    ((Prp50*(-(5.656854249492381D0*(Cpl105*dZW1*GS*tmp2)) - 
     -           5.656854249492381D0*
     -            (Cpl1*GS*tmp2*(dMWsq1 - dZW1*WProd433))))/
     -       WProd448 + (Prp49*
     -         (-(5.656854249492381D0*
     -              (Cpl105*dZW1*GS*(WProd449 - WProd450))) - 
     -           5.656854249492381D0*
     -            (Cpl1*GS*(dMWsq1 - dZW1*WProd433)*
     -              (WProd449 - WProd450))))/WProd448) + 
     -   Prp48*(-(5.656854249492381D0*
     -         (Cpl1*GS*Prp50**2*tmp2*
     -            (-WProd433 + WProd445 + WProd446)*dZfL1(3,1,1))/
     -          WProd448) - 
     -      5.656854249492381D0*
     -       (Cpl1*GS*Prp49**2*(-WProd445 - WProd446 + WProd447)*
     -          (WProd449 - WProd450)*dZfL1(4,1,1))/WProd448 + 
     -      (Prp50*(11.313708498984761D0*(Cpl106*dSW1*GS*tmp2) - 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*tmp2*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + 3*dZfL1(3,1,1) + dZfL1(4,1,1))
     -              )))/WProd448 + 
     -      (Prp49*(11.313708498984761D0*
     -            (Cpl106*dSW1*GS*(WProd449 - WProd450)) - 
     -           2.8284271247461903D0*
     -            (Cpl1*GS*(WProd449 - WProd450)*
     -              (4*dZe1 + 2*dZW1 + dZfL1(1,1,1) + 
     -                dZfL1(2,1,1) + dZfL1(3,1,1) + 3*dZfL1(4,1,1))
     -              )))/WProd448)

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

      endif

      end


