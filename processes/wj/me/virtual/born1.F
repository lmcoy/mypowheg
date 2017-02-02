      subroutine born1(amp,resnr)
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

        amp(1,1,1) = Prp1*
     -   (-(5.656854249492381D0*
     -        (Cpl1*GS*Prp3*(WProd11 - WProd7))/WProd9) - 
     -     5.656854249492381D0*
     -      (Cpl1*GS*Prp2*(WProd2*WProd3 - WProd4 + 
     -           WProd2*WProd5 - WProd6 - WProd7 - WProd8))/WProd9)

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        amp(1,1,2) = Prp1*
     -   (-(5.656854249492381D0*
     -        (Cpl1*GS*Prp2*
     -           (-WProd15 - WProd17 - WProd18 - WProd19 + 
     -             WProd14*WProd2 + WProd16*WProd2))/WProd13) - 
     -     5.656854249492381D0*
     -      (Cpl1*GS*Prp3*(-WProd18 + WProd20))/WProd13)

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

      endif

      end


