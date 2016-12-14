      subroutine realQED4(amp,resnr)
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
*     procnr=1:   { up  , ~do  ,  ve  , ~el  ,  ga  }
*     procnr=2:   { up  , ~do  ,  ve  , ~el  ,  gl  }
*     procnr=3:   { up  ,  ga  ,  ve  , ~el  ,  do  }
*     procnr=4:   { up  ,  gl  ,  ve  , ~el  ,  do  }
*     procnr=5:   {~do  ,  ga  ,  ve  , ~el  , ~up  }
*     procnr=6:   {~do  ,  gl  ,  ve  , ~el  , ~up  }
*
*   Helicities:
*     helnr=1:   {-1, 1, -1, 1, -1}
*     helnr=2:   {-1, 1, -1, 1, 1}
*     helnr=3:   {-1, -1, -1, 1, -1}
*     helnr=4:   {1, -1, -1, 1, 1}
*     helnr=5:   {1, 1, -1, 1, 1}
*
*   Color Structures:
*     colnr=1:   SUNT[Col1, Col2]
*     colnr=2:   SUNT[Glu1, Col2, Col1]
*     colnr=3:   SUNT[Glu1, Col1, Col2]


      if(resnr.eq.0)then

        amp(4,2,3) = -(5.656854249492381D0*
     -      (Cpl1*GS*Prp19*Prp20*
     -         (WProd86*WProd87 - WProd88 + WProd89))/WProd85) - 
     -   5.656854249492381D0*
     -    (Cpl1*GS*Prp20*Prp21*WProd87*
     -       (WProd86 - WProd92 - WProd93))/WProd85

#ifdef DEBUG
	print *, 'amp(4,2,3) =', amp(4,2,3)
#endif

        amp(4,2,1) = -(5.656854249492381D0*
     -      (Cpl1*GS*Prp20*Prp21*
     -         (-(WProd100*WProd87) - WProd101*WProd87 + 
     -           WProd87*WProd96 - WProd97))/WProd95) - 
     -   5.656854249492381D0*
     -    (Cpl1*GS*Prp19*Prp20*
     -       (WProd87*WProd96 - WProd97 + WProd98 - WProd99))/
     -     WProd95

#ifdef DEBUG
	print *, 'amp(4,2,1) =', amp(4,2,1)
#endif

      endif

      end


