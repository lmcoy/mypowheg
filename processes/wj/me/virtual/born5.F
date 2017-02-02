      subroutine born5(amp,resnr)
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

        amp(5,3,4) = -(5.656854249492381D0*
     -      (Cpl1*GS*Prp17*Prp18*WProd87*WProd88)/WProd86) + 
     -   5.656854249492381D0*
     -    (Cpl1*GS*Prp18*Prp19*
     -       (-(WProd87*WProd88) + WProd91 - WProd92 + 
     -         WProd88*WProd93 + WProd88*WProd94))/WProd86

#ifdef DEBUG
	print *, 'amp(5,3,4) =', amp(5,3,4)
#endif

        amp(5,3,5) = -(5.656854249492381D0*
     -      (Cpl1*GS*Prp17*Prp18*(WProd88*WProd97 + WProd98))/
     -       WProd96) - 5.656854249492381D0*
     -    (Cpl1*GS*Prp18*Prp19*
     -       (-WProd100 - WProd101*WProd88 - WProd102*WProd88 + 
     -         WProd88*WProd97 + WProd98 + WProd99))/WProd96

#ifdef DEBUG
	print *, 'amp(5,3,5) =', amp(5,3,5)
#endif

      endif

      end


