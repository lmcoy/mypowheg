      subroutine realQED6(amp,resnr)
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

        amp(6,3,4) = -(5.656854249492381D0*
     -      (Cpl1*GS*Prp27*Prp28*WProd127*WProd128)/WProd126) + 
     -   5.656854249492381D0*
     -    (Cpl1*GS*Prp28*Prp29*
     -       (-(WProd127*WProd128) + WProd131 - WProd132 + 
     -         WProd128*WProd133 + WProd128*WProd134))/WProd126

#ifdef DEBUG
	print *, 'amp(6,3,4) =', amp(6,3,4)
#endif

        amp(6,3,5) = -(5.656854249492381D0*
     -      (Cpl1*GS*Prp27*Prp28*(WProd128*WProd137 + WProd138))/
     -       WProd136) - 5.656854249492381D0*
     -    (Cpl1*GS*Prp28*Prp29*
     -       (WProd128*WProd137 + WProd138 + WProd139 - WProd140 - 
     -         WProd128*WProd141 - WProd128*WProd142))/WProd136

#ifdef DEBUG
	print *, 'amp(6,3,5) =', amp(6,3,5)
#endif

      endif

      end


