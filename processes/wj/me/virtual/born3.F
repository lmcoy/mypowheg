      subroutine born3(amp,resnr)
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

        amp(3,1,3) = -(5.656854249492381D0*
     -      (Cpl1*GS*Prp10*Prp9*
     -         (WProd46*WProd47 - WProd48 + WProd49))/WProd45) - 
     -   5.656854249492381D0*
     -    (Cpl1*GS*Prp10*Prp11*WProd47*
     -       (WProd46 - WProd52 - WProd53))/WProd45

#ifdef DEBUG
	print *, 'amp(3,1,3) =', amp(3,1,3)
#endif

        amp(3,1,1) = -(5.656854249492381D0*
     -      (Cpl1*GS*Prp10*Prp9*
     -         (WProd47*WProd56 - WProd57 + WProd58 - WProd59))/
     -       WProd55) - 5.656854249492381D0*
     -    (Cpl1*GS*Prp10*Prp11*
     -       (WProd47*WProd56 - WProd57 - WProd47*WProd60 - 
     -         WProd47*WProd61))/WProd55

#ifdef DEBUG
	print *, 'amp(3,1,1) =', amp(3,1,1)
#endif

      endif

      end


