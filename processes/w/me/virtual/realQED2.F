      subroutine realQED2(amp,resnr)
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

        amp(2,2,1) = Prp11*
     -   (-(5.656854249492381D0*
     -        (Cpl1*GS*Prp12*
     -           (WProd42*WProd43 - WProd44 + WProd42*WProd45 - 
     -             WProd46 - WProd47 - WProd48))/WProd49) - 
     -     5.656854249492381D0*
     -      (Cpl1*GS*Prp13*(-WProd47 + WProd51))/WProd49)

#ifdef DEBUG
	print *, 'amp(2,2,1) =', amp(2,2,1)
#endif

        amp(2,2,2) = Prp11*
     -   (-(5.656854249492381D0*
     -        (Cpl1*GS*Prp12*
     -           (WProd42*WProd54 - WProd55 + WProd42*WProd56 - 
     -             WProd57 - WProd58 - WProd59))/WProd53) - 
     -     5.656854249492381D0*
     -      (Cpl1*GS*Prp13*(-WProd58 + WProd60))/WProd53)

#ifdef DEBUG
	print *, 'amp(2,2,2) =', amp(2,2,2)
#endif

      endif

      end


