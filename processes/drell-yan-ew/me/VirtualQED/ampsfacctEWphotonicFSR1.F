      subroutine ampsfacctEWphotonicFSR1(amp,resnr)
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
*     procnr=1:   {~up  ,  up  , ~el  ,  el  }
*     procnr=2:   {~do  ,  do  , ~el  ,  el  }
*
*   Helicities:
*     helnr=1:   {-1, 1, -1, 1}
*     helnr=2:   {-1, 1, 1, -1}
*     helnr=3:   {1, -1, -1, 1}
*     helnr=4:   {1, -1, 1, -1}
*
*   Color Structures:
*     colnr=1:   SUNT[Col1, Col2]


      if(resnr.eq.0)then

        amp(1,1,1) = 5.333333333333333D0*(Cpl148*Prp25*WProd127) + 
     -   5.333333333333333D0*(Cpl149*Prp26*WProd127)

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        amp(1,1,2) = Prp26*
     -    (-(2.6666666666666665D0*(Cpl151*WProd129)) + 
     -      5.333333333333333D0*(Cpl152*WProd129)) + 
     -   5.333333333333333D0*(Cpl150*Prp25*WProd129)

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

        amp(1,1,3) = Prp26*
     -    (-4*Cpl153*WProd130 + 
     -      5.333333333333333D0*(Cpl149*WProd130)) + 
     -   5.333333333333333D0*(Cpl148*Prp25*WProd130)

#ifdef DEBUG
	print *, 'amp(1,1,3) =', amp(1,1,3)
#endif

        amp(1,1,4) = Prp26*
     -    (2*Cpl154*WProd131 - 
     -      6.666666666666667D0*(Cpl151*WProd131) + 
     -      5.333333333333333D0*(Cpl152*WProd131)) + 
     -   5.333333333333333D0*(Cpl150*Prp25*WProd131)

#ifdef DEBUG
	print *, 'amp(1,1,4) =', amp(1,1,4)
#endif

      endif

      end


