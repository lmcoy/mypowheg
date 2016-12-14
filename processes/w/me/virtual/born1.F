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
*     procnr=1:   { up  , ~do  ,  ve  , ~el  }
*
*   Helicities:
*     helnr=1:   {-1, 1, -1, 1}
*
*   Color Structures:
*     colnr=1:   SUNT[Col1, Col2]


      if(resnr.eq.0)then

          amp(1,1,1) = -4*Cpl1*Prp1*WProd1

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

      endif

      end


