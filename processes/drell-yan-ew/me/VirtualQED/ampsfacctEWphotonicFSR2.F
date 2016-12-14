      subroutine ampsfacctEWphotonicFSR2(amp,resnr)
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

        amp(2,1,1) = -(2.6666666666666665D0*
     -      (Cpl148*Prp27*WProd132)) - 
     -   2.6666666666666665D0*(Cpl149*Prp28*WProd132)

#ifdef DEBUG
	print *, 'amp(2,1,1) =', amp(2,1,1)
#endif

        amp(2,1,2) = Prp28*
     -    (1.3333333333333333D0*(Cpl151*WProd134) - 
     -      2.6666666666666665D0*(Cpl152*WProd134)) - 
     -   2.6666666666666665D0*(Cpl150*Prp27*WProd134)

#ifdef DEBUG
	print *, 'amp(2,1,2) =', amp(2,1,2)
#endif

        amp(2,1,3) = Prp28*
     -    (4*Cpl153*WProd135 - 
     -      2.6666666666666665D0*(Cpl149*WProd135)) - 
     -   2.6666666666666665D0*(Cpl148*Prp27*WProd135)

#ifdef DEBUG
	print *, 'amp(2,1,3) =', amp(2,1,3)
#endif

        amp(2,1,4) = Prp28*
     -    (-2*Cpl154*WProd136 + 
     -      5.333333333333333D0*(Cpl151*WProd136) - 
     -      2.6666666666666665D0*(Cpl152*WProd136)) - 
     -   2.6666666666666665D0*(Cpl150*Prp27*WProd136)

#ifdef DEBUG
	print *, 'amp(2,1,4) =', amp(2,1,4)
#endif

      endif

      end


