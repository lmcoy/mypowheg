      subroutine ampsfacctEWphotonicISR2(amp,resnr)
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
     -      (Cpl141*Prp23*WProd122)) - 
     -   2.6666666666666665D0*(Cpl142*Prp24*WProd122)

#ifdef DEBUG
	print *, 'amp(2,1,1) =', amp(2,1,1)
#endif

        amp(2,1,2) = Prp24*
     -    (-(2.6666666666666665D0*(Cpl142*WProd124)) + 
     -      1.3333333333333333D0*(Cpl143*WProd124)) - 
     -   2.6666666666666665D0*(Cpl141*Prp23*WProd124)

#ifdef DEBUG
	print *, 'amp(2,1,2) =', amp(2,1,2)
#endif

        amp(2,1,3) = Prp24*
     -    (4*Cpl145*WProd125 - 
     -      2.6666666666666665D0*(Cpl146*WProd125)) - 
     -   2.6666666666666665D0*(Cpl144*Prp23*WProd125)

#ifdef DEBUG
	print *, 'amp(2,1,3) =', amp(2,1,3)
#endif

        amp(2,1,4) = Prp24*
     -    (-2*Cpl147*WProd126 + 
     -      5.333333333333333D0*(Cpl145*WProd126) - 
     -      2.6666666666666665D0*(Cpl146*WProd126)) - 
     -   2.6666666666666665D0*(Cpl144*Prp23*WProd126)

#ifdef DEBUG
	print *, 'amp(2,1,4) =', amp(2,1,4)
#endif

      endif

      end


