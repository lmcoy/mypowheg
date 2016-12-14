      subroutine born2(amp,resnr)
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

        amp(2,1,1) = -(2.6666666666666665D0*(Cpl1*Prp3*WProd6)) - 
     -   2.6666666666666665D0*(Cpl2*Prp4*WProd6)

#ifdef DEBUG
	print *, 'amp(2,1,1) =', amp(2,1,1)
#endif

        amp(2,1,2) = Prp4*
     -    (-(2.6666666666666665D0*(Cpl2*WProd8)) + 
     -      1.3333333333333333D0*(Cpl3*WProd8)) - 
     -   2.6666666666666665D0*(Cpl1*Prp3*WProd8)

#ifdef DEBUG
	print *, 'amp(2,1,2) =', amp(2,1,2)
#endif

        amp(2,1,3) = Prp4*
     -    (4*Cpl3*WProd9 - 2.6666666666666665D0*(Cpl2*WProd9)) - 
     -   2.6666666666666665D0*(Cpl1*Prp3*WProd9)

#ifdef DEBUG
	print *, 'amp(2,1,3) =', amp(2,1,3)
#endif

        amp(2,1,4) = Prp4*
     -    (-2*Cpl4*WProd10 - 2.6666666666666665D0*(Cpl2*WProd10) + 
     -      5.333333333333333D0*(Cpl3*WProd10)) - 
     -   2.6666666666666665D0*(Cpl1*Prp3*WProd10)

#ifdef DEBUG
	print *, 'amp(2,1,4) =', amp(2,1,4)
#endif

      endif

      end


