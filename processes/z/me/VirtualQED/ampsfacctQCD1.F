      subroutine ampsfacctQCD1(amp,resnr)
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

        amp(1,1,1) = 5.333333333333333D0*
     -    (Alfas*Cpl1*Prp33*WProd151*dZfR1dimreg(3,1,1)) + 
     -   5.333333333333333D0*
     -    (Alfas*Cpl2*Prp34*WProd151*dZfR1dimreg(3,1,1))

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        amp(1,1,2) = Prp34*
     -    (5.333333333333333D0*
     -       (Alfas*Cpl2*WProd153*dZfR1dimreg(3,1,1)) - 
     -      2.6666666666666665D0*
     -       (Alfas*Cpl3*WProd153*dZfR1dimreg(3,1,1))) + 
     -   5.333333333333333D0*
     -    (Alfas*Cpl1*Prp33*WProd153*dZfR1dimreg(3,1,1))

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

        amp(1,1,3) = Prp34*
     -    (-4*Alfas*Cpl3*WProd154*dZfL1dimreg(3,1,1) + 
     -      5.333333333333333D0*
     -       (Alfas*Cpl2*WProd154*dZfL1dimreg(3,1,1))) + 
     -   5.333333333333333D0*
     -    (Alfas*Cpl1*Prp33*WProd154*dZfL1dimreg(3,1,1))

#ifdef DEBUG
	print *, 'amp(1,1,3) =', amp(1,1,3)
#endif

        amp(1,1,4) = Prp34*
     -    (2*Alfas*Cpl4*WProd155*dZfL1dimreg(3,1,1) + 
     -      5.333333333333333D0*
     -       (Alfas*Cpl2*WProd155*dZfL1dimreg(3,1,1)) - 
     -      6.666666666666667D0*
     -       (Alfas*Cpl3*WProd155*dZfL1dimreg(3,1,1))) + 
     -   5.333333333333333D0*
     -    (Alfas*Cpl1*Prp33*WProd155*dZfL1dimreg(3,1,1))

#ifdef DEBUG
	print *, 'amp(1,1,4) =', amp(1,1,4)
#endif

      endif

      end


