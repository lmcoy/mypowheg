      subroutine virtQCD1(amp,resnr)
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

        amp(1,1,1) = Prp5*
     -    (-(1.7777777777777777D0*(Alfas*Cpl5*WProd11)) + 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt1*WProd11) - 
     -      7.111111111111111D0*(Alfas*Cpl5*LInt3*WProd11) - 
     -      3.5555555555555554D0*
     -       (Alfas*Cpl5*LInt2*WProd11*WProd12) - 
     -      3.5555555555555554D0*
     -       (Alfas*Cpl5*LInt4*WProd11*WProd12) - 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt5*WProd11*WProd12))+
     -     Prp6*(-(1.7777777777777777D0*(Alfas*Cpl6*WProd11)) + 
     -      3.5555555555555554D0*(Alfas*Cpl6*LInt1*WProd11) - 
     -      7.111111111111111D0*(Alfas*Cpl6*LInt3*WProd11) - 
     -      3.5555555555555554D0*
     -       (Alfas*Cpl6*LInt2*WProd11*WProd12) - 
     -      3.5555555555555554D0*
     -       (Alfas*Cpl6*LInt4*WProd11*WProd12) - 
     -      3.5555555555555554D0*(Alfas*Cpl6*LInt5*WProd11*WProd12))

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        tmp167 = -(3.5555555555555554D0*
     -      (Alfas*Cpl6*WProd12*WProd13)) + 
     -   1.7777777777777777D0*(Alfas*Cpl7*WProd12*WProd13)
        amp(1,1,2) = Prp6*
     -    (LInt2*tmp167 + LInt4*tmp167 + LInt5*tmp167 - 
     -      1.7777777777777777D0*(Alfas*Cpl6*WProd13) + 
     -      0.8888888888888888D0*(Alfas*Cpl7*WProd13) + 
     -      LInt1*(3.5555555555555554D0*(Alfas*Cpl6*WProd13) - 
     -         1.7777777777777777D0*(Alfas*Cpl7*WProd13)) + 
     -      LInt3*(-(7.111111111111111D0*(Alfas*Cpl6*WProd13)) + 
     -         3.5555555555555554D0*(Alfas*Cpl7*WProd13))) + 
     -   Prp5*(-(1.7777777777777777D0*(Alfas*Cpl5*WProd13)) + 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt1*WProd13) - 
     -      7.111111111111111D0*(Alfas*Cpl5*LInt3*WProd13) - 
     -      3.5555555555555554D0*
     -       (Alfas*Cpl5*LInt2*WProd12*WProd13) - 
     -      3.5555555555555554D0*
     -       (Alfas*Cpl5*LInt4*WProd12*WProd13) - 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt5*WProd12*WProd13))

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

        tmp168 = -(3.5555555555555554D0*
     -      (Alfas*Cpl6*WProd12*WProd14)) + 
     -   2.6666666666666665D0*(Alfas*Cpl7*WProd12*WProd14)
        amp(1,1,3) = Prp6*
     -    (LInt2*tmp168 + LInt4*tmp168 + LInt5*tmp168 - 
     -      1.7777777777777777D0*(Alfas*Cpl6*WProd14) + 
     -      1.3333333333333333D0*(Alfas*Cpl7*WProd14) + 
     -      LInt1*(3.5555555555555554D0*(Alfas*Cpl6*WProd14) - 
     -         2.6666666666666665D0*(Alfas*Cpl7*WProd14)) + 
     -      LInt3*(-(7.111111111111111D0*(Alfas*Cpl6*WProd14)) + 
     -         5.333333333333333D0*(Alfas*Cpl7*WProd14))) + 
     -   Prp5*(-(1.7777777777777777D0*(Alfas*Cpl5*WProd14)) + 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt1*WProd14) - 
     -      7.111111111111111D0*(Alfas*Cpl5*LInt3*WProd14) - 
     -      3.5555555555555554D0*
     -       (Alfas*Cpl5*LInt2*WProd12*WProd14) - 
     -      3.5555555555555554D0*
     -       (Alfas*Cpl5*LInt4*WProd12*WProd14) - 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt5*WProd12*WProd14))

#ifdef DEBUG
	print *, 'amp(1,1,3) =', amp(1,1,3)
#endif

        tmp169 = -(3.5555555555555554D0*
     -      (Alfas*Cpl6*WProd12*WProd15)) + 
     -   4.444444444444445D0*(Alfas*Cpl7*WProd12*WProd15) - 
     -   1.3333333333333333D0*(Alfas*Cpl8*WProd12*WProd15)
        amp(1,1,4) = Prp6*
     -    (LInt2*tmp169 + LInt4*tmp169 + LInt5*tmp169 - 
     -      1.7777777777777777D0*(Alfas*Cpl6*WProd15) + 
     -      2.2222222222222223D0*(Alfas*Cpl7*WProd15) - 
     -      0.6666666666666666D0*(Alfas*Cpl8*WProd15) + 
     -      LInt1*(3.5555555555555554D0*(Alfas*Cpl6*WProd15) - 
     -         4.444444444444445D0*(Alfas*Cpl7*WProd15) + 
     -         1.3333333333333333D0*(Alfas*Cpl8*WProd15)) + 
     -      LInt3*(-(7.111111111111111D0*(Alfas*Cpl6*WProd15)) + 
     -         8.88888888888889D0*(Alfas*Cpl7*WProd15) - 
     -         2.6666666666666665D0*(Alfas*Cpl8*WProd15))) + 
     -   Prp5*(-(1.7777777777777777D0*(Alfas*Cpl5*WProd15)) + 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt1*WProd15) - 
     -      7.111111111111111D0*(Alfas*Cpl5*LInt3*WProd15) - 
     -      3.5555555555555554D0*
     -       (Alfas*Cpl5*LInt2*WProd12*WProd15) - 
     -      3.5555555555555554D0*
     -       (Alfas*Cpl5*LInt4*WProd12*WProd15) - 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt5*WProd12*WProd15))

#ifdef DEBUG
	print *, 'amp(1,1,4) =', amp(1,1,4)
#endif

      endif
      end


