      subroutine virtQCD2(amp,resnr)
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

        amp(2,1,1) = Prp7*
     -    (0.8888888888888888D0*(Alfas*Cpl5*WProd16) - 
     -      1.7777777777777777D0*(Alfas*Cpl5*LInt6*WProd16) + 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt8*WProd16) + 
     -      1.7777777777777777D0*
     -       (Alfas*Cpl5*LInt10*WProd16*WProd17) + 
     -      1.7777777777777777D0*
     -       (Alfas*Cpl5*LInt7*WProd16*WProd17) + 
     -      1.7777777777777777D0*(Alfas*Cpl5*LInt9*WProd16*WProd17))+
     -     Prp8*(0.8888888888888888D0*(Alfas*Cpl6*WProd16) - 
     -      1.7777777777777777D0*(Alfas*Cpl6*LInt6*WProd16) + 
     -      3.5555555555555554D0*(Alfas*Cpl6*LInt8*WProd16) + 
     -      1.7777777777777777D0*
     -       (Alfas*Cpl6*LInt10*WProd16*WProd17) + 
     -      1.7777777777777777D0*
     -       (Alfas*Cpl6*LInt7*WProd16*WProd17) + 
     -      1.7777777777777777D0*(Alfas*Cpl6*LInt9*WProd16*WProd17))

#ifdef DEBUG
	print *, 'amp(2,1,1) =', amp(2,1,1)
#endif

        tmp170 = 1.7777777777777777D0*
     -    (Alfas*Cpl6*WProd17*WProd18) - 
     -   0.8888888888888888D0*(Alfas*Cpl7*WProd17*WProd18)
        amp(2,1,2) = Prp8*
     -    (LInt10*tmp170 + LInt7*tmp170 + LInt9*tmp170 + 
     -      0.8888888888888888D0*(Alfas*Cpl6*WProd18) - 
     -      0.4444444444444444D0*(Alfas*Cpl7*WProd18) + 
     -      LInt6*(-(1.7777777777777777D0*(Alfas*Cpl6*WProd18)) + 
     -         0.8888888888888888D0*(Alfas*Cpl7*WProd18)) + 
     -      LInt8*(3.5555555555555554D0*(Alfas*Cpl6*WProd18) - 
     -         1.7777777777777777D0*(Alfas*Cpl7*WProd18))) + 
     -   Prp7*(0.8888888888888888D0*(Alfas*Cpl5*WProd18) - 
     -      1.7777777777777777D0*(Alfas*Cpl5*LInt6*WProd18) + 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt8*WProd18) + 
     -      1.7777777777777777D0*
     -       (Alfas*Cpl5*LInt10*WProd17*WProd18) + 
     -      1.7777777777777777D0*
     -       (Alfas*Cpl5*LInt7*WProd17*WProd18) + 
     -      1.7777777777777777D0*(Alfas*Cpl5*LInt9*WProd17*WProd18))

#ifdef DEBUG
	print *, 'amp(2,1,2) =', amp(2,1,2)
#endif

        tmp171 = 1.7777777777777777D0*
     -    (Alfas*Cpl6*WProd17*WProd19) - 
     -   2.6666666666666665D0*(Alfas*Cpl7*WProd17*WProd19)
        amp(2,1,3) = Prp8*
     -    (LInt10*tmp171 + LInt7*tmp171 + LInt9*tmp171 + 
     -      0.8888888888888888D0*(Alfas*Cpl6*WProd19) - 
     -      1.3333333333333333D0*(Alfas*Cpl7*WProd19) + 
     -      LInt6*(-(1.7777777777777777D0*(Alfas*Cpl6*WProd19)) + 
     -         2.6666666666666665D0*(Alfas*Cpl7*WProd19)) + 
     -      LInt8*(3.5555555555555554D0*(Alfas*Cpl6*WProd19) - 
     -         5.333333333333333D0*(Alfas*Cpl7*WProd19))) + 
     -   Prp7*(0.8888888888888888D0*(Alfas*Cpl5*WProd19) - 
     -      1.7777777777777777D0*(Alfas*Cpl5*LInt6*WProd19) + 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt8*WProd19) + 
     -      1.7777777777777777D0*
     -       (Alfas*Cpl5*LInt10*WProd17*WProd19) + 
     -      1.7777777777777777D0*
     -       (Alfas*Cpl5*LInt7*WProd17*WProd19) + 
     -      1.7777777777777777D0*(Alfas*Cpl5*LInt9*WProd17*WProd19))

#ifdef DEBUG
	print *, 'amp(2,1,3) =', amp(2,1,3)
#endif

        tmp172 = 1.7777777777777777D0*
     -    (Alfas*Cpl6*WProd17*WProd20) - 
     -   3.5555555555555554D0*(Alfas*Cpl7*WProd17*WProd20) + 
     -   1.3333333333333333D0*(Alfas*Cpl8*WProd17*WProd20)
        amp(2,1,4) = Prp8*
     -    (LInt10*tmp172 + LInt7*tmp172 + LInt9*tmp172 + 
     -      0.8888888888888888D0*(Alfas*Cpl6*WProd20) - 
     -      1.7777777777777777D0*(Alfas*Cpl7*WProd20) + 
     -      0.6666666666666666D0*(Alfas*Cpl8*WProd20) + 
     -      LInt6*(-(1.7777777777777777D0*(Alfas*Cpl6*WProd20)) + 
     -         3.5555555555555554D0*(Alfas*Cpl7*WProd20) - 
     -         1.3333333333333333D0*(Alfas*Cpl8*WProd20)) + 
     -      LInt8*(3.5555555555555554D0*(Alfas*Cpl6*WProd20) - 
     -         7.111111111111111D0*(Alfas*Cpl7*WProd20) + 
     -         2.6666666666666665D0*(Alfas*Cpl8*WProd20))) + 
     -   Prp7*(0.8888888888888888D0*(Alfas*Cpl5*WProd20) - 
     -      1.7777777777777777D0*(Alfas*Cpl5*LInt6*WProd20) + 
     -      3.5555555555555554D0*(Alfas*Cpl5*LInt8*WProd20) + 
     -      1.7777777777777777D0*
     -       (Alfas*Cpl5*LInt10*WProd17*WProd20) + 
     -      1.7777777777777777D0*
     -       (Alfas*Cpl5*LInt7*WProd17*WProd20) + 
     -      1.7777777777777777D0*(Alfas*Cpl5*LInt9*WProd17*WProd20))

#ifdef DEBUG
	print *, 'amp(2,1,4) =', amp(2,1,4)
#endif

      endif

      end


