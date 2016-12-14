      subroutine virtEWphotonicISR1(amp,resnr)
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

        amp(1,1,1) = Prp9*
     -    (-(0.5925925925925926D0*(Cpl92*WProd21)) + 
     -      1.1851851851851851D0*(Cpl92*LInt81*WProd21) - 
     -      2.3703703703703702D0*(Cpl92*LInt83*WProd21) - 
     -      1.1851851851851851D0*(Cpl92*LInt82*WProd21*WProd22) - 
     -      1.1851851851851851D0*(Cpl92*LInt84*WProd21*WProd22) - 
     -      1.1851851851851851D0*(Cpl92*LInt85*WProd21*WProd22)) + 
     -   Prp10*(-(0.5925925925925926D0*(Cpl93*WProd21)) + 
     -      1.1851851851851851D0*(Cpl93*LInt81*WProd21) - 
     -      2.3703703703703702D0*(Cpl93*LInt83*WProd21) - 
     -      1.1851851851851851D0*(Cpl93*LInt82*WProd21*WProd22) - 
     -      1.1851851851851851D0*(Cpl93*LInt84*WProd21*WProd22) - 
     -      1.1851851851851851D0*(Cpl93*LInt85*WProd21*WProd22))

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        tmp161 = -(1.1851851851851851D0*(Cpl93*WProd22*WProd23)) + 
     -   0.5925925925925926D0*(Cpl94*WProd22*WProd23)
        amp(1,1,2) = Prp10*
     -    (LInt82*tmp161 + LInt84*tmp161 + LInt85*tmp161 - 
     -      0.5925925925925926D0*(Cpl93*WProd23) + 
     -      0.2962962962962963D0*(Cpl94*WProd23) + 
     -      LInt81*(1.1851851851851851D0*(Cpl93*WProd23) - 
     -         0.5925925925925926D0*(Cpl94*WProd23)) + 
     -      LInt83*(-(2.3703703703703702D0*(Cpl93*WProd23)) + 
     -         1.1851851851851851D0*(Cpl94*WProd23))) + 
     -   Prp9*(-(0.5925925925925926D0*(Cpl92*WProd23)) + 
     -      1.1851851851851851D0*(Cpl92*LInt81*WProd23) - 
     -      2.3703703703703702D0*(Cpl92*LInt83*WProd23) - 
     -      1.1851851851851851D0*(Cpl92*LInt82*WProd22*WProd23) - 
     -      1.1851851851851851D0*(Cpl92*LInt84*WProd22*WProd23) - 
     -      1.1851851851851851D0*(Cpl92*LInt85*WProd22*WProd23))

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

        tmp162 = -(1.1851851851851851D0*(Cpl93*WProd22*WProd24)) + 
     -   0.8888888888888888D0*(Cpl94*WProd22*WProd24)
        amp(1,1,3) = Prp10*
     -    (LInt82*tmp162 + LInt84*tmp162 + LInt85*tmp162 - 
     -      0.5925925925925926D0*(Cpl93*WProd24) + 
     -      0.4444444444444444D0*(Cpl94*WProd24) + 
     -      LInt81*(1.1851851851851851D0*(Cpl93*WProd24) - 
     -         0.8888888888888888D0*(Cpl94*WProd24)) + 
     -      LInt83*(-(2.3703703703703702D0*(Cpl93*WProd24)) + 
     -         1.7777777777777777D0*(Cpl94*WProd24))) + 
     -   Prp9*(-(0.5925925925925926D0*(Cpl92*WProd24)) + 
     -      1.1851851851851851D0*(Cpl92*LInt81*WProd24) - 
     -      2.3703703703703702D0*(Cpl92*LInt83*WProd24) - 
     -      1.1851851851851851D0*(Cpl92*LInt82*WProd22*WProd24) - 
     -      1.1851851851851851D0*(Cpl92*LInt84*WProd22*WProd24) - 
     -      1.1851851851851851D0*(Cpl92*LInt85*WProd22*WProd24))

#ifdef DEBUG
	print *, 'amp(1,1,3) =', amp(1,1,3)
#endif

        tmp163 = -(1.1851851851851851D0*(Cpl93*WProd22*WProd25)) + 
     -   1.4814814814814814D0*(Cpl94*WProd22*WProd25) - 
     -   0.4444444444444444D0*(Cpl95*WProd22*WProd25)
        amp(1,1,4) = Prp10*
     -    (LInt82*tmp163 + LInt84*tmp163 + LInt85*tmp163 - 
     -      0.5925925925925926D0*(Cpl93*WProd25) + 
     -      0.7407407407407407D0*(Cpl94*WProd25) - 
     -      0.2222222222222222D0*(Cpl95*WProd25) + 
     -      LInt81*(1.1851851851851851D0*(Cpl93*WProd25) - 
     -         1.4814814814814814D0*(Cpl94*WProd25) + 
     -         0.4444444444444444D0*(Cpl95*WProd25)) + 
     -      LInt83*(-(2.3703703703703702D0*(Cpl93*WProd25)) + 
     -         2.962962962962963D0*(Cpl94*WProd25) - 
     -         0.8888888888888888D0*(Cpl95*WProd25))) + 
     -   Prp9*(-(0.5925925925925926D0*(Cpl92*WProd25)) + 
     -      1.1851851851851851D0*(Cpl92*LInt81*WProd25) - 
     -      2.3703703703703702D0*(Cpl92*LInt83*WProd25) - 
     -      1.1851851851851851D0*(Cpl92*LInt82*WProd22*WProd25) - 
     -      1.1851851851851851D0*(Cpl92*LInt84*WProd22*WProd25) - 
     -      1.1851851851851851D0*(Cpl92*LInt85*WProd22*WProd25))

#ifdef DEBUG
	print *, 'amp(1,1,4) =', amp(1,1,4)
#endif

      endif

      end


