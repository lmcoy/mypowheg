      subroutine virtEWphotonicISR2(amp,resnr)
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

        amp(2,1,1) = Prp11*
     -    (0.07407407407407407D0*(Cpl92*WProd26) - 
     -      0.14814814814814814D0*(Cpl92*LInt86*WProd26) + 
     -      0.2962962962962963D0*(Cpl92*LInt88*WProd26) + 
     -      0.14814814814814814D0*(Cpl92*LInt87*WProd26*WProd27) + 
     -      0.14814814814814814D0*(Cpl92*LInt89*WProd26*WProd27) + 
     -      0.14814814814814814D0*(Cpl92*LInt90*WProd26*WProd27)) + 
     -   Prp12*(0.07407407407407407D0*(Cpl93*WProd26) - 
     -      0.14814814814814814D0*(Cpl93*LInt86*WProd26) + 
     -      0.2962962962962963D0*(Cpl93*LInt88*WProd26) + 
     -      0.14814814814814814D0*(Cpl93*LInt87*WProd26*WProd27) + 
     -      0.14814814814814814D0*(Cpl93*LInt89*WProd26*WProd27) + 
     -      0.14814814814814814D0*(Cpl93*LInt90*WProd26*WProd27))

#ifdef DEBUG
	print *, 'amp(2,1,1) =', amp(2,1,1)
#endif

        tmp164 = 0.14814814814814814D0*(Cpl93*WProd27*WProd28) - 
     -   0.07407407407407407D0*(Cpl94*WProd27*WProd28)
        amp(2,1,2) = Prp12*
     -    (LInt87*tmp164 + LInt89*tmp164 + LInt90*tmp164 + 
     -      0.07407407407407407D0*(Cpl93*WProd28) - 
     -      0.037037037037037035D0*(Cpl94*WProd28) + 
     -      LInt86*(-(0.14814814814814814D0*(Cpl93*WProd28)) + 
     -         0.07407407407407407D0*(Cpl94*WProd28)) + 
     -      LInt88*(0.2962962962962963D0*(Cpl93*WProd28) - 
     -         0.14814814814814814D0*(Cpl94*WProd28))) + 
     -   Prp11*(0.07407407407407407D0*(Cpl92*WProd28) - 
     -      0.14814814814814814D0*(Cpl92*LInt86*WProd28) + 
     -      0.2962962962962963D0*(Cpl92*LInt88*WProd28) + 
     -      0.14814814814814814D0*(Cpl92*LInt87*WProd27*WProd28) + 
     -      0.14814814814814814D0*(Cpl92*LInt89*WProd27*WProd28) + 
     -      0.14814814814814814D0*(Cpl92*LInt90*WProd27*WProd28))

#ifdef DEBUG
	print *, 'amp(2,1,2) =', amp(2,1,2)
#endif

        tmp165 = 0.14814814814814814D0*(Cpl93*WProd27*WProd29) - 
     -   0.2222222222222222D0*(Cpl94*WProd27*WProd29)
        amp(2,1,3) = Prp12*
     -    (LInt87*tmp165 + LInt89*tmp165 + LInt90*tmp165 + 
     -      0.07407407407407407D0*(Cpl93*WProd29) - 
     -      0.1111111111111111D0*(Cpl94*WProd29) + 
     -      LInt86*(-(0.14814814814814814D0*(Cpl93*WProd29)) + 
     -         0.2222222222222222D0*(Cpl94*WProd29)) + 
     -      LInt88*(0.2962962962962963D0*(Cpl93*WProd29) - 
     -         0.4444444444444444D0*(Cpl94*WProd29))) + 
     -   Prp11*(0.07407407407407407D0*(Cpl92*WProd29) - 
     -      0.14814814814814814D0*(Cpl92*LInt86*WProd29) + 
     -      0.2962962962962963D0*(Cpl92*LInt88*WProd29) + 
     -      0.14814814814814814D0*(Cpl92*LInt87*WProd27*WProd29) + 
     -      0.14814814814814814D0*(Cpl92*LInt89*WProd27*WProd29) + 
     -      0.14814814814814814D0*(Cpl92*LInt90*WProd27*WProd29))

#ifdef DEBUG
	print *, 'amp(2,1,3) =', amp(2,1,3)
#endif

        tmp166 = 0.14814814814814814D0*(Cpl93*WProd27*WProd30) - 
     -   0.2962962962962963D0*(Cpl94*WProd27*WProd30) + 
     -   0.1111111111111111D0*(Cpl95*WProd27*WProd30)
        amp(2,1,4) = Prp12*
     -    (LInt87*tmp166 + LInt89*tmp166 + LInt90*tmp166 + 
     -      0.07407407407407407D0*(Cpl93*WProd30) - 
     -      0.14814814814814814D0*(Cpl94*WProd30) + 
     -      0.05555555555555555D0*(Cpl95*WProd30) + 
     -      LInt86*(-(0.14814814814814814D0*(Cpl93*WProd30)) + 
     -         0.2962962962962963D0*(Cpl94*WProd30) - 
     -         0.1111111111111111D0*(Cpl95*WProd30)) + 
     -      LInt88*(0.2962962962962963D0*(Cpl93*WProd30) - 
     -         0.5925925925925926D0*(Cpl94*WProd30) + 
     -         0.2222222222222222D0*(Cpl95*WProd30))) + 
     -   Prp11*(0.07407407407407407D0*(Cpl92*WProd30) - 
     -      0.14814814814814814D0*(Cpl92*LInt86*WProd30) + 
     -      0.2962962962962963D0*(Cpl92*LInt88*WProd30) + 
     -      0.14814814814814814D0*(Cpl92*LInt87*WProd27*WProd30) + 
     -      0.14814814814814814D0*(Cpl92*LInt89*WProd27*WProd30) + 
     -      0.14814814814814814D0*(Cpl92*LInt90*WProd27*WProd30))

#ifdef DEBUG
	print *, 'amp(2,1,4) =', amp(2,1,4)
#endif

      endif

      end


