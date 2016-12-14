      subroutine virtEWphotonicFSR1(amp,resnr)
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

        amp(1,1,1) = Prp13*
     -    (-(1.3333333333333333D0*(Cpl92*WProd31)) + 
     -      2.6666666666666665D0*(Cpl92*LInt91*WProd31) - 
     -      5.333333333333333D0*(Cpl92*LInt93*WProd31) + 
     -      2.6666666666666665D0*
     -       (Cpl92*LInt95*WProd31*
     -         (2*WProd32 - 3*WProd33 - 3*WProd34)) - 
     -      2.6666666666666665D0*
     -       (Cpl92*LInt92*WProd31*(WProd33 + WProd34)) - 
     -      2.6666666666666665D0*
     -       (Cpl92*LInt94*WProd31*(WProd33 + WProd34))) + 
     -   Prp14*(-(1.3333333333333333D0*(Cpl93*WProd31)) + 
     -      2.6666666666666665D0*(Cpl93*LInt91*WProd31) - 
     -      5.333333333333333D0*(Cpl93*LInt93*WProd31) + 
     -      2.6666666666666665D0*
     -       (Cpl93*LInt95*WProd31*
     -         (2*WProd32 - 3*WProd33 - 3*WProd34)) - 
     -      2.6666666666666665D0*
     -       (Cpl93*LInt92*WProd31*(WProd33 + WProd34)) - 
     -      2.6666666666666665D0*
     -       (Cpl93*LInt94*WProd31*(WProd33 + WProd34)))

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        tmp74 = -(2.6666666666666665D0*
     -      (Cpl93*(WProd33 + WProd34)*WProd35)) + 
     -   1.3333333333333333D0*(Cpl94*(WProd33 + WProd34)*WProd35)
        amp(1,1,2) = Prp14*
     -    (LInt92*tmp74 + LInt94*tmp74 - 
     -      1.3333333333333333D0*(Cpl93*WProd35) + 
     -      0.6666666666666666D0*(Cpl94*WProd35) + 
     -      LInt91*(2.6666666666666665D0*(Cpl93*WProd35) - 
     -         1.3333333333333333D0*(Cpl94*WProd35)) + 
     -      LInt93*(-(5.333333333333333D0*(Cpl93*WProd35)) + 
     -         2.6666666666666665D0*(Cpl94*WProd35)) + 
     -      LInt95*(2.6666666666666665D0*
     -          (Cpl93*(2*WProd32 - 3*WProd33 - 3*WProd34)*WProd35)
     -           - 1.3333333333333333D0*
     -          (Cpl94*(2*WProd32 - 3*WProd33 - 3*WProd34)*WProd35)
     -         )) + Prp13*(-(1.3333333333333333D0*(Cpl92*WProd35)) + 
     -      2.6666666666666665D0*(Cpl92*LInt91*WProd35) - 
     -      5.333333333333333D0*(Cpl92*LInt93*WProd35) + 
     -      2.6666666666666665D0*
     -       (Cpl92*LInt95*(2*WProd32 - 3*WProd33 - 3*WProd34)*
     -         WProd35) - 2.6666666666666665D0*
     -       (Cpl92*LInt92*(WProd33 + WProd34)*WProd35) - 
     -      2.6666666666666665D0*
     -       (Cpl92*LInt94*(WProd33 + WProd34)*WProd35))

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

        tmp75 = 2*Cpl94*(WProd33 + WProd34)*WProd36 - 
     -   2.6666666666666665D0*(Cpl93*(WProd33 + WProd34)*WProd36)
        amp(1,1,3) = Prp14*
     -    (LInt92*tmp75 + LInt94*tmp75 + Cpl94*WProd36 - 
     -      1.3333333333333333D0*(Cpl93*WProd36) + 
     -      LInt91*(-2*Cpl94*WProd36 + 
     -         2.6666666666666665D0*(Cpl93*WProd36)) + 
     -      LInt93*(4*Cpl94*WProd36 - 
     -         5.333333333333333D0*(Cpl93*WProd36)) + 
     -      LInt95*(-2*Cpl94*(2*WProd32 - 3*WProd33 - 3*WProd34)*
     -          WProd36 + 2.6666666666666665D0*
     -          (Cpl93*(2*WProd32 - 3*WProd33 - 3*WProd34)*WProd36)
     -         )) + Prp13*(-(1.3333333333333333D0*(Cpl92*WProd36)) + 
     -      2.6666666666666665D0*(Cpl92*LInt91*WProd36) - 
     -      5.333333333333333D0*(Cpl92*LInt93*WProd36) + 
     -      2.6666666666666665D0*
     -       (Cpl92*LInt95*(2*WProd32 - 3*WProd33 - 3*WProd34)*
     -         WProd36) - 2.6666666666666665D0*
     -       (Cpl92*LInt92*(WProd33 + WProd34)*WProd36) - 
     -      2.6666666666666665D0*
     -       (Cpl92*LInt94*(WProd33 + WProd34)*WProd36))

#ifdef DEBUG
	print *, 'amp(1,1,3) =', amp(1,1,3)
#endif

        tmp76 = -(Cpl95*(WProd33 + WProd34)*WProd37) - 
     -   2.6666666666666665D0*(Cpl93*(WProd33 + WProd34)*WProd37) + 
     -   3.3333333333333335D0*(Cpl94*(WProd33 + WProd34)*WProd37)
        amp(1,1,4) = Prp14*
     -    (LInt92*tmp76 + LInt94*tmp76 - 
     -      1.3333333333333333D0*(Cpl93*WProd37) + 
     -      1.6666666666666667D0*(Cpl94*WProd37) + 
     -      LInt91*(Cpl95*WProd37 + 
     -         2.6666666666666665D0*(Cpl93*WProd37) - 
     -         3.3333333333333335D0*(Cpl94*WProd37)) + 
     -      LInt93*(-2*Cpl95*WProd37 - 
     -         5.333333333333333D0*(Cpl93*WProd37) + 
     -         6.666666666666667D0*(Cpl94*WProd37)) - 
     -      0.5D0*(Cpl95*WProd37) + 
     -      LInt95*(Cpl95*(2*WProd32 - 3*WProd33 - 3*WProd34)*
     -          WProd37 + 2.6666666666666665D0*
     -          (Cpl93*(2*WProd32 - 3*WProd33 - 3*WProd34)*WProd37)
     -           - 3.3333333333333335D0*
     -          (Cpl94*(2*WProd32 - 3*WProd33 - 3*WProd34)*WProd37)
     -         )) + Prp13*(-(1.3333333333333333D0*(Cpl92*WProd37)) + 
     -      2.6666666666666665D0*(Cpl92*LInt91*WProd37) - 
     -      5.333333333333333D0*(Cpl92*LInt93*WProd37) + 
     -      2.6666666666666665D0*
     -       (Cpl92*LInt95*(2*WProd32 - 3*WProd33 - 3*WProd34)*
     -         WProd37) - 2.6666666666666665D0*
     -       (Cpl92*LInt92*(WProd33 + WProd34)*WProd37) - 
     -      2.6666666666666665D0*
     -       (Cpl92*LInt94*(WProd33 + WProd34)*WProd37))

#ifdef DEBUG
	print *, 'amp(1,1,4) =', amp(1,1,4)
#endif

      endif

      end


