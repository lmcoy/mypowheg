      subroutine virtEWphotonicFSR2(amp,resnr)
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

        amp(2,1,1) = Prp15*
     -    (0.6666666666666666D0*(Cpl92*WProd38) - 
     -      1.3333333333333333D0*(Cpl92*LInt96*WProd38) + 
     -      2.6666666666666665D0*(Cpl92*LInt98*WProd38) - 
     -      1.3333333333333333D0*
     -       (Cpl92*LInt100*WProd38*
     -         (2*WProd39 - 3*WProd40 - 3*WProd41)) + 
     -      1.3333333333333333D0*
     -       (Cpl92*LInt97*WProd38*(WProd40 + WProd41)) + 
     -      1.3333333333333333D0*
     -       (Cpl92*LInt99*WProd38*(WProd40 + WProd41))) + 
     -   Prp16*(0.6666666666666666D0*(Cpl93*WProd38) - 
     -      1.3333333333333333D0*(Cpl93*LInt96*WProd38) + 
     -      2.6666666666666665D0*(Cpl93*LInt98*WProd38) - 
     -      1.3333333333333333D0*
     -       (Cpl93*LInt100*WProd38*
     -         (2*WProd39 - 3*WProd40 - 3*WProd41)) + 
     -      1.3333333333333333D0*
     -       (Cpl93*LInt97*WProd38*(WProd40 + WProd41)) + 
     -      1.3333333333333333D0*
     -       (Cpl93*LInt99*WProd38*(WProd40 + WProd41)))

#ifdef DEBUG
	print *, 'amp(2,1,1) =', amp(2,1,1)
#endif

        tmp77 = 1.3333333333333333D0*
     -    (Cpl93*(WProd40 + WProd41)*WProd42) - 
     -   0.6666666666666666D0*(Cpl94*(WProd40 + WProd41)*WProd42)
        amp(2,1,2) = Prp16*
     -    (LInt97*tmp77 + LInt99*tmp77 + 
     -      0.6666666666666666D0*(Cpl93*WProd42) - 
     -      0.3333333333333333D0*(Cpl94*WProd42) + 
     -      LInt96*(-(1.3333333333333333D0*(Cpl93*WProd42)) + 
     -         0.6666666666666666D0*(Cpl94*WProd42)) + 
     -      LInt98*(2.6666666666666665D0*(Cpl93*WProd42) - 
     -         1.3333333333333333D0*(Cpl94*WProd42)) + 
     -      LInt100*(-(1.3333333333333333D0*
     -            (Cpl93*(2*WProd39 - 3*WProd40 - 3*WProd41)*
     -              WProd42)) + 
     -         0.6666666666666666D0*
     -          (Cpl94*(2*WProd39 - 3*WProd40 - 3*WProd41)*WProd42)
     -         )) + Prp15*(0.6666666666666666D0*(Cpl92*WProd42) - 
     -      1.3333333333333333D0*(Cpl92*LInt96*WProd42) + 
     -      2.6666666666666665D0*(Cpl92*LInt98*WProd42) - 
     -      1.3333333333333333D0*
     -       (Cpl92*LInt100*(2*WProd39 - 3*WProd40 - 3*WProd41)*
     -         WProd42) + 1.3333333333333333D0*
     -       (Cpl92*LInt97*(WProd40 + WProd41)*WProd42) + 
     -      1.3333333333333333D0*
     -       (Cpl92*LInt99*(WProd40 + WProd41)*WProd42))

#ifdef DEBUG
	print *, 'amp(2,1,2) =', amp(2,1,2)
#endif

        tmp78 = -2*Cpl94*(WProd40 + WProd41)*WProd43 + 
     -   1.3333333333333333D0*(Cpl93*(WProd40 + WProd41)*WProd43)
        amp(2,1,3) = Prp16*
     -    (LInt97*tmp78 + LInt99*tmp78 - Cpl94*WProd43 + 
     -      0.6666666666666666D0*(Cpl93*WProd43) + 
     -      LInt96*(2*Cpl94*WProd43 - 
     -         1.3333333333333333D0*(Cpl93*WProd43)) + 
     -      LInt98*(-4*Cpl94*WProd43 + 
     -         2.6666666666666665D0*(Cpl93*WProd43)) + 
     -      LInt100*(2*Cpl94*(2*WProd39 - 3*WProd40 - 3*WProd41)*
     -          WProd43 - 1.3333333333333333D0*
     -          (Cpl93*(2*WProd39 - 3*WProd40 - 3*WProd41)*WProd43)
     -         )) + Prp15*(0.6666666666666666D0*(Cpl92*WProd43) - 
     -      1.3333333333333333D0*(Cpl92*LInt96*WProd43) + 
     -      2.6666666666666665D0*(Cpl92*LInt98*WProd43) - 
     -      1.3333333333333333D0*
     -       (Cpl92*LInt100*(2*WProd39 - 3*WProd40 - 3*WProd41)*
     -         WProd43) + 1.3333333333333333D0*
     -       (Cpl92*LInt97*(WProd40 + WProd41)*WProd43) + 
     -      1.3333333333333333D0*
     -       (Cpl92*LInt99*(WProd40 + WProd41)*WProd43))

#ifdef DEBUG
	print *, 'amp(2,1,3) =', amp(2,1,3)
#endif

        tmp79 = Cpl95*(WProd40 + WProd41)*WProd44 + 
     -   1.3333333333333333D0*(Cpl93*(WProd40 + WProd41)*WProd44) - 
     -   2.6666666666666665D0*(Cpl94*(WProd40 + WProd41)*WProd44)
        amp(2,1,4) = Prp16*
     -    (LInt97*tmp79 + LInt99*tmp79 + 
     -      0.6666666666666666D0*(Cpl93*WProd44) - 
     -      1.3333333333333333D0*(Cpl94*WProd44) + 
     -      LInt96*(-(Cpl95*WProd44) - 
     -         1.3333333333333333D0*(Cpl93*WProd44) + 
     -         2.6666666666666665D0*(Cpl94*WProd44)) + 
     -      LInt98*(2*Cpl95*WProd44 + 
     -         2.6666666666666665D0*(Cpl93*WProd44) - 
     -         5.333333333333333D0*(Cpl94*WProd44)) + 
     -      0.5D0*(Cpl95*WProd44) + 
     -      LInt100*(-(Cpl95*(2*WProd39 - 3*WProd40 - 3*WProd41)*
     -            WProd44) - 
     -         1.3333333333333333D0*
     -          (Cpl93*(2*WProd39 - 3*WProd40 - 3*WProd41)*WProd44)
     -           + 2.6666666666666665D0*
     -          (Cpl94*(2*WProd39 - 3*WProd40 - 3*WProd41)*WProd44)
     -         )) + Prp15*(0.6666666666666666D0*(Cpl92*WProd44) - 
     -      1.3333333333333333D0*(Cpl92*LInt96*WProd44) + 
     -      2.6666666666666665D0*(Cpl92*LInt98*WProd44) - 
     -      1.3333333333333333D0*
     -       (Cpl92*LInt100*(2*WProd39 - 3*WProd40 - 3*WProd41)*
     -         WProd44) + 1.3333333333333333D0*
     -       (Cpl92*LInt97*(WProd40 + WProd41)*WProd44) + 
     -      1.3333333333333333D0*
     -       (Cpl92*LInt99*(WProd40 + WProd41)*WProd44))

#ifdef DEBUG
	print *, 'amp(2,1,4) =', amp(2,1,4)
#endif

      endif

      end


