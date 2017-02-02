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
*     procnr=1:   { up  , ~do  ,  ve  , ~el  ,  gl  }
*     procnr=2:   { up  ,  ga  ,  ve  , ~el  ,  do  }
*     procnr=3:   { up  ,  gl  ,  ve  , ~el  ,  do  }
*     procnr=4:   {~do  ,  ga  ,  ve  , ~el  , ~up  }
*     procnr=5:   {~do  ,  gl  ,  ve  , ~el  , ~up  }
*
*   Helicities:
*     helnr=1:   {-1, 1, -1, 1, -1}
*     helnr=2:   {-1, 1, -1, 1, 1}
*     helnr=3:   {-1, -1, -1, 1, -1}
*     helnr=4:   {1, -1, -1, 1, 1}
*     helnr=5:   {1, 1, -1, 1, 1}
*
*   Color Structures:
*     colnr=1:   SUNT[Glu1, Col2, Col1]
*     colnr=2:   SUNT[Col1, Col2]
*     colnr=3:   SUNT[Glu1, Col1, Col2]


      if(resnr.eq.0)then

        amp(2,2,3) = -(3.7712361663282534D0*
     -      (Cpl2*Prp4*Prp5*(WProd22*WProd23 - WProd24 + WProd25))/
     -       WProd21) - 5.656854249492381D0*
     -    (Cpl2*Prp6*Prp7*(WProd25 + WProd23*WProd28 - WProd29 - 
     -         WProd30))/WProd21 + 
     -   Prp5*(1.8856180831641267D0*
     -       (Cpl2*Prp8*WProd23*(WProd22 - WProd28 - WProd33))/
     -        WProd21 - 5.656854249492381D0*
     -       (Cpl2*Prp7*(-WProd24 + WProd25 + WProd23*WProd28 + 
     -            WProd23*WProd33))/WProd21)

#ifdef DEBUG
	print *, 'amp(2,2,3) =', amp(2,2,3)
#endif

        amp(2,2,1) = -(3.7712361663282534D0*
     -      (Cpl2*Prp4*Prp5*
     -         (WProd23*WProd36 - WProd37 + WProd38 - WProd39))/
     -       WProd35) - 5.656854249492381D0*
     -    (Cpl2*Prp6*Prp7*(WProd38 + WProd40 + WProd23*WProd41 - 
     -         WProd42 - WProd43))/WProd35 + 
     -   Prp5*(1.8856180831641267D0*
     -       (Cpl2*Prp8*(WProd23*WProd36 - WProd37 - 
     -            WProd23*WProd41 - WProd23*WProd44))/WProd35 - 
     -      5.656854249492381D0*
     -       (Cpl2*Prp7*(WProd38 - WProd39 + WProd23*WProd41 + 
     -            WProd23*WProd44))/WProd35)

#ifdef DEBUG
	print *, 'amp(2,2,1) =', amp(2,2,1)
#endif

      endif

      end


