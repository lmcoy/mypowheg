      subroutine realQED1(amp,resnr)
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
*     procnr=1:   { up  , ~do  ,  ve  , ~el  ,  ga  }
*     procnr=2:   { up  , ~do  ,  ve  , ~el  ,  gl  }
*     procnr=3:   { up  ,  ga  ,  ve  , ~el  ,  do  }
*     procnr=4:   { up  ,  gl  ,  ve  , ~el  ,  do  }
*     procnr=5:   {~do  ,  ga  ,  ve  , ~el  , ~up  }
*     procnr=6:   {~do  ,  gl  ,  ve  , ~el  , ~up  }
*
*   Helicities:
*     helnr=1:   {-1, 1, -1, 1, -1}
*     helnr=2:   {-1, 1, -1, 1, 1}
*     helnr=3:   {-1, -1, -1, 1, -1}
*     helnr=4:   {1, -1, -1, 1, 1}
*     helnr=5:   {1, 1, -1, 1, 1}
*
*   Color Structures:
*     colnr=1:   SUNT[Col1, Col2]
*     colnr=2:   SUNT[Glu1, Col2, Col1]
*     colnr=3:   SUNT[Glu1, Col1, Col2]


      if(resnr.eq.0)then

        amp(1,1,1) = Prp6*
     -    (-(3.7712361663282534D0*
     -         (Cpl109*Prp7*
     -            (WProd18*WProd19 - WProd20 + WProd18*WProd21 - 
     -              WProd22 - WProd23 - WProd24))/WProd25) + 
     -      1.8856180831641267D0*
     -       (Cpl109*Prp8*(-WProd23 + WProd27))/WProd25) + 
     -   Prp9*(-(5.656854249492381D0*
     -         (Cpl109*Prp6*
     -            (WProd18*WProd19 - WProd20 + WProd18*WProd21 - 
     -              WProd22 - WProd24 - WProd27))/WProd25) - 
     -      5.656854249492381D0*
     -       (Cpl109*Prp10*
     -          (WProd18*WProd19 - WProd20 + WProd18*WProd21 - 
     -            WProd22 - WProd27 - WProd30))/WProd25)

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        amp(1,1,2) = Prp6*
     -    (-(3.7712361663282534D0*
     -         (Cpl109*Prp7*
     -            (WProd18*WProd33 - WProd34 + WProd18*WProd35 - 
     -              WProd36 - WProd37 - WProd38))/WProd32) + 
     -      1.8856180831641267D0*
     -       (Cpl109*Prp8*(-WProd37 + WProd39))/WProd32) + 
     -   Prp9*(-(5.656854249492381D0*
     -         (Cpl109*Prp6*
     -            (WProd18*WProd33 - WProd34 + WProd18*WProd35 - 
     -              WProd36 - WProd38 - WProd39))/WProd32) - 
     -      5.656854249492381D0*
     -       (Cpl109*Prp10*
     -          (WProd18*WProd33 - WProd34 + WProd18*WProd35 - 
     -            WProd36 - WProd39 - WProd40))/WProd32)

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

      endif

      end


