      subroutine realQED3(amp,resnr)
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

        amp(3,1,3) = -(3.7712361663282534D0*
     -      (Cpl109*Prp14*Prp15*
     -         (WProd62*WProd63 - WProd64 + WProd65))/WProd61) - 
     -   5.656854249492381D0*
     -    (Cpl109*Prp16*Prp17*
     -       (WProd65 + WProd63*WProd68 - WProd69 - WProd70))/
     -     WProd61 + Prp15*
     -    (1.8856180831641267D0*
     -       (Cpl109*Prp18*WProd63*(WProd62 - WProd68 - WProd73))/
     -        WProd61 - 5.656854249492381D0*
     -       (Cpl109*Prp17*
     -          (-WProd64 + WProd65 + WProd63*WProd68 + 
     -            WProd63*WProd73))/WProd61)

#ifdef DEBUG
	print *, 'amp(3,1,3) =', amp(3,1,3)
#endif

        amp(3,1,1) = -(3.7712361663282534D0*
     -      (Cpl109*Prp14*Prp15*
     -         (WProd63*WProd76 - WProd77 + WProd78 - WProd79))/
     -       WProd75) - 5.656854249492381D0*
     -    (Cpl109*Prp16*Prp17*
     -       (WProd78 + WProd80 + WProd63*WProd81 - WProd82 - 
     -         WProd83))/WProd75 + 
     -   Prp15*(1.8856180831641267D0*
     -       (Cpl109*Prp18*
     -          (WProd63*WProd76 - WProd77 - WProd63*WProd81 - 
     -            WProd63*WProd84))/WProd75 - 
     -      5.656854249492381D0*
     -       (Cpl109*Prp17*
     -          (WProd78 - WProd79 + WProd63*WProd81 + 
     -            WProd63*WProd84))/WProd75)

#ifdef DEBUG
	print *, 'amp(3,1,1) =', amp(3,1,1)
#endif

      endif

      end


