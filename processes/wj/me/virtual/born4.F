      subroutine born4(amp,resnr)
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

        amp(4,2,4) = 1.8856180831641267D0*
     -    (Cpl2*Prp12*Prp13*WProd63*WProd64)/WProd62 + 
     -   5.656854249492381D0*
     -    (Cpl2*Prp14*Prp15*
     -       (WProd67 - WProd68 + WProd64*WProd69 - WProd70))/
     -     WProd62 + Prp13*
     -    (5.656854249492381D0*
     -       (Cpl2*Prp15*(WProd67 + WProd64*WProd69 - WProd73 + 
     -            WProd64*WProd74))/WProd62 + 
     -      3.7712361663282534D0*
     -       (Cpl2*Prp16*(-(WProd63*WProd64) + WProd67 + 
     -            WProd64*WProd69 - WProd73 + WProd64*WProd74))/
     -        WProd62)

#ifdef DEBUG
	print *, 'amp(4,2,4) =', amp(4,2,4)
#endif

        amp(4,2,5) = 1.8856180831641267D0*
     -    (Cpl2*Prp12*Prp13*(WProd64*WProd77 + WProd78))/WProd76 + 
     -   5.656854249492381D0*
     -    (Cpl2*Prp14*Prp15*
     -       (WProd79 - WProd80 + WProd64*WProd81 + WProd82 - 
     -         WProd83))/WProd76 + 
     -   Prp13*(-(5.656854249492381D0*
     -         (Cpl2*Prp15*
     -            (-WProd79 - WProd64*WProd81 + WProd84 - 
     -              WProd64*WProd85))/WProd76) - 
     -      3.7712361663282534D0*
     -       (Cpl2*Prp16*(WProd64*WProd77 + WProd78 - WProd79 - 
     -            WProd64*WProd81 + WProd84 - WProd64*WProd85))/
     -        WProd76)

#ifdef DEBUG
	print *, 'amp(4,2,5) =', amp(4,2,5)
#endif

      endif

      end


