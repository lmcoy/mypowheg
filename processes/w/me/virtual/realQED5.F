      subroutine realQED5(amp,resnr)
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

        amp(5,1,4) = 1.8856180831641267D0*
     -    (Cpl109*Prp22*Prp23*WProd103*WProd104)/WProd102 + 
     -   5.656854249492381D0*
     -    (Cpl109*Prp24*Prp25*
     -       (WProd107 - WProd108 + WProd104*WProd109 - WProd110))/
     -     WProd102 + Prp23*
     -    (5.656854249492381D0*
     -       (Cpl109*Prp25*
     -          (WProd107 + WProd104*WProd109 - WProd113 + 
     -            WProd104*WProd114))/WProd102 + 
     -      3.7712361663282534D0*
     -       (Cpl109*Prp26*
     -          (-(WProd103*WProd104) + WProd107 + 
     -            WProd104*WProd109 - WProd113 + WProd104*WProd114)
     -          )/WProd102)

#ifdef DEBUG
	print *, 'amp(5,1,4) =', amp(5,1,4)
#endif

        amp(5,1,5) = 1.8856180831641267D0*
     -    (Cpl109*Prp22*Prp23*(WProd104*WProd117 + WProd118))/
     -     WProd116 + 5.656854249492381D0*
     -    (Cpl109*Prp24*Prp25*
     -       (WProd119 - WProd120 + WProd104*WProd121 + WProd122 - 
     -         WProd123))/WProd116 + 
     -   Prp23*(-(5.656854249492381D0*
     -         (Cpl109*Prp25*
     -            (-WProd119 - WProd104*WProd121 + WProd124 - 
     -              WProd104*WProd125))/WProd116) - 
     -      3.7712361663282534D0*
     -       (Cpl109*Prp26*
     -          (WProd104*WProd117 + WProd118 - WProd119 - 
     -            WProd104*WProd121 + WProd124 - WProd104*WProd125)
     -          )/WProd116)

#ifdef DEBUG
	print *, 'amp(5,1,5) =', amp(5,1,5)
#endif

      endif

      end


