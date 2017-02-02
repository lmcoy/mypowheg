      subroutine realQED21(amp,resnr)
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
*     procnr=1:   {~up  ,  up  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=2:   {~up  , ~do  ,  ve  , ~el  , ~up  , ~up  }
*     procnr=3:   { up  ,  up  ,  ve  , ~el  ,  up  ,  do  }
*     procnr=4:   { up  , ~ch  ,  ve  , ~el  , ~ch  ,  do  }
*     procnr=5:   { up  ,  ch  ,  ve  , ~el  ,  ch  ,  do  }
*     procnr=6:   { up  , ~do  ,  ve  , ~el  , ~up  ,  up  }
*     procnr=7:   { up  , ~do  ,  ve  , ~el  , ~ch  ,  ch  }
*     procnr=8:   { up  , ~do  ,  ve  , ~el  , ~do  ,  do  }
*     procnr=9:   { up  , ~do  ,  ve  , ~el  , ~st  ,  st  }
*     procnr=10:   { up  , ~do  ,  ve  , ~el  , ~bo  ,  bo  }
*     procnr=11:   { up  , ~do  ,  ve  , ~el  ,  gl  ,  ga  }
*     procnr=12:   { up  , ~do  ,  ve  , ~el  ,  gl  ,  gl  }
*     procnr=13:   { up  ,  do  ,  ve  , ~el  ,  do  ,  do  }
*     procnr=14:   { up  , ~st  ,  ve  , ~el  ,  do  , ~st  }
*     procnr=15:   { up  ,  st  ,  ve  , ~el  ,  do  ,  st  }
*     procnr=16:   { up  ,  ga  ,  ve  , ~el  ,  do  ,  gl  }
*     procnr=17:   { up  ,  gl  ,  ve  , ~el  ,  do  ,  ga  }
*     procnr=18:   { up  ,  gl  ,  ve  , ~el  ,  do  ,  gl  }
*     procnr=19:   {~ch  ,  ch  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=20:   {~ch  , ~do  ,  ve  , ~el  , ~up  , ~ch  }
*     procnr=21:   { ch  , ~do  ,  ve  , ~el  , ~up  ,  ch  }
*     procnr=22:   {~do  , ~do  ,  ve  , ~el  , ~up  , ~do  }
*     procnr=23:   {~do  ,  do  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=24:   {~do  , ~st  ,  ve  , ~el  , ~up  , ~st  }
*     procnr=25:   {~do  ,  st  ,  ve  , ~el  , ~up  ,  st  }
*     procnr=26:   {~do  ,  ga  ,  ve  , ~el  , ~up  ,  gl  }
*     procnr=27:   {~do  ,  gl  ,  ve  , ~el  , ~up  ,  ga  }
*     procnr=28:   {~do  ,  gl  ,  ve  , ~el  , ~up  ,  gl  }
*     procnr=29:   {~st  ,  st  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=30:   { ga  ,  gl  ,  ve  , ~el  , ~up  ,  do  }
*     procnr=31:   { gl  ,  gl  ,  ve  , ~el  , ~up  ,  do  }
*
*   Helicities:
*     helnr=1:   {-1, -1, -1, 1, -1, -1}
*     helnr=2:   {-1, 1, -1, 1, 1, -1}
*     helnr=3:   {1, -1, -1, 1, 1, -1}
*     helnr=4:   {-1, 1, -1, 1, -1, 1}
*     helnr=5:   {1, 1, -1, 1, 1, 1}
*     helnr=6:   {-1, 1, -1, 1, -1, -1}
*     helnr=7:   {-1, 1, -1, 1, 1, 1}
*     helnr=8:   {-1, -1, -1, 1, -1, 1}
*     helnr=9:   {1, -1, -1, 1, 1, 1}
*     helnr=10:   {1, 1, -1, 1, 1, -1}
*     helnr=11:   {-1, -1, -1, 1, 1, -1}
*
*   Color Structures:
*     colnr=1:   SUNT[Col1, Col3] SUNT[Col2, Col4]
*     colnr=2:   SUNT[Col1, Col2] SUNT[Col3, Col4]
*     colnr=3:   SUNT[Col1, Col4] SUNT[Col2, Col3]
*     colnr=4:   SUNT[Glu1, Col2, Col1]
*     colnr=5:   SUNT[Glu1, Glu2, Col2, Col1]
*     colnr=6:   SUNT[Glu2, Glu1, Col2, Col1]
*     colnr=7:   SUNT[Glu1, Col1, Col2]
*     colnr=8:   SUNT[Glu1, Glu2, Col1, Col2]
*     colnr=9:   SUNT[Glu2, Glu1, Col1, Col2]


      if(resnr.eq.0)then

        amp(21,3,2) = Prp191*Prp192*
     -   (-(5.333333333333333D0*
     -        (Alfas*Cpl122*Prp193*
     -          (WProd1742 - WProd1743 - WProd1744 + WProd1745 - 
     -            WProd1746))) - 
     -     5.333333333333333D0*
     -      (Alfas*Cpl122*Prp194*
     -        (WProd1742 - 2*WProd1743 - WProd1746 + WProd1749 + 
     -          WProd1750 - WProd1751)))

#ifdef DEBUG
	print *, 'amp(21,3,2) =', amp(21,3,2)
#endif

        amp(21,3,5) = Prp191*Prp192*
     -   (-(5.333333333333333D0*
     -        (Alfas*Cpl122*Prp193*
     -          (WProd1754 + WProd1755 + WProd1756 - WProd1757 + 
     -            WProd1758))) - 
     -     5.333333333333333D0*
     -      (Alfas*Cpl122*Prp194*
     -        (WProd1754 + WProd1755 + WProd1756 - WProd1757 - 
     -          WProd1759 + WProd1760 - WProd1761 + WProd1762 + 
     -          WProd1763)))

#ifdef DEBUG
	print *, 'amp(21,3,5) =', amp(21,3,5)
#endif

        amp(21,2,2) = Prp191*Prp192*
     -   (16*Alfas*Cpl122*Prp193*
     -      (WProd1742 - WProd1743 - WProd1744 + WProd1745 - 
     -        WProd1746) + 
     -     16*Alfas*Cpl122*Prp194*
     -      (WProd1742 - 2*WProd1743 - WProd1746 + WProd1749 + 
     -        WProd1750 - WProd1751))

#ifdef DEBUG
	print *, 'amp(21,2,2) =', amp(21,2,2)
#endif

        amp(21,2,5) = Prp191*Prp192*
     -   (16*Alfas*Cpl122*Prp193*
     -      (WProd1754 + WProd1755 + WProd1756 - WProd1757 + 
     -        WProd1758) + 
     -     16*Alfas*Cpl122*Prp194*
     -      (WProd1754 + WProd1755 + WProd1756 - WProd1757 - 
     -        WProd1759 + WProd1760 - WProd1761 + WProd1762 + 
     -        WProd1763))

#ifdef DEBUG
	print *, 'amp(21,2,5) =', amp(21,2,5)
#endif

      endif

      end


