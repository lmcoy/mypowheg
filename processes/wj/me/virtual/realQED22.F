      subroutine realQED22(amp,resnr)
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

        amp(22,3,2) = Prp195*Prp196*
     -   (5.333333333333333D0*
     -      (Alfas*Cpl122*Prp197*
     -        (WProd1766 + WProd1767 + WProd1768 + WProd1769 - 
     -          WProd1770)) + 
     -     5.333333333333333D0*
     -      (Alfas*Cpl122*Prp198*
     -        (WProd1766 + WProd1767 + WProd1768 - WProd1770 + 
     -          WProd1773 + WProd1774 - WProd1775 + WProd1776 - 
     -          WProd1777)))

#ifdef DEBUG
	print *, 'amp(22,3,2) =', amp(22,3,2)
#endif

        amp(22,3,3) = Prp195*
     -   (16*Alfas*Cpl122*Prp199*Prp200*
     -      (WProd1780 + WProd1781 + WProd1782 - WProd1783 - 
     -        WProd1784 - WProd1785) + 
     -     16*Alfas*Cpl122*Prp197*Prp199*
     -      (WProd1782 - WProd1785 + WProd1788 + WProd1789 + 
     -        WProd1790))

#ifdef DEBUG
	print *, 'amp(22,3,3) =', amp(22,3,3)
#endif

        amp(22,3,5) = Prp195*
     -   (16*Alfas*Cpl122*Prp199*Prp200*
     -      (WProd1791 - WProd1792 - WProd1793 + WProd1794 + 
     -        WProd1795 - WProd1796) + 
     -     16*Alfas*Cpl122*Prp197*Prp199*
     -      (WProd1791 - WProd1792 - WProd1796 - WProd1797 + 
     -        WProd1798) + 
     -     Prp196*(5.333333333333333D0*
     -         (Alfas*Cpl122*Prp197*
     -           (-WProd1792 - WProd1793 + WProd1799 - WProd1800 + 
     -             WProd1801)) + 
     -        5.333333333333333D0*
     -         (Alfas*Cpl122*Prp198*
     -           (-2*WProd1792 - WProd1797 + WProd1799 - 
     -             WProd1800 + WProd1802 + WProd1803))))

#ifdef DEBUG
	print *, 'amp(22,3,5) =', amp(22,3,5)
#endif

        amp(22,1,2) = Prp195*Prp196*
     -   (-16*Alfas*Cpl122*Prp197*
     -      (WProd1766 + WProd1767 + WProd1768 + WProd1769 - 
     -        WProd1770) - 
     -     16*Alfas*Cpl122*Prp198*
     -      (WProd1766 + WProd1767 + WProd1768 - WProd1770 + 
     -        WProd1773 + WProd1774 - WProd1775 + WProd1776 - 
     -        WProd1777))

#ifdef DEBUG
	print *, 'amp(22,1,2) =', amp(22,1,2)
#endif

        amp(22,1,3) = Prp195*
     -   (-(5.333333333333333D0*
     -        (Alfas*Cpl122*Prp199*Prp200*
     -          (WProd1780 + WProd1781 + WProd1782 - WProd1783 - 
     -            WProd1784 - WProd1785))) - 
     -     5.333333333333333D0*
     -      (Alfas*Cpl122*Prp197*Prp199*
     -        (WProd1782 - WProd1785 + WProd1788 + WProd1789 + 
     -          WProd1790)))

#ifdef DEBUG
	print *, 'amp(22,1,3) =', amp(22,1,3)
#endif

        amp(22,1,5) = Prp195*
     -   (Prp196*(-16*Alfas*Cpl122*Prp197*
     -         (-WProd1792 - WProd1793 + WProd1799 - WProd1800 + 
     -           WProd1801) - 
     -        16*Alfas*Cpl122*Prp198*
     -         (-2*WProd1792 - WProd1797 + WProd1799 - WProd1800 + 
     -           WProd1802 + WProd1803)) - 
     -     5.333333333333333D0*
     -      (Alfas*Cpl122*Prp199*Prp200*
     -        (WProd1791 - WProd1792 - WProd1793 + WProd1794 + 
     -          WProd1795 - WProd1796)) - 
     -     5.333333333333333D0*
     -      (Alfas*Cpl122*Prp197*Prp199*
     -        (WProd1791 - WProd1792 - WProd1796 - WProd1797 + 
     -          WProd1798)))

#ifdef DEBUG
	print *, 'amp(22,1,5) =', amp(22,1,5)
#endif

      endif

      end


