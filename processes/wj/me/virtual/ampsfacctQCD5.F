      subroutine ampsfacctQCD5(amp,resnr)
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "params.h" 

      integer helnr,colnr,resnr,procnr

      double complex amp(nproc,ncol,nhel)

      double precision Alfas,Alfas2,GS

      double precision nflavour,mureg,Qren2

      double precision mydZg1dimreg

      common /my_alfas_common/Alfas,Alfas2,GS
      common /ct_common/nflavour,mureg,Qren2 

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


      mydZg1dimreg=dZg1dimreg+1D0/4D0/3.1415926535897932384626433832795029D0*
     &                         (nflavour/3D0 - 5.5D0)*(log(mureg)-log(Qren2))

      if(resnr.eq.0)then

        tmp19 = -(WProd597*WProd603) + WProd605 - WProd606 + 
     -   WProd603*WProd607 + WProd603*WProd608
        amp(5,3,4) = 5.656854249492381D0*
     -    (Alfas*Cpl1*GS*Prp73*Prp74**2*WProd597*
     -       (WProd598 + WProd599 + WProd600 + WProd601 - 
     -         WProd602)*WProd603*dZfL1dimreg(4,1,1))/WProd596 + 
     -   Prp73*(5.656854249492381D0*
     -       (Alfas*Cpl1*GS*Prp75**2*tmp19*
     -          (-WProd599 - WProd601 + WProd604)*
     -          dZfL1dimreg(3,1,1))/WProd596 + 
     -      2.8284271247461903D0*
     -       (Cpl1*GS*Prp75*tmp19*
     -          (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -            3*Alfas*dZfL1dimreg(3,1,1) + 
     -            Alfas*dZfL1dimreg(4,1,1)))/WProd596) - 
     -   2.8284271247461903D0*
     -    (Cpl1*GS*Prp73*Prp74*WProd597*WProd603*
     -       (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -         Alfas*dZfL1dimreg(3,1,1) + 
     -         3*Alfas*dZfL1dimreg(4,1,1)))/WProd596

#ifdef DEBUG
	print *, 'amp(5,3,4) =', amp(5,3,4)
#endif

        tmp20 = WProd603*WProd611 + WProd612 + WProd613 - 
     -   WProd614 - WProd603*WProd615 - WProd603*WProd616
        amp(5,3,5) = 5.656854249492381D0*
     -    (Alfas*Cpl1*GS*Prp73*Prp74**2*
     -       (WProd598 + WProd599 + WProd600 + WProd601 - 
     -         WProd602)*(WProd603*WProd611 + WProd612)*
     -       dZfL1dimreg(4,1,1))/WProd610 + 
     -   Prp73*(-(5.656854249492381D0*
     -         (Alfas*Cpl1*GS*Prp75**2*tmp20*
     -            (-WProd599 - WProd601 + WProd604)*
     -            dZfL1dimreg(3,1,1))/WProd610) - 
     -      2.8284271247461903D0*
     -       (Cpl1*GS*Prp75*tmp20*
     -          (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -            3*Alfas*dZfL1dimreg(3,1,1) + 
     -            Alfas*dZfL1dimreg(4,1,1)))/WProd610) - 
     -   2.8284271247461903D0*
     -    (Cpl1*GS*Prp73*Prp74*(WProd603*WProd611 + WProd612)*
     -       (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -         Alfas*dZfL1dimreg(3,1,1) + 
     -         3*Alfas*dZfL1dimreg(4,1,1)))/WProd610

#ifdef DEBUG
	print *, 'amp(5,3,5) =', amp(5,3,5)
#endif

      endif

      end


