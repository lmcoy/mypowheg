      subroutine ampsfacctQCD3(amp,resnr)
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

        amp(3,1,3) = 5.656854249492381D0*
     -    (Alfas*Cpl1*GS*Prp65*Prp66**2*
     -       (WProd549*WProd550 - WProd551 + WProd552)*
     -       (WProd553 + WProd554 + WProd555 + WProd556 - 
     -         WProd557)*dZfL1dimreg(3,1,1))/WProd548 - 
     -   2.8284271247461903D0*
     -    (Cpl1*GS*Prp65*Prp66*
     -       (WProd549*WProd550 - WProd551 + WProd552)*
     -       (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -         3*Alfas*dZfL1dimreg(3,1,1) + 
     -         Alfas*dZfL1dimreg(4,1,1)))/WProd548 + 
     -   Prp65*(-(5.656854249492381D0*
     -         (Alfas*Cpl1*GS*Prp67**2*WProd550*
     -            (-WProd554 - WProd556 + WProd558)*
     -            (WProd549 - WProd559 - WProd560)*
     -            dZfL1dimreg(4,1,1))/WProd548) - 
     -      2.8284271247461903D0*
     -       (Cpl1*GS*Prp67*WProd550*
     -          (WProd549 - WProd559 - WProd560)*
     -          (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -            Alfas*dZfL1dimreg(3,1,1) + 
     -            3*Alfas*dZfL1dimreg(4,1,1)))/WProd548)

#ifdef DEBUG
	print *, 'amp(3,1,3) =', amp(3,1,3)
#endif

	tmp13 = WProd550*WProd563 - WProd564 + WProd565 - WProd566
        tmp14 = WProd550*WProd563 - WProd564 - 
     -   WProd550*WProd567 - WProd550*WProd568
        amp(3,1,1) = 5.656854249492381D0*
     -    (Alfas*Cpl1*GS*Prp65*Prp66**2*tmp13*
     -       (WProd553 + WProd554 + WProd555 + WProd556 - 
     -         WProd557)*dZfL1dimreg(3,1,1))/WProd562 - 
     -   2.8284271247461903D0*
     -    (Cpl1*GS*Prp65*Prp66*tmp13*
     -       (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -         3*Alfas*dZfL1dimreg(3,1,1) + 
     -         Alfas*dZfL1dimreg(4,1,1)))/WProd562 + 
     -   Prp65*(-(5.656854249492381D0*
     -         (Alfas*Cpl1*GS*Prp67**2*tmp14*
     -            (-WProd554 - WProd556 + WProd558)*
     -            dZfL1dimreg(4,1,1))/WProd562) - 
     -      2.8284271247461903D0*
     -       (Cpl1*GS*Prp67*tmp14*
     -          (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -            Alfas*dZfL1dimreg(3,1,1) + 
     -            3*Alfas*dZfL1dimreg(4,1,1)))/WProd562)

#ifdef DEBUG
	print *, 'amp(3,1,1) =', amp(3,1,1)
#endif

      endif

      end


