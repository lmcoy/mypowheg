      subroutine ampsfacctQCD1(amp,resnr)
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

        tmp7 = WProd501*WProd502 - WProd503 + 
     -   WProd501*WProd504 - WProd505 - WProd506 - WProd507
        amp(1,1,1) = Prp57*
     -   (-(5.656854249492381D0*
     -        (Alfas*Cpl1*GS*Prp58**2*tmp7*
     -           (-WProd498 + WProd499 + WProd500)*
     -           dZfL1dimreg(3,1,1))/WProd508) - 
     -     5.656854249492381D0*
     -      (Alfas*Cpl1*GS*Prp59**2*
     -         (-WProd499 - WProd500 + WProd510)*
     -         (-WProd506 + WProd511)*dZfL1dimreg(4,1,1))/WProd508-
     -       2.8284271247461903D0*
     -      (Cpl1*GS*Prp58*tmp7*
     -         (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -           3*Alfas*dZfL1dimreg(3,1,1) + 
     -           Alfas*dZfL1dimreg(4,1,1)))/WProd508 - 
     -     2.8284271247461903D0*
     -      (Cpl1*GS*Prp59*(-WProd506 + WProd511)*
     -         (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -           Alfas*dZfL1dimreg(3,1,1) + 
     -           3*Alfas*dZfL1dimreg(4,1,1)))/WProd508)

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        tmp8 = WProd501*WProd514 - WProd515 + 
     -   WProd501*WProd516 - WProd517 - WProd518 - WProd519
        amp(1,1,2) = Prp57*
     -   (-(5.656854249492381D0*
     -        (Alfas*Cpl1*GS*Prp58**2*tmp8*
     -           (-WProd498 + WProd499 + WProd500)*
     -           dZfL1dimreg(3,1,1))/WProd513) - 
     -     5.656854249492381D0*
     -      (Alfas*Cpl1*GS*Prp59**2*
     -         (-WProd499 - WProd500 + WProd510)*
     -         (-WProd518 + WProd520)*dZfL1dimreg(4,1,1))/WProd513-
     -       2.8284271247461903D0*
     -      (Cpl1*GS*Prp58*tmp8*
     -         (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -           3*Alfas*dZfL1dimreg(3,1,1) + 
     -           Alfas*dZfL1dimreg(4,1,1)))/WProd513 - 
     -     2.8284271247461903D0*
     -      (Cpl1*GS*Prp59*(-WProd518 + WProd520)*
     -         (Alfas*dZGG1dimreg + 2*Alfas*mydZg1dimreg + 
     -           Alfas*dZfL1dimreg(3,1,1) + 
     -           3*Alfas*dZfL1dimreg(4,1,1)))/WProd513)

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

      endif

      end


