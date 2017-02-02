      subroutine ampsfacctQCD4(amp,resnr)
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

        tmp15 = -(WProd570*WProd576) + WProd578 + 
     -   WProd576*WProd580 - WProd583 + WProd576*WProd584
	tmp16 = Alfas*dZfL1dimreg(3,1,1) + Alfas*dZfL1dimreg(4,1,1)
        amp(4,2,4) = 2.8284271247461903D0*
     -    (Cpl2*Prp70*Prp71*tmp16*
     -       (WProd578 - WProd579 + WProd576*WProd580 - WProd581))/
     -     WProd569 - 1.8856180831641267D0*
     -    (Alfas*Cpl2*Prp68*Prp69**2*WProd570*
     -       (WProd571 + WProd572 + WProd573 + WProd574 - 
     -         WProd575)*WProd576*dZfL1dimreg(4,1,1))/WProd569 + 
     -   Prp68*(2.8284271247461903D0*
     -       (Cpl2*Prp71*tmp16*
     -          (WProd578 + WProd576*WProd580 - WProd583 + 
     -            WProd576*WProd584))/WProd569 + 
     -      3.7712361663282534D0*
     -       (Alfas*Cpl2*Prp72**2*tmp15*
     -          (-WProd572 - WProd574 + WProd577)*
     -          dZfL1dimreg(3,1,1))/WProd569 + 
     -      1.8856180831641267D0*
     -       (Cpl2*Prp72*tmp15*
     -          (3*Alfas*dZfL1dimreg(3,1,1) + 
     -            Alfas*dZfL1dimreg(4,1,1)))/WProd569) + 
     -   0.9428090415820634D0*
     -    (Cpl2*Prp68*Prp69*WProd570*WProd576*
     -       (Alfas*dZfL1dimreg(3,1,1) + 
     -         3*Alfas*dZfL1dimreg(4,1,1)))/WProd569

#ifdef DEBUG
	print *, 'amp(4,2,4) =', amp(4,2,4)
#endif

        tmp17 = WProd576*WProd587 + WProd588 - WProd589 - 
     -   WProd576*WProd591 + WProd594 - WProd576*WProd595
	tmp18 = Alfas*dZfL1dimreg(3,1,1) + Alfas*dZfL1dimreg(4,1,1)
        amp(4,2,5) = 2.8284271247461903D0*
     -    (Cpl2*Prp70*Prp71*tmp18*
     -       (WProd589 - WProd590 + WProd576*WProd591 + WProd592 - 
     -         WProd593))/WProd586 - 
     -   1.8856180831641267D0*
     -    (Alfas*Cpl2*Prp68*Prp69**2*
     -       (WProd571 + WProd572 + WProd573 + WProd574 - 
     -         WProd575)*(WProd576*WProd587 + WProd588)*
     -       dZfL1dimreg(4,1,1))/WProd586 + 
     -   Prp68*(-(2.8284271247461903D0*
     -         (Cpl2*Prp71*tmp18*
     -            (-WProd589 - WProd576*WProd591 + WProd594 - 
     -              WProd576*WProd595))/WProd586) - 
     -      3.7712361663282534D0*
     -       (Alfas*Cpl2*Prp72**2*tmp17*
     -          (-WProd572 - WProd574 + WProd577)*
     -          dZfL1dimreg(3,1,1))/WProd586 - 
     -      1.8856180831641267D0*
     -       (Cpl2*Prp72*tmp17*
     -          (3*Alfas*dZfL1dimreg(3,1,1) + 
     -            Alfas*dZfL1dimreg(4,1,1)))/WProd586) + 
     -   0.9428090415820634D0*
     -    (Cpl2*Prp68*Prp69*(WProd576*WProd587 + WProd588)*
     -       (Alfas*dZfL1dimreg(3,1,1) + 
     -         3*Alfas*dZfL1dimreg(4,1,1)))/WProd586

#ifdef DEBUG
	print *, 'amp(4,2,5) =', amp(4,2,5)
#endif

      endif

      end


