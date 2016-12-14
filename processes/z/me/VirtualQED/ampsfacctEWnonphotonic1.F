      subroutine ampsfacctEWnonphotonic1(amp,resnr)
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
*     procnr=1:   {~up  ,  up  , ~el  ,  el  }
*     procnr=2:   {~do  ,  do  , ~el  ,  el  }
*
*   Helicities:
*     helnr=1:   {-1, 1, -1, 1}
*     helnr=2:   {-1, 1, 1, -1}
*     helnr=3:   {1, -1, -1, 1}
*     helnr=4:   {1, -1, 1, -1}
*
*   Color Structures:
*     colnr=1:   SUNT[Col1, Col2]


      if(resnr.eq.0)then

        amp(1,1,1) = -(5.333333333333333D0*
     -      (Cpl1*dZAA1*Prp29**2*(WProd137 + WProd138)*WProd139))+
     -     Prp30**2*(5.333333333333333D0*(Cpl155*dZZZ1*WProd139) + 
     -      5.333333333333333D0*
     -       (Cpl2*(dMZsq1 - dZZZ1*WProd137 - dZZZ1*WProd138)*
     -         WProd139)) + 
     -   Prp29*(5.333333333333333D0*(Cpl156*dZZA1*WProd139) + 
     -      5.333333333333333D0*
     -       (Cpl1*WProd139*
     -         (dZAA1 + 2*dZe1 + dZfR1(2,1,1) + dZfR1(3,1,1)))) + 
     -   Prp30*(10.666666666666666D0*(Cpl157*dSW1*WProd139) + 
     -      5.333333333333333D0*(Cpl156*dZAZ1*WProd139) + 
     -      Prp29*(5.333333333333333D0*(Cpl158*dZZA1*WProd139) - 
     -         5.333333333333333D0*
     -          (Cpl156*(dZAZ1 + dZZA1)*(WProd137 + WProd138)*
     -            WProd139)) + 
     -      5.333333333333333D0*
     -       (Cpl2*WProd139*
     -         (2*dZe1 + dZZZ1 + dZfR1(2,1,1) + dZfR1(3,1,1))))

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

	tmp1 = 2*dZe1 + dZZZ1 + dZfL1(2,1,1) + dZfR1(3,1,1)
        amp(1,1,2) = Prp30*
     -    (5.333333333333333D0*(Cpl157*dSW1*WProd141) - 
     -      2.6666666666666665D0*(Cpl162*dSW1*WProd141) + 
     -      2.6666666666666665D0*(Cpl163*dSW1*WProd141) + 
     -      2.6666666666666665D0*(Cpl164*dSW1*WProd141) + 
     -      2.6666666666666665D0*(Cpl156*dZAZ1*WProd141) - 
     -      1.3333333333333333D0*(Cpl160*dZAZ1*WProd141) + 
     -      2.6666666666666665D0*(Cpl161*dZAZ1*WProd141) + 
     -      5.333333333333333D0*(Cpl2*tmp1*WProd141) - 
     -      2.6666666666666665D0*(Cpl3*tmp1*WProd141) + 
     -      Prp29*(2.6666666666666665D0*(Cpl158*dZZA1*WProd141) - 
     -         1.3333333333333333D0*(Cpl165*dZZA1*WProd141) + 
     -         2.6666666666666665D0*(Cpl166*dZZA1*WProd141) - 
     -         2.6666666666666665D0*
     -          (Cpl156*(dZAZ1 + dZZA1)*(WProd137 + WProd138)*
     -            WProd141) + 
     -         1.3333333333333333D0*
     -          (Cpl160*(dZAZ1 + dZZA1)*(WProd137 + WProd138)*
     -            WProd141) - 
     -         2.6666666666666665D0*
     -          (Cpl161*(dZAZ1 + dZZA1)*(WProd137 + WProd138)*
     -            WProd141))) - 
     -   5.333333333333333D0*
     -    (Cpl1*dZAA1*Prp29**2*(WProd137 + WProd138)*WProd141) + 
     -   Prp30**2*(5.333333333333333D0*(Cpl155*dZZZ1*WProd141) - 
     -      2.6666666666666665D0*(Cpl159*dZZZ1*WProd141) + 
     -      5.333333333333333D0*
     -       (Cpl2*(dMZsq1 - dZZZ1*WProd137 - dZZZ1*WProd138)*
     -         WProd141) - 
     -      2.6666666666666665D0*
     -       (Cpl3*(dMZsq1 - dZZZ1*WProd137 - dZZZ1*WProd138)*
     -         WProd141)) + 
     -   Prp29*(2.6666666666666665D0*(Cpl156*dZZA1*WProd141) - 
     -      1.3333333333333333D0*(Cpl160*dZZA1*WProd141) + 
     -      2.6666666666666665D0*(Cpl161*dZZA1*WProd141) + 
     -      5.333333333333333D0*
     -       (Cpl1*WProd141*
     -         (dZAA1 + 2*dZe1 + dZfL1(2,1,1) + dZfR1(3,1,1))))

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

	tmp2 = 2*dZe1 + dZZZ1 + dZfL1(3,1,1) + dZfR1(2,1,1)
        amp(1,1,3) = Prp30*
     -    (-4*Cpl162*dSW1*WProd142 + 4*Cpl163*dSW1*WProd142 - 
     -      2*Cpl160*dZAZ1*WProd142 - 4*Cpl3*tmp2*WProd142 + 
     -      5.333333333333333D0*(Cpl157*dSW1*WProd142) + 
     -      1.3333333333333333D0*(Cpl164*dSW1*WProd142) + 
     -      2.6666666666666665D0*(Cpl156*dZAZ1*WProd142) + 
     -      2.6666666666666665D0*(Cpl161*dZAZ1*WProd142) + 
     -      5.333333333333333D0*(Cpl2*tmp2*WProd142) + 
     -      Prp29*(-2*Cpl165*dZZA1*WProd142 + 
     -         2*Cpl160*(dZAZ1 + dZZA1)*(WProd137 + WProd138)*
     -          WProd142 + 
     -         2.6666666666666665D0*(Cpl158*dZZA1*WProd142) + 
     -         2.6666666666666665D0*(Cpl166*dZZA1*WProd142) - 
     -         2.6666666666666665D0*
     -          (Cpl156*(dZAZ1 + dZZA1)*(WProd137 + WProd138)*
     -            WProd142) - 
     -         2.6666666666666665D0*
     -          (Cpl161*(dZAZ1 + dZZA1)*(WProd137 + WProd138)*
     -            WProd142))) - 
     -   5.333333333333333D0*
     -    (Cpl1*dZAA1*Prp29**2*(WProd137 + WProd138)*WProd142) + 
     -   Prp30**2*(-4*Cpl159*dZZZ1*WProd142 - 
     -      4*Cpl3*(dMZsq1 - dZZZ1*WProd137 - dZZZ1*WProd138)*
     -       WProd142 + 5.333333333333333D0*
     -       (Cpl155*dZZZ1*WProd142) + 
     -      5.333333333333333D0*
     -       (Cpl2*(dMZsq1 - dZZZ1*WProd137 - dZZZ1*WProd138)*
     -         WProd142)) + 
     -   Prp29*(-2*Cpl160*dZZA1*WProd142 + 
     -      2.6666666666666665D0*(Cpl156*dZZA1*WProd142) + 
     -      2.6666666666666665D0*(Cpl161*dZZA1*WProd142) + 
     -      5.333333333333333D0*
     -       (Cpl1*WProd142*
     -         (dZAA1 + 2*dZe1 + dZfL1(3,1,1) + dZfR1(2,1,1))))

#ifdef DEBUG
	print *, 'amp(1,1,3) =', amp(1,1,3)
#endif

	tmp3 = 2*dZe1 + dZZZ1 + dZfL1(2,1,1) + dZfL1(3,1,1)
        amp(1,1,4) = Prp30*
     -    (4*Cpl164*dSW1*WProd143 - 4*Cpl168*dSW1*WProd143 + 
     -      2*Cpl4*tmp3*WProd143 - 
     -      2.6666666666666665D0*(Cpl162*dSW1*WProd143) + 
     -      6.666666666666667D0*(Cpl163*dSW1*WProd143) - 
     -      3.3333333333333335D0*(Cpl160*dZAZ1*WProd143) + 
     -      5.333333333333333D0*(Cpl161*dZAZ1*WProd143) + 
     -      5.333333333333333D0*(Cpl2*tmp3*WProd143) - 
     -      6.666666666666667D0*(Cpl3*tmp3*WProd143) + 
     -      Prp29*(-(3.3333333333333335D0*(Cpl165*dZZA1*WProd143)) + 
     -         5.333333333333333D0*(Cpl166*dZZA1*WProd143) + 
     -         3.3333333333333335D0*
     -          (Cpl160*(dZAZ1 + dZZA1)*(WProd137 + WProd138)*
     -            WProd143) - 
     -         5.333333333333333D0*
     -          (Cpl161*(dZAZ1 + dZZA1)*(WProd137 + WProd138)*
     -            WProd143))) - 
     -   5.333333333333333D0*
     -    (Cpl1*dZAA1*Prp29**2*(WProd137 + WProd138)*WProd143) + 
     -   Prp30**2*(2*Cpl167*dZZZ1*WProd143 + 
     -      2*Cpl4*(dMZsq1 - dZZZ1*WProd137 - dZZZ1*WProd138)*
     -       WProd143 + 5.333333333333333D0*
     -       (Cpl155*dZZZ1*WProd143) - 
     -      6.666666666666667D0*(Cpl159*dZZZ1*WProd143) + 
     -      5.333333333333333D0*
     -       (Cpl2*(dMZsq1 - dZZZ1*WProd137 - dZZZ1*WProd138)*
     -         WProd143) - 
     -      6.666666666666667D0*
     -       (Cpl3*(dMZsq1 - dZZZ1*WProd137 - dZZZ1*WProd138)*
     -         WProd143)) + 
     -   Prp29*(-(3.3333333333333335D0*(Cpl160*dZZA1*WProd143)) + 
     -      5.333333333333333D0*(Cpl161*dZZA1*WProd143) + 
     -      5.333333333333333D0*
     -       (Cpl1*WProd143*
     -         (dZAA1 + 2*dZe1 + dZfL1(2,1,1) + dZfL1(3,1,1))))

#ifdef DEBUG
	print *, 'amp(1,1,4) =', amp(1,1,4)
#endif

      endif

      end


