      subroutine ampsfacctEWnonphotonic2(amp,resnr)
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

        amp(2,1,1) = 2.6666666666666665D0*
     -    (Cpl1*dZAA1*Prp31**2*(WProd144 + WProd145)*WProd146) + 
     -   Prp32**2*(-(2.6666666666666665D0*(Cpl155*dZZZ1*WProd146)) - 
     -      2.6666666666666665D0*
     -       (Cpl2*(dMZsq1 - dZZZ1*WProd144 - dZZZ1*WProd145)*
     -         WProd146)) + 
     -   Prp31*(-(2.6666666666666665D0*(Cpl156*dZZA1*WProd146)) - 
     -      2.6666666666666665D0*
     -       (Cpl1*WProd146*
     -         (dZAA1 + 2*dZe1 + dZfR1(2,1,1) + dZfR1(4,1,1)))) + 
     -   Prp32*(-(5.333333333333333D0*(Cpl157*dSW1*WProd146)) - 
     -      2.6666666666666665D0*(Cpl156*dZAZ1*WProd146) + 
     -      Prp31*(-(2.6666666666666665D0*(Cpl158*dZZA1*WProd146)) + 
     -         2.6666666666666665D0*
     -          (Cpl156*(dZAZ1 + dZZA1)*(WProd144 + WProd145)*
     -            WProd146)) - 
     -      2.6666666666666665D0*
     -       (Cpl2*WProd146*
     -         (2*dZe1 + dZZZ1 + dZfR1(2,1,1) + dZfR1(4,1,1))))

#ifdef DEBUG
	print *, 'amp(2,1,1) =', amp(2,1,1)
#endif

	tmp4 = 2*dZe1 + dZZZ1 + dZfL1(2,1,1) + dZfR1(4,1,1)
        amp(2,1,2) = Prp32*
     -    (-(2.6666666666666665D0*(Cpl157*dSW1*WProd148)) + 
     -      1.3333333333333333D0*(Cpl162*dSW1*WProd148) - 
     -      1.3333333333333333D0*(Cpl163*dSW1*WProd148) - 
     -      1.3333333333333333D0*(Cpl164*dSW1*WProd148) - 
     -      1.3333333333333333D0*(Cpl156*dZAZ1*WProd148) + 
     -      0.6666666666666666D0*(Cpl160*dZAZ1*WProd148) - 
     -      1.3333333333333333D0*(Cpl161*dZAZ1*WProd148) - 
     -      2.6666666666666665D0*(Cpl2*tmp4*WProd148) + 
     -      1.3333333333333333D0*(Cpl3*tmp4*WProd148) + 
     -      Prp31*(-(1.3333333333333333D0*(Cpl158*dZZA1*WProd148)) + 
     -         0.6666666666666666D0*(Cpl165*dZZA1*WProd148) - 
     -         1.3333333333333333D0*(Cpl166*dZZA1*WProd148) + 
     -         1.3333333333333333D0*
     -          (Cpl156*(dZAZ1 + dZZA1)*(WProd144 + WProd145)*
     -            WProd148) - 
     -         0.6666666666666666D0*
     -          (Cpl160*(dZAZ1 + dZZA1)*(WProd144 + WProd145)*
     -            WProd148) + 
     -         1.3333333333333333D0*
     -          (Cpl161*(dZAZ1 + dZZA1)*(WProd144 + WProd145)*
     -            WProd148))) + 
     -   2.6666666666666665D0*
     -    (Cpl1*dZAA1*Prp31**2*(WProd144 + WProd145)*WProd148) + 
     -   Prp32**2*(-(2.6666666666666665D0*(Cpl155*dZZZ1*WProd148)) + 
     -      1.3333333333333333D0*(Cpl159*dZZZ1*WProd148) - 
     -      2.6666666666666665D0*
     -       (Cpl2*(dMZsq1 - dZZZ1*WProd144 - dZZZ1*WProd145)*
     -         WProd148) + 
     -      1.3333333333333333D0*
     -       (Cpl3*(dMZsq1 - dZZZ1*WProd144 - dZZZ1*WProd145)*
     -         WProd148)) + 
     -   Prp31*(-(1.3333333333333333D0*(Cpl156*dZZA1*WProd148)) + 
     -      0.6666666666666666D0*(Cpl160*dZZA1*WProd148) - 
     -      1.3333333333333333D0*(Cpl161*dZZA1*WProd148) - 
     -      2.6666666666666665D0*
     -       (Cpl1*WProd148*
     -         (dZAA1 + 2*dZe1 + dZfL1(2,1,1) + dZfR1(4,1,1))))

#ifdef DEBUG
	print *, 'amp(2,1,2) =', amp(2,1,2)
#endif

	tmp5 = 2*dZe1 + dZZZ1 + dZfL1(4,1,1) + dZfR1(2,1,1)
        amp(2,1,3) = Prp32*
     -    (4*Cpl162*dSW1*WProd149 - 4*Cpl163*dSW1*WProd149 + 
     -      2*Cpl160*dZAZ1*WProd149 + 4*Cpl3*tmp5*WProd149 - 
     -      2.6666666666666665D0*(Cpl157*dSW1*WProd149) + 
     -      1.3333333333333333D0*(Cpl164*dSW1*WProd149) - 
     -      1.3333333333333333D0*(Cpl156*dZAZ1*WProd149) - 
     -      1.3333333333333333D0*(Cpl161*dZAZ1*WProd149) - 
     -      2.6666666666666665D0*(Cpl2*tmp5*WProd149) + 
     -      Prp31*(2*Cpl165*dZZA1*WProd149 - 
     -         2*Cpl160*(dZAZ1 + dZZA1)*(WProd144 + WProd145)*
     -          WProd149 - 
     -         1.3333333333333333D0*(Cpl158*dZZA1*WProd149) - 
     -         1.3333333333333333D0*(Cpl166*dZZA1*WProd149) + 
     -         1.3333333333333333D0*
     -          (Cpl156*(dZAZ1 + dZZA1)*(WProd144 + WProd145)*
     -            WProd149) + 
     -         1.3333333333333333D0*
     -          (Cpl161*(dZAZ1 + dZZA1)*(WProd144 + WProd145)*
     -            WProd149))) + 
     -   2.6666666666666665D0*
     -    (Cpl1*dZAA1*Prp31**2*(WProd144 + WProd145)*WProd149) + 
     -   Prp32**2*(4*Cpl159*dZZZ1*WProd149 + 
     -      4*Cpl3*(dMZsq1 - dZZZ1*WProd144 - dZZZ1*WProd145)*
     -       WProd149 - 2.6666666666666665D0*
     -       (Cpl155*dZZZ1*WProd149) - 
     -      2.6666666666666665D0*
     -       (Cpl2*(dMZsq1 - dZZZ1*WProd144 - dZZZ1*WProd145)*
     -         WProd149)) + 
     -   Prp31*(2*Cpl160*dZZA1*WProd149 - 
     -      1.3333333333333333D0*(Cpl156*dZZA1*WProd149) - 
     -      1.3333333333333333D0*(Cpl161*dZZA1*WProd149) - 
     -      2.6666666666666665D0*
     -       (Cpl1*WProd149*
     -         (dZAA1 + 2*dZe1 + dZfL1(4,1,1) + dZfR1(2,1,1))))

#ifdef DEBUG
	print *, 'amp(2,1,3) =', amp(2,1,3)
#endif

	tmp6 = 2*dZe1 + dZZZ1 + dZfL1(2,1,1) + dZfL1(4,1,1)
        amp(2,1,4) = Prp32*
     -    (4*Cpl168*dSW1*WProd150 - 2*Cpl4*tmp6*WProd150 + 
     -      1.3333333333333333D0*(Cpl162*dSW1*WProd150) - 
     -      5.333333333333333D0*(Cpl163*dSW1*WProd150) + 
     -      2.6666666666666665D0*(Cpl160*dZAZ1*WProd150) - 
     -      2.6666666666666665D0*(Cpl161*dZAZ1*WProd150) - 
     -      2.6666666666666665D0*(Cpl2*tmp6*WProd150) + 
     -      5.333333333333333D0*(Cpl3*tmp6*WProd150) + 
     -      Prp31*(2.6666666666666665D0*(Cpl165*dZZA1*WProd150) - 
     -         2.6666666666666665D0*(Cpl166*dZZA1*WProd150) - 
     -         2.6666666666666665D0*
     -          (Cpl160*(dZAZ1 + dZZA1)*(WProd144 + WProd145)*
     -            WProd150) + 
     -         2.6666666666666665D0*
     -          (Cpl161*(dZAZ1 + dZZA1)*(WProd144 + WProd145)*
     -            WProd150))) + 
     -   2.6666666666666665D0*
     -    (Cpl1*dZAA1*Prp31**2*(WProd144 + WProd145)*WProd150) + 
     -   Prp32**2*(-2*Cpl167*dZZZ1*WProd150 - 
     -      2*Cpl4*(dMZsq1 - dZZZ1*WProd144 - dZZZ1*WProd145)*
     -       WProd150 - 2.6666666666666665D0*
     -       (Cpl155*dZZZ1*WProd150) + 
     -      5.333333333333333D0*(Cpl159*dZZZ1*WProd150) - 
     -      2.6666666666666665D0*
     -       (Cpl2*(dMZsq1 - dZZZ1*WProd144 - dZZZ1*WProd145)*
     -         WProd150) + 
     -      5.333333333333333D0*
     -       (Cpl3*(dMZsq1 - dZZZ1*WProd144 - dZZZ1*WProd145)*
     -         WProd150)) + 
     -   Prp31*(2.6666666666666665D0*(Cpl160*dZZA1*WProd150) - 
     -      2.6666666666666665D0*(Cpl161*dZZA1*WProd150) - 
     -      2.6666666666666665D0*
     -       (Cpl1*WProd150*
     -         (dZAA1 + 2*dZe1 + dZfL1(2,1,1) + dZfL1(4,1,1))))

#ifdef DEBUG
	print *, 'amp(2,1,4) =', amp(2,1,4)
#endif

      endif

      end


