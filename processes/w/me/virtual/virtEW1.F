      subroutine virtEW1(amp,resnr)
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
*     procnr=1:   { up  , ~do  ,  ve  , ~el  }
*
*   Helicities:
*     helnr=1:   {-1, 1, -1, 1}
*
*   Color Structures:
*     colnr=1:   SUNT[Col1, Col2]


      if(resnr.eq.0)then

        tmp1 = -(0.5D0*(Cpl87*WProd10)) + 
     -   1.3333333333333333D0*(Cpl88*WProd10) - 
     -   0.6666666666666666D0*(Cpl89*WProd10)
        tmp2 = Cpl88*WProd5*WProd6 - 0.5D0*(Cpl87*WProd5*WProd6) - 
     -   0.4444444444444444D0*(Cpl89*WProd5*WProd6)
	tmp3 = -WProd10 - WProd5*WProd6 + WProd5*WProd8
        tmp4 = WProd10 + 2*WProd5*WProd6 - 4*WProd5*WProd8 - 
     -   2*WProd5*WProd9
        tmp5 = -WProd10 - 2*WProd5*WProd6 + 4*WProd5*WProd8 + 
     -   2*WProd5*WProd9
        tmp6 = Cpl87*WProd5*(WProd6 - WProd8) - 
     -   0.6666666666666666D0*(Cpl88*WProd5*(WProd6 - WProd8))
        tmp7 = Cpl87*WProd5*(WProd6 - WProd8) - 
     -   3.3333333333333335D0*(Cpl88*WProd5*(WProd6 - WProd8)) + 
     -   2.6666666666666665D0*(Cpl89*WProd5*(WProd6 - WProd8))
        amp(1,1,1) = LInt101*tmp1 + LInt113*tmp1 + 
     -   LInt125*tmp1 + LInt109*tmp6 + LInt121*tmp6 + 
     -   LInt97*tmp6 + LInt111*tmp7 + LInt123*tmp7 + LInt99*tmp7 + 
     -   4*Cpl86*LInt82*WProd5 + 
     -   LInt84*(3*Cpl87*WProd5 - 4*Cpl88*WProd5) + 
     -   LInt83*(3*Cpl87*WProd5 - 8*Cpl88*WProd5 + 
     -      4*Cpl89*WProd5) - 
     -   0.6666666666666666D0*(Cpl86*LInt118*tmp4) + 
     -   LInt120*(0.6666666666666666D0*(Cpl88*tmp4) + 
     -      0.5D0*(Cpl87*tmp5)) + 
     -   LInt119*(1.3333333333333333D0*(Cpl88*tmp4) - 
     -      0.6666666666666666D0*(Cpl89*tmp4) + 0.5D0*(Cpl87*tmp5)) + 
     -   LInt102*(-(0.5D0*(Cpl87*WProd10)) + 
     -      0.6666666666666666D0*(Cpl88*WProd10)) + 
     -   LInt114*(-(0.5D0*(Cpl87*WProd10)) + 
     -      0.6666666666666666D0*(Cpl88*WProd10)) + 
     -   LInt126*(-(0.5D0*(Cpl87*WProd10)) + 
     -      0.6666666666666666D0*(Cpl88*WProd10)) - 
     -   0.6666666666666666D0*(Cpl86*LInt100*WProd10) - 
     -   0.6666666666666666D0*(Cpl86*LInt112*WProd10) - 
     -   0.6666666666666666D0*(Cpl86*LInt124*WProd10) + 
     -   LInt78*(-(Cpl87*WProd5) + 
     -      1.3333333333333333D0*(Cpl88*WProd5)) + 
     -   LInt79*(-4*Cpl87*WProd5 + 
     -      2.6666666666666665D0*(Cpl88*WProd5)) + 
     -   LInt77*(-(Cpl87*WProd5) + 
     -      2.6666666666666665D0*(Cpl88*WProd5) - 
     -      1.3333333333333333D0*(Cpl89*WProd5)) + 
     -   LInt81*(-4*Cpl87*WProd5 + 
     -      13.333333333333334D0*(Cpl88*WProd5) - 
     -      10.666666666666666D0*(Cpl89*WProd5)) - 
     -   1.3333333333333333D0*(Cpl86*LInt76*WProd5) - 
     -   10.666666666666666D0*(Cpl86*LInt80*WProd5) + 
     -   2.6666666666666665D0*
     -    (Cpl86*LInt110*WProd5*(WProd6 - WProd8)) + 
     -   2.6666666666666665D0*
     -    (Cpl86*LInt122*WProd5*(WProd6 - WProd8)) + 
     -   2.6666666666666665D0*
     -    (Cpl86*LInt98*WProd5*(WProd6 - WProd8)) + 
     -   LInt91*(Cpl87*(WProd10 + WProd5*WProd6 - WProd5*WProd8) - 
     -      0.6666666666666666D0*
     -       (Cpl88*(WProd10 + WProd5*WProd6 - WProd5*WProd8))) + 
     -   LInt96*(Cpl87*tmp3 + 
     -      1.3333333333333333D0*
     -       (Cpl88*(WProd10 + WProd5*WProd6 - WProd5*WProd8))) + 
     -   LInt95*(Cpl87*tmp3 + 
     -      2.6666666666666665D0*
     -       (Cpl88*(WProd10 + WProd5*WProd6 - WProd5*WProd8)) - 
     -      1.3333333333333333D0*
     -       (Cpl89*(WProd10 + WProd5*WProd6 - WProd5*WProd8))) + 
     -   LInt93*(Cpl87*(WProd10 + WProd5*WProd6 - WProd5*WProd8) - 
     -      3.3333333333333335D0*
     -       (Cpl88*(WProd10 + WProd5*WProd6 - WProd5*WProd8)) + 
     -      2.6666666666666665D0*
     -       (Cpl89*(WProd10 + WProd5*WProd6 - WProd5*WProd8))) + 
     -   2.6666666666666665D0*
     -    (Cpl86*LInt92*(WProd10 + WProd5*WProd6 - WProd5*WProd8))-
     -     1.3333333333333333D0*
     -    (Cpl86*LInt94*(WProd10 + WProd5*WProd6 - WProd5*WProd8))+
     -     LInt90*(-(0.5D0*(Cpl87*WProd5*
     -           (2*WProd6 - 3*WProd8 - 2*WProd9))) + 
     -      0.6666666666666666D0*
     -       (Cpl88*WProd5*(2*WProd6 - 3*WProd8 - 2*WProd9))) + 
     -   LInt89*(-(0.5D0*(Cpl87*WProd5*
     -           (2*WProd6 - 3*WProd8 - 2*WProd9))) + 
     -      1.3333333333333333D0*
     -       (Cpl88*WProd5*(2*WProd6 - 3*WProd8 - 2*WProd9)) - 
     -      0.6666666666666666D0*
     -       (Cpl89*WProd5*(2*WProd6 - 3*WProd8 - 2*WProd9))) + 
     -   LInt85*(0.5D0*(Cpl87*WProd5*
     -         (3*WProd6 - 3*WProd8 - 2*WProd9)) - 
     -      0.3333333333333333D0*
     -       (Cpl88*WProd5*(3*WProd6 - 3*WProd8 - 2*WProd9))) + 
     -   LInt87*(0.5D0*(Cpl87*WProd5*
     -         (3*WProd6 - 3*WProd8 - 2*WProd9)) - 
     -      1.6666666666666667D0*
     -       (Cpl88*WProd5*(3*WProd6 - 3*WProd8 - 2*WProd9)) + 
     -      1.3333333333333333D0*
     -       (Cpl89*WProd5*(3*WProd6 - 3*WProd8 - 2*WProd9))) + 
     -   LInt108*(0.5D0*(Cpl87*WProd5*(WProd6 - WProd9)) - 
     -      0.6666666666666666D0*(Cpl88*WProd5*(WProd6 - WProd9)))
        amp(1,1,1) = amp(1,1,1) - 
     -   0.6666666666666666D0*
     -    (Cpl86*LInt88*WProd5*(2*WProd6 - 3*WProd8 - 2*WProd9)) + 
     -   1.3333333333333333D0*
     -    (Cpl86*LInt86*WProd5*(3*WProd6 - 3*WProd8 - 2*WProd9)) + 
     -   LInt107*(0.5D0*(Cpl87*WProd5*(WProd6 - WProd9)) - 
     -      1.3333333333333333D0*(Cpl88*WProd5*(WProd6 - WProd9)) + 
     -      0.6666666666666666D0*(Cpl89*WProd5*(WProd6 - WProd9))) + 
     -   0.6666666666666666D0*
     -    (Cpl86*LInt106*WProd5*(WProd6 - WProd9)) + 
     -   LInt103*(0.5D0*(Cpl87*WProd5*WProd9) - 
     -      0.3333333333333333D0*(Cpl88*WProd5*WProd9)) + 
     -   LInt115*(Cpl87*WProd5*WProd9 - 
     -      0.6666666666666666D0*(Cpl88*WProd5*WProd9)) + 
     -   LInt105*(0.5D0*(Cpl87*WProd5*WProd9) - 
     -      1.6666666666666667D0*(Cpl88*WProd5*WProd9) + 
     -      1.3333333333333333D0*(Cpl89*WProd5*WProd9)) + 
     -   LInt117*(Cpl87*WProd5*WProd9 - 
     -      3.3333333333333335D0*(Cpl88*WProd5*WProd9) + 
     -      2.6666666666666665D0*(Cpl89*WProd5*WProd9)) + 
     -   1.3333333333333333D0*(Cpl86*LInt104*WProd5*WProd9) + 
     -   2.6666666666666665D0*(Cpl86*LInt116*WProd5*WProd9) + 
     -   Prp3*(LInt153*tmp2 + LInt163*tmp2 + Cpl88*WProd5 + 
     -      2*Cpl90*WProd5 - 4*Cpl86*LInt128*WProd5 - 
     -      2*Cpl90*LInt134*WProd5 - 4*Cpl86*LInt135*WProd5 + 
     -      LInt137*(-(Cpl87*WProd5) + 2*Cpl88*WProd5) + 
     -      LInt130*(4*Cpl86*WProd5 - 4*Cpl90*WProd5) + 
     -      LInt136*(4*Cpl86*WProd5 - 2*Cpl90*WProd5) + 
     -      LInt156*(2*Cpl86*WProd5*
     -          (2*WProd6 - 3*WProd8 - 3*WProd9) - 
     -         Cpl90*WProd5*(2*WProd6 - 3*WProd8 - 3*WProd9)) - 
     -      Cpl90*LInt154*WProd5*
     -       (2*WProd6 - 3*WProd8 - 3*WProd9) - 
     -      2*Cpl86*LInt155*WProd5*
     -       (2*WProd6 - 3*WProd8 - 3*WProd9) + 
     -      Cpl90*LInt144*WProd5*(WProd8 + WProd9) + 
     -      2*Cpl86*LInt145*WProd5*(WProd8 + WProd9) + 
     -      LInt146*(-2*Cpl86*WProd5*(WProd8 + WProd9) + 
     -         Cpl90*WProd5*(WProd8 + WProd9)) - 
     -      0.2222222222222222D0*(Cpl86*WProd5) + 
     -      LInt140*(-2*Cpl90*WProd5 + 
     -         1.3333333333333333D0*(Cpl86*WProd5)) + 
     -      LInt142*(-2*Cpl90*WProd5 + 
     -         2.6666666666666665D0*(Cpl86*WProd5)) - 
     -      0.5D0*(Cpl87*WProd5) + 
     -      LInt127*(-(Cpl88*WProd5) + 0.5D0*(Cpl87*WProd5)) - 
     -      0.2222222222222222D0*(Cpl89*WProd5) + 
     -      LInt129*(-(Cpl88*WProd5) + 
     -         0.4444444444444444D0*(Cpl86*WProd5) + 
     -         0.5D0*(Cpl87*WProd5) + 
     -         0.4444444444444444D0*(Cpl89*WProd5)) + 
     -      LInt143*(-(Cpl87*WProd5) + 2*Cpl88*WProd5 - 
     -         0.8888888888888888D0*(Cpl89*WProd5)) - 
     -      0.8888888888888888D0*(Cpl86*LInt138*WProd5) - 
     -      1.3333333333333333D0*(Cpl86*LInt139*WProd5) - 
     -      2.6666666666666665D0*(Cpl86*LInt141*WProd5) + 
     -      LInt150*(Cpl90*WProd5*WProd6 - 
     -         0.6666666666666666D0*(Cpl86*WProd5*WProd6)) + 
     -      LInt160*(Cpl90*WProd5*WProd6 - 
     -         0.6666666666666666D0*(Cpl86*WProd5*WProd6)) + 
     -      LInt152*(Cpl90*WProd5*WProd6 - 
     -         1.3333333333333333D0*(Cpl86*WProd5*WProd6)) + 
     -      LInt162*(Cpl90*WProd5*WProd6 - 
     -         1.3333333333333333D0*(Cpl86*WProd5*WProd6)) + 
     -      LInt133*(-(Cpl93*WProd5) + Cpl88*WProd5*WProd6 + 
     -         0.4444444444444444D0*(Cpl91*WProd5) + 
     -         0.5D0*(Cpl92*WProd5) - 0.5D0*(Cpl87*WProd5*WProd6) - 
     -         0.4444444444444444D0*(Cpl89*WProd5*WProd6)) - 
     -      0.4444444444444444D0*(Cpl86*LInt132*WProd5*WProd6) - 
     -      0.4444444444444444D0*(Cpl86*LInt148*WProd5*WProd6) + 
     -      0.6666666666666666D0*(Cpl86*LInt149*WProd5*WProd6) + 
     -      1.3333333333333333D0*(Cpl86*LInt151*WProd5*WProd6) - 
     -      0.4444444444444444D0*(Cpl86*LInt158*WProd5*WProd6) + 
     -      0.6666666666666666D0*(Cpl86*LInt159*WProd5*WProd6) + 
     -      1.3333333333333333D0*(Cpl86*LInt161*WProd5*WProd6) + 
     -      LInt157*(-(Cpl88*WProd5*
     -            (2*WProd6 - 3*WProd8 - 3*WProd9)) + 
     -         0.5D0*(Cpl87*WProd5*(2*WProd6 - 3*WProd8 - 3*WProd9)))
     -        + LInt147*(Cpl88*WProd5*(WProd8 + WProd9) - 
     -         0.5D0*(Cpl87*WProd5*(WProd8 + WProd9))) + 
     -      LInt131*(-(Cpl93*WProd5) + 
     -         Cpl88*WProd5*(WProd8 + WProd9) + 
     -         0.5D0*(Cpl92*WProd5) - 
     -         0.5D0*(Cpl87*WProd5*(WProd8 + WProd9))))
        amp(1,1,1) = amp(1,1,1) + 
     -   Prp3**2*(-2*Cpl95*WProd5 + 2*Cpl96*WProd5 + 
     -      Cpl95*LInt168*WProd5 - Cpl90*LInt169*WProd5 + 
     -      LInt170*(-(Cpl90*WProd5) - 8*Cpl94*WProd5) - 
     -      2*Cpl94*LInt172*WProd5*WProd6 + 
     -      LInt128*(Cpl99*WProd5 - 
     -         5*Cpl86*WProd5*(WProd8 + WProd9)) + 
     -      LInt130*(Cpl97*WProd5 - 2*Cpl98*WProd5 - 
     -         5*Cpl94*WProd5*(WProd8 + WProd9)) + 
     -      LInt178*(3*Cpl94*WProd5 + 0.25D0*(Cpl90*WProd5)) + 
     -      LInt177*(-2*Cpl86*WProd5 - 2*Cpl94*WProd5 + 
     -         3.5D0*(Cpl90*WProd5)) + 
     -      0.6666666666666666D0*(Cpl99*WProd5) + 
     -      0.25D0*(Cpl90*LInt176*WProd5) - 
     -      1.3333333333333333D0*(Cpl99*LInt179*WProd5) - 
     -      0.2222222222222222D0*(Cpl86*WProd5*WProd6) + 
     -      LInt171*(1.3333333333333333D0*(Cpl99*WProd5) - 
     -         3.3333333333333335D0*(Cpl86*WProd5*WProd6)) - 
     -      0.6666666666666666D0*(Cpl94*WProd5*WProd6))

	do Gen1 = 1, 3
        amp(1,1,1) = amp(1,1,1) + 
     -   Prp3**2*(WProd5*Cpl100(Gen1) - 
     -      0.4444444444444444D0*(Cpl90*WProd5*WProd6) - 
     -      0.6666666666666666D0*
     -       (Cpl90*WProd5*WProd6*LInt174(Gen1)) - 
     -      WProd5*Cpl100(Gen1)*LInt173(Gen1) + 
     -      (-2*Cpl90*WProd5*WProd6 + WProd5*Cpl100(Gen1))*
     -       LInt175(Gen1))
	enddo

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

      endif

      end


