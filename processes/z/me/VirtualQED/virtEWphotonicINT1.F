      subroutine virtEWphotonicINT1(amp,resnr)
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

        tmp80 = WProd45*WProd46 - WProd45*WProd47 + 
     -   2*WProd45*WProd49 - 2*WProd50 + WProd51 - WProd52
        tmp81 = 2*WProd45*WProd47 - 2*WProd45*WProd49 + 
     -   2*WProd50 - 2*WProd51 + WProd52
        tmp82 = WProd45*WProd46 - WProd45*WProd47 - 
     -   3*WProd45*WProd49 + WProd50 - WProd51 + WProd52
        tmp83 = WProd45*WProd47 - 2*WProd45*WProd49 + 
     -   2*WProd50 - WProd51 + WProd52
        tmp84 = 2*WProd45*WProd46 - 2*WProd45*WProd47 - 
     -   6*WProd45*WProd49 + 4*WProd50 - 3*WProd51 + 2*WProd52
        tmp85 = 2*WProd45*WProd46 - 2*WProd45*WProd47 - 
     -   4*WProd45*WProd49 + 2*WProd50 - WProd51 + 2*WProd52
        amp(1,1,1) = -(1.7777777777777777D0*
     -      (Cpl92*LInt113*tmp80)) - 
     -   1.7777777777777777D0*(Cpl93*LInt114*tmp80) - 
     -   1.7777777777777777D0*(Cpl93*LInt115*tmp80) + 
     -   1.7777777777777777D0*(Cpl92*LInt146*tmp81) + 
     -   1.7777777777777777D0*(Cpl93*LInt147*tmp81) + 
     -   1.7777777777777777D0*(Cpl93*LInt148*tmp81) + 
     -   3.5555555555555554D0*(Cpl92*LInt122*tmp82) + 
     -   3.5555555555555554D0*(Cpl93*LInt123*tmp82) + 
     -   3.5555555555555554D0*(Cpl93*LInt124*tmp82) + 
     -   3.5555555555555554D0*(Cpl92*LInt125*tmp82) + 
     -   3.5555555555555554D0*(Cpl93*LInt126*tmp82) + 
     -   3.5555555555555554D0*(Cpl93*LInt127*tmp82) + 
     -   1.7777777777777777D0*(Cpl92*LInt110*tmp83) + 
     -   1.7777777777777777D0*(Cpl93*LInt111*tmp83) + 
     -   1.7777777777777777D0*(Cpl93*LInt112*tmp83) + 
     -   1.7777777777777777D0*(Cpl92*LInt128*tmp84) + 
     -   1.7777777777777777D0*(Cpl93*LInt129*tmp84) + 
     -   1.7777777777777777D0*(Cpl93*LInt130*tmp84) + 
     -   1.7777777777777777D0*(Cpl92*LInt140*tmp85) + 
     -   1.7777777777777777D0*(Cpl93*LInt141*tmp85) + 
     -   1.7777777777777777D0*(Cpl93*LInt142*tmp85) + 
     -   1.7777777777777777D0*(Cpl92*LInt152*tmp85) + 
     -   1.7777777777777777D0*(Cpl93*LInt153*tmp85) + 
     -   1.7777777777777777D0*(Cpl93*LInt154*tmp85) - 
     -   3.5555555555555554D0*(Cpl92*LInt101*WProd45) - 
     -   3.5555555555555554D0*(Cpl93*LInt102*WProd45) - 
     -   3.5555555555555554D0*(Cpl93*LInt103*WProd45) + 
     -   10.666666666666666D0*(Cpl92*LInt104*WProd45) + 
     -   10.666666666666666D0*(Cpl93*LInt105*WProd45) + 
     -   10.666666666666666D0*(Cpl93*LInt106*WProd45) + 
     -   14.222222222222221D0*(Cpl92*LInt107*WProd45) + 
     -   14.222222222222221D0*(Cpl93*LInt108*WProd45) + 
     -   14.222222222222221D0*(Cpl93*LInt109*WProd45) + 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt134*WProd45*(WProd46 - WProd49)) + 
     -   1.7777777777777777D0*
     -    (Cpl93*LInt135*WProd45*(WProd46 - WProd49)) + 
     -   1.7777777777777777D0*
     -    (Cpl93*LInt136*WProd45*(WProd46 - WProd49)) - 
     -   1.7777777777777777D0*(Cpl92*LInt137*WProd45*WProd49) - 
     -   1.7777777777777777D0*(Cpl93*LInt138*WProd45*WProd49) - 
     -   1.7777777777777777D0*(Cpl93*LInt139*WProd45*WProd49) - 
     -   3.5555555555555554D0*(Cpl92*LInt143*WProd45*WProd49) - 
     -   3.5555555555555554D0*(Cpl93*LInt144*WProd45*WProd49) - 
     -   3.5555555555555554D0*(Cpl93*LInt145*WProd45*WProd49) - 
     -   3.5555555555555554D0*(Cpl92*LInt149*WProd45*WProd49) - 
     -   3.5555555555555554D0*(Cpl93*LInt150*WProd45*WProd49) - 
     -   3.5555555555555554D0*(Cpl93*LInt151*WProd45*WProd49) - 
     -   3.5555555555555554D0*(Cpl92*LInt155*WProd45*WProd49) - 
     -   3.5555555555555554D0*(Cpl93*LInt156*WProd45*WProd49) - 
     -   3.5555555555555554D0*(Cpl93*LInt157*WProd45*WProd49) - 
     -   3.5555555555555554D0*
     -    (Cpl92*LInt116*(WProd45*WProd49 - WProd50 + WProd51)) - 
     -   3.5555555555555554D0*
     -    (Cpl93*LInt117*(WProd45*WProd49 - WProd50 + WProd51)) - 
     -   3.5555555555555554D0*
     -    (Cpl93*LInt118*(WProd45*WProd49 - WProd50 + WProd51)) - 
     -   3.5555555555555554D0*
     -    (Cpl92*LInt119*(WProd45*WProd49 - WProd50 + WProd51)) - 
     -   3.5555555555555554D0*
     -    (Cpl93*LInt120*(WProd45*WProd49 - WProd50 + WProd51)) - 
     -   3.5555555555555554D0*
     -    (Cpl93*LInt121*(WProd45*WProd49 - WProd50 + WProd51)) - 
     -   3.5555555555555554D0*
     -    (Cpl92*LInt131*(2*WProd45*WProd49 - WProd50 + WProd51))-
     -     3.5555555555555554D0*
     -    (Cpl93*LInt132*(2*WProd45*WProd49 - WProd50 + WProd51))-
     -     3.5555555555555554D0*
     -    (Cpl93*LInt133*(2*WProd45*WProd49 - WProd50 + WProd51))

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        tmp86 = -(1.7777777777777777D0*(Cpl93*WProd49*WProd53)) + 
     -   0.8888888888888888D0*(Cpl94*WProd49*WProd53)
        tmp87 = WProd46*WProd53 - WProd47*WProd53 + 
     -   2*WProd49*WProd53 + WProd54 - WProd55
        tmp88 = WProd46*WProd53 - WProd47*WProd53 - 
     -   2*WProd49*WProd53 + WProd55
        tmp89 = WProd46*WProd53 - WProd47*WProd53 - 
     -   3*WProd49*WProd53 - WProd54 + WProd55
	tmp90 = WProd47*WProd53 - 2*WProd49*WProd53 - WProd54 + WProd55
        tmp91 = 2*WProd47*WProd53 - 2*WProd49*WProd53 - 
     -   WProd54 + WProd55
        tmp92 = 1.7777777777777777D0*
     -    (Cpl93*(WProd46 - WProd49)*WProd53) - 
     -   0.8888888888888888D0*(Cpl94*(WProd46 - WProd49)*WProd53)
        tmp93 = -(3.5555555555555554D0*
     -      (Cpl93*(WProd49*WProd53 + WProd54))) + 
     -   1.7777777777777777D0*(Cpl94*(WProd49*WProd53 + WProd54))
        tmp94 = -(1.7777777777777777D0*
     -      (Cpl93*(2*WProd49*WProd53 + WProd54))) + 
     -   0.8888888888888888D0*(Cpl94*(2*WProd49*WProd53 + WProd54))
        tmp95 = -(1.7777777777777777D0*
     -      (Cpl93*(4*WProd49*WProd53 + 3*WProd54))) + 
     -   0.8888888888888888D0*(Cpl94*(4*WProd49*WProd53 + 3*WProd54))
        amp(1,1,2) = LInt138*tmp86 + LInt139*tmp86 + 
     -   LInt135*tmp92 + LInt136*tmp92 + LInt117*tmp93 + 
     -   LInt118*tmp93 + LInt120*tmp93 + LInt121*tmp93 + 
     -   LInt144*tmp94 + LInt145*tmp94 + LInt150*tmp94 + 
     -   LInt151*tmp94 + LInt156*tmp94 + LInt157*tmp94 + 
     -   LInt132*tmp95 + LInt133*tmp95 + 
     -   LInt114*(-(1.7777777777777777D0*(Cpl93*tmp87)) + 
     -      0.8888888888888888D0*(Cpl94*tmp87)) + 
     -   LInt115*(-(1.7777777777777777D0*(Cpl93*tmp87)) + 
     -      0.8888888888888888D0*(Cpl94*tmp87)) - 
     -   1.7777777777777777D0*(Cpl92*LInt113*tmp87) + 
     -   LInt141*(3.5555555555555554D0*(Cpl93*tmp88) - 
     -      1.7777777777777777D0*(Cpl94*tmp88)) + 
     -   LInt142*(3.5555555555555554D0*(Cpl93*tmp88) - 
     -      1.7777777777777777D0*(Cpl94*tmp88)) + 
     -   LInt153*(3.5555555555555554D0*(Cpl93*tmp88) - 
     -      1.7777777777777777D0*(Cpl94*tmp88)) + 
     -   LInt154*(3.5555555555555554D0*(Cpl93*tmp88) - 
     -      1.7777777777777777D0*(Cpl94*tmp88)) + 
     -   3.5555555555555554D0*(Cpl92*LInt140*tmp88) + 
     -   3.5555555555555554D0*(Cpl92*LInt152*tmp88) + 
     -   LInt123*(3.5555555555555554D0*(Cpl93*tmp89) - 
     -      1.7777777777777777D0*(Cpl94*tmp89)) + 
     -   LInt124*(3.5555555555555554D0*(Cpl93*tmp89) - 
     -      1.7777777777777777D0*(Cpl94*tmp89)) + 
     -   LInt126*(3.5555555555555554D0*(Cpl93*tmp89) - 
     -      1.7777777777777777D0*(Cpl94*tmp89)) + 
     -   LInt127*(3.5555555555555554D0*(Cpl93*tmp89) - 
     -      1.7777777777777777D0*(Cpl94*tmp89)) + 
     -   LInt129*(3.5555555555555554D0*(Cpl93*tmp89) - 
     -      1.7777777777777777D0*(Cpl94*tmp89)) + 
     -   LInt130*(3.5555555555555554D0*(Cpl93*tmp89) - 
     -      1.7777777777777777D0*(Cpl94*tmp89)) + 
     -   3.5555555555555554D0*(Cpl92*LInt122*tmp89) + 
     -   3.5555555555555554D0*(Cpl92*LInt125*tmp89) + 
     -   3.5555555555555554D0*(Cpl92*LInt128*tmp89) + 
     -   LInt111*(1.7777777777777777D0*(Cpl93*tmp90) - 
     -      0.8888888888888888D0*(Cpl94*tmp90)) + 
     -   LInt112*(1.7777777777777777D0*(Cpl93*tmp90) - 
     -      0.8888888888888888D0*(Cpl94*tmp90)) + 
     -   1.7777777777777777D0*(Cpl92*LInt110*tmp90) + 
     -   LInt147*(1.7777777777777777D0*(Cpl93*tmp91) - 
     -      0.8888888888888888D0*(Cpl94*tmp91)) + 
     -   LInt148*(1.7777777777777777D0*(Cpl93*tmp91) - 
     -      0.8888888888888888D0*(Cpl94*tmp91)) + 
     -   1.7777777777777777D0*(Cpl92*LInt146*tmp91) + 
     -   LInt108*(3.5555555555555554D0*(Cpl93*WProd53) - 
     -      1.7777777777777777D0*(Cpl94*WProd53)) + 
     -   LInt109*(3.5555555555555554D0*(Cpl93*WProd53) - 
     -      1.7777777777777777D0*(Cpl94*WProd53)) + 
     -   LInt102*(-(3.5555555555555554D0*(Cpl93*WProd53)) + 
     -      1.7777777777777777D0*(Cpl94*WProd53)) + 
     -   LInt103*(-(3.5555555555555554D0*(Cpl93*WProd53)) + 
     -      1.7777777777777777D0*(Cpl94*WProd53)) - 
     -   3.5555555555555554D0*(Cpl92*LInt101*WProd53) + 
     -   3.5555555555555554D0*(Cpl92*LInt107*WProd53) + 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt134*(WProd46 - WProd49)*WProd53) - 
     -   1.7777777777777777D0*(Cpl92*LInt137*WProd49*WProd53) - 
     -   3.5555555555555554D0*
     -    (Cpl92*LInt116*(WProd49*WProd53 + WProd54)) - 
     -   3.5555555555555554D0*
     -    (Cpl92*LInt119*(WProd49*WProd53 + WProd54)) - 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt143*(2*WProd49*WProd53 + WProd54)) - 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt149*(2*WProd49*WProd53 + WProd54)) - 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt155*(2*WProd49*WProd53 + WProd54)) - 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt131*(4*WProd49*WProd53 + 3*WProd54))

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

        tmp96 = -(1.7777777777777777D0*(Cpl93*WProd49*WProd56)) + 
     -   1.3333333333333333D0*(Cpl94*WProd49*WProd56)
        tmp97 = -2*WProd47*WProd56 + 2*WProd49*WProd56 + 
     -   WProd57 - WProd58
        tmp98 = -(WProd47*WProd56) + 2*WProd49*WProd56 + 
     -   WProd57 - WProd58
        tmp99 = WProd46*WProd56 - WProd47*WProd56 + 
     -   2*WProd49*WProd56 + WProd57 - WProd58
        tmp100 = -(WProd46*WProd56) + WProd47*WProd56 + 
     -   3*WProd49*WProd56 + WProd57 - WProd58
        tmp101 = WProd46*WProd56 - WProd47*WProd56 - 
     -   2*WProd49*WProd56 + WProd58
        tmp102 = 1.7777777777777777D0*
     -    (Cpl93*(WProd46 - WProd49)*WProd56) - 
     -   1.3333333333333333D0*(Cpl94*(WProd46 - WProd49)*WProd56)
        tmp103 = -(3.5555555555555554D0*
     -      (Cpl93*(WProd49*WProd56 + WProd57))) + 
     -   2.6666666666666665D0*(Cpl94*(WProd49*WProd56 + WProd57))
        tmp104 = -(1.7777777777777777D0*
     -      (Cpl93*(2*WProd49*WProd56 + WProd57))) + 
     -   1.3333333333333333D0*(Cpl94*(2*WProd49*WProd56 + WProd57))
        tmp105 = -(1.7777777777777777D0*
     -      (Cpl93*(4*WProd49*WProd56 + 3*WProd57))) + 
     -   1.3333333333333333D0*(Cpl94*(4*WProd49*WProd56 + 3*WProd57))
        amp(1,1,3) = LInt135*tmp102 + LInt136*tmp102 + 
     -   LInt117*tmp103 + LInt118*tmp103 + LInt120*tmp103 + 
     -   LInt121*tmp103 + LInt144*tmp104 + LInt145*tmp104 + 
     -   LInt150*tmp104 + LInt151*tmp104 + LInt156*tmp104 + 
     -   LInt157*tmp104 + LInt132*tmp105 + LInt133*tmp105 + 
     -   LInt138*tmp96 + LInt139*tmp96 + 
     -   LInt123*(-(3.5555555555555554D0*(Cpl93*tmp100)) + 
     -      2.6666666666666665D0*(Cpl94*tmp100)) + 
     -   LInt124*(-(3.5555555555555554D0*(Cpl93*tmp100)) + 
     -      2.6666666666666665D0*(Cpl94*tmp100)) + 
     -   LInt126*(-(3.5555555555555554D0*(Cpl93*tmp100)) + 
     -      2.6666666666666665D0*(Cpl94*tmp100)) + 
     -   LInt127*(-(3.5555555555555554D0*(Cpl93*tmp100)) + 
     -      2.6666666666666665D0*(Cpl94*tmp100)) + 
     -   LInt129*(-(3.5555555555555554D0*(Cpl93*tmp100)) + 
     -      2.6666666666666665D0*(Cpl94*tmp100)) + 
     -   LInt130*(-(3.5555555555555554D0*(Cpl93*tmp100)) + 
     -      2.6666666666666665D0*(Cpl94*tmp100)) - 
     -   3.5555555555555554D0*(Cpl92*LInt122*tmp100) - 
     -   3.5555555555555554D0*(Cpl92*LInt125*tmp100) - 
     -   3.5555555555555554D0*(Cpl92*LInt128*tmp100) + 
     -   LInt141*(3.5555555555555554D0*(Cpl93*tmp101) - 
     -      2.6666666666666665D0*(Cpl94*tmp101)) + 
     -   LInt142*(3.5555555555555554D0*(Cpl93*tmp101) - 
     -      2.6666666666666665D0*(Cpl94*tmp101)) + 
     -   LInt153*(3.5555555555555554D0*(Cpl93*tmp101) - 
     -      2.6666666666666665D0*(Cpl94*tmp101)) + 
     -   LInt154*(3.5555555555555554D0*(Cpl93*tmp101) - 
     -      2.6666666666666665D0*(Cpl94*tmp101)) + 
     -   3.5555555555555554D0*(Cpl92*LInt140*tmp101) + 
     -   3.5555555555555554D0*(Cpl92*LInt152*tmp101) + 
     -   LInt147*(-(1.7777777777777777D0*(Cpl93*tmp97)) + 
     -      1.3333333333333333D0*(Cpl94*tmp97)) + 
     -   LInt148*(-(1.7777777777777777D0*(Cpl93*tmp97)) + 
     -      1.3333333333333333D0*(Cpl94*tmp97)) - 
     -   1.7777777777777777D0*(Cpl92*LInt146*tmp97) + 
     -   LInt111*(-(1.7777777777777777D0*(Cpl93*tmp98)) + 
     -      1.3333333333333333D0*(Cpl94*tmp98)) + 
     -   LInt112*(-(1.7777777777777777D0*(Cpl93*tmp98)) + 
     -      1.3333333333333333D0*(Cpl94*tmp98)) - 
     -   1.7777777777777777D0*(Cpl92*LInt110*tmp98) + 
     -   LInt114*(-(1.7777777777777777D0*(Cpl93*tmp99)) + 
     -      1.3333333333333333D0*(Cpl94*tmp99)) + 
     -   LInt115*(-(1.7777777777777777D0*(Cpl93*tmp99)) + 
     -      1.3333333333333333D0*(Cpl94*tmp99)) - 
     -   1.7777777777777777D0*(Cpl92*LInt113*tmp99) + 
     -   LInt108*(3.5555555555555554D0*(Cpl93*WProd56) - 
     -      2.6666666666666665D0*(Cpl94*WProd56)) + 
     -   LInt109*(3.5555555555555554D0*(Cpl93*WProd56) - 
     -      2.6666666666666665D0*(Cpl94*WProd56)) + 
     -   LInt102*(-(3.5555555555555554D0*(Cpl93*WProd56)) + 
     -      2.6666666666666665D0*(Cpl94*WProd56)) + 
     -   LInt103*(-(3.5555555555555554D0*(Cpl93*WProd56)) + 
     -      2.6666666666666665D0*(Cpl94*WProd56)) - 
     -   3.5555555555555554D0*(Cpl92*LInt101*WProd56) + 
     -   3.5555555555555554D0*(Cpl92*LInt107*WProd56) + 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt134*(WProd46 - WProd49)*WProd56) - 
     -   1.7777777777777777D0*(Cpl92*LInt137*WProd49*WProd56) - 
     -   3.5555555555555554D0*
     -    (Cpl92*LInt116*(WProd49*WProd56 + WProd57)) - 
     -   3.5555555555555554D0*
     -    (Cpl92*LInt119*(WProd49*WProd56 + WProd57)) - 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt143*(2*WProd49*WProd56 + WProd57)) - 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt149*(2*WProd49*WProd56 + WProd57)) - 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt155*(2*WProd49*WProd56 + WProd57)) - 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt131*(4*WProd49*WProd56 + 3*WProd57))

#ifdef DEBUG
	print *, 'amp(1,1,3) =', amp(1,1,3)
#endif

        tmp106 = -(3.5555555555555554D0*(Cpl93*WProd59)) + 
     -   4.444444444444445D0*(Cpl94*WProd59) - 
     -   1.3333333333333333D0*(Cpl95*WProd59)
        tmp107 = 4*Cpl95*WProd59 + 
     -   10.666666666666666D0*(Cpl93*WProd59) - 
     -   13.333333333333334D0*(Cpl94*WProd59)
        tmp108 = 14.222222222222221D0*(Cpl93*WProd59) - 
     -   17.77777777777778D0*(Cpl94*WProd59) + 
     -   5.333333333333333D0*(Cpl95*WProd59)
        tmp109 = -(3.5555555555555554D0*(Cpl93*WProd49*WProd59)) + 
     -   4.444444444444445D0*(Cpl94*WProd49*WProd59) - 
     -   1.3333333333333333D0*(Cpl95*WProd49*WProd59)
        tmp110 = -(1.7777777777777777D0*(Cpl93*WProd49*WProd59)) + 
     -   2.2222222222222223D0*(Cpl94*WProd49*WProd59) - 
     -   0.6666666666666666D0*(Cpl95*WProd49*WProd59)
        tmp111 = WProd46*WProd59 - WProd47*WProd59 + 
     -   2*WProd49*WProd59 - 2*WProd60 + WProd61 - WProd62
        tmp112 = 2*WProd47*WProd59 - 2*WProd49*WProd59 + 
     -   2*WProd60 - 2*WProd61 + WProd62
        tmp113 = WProd46*WProd59 - WProd47*WProd59 - 
     -   3*WProd49*WProd59 + WProd60 - WProd61 + WProd62
        tmp114 = WProd47*WProd59 - 2*WProd49*WProd59 + 
     -   2*WProd60 - WProd61 + WProd62
        tmp115 = 2*WProd46*WProd59 - 2*WProd47*WProd59 - 
     -   6*WProd49*WProd59 + 4*WProd60 - 3*WProd61 + 2*WProd62
        tmp116 = 2*WProd46*WProd59 - 2*WProd47*WProd59 - 
     -   4*WProd49*WProd59 + 2*WProd60 - WProd61 + 2*WProd62
        tmp117 = 1.7777777777777777D0*
     -    (Cpl93*(WProd46 - WProd49)*WProd59) - 
     -   2.2222222222222223D0*(Cpl94*(WProd46 - WProd49)*WProd59) + 
     -   0.6666666666666666D0*(Cpl95*(WProd46 - WProd49)*WProd59)
        tmp118 = -(3.5555555555555554D0*
     -      (Cpl93*(WProd49*WProd59 - WProd60 + WProd61))) + 
     -   4.444444444444445D0*
     -    (Cpl94*(WProd49*WProd59 - WProd60 + WProd61)) - 
     -   1.3333333333333333D0*
     -    (Cpl95*(WProd49*WProd59 - WProd60 + WProd61))
        tmp119 = -(3.5555555555555554D0*
     -      (Cpl93*(2*WProd49*WProd59 - WProd60 + WProd61))) + 
     -   4.444444444444445D0*
     -    (Cpl94*(2*WProd49*WProd59 - WProd60 + WProd61)) - 
     -   1.3333333333333333D0*
     -    (Cpl95*(2*WProd49*WProd59 - WProd60 + WProd61))
        amp(1,1,4) = LInt102*tmp106 + LInt103*tmp106 + 
     -   LInt105*tmp107 + LInt106*tmp107 + LInt108*tmp108 + 
     -   LInt109*tmp108 + LInt144*tmp109 + LInt145*tmp109 + 
     -   LInt150*tmp109 + LInt151*tmp109 + LInt156*tmp109 + 
     -   LInt157*tmp109 + LInt138*tmp110 + LInt139*tmp110 + 
     -   LInt135*tmp117 + LInt136*tmp117 + LInt117*tmp118 + 
     -   LInt118*tmp118 + LInt120*tmp118 + LInt121*tmp118 + 
     -   LInt132*tmp119 + LInt133*tmp119 + 
     -   LInt114*(-(1.7777777777777777D0*(Cpl93*tmp111)) + 
     -      2.2222222222222223D0*(Cpl94*tmp111) - 
     -      0.6666666666666666D0*(Cpl95*tmp111)) + 
     -   LInt115*(-(1.7777777777777777D0*(Cpl93*tmp111)) + 
     -      2.2222222222222223D0*(Cpl94*tmp111) - 
     -      0.6666666666666666D0*(Cpl95*tmp111)) - 
     -   1.7777777777777777D0*(Cpl92*LInt113*tmp111) + 
     -   LInt147*(1.7777777777777777D0*(Cpl93*tmp112) - 
     -      2.2222222222222223D0*(Cpl94*tmp112) + 
     -      0.6666666666666666D0*(Cpl95*tmp112)) + 
     -   LInt148*(1.7777777777777777D0*(Cpl93*tmp112) - 
     -      2.2222222222222223D0*(Cpl94*tmp112) + 
     -      0.6666666666666666D0*(Cpl95*tmp112)) + 
     -   1.7777777777777777D0*(Cpl92*LInt146*tmp112) + 
     -   LInt123*(3.5555555555555554D0*(Cpl93*tmp113) - 
     -      4.444444444444445D0*(Cpl94*tmp113) + 
     -      1.3333333333333333D0*(Cpl95*tmp113)) + 
     -   LInt124*(3.5555555555555554D0*(Cpl93*tmp113) - 
     -      4.444444444444445D0*(Cpl94*tmp113) + 
     -      1.3333333333333333D0*(Cpl95*tmp113)) + 
     -   LInt126*(3.5555555555555554D0*(Cpl93*tmp113) - 
     -      4.444444444444445D0*(Cpl94*tmp113) + 
     -      1.3333333333333333D0*(Cpl95*tmp113)) + 
     -   LInt127*(3.5555555555555554D0*(Cpl93*tmp113) - 
     -      4.444444444444445D0*(Cpl94*tmp113) + 
     -      1.3333333333333333D0*(Cpl95*tmp113)) + 
     -   3.5555555555555554D0*(Cpl92*LInt122*tmp113) + 
     -   3.5555555555555554D0*(Cpl92*LInt125*tmp113) + 
     -   LInt111*(1.7777777777777777D0*(Cpl93*tmp114) - 
     -      2.2222222222222223D0*(Cpl94*tmp114) + 
     -      0.6666666666666666D0*(Cpl95*tmp114)) + 
     -   LInt112*(1.7777777777777777D0*(Cpl93*tmp114) - 
     -      2.2222222222222223D0*(Cpl94*tmp114) + 
     -      0.6666666666666666D0*(Cpl95*tmp114)) + 
     -   1.7777777777777777D0*(Cpl92*LInt110*tmp114) + 
     -   LInt129*(1.7777777777777777D0*(Cpl93*tmp115) - 
     -      2.2222222222222223D0*(Cpl94*tmp115) + 
     -      0.6666666666666666D0*(Cpl95*tmp115)) + 
     -   LInt130*(1.7777777777777777D0*(Cpl93*tmp115) - 
     -      2.2222222222222223D0*(Cpl94*tmp115) + 
     -      0.6666666666666666D0*(Cpl95*tmp115)) + 
     -   1.7777777777777777D0*(Cpl92*LInt128*tmp115) + 
     -   LInt141*(1.7777777777777777D0*(Cpl93*tmp116) - 
     -      2.2222222222222223D0*(Cpl94*tmp116) + 
     -      0.6666666666666666D0*(Cpl95*tmp116)) + 
     -   LInt142*(1.7777777777777777D0*(Cpl93*tmp116) - 
     -      2.2222222222222223D0*(Cpl94*tmp116) + 
     -      0.6666666666666666D0*(Cpl95*tmp116)) + 
     -   LInt153*(1.7777777777777777D0*(Cpl93*tmp116) - 
     -      2.2222222222222223D0*(Cpl94*tmp116) + 
     -      0.6666666666666666D0*(Cpl95*tmp116)) + 
     -   LInt154*(1.7777777777777777D0*(Cpl93*tmp116) - 
     -      2.2222222222222223D0*(Cpl94*tmp116) + 
     -      0.6666666666666666D0*(Cpl95*tmp116)) + 
     -   1.7777777777777777D0*(Cpl92*LInt140*tmp116) + 
     -   1.7777777777777777D0*(Cpl92*LInt152*tmp116) - 
     -   3.5555555555555554D0*(Cpl92*LInt101*WProd59) + 
     -   10.666666666666666D0*(Cpl92*LInt104*WProd59) + 
     -   14.222222222222221D0*(Cpl92*LInt107*WProd59) + 
     -   1.7777777777777777D0*
     -    (Cpl92*LInt134*(WProd46 - WProd49)*WProd59) - 
     -   1.7777777777777777D0*(Cpl92*LInt137*WProd49*WProd59) - 
     -   3.5555555555555554D0*(Cpl92*LInt143*WProd49*WProd59) - 
     -   3.5555555555555554D0*(Cpl92*LInt149*WProd49*WProd59) - 
     -   3.5555555555555554D0*(Cpl92*LInt155*WProd49*WProd59) - 
     -   3.5555555555555554D0*
     -    (Cpl92*LInt116*(WProd49*WProd59 - WProd60 + WProd61)) - 
     -   3.5555555555555554D0*
     -    (Cpl92*LInt119*(WProd49*WProd59 - WProd60 + WProd61)) - 
     -   3.5555555555555554D0*
     -    (Cpl92*LInt131*(2*WProd49*WProd59 - WProd60 + WProd61))

#ifdef DEBUG
	print *, 'amp(1,1,4) =', amp(1,1,4)
#endif

      endif

      end


