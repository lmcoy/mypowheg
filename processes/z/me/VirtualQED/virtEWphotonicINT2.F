      subroutine virtEWphotonicINT2(amp,resnr)
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

        tmp120 = WProd63*WProd64 - WProd63*WProd65 + 
     -   2*WProd63*WProd67 - 2*WProd68 + WProd69 - WProd70
        tmp121 = 2*WProd63*WProd65 - 2*WProd63*WProd67 + 
     -   2*WProd68 - 2*WProd69 + WProd70
        tmp122 = WProd63*WProd64 - WProd63*WProd65 - 
     -   3*WProd63*WProd67 + WProd68 - WProd69 + WProd70
        tmp123 = WProd63*WProd65 - 2*WProd63*WProd67 + 
     -   2*WProd68 - WProd69 + WProd70
        tmp124 = 2*WProd63*WProd64 - 2*WProd63*WProd65 - 
     -   6*WProd63*WProd67 + 4*WProd68 - 3*WProd69 + 2*WProd70
        tmp125 = 2*WProd63*WProd64 - 2*WProd63*WProd65 - 
     -   4*WProd63*WProd67 + 2*WProd68 - WProd69 + 2*WProd70
        amp(2,1,1) = -(0.4444444444444444D0*
     -      (Cpl92*LInt170*tmp120)) - 
     -   0.4444444444444444D0*(Cpl93*LInt171*tmp120) - 
     -   0.4444444444444444D0*(Cpl93*LInt172*tmp120) + 
     -   0.4444444444444444D0*(Cpl92*LInt203*tmp121) + 
     -   0.4444444444444444D0*(Cpl93*LInt204*tmp121) + 
     -   0.4444444444444444D0*(Cpl93*LInt205*tmp121) + 
     -   0.8888888888888888D0*(Cpl92*LInt179*tmp122) + 
     -   0.8888888888888888D0*(Cpl93*LInt180*tmp122) + 
     -   0.8888888888888888D0*(Cpl93*LInt181*tmp122) + 
     -   0.8888888888888888D0*(Cpl92*LInt182*tmp122) + 
     -   0.8888888888888888D0*(Cpl93*LInt183*tmp122) + 
     -   0.8888888888888888D0*(Cpl93*LInt184*tmp122) + 
     -   0.4444444444444444D0*(Cpl92*LInt167*tmp123) + 
     -   0.4444444444444444D0*(Cpl93*LInt168*tmp123) + 
     -   0.4444444444444444D0*(Cpl93*LInt169*tmp123) + 
     -   0.4444444444444444D0*(Cpl92*LInt185*tmp124) + 
     -   0.4444444444444444D0*(Cpl93*LInt186*tmp124) + 
     -   0.4444444444444444D0*(Cpl93*LInt187*tmp124) + 
     -   0.4444444444444444D0*(Cpl92*LInt197*tmp125) + 
     -   0.4444444444444444D0*(Cpl93*LInt198*tmp125) + 
     -   0.4444444444444444D0*(Cpl93*LInt199*tmp125) + 
     -   0.4444444444444444D0*(Cpl92*LInt209*tmp125) + 
     -   0.4444444444444444D0*(Cpl93*LInt210*tmp125) + 
     -   0.4444444444444444D0*(Cpl93*LInt211*tmp125) - 
     -   0.8888888888888888D0*(Cpl92*LInt158*WProd63) - 
     -   0.8888888888888888D0*(Cpl93*LInt159*WProd63) - 
     -   0.8888888888888888D0*(Cpl93*LInt160*WProd63) + 
     -   2.6666666666666665D0*(Cpl92*LInt161*WProd63) + 
     -   2.6666666666666665D0*(Cpl93*LInt162*WProd63) + 
     -   2.6666666666666665D0*(Cpl93*LInt163*WProd63) + 
     -   3.5555555555555554D0*(Cpl92*LInt164*WProd63) + 
     -   3.5555555555555554D0*(Cpl93*LInt165*WProd63) + 
     -   3.5555555555555554D0*(Cpl93*LInt166*WProd63) + 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt191*WProd63*(WProd64 - WProd67)) + 
     -   0.4444444444444444D0*
     -    (Cpl93*LInt192*WProd63*(WProd64 - WProd67)) + 
     -   0.4444444444444444D0*
     -    (Cpl93*LInt193*WProd63*(WProd64 - WProd67)) - 
     -   0.4444444444444444D0*(Cpl92*LInt194*WProd63*WProd67) - 
     -   0.4444444444444444D0*(Cpl93*LInt195*WProd63*WProd67) - 
     -   0.4444444444444444D0*(Cpl93*LInt196*WProd63*WProd67) - 
     -   0.8888888888888888D0*(Cpl92*LInt200*WProd63*WProd67) - 
     -   0.8888888888888888D0*(Cpl93*LInt201*WProd63*WProd67) - 
     -   0.8888888888888888D0*(Cpl93*LInt202*WProd63*WProd67) - 
     -   0.8888888888888888D0*(Cpl92*LInt206*WProd63*WProd67) - 
     -   0.8888888888888888D0*(Cpl93*LInt207*WProd63*WProd67) - 
     -   0.8888888888888888D0*(Cpl93*LInt208*WProd63*WProd67) - 
     -   0.8888888888888888D0*(Cpl92*LInt212*WProd63*WProd67) - 
     -   0.8888888888888888D0*(Cpl93*LInt213*WProd63*WProd67) - 
     -   0.8888888888888888D0*(Cpl93*LInt214*WProd63*WProd67) - 
     -   0.8888888888888888D0*
     -    (Cpl92*LInt173*(WProd63*WProd67 - WProd68 + WProd69)) - 
     -   0.8888888888888888D0*
     -    (Cpl93*LInt174*(WProd63*WProd67 - WProd68 + WProd69)) - 
     -   0.8888888888888888D0*
     -    (Cpl93*LInt175*(WProd63*WProd67 - WProd68 + WProd69)) - 
     -   0.8888888888888888D0*
     -    (Cpl92*LInt176*(WProd63*WProd67 - WProd68 + WProd69)) - 
     -   0.8888888888888888D0*
     -    (Cpl93*LInt177*(WProd63*WProd67 - WProd68 + WProd69)) - 
     -   0.8888888888888888D0*
     -    (Cpl93*LInt178*(WProd63*WProd67 - WProd68 + WProd69)) - 
     -   0.8888888888888888D0*
     -    (Cpl92*LInt188*(2*WProd63*WProd67 - WProd68 + WProd69))-
     -     0.8888888888888888D0*
     -    (Cpl93*LInt189*(2*WProd63*WProd67 - WProd68 + WProd69))-
     -     0.8888888888888888D0*
     -    (Cpl93*LInt190*(2*WProd63*WProd67 - WProd68 + WProd69))

#ifdef DEBUG
	print *, 'amp(2,1,1) =', amp(2,1,1)
c	print *, Cpl92, Lint170, tmp120
#endif

        tmp126 = -(0.4444444444444444D0*(Cpl93*WProd67*WProd71)) + 
     -   0.2222222222222222D0*(Cpl94*WProd67*WProd71)
        tmp127 = WProd64*WProd71 - WProd65*WProd71 + 
     -   2*WProd67*WProd71 + WProd72 - WProd73
        tmp128 = WProd64*WProd71 - WProd65*WProd71 - 
     -   2*WProd67*WProd71 + WProd73
        tmp129 = WProd64*WProd71 - WProd65*WProd71 - 
     -   3*WProd67*WProd71 - WProd72 + WProd73
	tmp130 = WProd65*WProd71 - 2*WProd67*WProd71 - WProd72 + WProd73
        tmp131 = 2*WProd65*WProd71 - 2*WProd67*WProd71 - 
     -   WProd72 + WProd73
        tmp132 = 0.4444444444444444D0*
     -    (Cpl93*(WProd64 - WProd67)*WProd71) - 
     -   0.2222222222222222D0*(Cpl94*(WProd64 - WProd67)*WProd71)
        tmp133 = -(0.8888888888888888D0*
     -      (Cpl93*(WProd67*WProd71 + WProd72))) + 
     -   0.4444444444444444D0*(Cpl94*(WProd67*WProd71 + WProd72))
        tmp134 = -(0.4444444444444444D0*
     -      (Cpl93*(2*WProd67*WProd71 + WProd72))) + 
     -   0.2222222222222222D0*(Cpl94*(2*WProd67*WProd71 + WProd72))
        tmp135 = -(0.4444444444444444D0*
     -      (Cpl93*(4*WProd67*WProd71 + 3*WProd72))) + 
     -   0.2222222222222222D0*(Cpl94*(4*WProd67*WProd71 + 3*WProd72))
        amp(2,1,2) = LInt195*tmp126 + LInt196*tmp126 + 
     -   LInt192*tmp132 + LInt193*tmp132 + LInt174*tmp133 + 
     -   LInt175*tmp133 + LInt177*tmp133 + LInt178*tmp133 + 
     -   LInt201*tmp134 + LInt202*tmp134 + LInt207*tmp134 + 
     -   LInt208*tmp134 + LInt213*tmp134 + LInt214*tmp134 + 
     -   LInt189*tmp135 + LInt190*tmp135 + 
     -   LInt171*(-(0.4444444444444444D0*(Cpl93*tmp127)) + 
     -      0.2222222222222222D0*(Cpl94*tmp127)) + 
     -   LInt172*(-(0.4444444444444444D0*(Cpl93*tmp127)) + 
     -      0.2222222222222222D0*(Cpl94*tmp127)) - 
     -   0.4444444444444444D0*(Cpl92*LInt170*tmp127) + 
     -   LInt198*(0.8888888888888888D0*(Cpl93*tmp128) - 
     -      0.4444444444444444D0*(Cpl94*tmp128)) + 
     -   LInt199*(0.8888888888888888D0*(Cpl93*tmp128) - 
     -      0.4444444444444444D0*(Cpl94*tmp128)) + 
     -   LInt210*(0.8888888888888888D0*(Cpl93*tmp128) - 
     -      0.4444444444444444D0*(Cpl94*tmp128)) + 
     -   LInt211*(0.8888888888888888D0*(Cpl93*tmp128) - 
     -      0.4444444444444444D0*(Cpl94*tmp128)) + 
     -   0.8888888888888888D0*(Cpl92*LInt197*tmp128) + 
     -   0.8888888888888888D0*(Cpl92*LInt209*tmp128) + 
     -   LInt180*(0.8888888888888888D0*(Cpl93*tmp129) - 
     -      0.4444444444444444D0*(Cpl94*tmp129)) + 
     -   LInt181*(0.8888888888888888D0*(Cpl93*tmp129) - 
     -      0.4444444444444444D0*(Cpl94*tmp129)) + 
     -   LInt183*(0.8888888888888888D0*(Cpl93*tmp129) - 
     -      0.4444444444444444D0*(Cpl94*tmp129)) + 
     -   LInt184*(0.8888888888888888D0*(Cpl93*tmp129) - 
     -      0.4444444444444444D0*(Cpl94*tmp129)) + 
     -   LInt186*(0.8888888888888888D0*(Cpl93*tmp129) - 
     -      0.4444444444444444D0*(Cpl94*tmp129)) + 
     -   LInt187*(0.8888888888888888D0*(Cpl93*tmp129) - 
     -      0.4444444444444444D0*(Cpl94*tmp129)) + 
     -   0.8888888888888888D0*(Cpl92*LInt179*tmp129) + 
     -   0.8888888888888888D0*(Cpl92*LInt182*tmp129) + 
     -   0.8888888888888888D0*(Cpl92*LInt185*tmp129) + 
     -   LInt168*(0.4444444444444444D0*(Cpl93*tmp130) - 
     -      0.2222222222222222D0*(Cpl94*tmp130)) + 
     -   LInt169*(0.4444444444444444D0*(Cpl93*tmp130) - 
     -      0.2222222222222222D0*(Cpl94*tmp130)) + 
     -   0.4444444444444444D0*(Cpl92*LInt167*tmp130) + 
     -   LInt204*(0.4444444444444444D0*(Cpl93*tmp131) - 
     -      0.2222222222222222D0*(Cpl94*tmp131)) + 
     -   LInt205*(0.4444444444444444D0*(Cpl93*tmp131) - 
     -      0.2222222222222222D0*(Cpl94*tmp131)) + 
     -   0.4444444444444444D0*(Cpl92*LInt203*tmp131) + 
     -   LInt165*(0.8888888888888888D0*(Cpl93*WProd71) - 
     -      0.4444444444444444D0*(Cpl94*WProd71)) + 
     -   LInt166*(0.8888888888888888D0*(Cpl93*WProd71) - 
     -      0.4444444444444444D0*(Cpl94*WProd71)) + 
     -   LInt159*(-(0.8888888888888888D0*(Cpl93*WProd71)) + 
     -      0.4444444444444444D0*(Cpl94*WProd71)) + 
     -   LInt160*(-(0.8888888888888888D0*(Cpl93*WProd71)) + 
     -      0.4444444444444444D0*(Cpl94*WProd71)) - 
     -   0.8888888888888888D0*(Cpl92*LInt158*WProd71) + 
     -   0.8888888888888888D0*(Cpl92*LInt164*WProd71) + 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt191*(WProd64 - WProd67)*WProd71) - 
     -   0.4444444444444444D0*(Cpl92*LInt194*WProd67*WProd71) - 
     -   0.8888888888888888D0*
     -    (Cpl92*LInt173*(WProd67*WProd71 + WProd72)) - 
     -   0.8888888888888888D0*
     -    (Cpl92*LInt176*(WProd67*WProd71 + WProd72)) - 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt200*(2*WProd67*WProd71 + WProd72)) - 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt206*(2*WProd67*WProd71 + WProd72)) - 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt212*(2*WProd67*WProd71 + WProd72)) - 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt188*(4*WProd67*WProd71 + 3*WProd72))

#ifdef DEBUG
	print *, 'amp(2,1,2) =', amp(2,1,2)
#endif

        tmp136 = -(0.4444444444444444D0*(Cpl93*WProd67*WProd74)) + 
     -   0.6666666666666666D0*(Cpl94*WProd67*WProd74)
        tmp137 = -2*WProd65*WProd74 + 2*WProd67*WProd74 + 
     -   WProd75 - WProd76
        tmp138 = -(WProd65*WProd74) + 2*WProd67*WProd74 + 
     -   WProd75 - WProd76
        tmp139 = WProd64*WProd74 - WProd65*WProd74 + 
     -   2*WProd67*WProd74 + WProd75 - WProd76
        tmp140 = -(WProd64*WProd74) + WProd65*WProd74 + 
     -   3*WProd67*WProd74 + WProd75 - WProd76
        tmp141 = WProd64*WProd74 - WProd65*WProd74 - 
     -   2*WProd67*WProd74 + WProd76
        tmp142 = 0.4444444444444444D0*
     -    (Cpl93*(WProd64 - WProd67)*WProd74) - 
     -   0.6666666666666666D0*(Cpl94*(WProd64 - WProd67)*WProd74)
        tmp143 = -(0.8888888888888888D0*
     -      (Cpl93*(WProd67*WProd74 + WProd75))) + 
     -   1.3333333333333333D0*(Cpl94*(WProd67*WProd74 + WProd75))
        tmp144 = -(0.4444444444444444D0*
     -      (Cpl93*(2*WProd67*WProd74 + WProd75))) + 
     -   0.6666666666666666D0*(Cpl94*(2*WProd67*WProd74 + WProd75))
        tmp145 = -(0.4444444444444444D0*
     -      (Cpl93*(4*WProd67*WProd74 + 3*WProd75))) + 
     -   0.6666666666666666D0*(Cpl94*(4*WProd67*WProd74 + 3*WProd75))
        amp(2,1,3) = LInt195*tmp136 + LInt196*tmp136 + 
     -   LInt192*tmp142 + LInt193*tmp142 + LInt174*tmp143 + 
     -   LInt175*tmp143 + LInt177*tmp143 + LInt178*tmp143 + 
     -   LInt201*tmp144 + LInt202*tmp144 + LInt207*tmp144 + 
     -   LInt208*tmp144 + LInt213*tmp144 + LInt214*tmp144 + 
     -   LInt189*tmp145 + LInt190*tmp145 + 
     -   LInt204*(-(0.4444444444444444D0*(Cpl93*tmp137)) + 
     -      0.6666666666666666D0*(Cpl94*tmp137)) + 
     -   LInt205*(-(0.4444444444444444D0*(Cpl93*tmp137)) + 
     -      0.6666666666666666D0*(Cpl94*tmp137)) - 
     -   0.4444444444444444D0*(Cpl92*LInt203*tmp137) + 
     -   LInt168*(-(0.4444444444444444D0*(Cpl93*tmp138)) + 
     -      0.6666666666666666D0*(Cpl94*tmp138)) + 
     -   LInt169*(-(0.4444444444444444D0*(Cpl93*tmp138)) + 
     -      0.6666666666666666D0*(Cpl94*tmp138)) - 
     -   0.4444444444444444D0*(Cpl92*LInt167*tmp138) + 
     -   LInt171*(-(0.4444444444444444D0*(Cpl93*tmp139)) + 
     -      0.6666666666666666D0*(Cpl94*tmp139)) + 
     -   LInt172*(-(0.4444444444444444D0*(Cpl93*tmp139)) + 
     -      0.6666666666666666D0*(Cpl94*tmp139)) - 
     -   0.4444444444444444D0*(Cpl92*LInt170*tmp139) + 
     -   LInt180*(-(0.8888888888888888D0*(Cpl93*tmp140)) + 
     -      1.3333333333333333D0*(Cpl94*tmp140)) + 
     -   LInt181*(-(0.8888888888888888D0*(Cpl93*tmp140)) + 
     -      1.3333333333333333D0*(Cpl94*tmp140)) + 
     -   LInt183*(-(0.8888888888888888D0*(Cpl93*tmp140)) + 
     -      1.3333333333333333D0*(Cpl94*tmp140)) + 
     -   LInt184*(-(0.8888888888888888D0*(Cpl93*tmp140)) + 
     -      1.3333333333333333D0*(Cpl94*tmp140)) + 
     -   LInt186*(-(0.8888888888888888D0*(Cpl93*tmp140)) + 
     -      1.3333333333333333D0*(Cpl94*tmp140)) + 
     -   LInt187*(-(0.8888888888888888D0*(Cpl93*tmp140)) + 
     -      1.3333333333333333D0*(Cpl94*tmp140)) - 
     -   0.8888888888888888D0*(Cpl92*LInt179*tmp140) - 
     -   0.8888888888888888D0*(Cpl92*LInt182*tmp140) - 
     -   0.8888888888888888D0*(Cpl92*LInt185*tmp140) + 
     -   LInt198*(0.8888888888888888D0*(Cpl93*tmp141) - 
     -      1.3333333333333333D0*(Cpl94*tmp141)) + 
     -   LInt199*(0.8888888888888888D0*(Cpl93*tmp141) - 
     -      1.3333333333333333D0*(Cpl94*tmp141)) + 
     -   LInt210*(0.8888888888888888D0*(Cpl93*tmp141) - 
     -      1.3333333333333333D0*(Cpl94*tmp141)) + 
     -   LInt211*(0.8888888888888888D0*(Cpl93*tmp141) - 
     -      1.3333333333333333D0*(Cpl94*tmp141)) + 
     -   0.8888888888888888D0*(Cpl92*LInt197*tmp141) + 
     -   0.8888888888888888D0*(Cpl92*LInt209*tmp141) + 
     -   LInt165*(0.8888888888888888D0*(Cpl93*WProd74) - 
     -      1.3333333333333333D0*(Cpl94*WProd74)) + 
     -   LInt166*(0.8888888888888888D0*(Cpl93*WProd74) - 
     -      1.3333333333333333D0*(Cpl94*WProd74)) + 
     -   LInt159*(-(0.8888888888888888D0*(Cpl93*WProd74)) + 
     -      1.3333333333333333D0*(Cpl94*WProd74)) + 
     -   LInt160*(-(0.8888888888888888D0*(Cpl93*WProd74)) + 
     -      1.3333333333333333D0*(Cpl94*WProd74)) - 
     -   0.8888888888888888D0*(Cpl92*LInt158*WProd74) + 
     -   0.8888888888888888D0*(Cpl92*LInt164*WProd74) + 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt191*(WProd64 - WProd67)*WProd74) - 
     -   0.4444444444444444D0*(Cpl92*LInt194*WProd67*WProd74) - 
     -   0.8888888888888888D0*
     -    (Cpl92*LInt173*(WProd67*WProd74 + WProd75)) - 
     -   0.8888888888888888D0*
     -    (Cpl92*LInt176*(WProd67*WProd74 + WProd75)) - 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt200*(2*WProd67*WProd74 + WProd75)) - 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt206*(2*WProd67*WProd74 + WProd75)) - 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt212*(2*WProd67*WProd74 + WProd75)) - 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt188*(4*WProd67*WProd74 + 3*WProd75))

#ifdef DEBUG
	print *, 'amp(2,1,3) =', amp(2,1,3)
#endif

        tmp146 = -(0.8888888888888888D0*(Cpl93*WProd77)) + 
     -   1.7777777777777777D0*(Cpl94*WProd77) - 
     -   0.6666666666666666D0*(Cpl95*WProd77)
        tmp147 = 2*Cpl95*WProd77 + 
     -   2.6666666666666665D0*(Cpl93*WProd77) - 
     -   5.333333333333333D0*(Cpl94*WProd77)
        tmp148 = 3.5555555555555554D0*(Cpl93*WProd77) - 
     -   7.111111111111111D0*(Cpl94*WProd77) + 
     -   2.6666666666666665D0*(Cpl95*WProd77)
        tmp149 = -(0.8888888888888888D0*(Cpl93*WProd67*WProd77)) + 
     -   1.7777777777777777D0*(Cpl94*WProd67*WProd77) - 
     -   0.6666666666666666D0*(Cpl95*WProd67*WProd77)
        tmp150 = -(0.4444444444444444D0*(Cpl93*WProd67*WProd77)) + 
     -   0.8888888888888888D0*(Cpl94*WProd67*WProd77) - 
     -   0.3333333333333333D0*(Cpl95*WProd67*WProd77)
        tmp151 = WProd64*WProd77 - WProd65*WProd77 + 
     -   2*WProd67*WProd77 - 2*WProd78 + WProd79 - WProd80
        tmp152 = 2*WProd65*WProd77 - 2*WProd67*WProd77 + 
     -   2*WProd78 - 2*WProd79 + WProd80
        tmp153 = WProd64*WProd77 - WProd65*WProd77 - 
     -   3*WProd67*WProd77 + WProd78 - WProd79 + WProd80
        tmp154 = WProd65*WProd77 - 2*WProd67*WProd77 + 
     -   2*WProd78 - WProd79 + WProd80
        tmp155 = -(WProd64*WProd77) + WProd65*WProd77 - 
     -   2*WProd67*WProd77 + 2*WProd78 - WProd79 + WProd80
        tmp156 = 2*WProd64*WProd77 - 2*WProd65*WProd77 - 
     -   6*WProd67*WProd77 + 4*WProd78 - 3*WProd79 + 2*WProd80
        tmp157 = 2*WProd64*WProd77 - 2*WProd65*WProd77 - 
     -   4*WProd67*WProd77 + 2*WProd78 - WProd79 + 2*WProd80
        tmp158 = 0.4444444444444444D0*
     -    (Cpl93*(WProd64 - WProd67)*WProd77) - 
     -   0.8888888888888888D0*(Cpl94*(WProd64 - WProd67)*WProd77) + 
     -   0.3333333333333333D0*(Cpl95*(WProd64 - WProd67)*WProd77)
        tmp159 = -(0.8888888888888888D0*
     -      (Cpl93*(WProd67*WProd77 - WProd78 + WProd79))) + 
     -   1.7777777777777777D0*
     -    (Cpl94*(WProd67*WProd77 - WProd78 + WProd79)) - 
     -   0.6666666666666666D0*
     -    (Cpl95*(WProd67*WProd77 - WProd78 + WProd79))
        tmp160 = -(0.8888888888888888D0*
     -      (Cpl93*(2*WProd67*WProd77 - WProd78 + WProd79))) + 
     -   1.7777777777777777D0*
     -    (Cpl94*(2*WProd67*WProd77 - WProd78 + WProd79)) - 
     -   0.6666666666666666D0*
     -    (Cpl95*(2*WProd67*WProd77 - WProd78 + WProd79))
        amp(2,1,4) = LInt159*tmp146 + LInt160*tmp146 + 
     -   LInt162*tmp147 + LInt163*tmp147 + LInt165*tmp148 + 
     -   LInt166*tmp148 + LInt201*tmp149 + LInt202*tmp149 + 
     -   LInt207*tmp149 + LInt208*tmp149 + LInt213*tmp149 + 
     -   LInt214*tmp149 + LInt195*tmp150 + LInt196*tmp150 + 
     -   LInt192*tmp158 + LInt193*tmp158 + LInt174*tmp159 + 
     -   LInt175*tmp159 + LInt177*tmp159 + LInt178*tmp159 + 
     -   LInt189*tmp160 + LInt190*tmp160 - 
     -   0.4444444444444444D0*(Cpl92*LInt170*tmp151) + 
     -   LInt204*(0.4444444444444444D0*(Cpl93*tmp152) - 
     -      0.8888888888888888D0*(Cpl94*tmp152) + 
     -      0.3333333333333333D0*(Cpl95*tmp152)) + 
     -   LInt205*(0.4444444444444444D0*(Cpl93*tmp152) - 
     -      0.8888888888888888D0*(Cpl94*tmp152) + 
     -      0.3333333333333333D0*(Cpl95*tmp152)) + 
     -   0.4444444444444444D0*(Cpl92*LInt203*tmp152) + 
     -   LInt180*(0.8888888888888888D0*(Cpl93*tmp153) - 
     -      1.7777777777777777D0*(Cpl94*tmp153) + 
     -      0.6666666666666666D0*(Cpl95*tmp153)) + 
     -   LInt181*(0.8888888888888888D0*(Cpl93*tmp153) - 
     -      1.7777777777777777D0*(Cpl94*tmp153) + 
     -      0.6666666666666666D0*(Cpl95*tmp153)) + 
     -   LInt183*(0.8888888888888888D0*(Cpl93*tmp153) - 
     -      1.7777777777777777D0*(Cpl94*tmp153) + 
     -      0.6666666666666666D0*(Cpl95*tmp153)) + 
     -   LInt184*(0.8888888888888888D0*(Cpl93*tmp153) - 
     -      1.7777777777777777D0*(Cpl94*tmp153) + 
     -      0.6666666666666666D0*(Cpl95*tmp153)) + 
     -   0.8888888888888888D0*(Cpl92*LInt179*tmp153) + 
     -   0.8888888888888888D0*(Cpl92*LInt182*tmp153) + 
     -   LInt168*(0.4444444444444444D0*(Cpl93*tmp154) - 
     -      0.8888888888888888D0*(Cpl94*tmp154) + 
     -      0.3333333333333333D0*(Cpl95*tmp154)) + 
     -   LInt169*(0.4444444444444444D0*(Cpl93*tmp154) - 
     -      0.8888888888888888D0*(Cpl94*tmp154) + 
     -      0.3333333333333333D0*(Cpl95*tmp154)) + 
     -   0.4444444444444444D0*(Cpl92*LInt167*tmp154) + 
     -   LInt171*(-(0.4444444444444444D0*(Cpl93*tmp151)) + 
     -      0.8888888888888888D0*(Cpl94*tmp151) + 
     -      0.3333333333333333D0*(Cpl95*tmp155)) + 
     -   LInt172*(-(0.4444444444444444D0*(Cpl93*tmp151)) + 
     -      0.8888888888888888D0*(Cpl94*tmp151) + 
     -      0.3333333333333333D0*(Cpl95*tmp155)) + 
     -   LInt186*(0.4444444444444444D0*(Cpl93*tmp156) - 
     -      0.8888888888888888D0*(Cpl94*tmp156) + 
     -      0.3333333333333333D0*(Cpl95*tmp156)) + 
     -   LInt187*(0.4444444444444444D0*(Cpl93*tmp156) - 
     -      0.8888888888888888D0*(Cpl94*tmp156) + 
     -      0.3333333333333333D0*(Cpl95*tmp156)) + 
     -   0.4444444444444444D0*(Cpl92*LInt185*tmp156) + 
     -   LInt198*(0.4444444444444444D0*(Cpl93*tmp157) - 
     -      0.8888888888888888D0*(Cpl94*tmp157) + 
     -      0.3333333333333333D0*(Cpl95*tmp157)) + 
     -   LInt199*(0.4444444444444444D0*(Cpl93*tmp157) - 
     -      0.8888888888888888D0*(Cpl94*tmp157) + 
     -      0.3333333333333333D0*(Cpl95*tmp157)) + 
     -   LInt210*(0.4444444444444444D0*(Cpl93*tmp157) - 
     -      0.8888888888888888D0*(Cpl94*tmp157) + 
     -      0.3333333333333333D0*(Cpl95*tmp157)) + 
     -   LInt211*(0.4444444444444444D0*(Cpl93*tmp157) - 
     -      0.8888888888888888D0*(Cpl94*tmp157) + 
     -      0.3333333333333333D0*(Cpl95*tmp157)) + 
     -   0.4444444444444444D0*(Cpl92*LInt197*tmp157) + 
     -   0.4444444444444444D0*(Cpl92*LInt209*tmp157) - 
     -   0.8888888888888888D0*(Cpl92*LInt158*WProd77) + 
     -   2.6666666666666665D0*(Cpl92*LInt161*WProd77) + 
     -   3.5555555555555554D0*(Cpl92*LInt164*WProd77) + 
     -   0.4444444444444444D0*
     -    (Cpl92*LInt191*(WProd64 - WProd67)*WProd77) - 
     -   0.4444444444444444D0*(Cpl92*LInt194*WProd67*WProd77) - 
     -   0.8888888888888888D0*(Cpl92*LInt200*WProd67*WProd77) - 
     -   0.8888888888888888D0*(Cpl92*LInt206*WProd67*WProd77) - 
     -   0.8888888888888888D0*(Cpl92*LInt212*WProd67*WProd77) - 
     -   0.8888888888888888D0*
     -    (Cpl92*LInt173*(WProd67*WProd77 - WProd78 + WProd79)) - 
     -   0.8888888888888888D0*
     -    (Cpl92*LInt176*(WProd67*WProd77 - WProd78 + WProd79)) - 
     -   0.8888888888888888D0*
     -    (Cpl92*LInt188*(2*WProd67*WProd77 - WProd78 + WProd79))

#ifdef DEBUG
	print *, 'amp(2,1,4) =', amp(2,1,4)
#endif

      endif

      end


