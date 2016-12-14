      subroutine virtEWnonphotonic1(amp,resnr)
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

        tmp7 = WProd81*WProd82 - WProd81*WProd83 - 
     -   3*WProd81*WProd85 + WProd86 - WProd87 + WProd88
        tmp8 = 2*WProd81*WProd82 - 2*WProd81*WProd83 - 
     -   4*WProd81*WProd85 + 2*WProd86 - WProd87 + 2*WProd88
        amp(1,1,1) = 3.5555555555555554D0*(Cpl96*LInt222*tmp7) + 
     -   3.5555555555555554D0*(Cpl96*LInt223*tmp7) + 
     -   1.7777777777777777D0*(Cpl96*LInt228*tmp8) + 
     -   1.7777777777777777D0*(Cpl96*LInt232*tmp8) - 
     -   3.5555555555555554D0*(Cpl96*LInt215*WProd81) + 
     -   10.666666666666666D0*(Cpl96*LInt216*WProd81) + 
     -   14.222222222222221D0*(Cpl96*LInt217*WProd81) + 
     -   1.7777777777777777D0*
     -    (Cpl96*LInt226*WProd81*(WProd82 - WProd85)) - 
     -   1.7777777777777777D0*(Cpl96*LInt227*WProd81*WProd85) - 
     -   3.5555555555555554D0*(Cpl96*LInt229*WProd81*WProd85) - 
     -   3.5555555555555554D0*(Cpl96*LInt231*WProd81*WProd85) - 
     -   3.5555555555555554D0*(Cpl96*LInt233*WProd81*WProd85) + 
     -   Prp18**2*(LInt352*
     -       (-6*Cpl92*WProd81 + 
     -         1.3333333333333333D0*(Cpl93*WProd81) - 
     -         0.6666666666666666D0*(Cpl96*WProd81)) + 
     -      LInt236*(12*Cpl92*WProd81 - 
     -         2.6666666666666665D0*(Cpl93*WProd81) + 
     -         1.3333333333333333D0*(Cpl96*WProd81)) - 
     -      1.3333333333333333D0*(Cpl102*LInt250*WProd81) + 
     -      1.3333333333333333D0*(Cpl101*LInt251*WProd81) - 
     -      0.3333333333333333D0*(Cpl101*LInt353*WProd81) - 
     -      0.3333333333333333D0*(Cpl101*LInt354*WProd81) + 
     -      0.2222222222222222D0*(Cpl101*WProd81*WProd82) + 
     -      0.8888888888888888D0*(Cpl92*WProd81*WProd82) + 
     -      2.6666666666666665D0*(Cpl92*LInt237*WProd81*WProd82) + 
     -      1.3333333333333333D0*(Cpl101*LInt252*WProd81*WProd82) + 
     -      LInt235*(-(2.6666666666666665D0*(Cpl103*WProd81)) + 
     -         2.6666666666666665D0*(Cpl99*WProd81) + 
     -         6.666666666666667D0*
     -          (Cpl92*WProd81*(WProd83 + WProd85)))) + 
     -   Prp17**2*(16*Cpl92*LInt236*WProd81 - 
     -      8*Cpl92*LInt352*WProd81 + 
     -      0.8888888888888888D0*(Cpl92*WProd81*WProd82) + 
     -      2.6666666666666665D0*(Cpl92*LInt237*WProd81*WProd82) + 
     -      6.666666666666667D0*
     -       (Cpl92*LInt235*WProd81*(WProd83 + WProd85))) + 
     -   Prp17*(-(1.9259259259259258D0*(Cpl93*WProd81)) + 
     -      2.6666666666666665D0*(Cpl93*LInt238*WProd81) + 
     -      1.1851851851851851D0*(Cpl93*LInt239*WProd81) - 
     -      5.333333333333333D0*(Cpl93*LInt242*WProd81) - 
     -      2.3703703703703702D0*(Cpl93*LInt243*WProd81) + 
     -      LInt241*(1.1851851851851851D0*(Cpl98*WProd81) - 
     -         1.1851851851851851D0*(Cpl93*WProd81*WProd82)) - 
     -      1.1851851851851851D0*(Cpl93*LInt245*WProd81*WProd82) - 
     -      1.1851851851851851D0*(Cpl93*LInt247*WProd81*WProd82) + 
     -      2.6666666666666665D0*
     -       (Cpl93*LInt246*WProd81*
     -         (2*WProd82 - 3*WProd83 - 3*WProd85)) + 
     -      Prp18*(LInt352*
     -          (13.333333333333334D0*(Cpl92*WProd81) - 
     -            2.6666666666666665D0*(Cpl93*WProd81)) + 
     -         LInt236*(-(26.666666666666668D0*(Cpl92*WProd81)) + 
     -            5.333333333333333D0*(Cpl93*WProd81)) - 
     -         1.7777777777777777D0*(Cpl92*WProd81*WProd82) - 
     -         5.333333333333333D0*(Cpl92*LInt237*WProd81*WProd82) + 
     -         LInt235*(-(5.333333333333333D0*(Cpl100*WProd81)) - 
     -            5.333333333333333D0*(Cpl99*WProd81) - 
     -            13.333333333333334D0*
     -             (Cpl92*WProd81*(WProd83 + WProd85)))) + 
     -      LInt240*(2.6666666666666665D0*(Cpl98*WProd81) - 
     -         2.6666666666666665D0*
     -          (Cpl93*WProd81*(WProd83 + WProd85))) - 
     -      2.6666666666666665D0*
     -       (Cpl93*LInt244*WProd81*(WProd83 + WProd85))) + 
     -   Prp18*(-(1.9259259259259258D0*(Cpl96*WProd81)) + 
     -      2.6666666666666665D0*(Cpl96*LInt238*WProd81) + 
     -      1.1851851851851851D0*(Cpl96*LInt239*WProd81) - 
     -      5.333333333333333D0*(Cpl96*LInt242*WProd81) - 
     -      2.3703703703703702D0*(Cpl96*LInt243*WProd81) + 
     -      LInt241*(1.1851851851851851D0*(Cpl97*WProd81) - 
     -         1.1851851851851851D0*(Cpl96*WProd81*WProd82)) - 
     -      1.1851851851851851D0*(Cpl96*LInt245*WProd81*WProd82) - 
     -      1.1851851851851851D0*(Cpl96*LInt247*WProd81*WProd82) + 
     -      2.6666666666666665D0*
     -       (Cpl96*LInt246*WProd81*
     -         (2*WProd82 - 3*WProd83 - 3*WProd85)) + 
     -      LInt240*(2.6666666666666665D0*(Cpl97*WProd81) - 
     -         2.6666666666666665D0*
     -          (Cpl96*WProd81*(WProd83 + WProd85))) - 
     -      2.6666666666666665D0*
     -       (Cpl96*LInt244*WProd81*(WProd83 + WProd85))) - 
     -   3.5555555555555554D0*
     -    (Cpl96*LInt220*(WProd81*WProd85 - WProd86 + WProd87)) - 
     -   3.5555555555555554D0*
     -    (Cpl96*LInt221*(WProd81*WProd85 - WProd86 + WProd87)) - 
     -   3.5555555555555554D0*
     -    (Cpl96*LInt225*(2*WProd81*WProd85 - WProd86 + WProd87))-
     -     1.7777777777777777D0*
     -    (Cpl96*LInt219*(WProd81*WProd82 - WProd81*WProd83 + 
     -        2*WProd81*WProd85 - 2*WProd86 + WProd87 - WProd88))+
     -     1.7777777777777777D0*
     -    (Cpl96*LInt230*(2*WProd81*WProd83 - 2*WProd81*WProd85 + 
     -        2*WProd86 - 2*WProd87 + WProd88)) + 
     -   1.7777777777777777D0*
     -    (Cpl96*LInt218*(WProd81*WProd83 - 2*WProd81*WProd85 + 
     -        2*WProd86 - WProd87 + WProd88))
        amp(1,1,1) = amp(1,1,1) + 
     -   1.7777777777777777D0*
     -    (Cpl96*LInt224*(2*WProd81*WProd82 - 2*WProd81*WProd83 - 
     -        6*WProd81*WProd85 + 4*WProd86 - 3*WProd87 + 2*WProd88
     -        ))

	do Gen1 = 1, 3
        amp(1,1,1) = amp(1,1,1) + 
     -   Prp17**2*(1.5802469135802468D0*(Cpl92*WProd81*WProd82) - 
     -      4.7407407407407405D0*(WProd81*Cpl104(Gen1)) - 
     -      4.7407407407407405D0*
     -       (WProd81*Cpl104(Gen1)*LInt254(Gen1)) + 
     -      3.5555555555555554D0*
     -       (Cpl92*WProd81*WProd82*LInt255(Gen1)) + 
     -      1.1851851851851851D0*
     -       (Cpl92*WProd81*WProd82*LInt256(Gen1)) + 
     -      4.7407407407407405D0*
     -       (Cpl92*WProd81*WProd82*LInt257(Gen1)) + 
     -      4.7407407407407405D0*(WProd81*Cpl104(Gen1)*LInt355(Gen1))
     -      ) + Prp17*Prp18*
     -    (3.1604938271604937D0*(Cpl93*WProd81*WProd82) - 
     -      1.1851851851851851D0*(Cpl94*WProd81*WProd82) + 
     -      3.5555555555555554D0*(WProd81*Cpl105(Gen1)) - 
     -      9.481481481481481D0*(WProd81*Cpl106(Gen1)) + 
     -      (3.5555555555555554D0*(WProd81*Cpl105(Gen1)) - 
     -         9.481481481481481D0*(WProd81*Cpl106(Gen1)))*
     -       LInt254(Gen1) + 
     -      (7.111111111111111D0*(Cpl93*WProd81*WProd82) - 
     -         1.7777777777777777D0*(Cpl94*WProd81*WProd82))*
     -       LInt255(Gen1) + 
     -      (2.3703703703703702D0*(Cpl93*WProd81*WProd82) - 
     -         1.7777777777777777D0*(Cpl94*WProd81*WProd82))*
     -       LInt256(Gen1) + 
     -      (9.481481481481481D0*(Cpl93*WProd81*WProd82) - 
     -         3.5555555555555554D0*(Cpl94*WProd81*WProd82))*
     -       LInt257(Gen1) + 
     -      (-(3.5555555555555554D0*(WProd81*Cpl105(Gen1))) + 
     -         9.481481481481481D0*(WProd81*Cpl106(Gen1)))*
     -       LInt355(Gen1)) + 
     -   Prp18**2*(0.5185185185185185D0*(Cpl101*WProd81*WProd82) - 
     -      1.1851851851851851D0*(Cpl107*WProd81*WProd82) + 
     -      1.5802469135802468D0*(Cpl96*WProd81*WProd82) - 
     -      1.3333333333333333D0*(WProd81*Cpl108(Gen1)) + 
     -      3.5555555555555554D0*(WProd81*Cpl109(Gen1)) - 
     -      4.7407407407407405D0*(WProd81*Cpl110(Gen1)) + 
     -      (0.6666666666666666D0*(WProd81*Cpl108(Gen1)) + 
     -         3.5555555555555554D0*(WProd81*Cpl109(Gen1)) - 
     -         4.7407407407407405D0*(WProd81*Cpl110(Gen1)))*
     -       LInt254(Gen1) + 
     -      (0.4444444444444444D0*(Cpl101*WProd81*WProd82) - 
     -         1.7777777777777777D0*(Cpl107*WProd81*WProd82) + 
     -         3.5555555555555554D0*(Cpl96*WProd81*WProd82))*
     -       LInt255(Gen1) + 
     -      (1.3333333333333333D0*(Cpl101*WProd81*WProd82) - 
     -         1.7777777777777777D0*(Cpl107*WProd81*WProd82) + 
     -         1.1851851851851851D0*(Cpl96*WProd81*WProd82))*
     -       LInt256(Gen1) + 
     -      (1.3333333333333333D0*(Cpl101*WProd81*WProd82) - 
     -         3.5555555555555554D0*(Cpl107*WProd81*WProd82) + 
     -         4.7407407407407405D0*(Cpl96*WProd81*WProd82))*
     -       LInt257(Gen1) + 
     -      (1.3333333333333333D0*(WProd81*Cpl108(Gen1)) - 
     -         3.5555555555555554D0*(WProd81*Cpl109(Gen1)) + 
     -         4.7407407407407405D0*(WProd81*Cpl110(Gen1)))*
     -       LInt355(Gen1))
	enddo

#ifdef DEBUG
	print *, 'amp(1,1,1) =', amp(1,1,1)
#endif

        tmp9 = 0.5925925925925926D0*(Cpl107*WProd82*WProd89) - 
     -   1.1851851851851851D0*(Cpl96*WProd82*WProd89)
        tmp10 = WProd82*WProd89 - WProd83*WProd89 + 
     -   2*WProd85*WProd89 + WProd90 - WProd91
        tmp11 = WProd82*WProd89 - WProd83*WProd89 - 
     -   2*WProd85*WProd89 + WProd91
        tmp12 = WProd82*WProd89 - WProd83*WProd89 - 
     -   3*WProd85*WProd89 - WProd90 + WProd91
	tmp13 = WProd83*WProd89 - 2*WProd85*WProd89 - WProd90 + WProd91
        tmp14 = 2*WProd83*WProd89 - 2*WProd85*WProd89 - 
     -   WProd90 + WProd91
        tmp15 = -(0.8888888888888888D0*
     -      (Cpl101*(WProd85*WProd89 + WProd90))) + 
     -   3.5555555555555554D0*(Cpl107*(WProd85*WProd89 + WProd90)) - 
     -   3.5555555555555554D0*(Cpl96*(WProd85*WProd89 + WProd90))
        tmp16 = -(0.4444444444444444D0*
     -      (Cpl101*(2*WProd85*WProd89 + WProd90))) + 
     -   1.7777777777777777D0*
     -    (Cpl107*(2*WProd85*WProd89 + WProd90)) - 
     -   1.7777777777777777D0*(Cpl96*(2*WProd85*WProd89 + WProd90))
        amp(1,1,2) = LInt220*tmp15 + LInt221*tmp15 + 
     -   LInt229*tmp16 + LInt231*tmp16 + LInt233*tmp16 + 
     -   LInt219*(-(0.4444444444444444D0*(Cpl101*tmp10)) + 
     -      1.7777777777777777D0*(Cpl107*tmp10) - 
     -      1.7777777777777777D0*(Cpl96*tmp10)) + 
     -   LInt228*(0.8888888888888888D0*(Cpl101*tmp11) - 
     -      3.5555555555555554D0*(Cpl107*tmp11) + 
     -      3.5555555555555554D0*(Cpl96*tmp11)) + 
     -   LInt232*(0.8888888888888888D0*(Cpl101*tmp11) - 
     -      3.5555555555555554D0*(Cpl107*tmp11) + 
     -      3.5555555555555554D0*(Cpl96*tmp11)) + 
     -   LInt222*(0.8888888888888888D0*(Cpl101*tmp12) - 
     -      3.5555555555555554D0*(Cpl107*tmp12) + 
     -      3.5555555555555554D0*(Cpl96*tmp12)) + 
     -   LInt223*(0.8888888888888888D0*(Cpl101*tmp12) - 
     -      3.5555555555555554D0*(Cpl107*tmp12) + 
     -      3.5555555555555554D0*(Cpl96*tmp12)) + 
     -   LInt224*(0.8888888888888888D0*(Cpl101*tmp12) - 
     -      3.5555555555555554D0*(Cpl107*tmp12) + 
     -      3.5555555555555554D0*(Cpl96*tmp12)) + 
     -   LInt218*(0.4444444444444444D0*(Cpl101*tmp13) - 
     -      1.7777777777777777D0*(Cpl107*tmp13) + 
     -      1.7777777777777777D0*(Cpl96*tmp13)) + 
     -   LInt230*(0.4444444444444444D0*(Cpl101*tmp14) - 
     -      1.7777777777777777D0*(Cpl107*tmp14) + 
     -      1.7777777777777777D0*(Cpl96*tmp14)) + 
     -   LInt215*(-(0.8888888888888888D0*(Cpl101*WProd89)) + 
     -      3.5555555555555554D0*(Cpl107*WProd89) - 
     -      3.5555555555555554D0*(Cpl96*WProd89)) + 
     -   LInt217*(0.8888888888888888D0*(Cpl101*WProd89) - 
     -      3.5555555555555554D0*(Cpl107*WProd89) + 
     -      3.5555555555555554D0*(Cpl96*WProd89)) + 
     -   LInt226*(0.4444444444444444D0*
     -       (Cpl101*(WProd82 - WProd85)*WProd89) - 
     -      1.7777777777777777D0*
     -       (Cpl107*(WProd82 - WProd85)*WProd89) + 
     -      1.7777777777777777D0*(Cpl96*(WProd82 - WProd85)*WProd89))
     -     + LInt227*(-(0.4444444444444444D0*
     -         (Cpl101*WProd85*WProd89)) + 
     -      1.7777777777777777D0*(Cpl107*WProd85*WProd89) - 
     -      1.7777777777777777D0*(Cpl96*WProd85*WProd89)) + 
     -   Prp18**2*(LInt353*
     -       (-(0.3333333333333333D0*(Cpl101*WProd89)) + 
     -         0.16666666666666666D0*(Cpl112*WProd89)) + 
     -      LInt354*(-(0.3333333333333333D0*(Cpl101*WProd89)) + 
     -         0.16666666666666666D0*(Cpl112*WProd89)) + 
     -      LInt251*(1.3333333333333333D0*(Cpl101*WProd89) - 
     -         0.6666666666666666D0*(Cpl112*WProd89)) + 
     -      LInt250*(-(1.3333333333333333D0*(Cpl102*WProd89)) + 
     -         0.6666666666666666D0*(Cpl121*WProd89)) + 
     -      LInt352*(3*Cpl111*WProd89 - 6*Cpl92*WProd89 + 
     -         0.3333333333333333D0*(Cpl107*WProd89) + 
     -         1.3333333333333333D0*(Cpl93*WProd89) - 
     -         0.6666666666666666D0*(Cpl94*WProd89) - 
     -         0.6666666666666666D0*(Cpl96*WProd89)) + 
     -      LInt236*(-6*Cpl111*WProd89 + 12*Cpl92*WProd89 - 
     -         0.6666666666666666D0*(Cpl107*WProd89) - 
     -         2.6666666666666665D0*(Cpl93*WProd89) + 
     -         1.3333333333333333D0*(Cpl94*WProd89) + 
     -         1.3333333333333333D0*(Cpl96*WProd89)) + 
     -      0.2222222222222222D0*(Cpl101*WProd82*WProd89) - 
     -      0.4444444444444444D0*(Cpl111*WProd82*WProd89) - 
     -      0.1111111111111111D0*(Cpl112*WProd82*WProd89) + 
     -      LInt252*(1.3333333333333333D0*(Cpl101*WProd82*WProd89) - 
     -         0.6666666666666666D0*(Cpl112*WProd82*WProd89)) + 
     -      0.8888888888888888D0*(Cpl92*WProd82*WProd89) + 
     -      LInt237*(-(1.3333333333333333D0*
     -            (Cpl111*WProd82*WProd89)) + 
     -         2.6666666666666665D0*(Cpl92*WProd82*WProd89)) + 
     -      LInt235*(-(2.6666666666666665D0*(Cpl103*WProd89)) - 
     -         1.3333333333333333D0*(Cpl120*WProd89) + 
     -         1.3333333333333333D0*(Cpl122*WProd89) + 
     -         2.6666666666666665D0*(Cpl99*WProd89) - 
     -         3.3333333333333335D0*
     -          (Cpl111*(WProd83 + WProd85)*WProd89) + 
     -         6.666666666666667D0*
     -          (Cpl92*(WProd83 + WProd85)*WProd89))) + 
     -   Prp17**2*(16*Cpl92*LInt236*WProd89 - 
     -      8*Cpl92*LInt352*WProd89 + 
     -      0.8888888888888888D0*(Cpl92*WProd82*WProd89) + 
     -      2.6666666666666665D0*(Cpl92*LInt237*WProd82*WProd89) + 
     -      6.666666666666667D0*
     -       (Cpl92*LInt235*(WProd83 + WProd85)*WProd89))
        amp(1,1,2) = amp(1,1,2) + 
     -   Prp18*(LInt245*tmp9 + LInt247*tmp9 - Cpl101*WProd89 + 
     -      2.2962962962962963D0*(Cpl107*WProd89) + 
     -      0.6666666666666666D0*(Cpl111*WProd89) + 
     -      0.16666666666666666D0*(Cpl112*WProd89) - 
     -      0.3333333333333333D0*(Cpl95*WProd89) + 
     -      LInt239*(-(0.5925925925925926D0*(Cpl107*WProd89)) + 
     -         1.1851851851851851D0*(Cpl96*WProd89)) - 
     -      1.9259259259259258D0*(Cpl96*WProd89) + 
     -      LInt243*(1.1851851851851851D0*(Cpl107*WProd89) - 
     -         2.3703703703703702D0*(Cpl96*WProd89)) + 
     -      LInt238*(2*Cpl101*WProd89 - 4*Cpl107*WProd89 - 
     -         0.3333333333333333D0*(Cpl112*WProd89) + 
     -         2.6666666666666665D0*(Cpl96*WProd89)) + 
     -      LInt242*(-4*Cpl101*WProd89 + 8*Cpl107*WProd89 + 
     -         0.6666666666666666D0*(Cpl112*WProd89) - 
     -         5.333333333333333D0*(Cpl96*WProd89)) - 
     -      1.3333333333333333D0*(Cpl111*LInt235*WProd89) + 
     -      0.6666666666666666D0*(Cpl95*LInt258*WProd89) - 
     -      2.6666666666666665D0*(Cpl111*LInt260*WProd89) - 
     -      1.3333333333333333D0*(Cpl95*LInt261*WProd89) + 
     -      LInt241*(-(0.5925925925925926D0*(Cpl114*WProd89)) + 
     -         1.1851851851851851D0*(Cpl97*WProd89) + 
     -         0.5925925925925926D0*(Cpl107*WProd82*WProd89) - 
     -         1.1851851851851851D0*(Cpl96*WProd82*WProd89)) + 
     -      LInt246*(2*Cpl101*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -          WProd89 - 4*Cpl107*
     -          (2*WProd82 - 3*WProd83 - 3*WProd85)*WProd89 - 
     -         0.3333333333333333D0*
     -          (Cpl112*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -            WProd89) + 
     -         2.6666666666666665D0*
     -          (Cpl96*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd89)
     -         ) - 1.3333333333333333D0*
     -       (Cpl111*LInt264*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -         WProd89) + 0.6666666666666666D0*
     -       (Cpl95*LInt265*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -         WProd89) + LInt259*
     -       (0.6666666666666666D0*(Cpl113*WProd89) - 
     -         0.6666666666666666D0*
     -          (Cpl95*(WProd83 + WProd85)*WProd89)) + 
     -      LInt244*(-2*Cpl101*(WProd83 + WProd85)*WProd89 + 
     -         4*Cpl107*(WProd83 + WProd85)*WProd89 + 
     -         0.3333333333333333D0*
     -          (Cpl112*(WProd83 + WProd85)*WProd89) - 
     -         2.6666666666666665D0*
     -          (Cpl96*(WProd83 + WProd85)*WProd89)) + 
     -      LInt240*(-4*Cpl114*WProd89 + 2*Cpl115*WProd89 - 
     -         2*Cpl101*(WProd83 + WProd85)*WProd89 + 
     -         4*Cpl107*(WProd83 + WProd85)*WProd89 - 
     -         0.3333333333333333D0*(Cpl116*WProd89) + 
     -         2.6666666666666665D0*(Cpl97*WProd89) + 
     -         0.3333333333333333D0*
     -          (Cpl112*(WProd83 + WProd85)*WProd89) - 
     -         2.6666666666666665D0*
     -          (Cpl96*(WProd83 + WProd85)*WProd89)) + 
     -      1.3333333333333333D0*
     -       (Cpl111*LInt262*(WProd83 + WProd85)*WProd89) - 
     -      0.6666666666666666D0*
     -       (Cpl95*LInt263*(WProd83 + WProd85)*WProd89))
        amp(1,1,2) = amp(1,1,2) + 
     -   Prp17*(-(0.6666666666666666D0*(Cpl111*WProd89)) - 
     -      1.9259259259259258D0*(Cpl93*WProd89) + 
     -      1.3333333333333333D0*(Cpl94*WProd89) - 
     -      0.3333333333333333D0*(Cpl95*WProd89) + 
     -      LInt238*(2.6666666666666665D0*(Cpl93*WProd89) - 
     -         2.6666666666666665D0*(Cpl94*WProd89) + 
     -         0.6666666666666666D0*(Cpl95*WProd89)) + 
     -      LInt242*(-(5.333333333333333D0*(Cpl93*WProd89)) + 
     -         5.333333333333333D0*(Cpl94*WProd89) - 
     -         1.3333333333333333D0*(Cpl95*WProd89)) + 
     -      1.3333333333333333D0*(Cpl111*LInt235*WProd89) + 
     -      1.1851851851851851D0*(Cpl93*LInt239*WProd89) - 
     -      2.3703703703703702D0*(Cpl93*LInt243*WProd89) + 
     -      2.6666666666666665D0*(Cpl111*LInt260*WProd89) + 
     -      LInt241*(1.1851851851851851D0*(Cpl98*WProd89) - 
     -         1.1851851851851851D0*(Cpl93*WProd82*WProd89)) - 
     -      1.1851851851851851D0*(Cpl93*LInt245*WProd82*WProd89) - 
     -      1.1851851851851851D0*(Cpl93*LInt247*WProd82*WProd89) + 
     -      LInt246*(2.6666666666666665D0*
     -          (Cpl93*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd89)
     -           - 2.6666666666666665D0*
     -          (Cpl94*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd89)
     -           + 0.6666666666666666D0*
     -          (Cpl95*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd89)
     -         ) + 1.3333333333333333D0*
     -       (Cpl111*LInt264*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -         WProd89) + Prp18*
     -       (LInt352*(-(3.3333333333333335D0*(Cpl111*WProd89)) + 
     -            13.333333333333334D0*(Cpl92*WProd89) - 
     -            2.6666666666666665D0*(Cpl93*WProd89) + 
     -            0.6666666666666666D0*(Cpl94*WProd89)) + 
     -         LInt236*(6.666666666666667D0*(Cpl111*WProd89) - 
     -            26.666666666666668D0*(Cpl92*WProd89) + 
     -            5.333333333333333D0*(Cpl93*WProd89) - 
     -            1.3333333333333333D0*(Cpl94*WProd89)) + 
     -         0.4444444444444444D0*(Cpl111*WProd82*WProd89) - 
     -         1.7777777777777777D0*(Cpl92*WProd82*WProd89) + 
     -         LInt237*(1.3333333333333333D0*
     -             (Cpl111*WProd82*WProd89) - 
     -            5.333333333333333D0*(Cpl92*WProd82*WProd89)) + 
     -         LInt235*(-(5.333333333333333D0*(Cpl100*WProd89)) + 
     -            1.3333333333333333D0*(Cpl119*WProd89) + 
     -            1.3333333333333333D0*(Cpl120*WProd89) - 
     -            5.333333333333333D0*(Cpl99*WProd89) + 
     -            3.3333333333333335D0*
     -             (Cpl111*(WProd83 + WProd85)*WProd89) - 
     -            13.333333333333334D0*
     -             (Cpl92*(WProd83 + WProd85)*WProd89))) + 
     -      LInt244*(-(2.6666666666666665D0*
     -            (Cpl93*(WProd83 + WProd85)*WProd89)) + 
     -         2.6666666666666665D0*
     -          (Cpl94*(WProd83 + WProd85)*WProd89) - 
     -         0.6666666666666666D0*
     -          (Cpl95*(WProd83 + WProd85)*WProd89)) + 
     -      LInt240*(-(2.6666666666666665D0*(Cpl117*WProd89)) + 
     -         0.6666666666666666D0*(Cpl118*WProd89) + 
     -         2.6666666666666665D0*(Cpl98*WProd89) - 
     -         2.6666666666666665D0*
     -          (Cpl93*(WProd83 + WProd85)*WProd89) + 
     -         2.6666666666666665D0*
     -          (Cpl94*(WProd83 + WProd85)*WProd89) - 
     -         0.6666666666666666D0*
     -          (Cpl95*(WProd83 + WProd85)*WProd89)) - 
     -      1.3333333333333333D0*
     -       (Cpl111*LInt262*(WProd83 + WProd85)*WProd89)) + 
     -   LInt225*(-(0.4444444444444444D0*
     -         (Cpl101*(4*WProd85*WProd89 + 3*WProd90))) + 
     -      1.7777777777777777D0*
     -       (Cpl107*(4*WProd85*WProd89 + 3*WProd90)) - 
     -      1.7777777777777777D0*
     -       (Cpl96*(4*WProd85*WProd89 + 3*WProd90)))

	do Gen1 = 1, 3
        amp(1,1,2) = amp(1,1,2) + 
     -   Prp17**2*(1.5802469135802468D0*(Cpl92*WProd82*WProd89) - 
     -      4.7407407407407405D0*(WProd89*Cpl104(Gen1)) - 
     -      4.7407407407407405D0*
     -       (WProd89*Cpl104(Gen1)*LInt254(Gen1)) + 
     -      3.5555555555555554D0*
     -       (Cpl92*WProd82*WProd89*LInt255(Gen1)) + 
     -      1.1851851851851851D0*
     -       (Cpl92*WProd82*WProd89*LInt256(Gen1)) + 
     -      4.7407407407407405D0*
     -       (Cpl92*WProd82*WProd89*LInt257(Gen1)) + 
     -      4.7407407407407405D0*(WProd89*Cpl104(Gen1)*LInt355(Gen1))
     -      ) + Prp17*Prp18*
     -    (3.1604938271604937D0*(Cpl93*WProd82*WProd89) - 
     -      1.9753086419753085D0*(Cpl94*WProd82*WProd89) + 
     -      0.2962962962962963D0*(Cpl95*WProd82*WProd89) + 
     -      5.925925925925926D0*(WProd89*Cpl105(Gen1)) - 
     -      9.481481481481481D0*(WProd89*Cpl106(Gen1)) - 
     -      0.8888888888888888D0*(WProd89*Cpl123(Gen1)) + 
     -      (5.925925925925926D0*(WProd89*Cpl105(Gen1)) - 
     -         9.481481481481481D0*(WProd89*Cpl106(Gen1)) - 
     -         0.8888888888888888D0*(WProd89*Cpl123(Gen1)))*
     -       LInt254(Gen1) + 
     -      (7.111111111111111D0*(Cpl93*WProd82*WProd89) - 
     -         3.5555555555555554D0*(Cpl94*WProd82*WProd89) + 
     -         0.4444444444444444D0*(Cpl95*WProd82*WProd89))*
     -       LInt255(Gen1) + 
     -      (2.3703703703703702D0*(Cpl93*WProd82*WProd89) - 
     -         2.3703703703703702D0*(Cpl94*WProd82*WProd89) + 
     -         0.4444444444444444D0*(Cpl95*WProd82*WProd89))*
     -       LInt256(Gen1) + 
     -      (9.481481481481481D0*(Cpl93*WProd82*WProd89) - 
     -         5.925925925925926D0*(Cpl94*WProd82*WProd89) + 
     -         0.8888888888888888D0*(Cpl95*WProd82*WProd89))*
     -       LInt257(Gen1) + 
     -      (-(5.925925925925926D0*(WProd89*Cpl105(Gen1))) + 
     -         9.481481481481481D0*(WProd89*Cpl106(Gen1)) + 
     -         0.8888888888888888D0*(WProd89*Cpl123(Gen1)))*
     -       LInt355(Gen1)) + 
     -   Prp18**2*(1.1111111111111112D0*(Cpl101*WProd82*WProd89) - 
     -      1.9753086419753085D0*(Cpl107*WProd82*WProd89) - 
     -      0.25925925925925924D0*(Cpl112*WProd82*WProd89) + 
     -      1.5802469135802468D0*(Cpl96*WProd82*WProd89) - 
     -      3.111111111111111D0*(WProd89*Cpl108(Gen1)) + 
     -      5.925925925925926D0*(WProd89*Cpl109(Gen1)) - 
     -      4.7407407407407405D0*(WProd89*Cpl110(Gen1)) + 
     -      0.6666666666666666D0*(WProd89*Cpl124(Gen1)) + 
     -      (-(1.1111111111111112D0*(WProd89*Cpl108(Gen1))) + 
     -         5.925925925925926D0*(WProd89*Cpl109(Gen1)) - 
     -         4.7407407407407405D0*(WProd89*Cpl110(Gen1)) - 
     -         0.3333333333333333D0*(WProd89*Cpl124(Gen1)))*
     -       LInt254(Gen1) + 
     -      (1.3333333333333333D0*(Cpl101*WProd82*WProd89) - 
     -         3.5555555555555554D0*(Cpl107*WProd82*WProd89) - 
     -         0.2222222222222222D0*(Cpl112*WProd82*WProd89) + 
     -         3.5555555555555554D0*(Cpl96*WProd82*WProd89))*
     -       LInt255(Gen1) + 
     -      (2.2222222222222223D0*(Cpl101*WProd82*WProd89) - 
     -         2.3703703703703702D0*(Cpl107*WProd82*WProd89) - 
     -         0.6666666666666666D0*(Cpl112*WProd82*WProd89) + 
     -         1.1851851851851851D0*(Cpl96*WProd82*WProd89))*
     -       LInt256(Gen1) + 
     -      (3.111111111111111D0*(Cpl101*WProd82*WProd89) - 
     -         5.925925925925926D0*(Cpl107*WProd82*WProd89) - 
     -         0.6666666666666666D0*(Cpl112*WProd82*WProd89) + 
     -         4.7407407407407405D0*(Cpl96*WProd82*WProd89))*
     -       LInt257(Gen1) + 
     -      (3.111111111111111D0*(WProd89*Cpl108(Gen1)) - 
     -         5.925925925925926D0*(WProd89*Cpl109(Gen1)) + 
     -         4.7407407407407405D0*(WProd89*Cpl110(Gen1)) - 
     -         0.6666666666666666D0*(WProd89*Cpl124(Gen1)))*
     -       LInt355(Gen1))
	enddo

#ifdef DEBUG
	print *, 'amp(1,1,2) =', amp(1,1,2)
#endif

        tmp17 = -(Cpl95*WProd82*WProd92) + 
     -   0.6666666666666666D0*(Cpl94*WProd82*WProd92)
        tmp18 = -(1.1851851851851851D0*(Cpl93*WProd82*WProd92)) + 
     -   1.7777777777777777D0*(Cpl94*WProd82*WProd92) - 
     -   0.6666666666666666D0*(Cpl95*WProd82*WProd92)
        tmp19 = -2*Cpl101*WProd82*WProd92 + 
     -   2.6666666666666665D0*(Cpl107*WProd82*WProd92) + 
     -   0.5D0*(Cpl112*WProd82*WProd92) - 
     -   1.1851851851851851D0*(Cpl96*WProd82*WProd92)
        tmp20 = -2*WProd83*WProd92 + 2*WProd85*WProd92 + 
     -   WProd93 - WProd94
        tmp21 = -(WProd83*WProd92) + 2*WProd85*WProd92 + 
     -   WProd93 - WProd94
        tmp22 = WProd82*WProd92 - WProd83*WProd92 + 
     -   2*WProd85*WProd92 + WProd93 - WProd94
        tmp23 = -(WProd82*WProd92) + WProd83*WProd92 + 
     -   3*WProd85*WProd92 + WProd93 - WProd94
        tmp24 = WProd82*WProd92 - WProd83*WProd92 - 
     -   2*WProd85*WProd92 + WProd94
        tmp25 = -2*Cpl101*(WProd85*WProd92 + WProd93) + 
     -   5.333333333333333D0*(Cpl107*(WProd85*WProd92 + WProd93)) - 
     -   3.5555555555555554D0*(Cpl96*(WProd85*WProd92 + WProd93))
        tmp26 = Cpl101*(-2*WProd85*WProd92 - WProd93) + 
     -   2.6666666666666665D0*
     -    (Cpl107*(2*WProd85*WProd92 + WProd93)) - 
     -   1.7777777777777777D0*(Cpl96*(2*WProd85*WProd92 + WProd93))
        amp(1,1,3) = LInt220*tmp25 + LInt221*tmp25 + 
     -   LInt229*tmp26 + LInt231*tmp26 + LInt233*tmp26 + 
     -   LInt222*(-2*Cpl101*tmp23 + 
     -      5.333333333333333D0*(Cpl107*tmp23) - 
     -      3.5555555555555554D0*(Cpl96*tmp23)) + 
     -   LInt223*(-2*Cpl101*tmp23 + 
     -      5.333333333333333D0*(Cpl107*tmp23) - 
     -      3.5555555555555554D0*(Cpl96*tmp23)) + 
     -   LInt224*(-2*Cpl101*tmp23 + 
     -      5.333333333333333D0*(Cpl107*tmp23) - 
     -      3.5555555555555554D0*(Cpl96*tmp23)) + 
     -   LInt228*(2*Cpl101*tmp24 - 
     -      5.333333333333333D0*(Cpl107*tmp24) + 
     -      3.5555555555555554D0*(Cpl96*tmp24)) + 
     -   LInt232*(2*Cpl101*tmp24 - 
     -      5.333333333333333D0*(Cpl107*tmp24) + 
     -      3.5555555555555554D0*(Cpl96*tmp24)) + 
     -   LInt215*(-2*Cpl101*WProd92 + 
     -      5.333333333333333D0*(Cpl107*WProd92) - 
     -      3.5555555555555554D0*(Cpl96*WProd92)) + 
     -   LInt217*(2*Cpl101*WProd92 - 
     -      5.333333333333333D0*(Cpl107*WProd92) + 
     -      3.5555555555555554D0*(Cpl96*WProd92)) + 
     -   LInt226*(Cpl101*(WProd82 - WProd85)*WProd92 - 
     -      2.6666666666666665D0*
     -       (Cpl107*(WProd82 - WProd85)*WProd92) + 
     -      1.7777777777777777D0*(Cpl96*(WProd82 - WProd85)*WProd92))
     -     + LInt227*(-(Cpl101*WProd85*WProd92) + 
     -      2.6666666666666665D0*(Cpl107*WProd85*WProd92) - 
     -      1.7777777777777777D0*(Cpl96*WProd85*WProd92)) + 
     -   Prp18**2*(LInt251*
     -       (-(Cpl112*WProd92) + 
     -         1.3333333333333333D0*(Cpl101*WProd92)) + 
     -      LInt250*(Cpl121*WProd92 - 
     -         1.3333333333333333D0*(Cpl102*WProd92)) + 
     -      LInt353*(-(0.3333333333333333D0*(Cpl101*WProd92)) + 
     -         0.25D0*(Cpl112*WProd92)) + 
     -      LInt354*(-(0.3333333333333333D0*(Cpl101*WProd92)) + 
     -         0.25D0*(Cpl112*WProd92)) + 
     -      LInt352*(-6*Cpl92*WProd92 - Cpl94*WProd92 + 
     -         0.5D0*(Cpl107*WProd92) + 4.5D0*(Cpl111*WProd92) + 
     -         1.3333333333333333D0*(Cpl93*WProd92) - 
     -         0.6666666666666666D0*(Cpl96*WProd92)) + 
     -      LInt236*(-(Cpl107*WProd92) - 9*Cpl111*WProd92 + 
     -         12*Cpl92*WProd92 + 2*Cpl94*WProd92 - 
     -         2.6666666666666665D0*(Cpl93*WProd92) + 
     -         1.3333333333333333D0*(Cpl96*WProd92)) + 
     -      0.2222222222222222D0*(Cpl101*WProd82*WProd92) + 
     -      LInt252*(-(Cpl112*WProd82*WProd92) + 
     -         1.3333333333333333D0*(Cpl101*WProd82*WProd92)) - 
     -      0.6666666666666666D0*(Cpl111*WProd82*WProd92) - 
     -      0.16666666666666666D0*(Cpl112*WProd82*WProd92) + 
     -      0.8888888888888888D0*(Cpl92*WProd82*WProd92) + 
     -      LInt237*(-2*Cpl111*WProd82*WProd92 + 
     -         2.6666666666666665D0*(Cpl92*WProd82*WProd92)) + 
     -      LInt235*(-2*Cpl120*WProd92 + 2*Cpl122*WProd92 - 
     -         5*Cpl111*(WProd83 + WProd85)*WProd92 - 
     -         2.6666666666666665D0*(Cpl103*WProd92) + 
     -         2.6666666666666665D0*(Cpl99*WProd92) + 
     -         6.666666666666667D0*
     -          (Cpl92*(WProd83 + WProd85)*WProd92))) + 
     -   Prp17**2*(16*Cpl92*LInt236*WProd92 - 
     -      8*Cpl92*LInt352*WProd92 + 
     -      0.8888888888888888D0*(Cpl92*WProd82*WProd92) + 
     -      2.6666666666666665D0*(Cpl92*LInt237*WProd82*WProd92) + 
     -      6.666666666666667D0*
     -       (Cpl92*LInt235*(WProd83 + WProd85)*WProd92))
        amp(1,1,3) = amp(1,1,3) + 
     -   Prp18*(LInt271*tmp17 + LInt273*tmp17 + LInt245*tmp19 + 
     -      LInt247*tmp19 - Cpl101*WProd92 + Cpl111*WProd92 - 
     -      2*Cpl111*LInt235*WProd92 - 4*Cpl111*LInt268*WProd92 + 
     -      2*Cpl111*LInt270*WProd82*WProd92 + 
     -      2*Cpl111*LInt272*WProd82*WProd92 + 
     -      2.3333333333333335D0*(Cpl107*WProd92) + 
     -      0.25D0*(Cpl112*WProd92) + 
     -      0.3333333333333333D0*(Cpl94*WProd92) + 
     -      LInt266*(Cpl95*WProd92 - 
     -         0.6666666666666666D0*(Cpl94*WProd92)) + 
     -      LInt269*(-2*Cpl95*WProd92 + 
     -         1.3333333333333333D0*(Cpl94*WProd92)) - 
     -      0.5D0*(Cpl95*WProd92) + 
     -      LInt239*(2*Cpl101*WProd92 - 
     -         2.6666666666666665D0*(Cpl107*WProd92) - 
     -         0.5D0*(Cpl112*WProd92) + 
     -         1.1851851851851851D0*(Cpl96*WProd92)) - 
     -      1.9259259259259258D0*(Cpl96*WProd92) + 
     -      LInt243*(-4*Cpl101*WProd92 + Cpl112*WProd92 + 
     -         5.333333333333333D0*(Cpl107*WProd92) - 
     -         2.3703703703703702D0*(Cpl96*WProd92)) + 
     -      LInt238*(-2*Cpl107*WProd92 + 
     -         2.6666666666666665D0*(Cpl96*WProd92)) + 
     -      LInt242*(4*Cpl107*WProd92 - 
     -         5.333333333333333D0*(Cpl96*WProd92)) + 
     -      LInt267*(Cpl113*WProd92 - Cpl95*WProd82*WProd92 - 
     -         0.6666666666666666D0*(Cpl119*WProd92) + 
     -         0.6666666666666666D0*(Cpl94*WProd82*WProd92)) + 
     -      LInt241*(2*Cpl115*WProd92 - 2*Cpl101*WProd82*WProd92 - 
     -         2.6666666666666665D0*(Cpl114*WProd92) - 
     -         0.5D0*(Cpl116*WProd92) + 
     -         1.1851851851851851D0*(Cpl97*WProd92) + 
     -         2.6666666666666665D0*(Cpl107*WProd82*WProd92) + 
     -         0.5D0*(Cpl112*WProd82*WProd92) - 
     -         1.1851851851851851D0*(Cpl96*WProd82*WProd92)) + 
     -      LInt246*(-2*Cpl107*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -          WProd92 + 2.6666666666666665D0*
     -          (Cpl96*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd92)
     -         ) + LInt244*
     -       (2*Cpl107*(WProd83 + WProd85)*WProd92 - 
     -         2.6666666666666665D0*
     -          (Cpl96*(WProd83 + WProd85)*WProd92)) + 
     -      LInt240*(-2*Cpl114*WProd92 + 
     -         2*Cpl107*(WProd83 + WProd85)*WProd92 + 
     -         2.6666666666666665D0*(Cpl97*WProd92) - 
     -         2.6666666666666665D0*
     -          (Cpl96*(WProd83 + WProd85)*WProd92))) + 
     -   Prp17*(LInt245*tmp18 + LInt247*tmp18 + 
     -      2*Cpl111*LInt235*WProd92 + 4*Cpl111*LInt268*WProd92 - 
     -      2*Cpl111*LInt270*WProd82*WProd92 - 
     -      2*Cpl111*LInt272*WProd82*WProd92 - 
     -      0.6666666666666666D0*(Cpl111*WProd92) - 
     -      1.9259259259259258D0*(Cpl93*WProd92) + 
     -      0.8888888888888888D0*(Cpl94*WProd92) - 
     -      0.3333333333333333D0*(Cpl95*WProd92) + 
     -      LInt239*(1.1851851851851851D0*(Cpl93*WProd92) - 
     -         1.7777777777777777D0*(Cpl94*WProd92) + 
     -         0.6666666666666666D0*(Cpl95*WProd92)) + 
     -      LInt243*(-(2.3703703703703702D0*(Cpl93*WProd92)) + 
     -         3.5555555555555554D0*(Cpl94*WProd92) - 
     -         1.3333333333333333D0*(Cpl95*WProd92)) + 
     -      2.6666666666666665D0*(Cpl93*LInt238*WProd92) - 
     -      5.333333333333333D0*(Cpl93*LInt242*WProd92) - 
     -      0.6666666666666666D0*(Cpl111*LInt266*WProd92) + 
     -      1.3333333333333333D0*(Cpl111*LInt269*WProd92) + 
     -      LInt267*(-(0.6666666666666666D0*(Cpl120*WProd92)) + 
     -         0.6666666666666666D0*(Cpl111*WProd82*WProd92)) + 
     -      LInt241*(-(1.7777777777777777D0*(Cpl117*WProd92)) + 
     -         0.6666666666666666D0*(Cpl118*WProd92) + 
     -         1.1851851851851851D0*(Cpl98*WProd92) - 
     -         1.1851851851851851D0*(Cpl93*WProd82*WProd92) + 
     -         1.7777777777777777D0*(Cpl94*WProd82*WProd92) - 
     -         0.6666666666666666D0*(Cpl95*WProd82*WProd92)) + 
     -      0.6666666666666666D0*(Cpl111*LInt271*WProd82*WProd92) + 
     -      0.6666666666666666D0*(Cpl111*LInt273*WProd82*WProd92) + 
     -      2.6666666666666665D0*
     -       (Cpl93*LInt246*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -         WProd92) + Prp18*
     -       (LInt352*(-5*Cpl111*WProd92 + Cpl94*WProd92 + 
     -            13.333333333333334D0*(Cpl92*WProd92) - 
     -            2.6666666666666665D0*(Cpl93*WProd92)) + 
     -         LInt236*(10*Cpl111*WProd92 - 2*Cpl94*WProd92 - 
     -            26.666666666666668D0*(Cpl92*WProd92) + 
     -            5.333333333333333D0*(Cpl93*WProd92)) + 
     -         0.6666666666666666D0*(Cpl111*WProd82*WProd92) - 
     -         1.7777777777777777D0*(Cpl92*WProd82*WProd92) + 
     -         LInt237*(2*Cpl111*WProd82*WProd92 - 
     -            5.333333333333333D0*(Cpl92*WProd82*WProd92)) + 
     -         LInt235*(2*Cpl119*WProd92 + 2*Cpl120*WProd92 + 
     -            5*Cpl111*(WProd83 + WProd85)*WProd92 - 
     -            5.333333333333333D0*(Cpl100*WProd92) - 
     -            5.333333333333333D0*(Cpl99*WProd92) - 
     -            13.333333333333334D0*
     -             (Cpl92*(WProd83 + WProd85)*WProd92))) + 
     -      LInt240*(2.6666666666666665D0*(Cpl98*WProd92) - 
     -         2.6666666666666665D0*
     -          (Cpl93*(WProd83 + WProd85)*WProd92)) - 
     -      2.6666666666666665D0*
     -       (Cpl93*LInt244*(WProd83 + WProd85)*WProd92))
        amp(1,1,3) = amp(1,1,3) + 
     -   LInt230*(Cpl101*(2*WProd83*WProd92 - 2*WProd85*WProd92 - 
     -         WProd93 + WProd94) + 
     -      2.6666666666666665D0*(Cpl107*tmp20) - 
     -      1.7777777777777777D0*(Cpl96*tmp20)) + 
     -   LInt218*(Cpl101*(WProd83*WProd92 - 2*WProd85*WProd92 - 
     -         WProd93 + WProd94) + 
     -      2.6666666666666665D0*(Cpl107*tmp21) - 
     -      1.7777777777777777D0*(Cpl96*tmp21)) + 
     -   LInt219*(Cpl101*(-(WProd82*WProd92) + WProd83*WProd92 - 
     -         2*WProd85*WProd92 - WProd93 + WProd94) + 
     -      2.6666666666666665D0*(Cpl107*tmp22) - 
     -      1.7777777777777777D0*(Cpl96*tmp22)) + 
     -   LInt225*(Cpl101*(-4*WProd85*WProd92 - 3*WProd93) + 
     -      2.6666666666666665D0*
     -       (Cpl107*(4*WProd85*WProd92 + 3*WProd93)) - 
     -      1.7777777777777777D0*
     -       (Cpl96*(4*WProd85*WProd92 + 3*WProd93)))

	do Gen1 = 1, 3
        amp(1,1,3) = amp(1,1,3) + 
     -   Prp17**2*(1.5802469135802468D0*(Cpl92*WProd82*WProd92) - 
     -      4.7407407407407405D0*(WProd92*Cpl104(Gen1)) - 
     -      4.7407407407407405D0*
     -       (WProd92*Cpl104(Gen1)*LInt254(Gen1)) + 
     -      3.5555555555555554D0*
     -       (Cpl92*WProd82*WProd92*LInt255(Gen1)) + 
     -      1.1851851851851851D0*
     -       (Cpl92*WProd82*WProd92*LInt256(Gen1)) + 
     -      4.7407407407407405D0*
     -       (Cpl92*WProd82*WProd92*LInt257(Gen1)) + 
     -      4.7407407407407405D0*(WProd92*Cpl104(Gen1)*LInt355(Gen1))
     -      ) + Prp18**2*(-4*WProd92*Cpl108(Gen1) + 
     -      WProd92*Cpl124(Gen1) + 
     -      1.4074074074074074D0*(Cpl101*WProd82*WProd92) - 
     -      2.3703703703703702D0*(Cpl107*WProd82*WProd92) - 
     -      0.3888888888888889D0*(Cpl112*WProd82*WProd92) + 
     -      1.5802469135802468D0*(Cpl96*WProd82*WProd92) + 
     -      7.111111111111111D0*(WProd92*Cpl109(Gen1)) - 
     -      4.7407407407407405D0*(WProd92*Cpl110(Gen1)) + 
     -      (-2*WProd92*Cpl108(Gen1) + 
     -         7.111111111111111D0*(WProd92*Cpl109(Gen1)) - 
     -         4.7407407407407405D0*(WProd92*Cpl110(Gen1)) - 
     -         0.5D0*(WProd92*Cpl124(Gen1)))*LInt254(Gen1) + 
     -      (1.7777777777777777D0*(Cpl101*WProd82*WProd92) - 
     -         4.444444444444445D0*(Cpl107*WProd82*WProd92) - 
     -         0.3333333333333333D0*(Cpl112*WProd82*WProd92) + 
     -         3.5555555555555554D0*(Cpl96*WProd82*WProd92))*
     -       LInt255(Gen1) + 
     -      (-(Cpl112*WProd82*WProd92) + 
     -         2.6666666666666665D0*(Cpl101*WProd82*WProd92) - 
     -         2.6666666666666665D0*(Cpl107*WProd82*WProd92) + 
     -         1.1851851851851851D0*(Cpl96*WProd82*WProd92))*
     -       LInt256(Gen1) + 
     -      (4*Cpl101*WProd82*WProd92 - Cpl112*WProd82*WProd92 - 
     -         7.111111111111111D0*(Cpl107*WProd82*WProd92) + 
     -         4.7407407407407405D0*(Cpl96*WProd82*WProd92))*
     -       LInt257(Gen1) + 
     -      (4*WProd92*Cpl108(Gen1) - WProd92*Cpl124(Gen1) - 
     -         7.111111111111111D0*(WProd92*Cpl109(Gen1)) + 
     -         4.7407407407407405D0*(WProd92*Cpl110(Gen1)))*
     -       LInt355(Gen1)) + 
     -   Prp17*Prp18*(3.1604938271604937D0*(Cpl93*WProd82*WProd92) - 
     -      2.3703703703703702D0*(Cpl94*WProd82*WProd92) + 
     -      0.4444444444444444D0*(Cpl95*WProd82*WProd92) + 
     -      7.111111111111111D0*(WProd92*Cpl105(Gen1)) - 
     -      9.481481481481481D0*(WProd92*Cpl106(Gen1)) - 
     -      1.3333333333333333D0*(WProd92*Cpl123(Gen1)) + 
     -      (7.111111111111111D0*(WProd92*Cpl105(Gen1)) - 
     -         9.481481481481481D0*(WProd92*Cpl106(Gen1)) - 
     -         1.3333333333333333D0*(WProd92*Cpl123(Gen1)))*
     -       LInt254(Gen1) + 
     -      (7.111111111111111D0*(Cpl93*WProd82*WProd92) - 
     -         4.444444444444445D0*(Cpl94*WProd82*WProd92) + 
     -         0.6666666666666666D0*(Cpl95*WProd82*WProd92))*
     -       LInt255(Gen1) + 
     -      (2.3703703703703702D0*(Cpl93*WProd82*WProd92) - 
     -         2.6666666666666665D0*(Cpl94*WProd82*WProd92) + 
     -         0.6666666666666666D0*(Cpl95*WProd82*WProd92))*
     -       LInt256(Gen1) + 
     -      (9.481481481481481D0*(Cpl93*WProd82*WProd92) - 
     -         7.111111111111111D0*(Cpl94*WProd82*WProd92) + 
     -         1.3333333333333333D0*(Cpl95*WProd82*WProd92))*
     -       LInt257(Gen1) + 
     -      (-(7.111111111111111D0*(WProd92*Cpl105(Gen1))) + 
     -         9.481481481481481D0*(WProd92*Cpl106(Gen1)) + 
     -         1.3333333333333333D0*(WProd92*Cpl123(Gen1)))*
     -       LInt355(Gen1))
	enddo

#ifdef DEBUG
	print *, 'amp(1,1,3) =', amp(1,1,3)
#endif

        tmp27 = -(0.3333333333333333D0*(Cpl101*WProd95)) + 
     -   0.4166666666666667D0*(Cpl112*WProd95) - 
     -   0.125D0*(Cpl126*WProd95)
	tmp28 = 2*Cpl111*WProd82*WProd95 - Cpl125*WProd82*WProd95
        tmp29 = 0.5D0*(Cpl127*WProd82*WProd95) + 
     -   0.6666666666666666D0*(Cpl94*WProd82*WProd95) - 
     -   1.3333333333333333D0*(Cpl95*WProd82*WProd95)
        tmp30 = -(1.1851851851851851D0*(Cpl93*WProd82*WProd95)) + 
     -   1.7777777777777777D0*(Cpl94*WProd82*WProd95) - 
     -   0.6666666666666666D0*(Cpl95*WProd82*WProd95)
        tmp31 = -(3.3333333333333335D0*(Cpl101*WProd82*WProd95)) + 
     -   3.259259259259259D0*(Cpl107*WProd82*WProd95) + 
     -   1.5D0*(Cpl112*WProd82*WProd95) - 
     -   0.25D0*(Cpl126*WProd82*WProd95) - 
     -   1.1851851851851851D0*(Cpl96*WProd82*WProd95)
        tmp32 = -(8.222222222222221D0*(Cpl101*WProd85*WProd95)) + 
     -   8.88888888888889D0*(Cpl107*WProd85*WProd95) + 
     -   3.3333333333333335D0*(Cpl112*WProd85*WProd95) - 
     -   0.5D0*(Cpl126*WProd85*WProd95) - 
     -   3.5555555555555554D0*(Cpl96*WProd85*WProd95)
        tmp33 = WProd82*WProd95 - WProd83*WProd95 + 
     -   2*WProd85*WProd95 - 2*WProd96 + WProd97 - WProd98
        tmp34 = 2*WProd83*WProd95 - 2*WProd85*WProd95 + 
     -   2*WProd96 - 2*WProd97 + WProd98
        tmp35 = WProd82*WProd95 - WProd83*WProd95 - 
     -   3*WProd85*WProd95 + WProd96 - WProd97 + WProd98
        tmp36 = WProd83*WProd95 - 2*WProd85*WProd95 + 
     -   2*WProd96 - WProd97 + WProd98
        tmp37 = -(WProd82*WProd95) + WProd83*WProd95 - 
     -   2*WProd85*WProd95 + 2*WProd96 - WProd97 + WProd98
        tmp38 = 2*WProd82*WProd95 - 2*WProd83*WProd95 - 
     -   6*WProd85*WProd95 + 4*WProd96 - 3*WProd97 + 2*WProd98
        tmp39 = 2*WProd82*WProd95 - 2*WProd83*WProd95 - 
     -   4*WProd85*WProd95 + 2*WProd96 - WProd97 + 2*WProd98
        tmp40 = 0.5D0*(Cpl126*
     -      (-(WProd85*WProd95) + WProd96 - WProd97)) - 
     -   8.222222222222221D0*
     -    (Cpl101*(WProd85*WProd95 - WProd96 + WProd97)) + 
     -   8.88888888888889D0*(Cpl107*
     -      (WProd85*WProd95 - WProd96 + WProd97)) + 
     -   3.3333333333333335D0*
     -    (Cpl112*(WProd85*WProd95 - WProd96 + WProd97)) - 
     -   3.5555555555555554D0*
     -    (Cpl96*(WProd85*WProd95 - WProd96 + WProd97))
        tmp41 = LInt270*tmp28 + LInt272*tmp28 + LInt271*tmp29 + 
     -   LInt273*tmp29 + LInt245*tmp31 + LInt247*tmp31 - 
     -   Cpl125*WProd95 - Cpl95*WProd95 + 
     -   LInt268*(-4*Cpl111*WProd95 + 2*Cpl125*WProd95) - 
     -   4.166666666666667D0*(Cpl101*WProd95) + 
     -   4.62962962962963D0*(Cpl107*WProd95) + 
     -   1.6666666666666667D0*(Cpl111*WProd95) + 
     -   LInt260*(2*Cpl125*WProd95 - 
     -      2.6666666666666665D0*(Cpl111*WProd95)) + 
     -   LInt235*(2*Cpl125*WProd95 - 
     -      3.3333333333333335D0*(Cpl111*WProd95)) + 
     -   1.6666666666666667D0*(Cpl112*WProd95) - 
     -   0.25D0*(Cpl126*WProd95) + 0.5D0*(Cpl127*WProd95) + 
     -   0.3333333333333333D0*(Cpl94*WProd95) + 
     -   LInt258*(-(0.5D0*(Cpl127*WProd95)) + 
     -      0.6666666666666666D0*(Cpl95*WProd95)) + 
     -   LInt261*(Cpl127*WProd95 - 
     -      1.3333333333333333D0*(Cpl95*WProd95)) + 
     -   LInt266*(-(0.5D0*(Cpl127*WProd95)) - 
     -      0.6666666666666666D0*(Cpl94*WProd95) + 
     -      1.3333333333333333D0*(Cpl95*WProd95)) + 
     -   LInt269*(Cpl127*WProd95 + 
     -      1.3333333333333333D0*(Cpl94*WProd95) - 
     -      2.6666666666666665D0*(Cpl95*WProd95)) + 
     -   LInt239*(3.3333333333333335D0*(Cpl101*WProd95) - 
     -      3.259259259259259D0*(Cpl107*WProd95) - 
     -      1.5D0*(Cpl112*WProd95) + 0.25D0*(Cpl126*WProd95) + 
     -      1.1851851851851851D0*(Cpl96*WProd95)) - 
     -   1.9259259259259258D0*(Cpl96*WProd95) + 
     -   LInt243*(3*Cpl112*WProd95 - 
     -      6.666666666666667D0*(Cpl101*WProd95) + 
     -      6.518518518518518D0*(Cpl107*WProd95) - 
     -      0.5D0*(Cpl126*WProd95) - 
     -      2.3703703703703702D0*(Cpl96*WProd95)) + 
     -   LInt238*(5*Cpl101*WProd95 - 6*Cpl107*WProd95 - 
     -      1.8333333333333333D0*(Cpl112*WProd95) + 
     -      0.25D0*(Cpl126*WProd95) + 
     -      2.6666666666666665D0*(Cpl96*WProd95)) + 
     -   LInt242*(-10*Cpl101*WProd95 + 12*Cpl107*WProd95 + 
     -      3.6666666666666665D0*(Cpl112*WProd95) - 
     -      0.5D0*(Cpl126*WProd95) - 
     -      5.333333333333333D0*(Cpl96*WProd95)) + 
     -   LInt267*(1.3333333333333333D0*(Cpl113*WProd95) - 
     -      0.6666666666666666D0*(Cpl119*WProd95) - 
     -      0.5D0*(Cpl128*WProd95) + 0.5D0*(Cpl127*WProd82*WProd95) + 
     -      0.6666666666666666D0*(Cpl94*WProd82*WProd95) - 
     -      1.3333333333333333D0*(Cpl95*WProd82*WProd95)) + 
     -   LInt241*(-(3.259259259259259D0*(Cpl114*WProd95)) + 
     -      3.3333333333333335D0*(Cpl115*WProd95) - 
     -      1.5D0*(Cpl116*WProd95) + 0.25D0*(Cpl129*WProd95) + 
     -      1.1851851851851851D0*(Cpl97*WProd95) - 
     -      3.3333333333333335D0*(Cpl101*WProd82*WProd95) + 
     -      3.259259259259259D0*(Cpl107*WProd82*WProd95) + 
     -      1.5D0*(Cpl112*WProd82*WProd95) - 
     -      0.25D0*(Cpl126*WProd82*WProd95) - 
     -      1.1851851851851851D0*(Cpl96*WProd82*WProd95)) + 
     -   LInt264*(Cpl125*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -       WProd95 - 1.3333333333333333D0*
     -       (Cpl111*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd95))+
     -     LInt265*(-(0.5D0*
     -         (Cpl127*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd95)
     -         ) + 0.6666666666666666D0*
     -       (Cpl95*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd95))+
     -     LInt246*(5*Cpl101*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -       WProd95 - 6*Cpl107*
     -       (2*WProd82 - 3*WProd83 - 3*WProd85)*WProd95 - 
     -      1.8333333333333333D0*
     -       (Cpl112*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd95)+
     -        0.25D0*(Cpl126*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -         WProd95) + 2.6666666666666665D0*
     -       (Cpl96*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd95))+
     -     LInt262*(-(Cpl125*(WProd83 + WProd85)*WProd95) + 
     -      1.3333333333333333D0*(Cpl111*(WProd83 + WProd85)*WProd95)
     -      ) + LInt263*(0.5D0*
     -       (Cpl127*(WProd83 + WProd85)*WProd95) - 
     -      0.6666666666666666D0*(Cpl95*(WProd83 + WProd85)*WProd95))
     -     + LInt259*(0.6666666666666666D0*(Cpl113*WProd95) - 
     -      0.5D0*(Cpl128*WProd95) + 
     -      0.5D0*(Cpl127*(WProd83 + WProd85)*WProd95) - 
     -      0.6666666666666666D0*(Cpl95*(WProd83 + WProd85)*WProd95))
        tmp41 = tmp41 + LInt244*
     -    (-5*Cpl101*(WProd83 + WProd85)*WProd95 + 
     -      6*Cpl107*(WProd83 + WProd85)*WProd95 + 
     -      1.8333333333333333D0*
     -       (Cpl112*(WProd83 + WProd85)*WProd95) - 
     -      0.25D0*(Cpl126*(WProd83 + WProd85)*WProd95) - 
     -      2.6666666666666665D0*(Cpl96*(WProd83 + WProd85)*WProd95))
     -     + LInt240*(-6*Cpl114*WProd95 + 5*Cpl115*WProd95 - 
     -      5*Cpl101*(WProd83 + WProd85)*WProd95 + 
     -      6*Cpl107*(WProd83 + WProd85)*WProd95 - 
     -      1.8333333333333333D0*(Cpl116*WProd95) + 
     -      0.25D0*(Cpl129*WProd95) + 
     -      2.6666666666666665D0*(Cpl97*WProd95) + 
     -      1.8333333333333333D0*
     -       (Cpl112*(WProd83 + WProd85)*WProd95) - 
     -      0.25D0*(Cpl126*(WProd83 + WProd85)*WProd95) - 
     -      2.6666666666666665D0*(Cpl96*(WProd83 + WProd85)*WProd95))
        amp(1,1,4) = LInt229*tmp32 + LInt231*tmp32 + 
     -   LInt233*tmp32 + 2*Cpl125*LInt277*tmp35 + 
     -   Cpl125*LInt275*tmp37 + LInt220*tmp40 + LInt221*tmp40 + 
     -   Prp18*tmp41 + 8*Cpl125*LInt274*WProd95 - 
     -   Cpl125*LInt279*WProd85*WProd95 - 
     -   2*Cpl125*LInt280*WProd85*WProd95 - 
     -   2*Cpl125*LInt281*WProd85*WProd95 - 
     -   2*Cpl125*LInt282*WProd85*WProd95 + 
     -   LInt230*(4.111111111111111D0*(Cpl101*tmp34) - 
     -      4.444444444444445D0*(Cpl107*tmp34) - 
     -      1.6666666666666667D0*(Cpl112*tmp34) + 
     -      0.25D0*(Cpl126*tmp34) + 1.7777777777777777D0*(Cpl96*tmp34))
     -     + LInt222*(8.222222222222221D0*(Cpl101*tmp35) - 
     -      8.88888888888889D0*(Cpl107*tmp35) - 
     -      3.3333333333333335D0*(Cpl112*tmp35) + 
     -      0.5D0*(Cpl126*tmp35) + 3.5555555555555554D0*(Cpl96*tmp35))+
     -     LInt223*(8.222222222222221D0*(Cpl101*tmp35) - 
     -      8.88888888888889D0*(Cpl107*tmp35) - 
     -      3.3333333333333335D0*(Cpl112*tmp35) + 
     -      0.5D0*(Cpl126*tmp35) + 3.5555555555555554D0*(Cpl96*tmp35))+
     -     LInt218*(4.111111111111111D0*(Cpl101*tmp36) - 
     -      4.444444444444445D0*(Cpl107*tmp36) - 
     -      1.6666666666666667D0*(Cpl112*tmp36) + 
     -      0.25D0*(Cpl126*tmp36) + 1.7777777777777777D0*(Cpl96*tmp36))
     -     + LInt219*(-(4.111111111111111D0*(Cpl101*tmp33)) + 
     -      4.444444444444445D0*(Cpl107*tmp33) + 
     -      1.6666666666666667D0*(Cpl112*tmp33) - 
     -      1.7777777777777777D0*(Cpl96*tmp33) + 0.25D0*(Cpl126*tmp37))
     -     + LInt224*(4.111111111111111D0*(Cpl101*tmp38) - 
     -      4.444444444444445D0*(Cpl107*tmp38) - 
     -      1.6666666666666667D0*(Cpl112*tmp38) + 
     -      0.25D0*(Cpl126*tmp38) + 1.7777777777777777D0*(Cpl96*tmp38))
     -     + LInt228*(4.111111111111111D0*(Cpl101*tmp39) - 
     -      4.444444444444445D0*(Cpl107*tmp39) - 
     -      1.6666666666666667D0*(Cpl112*tmp39) + 
     -      0.25D0*(Cpl126*tmp39) + 1.7777777777777777D0*(Cpl96*tmp39))
     -     + LInt232*(4.111111111111111D0*(Cpl101*tmp39) - 
     -      4.444444444444445D0*(Cpl107*tmp39) - 
     -      1.6666666666666667D0*(Cpl112*tmp39) + 
     -      0.25D0*(Cpl126*tmp39) + 1.7777777777777777D0*(Cpl96*tmp39))
     -     + LInt215*(-(8.222222222222221D0*(Cpl101*WProd95)) + 
     -      8.88888888888889D0*(Cpl107*WProd95) + 
     -      3.3333333333333335D0*(Cpl112*WProd95) - 
     -      0.5D0*(Cpl126*WProd95) - 
     -      3.5555555555555554D0*(Cpl96*WProd95)) + 
     -   LInt216*(-10*Cpl112*WProd95 + 
     -      24.666666666666668D0*(Cpl101*WProd95) - 
     -      26.666666666666668D0*(Cpl107*WProd95) + 
     -      1.5D0*(Cpl126*WProd95) + 
     -      10.666666666666666D0*(Cpl96*WProd95)) + 
     -   LInt217*(2*Cpl126*WProd95 + 
     -      32.888888888888886D0*(Cpl101*WProd95) - 
     -      35.55555555555556D0*(Cpl107*WProd95) - 
     -      13.333333333333334D0*(Cpl112*WProd95) + 
     -      14.222222222222221D0*(Cpl96*WProd95)) + 
     -   LInt226*(4.111111111111111D0*
     -       (Cpl101*(WProd82 - WProd85)*WProd95) - 
     -      4.444444444444445D0*
     -       (Cpl107*(WProd82 - WProd85)*WProd95) - 
     -      1.6666666666666667D0*
     -       (Cpl112*(WProd82 - WProd85)*WProd95) + 
     -      0.25D0*(Cpl126*(WProd82 - WProd85)*WProd95) + 
     -      1.7777777777777777D0*(Cpl96*(WProd82 - WProd85)*WProd95))
     -     + LInt227*(-(4.111111111111111D0*
     -         (Cpl101*WProd85*WProd95)) + 
     -      4.444444444444445D0*(Cpl107*WProd85*WProd95) + 
     -      1.6666666666666667D0*(Cpl112*WProd85*WProd95) - 
     -      0.25D0*(Cpl126*WProd85*WProd95) - 
     -      1.7777777777777777D0*(Cpl96*WProd85*WProd95)) + 
     -   Prp17**2*(16*Cpl92*LInt236*WProd95 - 
     -      8*Cpl92*LInt352*WProd95 + 
     -      0.8888888888888888D0*(Cpl92*WProd82*WProd95) + 
     -      2.6666666666666665D0*(Cpl92*LInt237*WProd82*WProd95) + 
     -      6.666666666666667D0*
     -       (Cpl92*LInt235*(WProd83 + WProd85)*WProd95))
        amp(1,1,4) = amp(1,1,4) + 
     -   Prp18**2*(LInt353*tmp27 + LInt354*tmp27 + 
     -      LInt251*(1.3333333333333333D0*(Cpl101*WProd95) - 
     -         1.6666666666666667D0*(Cpl112*WProd95) + 
     -         0.5D0*(Cpl126*WProd95)) + 
     -      LInt250*(-(1.3333333333333333D0*(Cpl102*WProd95)) + 
     -         1.6666666666666667D0*(Cpl121*WProd95) - 
     -         0.5D0*(Cpl130*WProd95)) + 
     -      LInt352*(-6*Cpl92*WProd95 - 0.25D0*(Cpl101*WProd95) + 
     -         0.8333333333333334D0*(Cpl107*WProd95) + 
     -         7.5D0*(Cpl111*WProd95) - 2.25D0*(Cpl125*WProd95) + 
     -         1.3333333333333333D0*(Cpl93*WProd95) - 
     -         1.6666666666666667D0*(Cpl94*WProd95) + 
     -         0.5D0*(Cpl95*WProd95) - 
     -         0.6666666666666666D0*(Cpl96*WProd95)) + 
     -      LInt236*(-15*Cpl111*WProd95 + 12*Cpl92*WProd95 - 
     -         Cpl95*WProd95 + 0.5D0*(Cpl101*WProd95) - 
     -         1.6666666666666667D0*(Cpl107*WProd95) + 
     -         4.5D0*(Cpl125*WProd95) - 
     -         2.6666666666666665D0*(Cpl93*WProd95) + 
     -         3.3333333333333335D0*(Cpl94*WProd95) + 
     -         1.3333333333333333D0*(Cpl96*WProd95)) + 
     -      0.2222222222222222D0*(Cpl101*WProd82*WProd95) - 
     -      1.1111111111111112D0*(Cpl111*WProd82*WProd95) - 
     -      0.2777777777777778D0*(Cpl112*WProd82*WProd95) + 
     -      0.3333333333333333D0*(Cpl125*WProd82*WProd95) + 
     -      0.08333333333333333D0*(Cpl126*WProd82*WProd95) + 
     -      LInt252*(1.3333333333333333D0*(Cpl101*WProd82*WProd95) - 
     -         1.6666666666666667D0*(Cpl112*WProd82*WProd95) + 
     -         0.5D0*(Cpl126*WProd82*WProd95)) + 
     -      0.8888888888888888D0*(Cpl92*WProd82*WProd95) + 
     -      LInt237*(Cpl125*WProd82*WProd95 - 
     -         3.3333333333333335D0*(Cpl111*WProd82*WProd95) + 
     -         2.6666666666666665D0*(Cpl92*WProd82*WProd95)) + 
     -      LInt235*(-(Cpl131*WProd95) + Cpl132*WProd95 - 
     -         2.6666666666666665D0*(Cpl103*WProd95) - 
     -         3.3333333333333335D0*(Cpl120*WProd95) + 
     -         3.3333333333333335D0*(Cpl122*WProd95) + 
     -         2.6666666666666665D0*(Cpl99*WProd95) - 
     -         8.333333333333334D0*
     -          (Cpl111*(WProd83 + WProd85)*WProd95) + 
     -         2.5D0*(Cpl125*(WProd83 + WProd85)*WProd95) + 
     -         6.666666666666667D0*
     -          (Cpl92*(WProd83 + WProd85)*WProd95)))
        amp(1,1,4) = -2*Cpl125*LInt276*
     -    (WProd85*WProd95 - WProd96 + WProd97) - 
     -   2*Cpl125*LInt278*
     -    (2*WProd85*WProd95 - WProd96 + WProd97) + amp(1,1,4) + 
     -   Prp17*(LInt245*tmp30 + LInt247*tmp30 + 
     -      4*Cpl111*LInt268*WProd95 - 
     -      2*Cpl111*LInt270*WProd82*WProd95 - 
     -      2*Cpl111*LInt272*WProd82*WProd95 - 
     -      1.3333333333333333D0*(Cpl111*WProd95) - 
     -      1.9259259259259258D0*(Cpl93*WProd95) + 
     -      2.2222222222222223D0*(Cpl94*WProd95) - 
     -      0.6666666666666666D0*(Cpl95*WProd95) + 
     -      LInt239*(1.1851851851851851D0*(Cpl93*WProd95) - 
     -         1.7777777777777777D0*(Cpl94*WProd95) + 
     -         0.6666666666666666D0*(Cpl95*WProd95)) + 
     -      LInt238*(2.6666666666666665D0*(Cpl93*WProd95) - 
     -         2.6666666666666665D0*(Cpl94*WProd95) + 
     -         0.6666666666666666D0*(Cpl95*WProd95)) + 
     -      LInt243*(-(2.3703703703703702D0*(Cpl93*WProd95)) + 
     -         3.5555555555555554D0*(Cpl94*WProd95) - 
     -         1.3333333333333333D0*(Cpl95*WProd95)) + 
     -      LInt242*(-(5.333333333333333D0*(Cpl93*WProd95)) + 
     -         5.333333333333333D0*(Cpl94*WProd95) - 
     -         1.3333333333333333D0*(Cpl95*WProd95)) + 
     -      3.3333333333333335D0*(Cpl111*LInt235*WProd95) + 
     -      2.6666666666666665D0*(Cpl111*LInt260*WProd95) - 
     -      0.6666666666666666D0*(Cpl111*LInt266*WProd95) + 
     -      1.3333333333333333D0*(Cpl111*LInt269*WProd95) + 
     -      LInt267*(-(0.6666666666666666D0*(Cpl120*WProd95)) + 
     -         0.6666666666666666D0*(Cpl111*WProd82*WProd95)) + 
     -      LInt241*(-(1.7777777777777777D0*(Cpl117*WProd95)) + 
     -         0.6666666666666666D0*(Cpl118*WProd95) + 
     -         1.1851851851851851D0*(Cpl98*WProd95) - 
     -         1.1851851851851851D0*(Cpl93*WProd82*WProd95) + 
     -         1.7777777777777777D0*(Cpl94*WProd82*WProd95) - 
     -         0.6666666666666666D0*(Cpl95*WProd82*WProd95)) + 
     -      0.6666666666666666D0*(Cpl111*LInt271*WProd82*WProd95) + 
     -      0.6666666666666666D0*(Cpl111*LInt273*WProd82*WProd95) + 
     -      LInt246*(2.6666666666666665D0*
     -          (Cpl93*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd95)
     -           - 2.6666666666666665D0*
     -          (Cpl94*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd95)
     -           + 0.6666666666666666D0*
     -          (Cpl95*(2*WProd82 - 3*WProd83 - 3*WProd85)*WProd95)
     -         ) + 1.3333333333333333D0*
     -       (Cpl111*LInt264*(2*WProd82 - 3*WProd83 - 3*WProd85)*
     -         WProd95) + Prp18*
     -       (LInt352*(-(8.333333333333334D0*(Cpl111*WProd95)) + 
     -            13.333333333333334D0*(Cpl92*WProd95) - 
     -            2.6666666666666665D0*(Cpl93*WProd95) + 
     -            1.6666666666666667D0*(Cpl94*WProd95)) + 
     -         LInt236*(16.666666666666668D0*(Cpl111*WProd95) - 
     -            26.666666666666668D0*(Cpl92*WProd95) + 
     -            5.333333333333333D0*(Cpl93*WProd95) - 
     -            3.3333333333333335D0*(Cpl94*WProd95)) + 
     -         1.1111111111111112D0*(Cpl111*WProd82*WProd95) - 
     -         1.7777777777777777D0*(Cpl92*WProd82*WProd95) + 
     -         LInt237*(3.3333333333333335D0*
     -             (Cpl111*WProd82*WProd95) - 
     -            5.333333333333333D0*(Cpl92*WProd82*WProd95)) + 
     -         LInt235*(-(5.333333333333333D0*(Cpl100*WProd95)) + 
     -            3.3333333333333335D0*(Cpl119*WProd95) + 
     -            3.3333333333333335D0*(Cpl120*WProd95) - 
     -            5.333333333333333D0*(Cpl99*WProd95) + 
     -            8.333333333333334D0*
     -             (Cpl111*(WProd83 + WProd85)*WProd95) - 
     -            13.333333333333334D0*
     -             (Cpl92*(WProd83 + WProd85)*WProd95))) + 
     -      LInt244*(-(2.6666666666666665D0*
     -            (Cpl93*(WProd83 + WProd85)*WProd95)) + 
     -         2.6666666666666665D0*
     -          (Cpl94*(WProd83 + WProd85)*WProd95) - 
     -         0.6666666666666666D0*
     -          (Cpl95*(WProd83 + WProd85)*WProd95)) + 
     -      LInt240*(-(2.6666666666666665D0*(Cpl117*WProd95)) + 
     -         0.6666666666666666D0*(Cpl118*WProd95) + 
     -         2.6666666666666665D0*(Cpl98*WProd95) - 
     -         2.6666666666666665D0*
     -          (Cpl93*(WProd83 + WProd85)*WProd95) + 
     -         2.6666666666666665D0*
     -          (Cpl94*(WProd83 + WProd85)*WProd95) - 
     -         0.6666666666666666D0*
     -          (Cpl95*(WProd83 + WProd85)*WProd95)) - 
     -      1.3333333333333333D0*
     -       (Cpl111*LInt262*(WProd83 + WProd85)*WProd95)) + 
     -   LInt225*(0.5D0*(Cpl126*
     -         (-2*WProd85*WProd95 + WProd96 - WProd97)) - 
     -      8.222222222222221D0*
     -       (Cpl101*(2*WProd85*WProd95 - WProd96 + WProd97)) + 
     -      8.88888888888889D0*
     -       (Cpl107*(2*WProd85*WProd95 - WProd96 + WProd97)) + 
     -      3.3333333333333335D0*
     -       (Cpl112*(2*WProd85*WProd95 - WProd96 + WProd97)) - 
     -      3.5555555555555554D0*
     -       (Cpl96*(2*WProd85*WProd95 - WProd96 + WProd97)))

	do Gen1 = 1, 3
        amp(1,1,4) = amp(1,1,4) + 
     -   Prp17**2*(1.5802469135802468D0*(Cpl92*WProd82*WProd95) - 
     -      4.7407407407407405D0*(WProd95*Cpl104(Gen1)) - 
     -      4.7407407407407405D0*
     -       (WProd95*Cpl104(Gen1)*LInt254(Gen1)) + 
     -      3.5555555555555554D0*
     -       (Cpl92*WProd82*WProd95*LInt255(Gen1)) + 
     -      1.1851851851851851D0*
     -       (Cpl92*WProd82*WProd95*LInt256(Gen1)) + 
     -      4.7407407407407405D0*
     -       (Cpl92*WProd82*WProd95*LInt257(Gen1)) + 
     -      4.7407407407407405D0*(WProd95*Cpl104(Gen1)*LInt355(Gen1))
     -      ) + Prp17*Prp18*
     -    (3.1604938271604937D0*(Cpl93*WProd82*WProd95) - 
     -      3.1604938271604937D0*(Cpl94*WProd82*WProd95) + 
     -      0.7407407407407407D0*(Cpl95*WProd82*WProd95) + 
     -      9.481481481481481D0*(WProd95*Cpl105(Gen1)) - 
     -      9.481481481481481D0*(WProd95*Cpl106(Gen1)) - 
     -      2.2222222222222223D0*(WProd95*Cpl123(Gen1)) + 
     -      (9.481481481481481D0*(WProd95*Cpl105(Gen1)) - 
     -         9.481481481481481D0*(WProd95*Cpl106(Gen1)) - 
     -         2.2222222222222223D0*(WProd95*Cpl123(Gen1)))*
     -       LInt254(Gen1) + 
     -      (7.111111111111111D0*(Cpl93*WProd82*WProd95) - 
     -         6.222222222222222D0*(Cpl94*WProd82*WProd95) + 
     -         1.1111111111111112D0*(Cpl95*WProd82*WProd95))*
     -       LInt255(Gen1) + 
     -      (2.3703703703703702D0*(Cpl93*WProd82*WProd95) - 
     -         3.259259259259259D0*(Cpl94*WProd82*WProd95) + 
     -         1.1111111111111112D0*(Cpl95*WProd82*WProd95))*
     -       LInt256(Gen1) + 
     -      (9.481481481481481D0*(Cpl93*WProd82*WProd95) - 
     -         9.481481481481481D0*(Cpl94*WProd82*WProd95) + 
     -         2.2222222222222223D0*(Cpl95*WProd82*WProd95))*
     -       LInt257(Gen1) + 
     -      (-(9.481481481481481D0*(WProd95*Cpl105(Gen1))) + 
     -         9.481481481481481D0*(WProd95*Cpl106(Gen1)) + 
     -         2.2222222222222223D0*(WProd95*Cpl123(Gen1)))*
     -       LInt355(Gen1)) + 
     -   Prp18**2*(3*WProd95*Cpl124(Gen1) + 
     -      2.5925925925925926D0*(Cpl101*WProd82*WProd95) - 
     -      3.1604938271604937D0*(Cpl107*WProd82*WProd95) - 
     -      1.0925925925925926D0*(Cpl112*WProd82*WProd95) + 
     -      0.19444444444444445D0*(Cpl126*WProd82*WProd95) + 
     -      1.5802469135802468D0*(Cpl96*WProd82*WProd95) - 
     -      7.555555555555555D0*(WProd95*Cpl108(Gen1)) + 
     -      9.481481481481481D0*(WProd95*Cpl109(Gen1)) - 
     -      4.7407407407407405D0*(WProd95*Cpl110(Gen1)) - 
     -      0.5D0*(WProd95*Cpl133(Gen1)) + 
     -      (-(5.555555555555555D0*(WProd95*Cpl108(Gen1))) + 
     -         9.481481481481481D0*(WProd95*Cpl109(Gen1)) - 
     -         4.7407407407407405D0*(WProd95*Cpl110(Gen1)) + 
     -         0.5D0*(WProd95*Cpl124(Gen1)) + 
     -         0.25D0*(WProd95*Cpl133(Gen1)))*LInt254(Gen1) + 
     -      (4*Cpl101*WProd82*WProd95 - 
     -         6.222222222222222D0*(Cpl107*WProd82*WProd95) - 
     -         1.2222222222222223D0*(Cpl112*WProd82*WProd95) + 
     -         0.16666666666666666D0*(Cpl126*WProd82*WProd95) + 
     -         3.5555555555555554D0*(Cpl96*WProd82*WProd95))*
     -       LInt255(Gen1) + 
     -      (4*Cpl101*WProd82*WProd95 - 
     -         3.259259259259259D0*(Cpl107*WProd82*WProd95) - 
     -         2.3333333333333335D0*(Cpl112*WProd82*WProd95) + 
     -         0.5D0*(Cpl126*WProd82*WProd95) + 
     -         1.1851851851851851D0*(Cpl96*WProd82*WProd95))*
     -       LInt256(Gen1) + 
     -      (-3*Cpl112*WProd82*WProd95 + 
     -         7.555555555555555D0*(Cpl101*WProd82*WProd95) - 
     -         9.481481481481481D0*(Cpl107*WProd82*WProd95) + 
     -         0.5D0*(Cpl126*WProd82*WProd95) + 
     -         4.7407407407407405D0*(Cpl96*WProd82*WProd95))*
     -       LInt257(Gen1) + 
     -      (-3*WProd95*Cpl124(Gen1) + 
     -         7.555555555555555D0*(WProd95*Cpl108(Gen1)) - 
     -         9.481481481481481D0*(WProd95*Cpl109(Gen1)) + 
     -         4.7407407407407405D0*(WProd95*Cpl110(Gen1)) + 
     -         0.5D0*(WProd95*Cpl133(Gen1)))*LInt355(Gen1))
	enddo

#ifdef DEBUG
	print *, 'amp(1,1,4) =', amp(1,1,4)
#endif

      endif

      end


