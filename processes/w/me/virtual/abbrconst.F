      subroutine abbr_const()
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "params.h" 

      integer resnr,delta
      double complex myA0,myB0,myB1,myB00,myB11,myB001,myB111,myDB0
      double complex myDB1,myDB00,myDB11,myBval,myCval,myDval,myEval
      double complex Prp,myDeltaIR,myDeltaUV
      integer dr

      integer myBget,myCget,myDget,myEget

#include "looptools.h" 
#include "amps/ampvars.h" 
#include "amps/ampcommon.h"
#include "amps/ampcom.h" 
#include "regcommon.h"
	Cpl1 = (Alfa*Pi)/SW2c
	Cpl10 = Alfa/Pi*CW2c/SW2c
	Cpl101 = Alfa*Pi*MW2c/SW2c
	Cpl102 = Alfa*Pi*1/(SW2c*SWc)
	Cpl103 = SWc/CWc
	Cpl104 = Alfa/Pi/CW2c
	Cpl105 = 1/Pi*Log(lambda**2)
	Cpl106 = 1/Pi*Log(mureg)
	Cpl107 = 1/Pi*myDeltaIR()
	Cpl108 = 1/Pi*myDeltaUV()
	Cpl109 = Alfa*Pi*EL/SW2c
	Cpl11 = Alfa/Pi*MW2r
	Cpl12 = Alfa/Pi*(CW2c*MW2c)/SW2c
	Cpl13 = Alfa/Pi*(CW2c*MW2r)/SW2c
	Cpl14 = Alfa/Pi*(CW2c*MWr**2)/SW2c
	Cpl15 = Alfa/Pi*(CW2c*MZ2c)/SW2c
	Cpl16 = Alfa/Pi*(MW2c*SW2c)/CW2c
	Cpl17 = Alfa/Pi*MW2c**2
	Cpl18 = Alfa/Pi*(MW2c*MW2r)
	Cpl19 = Alfa/Pi*(MW2c*MWr**2)
	Cpl2 = Alfa/SW2c
	Cpl20 = Alfa/Pi*(MW2r*MWr**2)
	Cpl21 = Alfa/Pi*MW2c**2/SW2c
	Cpl22 = Alfa/Pi*(MW2c*MW2r)/SW2c
	Cpl23 = Alfa/Pi*(CW2c*MW2c*MWr**2)/SW2c
	Cpl24 = Alfa/Pi*(CW2c*MW2r*MWr**2)/SW2c
	Cpl25 = Alfa/Pi*(CW2c*MW2c*MZ2c)/SW2c
	Cpl26 = Alfa/Pi*(CW2c*MW2r*MZ2c)/SW2c
	Cpl27 = Alfa/Pi*(MW2c**2*SW2c)/CW2c
	Cpl28 = Alfa/Pi*(MW2c*MW2r*SW2c)/CW2c
	Cpl3 = Alfa/Pi*MW2c
	Cpl32 = Alfa/Pi*(MW2c*MWr**2)/SW2c
	Cpl33 = Alfa/Pi*(MW2r*MWr**2)/SW2c
	Cpl34 = Alfa/Pi*(CW2c*MZ2r)/SW2c
	Cpl35 = Alfa/Pi*MZr**2/(CW2c*SW2c)
	Cpl36 = Alfa/Pi*(CW2c*MZr**2)/SW2c
	Cpl37 = Alfa/Pi*1/(CW2c*SW2c)
	Cpl38 = Alfa/Pi*SW2c/CW2c
	Cpl39 = Alfa/Pi*MW2c/(CW2c**2*SW2c)
	Cpl4 = Alfa/Pi*MWr**2
	Cpl40 = Alfa/Pi*MZ2c/(CW2c*SW2c)
	Cpl41 = Alfa/Pi*MZ2r/(CW2c*SW2c)
	Cpl42 = Alfa/Pi*(MW2c*MZ2c)/(CW2c**2*SW2c)
	Cpl43 = Alfa/Pi*(MW2c*MZ2r)/(CW2c**2*SW2c)
	Cpl44 = Alfa/Pi*(CW2c*MW2c*MZ2r)/SW2c
	Cpl45 = Alfa/Pi*(CW2c*MZ2c*MZr**2)/SW2c
	Cpl46 = Alfa/Pi*(CW2c*MZ2r*MZr**2)/SW2c
	Cpl47 = Alfa/Pi*(MW2c*MZ2c*SW2c)/CW2c
	Cpl48 = Alfa/Pi*(MW2c*MZ2r*SW2c)/CW2c
	Cpl49 = Alfa/Pi*MZ2c
	Cpl5 = Alfa/Pi*MW2c/SW2c
	Cpl50 = Alfa/Pi*MZ2r
	Cpl51 = Alfa/Pi*(MZ2c*SW2c)/CW2c
	Cpl52 = Alfa/Pi*(MZ2r*SW2c)/CW2c
	Cpl53 = Alfa/Pi*(MZ2c*MZr**2)/(CW2c*SW2c)
	Cpl54 = Alfa/Pi*(MZ2r*MZr**2)/(CW2c*SW2c)
	Cpl55 = Alfa/Pi*MZr**2/CW2c
	Cpl56 = Alfa/Pi*(MZr**2*SW2c)/CW2c
	Cpl6 = Alfa/Pi*MW2r/SW2c
	Cpl60 = Alfa/Pi*MZ2c/CW2c
	Cpl61 = Alfa/Pi*MZ2r/CW2c
	Cpl64 = Alfa/Pi*(MZ2c*MZr**2)/CW2c
	Cpl65 = Alfa/Pi*(MZ2r*MZr**2)/CW2c
	Cpl66 = Alfa/Pi*(MZ2c*MZr**2*SW2c)/CW2c
	Cpl67 = Alfa/Pi*(MZ2r*MZr**2*SW2c)/CW2c
	Cpl68 = CW2c/(MW2c*SWc)
	Cpl69 = CW2c/(MZ2c*SWc)
	Cpl7 = Alfa/Pi*MWr**2/SW2c
	Cpl70 = Alfa*1/(pi*SW2c)
	Cpl71 = Alfa*Log(CW2c)/(pi*SW2c**2)
	Cpl72 = Alfa*Log(CW2c)/(pi*SW2c)
	Cpl73 = Alfa/Pi*1/(MW2c*SW2c)
	Cpl74 = Alfa/Pi/MW2c
	Cpl75 = Alfa/Pi*CW2c/(MW2c*SW2c)
	Cpl76 = Alfa/Pi*(CW2c*MZ2c)/(MW2c*SW2c)
	Cpl77 = 1/MW2c
	Cpl78 = 1/SWc
	Cpl79 = CWc/SWc
	Cpl8 = Alfa/Pi/SW2c
	Cpl81 = Alfa/Pi*CW2c/(CWc*MZ2c*SWc)
	Cpl82 = Alfa/Pi*SW2c/(CWc*MZ2c*SWc)
	Cpl83 = Alfa/Pi*MW2c/(CWc*MZ2c*SWc)
	Cpl86 = Alfa2/SW2c
	Cpl87 = Alfa2*1/(CW2c*SW2c**2)
	Cpl88 = Alfa2*1/(CW2c*SW2c)
	Cpl89 = Alfa2/CW2c
	Cpl9 = Alfa/Pi
	Cpl90 = Alfa2/SW2c**2
	Cpl91 = Alfa2*MZ2c/CW2c
	Cpl92 = Alfa2*MZ2c/(CW2c*SW2c**2)
	Cpl93 = Alfa2*MZ2c/(CW2c*SW2c)
	Cpl94 = Alfa2*CW2c/SW2c**2
	Cpl95 = Alfa2*MW2c/SW2c**2
	Cpl96 = Alfa2*(CW2c*MW2c)/SW2c**2
	Cpl97 = Alfa2*MW2c/CW2c
	Cpl98 = Alfa2*(CW2c*MZ2c)/SW2c**2
	Cpl99 = Alfa2*MW2c/SW2c
      

	do Gen1=1,3
	Cpl100(Gen1) = Alfa2*MT2*delta(Gen1,3)/SW2c**2
	Cpl29(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/SW2c
	Cpl30(Gen1) = (Alfa*MT2)/Pi*(MW2c*delta(Gen1,3))/SW2c
	Cpl31(Gen1) = (Alfa*MT2)/Pi*(MW2r*delta(Gen1,3))/SW2c
	Cpl57(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/CW2c
	Cpl58(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/(CW2c*SW2c)
	Cpl59(Gen1) = (Alfa*MT2)/Pi*(SW2c*delta(Gen1,3))/CW2c
	Cpl62(Gen1) = (Alfa*MT2)/Pi*(MZ2c*delta(Gen1,3))/(CW2c*SW2c)
	Cpl63(Gen1) = (Alfa*MT2)/Pi*(MZ2r*delta(Gen1,3))/(CW2c*SW2c)
	Cpl80(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/(MW2c*SW2c)
	Cpl84(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/(CWc*MZ2c*SWc)
	Cpl85(Gen1) = (Alfa*MT2)/Pi*(SW2c*delta(Gen1,3))/(CWc*MZ2c*SWc)
      
        enddo
      
      
      
	IInt11 = myBget(0.D0,MW2r,MZ2r,0)
	LInt68 = myBval(bb0,IInt11)
	IInt12 = myBget(0.D0,MZ2r,MW2r,2)
	LInt70 = myBval(bb00,IInt12)
      

	do Gen1=1,3
	IInt16(Gen1) = myBget(0.D0,0.D0,MQU2(Gen1),0)
	LInt46(Gen1) = myBval(bb0,IInt16(Gen1))
	IInt23(Gen1) = myBget(0.D0,MQD2(Gen1),MQU2(Gen1),0)
	LInt72(Gen1) = myBval(bb0,IInt23(Gen1))
	IInt24(Gen1) = myBget(0.D0,MQU2(Gen1),MQD2(Gen1),1)
	LInt73(Gen1) = myBval(bb1,IInt24(Gen1))
        enddo

	IInt10 = myBget(0.D0,MH2,MW2r,2)
	LInt67 = myBval(bb0,IInt10)
	LInt69 = myBval(bb00,IInt10)
	IInt3 = myBget(MWr**2,MW2r,MZ2r,0)
	LInt12 = myBval(bb0,IInt3)
	LInt19 = myBval(dbb0,IInt3)
      

	do Gen1=1,3
	IInt13(Gen1) = myBget(MWr**2,MQD2(Gen1),MQU2(Gen1),0)
	LInt25(Gen1) = myBval(bb0,IInt13(Gen1))
	LInt28(Gen1) = myBval(dbb0,IInt13(Gen1))
	IInt20(Gen1) = myBget(0.D0,MLE2(Gen1),MLE2(Gen1),2)
	LInt61(Gen1) = myBval(bb1,IInt20(Gen1))
	LInt64(Gen1) = myBval(dbb00,IInt20(Gen1))
	IInt21(Gen1) = myBget(0.D0,MQD2(Gen1),MQD2(Gen1),2)
	LInt62(Gen1) = myBval(bb1,IInt21(Gen1))
	LInt65(Gen1) = myBval(dbb00,IInt21(Gen1))
        enddo

	IInt5 = myBget(0.D0,0.D0,MW2r,1)
	LInt179 = myBval(bb0,IInt5)
	LInt71 = myBval(bb1,IInt5)
	LInt9 = myBval(bb0,IInt5)
	IInt8 = myBget(MZr**2,0.D0,0.D0,2)
	LInt37 = myBval(bb1,IInt8)
	LInt41 = myBval(dbb00,IInt8)
	LInt44 = myBval(dbb1,IInt8)
      

	do Gen1=1,3
	IInt14(Gen1) = myBget(MWr**2,0.D0,MLE2(Gen1),2)
	LInt26(Gen1) = myBval(bb1,IInt14(Gen1))
	LInt29(Gen1) = myBval(dbb00,IInt14(Gen1))
	LInt31(Gen1) = myBval(dbb1,IInt14(Gen1))
	IInt15(Gen1) = myBget(MWr**2,MQU2(Gen1),MQD2(Gen1),2)
	LInt27(Gen1) = myBval(bb1,IInt15(Gen1))
	LInt30(Gen1) = myBval(dbb00,IInt15(Gen1))
	LInt32(Gen1) = myBval(dbb1,IInt15(Gen1))
	IInt18(Gen1) = myBget(MZr**2,MLE2(Gen1),MLE2(Gen1),2)
	LInt48(Gen1) = myBval(bb1,IInt18(Gen1))
	LInt52(Gen1) = myBval(dbb00,IInt18(Gen1))
	LInt55(Gen1) = myBval(dbb1,IInt18(Gen1))
	IInt19(Gen1) = myBget(MZr**2,MQD2(Gen1),MQD2(Gen1),2)
	LInt49(Gen1) = myBval(bb1,IInt19(Gen1))
	LInt53(Gen1) = myBval(dbb00,IInt19(Gen1))
	LInt56(Gen1) = myBval(dbb1,IInt19(Gen1))
	IInt22(Gen1) = myBget(0.D0,MQU2(Gen1),MQU2(Gen1),2)
	LInt63(Gen1) = myBval(bb1,IInt22(Gen1))
	LInt66(Gen1) = myBval(dbb00,IInt22(Gen1))
	LInt75(Gen1) = myBval(bb0,IInt22(Gen1))
        enddo

	IInt2 = myBget(MWr**2,MH2,MW2r,2)
	LInt11 = myBval(bb0,IInt2)
	LInt13 = myBval(bb00,IInt2)
	LInt18 = myBval(dbb0,IInt2)
	LInt21 = myBval(dbb00,IInt2)
	IInt4 = myBget(MWr**2,MZ2r,MW2r,2)
	LInt14 = myBval(bb00,IInt4)
	LInt16 = myBval(bb1,IInt4)
	LInt22 = myBval(dbb00,IInt4)
	LInt24 = myBval(dbb1,IInt4)
	IInt6 = myBget(MZr**2,MH2,MZ2r,2)
	LInt33 = myBval(bb0,IInt6)
	LInt35 = myBval(bb00,IInt6)
	LInt39 = myBval(dbb0,IInt6)
	LInt42 = myBval(dbb00,IInt6)
	IInt9 = myBget(0.D0,MW2r,MW2r,2)
	LInt58 = myBval(bb0,IInt9)
	LInt59 = myBval(bb1,IInt9)
	LInt60 = myBval(dbb00,IInt9)
	LInt74 = myBval(bb00,IInt9)
	IInt1 = myBget(MWr**2,0.D0,MW2r,2)
	LInt10 = myBval(bb0,IInt1)
	LInt15 = myBval(bb1,IInt1)
	LInt17 = myBval(dbb0,IInt1)
	LInt20 = myBval(dbb00,IInt1)
	LInt23 = myBval(dbb1,IInt1)
      

	do Gen1=1,3
	IInt17(Gen1) = myBget(MZr**2,MQU2(Gen1),MQU2(Gen1),2)
	LInt47(Gen1) = myBval(bb0,IInt17(Gen1))
	LInt50(Gen1) = myBval(bb1,IInt17(Gen1))
	LInt51(Gen1) = myBval(dbb0,IInt17(Gen1))
	LInt54(Gen1) = myBval(dbb00,IInt17(Gen1))
	LInt57(Gen1) = myBval(dbb1,IInt17(Gen1))
        enddo

	IInt7 = myBget(MZr**2,MW2r,MW2r,2)
	LInt34 = myBval(bb0,IInt7)
	LInt36 = myBval(bb00,IInt7)
	LInt38 = myBval(bb1,IInt7)
	LInt40 = myBval(dbb0,IInt7)
	LInt43 = myBval(dbb00,IInt7)
	LInt45 = myBval(dbb1,IInt7)
	LInt176 = myA0(MH2)
	LInt177 = myA0(MW2r)
	LInt178 = myA0(MZ2r)
	LInt180 = myB1(0.D0,0.D0,MZ2r)

	LInt181 = myB1(0.D0,ME2,MW2r)
	LInt182 = myB1(ME2,0.D0,MW2r)
	LInt183 = myB1(ME2,ME2,0.D0)
	LInt184 = myB1(ME2,ME2,MZ2r)
	LInt185 = ME2*myDB0(ME2,0.D0,ME2)
	LInt186 = ME2*myDB1(ME2,ME2,0.D0)
	LInt187 = myB1(MU2,MD2,MW2r)
	LInt188 = myB1(MU2,MU2,0.D0)
	LInt189 = myB1(MU2,MU2,MZ2r)
	LInt190 = MU2*myDB0(MU2,0.D0,MU2)
	LInt191 = MU2*myDB1(MU2,MU2,0.D0)
	LInt192 = myB1(MD2,MD2,0.D0)
	LInt193 = myB1(MD2,MD2,MZ2r)
	LInt194 = myB1(MD2,MU2,MW2r)
	LInt195 = MD2*myDB0(MD2,0.D0,MD2)
	LInt196 = MD2*myDB1(MD2,MD2,0.D0)

	LInt6 = myA0(MH2)
	LInt7 = myA0(MW2r)
	LInt8 = myA0(MZ2r)
      
        dMWsq1 = (2*Cpl11 - 2*Cpl3)*LInt20 + 
     -   LInt16*(-(0.5D0*Cpl12) + 0.5D0*Cpl13 - 0.5D0*Cpl14) + 
     -   LInt12*(-(1.25D0*Cpl12) + 1.25D0*Cpl13 - 1.25D0*Cpl14 - 
     -      0.5D0*Cpl15 + 0.25D0*Cpl16) + 
     -   LInt23*(-(0.5D0*Cpl19) + 0.5D0*Cpl20) + 
     -   LInt17*(0.25D0*Cpl17 - 0.25D0*Cpl18 - 1.25D0*Cpl19 + 
     -      1.25D0*Cpl20) + LInt18*(0.25D0*Cpl21 - 0.25D0*Cpl22) + 
     -   LInt24*(-(0.5D0*Cpl23) + 0.5D0*Cpl24) + 
     -   LInt19*(-(1.25D0*Cpl23) + 1.25D0*Cpl24 - 0.5D0*Cpl25 + 
     -      0.5D0*Cpl26 + 0.25D0*Cpl27 - 0.25D0*Cpl28) - 
     -   0.3333333333333333D0*Cpl3 + 0.1111111111111111D0*Cpl4 + 
     -   LInt15*(0.5D0*Cpl11 - 0.16666666666666666D0*Cpl3 - 
     -      0.8333333333333334D0*Cpl4) + 
     -   LInt10*(-Cpl3 + 1.25D0*Cpl11 - 1.25D0*Cpl4) - 
     -   0.16666666666666666D0*Cpl5 + 0.16666666666666666D0*Cpl6 + 
     -   LInt21*(-(0.25D0*Cpl5) + 0.25D0*Cpl6) + 
     -   LInt22*(-2*Cpl12 + 2*Cpl13 - 0.25D0*Cpl5 + 0.25D0*Cpl6) - 
     -   0.16666666666666666D0*Cpl7 + 
     -   LInt8*(0.75D0*Cpl10 + 0.0625D0*Cpl8) + 
     -   LInt14*(-2*Cpl10 - 0.25D0*Cpl8) + 
     -   LInt7*(-(0.5D0*Cpl10) + 0.875D0*Cpl8 - 0.5D0*Cpl9) + 
     -   0.25D0*(Cpl5*LInt11) - 0.25D0*(Cpl8*LInt13) + 
     -   0.0625D0*(Cpl8*LInt6) - 0.3333333333333333D0*(Cpl3*LInt9)

	do Gen1 = 1, 3
        dMWsq1 = dMWsq1 - 0.1111111111111111D0*Cpl7 + 
     -   0.25D0*Cpl29(Gen1) - 0.25D0*(Cpl29(Gen1)*LInt25(Gen1)) + 
     -   (-(0.25D0*Cpl5) + 0.25D0*Cpl6 - 0.16666666666666666D0*Cpl7)*
     -    LInt26(Gen1) + (-(0.75D0*Cpl5) + 0.75D0*Cpl6 - 0.5D0*Cpl7 + 
     -      0.25D0*Cpl29(Gen1))*LInt27(Gen1) + 
     -   (-(0.75D0*Cpl30(Gen1)) + 0.75D0*Cpl31(Gen1))*LInt28(Gen1) + 
     -   (0.5D0*Cpl5 - 0.5D0*Cpl6)*LInt29(Gen1) + 
     -   (1.5D0*Cpl5 - 1.5D0*Cpl6)*LInt30(Gen1) + 
     -   (-(0.25D0*Cpl32) + 0.25D0*Cpl33)*LInt31(Gen1) + 
     -   (-(0.75D0*Cpl32) + 0.75D0*Cpl33)*LInt32(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dMWsq1 =', dMWsq1
#endif

        dMZsq1 = -(0.16666666666666666D0*Cpl15) + 
     -   0.16666666666666666D0*Cpl34 - 0.041666666666666664D0*Cpl35 - 
     -   0.16666666666666666D0*Cpl36 + 
     -   LInt38*(-(0.5D0*Cpl15) + 0.5D0*Cpl34 - 0.5D0*Cpl36) + 
     -   LInt34*(-(0.5D0*Cpl12) - 1.25D0*Cpl15 + 0.5D0*Cpl16 + 
     -      1.25D0*Cpl34 - 1.25D0*Cpl36) + 
     -   LInt42*(-(0.25D0*Cpl40) + 0.25D0*Cpl41) + 
     -   LInt37*(-(0.25D0*Cpl35) - 0.375D0*Cpl40 + 0.375D0*Cpl41) + 
     -   LInt41*(0.75D0*Cpl40 - 0.75D0*Cpl41) + 
     -   LInt39*(0.25D0*Cpl42 - 0.25D0*Cpl43) + 
     -   LInt45*(-(0.5D0*Cpl45) + 0.5D0*Cpl46) + 
     -   LInt40*(-(0.5D0*Cpl25) + 0.5D0*Cpl44 - 1.25D0*Cpl45 + 
     -      1.25D0*Cpl46 + 0.5D0*Cpl47 - 0.5D0*Cpl48) + 
     -   LInt43*(-(2.25D0*Cpl15) + 2.25D0*Cpl34 + 0.5D0*Cpl49 - 
     -      0.5D0*Cpl50 - 0.25D0*Cpl51 + 0.25D0*Cpl52) + 
     -   LInt44*(-(0.375D0*Cpl53) + 0.375D0*Cpl54) + 
     -   LInt7*(1.125D0*Cpl10 + 0.125D0*Cpl38 - 0.25D0*Cpl9) + 
     -   LInt36*(-(2.25D0*Cpl10) - 0.25D0*Cpl38 + 0.5D0*Cpl9) + 
     -   0.25D0*(Cpl39*LInt33) - 0.25D0*(Cpl37*LInt35) + 
     -   0.0625D0*(Cpl37*LInt6) + 0.0625D0*(Cpl37*LInt8)

	do Gen1 = 1, 3
        dMZsq1 = dMZsq1 - 0.09722222222222222D0*Cpl35 + 
     -   0.2222222222222222D0*Cpl55 - 0.2962962962962963D0*Cpl56 - 
     -   0.6666666666666666D0*Cpl57(Gen1) + 0.25D0*Cpl58(Gen1) + 
     -   0.8888888888888888D0*Cpl59(Gen1) + 
     -   (0.6666666666666666D0*Cpl57(Gen1) - 0.25D0*Cpl58(Gen1) - 
     -      0.8888888888888888D0*Cpl59(Gen1))*LInt46(Gen1) + 
     -   (-(0.6666666666666666D0*Cpl57(Gen1)) - 0.125D0*Cpl58(Gen1) + 
     -      0.8888888888888888D0*Cpl59(Gen1))*LInt47(Gen1) + 
     -   (-Cpl51 + Cpl52 - 0.08333333333333333D0*Cpl35 - 
     -      0.125D0*Cpl40 + 0.125D0*Cpl41 + 0.3333333333333333D0*Cpl55 - 
     -      0.6666666666666666D0*Cpl56 + 0.5D0*Cpl60 - 0.5D0*Cpl61)*
     -    LInt48(Gen1) + (-(0.25D0*Cpl35) - 0.375D0*Cpl40 + 
     -      0.375D0*Cpl41 - 0.3333333333333333D0*Cpl51 + 
     -      0.3333333333333333D0*Cpl52 + 0.3333333333333333D0*Cpl55 - 
     -      0.2222222222222222D0*Cpl56 + 0.5D0*Cpl60 - 0.5D0*Cpl61)*
     -    LInt49(Gen1) + (Cpl60 - Cpl61 - 0.25D0*Cpl35 - 
     -      0.375D0*Cpl40 + 0.375D0*Cpl41 - 1.3333333333333333D0*Cpl51 + 
     -      1.3333333333333333D0*Cpl52 + 0.6666666666666666D0*Cpl55 - 
     -      0.8888888888888888D0*Cpl56)*LInt50(Gen1) + 
     -   (-(0.375D0*Cpl62(Gen1)) + 0.375D0*Cpl63(Gen1))*LInt51(Gen1) + 
     -   (2*Cpl51 - 2*Cpl52 - Cpl60 + Cpl61 + 0.25D0*Cpl40 - 
     -      0.25D0*Cpl41)*LInt52(Gen1) + 
     -   (-Cpl60 + Cpl61 + 0.75D0*Cpl40 - 0.75D0*Cpl41 + 
     -      0.6666666666666666D0*Cpl51 - 0.6666666666666666D0*Cpl52)*
     -    LInt53(Gen1) + (-2*Cpl60 + 2*Cpl61 + 0.75D0*Cpl40 - 
     -      0.75D0*Cpl41 + 2.6666666666666665D0*Cpl51 - 
     -      2.6666666666666665D0*Cpl52)*LInt54(Gen1) + 
     -   (-Cpl66 + Cpl67 - 0.125D0*Cpl53 + 0.125D0*Cpl54 + 0.5D0*Cpl64 - 
     -      0.5D0*Cpl65)*LInt55(Gen1) + 
     -   (-(0.375D0*Cpl53) + 0.375D0*Cpl54 + 0.5D0*Cpl64 - 0.5D0*Cpl65 - 
     -      0.3333333333333333D0*Cpl66 + 0.3333333333333333D0*Cpl67)*
     -    LInt56(Gen1) + (Cpl64 - Cpl65 - 0.375D0*Cpl53 + 
     -      0.375D0*Cpl54 - 1.3333333333333333D0*Cpl66 + 
     -      1.3333333333333333D0*Cpl67)*LInt57(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dMZsq1 =', dMZsq1
#endif

        dZAA1 = 3*Cpl9*LInt60 + 0.16666666666666666D0*Cpl9 + 
     -   1.25D0*(Cpl9*LInt58) + 0.5D0*(Cpl9*LInt59)

	do Gen1 = 1, 3
        dZAA1 = dZAA1 + 0.3333333333333333D0*(Cpl9*LInt62(Gen1)) + 
     -   1.3333333333333333D0*(Cpl9*LInt63(Gen1)) - 
     -   0.6666666666666666D0*(Cpl9*LInt65(Gen1)) - 
     -   2.6666666666666665D0*(Cpl9*LInt66(Gen1)) + 
     -   Cpl9*LInt61(Gen1) - 2*Cpl9*LInt64(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dZAA1 =', dZAA1
#endif

        dZW1 = 2*Cpl9*LInt20 + 
     -   LInt19*(1.25D0*Cpl14 + 0.5D0*Cpl15 - 0.25D0*Cpl16) + 
     -   LInt17*(-(0.25D0*Cpl3) + 1.25D0*Cpl4) + 
     -   0.16666666666666666D0*Cpl8 + LInt22*(2*Cpl10 + 0.25D0*Cpl8) + 
     -   1.25D0*(Cpl9*LInt10) + 1.25D0*(Cpl10*LInt12) + 
     -   0.5D0*(Cpl9*LInt15) + 0.5D0*(Cpl10*LInt16) - 
     -   0.25D0*(Cpl5*LInt18) + 0.25D0*(Cpl8*LInt21) + 
     -   0.5D0*(Cpl4*LInt23) + 0.5D0*(Cpl14*LInt24)

	do Gen1 = 1, 3
        dZW1 = dZW1 + 0.25D0*(Cpl8*LInt26(Gen1)) + 
     -   0.75D0*(Cpl8*LInt27(Gen1)) + 
     -   0.75D0*(Cpl29(Gen1)*LInt28(Gen1)) - 
     -   0.5D0*(Cpl8*LInt29(Gen1)) - 1.5D0*(Cpl8*LInt30(Gen1)) + 
     -   0.25D0*(Cpl7*LInt31(Gen1)) + 0.75D0*(Cpl7*LInt32(Gen1))
	enddo

#ifdef DEBUG
	print *, 'dZW1 =', dZW1
#endif

        dZZA1 = Cpl83*LInt58 + (5*Cpl81 - Cpl82)*LInt74 + 
     -   LInt7*(-(2.5D0*Cpl81) + 0.5D0*Cpl82)

	do Gen1 = 1, 3
        dZZA1 = dZZA1 - 0.6666666666666666D0*Cpl84(Gen1) + 
     -   1.7777777777777777D0*Cpl85(Gen1) + 
     -   (0.6666666666666666D0*Cpl84(Gen1) - 
     -      1.7777777777777777D0*Cpl85(Gen1))*LInt46(Gen1) + 
     -   (-(0.6666666666666666D0*Cpl84(Gen1)) + 
     -      1.7777777777777777D0*Cpl85(Gen1))*LInt75(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dZZA1 =', dZZA1
#endif

        dZfL1(1,1,1) = 0.0625D0*Cpl37 + 0.125D0*Cpl8 + 
     -   0.125D0*(Cpl37*LInt180) + 0.25D0*(Cpl8*LInt181)

#ifdef DEBUG
	print *, 'dZfL1(1,1,1) =', dZfL1(1,1,1)
#endif

        dZfL1(2,1,1) = 2*Cpl9*LInt185 + Cpl9*LInt186 - 
     -   0.25D0*Cpl104 + 0.0625D0*Cpl37 + 0.25D0*Cpl38 + 
     -   LInt184*(-(0.5D0*Cpl104) + 0.125D0*Cpl37 + 0.5D0*Cpl38) + 
     -   0.125D0*Cpl8 + 0.25D0*Cpl9 + 0.25D0*(Cpl8*LInt182) + 
     -   0.5D0*(Cpl9*LInt183)

	If(dr.Eq.1) Then
	 dZfL1(2,1,1) = dZfL1(2,1,1) -
     -  (2*Cpl9*LInt185 + Cpl9*LInt186 + 0.25D0*Cpl9 + 
     -   0.5D0*(Cpl9*LInt183)) +
     -   Alfa*(0.25D0*Cpl105 - 0.25D0*Cpl106 + 
     -        0.25D0*Cpl107 - 0.25D0*Cpl108)
	Endif
	

#ifdef DEBUG
	print *, 'dZfL1(2,1,1) =', dZfL1(2,1,1)
#endif

        dZfL1(3,1,1) = -(0.16666666666666666D0*Cpl104) + 
     -   0.0625D0*Cpl37 + 0.1111111111111111D0*Cpl38 + 
     -   LInt189*(-(0.3333333333333333D0*Cpl104) + 0.125D0*Cpl37 + 
     -      0.2222222222222222D0*Cpl38) + 0.125D0*Cpl8 + 
     -   0.1111111111111111D0*Cpl9 + 0.25D0*(Cpl8*LInt187) + 
     -   0.2222222222222222D0*(Cpl9*LInt188) + 
     -   0.8888888888888888D0*(Cpl9*LInt190) + 
     -   0.4444444444444444D0*(Cpl9*LInt191)
       If(dr.Eq.1) Then
        dZfL1(3,1,1) = dZfL1(3,1,1) -
     - (0.1111111111111111D0*Cpl9+ 
     -   0.2222222222222222D0*(Cpl9*LInt188) + 
     -   0.8888888888888888D0*(Cpl9*LInt190) + 
     -   0.4444444444444444D0*(Cpl9*LInt191)) +
     -   4D0/9D0*Alfa*(0.25D0*Cpl105 - 0.25D0*Cpl106 + 
     -        0.25D0*Cpl107 - 0.25D0*Cpl108)
       Endif

#ifdef DEBUG
	print *, 'dZfL1(3,1,1) =', dZfL1(3,1,1)
#endif

        dZfL1(4,1,1) = -(0.08333333333333333D0*Cpl104) + 
     -   0.0625D0*Cpl37 + 0.027777777777777776D0*Cpl38 + 
     -   LInt193*(-(0.16666666666666666D0*Cpl104) + 0.125D0*Cpl37 + 
     -      0.05555555555555555D0*Cpl38) + 0.125D0*Cpl8 + 
     -   0.027777777777777776D0*Cpl9 + 
     -   0.05555555555555555D0*(Cpl9*LInt192) + 
     -   0.25D0*(Cpl8*LInt194) + 0.2222222222222222D0*(Cpl9*LInt195) + 
     -   0.1111111111111111D0*(Cpl9*LInt196)
       If(dr.Eq.1) Then
        dZfL1(4,1,1) = dZfL1(4,1,1) -
     - (0.027777777777777776D0*Cpl9 + 
     -   0.05555555555555555D0*(Cpl9*LInt192) + 
     -   0.2222222222222222D0*(Cpl9*LInt195) + 
     -   0.1111111111111111D0*(Cpl9*LInt196)) + 
     -   1D0/9D0*Alfa*(0.25D0*Cpl105 - 0.25D0*Cpl106 + 
     -        0.25D0*Cpl107 - 0.25D0*Cpl108)
       Endif
       
#ifdef DEBUG
	print *, 'dZfL1(4,1,1) =', dZfL1(4,1,1)
#endif

        dZfL1dimreg(3,1,1) = 
     -  0.3333333333333333D0*Cpl105 - 0.3333333333333333D0*Cpl106 + 
     -   0.3333333333333333D0*Cpl107 - 0.3333333333333333D0*Cpl108

#ifdef DEBUG
	print *, 'dZfL1dimreg(3,1,1) =', dZfL1dimreg(3,1,1)
#endif

        dZfL1dimreg(4,1,1) = 
     -  0.3333333333333333D0*Cpl105 - 0.3333333333333333D0*Cpl106 + 
     -   0.3333333333333333D0*Cpl107 - 0.3333333333333333D0*Cpl108

#ifdef DEBUG
	print *, 'dZfL1dimreg(4,1,1) =', dZfL1dimreg(4,1,1)
#endif

      
	dSW1 = -(0.5D0*(Cpl68*dMWsq1)) + 0.5D0*(Cpl69*dMZsq1)

#ifdef DEBUG
	print *, 'dSW1 =', dSW1
#endif

	dZe1 = -(0.5D0*dZAA1) - 0.5D0*(Cpl103*dZZA1)

#ifdef DEBUG
	print *, 'dZe1 =', dZe1
#endif

        dZR1 = -(Cpl77*dMWsq1) - dZAA1 + Cpl79*dZZA1 + 
     -   1.5D0*Cpl70 + 0.875D0*Cpl71 - 0.5D0*Cpl72 + 
     -   LInt70*(-2*Cpl75 - 0.25D0*Cpl73) + 
     -   LInt7*(0.875D0*Cpl73 - 0.5D0*Cpl74 - 0.5D0*Cpl75) + 
     -   LInt8*(0.0625D0*Cpl73 + 0.75D0*Cpl75) + 
     -   LInt68*(0.25D0*Cpl38 - 0.5D0*Cpl76) - 
     -   0.3333333333333333D0*Cpl9 - 2.D0*(Cpl78*dSW1) + 
     -   0.0625D0*(Cpl73*LInt6) + 0.25D0*(Cpl8*LInt67) - 
     -   0.25D0*(Cpl73*LInt69) + 0.3333333333333333D0*(Cpl9*LInt71) - 
     -   0.08333333333333333D0*(Cpl9*LInt9)

	do Gen1 = 1, 3
        dZR1 = dZR1 + 0.25D0*Cpl80(Gen1) - 
     -   0.25D0*(Cpl80(Gen1)*LInt72(Gen1)) + 
     -   0.25D0*(Cpl80(Gen1)*LInt73(Gen1))
	enddo

#ifdef DEBUG
	print *, 'dZR1 =', dZR1
#endif

      
      end


