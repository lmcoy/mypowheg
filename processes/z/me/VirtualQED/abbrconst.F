      subroutine abbr_const()
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "params.h" 

      integer resnr,delta
      double complex myA0,myB0,myB1,myB00,myB11,myB001,myB111,myDB0
      double complex myDB1,myDB00,myDB11,myBval,myCval,myDval,myEval
      double complex Prp,myDeltaIR,myDeltaUV
      double precision dZfL1photonic(4,1,1),dZfR1photonic(4,1,1),dr
      
      integer myBget,myCget,myDget,myEget,Gen1

#include "amps/looptools.h" 
#include "amps/ampcommon.h"
#include "amps/ampvars.h" 
#include "amps/ampcom.h" 

	Cpl1 = Alfa*Pi
	Cpl10 = Alfa/Pi*MWr**2
	Cpl100 = Alfa2*(MW2c*SW2c)/CW2c
	Cpl101 = Alfa2/CW2c**2
	Cpl102 = Alfa2*MW2c/CW2c**3
	Cpl103 = Alfa2*(MW2c*SW2c**2)/CW2c**2
	Cpl107 = Alfa2*SW2c/CW2c**2
	Cpl11 = Alfa/Pi*MW2c/SW2c
	Cpl111 = Alfa2/SW2c
	Cpl112 = Alfa2*1/(CW2c**2*SW2c)
	Cpl113 = Alfa2*MW2c/(CW2c*SW2c)
	Cpl114 = Alfa2*(MZ2c*SW2c)/CW2c**2
	Cpl115 = Alfa2*MZ2c/CW2c**2
	Cpl116 = Alfa2*MZ2c/(CW2c**2*SW2c)
	Cpl117 = Alfa2*MZ2c/CW2c
	Cpl118 = Alfa2*MZ2c/(CW2c*SW2c)
	Cpl119 = Alfa2*MW2c/CW2c
	Cpl12 = Alfa/Pi*MW2r/SW2c
	Cpl120 = Alfa2*MW2c/SW2c
	Cpl121 = Alfa2*MW2c/(CW2c**3*SW2c)
	Cpl122 = Alfa2*(MW2c*SW2c)/CW2c**2
	Cpl125 = Alfa2/SW2c**2
	Cpl126 = Alfa2*1/(CW2c**2*SW2c**2)
	Cpl127 = Alfa2*1/(CW2c*SW2c**2)
	Cpl128 = Alfa2*MW2c/(CW2c*SW2c**2)
	Cpl129 = Alfa2*MZ2c/(CW2c**2*SW2c**2)
	Cpl13 = Alfa/Pi*MWr**2/SW2c
	Cpl130 = Alfa2*MW2c/(CW2c**3*SW2c**2)
	Cpl131 = Alfa2*MW2c/CW2c**2
	Cpl132 = Alfa2*MW2c/SW2c**2
	Cpl14 = Alfa/Pi/SW2c
	Cpl15 = Alfa/Pi
	Cpl155 = Alfa*Pi*(MZ2c*SW2c)/CW2c
	Cpl156 = Alfa*Pi*SWc/CWc
	Cpl157 = Alfa*Pi*SWc/CW2c**2
	Cpl158 = Alfa*Pi*(MZ2c*SWc)/CWc
	Cpl159 = Alfa*Pi*MZ2c/CW2c
	Cpl16 = Alfa/Pi*CW2c/SW2c
	Cpl160 = Alfa*Pi*1/(CWc*SWc)
	Cpl161 = Alfa*Pi*SW2c/(CWc*SWc)
	Cpl162 = Alfa*Pi*1/(CW2c**2*SWc)
	Cpl163 = Alfa*Pi*1/(CW2c*SWc)
	Cpl164 = Alfa*Pi*SW2c/(CW2c**2*SWc)
	Cpl165 = Alfa*Pi*MZ2c/(CWc*SWc)
	Cpl166 = Alfa*Pi*(MZ2c*SW2c)/(CWc*SWc)
	Cpl167 = Alfa*Pi*MZ2c/(CW2c*SW2c)
	Cpl168 = Alfa*Pi*1/(CW2c*SW2c*SWc)
	Cpl169 = Alfa/Pi*(CWc*MZr**2)/(MZ2r*SWc)
	Cpl17 = Alfa/Pi*MW2r
	Cpl170 = Alfa/Pi*CW2c/(CWc*MZ2r*SWc)
	Cpl171 = Alfa/Pi*SW2c/(CWc*MZ2r*SWc)
	Cpl172 = Alfa/Pi*MW2c/(CWc*MZ2r*SWc)
	Cpl173 = MZ2c/MZ2r
	Cpl174 = Alfa/Pi*MZr**2/(CWc*MZ2r*SWc)
	Cpl175 = Alfa/Pi*(MZr**2*SW2c)/(CWc*MZ2r*SWc)
	Cpl178 = SWc/CWc
	Cpl179 = Alfa/Pi/CW2c
	Cpl18 = Alfa/Pi*(CW2c*MW2c)/SW2c
	Cpl180 = 1/Pi*Log(lambda**2)
	Cpl181 = 1/Pi*Log(mureg)
	Cpl182 = 1/Pi*myDeltaIR()
	Cpl183 = 1/Pi*myDeltaUV()
	Cpl184 = Alfa*Pi*EL
	Cpl185 = Alfa*Pi*(EL*SW2c)/CW2c
	Cpl186 = Alfa*Pi*EL/CW2c
	Cpl187 = Alfa*Pi*EL/(CW2c*SW2c)
	Cpl19 = Alfa/Pi*(CW2c*MW2r)/SW2c
	Cpl2 = Alfa*Pi*SW2c/CW2c
	Cpl20 = Alfa/Pi*(CW2c*MWr**2)/SW2c
	Cpl21 = Alfa/Pi*(CW2c*MZ2c)/SW2c
	Cpl22 = Alfa/Pi*(MW2c*SW2c)/CW2c
	Cpl23 = Alfa/Pi*MW2c**2
	Cpl24 = Alfa/Pi*(MW2c*MW2r)
	Cpl25 = Alfa/Pi*(MW2c*MWr**2)
	Cpl26 = Alfa/Pi*(MW2r*MWr**2)
	Cpl27 = Alfa/Pi*MW2c**2/SW2c
	Cpl28 = Alfa/Pi*(MW2c*MW2r)/SW2c
	Cpl29 = Alfa/Pi*(CW2c*MW2c*MWr**2)/SW2c
	Cpl3 = (Alfa*Pi)/CW2c
	Cpl30 = Alfa/Pi*(CW2c*MW2r*MWr**2)/SW2c
	Cpl31 = Alfa/Pi*(CW2c*MW2c*MZ2c)/SW2c
	Cpl32 = Alfa/Pi*(CW2c*MW2r*MZ2c)/SW2c
	Cpl33 = Alfa/Pi*(MW2c**2*SW2c)/CW2c
	Cpl34 = Alfa/Pi*(MW2c*MW2r*SW2c)/CW2c
	Cpl38 = Alfa/Pi*(MW2c*MWr**2)/SW2c
	Cpl39 = Alfa/Pi*(MW2r*MWr**2)/SW2c
	Cpl4 = Alfa*Pi*1/(CW2c*SW2c)
	Cpl40 = Alfa/Pi*(CW2c*MZ2r)/SW2c
	Cpl41 = Alfa/Pi*MZr**2/(CW2c*SW2c)
	Cpl42 = Alfa/Pi*(CW2c*MZr**2)/SW2c
	Cpl43 = Alfa/Pi*1/(CW2c*SW2c)
	Cpl44 = Alfa/Pi*SW2c/CW2c
	Cpl45 = Alfa/Pi*MW2c/(CW2c**2*SW2c)
	Cpl46 = Alfa/Pi*MZ2c/(CW2c*SW2c)
	Cpl47 = Alfa/Pi*MZ2r/(CW2c*SW2c)
	Cpl48 = Alfa/Pi*(MW2c*MZ2c)/(CW2c**2*SW2c)
	Cpl49 = Alfa/Pi*(MW2c*MZ2r)/(CW2c**2*SW2c)
	Cpl5 = Alfa
	Cpl50 = Alfa/Pi*(CW2c*MW2c*MZ2r)/SW2c
	Cpl51 = Alfa/Pi*(CW2c*MZ2c*MZr**2)/SW2c
	Cpl52 = Alfa/Pi*(CW2c*MZ2r*MZr**2)/SW2c
	Cpl53 = Alfa/Pi*(MW2c*MZ2c*SW2c)/CW2c
	Cpl54 = Alfa/Pi*(MW2c*MZ2r*SW2c)/CW2c
	Cpl55 = Alfa/Pi*MZ2c
	Cpl56 = Alfa/Pi*MZ2r
	Cpl57 = Alfa/Pi*(MZ2c*SW2c)/CW2c
	Cpl58 = Alfa/Pi*(MZ2r*SW2c)/CW2c
	Cpl59 = Alfa/Pi*(MZ2c*MZr**2)/(CW2c*SW2c)
	Cpl6 = Alfa*SW2c/CW2c
	Cpl60 = Alfa/Pi*(MZ2r*MZr**2)/(CW2c*SW2c)
	Cpl61 = Alfa/Pi*MZr**2/CW2c
	Cpl62 = Alfa/Pi*(MZr**2*SW2c)/CW2c
	Cpl66 = Alfa/Pi*MZ2c/CW2c
	Cpl67 = Alfa/Pi*MZ2r/CW2c
	Cpl7 = Alfa/CW2c
	Cpl70 = Alfa/Pi*(MZ2c*MZr**2)/CW2c
	Cpl71 = Alfa/Pi*(MZ2r*MZr**2)/CW2c
	Cpl72 = Alfa/Pi*(MZ2c*MZr**2*SW2c)/CW2c
	Cpl73 = Alfa/Pi*(MZ2r*MZr**2*SW2c)/CW2c
	Cpl74 = CW2c/(MW2c*SWc)
	Cpl75 = CW2c/(MZ2c*SWc)
	Cpl76 = Alfa*1/(pi*SW2c)
	Cpl77 = Alfa*Log(CW2c)/(pi*SW2c**2)
	Cpl78 = Alfa*Log(CW2c)/(pi*SW2c)
	Cpl79 = Alfa/Pi*1/(MW2c*SW2c)
	Cpl8 = Alfa*1/(CW2c*SW2c)
	Cpl80 = Alfa/Pi/MW2c
	Cpl81 = Alfa/Pi*CW2c/(MW2c*SW2c)
	Cpl82 = Alfa/Pi*(CW2c*MZ2c)/(MW2c*SW2c)
	Cpl83 = 1/MW2c
	Cpl84 = 1/SWc
	Cpl85 = CWc/SWc
	Cpl87 = Alfa/Pi*CW2c/(CWc*MZ2c*SWc)
	Cpl88 = Alfa/Pi*SW2c/(CWc*MZ2c*SWc)
	Cpl89 = Alfa/Pi*MW2c/(CWc*MZ2c*SWc)
	Cpl9 = Alfa/Pi*MW2c
	Cpl92 = Alfa2
	Cpl93 = Alfa2*SW2c/CW2c
	Cpl94 = Alfa2/CW2c
	Cpl95 = Alfa2*1/(CW2c*SW2c)
	Cpl96 = Alfa2*SW2c**2/CW2c**2
	Cpl97 = Alfa2*(MZ2c*SW2c**2)/CW2c**2
	Cpl98 = Alfa2*(MZ2c*SW2c)/CW2c
	Cpl99 = Alfa2*MW2c
      

	do Gen1=1,3
	Cpl104(Gen1) = Alfa2*MT2*delta(Gen1,3)
	Cpl105(Gen1) = Alfa2*MT2*delta(Gen1,3)/CW2c
	Cpl106(Gen1) = Alfa2*MT2*(SW2c*delta(Gen1,3))/CW2c
	Cpl108(Gen1) = Alfa2*MT2*delta(Gen1,3)/CW2c**2
	Cpl109(Gen1) = Alfa2*MT2*(SW2c*delta(Gen1,3))/CW2c**2
	Cpl110(Gen1) = Alfa2*MT2*(SW2c**2*delta(Gen1,3))/CW2c**2
	Cpl123(Gen1) = Alfa2*MT2*delta(Gen1,3)/(CW2c*SW2c)
	Cpl124(Gen1) = Alfa2*MT2*delta(Gen1,3)/(CW2c**2*SW2c)
	Cpl133(Gen1) = Alfa2*MT2*delta(Gen1,3)/(CW2c**2*SW2c**2)
	Cpl176(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/(CWc*MZ2r*SWc)
	Cpl177(Gen1) = (Alfa*MT2)/Pi*(SW2c*delta(Gen1,3))/(CWc*MZ2r*SWc)
	Cpl35(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/SW2c
	Cpl36(Gen1) = (Alfa*MT2)/Pi*(MW2c*delta(Gen1,3))/SW2c
	Cpl37(Gen1) = (Alfa*MT2)/Pi*(MW2r*delta(Gen1,3))/SW2c
	Cpl63(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/CW2c
	Cpl64(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/(CW2c*SW2c)
	Cpl65(Gen1) = (Alfa*MT2)/Pi*(SW2c*delta(Gen1,3))/CW2c
	Cpl68(Gen1) = (Alfa*MT2)/Pi*(MZ2c*delta(Gen1,3))/(CW2c*SW2c)
	Cpl69(Gen1) = (Alfa*MT2)/Pi*(MZ2r*delta(Gen1,3))/(CW2c*SW2c)
	Cpl86(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/(MW2c*SW2c)
	Cpl90(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/(CWc*MZ2c*SWc)
	Cpl91(Gen1) = (Alfa*MT2)/Pi*(SW2c*delta(Gen1,3))/(CWc*MZ2c*SWc)
      
        enddo
      
      
      
	IInt11 = myBget(0.D0,MW2r,MZ2r,0)
	LInt73 = myBval(bb0,IInt11)
	IInt12 = myBget(0.D0,MZ2r,MW2r,2)
	LInt75 = myBval(bb00,IInt12)
      

	do Gen1=1,3
	IInt23(Gen1) = myBget(0.D0,MQD2(Gen1),MQU2(Gen1),0)
	LInt77(Gen1) = myBval(bb0,IInt23(Gen1))
	IInt24(Gen1) = myBget(0.D0,MQU2(Gen1),MQD2(Gen1),1)
	LInt78(Gen1) = myBval(bb1,IInt24(Gen1))
        enddo

	IInt1 = myBget(0.D0,0.D0,MW2r,1)
	LInt14 = myBval(bb0,IInt1)
	LInt76 = myBval(bb1,IInt1)
	IInt10 = myBget(0.D0,MH2,MW2r,2)
	LInt72 = myBval(bb0,IInt10)
	LInt74 = myBval(bb00,IInt10)
	IInt4 = myBget(MWr**2,MW2r,MZ2r,0)
	LInt17 = myBval(bb0,IInt4)
	LInt24 = myBval(dbb0,IInt4)
      

	do Gen1=1,3
	IInt13(Gen1) = myBget(MWr**2,MQD2(Gen1),MQU2(Gen1),0)
	LInt30(Gen1) = myBval(bb0,IInt13(Gen1))
	LInt33(Gen1) = myBval(dbb0,IInt13(Gen1))
	IInt20(Gen1) = myBget(0.D0,MLE2(Gen1),MLE2(Gen1),2)
	LInt66(Gen1) = myBval(bb1,IInt20(Gen1))
	LInt69(Gen1) = myBval(dbb00,IInt20(Gen1))
	IInt21(Gen1) = myBget(0.D0,MQD2(Gen1),MQD2(Gen1),2)
	LInt67(Gen1) = myBval(bb1,IInt21(Gen1))
	LInt70(Gen1) = myBval(dbb00,IInt21(Gen1))
        enddo

	IInt8 = myBget(MZr**2,0.D0,0.D0,2)
	LInt42 = myBval(bb1,IInt8)
	LInt46 = myBval(dbb00,IInt8)
	LInt49 = myBval(dbb1,IInt8)
      

	do Gen1=1,3
	IInt14(Gen1) = myBget(MWr**2,0.D0,MLE2(Gen1),2)
	LInt31(Gen1) = myBval(bb1,IInt14(Gen1))
	LInt34(Gen1) = myBval(dbb00,IInt14(Gen1))
	LInt36(Gen1) = myBval(dbb1,IInt14(Gen1))
	IInt15(Gen1) = myBget(MWr**2,MQU2(Gen1),MQD2(Gen1),2)
	LInt32(Gen1) = myBval(bb1,IInt15(Gen1))
	LInt35(Gen1) = myBval(dbb00,IInt15(Gen1))
	LInt37(Gen1) = myBval(dbb1,IInt15(Gen1))
	IInt16(Gen1) = myBget(0.D0,0.D0,MQU2(Gen1),0)
	LInt355(Gen1) = myBval(bb0,IInt16(Gen1))
	LInt359(Gen1) = myBval(bb0,IInt16(Gen1))
	LInt51(Gen1) = myBval(bb0,IInt16(Gen1))
	IInt18(Gen1) = myBget(MZr**2,MLE2(Gen1),MLE2(Gen1),2)
	LInt53(Gen1) = myBval(bb1,IInt18(Gen1))
	LInt57(Gen1) = myBval(dbb00,IInt18(Gen1))
	LInt60(Gen1) = myBval(dbb1,IInt18(Gen1))
	IInt19(Gen1) = myBget(MZr**2,MQD2(Gen1),MQD2(Gen1),2)
	LInt54(Gen1) = myBval(bb1,IInt19(Gen1))
	LInt58(Gen1) = myBval(dbb00,IInt19(Gen1))
	LInt61(Gen1) = myBval(dbb1,IInt19(Gen1))
	IInt22(Gen1) = myBget(0.D0,MQU2(Gen1),MQU2(Gen1),2)
	LInt68(Gen1) = myBval(bb1,IInt22(Gen1))
	LInt71(Gen1) = myBval(dbb00,IInt22(Gen1))
	LInt80(Gen1) = myBval(bb0,IInt22(Gen1))
        enddo

	IInt3 = myBget(MWr**2,MH2,MW2r,2)
	LInt16 = myBval(bb0,IInt3)
	LInt18 = myBval(bb00,IInt3)
	LInt23 = myBval(dbb0,IInt3)
	LInt26 = myBval(dbb00,IInt3)
	IInt5 = myBget(MWr**2,MZ2r,MW2r,2)
	LInt19 = myBval(bb00,IInt5)
	LInt21 = myBval(bb1,IInt5)
	LInt27 = myBval(dbb00,IInt5)
	LInt29 = myBval(dbb1,IInt5)
	IInt6 = myBget(MZr**2,MH2,MZ2r,2)
	LInt38 = myBval(bb0,IInt6)
	LInt40 = myBval(bb00,IInt6)
	LInt44 = myBval(dbb0,IInt6)
	LInt47 = myBval(dbb00,IInt6)
	IInt9 = myBget(0.D0,MW2r,MW2r,2)
	LInt63 = myBval(bb0,IInt9)
	LInt64 = myBval(bb1,IInt9)
	LInt65 = myBval(dbb00,IInt9)
	LInt79 = myBval(bb00,IInt9)
	IInt2 = myBget(MWr**2,0.D0,MW2r,2)
	LInt15 = myBval(bb0,IInt2)
	LInt20 = myBval(bb1,IInt2)
	LInt22 = myBval(dbb0,IInt2)
	LInt25 = myBval(dbb00,IInt2)
	LInt28 = myBval(dbb1,IInt2)
      

	do Gen1=1,3
	IInt17(Gen1) = myBget(MZr**2,MQU2(Gen1),MQU2(Gen1),2)
	LInt52(Gen1) = myBval(bb0,IInt17(Gen1))
	LInt55(Gen1) = myBval(bb1,IInt17(Gen1))
	LInt56(Gen1) = myBval(dbb0,IInt17(Gen1))
	LInt59(Gen1) = myBval(dbb00,IInt17(Gen1))
	LInt62(Gen1) = myBval(dbb1,IInt17(Gen1))
        enddo

	IInt7 = myBget(MZr**2,MW2r,MW2r,2)
	LInt39 = myBval(bb0,IInt7)
	LInt41 = myBval(bb00,IInt7)
	LInt43 = myBval(bb1,IInt7)
	LInt45 = myBval(dbb0,IInt7)
	LInt48 = myBval(dbb00,IInt7)
	LInt50 = myBval(dbb1,IInt7)
	LInt11 = myA0(MH2)
	LInt12 = myA0(MW2r)
	LInt13 = myA0(MZ2r)
	LInt352 = myA0(MW2r)
	LInt353 = myA0(MH2)
	LInt354 = myA0(MZ2r)
	LInt356 = myA0(MW2r)
	LInt357 = myA0(MH2)
	LInt358 = myA0(MZ2r)
	LInt360 = myB1(ME2,0.D0,MW2r)
	LInt361 = myB1(ME2,ME2,0.D0)
	LInt362 = myB1(ME2,ME2,MZ2r)
	LInt363 = ME2*myDB0(ME2,0.D0,ME2)
	LInt364 = ME2*myDB1(ME2,ME2,0.D0)
	LInt365 = myB1(MU2,MD2,MW2r)
	LInt366 = myB1(MU2,MU2,0.D0)
	LInt367 = myB1(MU2,MU2,MZ2r)
	LInt368 = MU2*myDB0(MU2,0.D0,MU2)
	LInt369 = MU2*myDB1(MU2,MU2,0.D0)
	LInt370 = myB1(MD2,MD2,0.D0)
	LInt371 = myB1(MD2,MD2,MZ2r)
	LInt372 = myB1(MD2,MU2,MW2r)
	LInt373 = MD2*myDB0(MD2,0.D0,MD2)
	LInt374 = MD2*myDB1(MD2,MD2,0.D0)
      
        dMWsq1 = (2*Cpl17 - 2*Cpl9)*LInt25 + 
     -   0.1111111111111111D0*Cpl10 - 0.16666666666666666D0*Cpl11 + 
     -   0.16666666666666666D0*Cpl12 + 
     -   LInt26*(-(0.25D0*Cpl11) + 0.25D0*Cpl12) + 
     -   LInt27*(-2*Cpl18 + 2*Cpl19 - 0.25D0*Cpl11 + 0.25D0*Cpl12) - 
     -   0.16666666666666666D0*Cpl13 + 
     -   LInt19*(-2*Cpl16 - 0.25D0*Cpl14) + 
     -   LInt12*(0.875D0*Cpl14 - 0.5D0*Cpl15 - 0.5D0*Cpl16) + 
     -   LInt13*(0.0625D0*Cpl14 + 0.75D0*Cpl16) + 
     -   LInt15*(-Cpl9 - 1.25D0*Cpl10 + 1.25D0*Cpl17) + 
     -   LInt21*(-(0.5D0*Cpl18) + 0.5D0*Cpl19 - 0.5D0*Cpl20) + 
     -   LInt17*(-(1.25D0*Cpl18) + 1.25D0*Cpl19 - 1.25D0*Cpl20 - 
     -      0.5D0*Cpl21 + 0.25D0*Cpl22) + 
     -   LInt28*(-(0.5D0*Cpl25) + 0.5D0*Cpl26) + 
     -   LInt22*(0.25D0*Cpl23 - 0.25D0*Cpl24 - 1.25D0*Cpl25 + 
     -      1.25D0*Cpl26) + LInt23*(0.25D0*Cpl27 - 0.25D0*Cpl28) + 
     -   LInt29*(-(0.5D0*Cpl29) + 0.5D0*Cpl30) + 
     -   LInt24*(-(1.25D0*Cpl29) + 1.25D0*Cpl30 - 0.5D0*Cpl31 + 
     -      0.5D0*Cpl32 + 0.25D0*Cpl33 - 0.25D0*Cpl34) + 
     -   LInt20*(-(0.8333333333333334D0*Cpl10) + 0.5D0*Cpl17 - 
     -      0.16666666666666666D0*Cpl9) - 0.3333333333333333D0*Cpl9 + 
     -   0.0625D0*(Cpl14*LInt11) - 
     -   0.3333333333333333D0*(Cpl9*LInt14) + 0.25D0*(Cpl11*LInt16) - 
     -   0.25D0*(Cpl14*LInt18)

	do Gen1 = 1, 3
        dMWsq1 = dMWsq1 - 0.1111111111111111D0*Cpl13 + 
     -   0.25D0*Cpl35(Gen1) - 0.25D0*(Cpl35(Gen1)*LInt30(Gen1)) + 
     -   (-(0.25D0*Cpl11) + 0.25D0*Cpl12 - 0.16666666666666666D0*Cpl13)*
     -    LInt31(Gen1) + (-(0.75D0*Cpl11) + 0.75D0*Cpl12 - 0.5D0*Cpl13 + 
     -      0.25D0*Cpl35(Gen1))*LInt32(Gen1) + 
     -   (-(0.75D0*Cpl36(Gen1)) + 0.75D0*Cpl37(Gen1))*LInt33(Gen1) + 
     -   (0.5D0*Cpl11 - 0.5D0*Cpl12)*LInt34(Gen1) + 
     -   (1.5D0*Cpl11 - 1.5D0*Cpl12)*LInt35(Gen1) + 
     -   (-(0.25D0*Cpl38) + 0.25D0*Cpl39)*LInt36(Gen1) + 
     -   (-(0.75D0*Cpl38) + 0.75D0*Cpl39)*LInt37(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dMWsq1 =', dMWsq1
#endif

        dMZsq1 = -(0.16666666666666666D0*Cpl21) + 
     -   0.16666666666666666D0*Cpl40 - 0.041666666666666664D0*Cpl41 - 
     -   0.16666666666666666D0*Cpl42 + 
     -   LInt43*(-(0.5D0*Cpl21) + 0.5D0*Cpl40 - 0.5D0*Cpl42) + 
     -   LInt39*(-(0.5D0*Cpl18) - 1.25D0*Cpl21 + 0.5D0*Cpl22 + 
     -      1.25D0*Cpl40 - 1.25D0*Cpl42) + 
     -   LInt12*(-(0.25D0*Cpl15) + 1.125D0*Cpl16 + 0.125D0*Cpl44) + 
     -   LInt41*(0.5D0*Cpl15 - 2.25D0*Cpl16 - 0.25D0*Cpl44) + 
     -   LInt47*(-(0.25D0*Cpl46) + 0.25D0*Cpl47) + 
     -   LInt42*(-(0.25D0*Cpl41) - 0.375D0*Cpl46 + 0.375D0*Cpl47) + 
     -   LInt46*(0.75D0*Cpl46 - 0.75D0*Cpl47) + 
     -   LInt44*(0.25D0*Cpl48 - 0.25D0*Cpl49) + 
     -   LInt50*(-(0.5D0*Cpl51) + 0.5D0*Cpl52) + 
     -   LInt45*(-(0.5D0*Cpl31) + 0.5D0*Cpl50 - 1.25D0*Cpl51 + 
     -      1.25D0*Cpl52 + 0.5D0*Cpl53 - 0.5D0*Cpl54) + 
     -   LInt48*(-(2.25D0*Cpl21) + 2.25D0*Cpl40 + 0.5D0*Cpl55 - 
     -      0.5D0*Cpl56 - 0.25D0*Cpl57 + 0.25D0*Cpl58) + 
     -   LInt49*(-(0.375D0*Cpl59) + 0.375D0*Cpl60) + 
     -   0.0625D0*(Cpl43*LInt11) + 0.0625D0*(Cpl43*LInt13) + 
     -   0.25D0*(Cpl45*LInt38) - 0.25D0*(Cpl43*LInt40)

	do Gen1 = 1, 3
        dMZsq1 = dMZsq1 - 0.09722222222222222D0*Cpl41 + 
     -   0.2222222222222222D0*Cpl61 - 0.2962962962962963D0*Cpl62 - 
     -   0.6666666666666666D0*Cpl63(Gen1) + 0.25D0*Cpl64(Gen1) + 
     -   0.8888888888888888D0*Cpl65(Gen1) + 
     -   (0.6666666666666666D0*Cpl63(Gen1) - 0.25D0*Cpl64(Gen1) - 
     -      0.8888888888888888D0*Cpl65(Gen1))*LInt51(Gen1) + 
     -   (-(0.6666666666666666D0*Cpl63(Gen1)) - 0.125D0*Cpl64(Gen1) + 
     -      0.8888888888888888D0*Cpl65(Gen1))*LInt52(Gen1) + 
     -   (-Cpl57 + Cpl58 - 0.08333333333333333D0*Cpl41 - 
     -      0.125D0*Cpl46 + 0.125D0*Cpl47 + 0.3333333333333333D0*Cpl61 - 
     -      0.6666666666666666D0*Cpl62 + 0.5D0*Cpl66 - 0.5D0*Cpl67)*
     -    LInt53(Gen1) + (-(0.25D0*Cpl41) - 0.375D0*Cpl46 + 
     -      0.375D0*Cpl47 - 0.3333333333333333D0*Cpl57 + 
     -      0.3333333333333333D0*Cpl58 + 0.3333333333333333D0*Cpl61 - 
     -      0.2222222222222222D0*Cpl62 + 0.5D0*Cpl66 - 0.5D0*Cpl67)*
     -    LInt54(Gen1) + (Cpl66 - Cpl67 - 0.25D0*Cpl41 - 
     -      0.375D0*Cpl46 + 0.375D0*Cpl47 - 1.3333333333333333D0*Cpl57 + 
     -      1.3333333333333333D0*Cpl58 + 0.6666666666666666D0*Cpl61 - 
     -      0.8888888888888888D0*Cpl62)*LInt55(Gen1) + 
     -   (-(0.375D0*Cpl68(Gen1)) + 0.375D0*Cpl69(Gen1))*LInt56(Gen1) + 
     -   (2*Cpl57 - 2*Cpl58 - Cpl66 + Cpl67 + 0.25D0*Cpl46 - 
     -      0.25D0*Cpl47)*LInt57(Gen1) + 
     -   (-Cpl66 + Cpl67 + 0.75D0*Cpl46 - 0.75D0*Cpl47 + 
     -      0.6666666666666666D0*Cpl57 - 0.6666666666666666D0*Cpl58)*
     -    LInt58(Gen1) + (-2*Cpl66 + 2*Cpl67 + 0.75D0*Cpl46 - 
     -      0.75D0*Cpl47 + 2.6666666666666665D0*Cpl57 - 
     -      2.6666666666666665D0*Cpl58)*LInt59(Gen1) + 
     -   (-Cpl72 + Cpl73 - 0.125D0*Cpl59 + 0.125D0*Cpl60 + 0.5D0*Cpl70 - 
     -      0.5D0*Cpl71)*LInt60(Gen1) + 
     -   (-(0.375D0*Cpl59) + 0.375D0*Cpl60 + 0.5D0*Cpl70 - 0.5D0*Cpl71 - 
     -      0.3333333333333333D0*Cpl72 + 0.3333333333333333D0*Cpl73)*
     -    LInt61(Gen1) + (Cpl70 - Cpl71 - 0.375D0*Cpl59 + 
     -      0.375D0*Cpl60 - 1.3333333333333333D0*Cpl72 + 
     -      1.3333333333333333D0*Cpl73)*LInt62(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dMZsq1 =', dMZsq1
#endif

        dZAA1 = 3*Cpl15*LInt65 + 0.16666666666666666D0*Cpl15 + 
     -   1.25D0*(Cpl15*LInt63) + 0.5D0*(Cpl15*LInt64)

	do Gen1 = 1, 3
        dZAA1 = dZAA1 + 0.3333333333333333D0*
     -    (Cpl15*LInt67(Gen1)) + 
     -   1.3333333333333333D0*(Cpl15*LInt68(Gen1)) - 
     -   0.6666666666666666D0*(Cpl15*LInt70(Gen1)) - 
     -   2.6666666666666665D0*(Cpl15*LInt71(Gen1)) + 
     -   Cpl15*LInt66(Gen1) - 2*Cpl15*LInt69(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dZAA1 =', dZAA1
#endif

        dZZA1 = Cpl89*LInt63 + (5*Cpl87 - Cpl88)*LInt79 + 
     -   LInt12*(-(2.5D0*Cpl87) + 0.5D0*Cpl88)

	do Gen1 = 1, 3
        dZZA1 = dZZA1 - 0.6666666666666666D0*Cpl90(Gen1) + 
     -   1.7777777777777777D0*Cpl91(Gen1) + 
     -   (0.6666666666666666D0*Cpl90(Gen1) - 
     -      1.7777777777777777D0*Cpl91(Gen1))*LInt51(Gen1) + 
     -   (-(0.6666666666666666D0*Cpl90(Gen1)) + 
     -      1.7777777777777777D0*Cpl91(Gen1))*LInt80(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dZZA1 =', dZZA1
#endif

        dZZZ1 = 0.16666666666666666D0*Cpl16 + 
     -   LInt45*(0.5D0*Cpl18 - 0.5D0*Cpl22 + 1.25D0*Cpl42) + 
     -   LInt48*(-(0.5D0*Cpl15) + 2.25D0*Cpl16 + 0.25D0*Cpl44) + 
     -   1.25D0*(Cpl16*LInt39) + 0.375D0*(Cpl43*LInt42) + 
     -   0.5D0*(Cpl16*LInt43) - 0.25D0*(Cpl45*LInt44) - 
     -   0.75D0*(Cpl43*LInt46) + 0.25D0*(Cpl43*LInt47) + 
     -   0.375D0*(Cpl41*LInt49) + 0.5D0*(Cpl42*LInt50)

	do Gen1 = 1, 3
        dZZZ1 = dZZZ1 + 0.375D0*(Cpl64(Gen1)*LInt56(Gen1)) + 
     -   (Cpl44 - 0.5D0*Cpl179 + 0.125D0*Cpl43)*LInt53(Gen1) + 
     -   (-(0.5D0*Cpl179) + 0.375D0*Cpl43 + 0.3333333333333333D0*Cpl44)*
     -    LInt54(Gen1) + (-Cpl179 + 0.375D0*Cpl43 + 
     -      1.3333333333333333D0*Cpl44)*LInt55(Gen1) + 
     -   (Cpl179 - 2*Cpl44 - 0.25D0*Cpl43)*LInt57(Gen1) + 
     -   (Cpl179 - 0.75D0*Cpl43 - 0.6666666666666666D0*Cpl44)*
     -    LInt58(Gen1) + (2*Cpl179 - 0.75D0*Cpl43 - 
     -      2.6666666666666665D0*Cpl44)*LInt59(Gen1) + 
     -   (Cpl62 + 0.125D0*Cpl41 - 0.5D0*Cpl61)*LInt60(Gen1) + 
     -   (0.375D0*Cpl41 - 0.5D0*Cpl61 + 0.3333333333333333D0*Cpl62)*
     -    LInt61(Gen1) + (-Cpl61 + 0.375D0*Cpl41 + 
     -      1.3333333333333333D0*Cpl62)*LInt62(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dZZZ1 =', dZZZ1
#endif

        dZfL1(2,1,1) = 2*Cpl15*LInt363 + Cpl15*LInt364 + 
     -   0.125D0*Cpl14 + 0.25D0*Cpl15 - 0.25D0*Cpl179 + 0.0625D0*Cpl43 + 
     -   0.25D0*Cpl44 + LInt362*
     -    (-(0.5D0*Cpl179) + 0.125D0*Cpl43 + 0.5D0*Cpl44) + 
     -   0.25D0*(Cpl14*LInt360) + 0.5D0*(Cpl15*LInt361)
     -   -
     -   (2*Cpl15*LInt363 + Cpl15*LInt364 + 
     -    0.25D0*Cpl15 + 0.5D0*(Cpl15*LInt361))

#ifdef DEBUG
	print *, 'dZfL1(2,1,1) =', dZfL1(2,1,1)
#endif

        dZfL1(3,1,1) = 0.125D0*Cpl14 + 0.1111111111111111D0*Cpl15 - 
     -   0.16666666666666666D0*Cpl179 + 0.0625D0*Cpl43 + 
     -   0.1111111111111111D0*Cpl44 + 
     -   LInt367*(-(0.3333333333333333D0*Cpl179) + 0.125D0*Cpl43 + 
     -      0.2222222222222222D0*Cpl44) + 0.25D0*(Cpl14*LInt365) + 
     -   0.2222222222222222D0*(Cpl15*LInt366) + 
     -   0.8888888888888888D0*(Cpl15*LInt368) + 
     -   0.4444444444444444D0*(Cpl15*LInt369)
     -   -
     -   (0.1111111111111111D0*Cpl15 + 
     -   0.2222222222222222D0*(Cpl15*LInt366) +  
     -   0.8888888888888888D0*(Cpl15*LInt368) + 
     -   0.4444444444444444D0*(Cpl15*LInt369))

#ifdef DEBUG
	print *, 'dZfL1(3,1,1) =', dZfL1(3,1,1)
#endif

        dZfL1(4,1,1) = 0.125D0*Cpl14 + 
     -   0.027777777777777776D0*Cpl15 - 0.08333333333333333D0*Cpl179 + 
     -   0.0625D0*Cpl43 + 0.027777777777777776D0*Cpl44 + 
     -   LInt371*(-(0.16666666666666666D0*Cpl179) + 0.125D0*Cpl43 + 
     -      0.05555555555555555D0*Cpl44) + 
     -   0.05555555555555555D0*(Cpl15*LInt370) + 
     -   0.25D0*(Cpl14*LInt372) + 
     -   0.2222222222222222D0*(Cpl15*LInt373) + 
     -   0.1111111111111111D0*(Cpl15*LInt374)
     -   -
     -   (0.027777777777777776D0*Cpl15 + 
     -   0.05555555555555555D0*(Cpl15*LInt370) + 
     -   0.2222222222222222D0*(Cpl15*LInt373) + 
     -   0.1111111111111111D0*(Cpl15*LInt374))

#ifdef DEBUG
	print *, 'dZfL1(4,1,1) =', dZfL1(4,1,1)
#endif

        dZfL1dimreg(3,1,1) = 
     -  0.3333333333333333D0*Cpl180 - 0.3333333333333333D0*Cpl181 + 
     -   0.3333333333333333D0*Cpl182 - 0.3333333333333333D0*Cpl183

#ifdef DEBUG
	print *, 'dZfL1dimreg(3,1,1) =', dZfL1dimreg(3,1,1)
#endif

        dZfL1dimreg(4,1,1) = 
     -  0.3333333333333333D0*Cpl180 - 0.3333333333333333D0*Cpl181 + 
     -   0.3333333333333333D0*Cpl182 - 0.3333333333333333D0*Cpl183

#ifdef DEBUG
	print *, 'dZfL1dimreg(4,1,1) =', dZfL1dimreg(4,1,1)
#endif

        dZfR1(2,1,1) =  0.5D0*(Cpl44*LInt362)+0.25D0*Cpl44

#ifdef DEBUG
	print *, 'dZfR1(2,1,1) =', dZfR1(2,1,1)
#endif

        dZfR1(3,1,1) = 0.1111111111111111D0*Cpl44 + 
     -   0.2222222222222222D0*(Cpl44*LInt367)

#ifdef DEBUG
	print *, 'dZfR1(3,1,1) =', dZfR1(3,1,1)
#endif

        dZfR1(4,1,1) = 0.027777777777777776D0*Cpl44 + 
     -   0.05555555555555555D0*(Cpl44*LInt371)

#ifdef DEBUG
	print *, 'dZfR1(4,1,1) =', dZfR1(4,1,1)
#endif

        dZfR1dimreg(3,1,1) = 
     -  0.3333333333333333D0*Cpl180 - 0.3333333333333333D0*Cpl181 + 
     -   0.3333333333333333D0*Cpl182 - 0.3333333333333333D0*Cpl183

#ifdef DEBUG
	print *, 'dZfR1dimreg(3,1,1) =', dZfR1dimreg(3,1,1)
#endif

        dZfR1dimreg(4,1,1) = 
     -  0.3333333333333333D0*Cpl180 - 0.3333333333333333D0*Cpl181 + 
     -   0.3333333333333333D0*Cpl182 - 0.3333333333333333D0*Cpl183

#ifdef DEBUG
	print *, 'dZfR1dimreg(4,1,1) =', dZfR1dimreg(4,1,1)
#endif

      
	dSW1 = -(0.5D0*(Cpl74*dMWsq1)) + 0.5D0*(Cpl75*dMZsq1)

#ifdef DEBUG
	print *, 'dSW1 =', dSW1
#endif

        dZAZ1 = dZZA1 - Cpl173*dZZA1 + 
     -   (-5*Cpl170 + Cpl171)*LInt41 - Cpl169*LInt43 - 
     -   0.3333333333333333D0*Cpl169 + 
     -   LInt39*(-Cpl172 - 2.5D0*Cpl169) + 
     -   LInt12*(2.5D0*Cpl170 - 0.5D0*Cpl171)

	do Gen1 = 1, 3
        dZAZ1 = dZAZ1 - 0.2222222222222222D0*Cpl174 + 
     -   0.5925925925925926D0*Cpl175 + 
     -   0.6666666666666666D0*Cpl176(Gen1) - 
     -   1.7777777777777777D0*Cpl177(Gen1) + 
     -   (-(0.6666666666666666D0*Cpl176(Gen1)) + 
     -      1.7777777777777777D0*Cpl177(Gen1))*LInt51(Gen1) + 
     -   (0.6666666666666666D0*Cpl176(Gen1) - 
     -      1.7777777777777777D0*Cpl177(Gen1))*LInt52(Gen1) + 
     -   (-(0.3333333333333333D0*Cpl174) + 
     -      1.3333333333333333D0*Cpl175)*LInt53(Gen1) + 
     -   (-(0.3333333333333333D0*Cpl174) + 
     -      0.4444444444444444D0*Cpl175)*LInt54(Gen1) + 
     -   (-(0.6666666666666666D0*Cpl174) + 
     -      1.7777777777777777D0*Cpl175)*LInt55(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dZAZ1 =', dZAZ1
#endif

	dZe1 = -(0.5D0*dZAA1) - 0.5D0*(Cpl178*dZZA1)

#ifdef DEBUG
	print *, 'dZe1 =', dZe1
#endif

        dZR1 = -(Cpl83*dMWsq1) - dZAA1 + Cpl85*dZZA1 - 
     -   0.3333333333333333D0*Cpl15 + 1.5D0*Cpl76 + 0.875D0*Cpl77 - 
     -   0.5D0*Cpl78 + LInt75*(-2*Cpl81 - 0.25D0*Cpl79) + 
     -   LInt12*(0.875D0*Cpl79 - 0.5D0*Cpl80 - 0.5D0*Cpl81) + 
     -   LInt13*(0.0625D0*Cpl79 + 0.75D0*Cpl81) + 
     -   LInt73*(0.25D0*Cpl44 - 0.5D0*Cpl82) - 2.D0*(Cpl84*dSW1) + 
     -   0.0625D0*(Cpl79*LInt11) - 
     -   0.08333333333333333D0*(Cpl15*LInt14) + 
     -   0.25D0*(Cpl14*LInt72) - 0.25D0*(Cpl79*LInt74) + 
     -   0.3333333333333333D0*(Cpl15*LInt76)

	do Gen1 = 1, 3
        dZR1 = dZR1 + 0.25D0*Cpl86(Gen1) - 
     -   0.25D0*(Cpl86(Gen1)*LInt77(Gen1)) + 
     -   0.25D0*(Cpl86(Gen1)*LInt78(Gen1))
	enddo

#ifdef DEBUG
	print *, 'dZR1 =', dZR1
#endif

        dZfL1photonic(1,1,1)=0D0

	dr=1
	If(dr.Eq.0) Then
	dZfL1photonic(2,1,1)=2*Cpl15*LInt363 + Cpl15*LInt364 + 
     -   0.25D0*Cpl15 + 0.5D0*(Cpl15*LInt361)

	dZfL1photonic(3,1,1)=0.1111111111111111D0*Cpl15 + 
     -   0.2222222222222222D0*(Cpl15*LInt366) +  
     -   0.8888888888888888D0*(Cpl15*LInt368) + 
     -   0.4444444444444444D0*(Cpl15*LInt369)

	dZfL1photonic(4,1,1)=0.027777777777777776D0*Cpl15 + 
     -   0.05555555555555555D0*(Cpl15*LInt370) + 
     -   0.2222222222222222D0*(Cpl15*LInt373) + 
     -   0.1111111111111111D0*(Cpl15*LInt374)

	Else
	
	dZfL1photonic(2,1,1)=
     -  Alfa*(0.25D0*Cpl180 - 0.25D0*Cpl181 + 
     -        0.25D0*Cpl182 - 0.25D0*Cpl183)

	dZfL1photonic(3,1,1)=
     -  4D0/9D0*Alfa*(0.25D0*Cpl180 - 0.25D0*Cpl181 + 
     -                0.25D0*Cpl182 - 0.25D0*Cpl183)

	dZfL1photonic(4,1,1)=
     -  1D0/9D0*Alfa*(0.25D0*Cpl180 - 0.25D0*Cpl181 + 
     -                0.25D0*Cpl182 - 0.25D0*Cpl183)
        Endif

        dZfR1photonic(1,1,1)=dZfL1photonic(1,1,1)	
	dZfR1photonic(2,1,1)=dZfL1photonic(2,1,1)     
        dZfR1photonic(3,1,1)=dZfL1photonic(3,1,1)
        dZfR1photonic(4,1,1)=dZfL1photonic(4,1,1)
     
	Cpl134 = Alfa*Pi*dZfR1photonic(3,1,1)
	Cpl135 = Alfa*Pi*(SW2c*dZfR1photonic(3,1,1))/CW2c
	Cpl136 = Alfa*Pi*dZfR1photonic(3,1,1)/CW2c
	Cpl137 = Alfa*Pi*dZfL1photonic(3,1,1)
	Cpl138 = Alfa*Pi*dZfL1photonic(3,1,1)/CW2c
	Cpl139 = Alfa*Pi*(SW2c*dZfL1photonic(3,1,1))/CW2c
	Cpl140 = Alfa*Pi*dZfL1photonic(3,1,1)/(CW2c*SW2c)
	Cpl141 = Alfa*Pi*dZfR1photonic(4,1,1)
	Cpl142 = Alfa*Pi*(SW2c*dZfR1photonic(4,1,1))/CW2c
	Cpl143 = Alfa*Pi*dZfR1photonic(4,1,1)/CW2c
	Cpl144 = Alfa*Pi*dZfL1photonic(4,1,1)
	Cpl145 = Alfa*Pi*dZfL1photonic(4,1,1)/CW2c
	Cpl146 = Alfa*Pi*(SW2c*dZfL1photonic(4,1,1))/CW2c
	Cpl147 = Alfa*Pi*dZfL1photonic(4,1,1)/(CW2c*SW2c)
	Cpl148 = Alfa*Pi*dZfR1photonic(2,1,1)
	Cpl149 = Alfa*Pi*(SW2c*dZfR1photonic(2,1,1))/CW2c
	Cpl150 = Alfa*Pi*dZfL1photonic(2,1,1)
	Cpl151 = Alfa*Pi*dZfL1photonic(2,1,1)/CW2c
	Cpl152 = Alfa*Pi*(SW2c*dZfL1photonic(2,1,1))/CW2c
	Cpl153 = Alfa*Pi*dZfR1photonic(2,1,1)/CW2c
	Cpl154 = Alfa*Pi*dZfL1photonic(2,1,1)/(CW2c*SW2c)
      
      end


