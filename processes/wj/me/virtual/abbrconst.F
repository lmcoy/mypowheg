      subroutine abbr_const()
      implicit character (a-f,h-s,u-z)
      implicit double complex (t)
      implicit integer (g)

#include "regcommon.h"
#include "params.h"

      integer resnr,delta
      double complex myA0,myB0,myB1,myB00,myB11,myB001,myB111,myDB0
      double complex myDB1,myDB00,myDB11,myBval,myCval,myDval,myEval
      double complex Prp,myDeltaIR,myDeltaUV
      integer dr

      integer myBget,myCget,myDget,myEget

#include "amps/looptools.h"
#include "amps/ampvars.h"
#include "amps/ampcommon.h"
#include "amps/ampcom.h"
	Cpl1 = (Alfa*Pi)/SW2c
	Cpl10 = Alfa/Pi/SW2c
	Cpl100 = Alfa2*CW2c/SW2c**2
	Cpl101 = Alfa2*MW2c/SW2c**2
	Cpl102 = Alfa2*(CW2c*MW2c)/SW2c**2
	Cpl103 = Alfa2*(CW2c*MZ2c)/SW2c**2
	Cpl105 = Alfa*Pi*MW2c/SW2c
	Cpl106 = Alfa*Pi*1/(SW2c*SWc)
	Cpl107 = SWc/CWc
	Cpl108 = Alfa/Pi/CW2c
	Cpl109 = Log(MT2)/Pi
	Cpl11 = Alfa/Pi
	Cpl110 = 1/Pi*Log(mureg)
	Cpl111 = myDeltaUV()
	Cpl112 = nflavour*myDeltaUV()
	Cpl113 = Log(lambda**2)
	Cpl114 = nflavour*Log(lambda**2)
	Cpl115 = Log(mureg)
	Cpl116 = nflavour*Log(mureg)
	Cpl117 = myDeltaIR()
	Cpl118 = nflavour*myDeltaIR()
	Cpl119 = 1/Pi*Log(lambda**2)
	Cpl12 = Alfa/Pi*CW2c/SW2c
	Cpl120 = 1/Pi*myDeltaIR()
	Cpl121 = 1/Pi*myDeltaUV()
	Cpl122 = (Alfa*Pi**2)/SW2c
	Cpl13 = Alfa/Pi*MW2r
	Cpl14 = Alfa/Pi*(CW2c*MW2c)/SW2c
	Cpl15 = Alfa/Pi*(CW2c*MW2r)/SW2c
	Cpl16 = Alfa/Pi*(CW2c*MWr**2)/SW2c
	Cpl17 = Alfa/Pi*(CW2c*MZ2c)/SW2c
	Cpl18 = Alfa/Pi*(MW2c*SW2c)/CW2c
	Cpl19 = Alfa/Pi*MW2c**2
	Cpl2 = Alfa*Pi*EL/SW2c
	Cpl20 = Alfa/Pi*(MW2c*MW2r)
	Cpl21 = Alfa/Pi*(MW2c*MWr**2)
	Cpl22 = Alfa/Pi*(MW2r*MWr**2)
	Cpl23 = Alfa/Pi*MW2c**2/SW2c
	Cpl24 = Alfa/Pi*(MW2c*MW2r)/SW2c
	Cpl25 = Alfa/Pi*(CW2c*MW2c*MWr**2)/SW2c
	Cpl26 = Alfa/Pi*(CW2c*MW2r*MWr**2)/SW2c
	Cpl27 = Alfa/Pi*(CW2c*MW2c*MZ2c)/SW2c
	Cpl28 = Alfa/Pi*(CW2c*MW2r*MZ2c)/SW2c
	Cpl29 = Alfa/Pi*(MW2c**2*SW2c)/CW2c
	Cpl3 = Alfa/SW2c
	Cpl30 = Alfa/Pi*(MW2c*MW2r*SW2c)/CW2c
	Cpl34 = Alfa/Pi*(MW2c*MWr**2)/SW2c
	Cpl35 = Alfa/Pi*(MW2r*MWr**2)/SW2c
	Cpl36 = Alfa/Pi*(CW2c*MZ2r)/SW2c
	Cpl37 = Alfa/Pi*MZr**2/(CW2c*SW2c)
	Cpl38 = Alfa/Pi*(CW2c*MZr**2)/SW2c
	Cpl39 = Alfa/Pi*1/(CW2c*SW2c)
	Cpl4 = Alfa*EL/SW2c
	Cpl40 = Alfa/Pi*SW2c/CW2c
	Cpl41 = Alfa/Pi*MW2c/(CW2c**2*SW2c)
	Cpl42 = Alfa/Pi*MZ2c/(CW2c*SW2c)
	Cpl43 = Alfa/Pi*MZ2r/(CW2c*SW2c)
	Cpl44 = Alfa/Pi*(MW2c*MZ2c)/(CW2c**2*SW2c)
	Cpl45 = Alfa/Pi*(MW2c*MZ2r)/(CW2c**2*SW2c)
	Cpl46 = Alfa/Pi*(CW2c*MW2c*MZ2r)/SW2c
	Cpl47 = Alfa/Pi*(CW2c*MZ2c*MZr**2)/SW2c
	Cpl48 = Alfa/Pi*(CW2c*MZ2r*MZr**2)/SW2c
	Cpl49 = Alfa/Pi*(MW2c*MZ2c*SW2c)/CW2c
	Cpl5 = Alfa/Pi*MW2c
	Cpl50 = Alfa/Pi*(MW2c*MZ2r*SW2c)/CW2c
	Cpl51 = Alfa/Pi*MZ2c
	Cpl52 = Alfa/Pi*MZ2r
	Cpl53 = Alfa/Pi*(MZ2c*SW2c)/CW2c
	Cpl54 = Alfa/Pi*(MZ2r*SW2c)/CW2c
	Cpl55 = Alfa/Pi*(MZ2c*MZr**2)/(CW2c*SW2c)
	Cpl56 = Alfa/Pi*(MZ2r*MZr**2)/(CW2c*SW2c)
	Cpl57 = Alfa/Pi*MZr**2/CW2c
	Cpl58 = Alfa/Pi*(MZr**2*SW2c)/CW2c
	Cpl6 = Alfa/Pi*MWr**2
	Cpl62 = Alfa/Pi*MZ2c/CW2c
	Cpl63 = Alfa/Pi*MZ2r/CW2c
	Cpl66 = Alfa/Pi*(MZ2c*MZr**2)/CW2c
	Cpl67 = Alfa/Pi*(MZ2r*MZr**2)/CW2c
	Cpl68 = Alfa/Pi*(MZ2c*MZr**2*SW2c)/CW2c
	Cpl69 = Alfa/Pi*(MZ2r*MZr**2*SW2c)/CW2c
	Cpl7 = Alfa/Pi*MW2c/SW2c
	Cpl70 = CW2c/(MW2c*SWc)
	Cpl71 = CW2c/(MZ2c*SWc)
	Cpl72 = Alfa*1/(pi*SW2c)
	Cpl73 = Alfa*Log(CW2c)/(pi*SW2c**2)
	Cpl74 = Alfa*Log(CW2c)/(pi*SW2c)
	Cpl75 = Alfa/Pi*1/(MW2c*SW2c)
	Cpl76 = Alfa/Pi/MW2c
	Cpl77 = Alfa/Pi*CW2c/(MW2c*SW2c)
	Cpl78 = Alfa/Pi*(CW2c*MZ2c)/(MW2c*SW2c)
	Cpl79 = 1/MW2c
	Cpl8 = Alfa/Pi*MW2r/SW2c
	Cpl80 = 1/SWc
	Cpl81 = CWc/SWc
	Cpl83 = Alfa/Pi*CW2c/(CWc*MZ2c*SWc)
	Cpl84 = Alfa/Pi*SW2c/(CWc*MZ2c*SWc)
	Cpl85 = Alfa/Pi*MW2c/(CWc*MZ2c*SWc)
	Cpl88 = Alfa2/SW2c
	Cpl89 = Alfa2*1/(CW2c*SW2c**2)
	Cpl9 = Alfa/Pi*MWr**2/SW2c
	Cpl90 = Alfa2*1/(CW2c*SW2c)
	Cpl91 = Alfa2/CW2c
	Cpl92 = Alfa2*MW2c/(CW2c*SW2c**2)
	Cpl93 = Alfa2*MW2c/(CW2c*SW2c)
	Cpl94 = Alfa2*MZ2c/CW2c
	Cpl95 = Alfa2*MZ2c/(CW2c*SW2c**2)
	Cpl96 = Alfa2*MZ2c/(CW2c*SW2c)
	Cpl97 = Alfa2*MW2c/SW2c
	Cpl98 = Alfa2*MW2c/CW2c
	Cpl99 = Alfa2/SW2c**2


	do Gen1=1,3
	Cpl104(Gen1) = Alfa2*MT2*delta(Gen1,3)/SW2c**2
	Cpl31(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/SW2c
	Cpl32(Gen1) = (Alfa*MT2)/Pi*(MW2c*delta(Gen1,3))/SW2c
	Cpl33(Gen1) = (Alfa*MT2)/Pi*(MW2r*delta(Gen1,3))/SW2c
	Cpl59(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/CW2c
	Cpl60(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/(CW2c*SW2c)
	Cpl61(Gen1) = (Alfa*MT2)/Pi*(SW2c*delta(Gen1,3))/CW2c
	Cpl64(Gen1) = (Alfa*MT2)/Pi*(MZ2c*delta(Gen1,3))/(CW2c*SW2c)
	Cpl65(Gen1) = (Alfa*MT2)/Pi*(MZ2r*delta(Gen1,3))/(CW2c*SW2c)
	Cpl82(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/(MW2c*SW2c)
	Cpl86(Gen1) = (Alfa*MT2)/Pi*delta(Gen1,3)/(CWc*MZ2c*SWc)
	Cpl87(Gen1) = (Alfa*MT2)/Pi*(SW2c*delta(Gen1,3))/(CWc*MZ2c*SWc)

        enddo



	IInt11 = myBget(0.D0,MW2r,MZ2r,0)
	LInt536 = myBval(bb0,IInt11)
	IInt12 = myBget(0.D0,MZ2r,MW2r,2)
	LInt538 = myBval(bb00,IInt12)


	do Gen1=1,3
	IInt16(Gen1) = myBget(0.D0,0.D0,MQU2(Gen1),0)
	LInt514(Gen1) = myBval(bb0,IInt16(Gen1))
	IInt23(Gen1) = myBget(0.D0,MQD2(Gen1),MQU2(Gen1),0)
	LInt540(Gen1) = myBval(bb0,IInt23(Gen1))
	IInt24(Gen1) = myBget(0.D0,MQU2(Gen1),MQD2(Gen1),1)
	LInt541(Gen1) = myBval(bb1,IInt24(Gen1))
        enddo

	IInt10 = myBget(0.D0,MH2,MW2r,2)
	LInt535 = myBval(bb0,IInt10)
	LInt537 = myBval(bb00,IInt10)
	IInt4 = myBget(MWr**2,MW2r,MZ2r,0)
	LInt480 = myBval(bb0,IInt4)
	LInt487 = myBval(dbb0,IInt4)


	do Gen1=1,3
	IInt13(Gen1) = myBget(MWr**2,MQD2(Gen1),MQU2(Gen1),0)
	LInt493(Gen1) = myBval(bb0,IInt13(Gen1))
	LInt496(Gen1) = myBval(dbb0,IInt13(Gen1))
	IInt20(Gen1) = myBget(0.D0,MLE2(Gen1),MLE2(Gen1),2)
	LInt529(Gen1) = myBval(bb1,IInt20(Gen1))
	LInt532(Gen1) = myBval(dbb00,IInt20(Gen1))
	IInt21(Gen1) = myBget(0.D0,MQD2(Gen1),MQD2(Gen1),2)
	LInt530(Gen1) = myBval(bb1,IInt21(Gen1))
	LInt533(Gen1) = myBval(dbb00,IInt21(Gen1))
        enddo

	IInt8 = myBget(MZr**2,0.D0,0.D0,2)
	LInt505 = myBval(bb1,IInt8)
	LInt509 = myBval(dbb00,IInt8)
	LInt512 = myBval(dbb1,IInt8)


	do Gen1=1,3
	IInt14(Gen1) = myBget(MWr**2,0.D0,MLE2(Gen1),2)
	LInt494(Gen1) = myBval(bb1,IInt14(Gen1))
	LInt497(Gen1) = myBval(dbb00,IInt14(Gen1))
	LInt499(Gen1) = myBval(dbb1,IInt14(Gen1))
	IInt15(Gen1) = myBget(MWr**2,MQU2(Gen1),MQD2(Gen1),2)
	LInt495(Gen1) = myBval(bb1,IInt15(Gen1))
	LInt498(Gen1) = myBval(dbb00,IInt15(Gen1))
	LInt500(Gen1) = myBval(dbb1,IInt15(Gen1))
	IInt18(Gen1) = myBget(MZr**2,MLE2(Gen1),MLE2(Gen1),2)
	LInt516(Gen1) = myBval(bb1,IInt18(Gen1))
	LInt520(Gen1) = myBval(dbb00,IInt18(Gen1))
	LInt523(Gen1) = myBval(dbb1,IInt18(Gen1))
	IInt19(Gen1) = myBget(MZr**2,MQD2(Gen1),MQD2(Gen1),2)
	LInt517(Gen1) = myBval(bb1,IInt19(Gen1))
	LInt521(Gen1) = myBval(dbb00,IInt19(Gen1))
	LInt524(Gen1) = myBval(dbb1,IInt19(Gen1))
	IInt22(Gen1) = myBget(0.D0,MQU2(Gen1),MQU2(Gen1),2)
	LInt531(Gen1) = myBval(bb1,IInt22(Gen1))
	LInt534(Gen1) = myBval(dbb00,IInt22(Gen1))
	LInt543(Gen1) = myBval(bb0,IInt22(Gen1))
        enddo

	IInt3 = myBget(MWr**2,MH2,MW2r,2)
	LInt479 = myBval(bb0,IInt3)
	LInt481 = myBval(bb00,IInt3)
	LInt486 = myBval(dbb0,IInt3)
	LInt489 = myBval(dbb00,IInt3)
	IInt5 = myBget(MWr**2,MZ2r,MW2r,2)
	LInt482 = myBval(bb00,IInt5)
	LInt484 = myBval(bb1,IInt5)
	LInt490 = myBval(dbb00,IInt5)
	LInt492 = myBval(dbb1,IInt5)
	IInt6 = myBget(MZr**2,MH2,MZ2r,2)
	LInt501 = myBval(bb0,IInt6)
	LInt503 = myBval(bb00,IInt6)
	LInt507 = myBval(dbb0,IInt6)
	LInt510 = myBval(dbb00,IInt6)
	IInt9 = myBget(0.D0,MW2r,MW2r,2)
	LInt526 = myBval(bb0,IInt9)
	LInt527 = myBval(bb1,IInt9)
	LInt528 = myBval(dbb00,IInt9)
	LInt542 = myBval(bb00,IInt9)
	IInt1 = myBget(0.D0,0.D0,MW2r,1)
	LInt2576 = myBval(bb0,IInt1)
	LInt2580 = myBval(bb0,IInt1)
	LInt2584 = myBval(bb0,IInt1)
	LInt477 = myBval(bb0,IInt1)
	LInt539 = myBval(bb1,IInt1)
	IInt2 = myBget(MWr**2,0.D0,MW2r,2)
	LInt478 = myBval(bb0,IInt2)
	LInt483 = myBval(bb1,IInt2)
	LInt485 = myBval(dbb0,IInt2)
	LInt488 = myBval(dbb00,IInt2)
	LInt491 = myBval(dbb1,IInt2)


	do Gen1=1,3
	IInt17(Gen1) = myBget(MZr**2,MQU2(Gen1),MQU2(Gen1),2)
	LInt515(Gen1) = myBval(bb0,IInt17(Gen1))
	LInt518(Gen1) = myBval(bb1,IInt17(Gen1))
	LInt519(Gen1) = myBval(dbb0,IInt17(Gen1))
	LInt522(Gen1) = myBval(dbb00,IInt17(Gen1))
	LInt525(Gen1) = myBval(dbb1,IInt17(Gen1))
        enddo

	IInt7 = myBget(MZr**2,MW2r,MW2r,2)
	LInt502 = myBval(bb0,IInt7)
	LInt504 = myBval(bb00,IInt7)
	LInt506 = myBval(bb1,IInt7)
	LInt508 = myBval(dbb0,IInt7)
	LInt511 = myBval(dbb00,IInt7)
	LInt513 = myBval(dbb1,IInt7)
	LInt2573 = myA0(MH2)
	LInt2574 = myA0(MW2r)
	LInt2575 = myA0(MZ2r)
	LInt2577 = myA0(MH2)
	LInt2578 = myA0(MW2r)
	LInt2579 = myA0(MZ2r)
	LInt2581 = myA0(MH2)
	LInt2582 = myA0(MW2r)
	LInt2583 = myA0(MZ2r)
	LInt2585 = myB1(0.D0,0.D0,MZ2r)
	LInt2586 = myB1(0.D0,ME2,MW2r)
	LInt2587 = myB1(ME2,0.D0,MW2r)
	LInt2588 = myB1(ME2,ME2,0.D0)
	LInt2589 = myB1(ME2,ME2,MZ2r)
	LInt2590 = ME2*myDB0(ME2,0.D0,ME2)
	LInt2591 = ME2*myDB1(ME2,ME2,0.D0)
	LInt2592 = myB1(MU2,MD2,MW2r)
	LInt2593 = myB1(MU2,MU2,0.D0)
	LInt2594 = myB1(MU2,MU2,MZ2r)
	LInt2595 = MU2*myDB0(MU2,0.D0,MU2)
	LInt2596 = MU2*myDB1(MU2,MU2,0.D0)
	LInt2597 = myB1(MD2,MD2,0.D0)
	LInt2598 = myB1(MD2,MD2,MZ2r)
	LInt2599 = myB1(MD2,MU2,MW2r)
	LInt2600 = MD2*myDB0(MD2,0.D0,MD2)
	LInt2601 = MD2*myDB1(MD2,MD2,0.D0)
	LInt474 = myA0(MH2)
	LInt475 = myA0(MW2r)
	LInt476 = myA0(MZ2r)

        dMWsq1 = (2*Cpl13 - 2*Cpl5)*LInt488 +
     -   LInt482*(-2*Cpl12 - 0.25D0*Cpl10) +
     -   LInt475*(0.875D0*Cpl10 - 0.5D0*Cpl11 - 0.5D0*Cpl12) +
     -   LInt476*(0.0625D0*Cpl10 + 0.75D0*Cpl12) +
     -   LInt484*(-(0.5D0*Cpl14) + 0.5D0*Cpl15 - 0.5D0*Cpl16) +
     -   LInt480*(-(1.25D0*Cpl14) + 1.25D0*Cpl15 - 1.25D0*Cpl16 -
     -      0.5D0*Cpl17 + 0.25D0*Cpl18) +
     -   LInt491*(-(0.5D0*Cpl21) + 0.5D0*Cpl22) +
     -   LInt485*(0.25D0*Cpl19 - 0.25D0*Cpl20 - 1.25D0*Cpl21 +
     -      1.25D0*Cpl22) + LInt486*(0.25D0*Cpl23 - 0.25D0*Cpl24) +
     -   LInt492*(-(0.5D0*Cpl25) + 0.5D0*Cpl26) +
     -   LInt487*(-(1.25D0*Cpl25) + 1.25D0*Cpl26 - 0.5D0*Cpl27 +
     -      0.5D0*Cpl28 + 0.25D0*Cpl29 - 0.25D0*Cpl30) -
     -   0.3333333333333333D0*Cpl5 + 0.1111111111111111D0*Cpl6 +
     -   LInt483*(0.5D0*Cpl13 - 0.16666666666666666D0*Cpl5 -
     -      0.8333333333333334D0*Cpl6) +
     -   LInt478*(-Cpl5 + 1.25D0*Cpl13 - 1.25D0*Cpl6) -
     -   0.16666666666666666D0*Cpl7 + 0.16666666666666666D0*Cpl8 +
     -   LInt489*(-(0.25D0*Cpl7) + 0.25D0*Cpl8) +
     -   LInt490*(-2*Cpl14 + 2*Cpl15 - 0.25D0*Cpl7 + 0.25D0*Cpl8) -
     -   0.16666666666666666D0*Cpl9 + 0.0625D0*(Cpl10*LInt474) -
     -   0.3333333333333333D0*(Cpl5*LInt477) + 0.25D0*(Cpl7*LInt479) -
     -   0.25D0*(Cpl10*LInt481)

	do Gen1 = 1, 3
        dMWsq1 = dMWsq1 - 0.1111111111111111D0*Cpl9 +
     -   0.25D0*Cpl31(Gen1) - 0.25D0*(Cpl31(Gen1)*LInt493(Gen1)) +
     -   (-(0.25D0*Cpl7) + 0.25D0*Cpl8 - 0.16666666666666666D0*Cpl9)*
     -    LInt494(Gen1) + (-(0.75D0*Cpl7) + 0.75D0*Cpl8 - 0.5D0*Cpl9 +
     -      0.25D0*Cpl31(Gen1))*LInt495(Gen1) +
     -   (-(0.75D0*Cpl32(Gen1)) + 0.75D0*Cpl33(Gen1))*LInt496(Gen1) +
     -   (0.5D0*Cpl7 - 0.5D0*Cpl8)*LInt497(Gen1) +
     -   (1.5D0*Cpl7 - 1.5D0*Cpl8)*LInt498(Gen1) +
     -   (-(0.25D0*Cpl34) + 0.25D0*Cpl35)*LInt499(Gen1) +
     -   (-(0.75D0*Cpl34) + 0.75D0*Cpl35)*LInt500(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dMWsq1 =', dMWsq1
#endif

        dMZsq1 = -(0.16666666666666666D0*Cpl17) +
     -   0.16666666666666666D0*Cpl36 - 0.041666666666666664D0*Cpl37 -
     -   0.16666666666666666D0*Cpl38 +
     -   LInt506*(-(0.5D0*Cpl17) + 0.5D0*Cpl36 - 0.5D0*Cpl38) +
     -   LInt502*(-(0.5D0*Cpl14) - 1.25D0*Cpl17 + 0.5D0*Cpl18 +
     -      1.25D0*Cpl36 - 1.25D0*Cpl38) +
     -   LInt475*(-(0.25D0*Cpl11) + 1.125D0*Cpl12 + 0.125D0*Cpl40) +
     -   LInt504*(0.5D0*Cpl11 - 2.25D0*Cpl12 - 0.25D0*Cpl40) +
     -   LInt510*(-(0.25D0*Cpl42) + 0.25D0*Cpl43) +
     -   LInt505*(-(0.25D0*Cpl37) - 0.375D0*Cpl42 + 0.375D0*Cpl43) +
     -   LInt509*(0.75D0*Cpl42 - 0.75D0*Cpl43) +
     -   LInt507*(0.25D0*Cpl44 - 0.25D0*Cpl45) +
     -   LInt513*(-(0.5D0*Cpl47) + 0.5D0*Cpl48) +
     -   LInt508*(-(0.5D0*Cpl27) + 0.5D0*Cpl46 - 1.25D0*Cpl47 +
     -      1.25D0*Cpl48 + 0.5D0*Cpl49 - 0.5D0*Cpl50) +
     -   LInt511*(-(2.25D0*Cpl17) + 2.25D0*Cpl36 + 0.5D0*Cpl51 -
     -      0.5D0*Cpl52 - 0.25D0*Cpl53 + 0.25D0*Cpl54) +
     -   LInt512*(-(0.375D0*Cpl55) + 0.375D0*Cpl56) +
     -   0.0625D0*(Cpl39*LInt474) + 0.0625D0*(Cpl39*LInt476) +
     -   0.25D0*(Cpl41*LInt501) - 0.25D0*(Cpl39*LInt503)

	do Gen1 = 1, 3
        dMZsq1 = dMZsq1 - 0.09722222222222222D0*Cpl37 +
     -   0.2222222222222222D0*Cpl57 - 0.2962962962962963D0*Cpl58 -
     -   0.6666666666666666D0*Cpl59(Gen1) + 0.25D0*Cpl60(Gen1) +
     -   0.8888888888888888D0*Cpl61(Gen1) +
     -   (0.6666666666666666D0*Cpl59(Gen1) - 0.25D0*Cpl60(Gen1) -
     -      0.8888888888888888D0*Cpl61(Gen1))*LInt514(Gen1) +
     -   (-(0.6666666666666666D0*Cpl59(Gen1)) - 0.125D0*Cpl60(Gen1) +
     -      0.8888888888888888D0*Cpl61(Gen1))*LInt515(Gen1) +
     -   (-Cpl53 + Cpl54 - 0.08333333333333333D0*Cpl37 -
     -      0.125D0*Cpl42 + 0.125D0*Cpl43 + 0.3333333333333333D0*Cpl57 -
     -      0.6666666666666666D0*Cpl58 + 0.5D0*Cpl62 - 0.5D0*Cpl63)*
     -    LInt516(Gen1) + (-(0.25D0*Cpl37) - 0.375D0*Cpl42 +
     -      0.375D0*Cpl43 - 0.3333333333333333D0*Cpl53 +
     -      0.3333333333333333D0*Cpl54 + 0.3333333333333333D0*Cpl57 -
     -      0.2222222222222222D0*Cpl58 + 0.5D0*Cpl62 - 0.5D0*Cpl63)*
     -    LInt517(Gen1) + (Cpl62 - Cpl63 - 0.25D0*Cpl37 -
     -      0.375D0*Cpl42 + 0.375D0*Cpl43 - 1.3333333333333333D0*Cpl53 +
     -      1.3333333333333333D0*Cpl54 + 0.6666666666666666D0*Cpl57 -
     -      0.8888888888888888D0*Cpl58)*LInt518(Gen1) +
     -   (-(0.375D0*Cpl64(Gen1)) + 0.375D0*Cpl65(Gen1))*
     -    LInt519(Gen1) + (2*Cpl53 - 2*Cpl54 - Cpl62 + Cpl63 +
     -      0.25D0*Cpl42 - 0.25D0*Cpl43)*LInt520(Gen1) +
     -   (-Cpl62 + Cpl63 + 0.75D0*Cpl42 - 0.75D0*Cpl43 +
     -      0.6666666666666666D0*Cpl53 - 0.6666666666666666D0*Cpl54)*
     -    LInt521(Gen1) + (-2*Cpl62 + 2*Cpl63 + 0.75D0*Cpl42 -
     -      0.75D0*Cpl43 + 2.6666666666666665D0*Cpl53 -
     -      2.6666666666666665D0*Cpl54)*LInt522(Gen1) +
     -   (-Cpl68 + Cpl69 - 0.125D0*Cpl55 + 0.125D0*Cpl56 + 0.5D0*Cpl66 -
     -      0.5D0*Cpl67)*LInt523(Gen1) +
     -   (-(0.375D0*Cpl55) + 0.375D0*Cpl56 + 0.5D0*Cpl66 - 0.5D0*Cpl67 -
     -      0.3333333333333333D0*Cpl68 + 0.3333333333333333D0*Cpl69)*
     -    LInt524(Gen1) + (Cpl66 - Cpl67 - 0.375D0*Cpl55 +
     -      0.375D0*Cpl56 - 1.3333333333333333D0*Cpl68 +
     -      1.3333333333333333D0*Cpl69)*LInt525(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dMZsq1 =', dMZsq1
#endif

        dZAA1 = 3*Cpl11*LInt528 + 0.16666666666666666D0*Cpl11 +
     -   1.25D0*(Cpl11*LInt526) + 0.5D0*(Cpl11*LInt527)

	do Gen1 = 1, 3
        dZAA1 = dZAA1 + 0.3333333333333333D0*
     -    (Cpl11*LInt530(Gen1)) +
     -   1.3333333333333333D0*(Cpl11*LInt531(Gen1)) -
     -   0.6666666666666666D0*(Cpl11*LInt533(Gen1)) -
     -   2.6666666666666665D0*(Cpl11*LInt534(Gen1)) +
     -   Cpl11*LInt529(Gen1) - 2*Cpl11*LInt532(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dZAA1 =', dZAA1
#endif

        dZg1dimreg = -(0.08333333333333333D0*Cpl109) +
     -   0.08333333333333333D0*Cpl110 - 0.411150269654063D0*Cpl111 +
     -   0.026525823848649224D0*Cpl112

#ifdef DEBUG
	print *, 'dZg1dimreg =', dZg1dimreg
#endif

        dZGG1dimreg = 0.16666666666666666D0*Cpl109 +
     -   0.3448357100324399D0*Cpl111 - 0.05305164769729845D0*Cpl112 -
     -   0.3978873577297384D0*Cpl113 + 0.05305164769729845D0*Cpl114 +
     -   0.3448357100324399D0*Cpl115 - 0.05305164769729845D0*Cpl116 -
     -   0.3978873577297384D0*Cpl117 + 0.05305164769729845D0*Cpl118

#ifdef DEBUG
	print *, 'dZGG1dimreg =', dZGG1dimreg
#endif

        dZW1 = 2*Cpl11*LInt488 + 0.16666666666666666D0*Cpl10 +
     -   LInt490*(2*Cpl12 + 0.25D0*Cpl10) +
     -   LInt487*(1.25D0*Cpl16 + 0.5D0*Cpl17 - 0.25D0*Cpl18) +
     -   LInt485*(-(0.25D0*Cpl5) + 1.25D0*Cpl6) +
     -   1.25D0*(Cpl11*LInt478) + 1.25D0*(Cpl12*LInt480) +
     -   0.5D0*(Cpl11*LInt483) + 0.5D0*(Cpl12*LInt484) -
     -   0.25D0*(Cpl7*LInt486) + 0.25D0*(Cpl10*LInt489) +
     -   0.5D0*(Cpl6*LInt491) + 0.5D0*(Cpl16*LInt492)

	do Gen1 = 1, 3
        dZW1 = dZW1 + 0.25D0*(Cpl10*LInt494(Gen1)) +
     -   0.75D0*(Cpl10*LInt495(Gen1)) +
     -   0.75D0*(Cpl31(Gen1)*LInt496(Gen1)) -
     -   0.5D0*(Cpl10*LInt497(Gen1)) - 1.5D0*(Cpl10*LInt498(Gen1)) +
     -   0.25D0*(Cpl9*LInt499(Gen1)) + 0.75D0*(Cpl9*LInt500(Gen1))
	enddo

#ifdef DEBUG
	print *, 'dZW1 =', dZW1
#endif

        dZZA1 = Cpl85*LInt526 + (5*Cpl83 - Cpl84)*LInt542 +
     -   LInt475*(-(2.5D0*Cpl83) + 0.5D0*Cpl84)

	do Gen1 = 1, 3
        dZZA1 = dZZA1 - 0.6666666666666666D0*Cpl86(Gen1) +
     -   1.7777777777777777D0*Cpl87(Gen1) +
     -   (0.6666666666666666D0*Cpl86(Gen1) -
     -      1.7777777777777777D0*Cpl87(Gen1))*LInt514(Gen1) +
     -   (-(0.6666666666666666D0*Cpl86(Gen1)) +
     -      1.7777777777777777D0*Cpl87(Gen1))*LInt543(Gen1)
	enddo

#ifdef DEBUG
	print *, 'dZZA1 =', dZZA1
#endif

        dZfL1(1,1,1) = 0.125D0*Cpl10 + 0.0625D0*Cpl39 +
     -   0.125D0*(Cpl39*LInt2585) + 0.25D0*(Cpl10*LInt2586)

#ifdef DEBUG
	print *, 'dZfL1(1,1,1) =', dZfL1(1,1,1)
#endif

        dZfL1(2,1,1) = 2*Cpl11*LInt2590 + Cpl11*LInt2591 +
     -   0.125D0*Cpl10 - 0.25D0*Cpl108 + 0.25D0*Cpl11 + 0.0625D0*Cpl39 +
     -   0.25D0*Cpl40 + LInt2589*
     -    (-(0.5D0*Cpl108) + 0.125D0*Cpl39 + 0.5D0*Cpl40) +
     -   0.25D0*(Cpl10*LInt2587) + 0.5D0*(Cpl11*LInt2588)

	If(dr.Eq.1) Then
	 dZfL1(2,1,1) = dZfL1(2,1,1) -
     -  (2*Cpl11*LInt2590 + Cpl11*LInt2591 + 0.25D0*Cpl11 +
     -   0.5D0*(Cpl11*LInt2588)) +
     -   Alfa*(0.25D0*Cpl119 - 0.25D0*Cpl110 +
     -        0.25D0*Cpl120 - 0.25D0*Cpl121)
	Endif

#ifdef DEBUG
	print *, 'dZfL1(2,1,1) =', dZfL1(2,1,1)
#endif

        dZfL1(3,1,1) = 0.125D0*Cpl10 -
     -   0.16666666666666666D0*Cpl108 + 0.1111111111111111D0*Cpl11 +
     -   0.0625D0*Cpl39 + 0.1111111111111111D0*Cpl40 +
     -   LInt2594*(-(0.3333333333333333D0*Cpl108) + 0.125D0*Cpl39 +
     -      0.2222222222222222D0*Cpl40) + 0.25D0*(Cpl10*LInt2592) +
     -   0.2222222222222222D0*(Cpl11*LInt2593) +
     -   0.8888888888888888D0*(Cpl11*LInt2595) +
     -   0.4444444444444444D0*(Cpl11*LInt2596)
       If(dr.Eq.1) Then
        dZfL1(3,1,1) = dZfL1(3,1,1) -
     - (0.1111111111111111D0*Cpl11+
     -   0.2222222222222222D0*(Cpl11*LInt2593) +
     -   0.8888888888888888D0*(Cpl11*LInt2595) +
     -   0.4444444444444444D0*(Cpl11*LInt2596)) +
     -   4D0/9D0*Alfa*(0.25D0*Cpl119 - 0.25D0*Cpl110 +
     -        0.25D0*Cpl120 - 0.25D0*Cpl121)
       Endif

#ifdef DEBUG
	print *, 'dZfL1(3,1,1) =', dZfL1(3,1,1)
#endif

        dZfL1(4,1,1) = 0.125D0*Cpl10 -
     -   0.08333333333333333D0*Cpl108 + 0.027777777777777776D0*Cpl11 +
     -   0.0625D0*Cpl39 + 0.027777777777777776D0*Cpl40 +
     -   LInt2598*(-(0.16666666666666666D0*Cpl108) + 0.125D0*Cpl39 +
     -      0.05555555555555555D0*Cpl40) +
     -   0.05555555555555555D0*(Cpl11*LInt2597) +
     -   0.25D0*(Cpl10*LInt2599) +
     -   0.2222222222222222D0*(Cpl11*LInt2600) +
     -   0.1111111111111111D0*(Cpl11*LInt2601)
       If(dr.Eq.1) Then
        dZfL1(4,1,1) = dZfL1(4,1,1) -
     - (0.027777777777777776D0*Cpl11 +
     -   0.05555555555555555D0*(Cpl11*LInt2597) +
     -   0.2222222222222222D0*(Cpl11*LInt2600) +
     -   0.1111111111111111D0*(Cpl11*LInt2601)) +
     -   1D0/9D0*Alfa*(0.25D0*Cpl119 - 0.25D0*Cpl110 +
     -        0.25D0*Cpl120 - 0.25D0*Cpl121)
       Endif

#ifdef DEBUG
	print *, 'dZfL1(4,1,1) =', dZfL1(4,1,1)
#endif

        dZfL1dimreg(3,1,1) =
     -  -(0.3333333333333333D0*Cpl110) + 0.3333333333333333D0*Cpl119 +
     -   0.3333333333333333D0*Cpl120 - 0.3333333333333333D0*Cpl121

#ifdef DEBUG
	print *, 'dZfL1dimreg(3,1,1) =', dZfL1dimreg(3,1,1)
#endif

        dZfL1dimreg(4,1,1) =
     -  -(0.3333333333333333D0*Cpl110) + 0.3333333333333333D0*Cpl119 +
     -   0.3333333333333333D0*Cpl120 - 0.3333333333333333D0*Cpl121

#ifdef DEBUG
	print *, 'dZfL1dimreg(4,1,1) =', dZfL1dimreg(4,1,1)
#endif


	dSW1 = -(0.5D0*(Cpl70*dMWsq1)) + 0.5D0*(Cpl71*dMZsq1)

#ifdef DEBUG
	print *, 'dSW1 =', dSW1
#endif

	dZe1 = -(0.5D0*dZAA1) - 0.5D0*(Cpl107*dZZA1)

#ifdef DEBUG
	print *, 'dZe1 =', dZe1
#endif

        dZR1 = -(Cpl79*dMWsq1) - dZAA1 + Cpl81*dZZA1 -
     -   0.3333333333333333D0*Cpl11 + 1.5D0*Cpl72 + 0.875D0*Cpl73 -
     -   0.5D0*Cpl74 + LInt538*(-2*Cpl77 - 0.25D0*Cpl75) +
     -   LInt475*(0.875D0*Cpl75 - 0.5D0*Cpl76 - 0.5D0*Cpl77) +
     -   LInt476*(0.0625D0*Cpl75 + 0.75D0*Cpl77) +
     -   LInt536*(0.25D0*Cpl40 - 0.5D0*Cpl78) - 2.D0*(Cpl80*dSW1) +
     -   0.0625D0*(Cpl75*LInt474) -
     -   0.08333333333333333D0*(Cpl11*LInt477) +
     -   0.25D0*(Cpl10*LInt535) - 0.25D0*(Cpl75*LInt537) +
     -   0.3333333333333333D0*(Cpl11*LInt539)

	do Gen1 = 1, 3
        dZR1 = dZR1 + 0.25D0*Cpl82(Gen1) -
     -   0.25D0*(Cpl82(Gen1)*LInt540(Gen1)) +
     -   0.25D0*(Cpl82(Gen1)*LInt541(Gen1))
	enddo

#ifdef DEBUG
	print *, 'dZR1 =', dZR1
#endif


      end
