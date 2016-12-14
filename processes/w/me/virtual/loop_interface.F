#define COLILib
C#define DDLib 

#define rmax1 16
#define rmax2 16
#define rmax3 16
#define rmax4 16
#define rmax5 5

**************************************************************************************

      double complex function myDeltaIR()
      implicit none

#include "params.h"
#include "amps/ampcommon.h"

      myDeltaIR=epsfac

      end

**************************************************************************************

      double complex function myDeltaUV()
      implicit none

#include "params.h"
#include "amps/ampcommon.h"

      myDeltaUV=deltareg

      end

**************************************************************************************

      double complex function com_mass(m1)
      implicit none

#include "params.h"
#include "amps/ampcommon.h"

      double precision m1
      double complex mass1
      
      If(m1.EQ.MW2r) Then
        com_mass=MW2c
      Else If(m1.EQ.MZ2r) Then
        com_mass=MZ2c
      Else
        com_mass=m1
      Endif

#ifdef LTLib
      If(dimag(com_mass).NE.0D0) Then
        Print *,"Looptools not compatible with complex masses"
	Stop
      Endif
#endif
      
      end

**************************************************************************************

#ifdef DDLib
      Subroutine small_mass1(mass1,mm1)
      Implicit None
      
      double complex mass1,mm1
      double precision smallm2(10)
      integer ii,iimax

      common/smallm2_common/smallm2,iimax
      
      do ii=1,iimax
        If(dreal(mass1).EQ.smallm2(ii)) Then
	   mm1=ii*1D-20
	   Return
	Endif
      Enddo
      
      mm1=mass1
      
      End

**************************************************************************************

      Subroutine small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      Implicit None
      
      double complex mass1,mm1,mass2,mm2
      double precision smallm2(10),p1,pp1
      integer ii,iimax

      common/smallm2_common/smallm2,iimax
      
      do ii=1,iimax
        If(dreal(mass1).EQ.smallm2(ii)) Then
	   mm1=ii*1D-20
	   goto 99
	Endif
      Enddo
      mm1=mass1
 99   Continue      
      do ii=1,iimax
        If(dreal(mass2).EQ.smallm2(ii)) Then
	   mm2=ii*1D-20
	   goto 98
	Endif
      Enddo
      mm2=mass2
 98   Continue      
      do ii=1,iimax
        If(p1.EQ.smallm2(ii)) Then
	   pp1=ii*1D-20
	   goto 89
	Endif
      Enddo
      pp1=p1
 89   Continue      
      
      End

**************************************************************************************

      Subroutine small_mass3(p1,p2,p3,mass1,mass2,mass3,pp1,pp2,pp3,mm1,mm2,mm3)
      Implicit None
      
      double complex mass1,mm1,mass2,mm2,mass3,mm3
      double precision smallm2(10),p1,pp1,p2,pp2,p3,pp3
      integer ii,iimax

      common/smallm2_common/smallm2,iimax
      
      do ii=1,iimax
        If(dreal(mass1).EQ.smallm2(ii)) Then
	   mm1=ii*1D-20
	   goto 99
	Endif
      Enddo
      mm1=mass1
 99   Continue      
      do ii=1,iimax
        If(dreal(mass2).EQ.smallm2(ii)) Then
	   mm2=ii*1D-20
	   goto 98
	Endif
      Enddo
      mm2=mass2
 98   Continue      
      do ii=1,iimax
        If(dreal(mass3).EQ.smallm2(ii)) Then
	   mm3=ii*1D-20
	   goto 97
	Endif
      Enddo
      mm3=mass3
 97   Continue      
      do ii=1,iimax
        If(p1.EQ.smallm2(ii)) Then
	   pp1=ii*1D-20
	   goto 89
	Endif
      Enddo
      pp1=p1
 89   Continue      
      do ii=1,iimax
        If(p2.EQ.smallm2(ii)) Then
	   pp2=ii*1D-20
	   goto 88
	Endif
      Enddo
      pp2=p2
 88   Continue      
      do ii=1,iimax
        If(p3.EQ.smallm2(ii)) Then
	   pp3=ii*1D-20
	   goto 87
	Endif
      Enddo
      pp3=p3
 87   Continue      
      
      End

**************************************************************************************

      Subroutine small_mass4(p1,p2,p3,p4,p5,p6,mass1,mass2,mass3,mass4,
     &                       pp1,pp2,pp3,pp4,pp5,pp6,mm1,mm2,mm3,mm4)
      Implicit None
      
      double complex mass1,mm1,mass2,mm2,mass3,mm3,mass4,mm4
      double precision smallm2(10),p1,pp1,p2,pp2,p3,pp3,p4,pp4,p5,pp5,p6,pp6
      integer ii,iimax

      common/smallm2_common/smallm2,iimax
      
      do ii=1,iimax
        If(dreal(mass1).EQ.smallm2(ii)) Then
	   mm1=ii*1D-20
	   goto 99
	Endif
      Enddo
      mm1=mass1
 99   Continue      
      do ii=1,iimax
        If(dreal(mass2).EQ.smallm2(ii)) Then
	   mm2=ii*1D-20
	   goto 98
	Endif
      Enddo
      mm2=mass2
 98   Continue      
      do ii=1,iimax
        If(dreal(mass3).EQ.smallm2(ii)) Then
	   mm3=ii*1D-20
	   goto 97
	Endif
      Enddo
      mm3=mass3
 97   Continue      
      do ii=1,iimax
        If(dreal(mass4).EQ.smallm2(ii)) Then
	   mm4=ii*1D-20
	   goto 96
	Endif
      Enddo
      mm4=mass4
 96   Continue      
      do ii=1,iimax
        If(p1.EQ.smallm2(ii)) Then
	   pp1=ii*1D-20
	   goto 89
	Endif
      Enddo
      pp1=p1
 89   Continue      
      do ii=1,iimax
        If(p2.EQ.smallm2(ii)) Then
	   pp2=ii*1D-20
	   goto 88
	Endif
      Enddo
      pp2=p2
 88   Continue      
      do ii=1,iimax
        If(p3.EQ.smallm2(ii)) Then
	   pp3=ii*1D-20
	   goto 87
	Endif
      Enddo
      pp3=p3
 87   Continue      
      do ii=1,iimax
        If(p4.EQ.smallm2(ii)) Then
	   pp4=ii*1D-20
	   goto 86
	Endif
      Enddo
      pp4=p4
 86   Continue      
      do ii=1,iimax
        If(p5.EQ.smallm2(ii)) Then
	   pp5=ii*1D-20
	   goto 85
	Endif
      Enddo
      pp5=p5
 85   Continue      
      do ii=1,iimax
        If(p6.EQ.smallm2(ii)) Then
	   pp6=ii*1D-20
	   goto 84
	Endif
      Enddo
      pp6=p6
 84   Continue      
      
      End

**************************************************************************************

      Subroutine small_mass5(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,mass1,mass2,mass3,mass4,mass5,
     &                       pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10,mm1,mm2,mm3,mm4,mm5)
      Implicit None
      
      double complex mass1,mm1,mass2,mm2,mass3,mm3,mass4,mm4,mass5,mm5
      double precision smallm2(10),p1,pp1,p2,pp2,p3,pp3,p4,pp4,p5,pp5,p6,pp6
      double precision p7,pp7,p8,pp8,p9,pp9,p10,pp10
      integer ii,iimax

      common/smallm2_common/smallm2,iimax
      
      do ii=1,iimax
        If(dreal(mass1).EQ.smallm2(ii)) Then
	   mm1=ii*1D-20
	   goto 99
	Endif
      Enddo
      mm1=mass1
 99   Continue      
      do ii=1,iimax
        If(dreal(mass2).EQ.smallm2(ii)) Then
	   mm2=ii*1D-20
	   goto 98
	Endif
      Enddo
      mm2=mass2
 98   Continue      
      do ii=1,iimax
        If(dreal(mass3).EQ.smallm2(ii)) Then
	   mm3=ii*1D-20
	   goto 97
	Endif
      Enddo
      mm3=mass3
 97   Continue      
      do ii=1,iimax
        If(dreal(mass4).EQ.smallm2(ii)) Then
	   mm4=ii*1D-20
	   goto 96
	Endif
      Enddo
      mm4=mass4
 96   Continue      
      do ii=1,iimax
        If(dreal(mass5).EQ.smallm2(ii)) Then
	   mm5=ii*1D-20
	   goto 95
	Endif
      Enddo
      mm5=mass5
 95   Continue      
      do ii=1,iimax
	If(p1.EQ.smallm2(ii)) Then
	   pp1=ii*1D-20
	   goto 89
	Endif
      Enddo
      pp1=p1
 89   Continue      
      do ii=1,iimax
        If(p2.EQ.smallm2(ii)) Then
	   pp2=ii*1D-20
	   goto 88
	Endif
      Enddo
      pp2=p2
 88   Continue      
      do ii=1,iimax
        If(p3.EQ.smallm2(ii)) Then
	   pp3=ii*1D-20
	   goto 87
	Endif
      Enddo
      pp3=p3
 87   Continue      
      do ii=1,iimax
        If(p4.EQ.smallm2(ii)) Then
	   pp4=ii*1D-20
	   goto 86
	Endif
      Enddo
      pp4=p4
 86   Continue      
      do ii=1,iimax
        If(p5.EQ.smallm2(ii)) Then
	   pp5=ii*1D-20
	   goto 85
	Endif
      Enddo
      pp5=p5
 85   Continue      
      do ii=1,iimax
        If(p6.EQ.smallm2(ii)) Then
	   pp6=ii*1D-20
	   goto 84
	Endif
      Enddo
      pp6=p6
 84   Continue      
      do ii=1,iimax
        If(p7.EQ.smallm2(ii)) Then
	   pp7=ii*1D-20
	   goto 83
	Endif
      Enddo
      pp7=p7
 83   Continue      
      do ii=1,iimax
        If(p8.EQ.smallm2(ii)) Then
	   pp8=ii*1D-20
	   goto 82
	Endif
      Enddo
      pp8=p8
 82   Continue      
      do ii=1,iimax
        If(p9.EQ.smallm2(ii)) Then
	   pp9=ii*1D-20
	   goto 81
	Endif
      Enddo
      pp9=p9
 81   Continue      
      do ii=1,iimax
        If(p10.EQ.smallm2(ii)) Then
	   pp10=ii*1D-20
	   goto 80
	Endif
      Enddo
      pp10=p10
 80   Continue      
      
      End
#endif

**************************************************************************************

      double complex function myA0(m1)
      implicit none

      double precision m1
      double complex mass1,cA0f,com_mass,A(0:rmax1),Auv(0:rmax1),mm1
      
#include "looptools.h" 

      mass1=com_mass(m1)
#ifdef LTLib
      myA0=A0(m1)
      Return
#else
#ifdef COLILib
#ifndef DDLib
      myA0=cA0f(mass1)
      Return
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass1(mass1,mm1)
      call A_dd(A,Auv,mm1,0,0)
      myA0=A(0)
      Return
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      myA0=cA0f(mass1)
      call small_mass1(mass1,mm1)
      call A_dd(A,Auv,mm1,0,0)
      call check_LI(myA0,A(0),1)
      myA0=A(0)
      Return
#endif
#endif
#endif

      Print *,"Hallo 1"
      Print *,"no loop library inluded"
      Stop
      					      
      End

**************************************************************************************

      double complex function myB0(p1,m1,m2)
      implicit none

      double precision p1,pp1,m1,m2
      double complex mom1,cB0f,B(0:rmax2,0:rmax2),Buv(0:rmax2,0:rmax2)
      double complex mass1,mass2,com_mass,mm1,mm2
      
#include "looptools.h" 

      mom1=p1
      mass1=com_mass(m1)
      mass2=com_mass(m2)
#ifdef LTLib
      myB0=B0(p1,m1,m2)
      Return
#else 
#ifdef COLILib
#ifndef DDLib
      myB0=cB0f(mom1,mass1,mass2)
      Return
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call B_dd(B,Buv,pp1,mm1,mm2,0,0,0,0)
      myB0=B(0,0)
      Print *,"From myB0: Check implementation of Stefans library first"
      Stop
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      myB0=cB0f(mom1,mass1,mass2)
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call B_dd(B,Buv,pp1,mm1,mm2,0,0,0,0)
      call check_LI(myB0,B(0,0),2)
      myB0=B(0,0)
      Print *,"From myB0: Check implementation of Stefans library first"
      Stop
#endif
#endif
#endif

      Print *,"Hallo 2"
      Print *,"no loop library inluded"
      Stop
      					      
      End

**************************************************************************************

      double complex function myB1(p1,m1,m2)
      implicit none

      double precision p1,pp1,m1,m2
      double complex mom1,cB1f,B(0:rmax2,0:rmax2),Buv(0:rmax2,0:rmax2)
      double complex mass1,mass2,com_mass,mm1,mm2
      
#include "looptools.h" 

      mom1=p1
      mass1=com_mass(m1)
      mass2=com_mass(m2)
#ifdef LTLib
      myB1=B1(p1,m1,m2)
      Return
#else 
#ifdef COLILib
#ifndef DDLib
      myB1=cB1f(mom1,mass1,mass2)
      Return
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call B_dd(B,Buv,pp1,mm1,mm2,0,1,0,0)
      myB1=B(0,1)
      Return
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      myB1=cB1f(mom1,mass1,mass2)
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call B_dd(B,Buv,pp1,mm1,mm2,0,1,0,0)
      call check_LI(myB1,B(0,1),3)
      myB1=B(0,1)
      Return
#endif
#endif
#endif

      Print *,"Hallo 3"
      Print *,"no loop library inluded"
      Stop
      					      
      End

**************************************************************************************

      double complex function myB00(p1,m1,m2)
      implicit none

      double precision p1,pp1,m1,m2
      double complex mom1,cB00f,B(0:rmax2,0:rmax2),Buv(0:rmax2,0:rmax2)
      double complex mass1,mass2,com_mass,mm1,mm2
      
#include "looptools.h" 
            
      mom1=p1
      mass1=com_mass(m1)
      mass2=com_mass(m2)
#ifdef LTLib
      myB00=B00(p1,m1,m2)
      Return
#else 
#ifdef COLILib
#ifndef DDLib
      myB00=cB00f(mom1,mass1,mass2)
      Return
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call B_dd(B,Buv,pp1,mm1,mm2,0,1,0,0)
      myB00=B(2,0)
      Print *,"From myB00: Check implementation of Stefans library first"
      Stop
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      myB00=cB00f(mom1,mass1,mass2)
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call B_dd(B,Buv,pp1,mm1,mm2,0,1,0,0)
      call check_LI(myB00,B(2,0),4)
      myB00=B(2,0)
      Print *,"From myB00: Check implementation of Stefans library first"
      Stop
#endif
#endif
#endif

      Print *,"Hallo 4"
      Print *,"no loop library inluded"
      Stop
      					      
      End

**************************************************************************************

      double complex function myB11(p1,m1,m2)
      implicit none

      double precision p1,pp1,m1,m2
      double complex mom1,cB11f,B(0:rmax2,0:rmax2),Buv(0:rmax2,0:rmax2)
      double complex mass1,mass2,com_mass,mm1,mm2
      
#include "looptools.h" 

      mom1=p1
      mass1=com_mass(m1)
      mass2=com_mass(m2)
#ifdef LTLib
      myB11=B11(p1,m1,m2)
      Return
#else 
#ifdef COLILib
#ifndef DDLib
      myB11=cB11f(mom1,mass1,mass2)
      Return
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call B_dd(B,Buv,pp1,mm1,mm2,0,2,0,0)
      myB11=B(0,2)
      Print *,"From myB11: Check implementation of Stefans library first"
      Stop
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      myB11=cB11f(mom1,mass1,mass2)
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call B_dd(B,Buv,pp1,mm1,mm2,0,2,0,0)
      call check_LI(myB11,B(0,2),5)
      myB11=B(0,2)
      Print *,"From myB11: Check implementation of Stefans library first"
      Stop
#endif
#endif
#endif

      Print *,"Hallo 5"
      Print *,"no loop library inluded"
      Stop
      					      
      End

**************************************************************************************

      double complex function myDB0(p1,m1,m2)
      implicit none

      double precision p1,pp1,m1,m2
      double complex mom1,cDB0f,DB0DD,DB1DD,myDB0test
      double complex mass1,mass2,com_mass,mm1,mm2
      
#include "looptools.h" 

      mom1=p1
      mass1=com_mass(m1)
      mass2=com_mass(m2)
#ifdef LTLib
      myDB0=DB0(p1,m1,m2)
      Return
#else 
#ifdef COLILib
#ifndef DDLib
      myDB0=cDB0f(mom1,mass1,mass2)
      Return
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call DB_dd(DB0DD,DB1DD,pp1,mm1,mm2,0)
      myDB0=DB0DD
      Return
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      myDB0=cDB0f(mom1,mass1,mass2)
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call DB_dd(DB0DD,DB1DD,pp1,mm1,mm2,0)
      call check_LI(myDB0,DB0DD,6)
      myDB0=DB0DD
      Return
#endif
#endif
#endif

      Print *,"Hallo 6"
      Print *,"no loop library inluded"
      Stop
            
      End

**************************************************************************************

      double complex function myDB1(p1,m1,m2)
      implicit none

      double precision p1,pp1,m1,m2
      double complex mom1,cDB1f,DB0DD,DB1DD
      double complex mass1,mass2,com_mass,mm1,mm2
      
#include "looptools.h" 

      mom1=p1
      mass1=com_mass(m1)
      mass2=com_mass(m2)
#ifdef LTLib
      myDB1=DB1(p1,m1,m2)
      Return
#else 
#ifdef COLILib
#ifndef DDLib
      myDB1=cDB1f(mom1,mass1,mass2)
      Return
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call DB_dd(DB0DD,DB1DD,pp1,mm1,mm2,0)
      myDB1=DB1DD
      Return
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      myDB1=cDB1f(mom1,mass1,mass2)
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call DB_dd(DB0DD,DB1DD,pp1,mm1,mm2,0)
      call check_LI(myDB1,DB1DD,7)
      myDB1=DB1DD
      Return
#endif
#endif
#endif

      Print *,"Hallo 7"
      Print *,"no loop library inluded"
      Stop
            
      End

**************************************************************************************

      double complex function myDB00(p1,m1,m2)
      implicit none

      double precision p1,pp1,m1,m2
      double complex mom1,cDB00f,B(0:rmax2,0:rmax2),Buv(0:rmax2,0:rmax2)
      double complex mass1,mass2,com_mass,DB0DD,DB1DD,mm1,mm2,myDB00keep
      
#include "looptools.h" 

      mom1=p1
      mass1=com_mass(m1)
      mass2=com_mass(m2)
#ifdef LTLib
      myDB00=DB00(p1,m1,m2)
      Return
#else 
#ifdef COLILib
#ifndef DDLib
      myDB00=cDB00f(mom1,mass1,mass2)
      Return
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)	      
      call B_dd(B,Buv,pp1,mm1,mm2,0,1,0,0)		      
      call DB_dd(DB0DD,DB1DD,pp1,mm1,mm2,0)		      
      myDB00=-1D0/18D0+1D0/6D0*B(0,1)+1D0/3D0*mass1*DB0DD+    
     &        1D0/6D0*(p1-mass2+mass1)*DB1DD		      
      Return
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      myDB00keep=cDB00f(mom1,mass1,mass2)		      
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)	      
      call B_dd(B,Buv,pp1,mm1,mm2,0,1,0,0)		      
      call DB_dd(DB0DD,DB1DD,pp1,mm1,mm2,0)		      
      myDB00=-1D0/18D0+1D0/6D0*B(0,1)+1D0/3D0*mass1*DB0DD+    
     &        1D0/6D0*(p1-mass2+mass1)*DB1DD		      
      call check_LI(myDB00keep,myDB00,8)		      
      Return
#endif
#endif
#endif

      Print *,"Hallo 8"
      Print *,"no loop library inluded"
      Stop
            
      End

**************************************************************************************

      double complex function myDB11(p1,m1,m2)
      implicit none

      double precision p1,m1,m2
      double complex mom1
      double complex mass1,mass2,com_mass
      
#include "looptools.h" 

#ifdef LTLib
      myDB11=DB11(p1,m1,m2)
      Return
#endif

      myDB11=0D0
C      Print *,"Hallo 9"
C      Print *,"no loop library inluded"
C      Stop
            
      End

**************************************************************************************

      integer function myBget(p1,m1,m2,rank)
      implicit none

      Integer rank
      double precision p1,pp1,m1,m2
      double complex mom1
      double complex mass1,mass2,com_mass
      double complex B0result(10),myDB0,myDB1,myDB00,myDB11
      double complex B(0:rmax2,0:rmax2),Buv(0:rmax2,0:rmax2)
      complex*16 aB0,aB1,aB2(0:1,0:1),aB3(0:1,0:1,0:1)
      complex*16 aB4(0:1,0:1,0:1,0:1),aB5(0:1,0:1,0:1,0:1,0:1)
      complex*16 aB6(0:1,0:1,0:1,0:1,0:1,0:1),mm1,mm2

      
#include "looptools.h" 

      common/B0_result_common/B0result

      mom1=p1
      mass1=com_mass(m1)
      mass2=com_mass(m2)
#ifdef LTLib
      myBget=Bget(p1,m1,m2)
      Return
#else 
#ifdef COLILib
#ifndef DDLib
      call cbp12345(mom1,mass1,mass2,		           
     &  	    aB0,aB1,aB2,aB3,aB4,aB5,aB6,rank,0) 	   
      B0result(1)=aB0				           
      B0result(7)=myDB0(p1,m1,m2)		           
      If(rank.GE.1) Then			           
      B0result(2)=aB1				           
      B0result(8)=myDB1(p1,m1,m2)		           
      If(rank.GE.2) Then			           
      B0result(3)=aB2(0,0)			           
      B0result(4)=aB2(1,1)			           
      B0result(9)=myDB00(p1,m1,m2)		           
      B0result(10)=myDB11(p1,m1,m2)		           
      If(rank.GE.3) Then			           
      B0result(5)=aB3(0,0,1)			           
      B0result(6)=aB3(1,1,1)			           
      Endif					           
      Endif					           
      Endif	
      Return				           
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)
      call B_dd(B,Buv,pp1,mm1,mm2,0,rank,0,0)
      B0result(1)=B(0,0)
      B0result(7)=myDB0(p1,m1,m2)
      If(rank.GE.1) Then
      B0result(2)=B(0,1)
      B0result(8)=myDB1(p1,m1,m2)
      If(rank.GE.2) Then
      B0result(3)=B(2,0)
      B0result(4)=B(0,2)
      B0result(9)=myDB00(p1,m1,m2)
      B0result(10)=myDB11(p1,m1,m2)
      If(rank.GE.3) Then
      B0result(5)=B(2,1)
      B0result(6)=B(0,3)
      Endif
      Endif
      Endif	  	  
      Return
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      call cbp12345(mom1,mass1,mass2,			   
     &  	    aB0,aB1,aB2,aB3,aB4,aB5,aB6,rank,0)    
      call small_mass2(p1,mass1,mass2,pp1,mm1,mm2)	   
      call B_dd(B,Buv,pp1,mm1,mm2,0,rank,0,0)		   
      call check_LI(aB0,B(0,0),9)			   
      B0result(1)=B(0,0)				   
      B0result(7)=myDB0(p1,m1,m2)			   
      If(rank.GE.1) Then				   
      call check_LI(aB1,B(0,1),10)			   
      B0result(2)=B(0,1)				   
      B0result(8)=myDB1(p1,m1,m2)			   
      If(rank.GE.2) Then				   
      call check_LI(aB2(0,0),B(2,0),11) 		   
      B0result(3)=B(2,0)				   
      call check_LI(aB2(1,1),B(0,2),12) 		   
      B0result(4)=B(0,2)				   
      B0result(9)=myDB00(p1,m1,m2)			   
      B0result(10)=myDB11(p1,m1,m2)			   
      If(rank.GE.3) Then				   
      call check_LI(aB3(0,0,1),B(2,1),13)		   
      B0result(5)=B(2,1)				   
      call check_LI(aB3(1,1,1),B(0,3),14)		   
      B0result(6)=B(0,3)				   
      Endif						   
      Endif						   
      Endif						   
      Return
#endif
#endif
#endif

      myBget=0
      Print *,"Hallo 10"
      Print *,"no loop library inluded"
      Stop
      					      
      End

**************************************************************************************

      double complex function myBval(i1,i2)
      implicit none 
      
      integer i1,i2
      double complex B0result(10)
      
#include "looptools.h" 

      common/B0_result_common/B0result
            
#ifdef LTLib
        myBval=Bval(i1,i2)
#else 
        myBval=B0result(i1)
#endif
      
      End


**************************************************************************************

      integer function myCget(p1,p2,p3,m1,m2,m3,rank)
      implicit none

      Integer rank
      double precision p1,p2,p3,m1,m2,m3,pp1,pp2,pp3
      double complex mom1,mom2,mom3,mm1,mm2,mm3
      double complex mass1,mass2,mass3,com_mass
      double complex C0result(22)
      complex*16 C(0:rmax3,0:rmax3,0:rmax3)
      complex*16 Cuv(0:rmax3,0:rmax3,0:rmax3)
      complex*16 aC0,aC1(2),aC2(0:2,0:2),aC3(0:2,0:2,0:2)
      complex*16 aC4(0:2,0:2,0:2,0:2),aC5(0:2,0:2,0:2,0:2,0:2)
      complex*16 aC6(0:2,0:2,0:2,0:2,0:2,0:2)
      
#include "looptools.h" 

      common/C0_result_common/C0result

      mom1=p1
      mom2=p2
      mom3=p3
      mass1=com_mass(m1)
      mass2=com_mass(m2)
      mass3=com_mass(m3)
#ifdef LTLib
      myCget=Cget(p1,p2,p3,m1,m2,m3)
      Return
#else 
#ifdef COLILib
#ifndef DDLib
      call cCp12345(mom1,mom2,mom3,mass1,mass2,mass3,	   
     &  	    aC0,aC1,aC2,aC3,aC4,aC5,aC6,rank,0)    
      C0result(1)=aC0					   
      If(rank.GE.1) Then				   
      C0result(2)=aC1(1)				   
      C0result(3)=aC1(2)				   
      If(rank.GE.2) Then				   
      C0result(4)=aC2(0,0)				   
      C0result(5)=aC2(1,1)				   
      C0result(6)=aC2(1,2)				   
      C0result(7)=aC2(2,2)				   
      If(rank.GE.3) Then				   
      C0result(8)=aC3(0,0,1)				   
      C0result(9)=aC3(0,0,2)				   
      C0result(10)=aC3(1,1,1)				   
      C0result(11)=aC3(1,1,2)				   
      C0result(12)=aC3(1,2,2)				   
      C0result(13)=aC3(2,2,2)				   
      If(rank.GE.4) Then				   
      C0result(14)=aC4(0,0,0,0) 			   
      C0result(15)=aC4(0,0,1,1) 			   
      C0result(16)=aC4(0,0,1,2) 			   
      C0result(17)=aC4(0,0,2,2) 			   
      C0result(18)=aC4(1,1,1,1) 			   
      C0result(19)=aC4(1,1,1,2) 			   
      C0result(20)=aC4(1,1,2,2) 			   
      C0result(21)=aC4(1,2,2,2) 			   
      C0result(22)=aC4(2,2,2,2) 			   
      If(rank.GE.5) Then				   
        Print *,"rank not implemented yet"		   
      Endif						   
      Endif						   
      Endif						   
      Endif						   
      Endif
      Return						   
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      Call small_mass3(p1,p2,p3,mass1,mass2,mass3,pp1,pp2,pp3,mm1,mm2,mm3)    
      Call C_dd(C,Cuv,pp1,pp2,pp3,mm1,mm2,mm3,0,rank,0,0)		      
      C0result(1)=C(0,0,0)						      
      If(rank.GE.1) Then						      
      C0result(2)=C(0,1,0)						      
      C0result(3)=C(0,0,1)						      
      If(rank.GE.2) Then						      
      C0result(4)=C(2,0,0)						      
      C0result(5)=C(0,2,0)						      
      C0result(6)=C(0,1,1)						      
      C0result(7)=C(0,0,2)						      
      If(rank.GE.3) Then						      
      C0result(8)=C(2,1,0)						      
      C0result(9)=C(2,0,1)						      
      C0result(10)=C(0,3,0)						      
      C0result(11)=C(0,2,1)						      
      C0result(12)=C(0,1,2)						      
      C0result(13)=C(0,0,3)						      
      If(rank.GE.4) Then						      
      C0result(14)=C(4,0,0)						      
      C0result(15)=C(2,2,0)						      
      C0result(16)=C(2,1,1)						      
      C0result(17)=C(2,0,2)						      
      C0result(18)=C(0,4,0)						      
      C0result(19)=C(0,3,1)						      
      C0result(20)=C(0,2,2)						      
      C0result(21)=C(0,1,3)						      
      C0result(22)=C(0,0,4)						      
      If(rank.GE.5) Then						      
        Print *,"rank not implemented yet"				      
      Endif								      
      Endif								      
      Endif								      
      Endif								      
      Endif								      
      Return
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      call cCp12345(mom1,mom2,mom3,mass1,mass2,mass3,			      
     &  	    aC0,aC1,aC2,aC3,aC4,aC5,aC6,rank,0) 		      
      Call small_mass3(p1,p2,p3,mass1,mass2,mass3,pp1,pp2,pp3,mm1,mm2,mm3)    
      Call C_dd(C,Cuv,pp1,pp2,pp3,mm1,mm2,mm3,0,rank,0,0)		      
      call check_LI(aC0,C(0,0,0),15)					      
      C0result(1)=C(0,0,0)						      
      If(rank.GE.1) Then						      
      call check_LI(aC1(1),C(0,1,0),16) 				      
      C0result(2)=C(0,1,0)						      
      call check_LI(aC1(2),C(0,0,1),17) 				      
      C0result(3)=C(0,0,1)						      
      If(rank.GE.2) Then						      
      call check_LI(aC2(0,0),C(2,0,0),18)				      
      C0result(4)=C(2,0,0)						      
      call check_LI(aC2(1,1),C(0,2,0),19)				      
      C0result(5)=C(0,2,0)						      
      call check_LI(aC2(1,2),C(0,1,1),20)				      
      C0result(6)=C(0,1,1)						      
      call check_LI(aC2(2,2),C(0,0,2),21)				      
      C0result(7)=C(0,0,2)						      
      If(rank.GE.3) Then						      
      call check_LI(aC3(0,0,1),C(2,1,0),22)				      
      C0result(8)=C(2,1,0)						      
      call check_LI(aC3(0,0,2),C(2,0,1),23)				      
      C0result(9)=C(2,0,1)						      
      call check_LI(aC3(1,1,1),C(0,3,0),24)				      
      C0result(10)=C(0,3,0)						      
      call check_LI(aC3(1,1,2),C(0,2,1),25)				      
      C0result(11)=C(0,2,1)						      
      call check_LI(aC3(1,2,2),C(0,1,2),26)				      
      C0result(12)=C(0,1,2)						      
      call check_LI(aC3(2,2,2),C(0,0,3),27)				      
      C0result(13)=C(0,0,3)						      
      If(rank.GE.4) Then						      
      call check_LI(aC4(0,0,0,0),C(4,0,0),28)				      
      C0result(14)=C(4,0,0)						      
      call check_LI(aC4(0,0,1,1),C(2,2,0),29)				      
      C0result(15)=C(2,2,0)						      
      call check_LI(aC4(0,0,1,2),C(2,1,1),30)				      
      C0result(16)=C(2,1,1)						      
      call check_LI(aC4(0,0,2,2),C(2,0,2),31)				      
      C0result(17)=C(2,0,2)						      
      call check_LI(aC4(1,1,1,1),C(0,4,0),32)				      
      C0result(18)=C(0,4,0)						      
      call check_LI(aC4(1,1,1,2),C(0,3,1),33)				      
      C0result(19)=C(0,3,1)						      
      call check_LI(aC4(1,1,2,2),C(0,2,2),34)				      
      C0result(20)=C(0,2,2)						      
      call check_LI(aC4(1,2,2,2),C(0,1,3),35)				      
      C0result(21)=C(0,1,3)						      
      call check_LI(aC4(2,2,2,2),C(0,0,4),36)				      
      C0result(22)=C(0,0,4)						      
      If(rank.GE.5) Then						      
        Print *,"rank not implemented yet"				      
      Endif								      
      Endif								      
      Endif								      
      Endif								      
      Endif	
      Return							      
#endif
#endif
#endif

      myCget=0
      Print *,"Hallo 11"
      Print *,"no loop library inluded"
      Stop
      					      
      End

**************************************************************************************

      double complex function myCval(i1,i2)
      implicit none 
      
      integer i1,i2
      double complex C0result(22)
      
#include "looptools.h" 

      common/C0_result_common/C0result
            
#ifdef LTLib
        myCval=Cval(i1,i2)
#else 
        myCval=C0result(i1)
#endif
      
      End

**************************************************************************************

      integer function myDget(p1,p2,p3,p4,p5,p6,m1,m2,m3,m4,rank)
      implicit none

      Integer rank
      double precision p1,p2,p3,p4,p5,p6,m1,m2,m3,m4,pp1,pp2,pp3,pp4,pp5,pp6
      double complex mom1,mom2,mom3,mom4,mom5,mom6,mm1,mm2,mm3,mm4
      double complex mass1,mass2,mass3,mass4,com_mass
      double complex D0result(80)
      complex*16 D(0:rmax4,0:rmax4,0:rmax4,0:rmax4)
      complex*16 Duv(0:rmax4,0:rmax4,0:rmax4,0:rmax4)
      complex*16 aD0,aD1(3),aD2(0:3,0:3),aD3(0:3,0:3,0:3)
      complex*16 aD4(0:3,0:3,0:3,0:3),aD5(0:3,0:3,0:3,0:3,0:3)
      
#include "looptools.h" 

      common/D0_result_common/D0result

      mom1=p1
      mom2=p2
      mom3=p3
      mom4=p4
      mom5=p5
      mom6=p6
      mass1=com_mass(m1)
      mass2=com_mass(m2)
      mass3=com_mass(m3)
      mass4=com_mass(m4)
#ifdef LTLib
      myDget=Dget(p1,p2,p3,p4,p5,p6,m1,m2,m3,m4)
      Return
#else
#ifdef COLILib
#ifndef DDLib
      call cDp12345(mom1,mom2,mom3,mom4,mom5,mom6,mass1,mass2,mass3,mass4,    
     &  	    aD0,aD1,aD2,aD3,aD4,aD5,rank,0)			      
      D0result(1)=aD0							      
      If(rank.GE.1) Then						      
      D0result(2)=aD1(1)						      
      D0result(3)=aD1(2)						      
      D0result(4)=aD1(3)						      
      If(rank.GE.2) Then						      
      D0result(5)=aD2(0,0)						      
      D0result(6)=aD2(1,1)						      
      D0result(7)=aD2(1,2)						      
      D0result(8)=aD2(1,3)						      
      D0result(9)=aD2(2,2)						      
      D0result(10)=aD2(2,3)						      
      D0result(11)=aD2(3,3)						      
      If(rank.GE.3) Then						      
      D0result(12)=aD3(0,0,1)						      
      D0result(13)=aD3(0,0,2)						      
      D0result(14)=aD3(0,0,3)						      
      D0result(15)=aD3(1,1,1)						      
      D0result(16)=aD3(1,1,2)						      
      D0result(17)=aD3(1,1,3)						      
      D0result(18)=aD3(1,2,2)						      
      D0result(19)=aD3(1,2,3)						      
      D0result(20)=aD3(1,3,3)						      
      D0result(21)=aD3(2,2,2)						      
      D0result(22)=aD3(2,2,3)						      
      D0result(23)=aD3(2,3,3)						      
      D0result(24)=aD3(3,3,3)						      
      If(rank.GE.4) Then						      
      D0result(25)=aD4(0,0,0,0) 					      
      D0result(26)=aD4(0,0,1,1) 					      
      D0result(27)=aD4(0,0,1,2) 					      
      D0result(28)=aD4(0,0,1,3) 					      
      D0result(29)=aD4(0,0,2,2) 					      
      D0result(30)=aD4(0,0,2,3) 					      
      D0result(31)=aD4(0,0,3,3) 					      
      D0result(32)=aD4(1,1,1,1) 					      
      D0result(33)=aD4(1,1,1,2) 					      
      D0result(34)=aD4(1,1,1,3) 					      
      D0result(35)=aD4(1,1,2,2) 					      
      D0result(36)=aD4(1,1,2,3) 					      
      D0result(37)=aD4(1,1,3,3) 					      
      D0result(38)=aD4(1,2,2,2) 					      
      D0result(39)=aD4(1,2,2,3) 					      
      D0result(40)=aD4(1,2,3,3) 					      
      D0result(41)=aD4(1,3,3,3) 					      
      D0result(42)=aD4(2,2,2,2) 					      
      D0result(43)=aD4(2,2,2,3) 					      
      D0result(44)=aD4(2,2,3,3) 					      
      D0result(45)=aD4(2,3,3,3) 					      
      D0result(46)=aD4(3,3,3,3) 					      
      If(rank.GE.5) Then						      
      D0result(47)=aD5(0,0,0,0,1)					      
      D0result(48)=aD5(0,0,0,0,2)					      
      D0result(49)=aD5(0,0,0,0,3)					      
      D0result(50)=aD5(0,0,1,1,1)					      
      D0result(51)=aD5(0,0,1,1,2)					      
      D0result(52)=aD5(0,0,1,1,3)					      
      D0result(53)=aD5(0,0,1,2,2)					      
      D0result(54)=aD5(0,0,1,2,3)					      
      D0result(55)=aD5(0,0,1,3,3)					      
      D0result(56)=aD5(0,0,2,2,2)					      
      D0result(57)=aD5(0,0,2,2,3)					      
      D0result(58)=aD5(0,0,2,3,3)					      
      D0result(59)=aD5(0,0,3,3,3)					      
      D0result(60)=aD5(1,1,1,1,1)					      
      D0result(61)=aD5(1,1,1,1,2)					      
      D0result(62)=aD5(1,1,1,1,3)					      
      D0result(63)=aD5(1,1,1,2,2)					      
      D0result(64)=aD5(1,1,1,2,3)					      
      D0result(65)=aD5(1,1,1,3,3)					      
      D0result(66)=aD5(1,1,2,2,2)					      
      D0result(67)=aD5(1,1,2,2,3)					      
      D0result(68)=aD5(1,1,2,3,3)					      
      D0result(69)=aD5(1,1,3,3,3)					      
      D0result(70)=aD5(1,2,2,2,2)					      
      D0result(71)=aD5(1,2,2,2,3)					      
      D0result(72)=aD5(1,2,2,3,3)					      
      D0result(73)=aD5(1,2,3,3,3)					      
      D0result(74)=aD5(1,3,3,3,3)					      
      D0result(75)=aD5(2,2,2,2,2)					      
      D0result(76)=aD5(2,2,2,2,3)					      
      D0result(77)=aD5(2,2,2,3,3)					      
      D0result(78)=aD5(2,2,3,3,3)					      
      D0result(79)=aD5(2,3,3,3,3)					      
      D0result(80)=aD5(3,3,3,3,3)					      
      If(rank.GE.6) Then						      
        Print *,"rank not implemented yet"				      
      Endif								      
      Endif								      
      Endif								      
      Endif								      
      Endif								      
      Endif								      
      Return
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass4(p1,p2,p3,p4,p5,p6,mass1,mass2,mass3,mass4,    
     &  		 pp1,pp2,pp3,pp4,pp5,pp6,mm1,mm2,mm3,mm4)    
      call D_dd(D,Duv,pp1,pp2,pp3,pp4,pp5,pp6,mm1,mm2,mm3,mm4,       
     &  	      rank,1,0) 				     
      D0result(1)=D(0,0,0,0)					     
      If(rank.GE.1) Then					     
      D0result(2)=D(0,1,0,0)					     
      D0result(3)=D(0,0,1,0)					     
      D0result(4)=D(0,0,0,1)					     
      If(rank.GE.2) Then					     
      D0result(5)=D(2,0,0,0)					     
      D0result(6)=D(0,2,0,0)					     
      D0result(7)=D(0,1,1,0)					     
      D0result(8)=D(0,1,0,1)					     
      D0result(9)=D(0,0,2,0)					     
      D0result(10)=D(0,0,1,1)					     
      D0result(11)=D(0,0,0,2)					     
      If(rank.GE.3) Then					     
      D0result(12)=D(2,1,0,0)					     
      D0result(13)=D(2,0,1,0)					     
      D0result(14)=D(2,0,0,1)					     
      D0result(15)=D(0,3,0,0)					     
      D0result(16)=D(0,2,1,0)					     
      D0result(17)=D(0,2,0,1)					     
      D0result(18)=D(0,1,2,0)					     
      D0result(19)=D(0,1,1,1)					     
      D0result(20)=D(0,1,0,2)					     
      D0result(21)=D(0,0,3,0)					     
      D0result(22)=D(0,0,2,1)					     
      D0result(23)=D(0,0,1,2)					     
      D0result(24)=D(0,0,0,3)					     
      If(rank.GE.4) Then					     
      D0result(25)=D(4,0,0,0)					     
      D0result(26)=D(2,2,0,0)					     
      D0result(27)=D(2,1,1,0)					     
      D0result(28)=D(2,1,0,1)					     
      D0result(29)=D(2,0,2,0)					     
      D0result(30)=D(2,0,1,1)					     
      D0result(31)=D(2,0,0,2)					     
      D0result(32)=D(0,4,0,0)					     
      D0result(33)=D(0,3,1,0)					     
      D0result(34)=D(0,3,0,1)					     
      D0result(35)=D(0,2,2,0)					     
      D0result(36)=D(0,2,1,1)					     
      D0result(37)=D(0,2,0,2)					     
      D0result(38)=D(0,1,3,0)					     
      D0result(39)=D(0,1,2,1)					     
      D0result(40)=D(0,1,1,2)					     
      D0result(41)=D(0,1,0,3)					     
      D0result(42)=D(0,0,4,0)					     
      D0result(43)=D(0,0,3,1)					     
      D0result(44)=D(0,0,2,2)					     
      D0result(45)=D(0,0,1,3)					     
      D0result(46)=D(0,0,0,4)					     
      If(rank.GE.5) Then					     
      D0result(47)=D(4,1,0,0)					     
      D0result(48)=D(4,0,1,0)					     
      D0result(49)=D(4,0,0,1)					     
      D0result(50)=D(2,3,0,0)					     
      D0result(51)=D(2,2,1,0)					     
      D0result(52)=D(2,2,0,1)					     
      D0result(53)=D(2,1,2,0)					     
      D0result(54)=D(2,1,1,1)					     
      D0result(55)=D(2,1,0,2)					     
      D0result(56)=D(2,0,3,0)					     
      D0result(57)=D(2,0,2,1)					     
      D0result(58)=D(2,0,1,2)					     
      D0result(59)=D(2,0,0,3)					     
      D0result(60)=D(0,5,0,0)					     
      D0result(61)=D(0,4,1,0)					     
      D0result(62)=D(0,4,0,1)					     
      D0result(63)=D(0,3,2,0)					     
      D0result(64)=D(0,3,1,1)					     
      D0result(65)=D(0,3,0,2)					     
      D0result(66)=D(0,2,3,0)					     
      D0result(67)=D(0,2,2,1)					     
      D0result(68)=D(0,2,1,2)					     
      D0result(69)=D(0,2,0,3)					     
      D0result(70)=D(0,1,4,0)					     
      D0result(71)=D(0,1,3,1)					     
      D0result(72)=D(0,1,2,2)					     
      D0result(73)=D(0,1,1,3)					     
      D0result(74)=D(0,1,0,4)					     
      D0result(75)=D(0,0,5,0)					     
      D0result(76)=D(0,0,4,1)					     
      D0result(77)=D(0,0,3,2)					     
      D0result(78)=D(0,0,2,3)					     
      D0result(79)=D(0,0,1,4)					     
      D0result(80)=D(0,0,0,5)					     
      If(rank.GE.6) Then					     
        Print *,"rank not implemented yet"			     
      Endif							     
      Endif							     
      Endif							     
      Endif							     
      Endif							     
      Endif							     
      Return
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      call cDp12345(mom1,mom2,mom3,mom4,mom5,mom6,mass1,mass2,mass3,mass4,    
     &  	    aD0,aD1,aD2,aD3,aD4,aD5,rank,0)			      
      call small_mass4(p1,p2,p3,p4,p5,p6,mass1,mass2,mass3,mass4,	      
     &  		 pp1,pp2,pp3,pp4,pp5,pp6,mm1,mm2,mm3,mm4)	      
      call D_dd(D,Duv,pp1,pp2,pp3,pp4,pp5,pp6,mm1,mm2,mm3,mm4,  	      
     &  	      rank,1,0) 					      
      call check_LI(aD0,D(0,0,0,0),37)  				      
      D0result(1)=D(0,0,0,0)						      
      If(rank.GE.1) Then						      
      call check_LI(aD1(1),D(0,1,0,0),37)				      
      D0result(2)=D(0,1,0,0)						      
      call check_LI(aD1(2),D(0,0,1,0),38)				      
      D0result(3)=D(0,0,1,0)						      
      call check_LI(aD1(3),D(0,0,0,1),39)				      
      D0result(4)=D(0,0,0,1)						      
      If(rank.GE.2) Then						      
      call check_LI(aD2(0,0),D(2,0,0,0),40)				      
      D0result(5)=D(2,0,0,0)						      
      call check_LI(aD2(1,1),D(0,2,0,0),41)				      
      D0result(6)=D(0,2,0,0)						      
      call check_LI(aD2(1,2),D(0,1,1,0),42)				      
      D0result(7)=D(0,1,1,0)						      
      call check_LI(aD2(1,3),D(0,1,0,1),43)				      
      D0result(8)=D(0,1,0,1)						      
      call check_LI(aD2(2,2),D(0,0,2,0),44)				      
      D0result(9)=D(0,0,2,0)						      
      call check_LI(aD2(2,3),D(0,0,1,1),45)				      
      D0result(10)=D(0,0,1,1)						      
      call check_LI(aD2(3,3),D(0,0,0,2),46)				      
      D0result(11)=D(0,0,0,2)						      
      If(rank.GE.3) Then						      
      call check_LI(aD3(0,0,1),D(2,1,0,0),47)				      
      D0result(12)=D(2,1,0,0)						      
      call check_LI(aD3(0,0,2),D(2,0,1,0),48)				      
      D0result(13)=D(2,0,1,0)						      
      call check_LI(aD3(0,0,3),D(2,0,0,1),49)				      
      D0result(14)=D(2,0,0,1)						      
      call check_LI(aD3(1,1,1),D(0,3,0,0),50)				      
      D0result(15)=D(0,3,0,0)						      
      call check_LI(aD3(1,1,2),D(0,2,1,0),51)				      
      D0result(16)=D(0,2,1,0)						      
      call check_LI(aD3(1,1,3),D(0,2,0,1),52)				      
      D0result(17)=D(0,2,0,1)						      
      call check_LI(aD3(1,2,2),D(0,1,2,0),53)				      
      D0result(18)=D(0,1,2,0)						      
      call check_LI(aD3(1,2,3),D(0,1,1,1),54)				      
      D0result(19)=D(0,1,1,1)						      
      call check_LI(aD3(1,3,3),D(0,1,0,2),55)				      
      D0result(20)=D(0,1,0,2)						      
      call check_LI(aD3(2,2,2),D(0,0,3,0),56)				      
      D0result(21)=D(0,0,3,0)						      
      call check_LI(aD3(2,2,3),D(0,0,2,1),57)				      
      D0result(22)=D(0,0,2,1)						      
      call check_LI(aD3(2,3,3),D(0,0,1,2),58)				      
      D0result(23)=D(0,0,1,2)						      
      call check_LI(aD3(3,3,3),D(0,0,0,3),59)				      
      D0result(24)=D(0,0,0,3)						      
      If(rank.GE.4) Then						      
      call check_LI(aD4(0,0,0,0),D(0,0,0,0),60) 			      
      D0result(25)=D(4,0,0,0)						      
      call check_LI(aD4(0,0,1,1),D(2,2,0,0),61) 			      
      D0result(26)=D(2,2,0,0)						      
      call check_LI(aD4(0,0,1,2),D(2,1,1,0),62) 			      
      D0result(27)=D(2,1,1,0)						      
      call check_LI(aD4(0,0,1,3),D(2,1,0,1),63) 			      
      D0result(28)=D(2,1,0,1)						      
      call check_LI(aD4(0,0,2,2),D(2,0,2,0),64) 			      
      D0result(29)=D(2,0,2,0)						      
      call check_LI(aD4(0,0,2,3),D(2,0,1,1),65) 			      
      D0result(30)=D(2,0,1,1)						      
      call check_LI(aD4(0,0,3,3),D(2,0,0,2),66) 			      
      D0result(31)=D(2,0,0,2)						      
      call check_LI(aD4(1,1,1,1),D(0,4,0,0),67) 			      
      D0result(32)=D(0,4,0,0)						      
      call check_LI(aD4(1,1,1,2),D(0,3,1,0),68) 			      
      D0result(33)=D(0,3,1,0)						      
      call check_LI(aD4(1,1,1,3),D(0,3,0,1),69) 			      
      D0result(34)=D(0,3,0,1)						      
      call check_LI(aD4(1,1,2,2),D(0,2,2,0),70) 			      
      D0result(35)=D(0,2,2,0)						      
      call check_LI(aD4(1,1,2,3),D(0,2,1,1),71) 			      
      D0result(36)=D(0,2,1,1)						      
      call check_LI(aD4(1,1,3,3),D(0,2,0,2),72) 			      
      D0result(37)=D(0,2,0,2)						      
      call check_LI(aD4(1,2,2,2),D(0,1,3,0),73) 			      
      D0result(38)=D(0,1,3,0)						      
      call check_LI(aD4(1,2,2,3),D(0,1,2,1),74) 			      
      D0result(39)=D(0,1,2,1)						      
      call check_LI(aD4(1,2,3,3),D(0,1,1,2),75) 			      
      D0result(40)=D(0,1,1,2)						      
      call check_LI(aD4(1,3,3,3),D(0,1,0,3),76) 			      
      D0result(41)=D(0,1,0,3)						      
      call check_LI(aD4(2,2,2,2),D(0,0,4,0),77) 			      
      D0result(42)=D(0,0,4,0)						      
      call check_LI(aD4(2,2,2,3),D(0,0,3,1),78) 			      
      D0result(43)=D(0,0,3,1)						      
      call check_LI(aD4(2,2,3,3),D(0,0,2,2),79) 			      
      D0result(44)=D(0,0,2,2)						      
      call check_LI(aD4(2,3,3,3),D(0,0,1,3),80) 			      
      D0result(45)=D(0,0,1,3)						      
      call check_LI(aD4(3,3,3,3),D(0,0,0,4),81) 			      
      D0result(46)=D(0,0,0,4)						      
      If(rank.GE.5) Then						      
      call check_LI(aD5(0,0,0,0,1),D(4,1,0,0),82)			      
      D0result(47)=D(4,1,0,0)						      
      call check_LI(aD5(0,0,0,0,2),D(4,0,1,0),83)			      
      D0result(48)=D(4,0,1,0)						      
      call check_LI(aD5(0,0,0,0,3),D(4,0,0,1),84)			      
      D0result(49)=D(4,0,0,1)						      
      call check_LI(aD5(0,0,1,1,1),D(2,3,0,0),85)			      
      D0result(50)=D(2,3,0,0)						      
      call check_LI(aD5(0,0,1,1,2),D(2,2,1,0),86)			      
      D0result(51)=D(2,2,1,0)						      
      call check_LI(aD5(0,0,1,1,3),D(2,2,0,1),87)			      
      D0result(52)=D(2,2,0,1)						      
      call check_LI(aD5(0,0,1,2,2),D(2,1,2,0),88)			      
      D0result(53)=D(2,1,2,0)						      
      call check_LI(aD5(0,0,1,2,3),D(2,1,1,1),89)			      
      D0result(54)=D(2,1,1,1)						      
      call check_LI(aD5(0,0,1,3,3),D(2,1,0,2),90)			      
      D0result(55)=D(2,1,0,2)						      
      call check_LI(aD5(0,0,2,2,2),D(2,0,3,0),91)			      
      D0result(56)=D(2,0,3,0)						      
      call check_LI(aD5(0,0,2,2,3),D(2,0,2,1),92)			      
      D0result(57)=D(2,0,2,1)						      
      call check_LI(aD5(0,0,2,3,3),D(2,0,1,2),93)			      
      D0result(58)=D(2,0,1,2)						      
      call check_LI(aD5(0,0,3,3,3),D(2,0,0,2),94)			      
      D0result(59)=D(2,0,0,3)						      
      call check_LI(aD5(1,1,1,1,1),D(0,5,0,0),95)			      
      D0result(60)=D(0,5,0,0)						      
      call check_LI(aD5(1,1,1,1,2),D(0,4,1,0),96)			      
      D0result(61)=D(0,4,1,0)						      
      call check_LI(aD5(1,1,1,1,3),D(0,4,0,1),97)			      
      D0result(62)=D(0,4,0,1)						      
      call check_LI(aD5(1,1,1,2,2),D(0,3,2,0),98)			      
      D0result(63)=D(0,3,2,0)						      
      call check_LI(aD5(1,1,1,2,3),D(0,3,1,1),99)			      
      D0result(64)=D(0,3,1,1)						      
      call check_LI(aD5(1,1,1,3,3),D(0,3,0,2),100)			      
      D0result(65)=D(0,3,0,2)						      
      call check_LI(aD5(1,1,2,2,2),D(0,2,3,0),101)			      
      D0result(66)=D(0,2,3,0)						      
      call check_LI(aD5(1,1,2,2,3),D(0,2,2,1),102)			      
      D0result(67)=D(0,2,2,1)						      
      call check_LI(aD5(1,1,2,3,3),D(0,2,1,2),103)			      
      D0result(68)=D(0,2,1,2)						      
      call check_LI(aD5(1,1,3,3,3),D(0,2,0,3),104)			      
      D0result(69)=D(0,2,0,3)						      
      call check_LI(aD5(1,2,2,2,2),D(0,1,4,0),105)			      
      D0result(70)=D(0,1,4,0)						      
      call check_LI(aD5(1,2,2,2,3),D(0,1,3,1),106)			      
      D0result(71)=D(0,1,3,1)						      
      call check_LI(aD5(1,2,2,3,3),D(0,1,2,2),107)			      
      D0result(72)=D(0,1,2,2)						      
      call check_LI(aD5(1,2,3,3,3),D(0,1,1,3),108)			      
      D0result(73)=D(0,1,1,3)						      
      call check_LI(aD5(1,3,3,3,3),D(0,1,0,4),109)			      
      D0result(74)=D(0,1,0,4)						      
      call check_LI(aD5(2,2,2,2,2),D(0,0,5,0),110)			      
      D0result(75)=D(0,0,5,0)						      
      call check_LI(aD5(2,2,2,2,3),D(0,0,4,1),111)			      
      D0result(76)=D(0,0,4,1)						      
      call check_LI(aD5(2,2,2,3,3),D(0,0,3,2),112)			      
      D0result(77)=D(0,0,3,2)						      
      call check_LI(aD5(2,2,3,3,3),D(0,0,2,3),113)			      
      D0result(78)=D(0,0,2,3)						      
      call check_LI(aD5(2,3,3,3,3),D(0,0,1,4),114)			      
      D0result(79)=D(0,0,1,4)						      
      call check_LI(aD5(3,3,3,3,3),D(0,0,0,5),115)			      
      D0result(80)=D(0,0,0,5)						      
      If(rank.GE.6) Then						      
        Print *,"rank not implemented yet"				      
      Endif								      
      Endif								      
      Endif								      
      Endif								      
      Endif								      
      Endif								      
      Return
#endif
#endif
#endif

      myDget=0
      Print *,"Hallo 12"
      Print *,"no loop library inluded"
      Stop
      					      
      End

**************************************************************************************

      double complex function myDval(i1,i2)
      implicit none 
      
      integer i1,i2
      double complex D0result(80)
      
#include "looptools.h" 

      common/D0_result_common/D0result
            
#ifdef LTLib
        myDval=Dval(i1,i2)
#else 
        myDval=D0result(i1)
#endif
      
      End

**************************************************************************************

      integer function myEget(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,m1,m2,m3,m4,m5,rank)
      implicit none

      Integer rank
      double precision p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,m1,m2,m3,m4,m5
      double precision pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10
      double complex mom1,mom2,mom3,mom4,mom5,mom6,mom7,mom8,mom9,mom10
      double complex mm1,mm2,mm3,mm4,mm5
      double complex mass1,mass2,mass3,mass4,mass5,com_mass
      double complex E0result(86)
      complex*16 E(0:rmax5,0:rmax5,0:rmax5,0:rmax5,0:rmax5)
      complex*16 aEn0,aEn1(4),aEn2(0:4,0:4),aEn3(0:4,0:4,0:4)
      complex*16 aEn4(0:4,0:4,0:4,0:4)
      
#include "looptools.h" 

      common/E0_result_common/E0result

      mom1=p1
      mom2=p2
      mom3=p3
      mom4=p4
      mom5=p5
      mom6=p6
      mom7=p7
      mom8=p8
      mom9=p9
      mom10=p10
      mass1=com_mass(m1)
      mass2=com_mass(m2)
      mass3=com_mass(m3)
      mass4=com_mass(m4)
      mass5=com_mass(m5)
#ifdef LTLib
      myEget=Eget(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,m1,m2,m3,m4,m5)
      Return
#else 
#ifdef COLILib
#ifndef DDLib
      call cEp1234(mom1,mom2,mom3,mom4,mom5,mom6,mom7,mom8,mom9,mom10,  		    
     &  	    mass1,mass2,mass3,mass4,mass5,aEn0,aEn1,aEn2,aEn3,aEn4,rank)	    
      E0result(1)=aEn0  								    
      If(rank.GE.1) Then								    
      E0result(2)=aEn1(1)								    
      E0result(3)=aEn1(2)								    
      E0result(4)=aEn1(3)								    
      E0result(5)=aEn1(4)								    
      If(rank.GE.2) Then								    
      E0result(6)=aEn2(0,0)								    
      E0result(7)=aEn2(1,1)								    
      E0result(8)=aEn2(1,2)								    
      E0result(9)=aEn2(1,3)								    
      E0result(10)=aEn2(1,4)								    
      E0result(11)=aEn2(2,2)								    
      E0result(12)=aEn2(2,3)								    
      E0result(13)=aEn2(2,4)								    
      E0result(14)=aEn2(3,3)								    
      E0result(15)=aEn2(3,4)								    
      E0result(16)=aEn2(4,4)								    
      If(rank.GE.3) Then								    
      E0result(17)=aEn3(0,0,1)  							    
      E0result(18)=aEn3(0,0,2)  							    
      E0result(19)=aEn3(0,0,3)  							    
      E0result(20)=aEn3(0,0,4)  							    
      E0result(21)=aEn3(1,1,1)  							    
      E0result(22)=aEn3(1,1,2)  							    
      E0result(23)=aEn3(1,1,3)  							    
      E0result(24)=aEn3(1,1,4)  							    
      E0result(25)=aEn3(1,2,2)  							    
      E0result(26)=aEn3(1,2,3)  							    
      E0result(27)=aEn3(1,2,4)  							    
      E0result(28)=aEn3(1,3,3)  							    
      E0result(29)=aEn3(1,3,4)  							    
      E0result(30)=aEn3(1,4,4)  							    
      E0result(31)=aEn3(2,2,2)  							    
      E0result(32)=aEn3(2,2,3)  							    
      E0result(33)=aEn3(2,2,4)  							    
      E0result(34)=aEn3(2,3,3)  							    
      E0result(35)=aEn3(2,3,4)  							    
      E0result(36)=aEn3(2,4,4)  							    
      E0result(37)=aEn3(3,3,3)  							    
      E0result(38)=aEn3(3,3,4)  							    
      E0result(39)=aEn3(3,4,4)  							    
      E0result(40)=aEn3(4,4,4)  							    
      If(rank.GE.4) Then								    
      E0result(41)=aEn4(0,0,0,0)							    
      E0result(42)=aEn4(0,0,1,1)							    
      E0result(43)=aEn4(0,0,1,2)							    
      E0result(44)=aEn4(0,0,1,3)							    
      E0result(45)=aEn4(0,0,1,4)							    
      E0result(46)=aEn4(0,0,2,2)							    
      E0result(47)=aEn4(0,0,2,3)							    
      E0result(48)=aEn4(0,0,2,4)							    
      E0result(49)=aEn4(0,0,3,3)							    
      E0result(50)=aEn4(0,0,3,4)							    
      E0result(51)=aEn4(0,0,4,4)							    
      E0result(52)=aEn4(1,1,1,1)							    
      E0result(53)=aEn4(1,1,1,2)							    
      E0result(54)=aEn4(1,1,1,3)							    
      E0result(55)=aEn4(1,1,1,4)							    
      E0result(56)=aEn4(1,1,2,2)							    
      E0result(57)=aEn4(1,1,2,3)							    
      E0result(58)=aEn4(1,1,2,4)							    
      E0result(59)=aEn4(1,1,3,3)							    
      E0result(60)=aEn4(1,1,3,4)							    
      E0result(61)=aEn4(1,1,4,4)							    
      E0result(62)=aEn4(1,2,2,2)							    
      E0result(63)=aEn4(1,2,2,3)							    
      E0result(64)=aEn4(1,2,2,4)							    
      E0result(65)=aEn4(1,2,3,3)							    
      E0result(66)=aEn4(1,2,3,4)							    
      E0result(67)=aEn4(1,2,4,4)							    
      E0result(68)=aEn4(1,3,3,3)							    
      E0result(69)=aEn4(1,3,3,4)							    
      E0result(70)=aEn4(1,3,4,4)							    
      E0result(71)=aEn4(1,4,4,4)							    
      E0result(72)=aEn4(2,2,2,2)							    
      E0result(73)=aEn4(2,2,2,3)							    
      E0result(74)=aEn4(2,2,2,4)							    
      E0result(75)=aEn4(2,2,3,3)							    
      E0result(76)=aEn4(2,2,3,4)							    
      E0result(77)=aEn4(2,2,4,4)							    
      E0result(78)=aEn4(2,3,3,3)							    
      E0result(79)=aEn4(2,3,3,4)							    
      E0result(80)=aEn4(2,3,4,4)							    
      E0result(81)=aEn4(2,4,4,4)							    
      E0result(82)=aEn4(3,3,3,3)							    
      E0result(83)=aEn4(3,3,3,4)							    
      E0result(84)=aEn4(3,3,4,4)							    
      E0result(85)=aEn4(3,4,4,4)							    
      E0result(86)=aEn4(4,4,4,4)							    
      If(rank.GE.5) Then								    
        Print *,"rank not implemented yet"						    
      Endif										    
      Endif										    
      Endif										    
      Endif										    
      Endif										    
      Return
#endif
#endif
#ifndef COLILib
#ifdef DDLib
      call small_mass5(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,  	    
     &  	       mass1,mass2,mass3,mass4,mass5,		    
     &  	       pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10,    
     &  	       mm1,mm2,mm3,mm4,mm5)			    
      call E_dd(E,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10,	    
     &  	    mm1,mm2,mm3,mm4,mm5,rank,0) 		    
      E0result(1)=E(0,0,0,0,0)  				    
      If(rank.GE.1) Then					    
      E0result(2)=E(0,1,0,0,0)  				    
      E0result(3)=E(0,0,1,0,0)  				    
      E0result(4)=E(0,0,0,1,0)  				    
      E0result(5)=E(0,0,0,0,1)  				    
      If(rank.GE.2) Then					    
      E0result(6)=E(2,0,0,0,0)  				    
      E0result(7)=E(0,2,0,0,0)  				    
      E0result(8)=E(0,1,1,0,0)  				    
      E0result(9)=E(0,1,0,1,0)  				    
      E0result(10)=E(0,1,0,0,1) 				    
      E0result(11)=E(0,0,2,0,0) 				    
      E0result(12)=E(0,0,1,1,0) 				    
      E0result(13)=E(0,0,1,0,1) 				    
      E0result(14)=E(0,0,0,2,0) 				    
      E0result(15)=E(0,0,0,1,1) 				    
      E0result(16)=E(0,0,0,0,2) 				    
      If(rank.GE.3) Then					    
      E0result(17)=E(2,1,0,0,0) 				    
      E0result(18)=E(2,0,1,0,0) 				    
      E0result(19)=E(2,0,0,1,0) 				    
      E0result(20)=E(2,0,0,0,1) 				    
      E0result(21)=E(0,3,0,0,0) 				    
      E0result(22)=E(0,2,1,0,0) 				    
      E0result(23)=E(0,2,0,1,0) 				    
      E0result(24)=E(0,2,0,0,1) 				    
      E0result(25)=E(0,1,2,0,0) 				    
      E0result(26)=E(0,1,1,1,0) 				    
      E0result(27)=E(0,1,1,0,1) 				    
      E0result(28)=E(0,1,0,2,0) 				    
      E0result(29)=E(0,1,0,1,1) 				    
      E0result(30)=E(0,1,0,0,2) 				    
      E0result(31)=E(0,0,3,0,0) 				    
      E0result(32)=E(0,0,2,1,0) 				    
      E0result(33)=E(0,0,2,0,1) 				    
      E0result(34)=E(0,0,1,2,0) 				    
      E0result(35)=E(0,0,1,1,1) 				    
      E0result(36)=E(0,0,1,0,2) 				    
      E0result(37)=E(0,0,0,3,0) 				    
      E0result(38)=E(0,0,0,2,1) 				    
      E0result(39)=E(0,0,0,1,2) 				    
      E0result(40)=E(0,0,0,0,3) 				    
      If(rank.GE.4) Then					    
      E0result(41)=E(4,0,0,0,0) 				    
      E0result(42)=E(2,2,0,0,0) 				    
      E0result(43)=E(2,1,1,0,0) 				    
      E0result(44)=E(2,1,0,1,0) 				    
      E0result(45)=E(2,1,0,0,1) 				    
      E0result(46)=E(2,0,2,0,0) 				    
      E0result(47)=E(2,0,1,1,0) 				    
      E0result(48)=E(2,0,1,0,1) 				    
      E0result(49)=E(2,0,0,2,0) 				    
      E0result(50)=E(2,0,0,1,1) 				    
      E0result(51)=E(2,0,0,0,2) 				    
      E0result(52)=E(0,4,0,0,0) 				    
      E0result(53)=E(0,3,1,0,0) 				    
      E0result(54)=E(0,3,0,1,0) 				    
      E0result(55)=E(0,3,0,0,1) 				    
      E0result(56)=E(0,2,2,0,0) 				    
      E0result(57)=E(0,2,1,1,0) 				    
      E0result(58)=E(0,2,1,0,1) 				    
      E0result(59)=E(0,2,0,2,0) 				    
      E0result(60)=E(0,2,0,1,1) 				    
      E0result(61)=E(0,2,0,0,2) 				    
      E0result(62)=E(0,1,3,0,0) 				    
      E0result(63)=E(0,1,2,1,0) 				    
      E0result(64)=E(0,1,2,0,1) 				    
      E0result(65)=E(0,1,1,2,0) 				    
      E0result(66)=E(0,1,1,1,1) 				    
      E0result(67)=E(0,1,1,0,2) 				    
      E0result(68)=E(0,1,0,3,0) 				    
      E0result(69)=E(0,1,0,2,1) 				    
      E0result(70)=E(0,1,0,1,2) 				    
      E0result(71)=E(0,1,0,0,3) 				    
      E0result(72)=E(0,0,4,0,0) 				    
      E0result(73)=E(0,0,3,1,0) 				    
      E0result(74)=E(0,0,3,0,1) 				    
      E0result(75)=E(0,0,2,2,0) 				    
      E0result(76)=E(0,0,2,1,1) 				    
      E0result(77)=E(0,0,2,0,2) 				    
      E0result(78)=E(0,0,1,3,0) 				    
      E0result(79)=E(0,0,1,2,1) 				    
      E0result(80)=E(0,0,1,1,2) 				    
      E0result(81)=E(0,0,1,0,3) 				    
      E0result(82)=E(0,0,0,4,0) 				    
      E0result(83)=E(0,0,0,3,1) 				    
      E0result(84)=E(0,0,0,2,2) 				    
      E0result(85)=E(0,0,0,1,3) 				    
      E0result(86)=E(0,0,0,0,4) 				    
      If(rank.GE.5) Then					    
        Print *,"rank not implemented yet"			    
      Endif							    
      Endif							    
      Endif							    
      Endif							    
      Endif							    
      Return
#endif
#endif
#ifdef COLILib
#ifdef DDLib
      call cEp1234(mom1,mom2,mom3,mom4,mom5,mom6,mom7,mom8,mom9,mom10,  	           
     &  	    mass1,mass2,mass3,mass4,mass5,aEn0,aEn1,aEn2,aEn3,aEn4,rank)            
      call small_mass5(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,  			           
     &  	       mass1,mass2,mass3,mass4,mass5,				           
     &  	       pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10,		           
     &  	       mm1,mm2,mm3,mm4,mm5)					           
      call E_dd(E,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10,			           
     &  	    mm1,mm2,mm3,mm4,mm5,rank,0) 				          
      call check_LI(aEn0,E(0,0,0,0,0),116)					           
      E0result(1)=E(0,0,0,0,0)  						           
      If(rank.GE.1) Then							           
      call check_LI(aEn1(1),E(0,1,0,0,0),117)					           
      E0result(2)=E(0,1,0,0,0)  						           
      call check_LI(aEn1(2),E(0,0,1,0,0),118)					           
      E0result(3)=E(0,0,1,0,0)  						           
      call check_LI(aEn1(3),E(0,0,0,1,0),119)					           
      E0result(4)=E(0,0,0,1,0)  						           
      call check_LI(aEn1(4),E(0,0,0,0,1),120)					           
      E0result(5)=E(0,0,0,0,1)  						           
      If(rank.GE.2) Then							           
      call check_LI(aEn2(0,0),E(2,0,0,0,0),121) 				           
      E0result(6)=E(2,0,0,0,0)  						           
      call check_LI(aEn2(1,1),E(0,2,0,0,0),122) 				           
      E0result(7)=E(0,2,0,0,0)  						           
      call check_LI(aEn2(1,2),E(0,1,1,0,0),123) 				           
      E0result(8)=E(0,1,1,0,0)  						           
      call check_LI(aEn2(1,3),E(0,1,0,1,0),124) 				           
      E0result(9)=E(0,1,0,1,0)  						           
      call check_LI(aEn2(1,4),E(0,1,0,0,1),125) 				           
      E0result(10)=E(0,1,0,0,1) 						           
      call check_LI(aEn2(2,2),E(0,0,2,0,0),126) 				           
      E0result(11)=E(0,0,2,0,0) 						           
      call check_LI(aEn2(2,3),E(0,0,1,1,0),127) 				           
      E0result(12)=E(0,0,1,1,0) 						           
      call check_LI(aEn2(2,4),E(0,0,1,0,1),128) 				           
      E0result(13)=E(0,0,1,0,1) 						           
      call check_LI(aEn2(3,3),E(0,0,0,2,0),129) 				           
      E0result(14)=E(0,0,0,2,0) 						           
      call check_LI(aEn2(3,4),E(0,0,0,1,1),130) 				           
      E0result(15)=E(0,0,0,1,1) 						           
      call check_LI(aEn2(4,4),E(0,0,0,0,2),131) 				           
      E0result(16)=E(0,0,0,0,2) 						           
      If(rank.GE.3) Then							           
      call check_LI(aEn3(0,0,1),E(2,1,0,0,0),132)				           
      E0result(17)=E(2,1,0,0,0) 						           
      call check_LI(aEn3(0,0,2),E(2,0,1,0,0),133)				           
      E0result(18)=E(2,0,1,0,0) 						           
      call check_LI(aEn3(0,0,3),E(2,0,0,1,0),134)				           
      E0result(19)=E(2,0,0,1,0) 						           
      call check_LI(aEn3(0,0,4),E(2,0,0,0,1),135)				           
      E0result(20)=E(2,0,0,0,1) 						           
      call check_LI(aEn3(1,1,1),E(0,3,0,0,0),136)				           
      E0result(21)=E(0,3,0,0,0) 						           
      call check_LI(aEn3(1,1,2),E(0,2,1,0,0),137)				           
      E0result(22)=E(0,2,1,0,0) 						           
      call check_LI(aEn3(1,1,3),E(0,2,0,1,0),138)				           
      E0result(23)=E(0,2,0,1,0) 						           
      call check_LI(aEn3(1,1,4),E(0,2,0,0,1),139)				           
      E0result(24)=E(0,2,0,0,1) 						           
      call check_LI(aEn3(1,2,2),E(0,1,2,0,0),140)				           
      E0result(25)=E(0,1,2,0,0) 						           
      call check_LI(aEn3(1,2,3),E(0,1,1,1,0),141)				           
      E0result(26)=E(0,1,1,1,0) 						           
      call check_LI(aEn3(1,2,4),E(0,1,1,0,1),142)				           
      E0result(27)=E(0,1,1,0,1) 						           
      call check_LI(aEn3(1,3,3),E(0,1,0,2,0),143)				           
      E0result(28)=E(0,1,0,2,0) 						           
      call check_LI(aEn3(1,3,4),E(0,1,0,1,1),144)				           
      E0result(29)=E(0,1,0,1,1) 						           
      call check_LI(aEn3(1,4,4),E(0,1,0,0,2),145)				           
      E0result(30)=E(0,1,0,0,2) 						           
      call check_LI(aEn3(2,2,2),E(0,0,3,0,0),146)				           
      E0result(31)=E(0,0,3,0,0) 						           
      call check_LI(aEn3(2,2,3),E(0,0,2,1,0),147)				           
      E0result(32)=E(0,0,2,1,0) 						           
      call check_LI(aEn3(2,2,4),E(0,0,2,0,1),148)				           
      E0result(33)=E(0,0,2,0,1) 						           
      call check_LI(aEn3(2,3,3),E(0,0,1,2,0),149)				           
      E0result(34)=E(0,0,1,2,0) 						           
      call check_LI(aEn3(2,3,4),E(0,0,1,1,1),150)				           
      E0result(35)=E(0,0,1,1,1) 						           
      call check_LI(aEn3(2,4,4),E(0,0,1,0,2),151)				           
      E0result(36)=E(0,0,1,0,2) 						           
      call check_LI(aEn3(3,3,3),E(0,0,0,3,0),152)				           
      E0result(37)=E(0,0,0,3,0) 						           
      call check_LI(aEn3(3,3,4),E(0,0,0,2,1),153)				           
      E0result(38)=E(0,0,0,2,1) 						           
      call check_LI(aEn3(3,4,4),E(0,0,0,1,2),154)				           
      E0result(39)=E(0,0,0,1,2) 						           
      call check_LI(aEn3(4,4,4),E(0,0,0,0,3),155)				           
      E0result(40)=E(0,0,0,0,3) 						           
      If(rank.GE.4) Then							           
      call check_LI(aEn4(0,0,0,0),E(4,0,0,0,0),156)				           
      E0result(41)=E(4,0,0,0,0) 						           
      call check_LI(aEn4(0,0,1,1),E(2,2,0,0,0),157)				           
      E0result(42)=E(2,2,0,0,0) 						           
      call check_LI(aEn4(0,0,1,2),E(2,1,1,0,0),158)				           
      E0result(43)=E(2,1,1,0,0) 						           
      call check_LI(aEn4(0,0,1,3),E(2,1,0,1,0),159)				           
      E0result(44)=E(2,1,0,1,0) 						           
      call check_LI(aEn4(0,0,1,4),E(2,1,0,0,1),10)				           
      E0result(45)=E(2,1,0,0,1) 						           
      call check_LI(aEn4(0,0,2,2),E(2,0,2,0,0),161)				           
      E0result(46)=E(2,0,2,0,0) 						           
      call check_LI(aEn4(0,0,2,3),E(2,0,1,1,0),162)				           
      E0result(47)=E(2,0,1,1,0) 						           
      call check_LI(aEn4(0,0,2,4),E(2,0,1,0,1),163)				           
      E0result(48)=E(2,0,1,0,1) 						           
      call check_LI(aEn4(0,0,3,3),E(2,0,0,2,0),164)				           
      E0result(49)=E(2,0,0,2,0) 						           
      call check_LI(aEn4(0,0,3,4),E(2,0,0,1,1),165)				           
      E0result(50)=E(2,0,0,1,1) 						           
      call check_LI(aEn4(0,0,4,4),E(2,0,0,0,2),166)				           
      E0result(51)=E(2,0,0,0,2) 						           
      call check_LI(aEn4(1,1,1,1),E(0,4,0,0,0),167)				           
      E0result(52)=E(0,4,0,0,0) 						           
      call check_LI(aEn4(1,1,1,2),E(0,3,1,0,0),168)				           
      E0result(53)=E(0,3,1,0,0) 						           
      call check_LI(aEn4(1,1,1,3),E(0,3,0,1,0),169)				           
      E0result(54)=E(0,3,0,1,0) 						           
      call check_LI(aEn4(1,1,1,4),E(0,3,0,0,1),170)				           
      E0result(55)=E(0,3,0,0,1) 						           
      call check_LI(aEn4(1,1,2,2),E(0,2,2,0,0),171)				           
      E0result(56)=E(0,2,2,0,0) 						           
      call check_LI(aEn4(1,1,2,3),E(0,2,1,1,0),172)				           
      E0result(57)=E(0,2,1,1,0) 						           
      call check_LI(aEn4(1,1,2,4),E(0,2,1,0,1),173)				           
      E0result(58)=E(0,2,1,0,1) 						           
      call check_LI(aEn4(1,1,3,3),E(0,2,0,2,0),174)				           
      E0result(59)=E(0,2,0,2,0) 						           
      call check_LI(aEn4(1,1,3,4),E(0,2,0,1,1),175)				           
      E0result(60)=E(0,2,0,1,1) 						           
      call check_LI(aEn4(1,1,4,4),E(0,2,0,0,2),176)				           
      E0result(61)=E(0,2,0,0,2) 						           
      call check_LI(aEn4(1,2,2,2),E(0,1,3,0,0),177)				           
      E0result(62)=E(0,1,3,0,0) 						           
      call check_LI(aEn4(1,2,2,3),E(0,1,2,1,0),178)				           
      E0result(63)=E(0,1,2,1,0) 						           
      call check_LI(aEn4(1,2,2,4),E(0,1,2,0,1),179)				           
      E0result(64)=E(0,1,2,0,1) 						           
      call check_LI(aEn4(1,2,3,3),E(0,1,1,2,0),180)				           
      E0result(65)=E(0,1,1,2,0) 						           
      call check_LI(aEn4(1,2,3,4),E(0,1,1,1,1),181)				           
      E0result(66)=E(0,1,1,1,1) 						           
      call check_LI(aEn4(1,2,4,4),E(0,1,1,0,2),182)				           
      E0result(67)=E(0,1,1,0,2) 						           
      call check_LI(aEn4(1,3,3,3),E(0,1,0,3,0),183)				           
      E0result(68)=E(0,1,0,3,0) 						           
      call check_LI(aEn4(1,3,3,4),E(0,1,0,2,1),184)				           
      E0result(69)=E(0,1,0,2,1) 						           
      call check_LI(aEn4(1,3,4,4),E(0,1,0,1,2),185)				           
      E0result(70)=E(0,1,0,1,2) 						           
      call check_LI(aEn4(1,4,4,4),E(0,1,0,0,3),186)				           
      E0result(71)=E(0,1,0,0,3) 						           
      call check_LI(aEn4(2,2,2,2),E(0,0,4,0,0),187)				           
      E0result(72)=E(0,0,4,0,0) 						           
      call check_LI(aEn4(2,2,2,3),E(0,0,3,1,0),188)				           
      E0result(73)=E(0,0,3,1,0) 						           
      call check_LI(aEn4(2,2,2,4),E(0,0,3,0,1),189)				           
      E0result(74)=E(0,0,3,0,1) 						           
      call check_LI(aEn4(2,2,3,3),E(0,0,2,2,0),190)				           
      E0result(75)=E(0,0,2,2,0) 						           
      call check_LI(aEn4(2,2,3,4),E(0,0,2,1,1),191)				           
      E0result(76)=E(0,0,2,1,1) 						           
      call check_LI(aEn4(2,2,4,4),E(0,0,2,0,2),192)				           
      E0result(77)=E(0,0,2,0,2) 						           
      call check_LI(aEn4(2,3,3,3),E(0,0,1,3,0),193)				           
      E0result(78)=E(0,0,1,3,0) 						           
      call check_LI(aEn4(2,3,3,4),E(0,0,1,2,1),194)				           
      E0result(79)=E(0,0,1,2,1) 						           
      call check_LI(aEn4(2,3,4,4),E(0,0,1,1,2),195)				           
      E0result(80)=E(0,0,1,1,2) 						           
      call check_LI(aEn4(2,4,4,4),E(0,0,1,0,3),196)				           
      E0result(81)=E(0,0,1,0,3) 						           
      call check_LI(aEn4(3,3,3,3),E(0,0,0,4,0),197)				           
      E0result(82)=E(0,0,0,4,0) 						           
      call check_LI(aEn4(3,3,3,4),E(0,0,0,3,1),198)				           
      E0result(83)=E(0,0,0,3,1) 						           
      call check_LI(aEn4(3,3,4,4),E(0,0,0,2,2),199)				           
      E0result(84)=E(0,0,0,2,2) 						           
      call check_LI(aEn4(3,4,4,4),E(0,0,0,1,3),200)				           
      E0result(85)=E(0,0,0,1,3) 						           
      call check_LI(aEn4(4,4,4,4),E(0,0,0,0,4),201)				           
      E0result(86)=E(0,0,0,0,4) 						           
      If(rank.GE.5) Then							           
        Print *,"rank not implemented yet"					           
      Endif									           
      Endif									           
      Endif									           
      Endif									           
      Endif
      Return
#endif
#endif
#endif

      myEget=0
      Print *,"Hallo 13"
      Print *,"no loop library inluded"
      Stop
           
      End	     

**************************************************************************************

      double complex function myEval(i1,i2)
      implicit none 
      
      integer i1,i2
      double complex E0result(86)
      
#include "looptools.h" 

      common/E0_result_common/E0result
            
#ifdef LTLib
        myEval=Eval(i1,i2)
#else
        myEval=E0result(i1)
#endif
      
      End

**************************************************************************************

      subroutine check_LI(LI1,LI2,id)
      Implicit none 
      
      double complex LI1,LI2
      double precision acc,rLI1,iLI1,rLI2,iLI2
      integer id,toprint
      
      rLI1=dreal(LI1)
      iLI1=dimag(LI1)
      rLI2=dreal(LI2)
      iLI2=dimag(LI2)
      
      If(id.EQ.7 .OR. id.EQ.6) Then
        acc=1D-1
      Else
        acc=1D-1
      Endif
      toprint=0
      If(abs(rLI1-rLI2)/abs(rLI1).GT.acc) Then
        toprint=1
      Endif
      If(iLI1/rLI1.GT.acc .OR. iLI2/rLI2.GT.acc) Then
        If(abs(iLI1-iLI2)/abs(iLI1).GT.acc) Then
          toprint=1
        Endif    
      Endif	
      
      If(toprint.EQ.1) Then
        Print *,"DD and coli differ id=",id
	Print *,LI1
	Print *,LI2
      Endif  
      
      End
