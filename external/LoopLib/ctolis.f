***********************************************************************
*                                                                     *
*     reduction of 3-, 4- and 5-point tensor functions                *
*     reduction formula according to                                  *
*     A. Denner, S. Dittmaier, Nucl. Phys. B658 (2003) 175-202        *
*     version that calculates  only coefficientes with                *
*     ascending indices                                               *
*                                                                     *
***********************************************************************
*                                                                     *
*     last changed  12.12.05  Ansgar Denner                           *
*                                                                     *
***********************************************************************
* subroutines:                                                        *
* cBp123,cCpv12345,cDpv12345,cEp1234,cFp1234,chinv                    *
* cE0f,cF0f,chdet                                                     *
***********************************************************************
      subroutine cBp123(p10,m02,m12,B0,B1,B2,B3,rank)
***********************************************************************
*     2-point tensor coefficient functions B0,B1,B2,B3                *
*     up to 3 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     definition of integrals:                                        *
*     {B0,B1,B2,B3}= 1/(i pi^2)*\int d^4 q {1,q,q.q,q.q.q} *          *
*        1/[(q**2 - m02)*[(q+k1)**2-m12] ]                            *
*     arguments:   p10=k1**2                                          *
*     definition of tensor coefficients:                              *
*          B1 = B1(i) *k1                                             *
*          B2 = B2(1,1) *k1*k1  + B2(0,0) * g                         *
*          B3 = B3(1,1,1) *k1*k1*k1 + B3(0,0,i) * (k1*g + 2*sym)      *
*---------------------------------------------------------------------*
*     09.12.02  Ansgar Denner     last changed  07.06.04              *
***********************************************************************
      implicit   none
      complex*16 p10
      complex*16 m02,m12
c     complex*16 q10, mm0, mm1
c     complex*16 f1,mm02,mm12
c     complex*16 A0w1,A0w0
      complex*16 B0,B1,B2(0:1,0:1),B3(0:1,0:1,0:1)
      complex*16 cB0f,cB1f,cB00f,cB11f,cB001f,cB111f
      real*8     mudim2
c     real*8     elimcminf2
      integer    rank
      integer    ltest

      common /uv/ mudim2

      common /ltest/ ltest
 
c      write(*,*) 'cBp123 in',p10,m02,m12,rank

      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cBp123'
        write(*,*) 'rank = ',rank
        stop
      end if

      B0   = cB0f(p10,m02,m12)

      if (rank.eq.0) return

CC      f1 = q10+mm02-mm12
CC
CC      A0w1 = mm02*(1d0-log(m02/mudim2))
CC      A0w0 = mm12*(1d0-log(m12/mudim2))
CC      
CC      B1 = (A0w1-A0w0-f1*B0)/(2d0*q10)
      B1   = cB1f(p10,m02,m12)

      if (rank.eq.1) return

CC      B2(0,0) = (A0w0+f1*B1+2*mm02*B0+mm02+mm12-q10/3d0)/6d0
CC      B2(1,1) = (A0w0-2d0*f1*B1-mm02*B0-(m02+mm12-q10/3d0)/2d0)/(3d0*q10)
      B2(0,0) = cB00f(p10,m02,m12)
      B2(1,1) = cB11f(p10,m02,m12)

      if (rank.eq.2) return

c      B3(0,0,1) = (-A0w0+f1*B2(1,1)+2*mm02*B1
c     &            -(2d0*mm02+4d0*mm12-q10)/6d0)/8d0
c      B3(1,1,1) = -(A0w0+3d0*f1*B2(1,1)+2d0*mm02*B1
c     &            -(2d0*mm02+4d0*mm12-q10)/6d0)/(4d0*q10)

      B3(0,0,1) = cB001f(p10,m02,m12)
      B3(1,1,1) = cB111f(p10,m02,m12)
     
      if (rank.gt.3) then
        write(*,*) 'rank > 3 not implemented in cBp123'
        write(*,*) 'rank = ',rank
      end if

      end
***********************************************************************
      subroutine cCpv12345(p10,p21,p20,m02,m12,m22,
     &                C0,C1,C2,C3,C4,C5,rank)
c      subroutine cCpv1234(p10,p21,p20,m02,m12,m22,
c     &                C0,C1,C2,C3,C4,rank)
c      entry Cp1234B(p10,p21,p20,m02,m12,m22,
c     &                C0,C1,C2,C3,C4,B0w,B1w,B2w,B3w,rank)
***********************************************************************
*     3-point tensor coefficient functions C0,C1,C2,C3,C4             *
*     up to 4 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     definition of integrals:                                        *
*     {C0,C1,C2,C3,C4}= 1/(i pi^2)*\int d^4 q {1,q,q.q,q.q.q,q.q.q.q} *
*        1/[(q**2 - m02)*[(q+k1)**2-m12]*[(q+k2)**2 - m22] ]          *
*     arguments:   p10=k1**2, p21=(k2-k1)**2, p20=k2**2               *
*     definition of tensor coefficients:                              *
*          C1 = \sum_{i=1}^2 C1(i) *ki                                *
*          C2 = \sum_{i,j=1}^2 C2(i,j) *ki*kj  + C2(0,0) * g          *
*          C3 = \sum_{i,j,k=1}^2 C3(i,j,k) *ki*kj*kk                  *
*              + \sum_{i=1}^2 C3(0,0,i) * (ki*g + 2*sym)              *
*          C4 = \sum_{i,j,k,l=1}^2 C4(i,j,k) *ki*kj*kk*kl             *
*              + \sum_{i,j=1}^2 C4(0,0,i,j) * (ki*kj*g + 5*sym)       *
*              +   C4(0,0,0,0) * (g*g + 2*sym)                        *
*          B1w0 =  B1w0(i) *( k2 - k1 )                               *
*     2-point functions resulting from 3-point functions with         *
*     i-th denominator cancelled: BxwI                                *
*---------------------------------------------------------------------*
*     06.12.02  Ansgar Denner     last changed  24.01.05              *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p20
      complex*16 q10,q21,q20
      complex*16 m02,m12,m22
      complex*16 mm02,mm12,mm22
c      real*8     k1k2,mm0,mm1,mm2
      real*8     k1k2
      complex*16 f(2)
      complex*16 zinv(2,2),detz
      complex*16 B0w2,B1w2,B2w2(0:1,0:1),B3w2(0:1,0:1,0:1)
      complex*16 B4w2(0:1,0:1,0:1,0:1),B5w2(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w2(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0w1,B1w1,B2w1(0:1,0:1),B3w1(0:1,0:1,0:1)
      complex*16 B4w1(0:1,0:1,0:1,0:1),B5w1(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w1(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0w0,B1w0,B2w0(0:1,0:1),B3w0(0:1,0:1,0:1)
      complex*16 B4w0(0:1,0:1,0:1,0:1),B5w0(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w0(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0m(0:2),B1m(0:2,2),B2m(0:2,0:2,0:2),
     &           B3m(0:2,0:2,0:2,2),B4m(0:2,0:2,0:2,0:2,0:2)
c     &           ,B5m(0:2,0:2,0:2,0:2,0:2,2)
      complex*16 S1(2),S2(2,2),S3(2,0:2,0:2),S4(2,0:2,0:2,2)
      complex*16 S5(2,0:2,0:2,0:2,0:2)
      complex*16 C0,C1(2),C2(0:2,0:2),C3(0:2,0:2,0:2)
      complex*16 C4(0:2,0:2,0:2,0:2),C5(0:2,0:2,0:2,0:2,0:2)
c      complex*16 C6(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 cC0f 
      complex*16 elimcminf2
      integer    rank,rankB
      integer    i1,i2,i3,i4,k,l
      integer    ltest,sym

      common /ltest/  ltest
      common /sym/    sym
                                                                               
      data       B1m /6*0d0/, B2m /27*0d0/, B3m/54*0d0/, B4m /243*0d0/
c      data       B5m /486*0d0/

c      write(*,*) 'Cp1234 in',p10,p21,p20,m02,m12,m22,rank
c      write(*,*) 'Cp1234 in',C0,C1,C2
c      write(*,*) 'Cp1234 in',C3
c      write(*,*) 'Cp1234 in',C4

      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cCpv12345'
        write(*,*) 'rank = ',rank
        stop
      end if

      C0   = cC0f(p10,p21,p20,m02,m12,m22)
c      write(*,*) 'C0 = ',C0

c  cBp must be called also for rank = 0 because of cache!

      rankB = max(rank-1,0)
c      write(*,*) 'Cp0'
      call cBp12345(p21,m12,m22,B0w0,B1w0,B2w0,B3w0,B4w0,B5w0,B6w0,
     &    rankB,0)
c      write(*,*) 'Cp1'
      call cBp12345(p20,m02,m22,B0w1,B1w1,B2w1,B3w1,B4w1,B5w1,B6w1,
     &    rankB,0)
c      write(*,*) 'Cp2'
      call cBp12345(p10,m02,m12,B0w2,B1w2,B2w2,B3w2,B4w2,B5w2,B6w2,
     &    rankB,0)

      if (rank.eq.0) return

      mm02 = elimcminf2(m02)
      mm12 = elimcminf2(m12)
      mm22 = elimcminf2(m22)
      q10  = elimcminf2(p10)
      q21  = elimcminf2(p21)
      q20  = elimcminf2(p20)
 
      k1k2 = (q10+q20-q21)/2D0
      detz = q10*q20-k1k2*k1k2
      zinv(1,1) = q20/(2d0*detz)
      zinv(1,2) = -k1k2/(2d0*detz)
      zinv(2,1) = zinv(1,2)
      zinv(2,2) = q10/(2d0*detz)
      f(1) = q10+mm02-mm12
      f(2) = q20+mm02-mm22
 
      B0m(0) = B0w0
      B0m(1) = B0w1
      B0m(2) = B0w2

      do 11 k =1,2
        S1(k) = B0m(k) - B0m(0) - f(k)*C0
11    continue
      do 16 k =1,2
        C1(k) = 0D0
      do 16 l =1,2
        C1(k) = C1(k) + zinv(k,l)*S1(l)
16    continue

      if (rank.eq.1) return

      C2(0,0) = mm02*C0/2d0+ (B0w0 +1d0 +f(1)*C1(1)+f(2)*C1(2))/4d0

      B1m(0,2) = B1w0
      B1m(0,1) = -B0m(0) - B1m(0,2)
      B1m(1,2) = B1w1
      B1m(2,1) = B1w2

      do 26 i1=1,2
      do 21 k =1,2
        S2(k,i1) = B1m(k,i1) - B1m(0,i1) - f(k)*C1(i1)
        if(k.eq.i1)  S2(k,i1) = S2(k,i1) - 2d0*C2(0,0)
21    continue
      do 26 k =1,i1
        C2(k,i1) = 0D0
      do 26 l =1,2
        C2(k,i1) = C2(k,i1) + zinv(k,l)*S2(l,i1)
26    continue
      if (sym.eq.1) then
        C2(2,1) = C2(1,2)
      end if

      if (rank.eq.2) return
 
      B2m(0,0,0) = B2w0(0,0)
      B2m(1,0,0) = B2w1(0,0)
      B2m(2,0,0) = B2w2(0,0)
      B2m(0,2,2) = B2w0(1,1)
      B2m(0,1,2) = -B1m(0,2) - B2m(0,2,2)
      B2m(0,1,1) = -B1m(0,1) - B2m(0,1,2)
      B2m(1,2,2) = B2w1(1,1)
      B2m(2,1,1) = B2w2(1,1)
 
c     do 35 i1 =1,2
c        C3(0,0,i1) = -1d0/18d0+mm02*C1(i1)/3d0 +(B1m(0,i1)+
c    &     +f(1)*C2(1,i1)+f(2)*C2(i1,2))/6d0      
c35   continue

      do 30 k =1,2
        S3(k,0,0) = B2m(k,0,0) - B2m(0,0,0) - f(k)*C2(0,0)
 30   continue

      do 34 k =1,2
        C3(0,0,k) = 0D0
      do 34 l =1,2
        C3(0 ,0,k) = C3(0 ,0,k) + zinv(k,l)*S3(l,0 ,0 )
 34   continue

      do 38 i1=1,2
      do 38 i2=i1,2
      do 31 k =1,2
        S3(k,i1,i2) = B2m(k,i1,i2) - B2m(0,i1,i2) - f(k)*C2(i1,i2)
        if(k.eq.i1)  S3(k,i1,i2) = S3(k,i1,i2) - 2d0*C3(0,0,i2)
        if(k.eq.i2)  S3(k,i1,i2) = S3(k,i1,i2) - 2d0*C3(0,0,i1)
31    continue
      do 36 k =1,i1
        C3(k,i1,i2) = 0D0
      do 36 l =1,2
        C3(k,i1,i2) = C3(k,i1,i2) + zinv(k,l)*S3(l,i1,i2)
36    continue
 38   continue
      if (sym.eq.1) then
        C3(2,2,1) = C3(1,2,2)
        C3(2,1,2) = C3(1,2,2)
        C3(2,1,1) = C3(1,1,2)
        C3(1,2,1) = C3(1,1,2)
      end if

      if (rank.eq.3) return
 
      C4(0,0,0,0) = (m02+mm12+mm22)/48d0 -(q10+q20+q21)/192d0
     &   + mm02*C2(0,0)/4d0 
     &   + (B2w0(0,0)+f(1)*C3(0,0,1)+f(2)*C3(0,0,2))/8d0
      
      B3m(0,0,0,1) = -B2w0(0,0)-B3w0(0,0,1)
      B3m(0,0,0,2) = B3w0(0,0,1)
      B3m(1,0,0,2) = B3w1(0,0,1)
      B3m(2,0,0,1) = B3w2(0,0,1)
      B3m(0,2,2,2) = B3w0(1,1,1)

      B3m(0,1,2,2) = -B2m(0,2,2) - B3m(0,2,2,2)
      B3m(0,1,1,2) = -B2m(0,1,2) - B3m(0,1,2,2)
      B3m(0,1,1,1) = -B2m(0,1,1) - B3m(0,1,1,2)

      B3m(1,2,2,2) = B3w1(1,1,1)
      B3m(2,1,1,1) = B3w2(1,1,1)

c     do 43 i1 =1,2
c       C4(0,0,i1,i1) = 1d0/48d0
c       C4(0,0,i1,3-i1) = 1d0/96d0
c     do 43 i2 =i1,2
c       C4(0,0,i1,i2) =  C4(0,0,i1,i2)+ mm02*C2(i1,i2)/4d0 + 
c    &     (B2m(0,i1,i2) +f(1)*C3(1,i1,i2)+f(2)*C3(i1,i2,2))/8d0      
c43   continue
c     C4(0,0,2,1) = C4(0,0,1,2)

      do 45 i1=1,2
      do 41 k =1,2
        S4(k,0,0,i1) = B3m(k,0,0,i1) - B3m(0,0,0,i1)
     &                   - f(k)*C3(0,0,i1)
        if(k.eq.i1)  S4(k,0,0,i1) = S4(k,0,0,i1) - 2d0*C4(0,0,0,0)
41    continue
      do 45 k =1,i1
         C4(0,0,k,i1) = 0D0
      do 45 l =1,2
         C4(0,0,k,i1) = C4(0,0,k,i1) + zinv(k,l)*S4(l,0,0,i1)
 45   continue
      if (sym.eq.1) then
        C4(0,0,2,1) = C4(0,0,1,2)
      end if

      do 49 i1=1,2
      do 49 i2=i1,2
      do 49 i3=i2,2
      do 47 k =1,2
        S4(k,i1,i2,i3) = B3m(k,i1,i2,i3) - B3m(0,i1,i2,i3)
     &                   - f(k)*C3(i1,i2,i3)
        if(k.eq.i1) S4(k,i1,i2,i3) = S4(k,i1,i2,i3) - 2d0*C4(0,0,i2,i3)
        if(k.eq.i2) S4(k,i1,i2,i3) = S4(k,i1,i2,i3) - 2d0*C4(0,0,i1,i3)
        if(k.eq.i3) S4(k,i1,i2,i3) = S4(k,i1,i2,i3) - 2d0*C4(0,0,i1,i2)
47    continue
      do 49 k =1,i1
         C4(k,i1,i2,i3) = 0D0
      do 49 l =1,2
         C4(k,i1,i2,i3) = C4(k,i1,i2,i3) + zinv(k,l)*S4(l,i1,i2,i3)
49    continue
      if (sym.eq.1) then
        C4(1,1,2,1) = C4(1,1,1,2)
        C4(1,2,1,1) = C4(1,1,1,2)
        C4(2,1,1,1) = C4(1,1,1,2)
        C4(1,2,2,1) = C4(1,1,2,2)
        C4(2,1,2,1) = C4(1,1,2,2)
        C4(2,2,1,1) = C4(1,1,2,2)
        C4(1,2,1,2) = C4(1,1,2,2)
        C4(2,1,1,2) = C4(1,1,2,2)
        C4(2,2,2,1) = C4(1,2,2,2)
        C4(2,2,1,2) = C4(1,2,2,2)
        C4(2,1,2,2) = C4(1,2,2,2)
      end if

      if (rank.eq.4) return
 
      do i1 = 1,2
      C5(0,0,0,0,i1) = -(mm02+mm12+mm22)/240d0 +(q10+q20+2d0*q21)/1200d0
     &   + mm02*C3(0,0,i1)/5d0 
     &   + (B3m(0,0,0,i1)+f(1)*C4(0,0,1,i1)+f(2)*C4(0,0,i1,2))/10d0
      end do      
      C5(0,0,0,0,1) = C5(0,0,0,0,1) - (5d0*mm12-q10)/1200d0
      C5(0,0,0,0,2) = C5(0,0,0,0,2) - (5d0*mm22-q20)/1200d0

      B4m(0,0,0,0,0) = B4w0(0,0,0,0)
      B4m(1,0,0,0,0) = B4w1(0,0,0,0)
      B4m(2,0,0,0,0) = B4w2(0,0,0,0)

      B4m(0,0,0,2,2) = B4w0(0,0,1,1)
      B4m(0,0,0,1,2) = -B3m(0,0,0,2)-B4m(0,0,0,2,2)
      B4m(0,0,0,1,1) = -B3m(0,0,0,1)-B4m(0,0,0,1,2)

      B4m(0,2,2,2,2) = B4w0(1,1,1,1)
      B4m(0,1,2,2,2) = -B3m(0,2,2,2)-B4m(0,2,2,2,2)
      B4m(0,1,1,2,2) = -B3m(0,1,2,2)-B4m(0,1,2,2,2)
      B4m(0,1,1,1,2) = -B3m(0,1,1,2)-B4m(0,1,1,2,2)
      B4m(0,1,1,1,1) = -B3m(0,1,1,1)-B4m(0,1,1,1,2)

      B4m(1,0,0,2,2) = B4w1(0,0,1,1)
      B4m(2,0,0,1,1) = B4w2(0,0,1,1)      
      
      B4m(1,2,2,2,2) = B4w1(1,1,1,1)
      B4m(2,1,1,1,1) = B4w2(1,1,1,1)
      
      B4m(0,0,0,2,1) = B4m(0,0,0,1,2) 

      B4m(0,2,1,1,1) = B4m(0,1,1,1,2) 
cs      B4m(0,1,2,1,1) = B4m(0,1,1,1,2) 
cs      B4m(0,1,1,2,1) = B4m(0,1,1,1,2) 
cs      B4m(0,1,2,1,2) = B4m(0,1,1,2,2) 
cs      B4m(0,1,2,2,1) = B4m(0,1,1,2,2) 
      B4m(0,2,1,1,2) = B4m(0,1,1,2,2) 
cs      B4m(0,2,1,2,1) = B4m(0,1,1,2,2) 
cs      B4m(0,2,2,1,1) = B4m(0,1,1,2,2) 
      B4m(0,2,1,2,2) = B4m(0,1,2,2,2) 
cs      B4m(0,2,2,1,2) = B4m(0,1,2,2,2) 
cs      B4m(0,2,2,2,1) = B4m(0,1,2,2,2) 


      do 53 i1 =1,2
        C5(0,0,i1,i1,i1) = -1d0/100d0
        C5(0,0,i1,i1,3-i1) = -1d0/300d0
      do 53 i2 =i1,2
      do 53 i3 =i1,2
        C5(0,0,i1,i2,i3) =  C5(0,0,i1,i2,i3)+ mm02*C3(i1,i2,i3)/5d0 + 
     &      (B3m(0,i1,i2,i3) +f(1)*C4(1,i1,i2,i3)+f(2)*C4(i1,i2,i3,2))
     &      /10d0      
 53   continue
      if (sym.eq.1) then
        C5(0,0,1,2,1) = C5(0,0,1,1,2)
        C5(0,0,2,1,1) = C5(0,0,1,1,2)
        C5(0,0,2,2,1) = C5(0,0,1,2,2)
        C5(0,0,2,1,2) = C5(0,0,1,2,2)
      end if

      do 55 i1=1,2
      do 55 i2=1,2
      do 51 k =1,2
        S5(k,0,0,i1,i2) = B4m(k,0,0,i1,i2) - B4m(0,0,0,i1,i2)
     &                   - f(k)*C4(0,0,i1,i2)
        if(k.eq.i1) S5(k,0,0,i1,i2) = S5(k,0,0,i1,i2)-2d0*C5(0,0,0,0,i2)
        if(k.eq.i2) S5(k,0,0,i1,i2) = S5(k,0,0,i1,i2)-2d0*C5(0,0,0,0,i1)
51    continue
      do 55 k =1,i1
         C5(0,0,k,i1,i2) = 0D0
      do 55 l =1,2
         C5(0,0,k,i1,i2) = C5(0,0,k,i1,i2) + zinv(k,l)*S5(l,0,0,i1,i2)
 55   continue
      if (sym.eq.1) then
        C5(0,0,1,2,1) = C5(0,0,1,1,2)
        C5(0,0,2,1,1) = C5(0,0,1,1,2)
        C5(0,0,2,2,1) = C5(0,0,1,2,2)
        C5(0,0,2,1,2) = C5(0,0,1,2,2)
      end if

      do 59 i1=1,2
      do 59 i2=i1,2
      do 59 i3=i2,2
      do 59 i4=i3,2
      do 57 k =1,2
        S5(k,i1,i2,i3,i4) = B4m(k,i1,i2,i3,i4) - B4m(0,i1,i2,i3,i4)
     &                   - f(k)*C4(i1,i2,i3,i4)
        if(k.eq.i1)  S5(k,i1,i2,i3,i4) = S5(k,i1,i2,i3,i4) 
     &      - 2d0*C5(0,0,i2,i3,i4)
        if(k.eq.i2)  S5(k,i1,i2,i3,i4) = S5(k,i1,i2,i3,i4)
     &      - 2d0*C5(0,0,i1,i3,i4)
        if(k.eq.i3)  S5(k,i1,i2,i3,i4) = S5(k,i1,i2,i3,i4) 
     &      - 2d0*C5(0,0,i1,i2,i4)
        if(k.eq.i4)  S5(k,i1,i2,i3,i4) = S5(k,i1,i2,i3,i4) 
     &      - 2d0*C5(0,0,i1,i2,i3)
57    continue
      do 59 k =1,i1
        C5(k,i1,i2,i3,i4) = 0D0
      do 59 l =1,2
        C5(k,i1,i2,i3,i4) = C5(k,i1,i2,i3,i4) 
     &      + zinv(k,l)*S5(l,i1,i2,i3,i4)
59    continue
      if (sym.eq.1) then
        C5(1,1,1,2,1) = C5(1,1,1,1,2)
        C5(1,1,2,1,1) = C5(1,1,1,1,2)
        C5(1,2,1,1,1) = C5(1,1,1,1,2)
        C5(2,1,1,1,1) = C5(1,1,1,1,2)

        C5(1,1,2,1,2) = C5(1,1,1,2,2)
        C5(1,2,1,1,2) = C5(1,1,1,2,2)
        C5(2,1,1,1,2) = C5(1,1,1,2,2)
        C5(1,1,2,2,1) = C5(1,1,1,2,2)
        C5(1,2,1,2,1) = C5(1,1,1,2,2)
        C5(2,1,1,2,1) = C5(1,1,1,2,2)
        C5(1,2,2,1,1) = C5(1,1,1,2,2)
        C5(2,1,2,1,1) = C5(1,1,1,2,2)
        C5(2,2,1,1,1) = C5(1,1,1,2,2)

        C5(1,2,1,2,2) = C5(1,1,2,2,2)
        C5(2,1,1,2,2) = C5(1,1,2,2,2)
        C5(1,2,2,1,2) = C5(1,1,2,2,2)
        C5(2,1,2,1,2) = C5(1,1,2,2,2)
        C5(2,2,1,1,2) = C5(1,1,2,2,2)
        C5(1,2,2,2,1) = C5(1,1,2,2,2)
        C5(2,1,2,2,1) = C5(1,1,2,2,2)
        C5(2,2,1,2,1) = C5(1,1,2,2,2)
        C5(2,2,2,1,1) = C5(1,1,2,2,2)

        C5(2,1,2,2,2) = C5(1,2,2,2,2)
        C5(2,2,1,2,2) = C5(1,2,2,2,2)
        C5(2,2,2,1,2) = C5(1,2,2,2,2)
        C5(2,2,2,2,1) = C5(1,2,2,2,2)
      end if

      if (rank.gt.5) then
        write(*,*) 'rank > 5 not implemented in cCpv12345'
        write(*,*) 'rank = ',rank
        stop
      end if

      end
***********************************************************************
      subroutine cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &              D0,D1,D2,D3,D4,D5,rank,eswitch)
c      entry Dp12345C(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
c     &        D0,D1,D2,D3,D4,D5,C0w0,C1w0,C2w0,C3w0,C4w0,rank,eswitch)
***********************************************************************
*     4-point tensor coefficient functions D0,D1,D2,D3,D4,D5          *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     eswitch should be set to 0                                      *
*     for eswitch=1 some coefficients of rank+1 are calculated        *
*                   as required for the 5-point functions             *
*     definition of integrals:                                        *
*     {D0,D1,D2,D3,D4,D5}=                                            *
*        1/(i pi^2)*\int d^4 q {1,q,q.q,q.q.q,q.q.q.q,q.q.q.q.q}      *
*        1/[q**2 - m02]*[(q+k1)**2-m12]*[(q+k2)**2 - m22]*            *
*        [(q+k3)**2-m32] ]                                            *
*     arguments:   p10=k1**2, p21=(k2-k1)**2, p32=(k3-k2)**2          *
*                      p30=k3**2, p31=(k3-k1)**2, p20=k2**2           *
*     definition of tensor coefficients:                              *
*          D1 = \sum_{i=1}^3 D1(i) *ki                                *
*          D2 = \sum_{i,j=1}^3 D2(i,j) *ki*kj  + D2(0,0) * g          *
*          D3 = \sum_{i,j,k=1}^3 D3(i,j,k) *ki*kj*kk                  *
*              + \sum_{i=1}^3 D3(0,0,i) * (ki*g + 2*sym)              *
*          D4 = \sum_{i,j,k,l=1}^3 D4(i,j,k) *ki*kj*kk*kl             *
*              + \sum_{i,j=1}^3 D4(0,0,i,j) * (ki*kj*g + 5*sym)       *
*              +   D4(0,0,0,0) * (g*g + 2*sym)                        *
*          D5 = \sum_{i,j,k,l,m=1}^3 D4(i,j,k) *ki*kj*kk*kl*km        *
*              + \sum_{i,j,k=1}^3 D4(0,0,i,j,k) * (ki*kj*kk*g +10*sym)*
*              +   \sum_{i=1}^3 D4(0,0,0,0,i) * (g*g*ki + 14*sym)     *
*          C1w0 =  \sum_{i=1}^2 C1w0(i) *( k(i+1)-k1 )                *
*     3-point functions resulting from 4-point functions with         *
*     i-th denominator cancelled: CxwI                                *
***********************************************************************
*    06.12.02 Ansgar Denner    last changed  28.05.04                 *
*                            cosmetic changes 01.06.06 Ansgar Denner  *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p32,p30,p20,p31
      complex*16 q10,q21,q32,q30,q20,q31
      complex*16 m02,m12,m22,m32
      complex*16 mm02,mm12,mm22,mm32,f(3)
      complex*16 k1k2,k1k3,k2k3,zinv(3,3),detz
      complex*16 C0w3,C1w3(2),C2w3(0:2,0:2),C3w3(0:2,0:2,0:2)
     &          ,C4w3(0:2,0:2,0:2,0:2),C5w3(0:2,0:2,0:2,0:2,0:2)
     &          ,C6w3(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w2,C1w2(2),C2w2(0:2,0:2),C3w2(0:2,0:2,0:2)
     &          ,C4w2(0:2,0:2,0:2,0:2),C5w2(0:2,0:2,0:2,0:2,0:2)
     &          ,C6w2(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w1,C1w1(2),C2w1(0:2,0:2),C3w1(0:2,0:2,0:2)
     &          ,C4w1(0:2,0:2,0:2,0:2),C5w1(0:2,0:2,0:2,0:2,0:2)
     &          ,C6w1(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w0,C1w0(2),C2w0(0:2,0:2),C3w0(0:2,0:2,0:2)
     &          ,C4w0(0:2,0:2,0:2,0:2),C5w0(0:2,0:2,0:2,0:2,0:2)
     &          ,C6w0(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0m(0:3),C1m(0:3,3),C2m(0:3,0:3,0:3)
     &          ,C3m(0:3,0:3,0:3,0:3),C4m(0:3,0:3,0:3,0:3,0:3)
      complex*16 S1(3),S2(3,3),S3(3,0:3,0:3),S4(3,0:3,0:3,3),
     &           S5(0:3,0:3,0:3,0:3,0:3)
      complex*16 D0,D1(3),D2(0:3,0:3),D3(0:3,0:3,0:3)
      complex*16 D4(0:3,0:3,0:3,0:3),D5(0:3,0:3,0:3,0:3,0:3)
      complex*16 cD0f
      complex*16 elimcminf2
      integer    rank,eswitch,rankC
      integer    i1,i2,i3,i4,k,l
      integer    ltest,sym

      common /ltest/  ltest
      common /sym/ sym

      data       C1m /12*0d0/,C2m /64*0d0/,C3m/256*0d0/,C4m/1024*0d0/

c      write(*,*)'Dp12345 in',p10,p21,p32,p30,p20,p31,m02,m12,m22,m32
c1     &    ,rank,eswitch


      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cDpv12345'
        write(*,*) 'rank = ',rank
        stop
      end if



c calculate coefficients
      
c      write(*,*) 'cDp0',p10,p21,p32,p30,p20,p31,m02,m12,m22,m32
      D0   = cD0f(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32)
c      write(*,*) 'cDp0',D0

      rankC = max(rank-1,0)

c  cCp must be called also for rank = 0 because of cache!

c      write(*,*) 'cDp1',p21,p32,p31,m12,m22,m32,rank-1
c     shift by k1
      call cCp12345(p21,p32,p31,m12,m22,m32,C0w0,C1w0,C2w0,C3w0,C4w0,
     &                    C5w0,C6w0,              rankC)
c      write(*,*) 'cDp2'
c     shift by k3
c     call cCp12345(p31,p21,p32,m32,m12,m22,C0w0,C1w0,C2w0,C3w0,C4w0,
c     &                   C5w0,C6w0,                rankC)
c     shift by k2
c      call cCp12345(p21,p31,p32,m22,m12,m32,C0w0,C1w0,C2w0,C3w0,C4w0,
c     &                    C5w0,C6w0,               rankC)
      call cCp12345(p20,p32,p30,m02,m22,m32,C0w1,C1w1,C2w1,C3w1,C4w1,
     &                      C5w1,C6w1,             rankC)
      call cCp12345(p10,p31,p30,m02,m12,m32,C0w2,C1w2,C2w2,C3w2,C4w2,
     &                     C5w2,C6w2,              rankC)
      call cCp12345(p10,p21,p20,m02,m12,m22,C0w3,C1w3,C2w3,C3w3,C4w3,
     &                     C5w3,C6w3,              rankC)

      if (rank.eq.0) goto 522

      mm02 = elimcminf2(m02)
      mm12 = elimcminf2(m12)
      mm22 = elimcminf2(m22)
      mm32 = elimcminf2(m32)
      q10  = elimcminf2(p10)
      q21  = elimcminf2(p21)
      q32  = elimcminf2(p32)
      q30  = elimcminf2(p30)
      q31  = elimcminf2(p31)
      q20  = elimcminf2(p20)
 
      k1k2 = (q10+q20-q21)/2D0
      k1k3 = (q10+q30-q31)/2D0
      k2k3 = (q20+q30-q32)/2D0
      detz  = q10*q30*q20+2D0*k1k2*k1k3*k2k3
     &      -q10*k2k3*k2k3-q20*k1k3*k1k3-q30*k1k2*k1k2
      zinv(1,1) = (q30*q20-k2k3*k2k3)/(2d0*detz)
      zinv(1,2) = (k1k3*k2k3-q30*k1k2)/(2d0*detz)
      zinv(1,3) = (k1k2*k2k3-q20*k1k3)/(2d0*detz)
      zinv(2,1) = zinv(1,2)
      zinv(2,2) = (q10*q30-k1k3*k1k3)/(2d0*detz)
      zinv(2,3) = (k1k2*k1k3-q10*k2k3)/(2d0*detz)
      zinv(3,1) = zinv(1,3)
      zinv(3,2) = zinv(2,3)
      zinv(3,3) = (q10*q20-k1k2*k1k2)/(2d0*detz)
      f(1) = q10+mm02-mm12
      f(2) = q20+mm02-mm22
      f(3) = q30+mm02-mm32
 
      C0m(0) = C0w0
      C0m(1) = C0w1
      C0m(2) = C0w2
      C0m(3) = C0w3

c      write(*,*) 'Dad Det ',detz
c      write(*,*) 'Dad inv ',zinv
c      write(*,*) 'Dad C0  ',C0w0,C0w1,C0w2,C0w3


      do 11 k =1,3
        S1(k) = (C0m(k) - C0m(0) - f(k)*D0)
11    continue
      do 16 k =1,3
         D1(k) = 0D0
      do 16 l =1,3
         D1(k) = D1(k) + zinv(k,l)*S1(l)
16    continue

c>      write(*,*) 'cDpv1 ',D1(1),zinv(1,1)*S1(1)
c>     &    ,zinv(1,2)*S1(2) ,zinv(1,3)*S1(3)
c>      write(*,*) 'cDpv1 ',zinv(1,1),S1(1)
c>     &    ,zinv(1,2),S1(2) ,zinv(1,3),S1(3)
c>      write(*,*) 'cDpv1 ',S1(2),C0m(2), - C0m(0),
c>     &    - f(2)*D0 
 
c      write(*,*) 'Dad D1(2)',zinv(2,1)*S1(1),zinv(2,2)*S1(2),
c     &    zinv(2,3)*S1(3)
c      write(*,*) 'Dad S1   ',S1

 
      if (rank.eq.1.and.eswitch.eq.0) goto 522

      D2(0,0) =  mm02*D0+ (C0w0
     &     +f(1)*D1(1)+f(2)*D1(2)+f(3)*D1(3))/2d0

      if (rank.eq.1) goto 522

      C1m(0,2) = C1w0(1)
      C1m(0,3) = C1w0(2)
      C1m(0,1) = -C0m(0) - C1m(0,2) - C1m(0,3)
      C1m(1,2) = C1w1(1)
      C1m(1,3) = C1w1(2)
      C1m(2,1) = C1w2(1)
      C1m(2,3) = C1w2(2)
      C1m(3,1) = C1w3(1)
      C1m(3,2) = C1w3(2)

      do 26 i1=1,3
        do 27 k =1,3
          S2(k,i1) = (C1m(k,i1) - C1m(0,i1) - f(k)*D1(i1))
          if(k.eq.i1)  S2(k,i1) = S2(k,i1) - 2d0*D2(0,0)
 27     continue
      do 26 k =1,i1
        D2(k,i1) = 0D0
      do 26 l =1,3
        D2(k,i1) = D2(k,i1) + zinv(k,l)*S2(l,i1)
26    continue
c      write(*,*) 'cDpv2 ',D2(1,1),zinv(1,1)*S2(1,1)
c     &    ,zinv(1,2)*S2(2,1) ,zinv(1,3)*S2(3,1)
c      write(*,*) 'cDpv2 ',zinv(1,1),S2(1,1)
c     &    ,zinv(1,2),S2(2,1) ,zinv(1,3),S2(3,1)
c      write(*,*) 'cDpv2 ',S2(1,1),C1m(1,1), - C1m(0,1),
c     &    - f(1)*D1(1), -2d0*D2(0,0) 

      if (sym.eq.1) then
        D2(2,1) = D2(1,2)
        D2(3,1) = D2(1,3)
        D2(3,2) = D2(2,3)
      end if

      if (rank.eq.2) then
        if (eswitch.eq.1) then
          D2(3,2) = D2(2,3)      
          do 35 i1 =1,3
            D3(0,0,i1) = mm02*D1(i1)/2 +(C1m(0,i1)
     &           +f(1)*D2(1,i1)+f(2)*D2(i1,2)+f(3)*D2(i1,3))/4d0      
 35       continue
        end if
        goto 522
      end if

      C2m(0,0,0) = C2w0(0,0)
      C2m(1,0,0) = C2w1(0,0)
      C2m(2,0,0) = C2w2(0,0)
      C2m(3,0,0) = C2w3(0,0)

      do 30 k =1,3
        S3(k,0,0) = (C2m(k,0,0) - C2m(0,0,0) - f(k)*D2(0,0))
 30   continue

      do 34 k =1,3
        D3(0,0,k) = 0D0
      do 34 l =1,3
        D3(0 ,0,k) = D3(0 ,0,k) + zinv(k,l)*S3(l,0 ,0 )
 34   continue

c>      do i1 =1,3
c>            D3(0,0,i1) = mm02*D1(i1)/2 +(C1m(0,i1)+
c>     &           +f(1)*D2(1,i1)+f(2)*D2(i1,2)+f(3)*D2(i1,3))/4d0      
c>      end do 

      C2m(0,2,2) = C2w0(1,1)
      C2m(0,2,3) = C2w0(1,2)
      C2m(0,3,3) = C2w0(2,2)
      C2m(0,1,2) = -C1m(0,2) - C2m(0,2,2) - C2m(0,2,3)
      C2m(0,1,3) = -C1m(0,3) - C2m(0,2,3) - C2m(0,3,3)
      C2m(0,1,1) = -C1m(0,1) - C2m(0,1,2) - C2m(0,1,3)
      C2m(1,2,2) = C2w1(1,1)
      C2m(1,2,3) = C2w1(1,2)
      C2m(1,3,3) = C2w1(2,2)
      C2m(2,1,1) = C2w2(1,1)
      C2m(2,1,3) = C2w2(1,2)
      C2m(2,3,3) = C2w2(2,2)
      C2m(3,1,1) = C2w3(1,1)
      C2m(3,1,2) = C2w3(1,2)
      C2m(3,2,2) = C2w3(2,2)
 
      do 38 i1=1,3
      do 38 i2=i1,3
        do 37 k =1,3
          S3(k,i1,i2) = (C2m(k,i1,i2) - C2m(0,i1,i2) - f(k)*D2(i1,i2))
          if(k.eq.i1)  S3(k,i1,i2) = S3(k,i1,i2) - 2d0*D3(0,0,i2)
          if(k.eq.i2)  S3(k,i1,i2) = S3(k,i1,i2) - 2d0*D3(0,0,i1)
 37     continue
      do 38 k =1,i1
        D3(k,i1,i2) = 0D0
      do 38 l =1,3
        D3(k,i1,i2) = D3(k,i1,i2) + zinv(k,l)*S3(l,i1,i2)
38    continue
c      write(*,*) 'cDpv3 ',D3(1,1,1),zinv(1,1)*S3(1,1,1)
c     &    ,zinv(1,2)*S3(2,1,1) ,zinv(1,3)*S3(3,1,1)
c      write(*,*) 'cDpv3 ',zinv(1,1),S3(1,1,1)
c     &    ,zinv(1,2),S3(2,1,1) ,zinv(1,3),S3(3,1,1)
c      write(*,*) 'cDpv3 ',S3(1,1,1),C2m(1,1,1), - C2m(0,1,1),
c     &    - f(1)*D2(1,1), -4d0*D3(0,0,1) 
      if (sym.eq.1) then
        D3(2,1,1) = D3(1,1,2)
        D3(1,2,1) = D3(1,1,2)
        D3(3,1,1) = D3(1,1,3)
        D3(1,3,1) = D3(1,1,3)
        D3(2,1,2) = D3(1,2,2)
        D3(2,2,1) = D3(1,2,2)
        D3(3,1,3) = D3(1,3,3)
        D3(3,3,1) = D3(1,3,3)
        D3(3,3,2) = D3(2,3,3)
        D3(3,2,3) = D3(2,3,3)
        D3(3,2,2) = D3(2,2,3)
        D3(2,3,2) = D3(2,2,3)
        D3(2,1,3) = D3(1,2,3)
        D3(1,3,2) = D3(1,2,3)
        D3(3,1,2) = D3(1,2,3)
        D3(2,3,1) = D3(1,2,3)
        D3(3,2,1) = D3(1,2,3)
      end if

      if (rank.eq.3) then
        if (eswitch.eq.1) then
          D4(0,0,0,0) = 1d0/36d0+ mm02*D2(0,0)/3d0
     &         + (C2w0(0,0) + f(1)*D3(0,0,1) + f(2)*D3(0,0,2)
     &         + f(3)*D3(0,0,3))/6d0
                                                                               
          D3(3,3,2) = D3(2,3,3)
          D3(2,3,2) = D3(2,2,3)
          D3(1,3,2) = D3(1,2,3)
          do 43 i1 =1,3
          do 43 i2 =i1,3
            D4(0,0,i1,i2) = mm02*D2(i1,i2)/3d0 + (C2m(0,i1,i2)
     &           + f(1)*D3(1,i1,i2) + f(2)*D3(i1,i2,2)
     &           + f(3)*D3(i1,i2,3))/6d0
            if (sym.eq.1) then
              D4(0,0,i2,i1) = D4(0,0,i1,i2)
            end if
 43       continue
        end if
        goto 522
      end if                                                      
 
      D4(0,0,0,0) = 1d0/36d0+ mm02*D2(0,0)/3d0 
     & + (C2w0(0,0)+f(1)*D3(0,0,1)+f(2)*D3(0,0,2)+f(3)*D3(0,0,3))/6d0

      C3m(0,0,0,1) = -C2w0(0,0)-C3w0(0,0,1)-C3w0(0,0,2)
      C3m(0,0,0,2) = C3w0(0,0,1)
      C3m(0,0,0,3) = C3w0(0,0,2)
      C3m(1,0,0,2) = C3w1(0,0,1)
      C3m(1,0,0,3) = C3w1(0,0,2)
      C3m(2,0,0,1) = C3w2(0,0,1)
      C3m(2,0,0,3) = C3w2(0,0,2)
      C3m(3,0,0,1) = C3w3(0,0,1)
      C3m(3,0,0,2) = C3w3(0,0,2)
      C3m(0,2,2,2) = C3w0(1,1,1)
      C3m(0,2,2,3) = C3w0(1,1,2)
      C3m(0,2,3,3) = C3w0(1,2,2)
      C3m(0,3,3,3) = C3w0(2,2,2)

      C3m(0,1,2,2) = -C2m(0,2,2) - C3m(0,2,2,2) - C3m(0,2,2,3)
      C3m(0,1,2,3) = -C2m(0,2,3) - C3m(0,2,2,3) - C3m(0,2,3,3)
      C3m(0,1,3,3) = -C2m(0,3,3) - C3m(0,2,3,3) - C3m(0,3,3,3)
      C3m(0,1,1,2) = -C2m(0,1,2) - C3m(0,1,2,2) - C3m(0,1,2,3)
      C3m(0,1,1,3) = -C2m(0,1,3) - C3m(0,1,2,3) - C3m(0,1,3,3)
      C3m(0,1,1,1) = -C2m(0,1,1) - C3m(0,1,1,2) - C3m(0,1,1,3)

      C3m(1,2,2,2) = C3w1(1,1,1)
      C3m(1,2,2,3) = C3w1(1,1,2)
      C3m(1,2,3,3) = C3w1(1,2,2)
      C3m(1,3,3,3) = C3w1(2,2,2)
      C3m(2,1,1,1) = C3w2(1,1,1)
      C3m(2,1,1,3) = C3w2(1,1,2)
      C3m(2,1,3,3) = C3w2(1,2,2)
      C3m(2,3,3,3) = C3w2(2,2,2)
      C3m(3,1,1,1) = C3w3(1,1,1)
      C3m(3,1,1,2) = C3w3(1,1,2)
      C3m(3,1,2,2) = C3w3(1,2,2)
      C3m(3,2,2,2) = C3w3(2,2,2)

c     do 43 i1 =1,3
c     do 43 i2 =i1,3
c       D4(0,0,i1,i2) = mm02*D2(i1,i2)/3d0 + (C2m(0,i1,i2)
c    &     +f(1)*D3(1,i1,i2)+f(2)*D3(i1,i2,2)+f(3)*D3(i1,i2,3))/6d0      
c43   continue
c     D4(0,0,2,1) = D4(0,0,1,2)
c     D4(0,0,3,1) = D4(0,0,1,3)
c     D4(0,0,3,2) = D4(0,0,2,3)

c version with  zinv(k,l) numerically more stable!!!!!
      do 45 i1=1,3
        do 44 k =1,3
          S4(k,0,0,i1) = C3m(k,0,0,i1) - C3m(0,0,0,i1)
     &         - f(k)*D3(0,0,i1)
          if(k.eq.i1)  S4(k,0,0,i1) = S4(k,0,0,i1) - 2d0*D4(0,0,0,0)
 44     continue
      do 45 k =1,i1
         D4(0,0,k,i1) = 0D0
      do 45 l =1,3
         D4(0,0,k,i1) = D4(0,0,k,i1) + zinv(k,l)*S4(l, 0, 0,i1)
 45   continue
      if (sym.eq.1) then
        D4(0,0,2,1) = D4(0,0,1,2)
        D4(0,0,3,1) = D4(0,0,1,3)
        D4(0,0,3,2) = D4(0,0,2,3)
      end if

      do 49 i1=1,3
      do 49 i2=i1,3
      do 49 i3=i2,3
        do 47 k =1,3
          S4(k,i1,i2,i3) = C3m(k,i1,i2,i3) - C3m(0,i1,i2,i3)
     &         - f(k)*D3(i1,i2,i3)
          if(k.eq.i1) S4(k,i1,i2,i3)= S4(k,i1,i2,i3) - 2d0*D4(0,0,i2,i3)
          if(k.eq.i2) S4(k,i1,i2,i3)= S4(k,i1,i2,i3) - 2d0*D4(0,0,i1,i3)
          if(k.eq.i3) S4(k,i1,i2,i3)= S4(k,i1,i2,i3) - 2d0*D4(0,0,i1,i2)
 47     continue
        do 46 k =1,i1
          D4(k,i1,i2,i3) = 0D0
        do 46 l =1,3
          D4(k,i1,i2,i3) = D4(k,i1,i2,i3) + zinv(k,l)*S4(l,i1,i2,i3)
 46     continue
        if (sym.eq.1) then
          do 149 k =1,i1
            D4(i1,k,i2,i3) = D4(k,i1,i2,i3)
            D4(i1,i2,k,i3) = D4(k,i1,i2,i3)
            D4(i1,i2,i3,k) = D4(k,i1,i2,i3)
            D4(k,i1,i3,i2) = D4(k,i1,i2,i3)
            D4(i1,k,i3,i2) = D4(k,i1,i2,i3)
            D4(i1,i3,k,i2) = D4(k,i1,i2,i3)
            D4(i1,i3,i2,k) = D4(k,i1,i2,i3)
            D4(k,i2,i3,i1) = D4(k,i1,i2,i3)
            D4(i2,k,i3,i1) = D4(k,i1,i2,i3)
            D4(i2,i3,k,i1) = D4(k,i1,i2,i3)
            D4(i2,i3,i1,k) = D4(k,i1,i2,i3)
            D4(k,i2,i1,i3) = D4(k,i1,i2,i3)
            D4(i2,k,i1,i3) = D4(k,i1,i2,i3)
            D4(i2,i1,k,i3) = D4(k,i1,i2,i3)
            D4(i2,i1,i3,k) = D4(k,i1,i2,i3)
            D4(k,i3,i1,i2) = D4(k,i1,i2,i3)
            D4(i3,k,i1,i2) = D4(k,i1,i2,i3)
            D4(i3,i1,k,i2) = D4(k,i1,i2,i3)
            D4(i3,i1,i2,k) = D4(k,i1,i2,i3)
            D4(k,i3,i2,i1) = D4(k,i1,i2,i3)
            D4(i3,k,i2,i1) = D4(k,i1,i2,i3)
            D4(i3,i2,k,i1) = D4(k,i1,i2,i3)
            D4(i3,i2,i1,k) = D4(k,i1,i2,i3)
 149      continue
        end if
 49   continue


      if (rank.eq.4) then
        if (eswitch.eq.1) then
          
          C3m(0,0,0,1) = -C2w0(0,0)-C3w0(0,0,1)-C3w0(0,0,2)
          C3m(0,0,0,2) = C3w0(0,0,1)
          C3m(0,0,0,3) = C3w0(0,0,2)
          
          D4(0,0,3,2)=D4(0,0,2,3)
          do 52 i1 =1,3
            D5(0,0,0,0,i1) = -1d0/192d0+mm02*D3(0,0,i1)/4d0 
     &           + (C3m(0,0,0,i1) + f(1)*D4(0,0,1,i1)
     &           + f(2)*D4(0,0,i1,2) + f(3)*D4(0,0,i1,3))/8d0
 52       continue
          
          D4(1,1,3,2) = D4(1,1,2,3)
          D4(1,2,3,2) = D4(1,2,2,3)
          D4(1,3,3,2) = D4(1,2,3,3)
          D4(2,2,3,2) = D4(2,2,2,3)
          D4(2,3,3,2) = D4(2,2,3,3)
          D4(3,3,3,2) = D4(2,3,3,3)
          do 55 i1 =1,3
          do 55 i2 =i1,3
          do 55 i3 =i2,3
            D5(0,0,i1,i2,i3) = mm02*D3(i1,i2,i3)/4d0 + (C3m(0,i1,i2,i3)
     &           +f(1)*D4(1,i1,i2,i3)+f(2)*D4(i1,i2,i3,2)
     &           +f(3)*D4(i1,i2,i3,3))/8d0
 55       continue
          if (sym.eq.1) then
            D5(0,0,2,1,1) = D5(0,0,1,1,2)
            D5(0,0,1,2,1) = D5(0,0,1,1,2)
            D5(0,0,3,1,1) = D5(0,0,1,1,3)
            D5(0,0,1,3,1) = D5(0,0,1,1,3)
            D5(0,0,2,1,2) = D5(0,0,1,2,2)
            D5(0,0,2,2,1) = D5(0,0,1,2,2)
            D5(0,0,3,1,3) = D5(0,0,1,3,3)
            D5(0,0,3,3,1) = D5(0,0,1,3,3)
            D5(0,0,3,3,2) = D5(0,0,2,3,3)
            D5(0,0,3,2,3) = D5(0,0,2,3,3)
            D5(0,0,3,2,2) = D5(0,0,2,2,3)
            D5(0,0,2,3,2) = D5(0,0,2,2,3)
            D5(0,0,2,1,3) = D5(0,0,1,2,3)
            D5(0,0,1,3,2) = D5(0,0,1,2,3)
            D5(0,0,3,1,2) = D5(0,0,1,2,3)
            D5(0,0,2,3,1) = D5(0,0,1,2,3)
            D5(0,0,3,2,1) = D5(0,0,1,2,3)
          end if
        end if
        goto 522
      end if

      C4m(0,0,0,0,0) = C4w0(0,0,0,0)
      C4m(1,0,0,0,0) = C4w1(0,0,0,0)
      C4m(2,0,0,0,0) = C4w2(0,0,0,0)
      C4m(3,0,0,0,0) = C4w3(0,0,0,0)

      C4m(0,0,0,2,2) = C4w0(0,0,1,1)
      C4m(0,0,0,2,3) = C4w0(0,0,1,2)
      C4m(0,0,0,3,3) = C4w0(0,0,2,2)
      C4m(0,0,0,1,2) = -C3w0(0,0,1)-C4w0(0,0,1,1)-C4w0(0,0,1,2)
      C4m(0,0,0,1,3) = -C3w0(0,0,2)-C4w0(0,0,1,2)-C4w0(0,0,2,2)
      C4m(0,0,0,1,1) = -C3m(0,0,0,1)-C4m(0,0,0,1,2)-C4m(0,0,0,1,3)

      C4m(1,0,0,2,2) = C4w1(0,0,1,1)
      C4m(1,0,0,2,3) = C4w1(0,0,1,2)
      C4m(1,0,0,3,3) = C4w1(0,0,2,2)
      C4m(2,0,0,1,1) = C4w2(0,0,1,1)
      C4m(2,0,0,1,3) = C4w2(0,0,1,2)
      C4m(2,0,0,3,3) = C4w2(0,0,2,2)
      C4m(3,0,0,1,1) = C4w3(0,0,1,1)
      C4m(3,0,0,1,2) = C4w3(0,0,1,2)
      C4m(3,0,0,2,2) = C4w3(0,0,2,2)

      C4m(0,2,2,2,2) = C4w0(1,1,1,1)
      C4m(0,2,2,2,3) = C4w0(1,1,1,2)
      C4m(0,2,2,3,3) = C4w0(1,1,2,2)
      C4m(0,2,3,3,3) = C4w0(1,2,2,2)
      C4m(0,3,3,3,3) = C4w0(2,2,2,2)

      C4m(0,1,2,2,2) = -C3m(0,2,2,2) - C4m(0,2,2,2,2) - C4m(0,2,2,2,3)
      C4m(0,1,2,2,3) = -C3m(0,2,2,3) - C4m(0,2,2,2,3) - C4m(0,2,2,3,3)
      C4m(0,1,2,3,3) = -C3m(0,2,3,3) - C4m(0,2,2,3,3) - C4m(0,2,3,3,3)
      C4m(0,1,3,3,3) = -C3m(0,3,3,3) - C4m(0,2,3,3,3) - C4m(0,3,3,3,3)
      C4m(0,1,1,2,2) = -C3m(0,1,2,2) - C4m(0,1,2,2,2) - C4m(0,1,2,2,3)
      C4m(0,1,1,2,3) = -C3m(0,1,2,3) - C4m(0,1,2,2,3) - C4m(0,1,2,3,3)
      C4m(0,1,1,3,3) = -C3m(0,1,3,3) - C4m(0,1,2,3,3) - C4m(0,1,3,3,3)
      C4m(0,1,1,1,2) = -C3m(0,1,1,2) - C4m(0,1,1,2,2) - C4m(0,1,1,2,3)
      C4m(0,1,1,1,3) = -C3m(0,1,1,3) - C4m(0,1,1,2,3) - C4m(0,1,1,3,3)
      C4m(0,1,1,1,1) = -C3m(0,1,1,1) - C4m(0,1,1,1,2) - C4m(0,1,1,1,3)

      C4m(1,2,2,2,2) = C4w1(1,1,1,1)
      C4m(1,2,2,2,3) = C4w1(1,1,1,2)
      C4m(1,2,2,3,3) = C4w1(1,1,2,2)
      C4m(1,2,3,3,3) = C4w1(1,2,2,2)
      C4m(1,3,3,3,3) = C4w1(2,2,2,2)
      C4m(2,1,1,1,1) = C4w2(1,1,1,1)
      C4m(2,1,1,1,3) = C4w2(1,1,1,2)
      C4m(2,1,1,3,3) = C4w2(1,1,2,2)
      C4m(2,1,3,3,3) = C4w2(1,2,2,2)
      C4m(2,3,3,3,3) = C4w2(2,2,2,2)
      C4m(3,1,1,1,1) = C4w3(1,1,1,1)
      C4m(3,1,1,1,2) = C4w3(1,1,1,2)
      C4m(3,1,1,2,2) = C4w3(1,1,2,2)
      C4m(3,1,2,2,2) = C4w3(1,2,2,2)
      C4m(3,2,2,2,2) = C4w3(2,2,2,2)

c     do 53 i1 =1,3
c       D5(0,0,0,0,i1) = -1d0/192d0+mm02*D3(0,0,i1)/4d0 + (C3m(0,0,0,i1)
c    &     +f(1)*D4(0,0,1,i1)+f(2)*D4(0,0,i1,2)+f(3)*D4(0,0,i1,3))/8d0 
c53   continue
      
      do 51 k =1,3
        S5(k,0,0,0,0) = C4m(k,0,0,0,0) - C4m(0,0,0,0,0)
     &                   - f(k)*D4(0,0,0,0)
51    continue

      do 54 k =1,3
        D5(0,0,0,0,k) = 0D0
      do 54 l =1,3
        D5(0,0,0,0,k) = D5(0,0,0,0,k) + zinv(k,l)*S5(l,0,0,0,0)
 54   continue

c     do 58 i1 =1,3
c     do 58 i2 =i1,3
c     do 58 i3 =i2,3
c       D5(0,0,i1,i2,i3) = mm02*D3(i1,i2,i3)/4d0 + (C3m(0,i1,i2,i3)
c    &      +f(1)*D4(1,i1,i2,i3)+f(2)*D4(i1,i2,i3,2)
c    &      +f(3)*D4(i1,i2,i3,3))/8d0      
c58   continue

      do 60 i1=1,3
      do 60 i2=i1,3
        do 57 k =1,3
          S5(k,0,0,i1,i2) = C4m(k,0,0,i1,i2) - C4m(0,0,0,i1,i2)
     &         - f(k)*D4(0,0,i1,i2)
          if(k.eq.i1)  S5(k,0,0,i1,i2) = S5(k,0,0,i1,i2)
     &         - 2d0*D5(0,0,0,0,i2)
          if(k.eq.i2)  S5(k,0,0,i1,i2) = S5(k,0,0,i1,i2)
     &         - 2d0*D5(0,0,0,0,i1)
 57     continue
        do 56 k =1,i1
          D5(0,0,k,i1,i2) = 0D0
          do 56 l =1,3
            D5(0,0,k,i1,i2)= D5(0,0,k,i1,i2) + zinv(k,l)*S5(l,0,0,i1,i2)
 56       continue
 60   continue
      if (sym.eq.1) then
        D5(0,0,2,1,1) = D5(0,0,1,1,2)
        D5(0,0,1,2,1) = D5(0,0,1,1,2)
        D5(0,0,3,1,1) = D5(0,0,1,1,3)
        D5(0,0,1,3,1) = D5(0,0,1,1,3)
        D5(0,0,2,1,2) = D5(0,0,1,2,2)
        D5(0,0,2,2,1) = D5(0,0,1,2,2)
        D5(0,0,3,1,3) = D5(0,0,1,3,3)
        D5(0,0,3,3,1) = D5(0,0,1,3,3)
        D5(0,0,3,3,2) = D5(0,0,2,3,3)
        D5(0,0,3,2,3) = D5(0,0,2,3,3)
        D5(0,0,3,2,2) = D5(0,0,2,2,3)
        D5(0,0,2,3,2) = D5(0,0,2,2,3)
        D5(0,0,2,1,3) = D5(0,0,1,2,3)
        D5(0,0,1,3,2) = D5(0,0,1,2,3)
        D5(0,0,3,1,2) = D5(0,0,1,2,3)
        D5(0,0,2,3,1) = D5(0,0,1,2,3)
        D5(0,0,3,2,1) = D5(0,0,1,2,3)
      end if

      do 61 i1=1,3
      do 61 i2=i1,3
      do 61 i3=i2,3
      do 61 i4=i3,3
        do 59 k =1,3
          S5(k,i1,i2,i3,i4) = C4m(k,i1,i2,i3,i4) - C4m(0,i1,i2,i3,i4)
     &         - f(k)*D4(i1,i2,i3,i4)
          if(k.eq.i1)  S5(k,i1,i2,i3,i4) = S5(k,i1,i2,i3,i4)
     &         - 2d0*D5(0,0,i2,i3,i4)
          if(k.eq.i2)  S5(k,i1,i2,i3,i4) = S5(k,i1,i2,i3,i4)
     &         - 2d0*D5(0,0,i1,i3,i4)
          if(k.eq.i3)  S5(k,i1,i2,i3,i4) = S5(k,i1,i2,i3,i4)
     &         - 2d0*D5(0,0,i1,i2,i4)
          if(k.eq.i4)  S5(k,i1,i2,i3,i4) = S5(k,i1,i2,i3,i4)
     &         - 2d0*D5(0,0,i1,i2,i3)
 59     continue
        do 63 k =1,i1
          D5(k,i1,i2,i3,i4) = 0D0
        do 63 l =1,3
          D5(k,i1,i2,i3,i4) = D5(k,i1,i2,i3,i4) 
     &         + zinv(k,l)*S5(l,i1,i2,i3,i4)
 63     continue

        if (sym.eq.1) then
          do 62 k =1,i1
c            D5(i1,i2,i3,i4,k) = D5(k,i1,i2,i3,i4)
            D5(k,i1,i2,i4,i3) = D5(k,i1,i2,i3,i4)
            D5(k,i1,i3,i2,i4) = D5(k,i1,i2,i3,i4)
            D5(k,i1,i3,i4,i2) = D5(k,i1,i2,i3,i4)
            D5(k,i1,i4,i2,i3) = D5(k,i1,i2,i3,i4)
            D5(k,i1,i4,i3,i2) = D5(k,i1,i2,i3,i4)
            D5(k,i2,i1,i3,i4) = D5(k,i1,i2,i3,i4)
            D5(k,i2,i1,i4,i3) = D5(k,i1,i2,i3,i4)
            D5(k,i2,i3,i1,i4) = D5(k,i1,i2,i3,i4)
            D5(k,i2,i3,i4,i1) = D5(k,i1,i2,i3,i4)
            D5(k,i2,i4,i1,i3) = D5(k,i1,i2,i3,i4)
            D5(k,i2,i4,i3,i1) = D5(k,i1,i2,i3,i4)
            D5(k,i3,i1,i2,i4) = D5(k,i1,i2,i3,i4)
            D5(k,i3,i1,i4,i2) = D5(k,i1,i2,i3,i4)
            D5(k,i3,i2,i1,i4) = D5(k,i1,i2,i3,i4)
            D5(k,i3,i2,i4,i1) = D5(k,i1,i2,i3,i4)
            D5(k,i3,i4,i1,i2) = D5(k,i1,i2,i3,i4)
            D5(k,i3,i4,i2,i1) = D5(k,i1,i2,i3,i4)
            D5(k,i4,i1,i3,i2) = D5(k,i1,i2,i3,i4)
            D5(k,i4,i1,i2,i3) = D5(k,i1,i2,i3,i4)
            D5(k,i4,i3,i1,i2) = D5(k,i1,i2,i3,i4)
            D5(k,i4,i3,i2,i1) = D5(k,i1,i2,i3,i4)
            D5(k,i4,i2,i1,i3) = D5(k,i1,i2,i3,i4)
            D5(k,i4,i2,i3,i1) = D5(k,i1,i2,i3,i4)
            
            D5(i1,k,i2,i3,i4) = D5(k,i1,i2,i3,i4)
            D5(i1,k,i2,i4,i3) = D5(k,i1,i2,i3,i4)
            D5(i1,k,i3,i2,i4) = D5(k,i1,i2,i3,i4)
            D5(i1,k,i3,i4,i2) = D5(k,i1,i2,i3,i4)
            D5(i1,k,i4,i2,i3) = D5(k,i1,i2,i3,i4)
            D5(i1,k,i4,i3,i2) = D5(k,i1,i2,i3,i4)
            D5(i2,k,i1,i3,i4) = D5(k,i1,i2,i3,i4)
            D5(i2,k,i1,i4,i3) = D5(k,i1,i2,i3,i4)
            D5(i2,k,i3,i1,i4) = D5(k,i1,i2,i3,i4)
            D5(i2,k,i3,i4,i1) = D5(k,i1,i2,i3,i4)
            D5(i2,k,i4,i1,i3) = D5(k,i1,i2,i3,i4)
            D5(i2,k,i4,i3,i1) = D5(k,i1,i2,i3,i4)
            D5(i3,k,i1,i2,i4) = D5(k,i1,i2,i3,i4)
            D5(i3,k,i1,i4,i2) = D5(k,i1,i2,i3,i4)
            D5(i3,k,i2,i1,i4) = D5(k,i1,i2,i3,i4)
            D5(i3,k,i2,i4,i1) = D5(k,i1,i2,i3,i4)
            D5(i3,k,i4,i1,i2) = D5(k,i1,i2,i3,i4)
            D5(i3,k,i4,i2,i1) = D5(k,i1,i2,i3,i4)
            D5(i4,k,i1,i3,i2) = D5(k,i1,i2,i3,i4)
            D5(i4,k,i1,i2,i3) = D5(k,i1,i2,i3,i4)
            D5(i4,k,i3,i1,i2) = D5(k,i1,i2,i3,i4)
            D5(i4,k,i3,i2,i1) = D5(k,i1,i2,i3,i4)
            D5(i4,k,i2,i1,i3) = D5(k,i1,i2,i3,i4)
            D5(i4,k,i2,i3,i1) = D5(k,i1,i2,i3,i4)
            
            D5(i1,i2,k,i3,i4) = D5(k,i1,i2,i3,i4)
            D5(i1,i2,k,i4,i3) = D5(k,i1,i2,i3,i4)
            D5(i1,i3,k,i2,i4) = D5(k,i1,i2,i3,i4)
            D5(i1,i3,k,i4,i2) = D5(k,i1,i2,i3,i4)
            D5(i1,i4,k,i2,i3) = D5(k,i1,i2,i3,i4)
            D5(i1,i4,k,i3,i2) = D5(k,i1,i2,i3,i4)
            D5(i2,i1,k,i3,i4) = D5(k,i1,i2,i3,i4)
            D5(i2,i1,k,i4,i3) = D5(k,i1,i2,i3,i4)
            D5(i2,i3,k,i1,i4) = D5(k,i1,i2,i3,i4)
            D5(i2,i3,k,i4,i1) = D5(k,i1,i2,i3,i4)
            D5(i2,i4,k,i1,i3) = D5(k,i1,i2,i3,i4)
            D5(i2,i4,k,i3,i1) = D5(k,i1,i2,i3,i4)
            D5(i3,i1,k,i2,i4) = D5(k,i1,i2,i3,i4)
            D5(i3,i1,k,i4,i2) = D5(k,i1,i2,i3,i4)
            D5(i3,i2,k,i1,i4) = D5(k,i1,i2,i3,i4)
            D5(i3,i2,k,i4,i1) = D5(k,i1,i2,i3,i4)
            D5(i3,i4,k,i1,i2) = D5(k,i1,i2,i3,i4)
            D5(i3,i4,k,i2,i1) = D5(k,i1,i2,i3,i4)
            D5(i4,i1,k,i3,i2) = D5(k,i1,i2,i3,i4)
            D5(i4,i1,k,i2,i3) = D5(k,i1,i2,i3,i4)
            D5(i4,i3,k,i1,i2) = D5(k,i1,i2,i3,i4)
            D5(i4,i3,k,i2,i1) = D5(k,i1,i2,i3,i4)
            D5(i4,i2,k,i1,i3) = D5(k,i1,i2,i3,i4)
            D5(i4,i2,k,i3,i1) = D5(k,i1,i2,i3,i4)
            
            D5(i1,i2,i3,k,i4) = D5(k,i1,i2,i3,i4)
            D5(i1,i2,i4,k,i3) = D5(k,i1,i2,i3,i4)
            D5(i1,i3,i2,k,i4) = D5(k,i1,i2,i3,i4)
            D5(i1,i3,i4,k,i2) = D5(k,i1,i2,i3,i4)
            D5(i1,i4,i2,k,i3) = D5(k,i1,i2,i3,i4)
            D5(i1,i4,i3,k,i2) = D5(k,i1,i2,i3,i4)
            D5(i2,i1,i3,k,i4) = D5(k,i1,i2,i3,i4)
            D5(i2,i1,i4,k,i3) = D5(k,i1,i2,i3,i4)
            D5(i2,i3,i1,k,i4) = D5(k,i1,i2,i3,i4)
            D5(i2,i3,i4,k,i1) = D5(k,i1,i2,i3,i4)
            D5(i2,i4,i1,k,i3) = D5(k,i1,i2,i3,i4)
            D5(i2,i4,i3,k,i1) = D5(k,i1,i2,i3,i4)
            D5(i3,i1,i2,k,i4) = D5(k,i1,i2,i3,i4)
            D5(i3,i1,i4,k,i2) = D5(k,i1,i2,i3,i4)
            D5(i3,i2,i1,k,i4) = D5(k,i1,i2,i3,i4)
            D5(i3,i2,i4,k,i1) = D5(k,i1,i2,i3,i4)
            D5(i3,i4,i1,k,i2) = D5(k,i1,i2,i3,i4)
            D5(i3,i4,i2,k,i1) = D5(k,i1,i2,i3,i4)
            D5(i4,i1,i3,k,i2) = D5(k,i1,i2,i3,i4)
            D5(i4,i1,i2,k,i3) = D5(k,i1,i2,i3,i4)
            D5(i4,i3,i1,k,i2) = D5(k,i1,i2,i3,i4)
            D5(i4,i3,i2,k,i1) = D5(k,i1,i2,i3,i4)
            D5(i4,i2,i1,k,i3) = D5(k,i1,i2,i3,i4)
            D5(i4,i2,i3,k,i1) = D5(k,i1,i2,i3,i4)

            D5(i1,i2,i3,i4,k) = D5(k,i1,i2,i3,i4)
            D5(i1,i2,i4,i3,k) = D5(k,i1,i2,i3,i4)
            D5(i1,i3,i2,i4,k) = D5(k,i1,i2,i3,i4)
            D5(i1,i3,i4,i2,k) = D5(k,i1,i2,i3,i4)
            D5(i1,i4,i2,i3,k) = D5(k,i1,i2,i3,i4)
            D5(i1,i4,i3,i2,k) = D5(k,i1,i2,i3,i4)
            D5(i2,i1,i3,i4,k) = D5(k,i1,i2,i3,i4)
            D5(i2,i1,i4,i3,k) = D5(k,i1,i2,i3,i4)
            D5(i2,i3,i1,i4,k) = D5(k,i1,i2,i3,i4)
            D5(i2,i3,i4,i1,k) = D5(k,i1,i2,i3,i4)
            D5(i2,i4,i1,i3,k) = D5(k,i1,i2,i3,i4)
            D5(i2,i4,i3,i1,k) = D5(k,i1,i2,i3,i4)
            D5(i3,i1,i2,i4,k) = D5(k,i1,i2,i3,i4)
            D5(i3,i1,i4,i2,k) = D5(k,i1,i2,i3,i4)
            D5(i3,i2,i1,i4,k) = D5(k,i1,i2,i3,i4)
            D5(i3,i2,i4,i1,k) = D5(k,i1,i2,i3,i4)
            D5(i3,i4,i1,i2,k) = D5(k,i1,i2,i3,i4)
            D5(i3,i4,i2,i1,k) = D5(k,i1,i2,i3,i4)
            D5(i4,i1,i3,i2,k) = D5(k,i1,i2,i3,i4)
            D5(i4,i1,i2,i3,k) = D5(k,i1,i2,i3,i4)
            D5(i4,i3,i1,i2,k) = D5(k,i1,i2,i3,i4)
            D5(i4,i3,i2,i1,k) = D5(k,i1,i2,i3,i4)
            D5(i4,i2,i1,i3,k) = D5(k,i1,i2,i3,i4)
            D5(i4,i2,i3,i1,k) = D5(k,i1,i2,i3,i4)
            
 62       continue
        end if
 61   continue


      if (rank.gt.5.or.rank.eq.5.and.eswitch.eq.1) then
        write(*,*) 'rank > 5 not implemented in cDpv12345'
        write(*,*) 'rank = ',rank
        write(*,*) 'eswitch = ',eswitch
        stop
      end if

 522  continue

      end
***********************************************************************
      subroutine cEp1234(p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
     &                                    m02,m12,m22,m32,m42,
     &                   En0,En1,En2,En3,En4,rankin)
c      entry Ep1234D(p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
c     &                                    m02,m12,m22,m32,m42,
c     &              En0,En1,En2,En3,En4,D0w0,D1w0,D2w0,D3w0,rankin)
***********************************************************************
*     5-point tensor coefficient functions  E0,E1,E2,E3,E4            *
*     up to 4 integration momenta in numerator                        *
*     g(mu,nu) terms in decomposition of En...                        *
*     version that needs only D's of rank-1                           *
*     coefficients calculated up to *rank* momenta in numerator       *
*     definition of integrals:                                        *
*     {E0,E1,E2,E3,E4}= 1/(i pi^2)*\int d^4 q {1,q,q.q,q.q.q,q.q.q.q} *
*        1/[(q**2 - m02)*[(q+k1)**2-m12]*[(q+k2)**2 - m22]*           *
*        [(q+k3)**2-m32]*[(q+k4)**2-m42] ]                            *
*     arguments:  p10=k1**2, p21=(k2-k1)**2, p32=(k3-k2)**2,          *
*               p43=(k4-k3)**2, p40=k4**2, p20=k2**2, p31=(k3-k1)**2, *
*               p42=(k4-k2)**2, p30=k3**2, p41=(k4-k1)**2             *
*     definition of tensor coefficients:                              *
*          E1 = \sum_{i=1}^4 En1(i) *ki                               *
*          E2 = \sum_{i,j=1}^4 En2(i,j) *ki*kj  + En2(0,0) * g        *
*          E3 = \sum_{i,j,k=1}^4 En3(i,j,k) *ki*kj*kk                 *
*              + \sum_{i=1}^4 En3(0,0,i) * (ki*g + 2*sym)             *
*          E4 = \sum_{i,j,k,l=1}^4 En4(i,j,k) *ki*kj*kk*kl            *
*              + \sum_{i,j=1}^4 En4(0,0,i,j) * (ki*kj*g + 5*sym)      *
*              +   En4(0,0,0,0) * (g*g + 2*sym)                       *
*          D1w0 =  \sum_{i=1}^3 D1w0(i) *( k(i+1)-k1 )                *
*     4-point functions resulting from 5-point functions with         *
*     i-th denominator cancelled: DXwI                                *
***********************************************************************
*    08.06.05 Ansgar Denner    last changed  12.12.05                 *
*                            cosmetic changes 01.06.06 Ansgar Denner  *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p32,p43,p40,p20,p31,p42,p30,p41
      complex*16 m02,m12,m22,m32,m42
      complex*16 q10,q21,q32,q43,q40,q20,q31,q42,q30,q41
c      complex*16 k1k1,k1k2,k1k3,k1k4,k2k2,k2k3,k2k4,k3k3,k3k4,k4k4
c      complex*16 detz
      complex*16 mm02,mm12,mm22,mm32,mm42
      complex*16 mx(0:4,0:4),mxinv(0:4,0:4)
c      complex*16 z(4,4),zinv(4,4)
      complex*16 D0w4,D1w4(3),D2w4(0:3,0:3),D3w4(0:3,0:3,0:3)
      complex*16 D0w3,D1w3(3),D2w3(0:3,0:3),D3w3(0:3,0:3,0:3)
      complex*16 D0w2,D1w2(3),D2w2(0:3,0:3),D3w2(0:3,0:3,0:3)
      complex*16 D0w1,D1w1(3),D2w1(0:3,0:3),D3w1(0:3,0:3,0:3)
      complex*16 D0w0,D1w0(3),D2w0(0:3,0:3),D3w0(0:3,0:3,0:3)
      complex*16 D4w4(0:3,0:3,0:3,0:3),D5w4(0:3,0:3,0:3,0:3,0:3)
      complex*16 D4w3(0:3,0:3,0:3,0:3),D5w3(0:3,0:3,0:3,0:3,0:3)
      complex*16 D4w2(0:3,0:3,0:3,0:3),D5w2(0:3,0:3,0:3,0:3,0:3)
      complex*16 D4w1(0:3,0:3,0:3,0:3),D5w1(0:3,0:3,0:3,0:3,0:3)
      complex*16 D4w0(0:3,0:3,0:3,0:3),D5w0(0:3,0:3,0:3,0:3,0:3)
      complex*16 D0m(0:4),D1m(0:4,4),D2m(0:4,1:4,1:4)
     &          ,D3m(0:4,1:4,1:4,1:4),D4m(0:4,1:4,1:4,1:4,1:4)
      complex*16 D2m00(0:4),D3m00(0:4,4),D4m00(0:4,1:4,1:4)
     &          ,D5m00(0:4,4,4,4)
     &          ,D4m0000(0:4),D5m0000(0:4,4)
      complex*16 En0,En1(4),En2(0:4,0:4),En3(0:4,0:4,0:4)
      complex*16 En4(0:4,0:4,0:4,0:4)
      complex*16 elimcminf2
      integer    rankin,switch,rank,rankD
      integer    i,j,i1,i2,i3,i4,k,l,m
      integer    ltest,sym

      common /ltest/  ltest
      common /sym/    sym
 
      data       D1m /20*0d0/,D2m /80*0d0/,D3m/320*0d0/,D4m/1280*0d0/
      data       D2m00/5*0d0/,D3m00/20*0d0/,D4m00/80*0d0/,D5m00/320*0d0/
      data       D4m0000/5*0d0/,D5m0000/20*0d0/

      complex*16 fct(86),x(15)
      logical    nocalc
      integer    type
      character  name*11

c>      integer maxt,maxf,maxn,maxfct,maxx
c>      parameter(maxt=5,maxf=500,maxn=1500,maxfct=86,maxx=15)
c>c fast
c>      complex*16 f(maxt,maxf,maxfct),arg(maxt,maxf,maxx)
c>      integer grank(maxt,maxf),gswitch(maxt,maxf)
c>      integer lrank(maxt,maxf),lswitch(maxt,maxf)
c>      integer ncalc(maxt),calc(maxt,maxn)
c>      integer pointer(maxt,maxn),ncall(maxt),n
c>      character*11 cname(maxt)
c>      integer usecache,usecachesave,cachelevel,cacheleveli
c>      common/cache/f,arg,grank,lrank,gswitch,lswitch,
c>     &             ncalc,calc,pointer,ncall,n,cname
c>      common/cachem/usecache,usecachesave,cachelevel,cacheleveli
c>
c>      integer countE
c>      data countE /0/
c>      common /countE/ countE
c>
c>      countE = countE+1

c      write(*,*)'Dp12345 in',p10,p21,p32,p30,p20,p31,m02,m12,m22,m32
c1     &    ,rank,eswitch


c      write(*,*) 'cEp in', p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
c     &                                    m02,m12,m22,m32,m42

      if (rankin.lt.0) then
        write(*,*) 'rank < 0 not implemented in cEp1234'
        write(*,*) 'rank = ',rank
      end if

      x(1)=p10
      x(2)=p21
      x(3)=p32
      x(4)=p43
      x(5)=p40
      x(6)=p20
      x(7)=p31
      x(8)=p42
      x(9)=p30
      x(10)=p41
      x(11)=m02
      x(12)=m12
      x(13)=m22
      x(14)=m32
      x(15)=m42
      type = 1
      rank = rankin
      switch = 0
      name='cEp1234'

c      rank = 3

      call cacheread(fct,x,86,15,type,rank,switch,name,nocalc)
      if(nocalc)then
        En0=fct(1)
        En1(1)=fct(2)
        En1(2)=fct(3)
        En1(3)=fct(4)
        En1(4)=fct(5)
        if (rank.eq.1) goto 511
        En2(0,0)=fct(6)
        En2(1,1)=fct(7)
        En2(1,2)=fct(8)
        En2(1,3)=fct(9)
        En2(1,4)=fct(10)
        En2(2,2)=fct(11)
        En2(2,3)=fct(12)
        En2(2,4)=fct(13)
        En2(3,3)=fct(14)
        En2(3,4)=fct(15)
        En2(4,4)=fct(16)
        if (rank.eq.2) goto 511
        En3(0,0,1)=fct(17)
        En3(0,0,2)=fct(18)
        En3(0,0,3)=fct(19)
        En3(0,0,4)=fct(20)
        En3(1,1,1)=fct(21)
        En3(1,1,2)=fct(22)
        En3(1,1,3)=fct(23)
        En3(1,1,4)=fct(24)
        En3(1,2,2)=fct(25)
        En3(1,2,3)=fct(26)
        En3(1,2,4)=fct(27)
        En3(1,3,3)=fct(28)
        En3(1,3,4)=fct(29)
        En3(1,4,4)=fct(30)
        En3(2,2,2)=fct(31)
        En3(2,2,3)=fct(32)
        En3(2,2,4)=fct(33)
        En3(2,3,3)=fct(34)
        En3(2,3,4)=fct(35)
        En3(2,4,4)=fct(36)
        En3(3,3,3)=fct(37)
        En3(3,3,4)=fct(38)
        En3(3,4,4)=fct(39)
        En3(4,4,4)=fct(40)
        if (rank.eq.3) goto 511
        En4(0,0,0,0)=fct(41)
        En4(0,0,1,1)=fct(42)
        En4(0,0,1,2)=fct(43)
        En4(0,0,1,3)=fct(44)
        En4(0,0,1,4)=fct(45)
        En4(0,0,2,2)=fct(46)
        En4(0,0,2,3)=fct(47)
        En4(0,0,2,4)=fct(48)
        En4(0,0,3,3)=fct(49)
        En4(0,0,3,4)=fct(50)
        En4(0,0,4,4)=fct(51)
        En4(1,1,1,1)=fct(52)
        En4(1,1,1,2)=fct(53)
        En4(1,1,1,3)=fct(54)
        En4(1,1,1,4)=fct(55)
        En4(1,1,2,2)=fct(56)
        En4(1,1,2,3)=fct(57)
        En4(1,1,2,4)=fct(58)
        En4(1,1,3,3)=fct(59)
        En4(1,1,3,4)=fct(60)
        En4(1,1,4,4)=fct(61)
        En4(1,2,2,2)=fct(62)
        En4(1,2,2,3)=fct(63)
        En4(1,2,2,4)=fct(64)
        En4(1,2,3,3)=fct(65)
        En4(1,2,3,4)=fct(66)
        En4(1,2,4,4)=fct(67)
        En4(1,3,3,3)=fct(68)
        En4(1,3,3,4)=fct(69)
        En4(1,3,4,4)=fct(70)
        En4(1,4,4,4)=fct(71)
        En4(2,2,2,2)=fct(72)
        En4(2,2,2,3)=fct(73)
        En4(2,2,2,4)=fct(74)
        En4(2,2,3,3)=fct(75)
        En4(2,2,3,4)=fct(76)
        En4(2,2,4,4)=fct(77)
        En4(2,3,3,3)=fct(78)
        En4(2,3,3,4)=fct(79)
        En4(2,3,4,4)=fct(80)
        En4(2,4,4,4)=fct(81)
        En4(3,3,3,3)=fct(82)
        En4(3,3,3,4)=fct(83)
        En4(3,3,4,4)=fct(84)
        En4(3,4,4,4)=fct(85)
        En4(4,4,4,4)=fct(86)

 511    continue
        if (sym.eq.1) then
        do i1=1,4
          do i2=i1+1,4
            En2(i2,i1)=En2(i1,i2) 
          enddo
        enddo
        if (rank.eq.2) goto 522
        do i1=1,4
          do i2=i1,4
            do i3=i2,4
              En3(i1,i3,i2)=En3(i1,i2,i3) 
              En3(i2,i1,i3)=En3(i1,i2,i3) 
              En3(i2,i3,i1)=En3(i1,i2,i3) 
              En3(i3,i1,i2)=En3(i1,i2,i3) 
              En3(i3,i2,i1)=En3(i1,i2,i3) 
            enddo
          enddo
        enddo
        if (rank.eq.3) goto 522
        do i1=1,4
          do i2=i1+1,4
            En4(0,0,i2,i1)=En4(0,0,i1,i2) 
          enddo
        enddo
        do i1=1,4
          do i2=i1,4
            do i3=i2,4
              do i4=i3,4
                En4(i1,i2,i4,i3)=En4(i1,i2,i3,i4) 
                En4(i1,i3,i2,i4)=En4(i1,i2,i3,i4) 
                En4(i1,i3,i4,i2)=En4(i1,i2,i3,i4) 
                En4(i1,i4,i2,i3)=En4(i1,i2,i3,i4) 
                En4(i1,i4,i3,i2)=En4(i1,i2,i3,i4) 
                En4(i2,i1,i3,i4)=En4(i1,i2,i3,i4) 
                En4(i2,i1,i4,i3)=En4(i1,i2,i3,i4) 
                En4(i2,i3,i1,i4)=En4(i1,i2,i3,i4) 
                En4(i2,i3,i4,i1)=En4(i1,i2,i3,i4) 
                En4(i2,i4,i1,i3)=En4(i1,i2,i3,i4) 
                En4(i2,i4,i3,i1)=En4(i1,i2,i3,i4) 
                En4(i3,i1,i2,i4)=En4(i1,i2,i3,i4) 
                En4(i3,i1,i4,i2)=En4(i1,i2,i3,i4) 
                En4(i3,i2,i1,i4)=En4(i1,i2,i3,i4) 
                En4(i3,i2,i4,i1)=En4(i1,i2,i3,i4) 
                En4(i3,i4,i1,i2)=En4(i1,i2,i3,i4) 
                En4(i3,i4,i2,i1)=En4(i1,i2,i3,i4) 
                En4(i4,i1,i2,i3)=En4(i1,i2,i3,i4) 
                En4(i4,i1,i3,i2)=En4(i1,i2,i3,i4) 
                En4(i4,i2,i1,i3)=En4(i1,i2,i3,i4) 
                En4(i4,i2,i3,i1)=En4(i1,i2,i3,i4) 
                En4(i4,i3,i1,i2)=En4(i1,i2,i3,i4) 
                En4(i4,i3,i2,i1)=En4(i1,i2,i3,i4) 
              enddo
            enddo
          enddo
        enddo
      endif
      goto 522
      end if

      mm02 = elimcminf2(m02)
      mm12 = elimcminf2(m12)
      mm22 = elimcminf2(m22)
      mm32 = elimcminf2(m32)
      mm42 = elimcminf2(m42)
      q10  = elimcminf2(p10)
      q21  = elimcminf2(p21)
      q32  = elimcminf2(p32)
      q43  = elimcminf2(p43)
      q40  = elimcminf2(p40)
      q31  = elimcminf2(p31)
      q42  = elimcminf2(p42)
      q30  = elimcminf2(p30)
      q41  = elimcminf2(p41)
      q20  = elimcminf2(p20)
  
      mx(0,0) = 2d0*mm02
      mx(1,0) = q10 - mm12 + mm02
      mx(2,0) = q20 - mm22 + mm02
      mx(3,0) = q30 - mm32 + mm02
      mx(4,0) = q40 - mm42 + mm02
      mx(0,1) = mx(1,0)
      mx(0,2) = mx(2,0)
      mx(0,3) = mx(3,0)
      mx(0,4) = mx(4,0)

      mx(1,1) = 2d0*q10
      mx(2,2) = 2d0*q20
      mx(3,3) = 2d0*q30
      mx(4,4) = 2d0*q40
      mx(1,2) = q10+q20-q21
      mx(1,3) = q10+q30-q31
      mx(1,4) = q10+q40-q41
      mx(2,3) = q20+q30-q32
      mx(2,4) = q20+q40-q42
      mx(3,4) = q30+q40-q43
      mx(2,1) = mx(1,2)
      mx(3,1) = mx(1,3)
      mx(4,1) = mx(1,4)
      mx(3,2) = mx(2,3)
      mx(4,2) = mx(2,4)
      mx(4,3) = mx(3,4)

      call chinv(5,mx,mxinv)

c      do i=1,4
c      do j=1,4
c        z(i,j) = mx(i,j)
c      end do
c      end do
 
c      call chinv(4,z,zinv)

c      write(55,*) 'x01',mxinv(0,1)
c      write(55,*) 'x02',mxinv(0,2)
c      write(55,*) 'x03',mxinv(0,3)
c      write(55,*) 'x04',mxinv(0,4)
c      write(55,*) 'x00',mxinv(0,0)
c      write(55,*) 'sum',
c     &    -mxinv(0,0)-mxinv(0,1)-mxinv(0,2)-mxinv(0,3)-mxinv(0,4)


      rankD = max(rank-1,0)
c      write(*,*) 'Ep1 '
      call cDp12345(p21,p32,p43,p41,p31,p42,m12,m22,m32,m42,
     &            D0w0,D1w0,D2w0,D3w0,D4w0,D5w0,rankD,1)
c      write(*,*) 'Ep2 '
      call cDp12345(p20,p32,p43,p40,p30,p42,m02,m22,m32,m42,
     &            D0w1,D1w1,D2w1,D3w1,D4w1,D5w1,rankD,1)
c      write(*,*) 'Ep3 '
      call cDp12345(p10,p31,p43,p40,p30,p41,m02,m12,m32,m42,
     &            D0w2,D1w2,D2w2,D3w2,D4w2,D5w2,rankD,1)
c      write(*,*) 'Ep4 '
      call cDp12345(p10,p21,p42,p40,p20,p41,m02,m12,m22,m42,
     &            D0w3,D1w3,D2w3,D3w3,D4w3,D5w3,rankD,1)
c      write(*,*) 'Ep5 '
      call cDp12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &            D0w4,D1w4,D2w4,D3w4,D4w4,D5w4,rankD,1)

      D0m(0) = D0w0
      D0m(1) = D0w1
      D0m(2) = D0w2
      D0m(3) = D0w3
      D0m(4) = D0w4

c>      write(*,*) 'Ep D ',D0w0,D0w1,D0w2,D0w3,D0w4
c>      write(*,*) 'Ep e ',eta

      En0 = - mxinv(0,0)*D0m(0)
      do 5 i=1,4
         En0 = En0 + mxinv(i,0)*(D0m(i)-D0m(0))
  5   continue

C>      write(*,*) 'Ep: dety ',dety
C>      write(*,*) 'Ep: eta0 ',eta(0)
C>      write(*,*) 'Ep: eta1 ',eta(1)
C>      write(*,*) 'Ep: eta2 ',eta(2)
C>      write(*,*) 'Ep: eta3 ',eta(3)
C>      write(*,*) 'Ep: eta4 ',eta(4)
C>      write(*,*) 'Ep: D0w0 ',D0w0
C>      write(*,*) 'Ep: D0w1 ',D0w1
C>      write(*,*) 'Ep: D0w2 ',D0w2
C>      write(*,*) 'Ep: D0w3 ',D0w3
C>      write(*,*) 'Ep: D0w4 ',D0w4
C>      write(*,*) 'Ep: E0   ',En0
C>      write(*,*) 'Ep '

      fct(1)=En0

      if(rank.gt.0) then
 
      do 15 j=1,4
      En1(j) = -mxinv(j,0)*D0m(0)
      do 15 i=1,4
         En1(j) = En1(j) + mxinv(j,i)*(D0m(i)- D0m(0))
C> 16   continue
C>      write(*,*) 'Ep: D1m0 ',D1m(0,j)
C>      write(*,*) 'Ep: D1m1 ',D1m(1,j)
C>      write(*,*) 'Ep: D1m2 ',D1m(2,j)
C>      write(*,*) 'Ep: D1m3 ',D1m(3,j)
C>      write(*,*) 'Ep: D1m4 ',D1m(4,j)
C>      write(*,*) 'Ep: zeta ',zeta(1,j)
C>      write(*,*) 'Ep: zeta ',zeta(2,j)
C>      write(*,*) 'Ep: zeta ',zeta(3,j)
C>      write(*,*) 'Ep: zeta ',zeta(4,j)
C>      write(*,*) 'Ep: D2m1 ',D2m00(1)-D2m00(0)
C>      write(*,*) 'Ep: D2m2 ',D2m00(2)-D2m00(0)
C>      write(*,*) 'Ep: D2m3 ',D2m00(3)-D2m00(0)
C>      write(*,*) 'Ep: D2m4 ',D2m00(4)-D2m00(0)
C>      write(*,*) 'Ep: E1   ',En1(j)
C>      write(*,*) 'Ep '
 15   continue
        
      fct(2)=En1(1)
      fct(3)=En1(2)
      fct(4)=En1(3)
      fct(5)=En1(4)

      if(rank.gt.1) then

      D1m(0,1) = -D0w0 - D1w0(1) - D1w0(2) - D1w0(3)
      D1m(0,2) = D1w0(1)
      D1m(0,3) = D1w0(2)
      D1m(0,4) = D1w0(3)
      D1m(1,2) = D1w1(1)
      D1m(1,3) = D1w1(2)
      D1m(1,4) = D1w1(3)
      D1m(2,1) = D1w2(1)
      D1m(2,3) = D1w2(2)
      D1m(2,4) = D1w2(3)
      D1m(3,1) = D1w3(1)
      D1m(3,2) = D1w3(2)
      D1m(3,4) = D1w3(3)
      D1m(4,1) = D1w4(1)
      D1m(4,2) = D1w4(2)
      D1m(4,3) = D1w4(3)

      D2m00(0) = D2w0(0,0)
      D2m00(1) = D2w1(0,0)
      D2m00(2) = D2w2(0,0)
      D2m00(3) = D2w3(0,0)
      D2m00(4) = D2w4(0,0)

      En2(0,0) = 0d0
c      En2(0,0) = - mxinv(0,0)*D2m00(0)
      do 25 i=1,4
        En2(0,0) = En2(0,0) + mxinv(i,0)*(D2m00(i)-D2m00(0))
 25   continue

      do 27 j=1,4
      do 27 k=j,4
        En2(j,k) = -mxinv(j,0)*D1m(0,k)
     &             -mxinv(k,0)*D1m(0,j)
c     &      - 2d0*(mxinv(0,j)*mxinv(0,k) - mxinv(0,0)*mxinv(j,k))
c     &      * D2m00(0)
c     &      - 2d0*(mxinv(0,k)*mxinv(0,j) - mxinv(0,0)*mxinv(j,k))
c     &      * D2m00(0)
      do 26 i=1,4
        En2(j,k) = En2(j,k) + mxinv(i,j)*(D1m(i,k)- D1m(0,k))
     &      + 2d0*(mxinv(0,j)*mxinv(i,k) - mxinv(0,i)*mxinv(j,k))
     &      * (D2m00(i)- D2m00(0))
     &                      + mxinv(i,k)*(D1m(i,j)- D1m(0,j))
     &      + 2d0*(mxinv(0,k)*mxinv(i,j) - mxinv(0,i)*mxinv(j,k))
     &      * (D2m00(i)- D2m00(0))
 26   continue
        En2(j,k) = En2(j,k)/2d0
        En2(k,j) = En2(j,k)
 27   continue

      fct(6) =En2(0,0)
      fct(7) =En2(1,1)
      fct(8) =En2(1,2)
      fct(9) =En2(1,3)
      fct(10)=En2(1,4)
      fct(11)=En2(2,2)
      fct(12)=En2(2,3)
      fct(13)=En2(2,4)
      fct(14)=En2(3,3)
      fct(15)=En2(3,4)
      fct(16)=En2(4,4)

      if(rank.gt.2) then

c      D2m(0,1,1) = D0w0 + 2*(D1w0(1) + D1w0(2) + D1w0(3))
c     &           + D2w0(1,1) + D2w0(2,2) + D2w0(3,3)  
c     &           + 2*(D2w0(1,2) + D2w0(1,3) + D2w0(2,3))  
      D2m(0,1,2) = -D1w0(1) - D2w0(1,1) - D2w0(1,2) - D2w0(1,3) 
      D2m(0,1,3) = -D1w0(2) - D2w0(1,2) - D2w0(2,2) - D2w0(2,3) 
      D2m(0,1,4) = -D1w0(3) - D2w0(1,3) - D2w0(2,3) - D2w0(3,3) 
      D2m(0,2,2) = D2w0(1,1) 
      D2m(0,2,3) = D2w0(1,2) 
      D2m(0,2,4) = D2w0(1,3) 
      D2m(0,3,3) = D2w0(2,2) 
      D2m(0,3,4) = D2w0(2,3) 
      D2m(0,4,4) = D2w0(3,3) 
      D2m(0,1,1) = -D1m(0,1)-D2m(0,1,2)-D2m(0,1,3)-D2m(0,1,4)


      D2m(1,2,2) = D2w1(1,1) 
      D2m(1,2,3) = D2w1(1,2) 
      D2m(1,2,4) = D2w1(1,3) 
      D2m(1,3,3) = D2w1(2,2) 
      D2m(1,3,4) = D2w1(2,3) 
      D2m(1,4,4) = D2w1(3,3) 
      D2m(2,1,1) = D2w2(1,1) 
      D2m(2,1,3) = D2w2(1,2) 
      D2m(2,1,4) = D2w2(1,3) 
      D2m(2,3,3) = D2w2(2,2) 
      D2m(2,3,4) = D2w2(2,3) 
      D2m(2,4,4) = D2w2(3,3) 
      D2m(3,1,1) = D2w3(1,1) 
      D2m(3,1,2) = D2w3(1,2) 
      D2m(3,1,4) = D2w3(1,3) 
      D2m(3,2,2) = D2w3(2,2) 
      D2m(3,2,4) = D2w3(2,3) 
      D2m(3,4,4) = D2w3(3,3) 
      D2m(4,1,1) = D2w4(1,1) 
      D2m(4,1,2) = D2w4(1,2) 
      D2m(4,1,3) = D2w4(1,3) 
      D2m(4,2,2) = D2w4(2,2) 
      D2m(4,2,3) = D2w4(2,3) 
      D2m(4,3,3) = D2w4(3,3) 

      D3m00(0,1)   =-D2w0(0,0) - D3w0(0,0,1) - D3w0(0,0,2) - D3w0(0,0,3) 
      D3m00(0,2)   = D3w0(0,0,1)
      D3m00(0,3)   = D3w0(0,0,2)
      D3m00(0,4)   = D3w0(0,0,3)
      D3m00(1,2)   = D3w1(0,0,1)
      D3m00(1,3)   = D3w1(0,0,2)
      D3m00(1,4)   = D3w1(0,0,3)
      D3m00(2,1)   = D3w2(0,0,1) 
      D3m00(2,3)   = D3w2(0,0,2)
      D3m00(2,4)   = D3w2(0,0,3)
      D3m00(3,1)   = D3w3(0,0,1) 
      D3m00(3,2)   = D3w3(0,0,2)
      D3m00(3,4)   = D3w3(0,0,3)
      D3m00(4,1)   = D3w4(0,0,1) 
      D3m00(4,2)   = D3w4(0,0,2)
      D3m00(4,3)   = D3w4(0,0,3)


      do 35 j=1,4
        En3(0,0,j) = - mxinv(j,0)*D2m00(0)
c     &               -2d0*mxinv(0,0)*D3m00(0,j)
      do 36 i=1,4
        En3(0,0,j) = En3(0,0,j) + mxinv(j,i)*(D2m00(i)-D2m00(0))
        En3(0,0,j) = En3(0,0,j) + 2d0*mxinv(0,i)*(D3m00(i,j)-D3m00(0,j))
 36   continue
        En3(0,0,j) = En3(0,0,j)/3d0
 35   continue

c>      do i=0,4
c>      do k=1,4
c>      do l=k,4
c>       D2m(i,l,k) = D2m(i,k,l)
c>      end do
c>      end do
c>      end do

      do 37 j=1,4
      do 37 k=j,4
      do 37 l=k,4        
        En3(j,k,l) = -mxinv(j,0)*D2m(0,k,l)
     &               -mxinv(k,0)*D2m(0,j,l)
     &               -mxinv(l,0)*D2m(0,j,k)
c>     &      - 4d0*(mxinv(0,j)*mxinv(0,k) - mxinv(0,0)*mxinv(j,k))
c>     &      * D3m00(0,l)
c>     &      - 4d0*(mxinv(0,j)*mxinv(0,l) - mxinv(0,0)*mxinv(j,l))
c>     &      * D3m00(0,k)
c>     &      - 4d0*(mxinv(0,k)*mxinv(0,l) - mxinv(0,0)*mxinv(k,l))
c>     &      * D3m00(0,j)
      do 38 i=1,4
        En3(j,k,l) = En3(j,k,l) + mxinv(i,j)*(D2m(i,k,l)- D2m(0,k,l))
     &      + 2d0*(mxinv(0,j)*mxinv(i,k) - mxinv(0,i)*mxinv(j,k))
     &           * (D3m00(i,l)- D3m00(0,l))
     &      + 2d0*(mxinv(0,j)*mxinv(i,l) - mxinv(0,i)*mxinv(j,l))
     &           * (D3m00(i,k)- D3m00(0,k))
        En3(j,k,l) = En3(j,k,l) + mxinv(i,k)*(D2m(i,j,l)- D2m(0,j,l))
     &      + 2d0*(mxinv(0,k)*mxinv(i,j) - mxinv(0,i)*mxinv(j,k))
     &           * (D3m00(i,l)- D3m00(0,l))
     &      + 2d0*(mxinv(0,k)*mxinv(i,l) - mxinv(0,i)*mxinv(k,l))
     &           * (D3m00(i,j)- D3m00(0,j))
        En3(j,k,l) = En3(j,k,l) + mxinv(i,l)*(D2m(i,j,k)- D2m(0,j,k))
     &      + 2d0*(mxinv(0,l)*mxinv(i,k) - mxinv(0,i)*mxinv(l,k))
     &           * (D3m00(i,j)- D3m00(0,j))
     &      + 2d0*(mxinv(0,l)*mxinv(i,j) - mxinv(0,i)*mxinv(j,l))
     &           * (D3m00(i,k)- D3m00(0,k))

c        write(55,*) i,k,l,En3(j,k,l),En3(j,l,k)
c        write(55,*) i,k,l,D2m(i,k,l),D2m(i,l,k)
 38   continue
        En3(j,k,l) = En3(j,k,l)/3d0   
        if (sym.eq.1) then
          En3(l,j,k) = En3(j,k,l)
          En3(k,l,j) = En3(j,k,l)
          En3(j,l,k) = En3(j,k,l)
          En3(l,k,j) = En3(j,k,l)
          En3(k,j,l) = En3(j,k,l)
        end if
 37   continue
        
      fct(17)=En3(0,0,1)
      fct(18)=En3(0,0,2)
      fct(19)=En3(0,0,3)
      fct(20)=En3(0,0,4)
      fct(21)=En3(1,1,1)
      fct(22)=En3(1,1,2)
      fct(23)=En3(1,1,3)
      fct(24)=En3(1,1,4)
      fct(25)=En3(1,2,2)
      fct(26)=En3(1,2,3)
      fct(27)=En3(1,2,4)
      fct(28)=En3(1,3,3)
      fct(29)=En3(1,3,4)
      fct(30)=En3(1,4,4)
      fct(31)=En3(2,2,2)
      fct(32)=En3(2,2,3)
      fct(33)=En3(2,2,4)
      fct(34)=En3(2,3,3)
      fct(35)=En3(2,3,4)
      fct(36)=En3(2,4,4)
      fct(37)=En3(3,3,3)
      fct(38)=En3(3,3,4)
      fct(39)=En3(3,4,4)
      fct(40)=En3(4,4,4)
      
      if(rank.gt.3) then

      do 31 i1=1,3
      do 31 i2=i1,3
          D4m00(0,i1+1,i2+1) = D4w0(0,0,i1,i2) 
      do 31 i3=i2,3
          D3m(0,i1+1,i2+1,i3+1) = D3w0(i1,i2,i3)
 31   continue    

      D3m(0,1,2,2) = -D2m(0,2,2)-D3m(0,2,2,2)-D3m(0,2,2,3)-D3m(0,2,2,4) 
      D3m(0,1,2,3) = -D2m(0,2,3)-D3m(0,2,2,3)-D3m(0,2,3,3)-D3m(0,2,3,4) 
      D3m(0,1,2,4) = -D2m(0,2,4)-D3m(0,2,2,4)-D3m(0,2,3,4)-D3m(0,2,4,4) 
      D3m(0,1,3,3) = -D2m(0,3,3)-D3m(0,2,3,3)-D3m(0,3,3,3)-D3m(0,3,3,4) 
      D3m(0,1,3,4) = -D2m(0,3,4)-D3m(0,2,3,4)-D3m(0,3,3,4)-D3m(0,3,4,4) 
      D3m(0,1,4,4) = -D2m(0,4,4)-D3m(0,2,4,4)-D3m(0,3,4,4)-D3m(0,4,4,4) 
      D3m(0,1,1,2) = -D2m(0,1,2)-D3m(0,1,2,2)-D3m(0,1,2,3)-D3m(0,1,2,4) 
      D3m(0,1,1,3) = -D2m(0,1,3)-D3m(0,1,2,3)-D3m(0,1,3,3)-D3m(0,1,3,4) 
      D3m(0,1,1,4) = -D2m(0,1,4)-D3m(0,1,2,4)-D3m(0,1,3,4)-D3m(0,1,4,4) 
      D3m(0,1,1,1) = -D2m(0,1,1)-D3m(0,1,1,2)-D3m(0,1,1,3)-D3m(0,1,1,4) 

      D4m00(0,1,2) = -D3m00(0,2)-D4m00(0,2,2)-D4m00(0,2,3)-D4m00(0,2,4)
      D4m00(0,1,3) = -D3m00(0,3)-D4m00(0,2,3)-D4m00(0,3,3)-D4m00(0,3,4)
      D4m00(0,1,4) = -D3m00(0,4)-D4m00(0,2,4)-D4m00(0,3,4)-D4m00(0,4,4)
      D4m00(0,1,1) = -D3m00(0,1)-D4m00(0,1,2)-D4m00(0,1,3)-D4m00(0,1,4)

      do 33 i1=1,3
      do 33 i2=i1,3
          D4m00(1,i1+1,i2+1) = D4w1(0,0,i1,i2) 
          D4m00(4,i1,i2) = D4w4(0,0,i1,i2)
      do 33 i3=i2,3
          D3m(1,i1+1,i2+1,i3+1) = D3w1(i1,i2,i3)
          D3m(4,i1,i2,i3) = D3w4(i1,i2,i3)
 33   continue    
                     
      D3m(2,1,1,1) = D3w2(1,1,1)
      D3m(2,1,1,3) = D3w2(1,1,2)
      D3m(2,1,1,4) = D3w2(1,1,3)
      D3m(2,1,3,3) = D3w2(1,2,2)
      D3m(2,1,3,4) = D3w2(1,2,3)
      D3m(2,1,4,4) = D3w2(1,3,3)
      D3m(2,3,3,3) = D3w2(2,2,2)
      D3m(2,3,3,4) = D3w2(2,2,3)
      D3m(2,3,4,4) = D3w2(2,3,3)
      D3m(2,4,4,4) = D3w2(3,3,3)

      D3m(3,1,1,1) = D3w3(1,1,1)
      D3m(3,1,1,2) = D3w3(1,1,2)
      D3m(3,1,1,4) = D3w3(1,1,3)
      D3m(3,1,2,2) = D3w3(1,2,2)
      D3m(3,1,2,4) = D3w3(1,2,3)
      D3m(3,1,4,4) = D3w3(1,3,3)
      D3m(3,2,2,2) = D3w3(2,2,2)
      D3m(3,2,2,4) = D3w3(2,2,3)
      D3m(3,2,4,4) = D3w3(2,3,3)
      D3m(3,4,4,4) = D3w3(3,3,3)
 
      D4m00(2,1,1) = D4w2(0,0,1,1) 
      D4m00(2,1,3) = D4w2(0,0,1,2) 
      D4m00(2,1,4) = D4w2(0,0,1,3) 
      D4m00(2,3,3) = D4w2(0,0,2,2) 
      D4m00(2,3,4) = D4w2(0,0,2,3) 
      D4m00(2,4,4) = D4w2(0,0,3,3) 

      D4m00(3,1,1) = D4w3(0,0,1,1) 
      D4m00(3,1,2) = D4w3(0,0,1,2) 
      D4m00(3,1,4) = D4w3(0,0,1,3) 
      D4m00(3,2,2) = D4w3(0,0,2,2) 
      D4m00(3,2,4) = D4w3(0,0,2,3) 
      D4m00(3,4,4) = D4w3(0,0,3,3) 

      D4m0000(0) = D4w0(0,0,0,0) 
      D4m0000(1) = D4w1(0,0,0,0) 
      D4m0000(2) = D4w2(0,0,0,0) 
      D4m0000(3) = D4w3(0,0,0,0) 
      D4m0000(4) = D4w4(0,0,0,0) 

      En4(0,0,0,0) = 0d0
c      En4(0,0,0,0) = - mxinv(0,0)*(D4m0000(0)+1d0/48d0)
      do 44 i=1,4
        En4(0,0,0,0) = En4(0,0,0,0) + mxinv(i,0)*(D4m0000(i)-D4m0000(0))
 44   continue

c>      do i=0,4
c>      do k=1,4
c>      do l=k+1,4
c>        D4m00(i,l,k) = D4m00(i,k,l)
c>      end do
c>      end do
c>      end do

      do 45 j=1,4
      do 45 k=j,4
        En4(0,0,j,k) = 
     &      - mxinv(j,0)*D3m00(0,k)
     &      - mxinv(k,0)*D3m00(0,j)
c     &      - 2d0*(mxinv(0,j)*mxinv(0,k) - mxinv(0,0)*mxinv(j,k)
c     &           + mxinv(0,k)*mxinv(0,j) - mxinv(0,0)*mxinv(j,k))
c     &           * (D4m0000(0)+1d0/48d0)
c     &      - 2d0*mxinv(0,0)*D4m00(0,j,k)
      do 46 i=1,4
        En4(0,0,j,k) = En4(0,0,j,k) 
     &      + mxinv(j,i)*(D3m00(i,k)-D3m00(0,k))
     &      + mxinv(k,i)*(D3m00(i,j)-D3m00(0,j))
     &      + 2d0*(mxinv(0,j)*mxinv(i,k) - mxinv(0,i)*mxinv(j,k)
     &           + mxinv(0,k)*mxinv(i,j) - mxinv(0,i)*mxinv(j,k))
     &           * (D4m0000(i) - D4m0000(0))
     &      + 2d0*mxinv(0,i)*(D4m00(i,j,k)-D4m00(0,j,k))
 46   continue
        En4(0,0,j,k) = En4(0,0,j,k)/4d0
        En4(0,0,k,j) = En4(0,0,j,k)
 45   continue

c>      do i=0,4
c>      do k=1,4
c>      do l=k,4
c>      do m=l,4
c>        D3m(i,l,k,m) = D3m(i,k,l,m)
c>        D3m(i,k,m,l) = D3m(i,k,l,m)
c>        D3m(i,m,k,l) = D3m(i,k,l,m)
c>        D3m(i,l,m,k) = D3m(i,k,l,m)
c>        D3m(i,m,l,k) = D3m(i,k,l,m)
c>      end do
c>      end do
c>      end do
c>      end do

      do 47 j=1,4
      do 47 k=j,4
      do 47 l=k,4        
      do 47 m=l,4        
        En4(j,k,l,m) = -mxinv(j,0)*D3m(0,k,l,m)
     &               -mxinv(k,0)*D3m(0,j,l,m)
     &               -mxinv(l,0)*D3m(0,j,k,m)
     &               -mxinv(m,0)*D3m(0,j,k,l)
c>     &      - 4d0*(mxinv(0,j)*mxinv(0,k) - mxinv(0,0)*mxinv(j,k))
c>     &      * D4m00(0,l,m)
c>     &      - 4d0*(mxinv(0,j)*mxinv(0,l) - mxinv(0,0)*mxinv(j,l))
c>     &      * D4m00(0,k,m)
c>     &      - 4d0*(mxinv(0,k)*mxinv(0,l) - mxinv(0,0)*mxinv(k,l))
c>     &      * D4m00(0,j,m)
c>     &      - 4d0*(mxinv(0,j)*mxinv(0,m) - mxinv(0,0)*mxinv(j,m))
c>     &      * D4m00(0,k,l)
c>     &      - 4d0*(mxinv(0,k)*mxinv(0,m) - mxinv(0,0)*mxinv(k,m))
c>     &      * D4m00(0,j,l)
c>     &      - 4d0*(mxinv(0,l)*mxinv(0,m) - mxinv(0,0)*mxinv(l,m))
c>     &      * D4m00(0,j,k)
      do 48 i=1,4
        En4(j,k,l,m) = En4(j,k,l,m) 
     &      + mxinv(i,j)*(D3m(i,k,l,m)- D3m(0,k,l,m))
     &      + mxinv(i,k)*(D3m(i,j,l,m)- D3m(0,j,l,m))
     &      + mxinv(i,l)*(D3m(i,j,k,m)- D3m(0,j,k,m))
     &      + mxinv(i,m)*(D3m(i,j,k,l)- D3m(0,j,k,l))
     &      + 2d0*(mxinv(0,j)*mxinv(i,k) - mxinv(0,i)*mxinv(j,k)
     &           + mxinv(0,k)*mxinv(i,j) - mxinv(0,i)*mxinv(j,k))
     &           * (D4m00(i,l,m)- D4m00(0,l,m))
     &      + 2d0*(mxinv(0,j)*mxinv(i,l) - mxinv(0,i)*mxinv(j,l)
     &           + mxinv(0,l)*mxinv(i,j) - mxinv(0,i)*mxinv(j,l))
     &           * (D4m00(i,k,m)- D4m00(0,k,m))
     &      + 2d0*(mxinv(0,k)*mxinv(i,l) - mxinv(0,i)*mxinv(k,l)
     &           + mxinv(0,l)*mxinv(i,k) - mxinv(0,i)*mxinv(k,l))
     &           * (D4m00(i,j,m)- D4m00(0,j,m))
     &      + 2d0*(mxinv(0,j)*mxinv(i,m) - mxinv(0,i)*mxinv(j,m)
     &           + mxinv(0,m)*mxinv(i,j) - mxinv(0,i)*mxinv(j,m))
     &           * (D4m00(i,k,l)- D4m00(0,k,l))
     &      + 2d0*(mxinv(0,k)*mxinv(i,m) - mxinv(0,i)*mxinv(k,m)
     &           + mxinv(0,m)*mxinv(i,k) - mxinv(0,i)*mxinv(k,m))
     &           * (D4m00(i,j,l)- D4m00(0,j,l))
     &      + 2d0*(mxinv(0,l)*mxinv(i,m) - mxinv(0,i)*mxinv(l,m)
     &           + mxinv(0,m)*mxinv(i,l) - mxinv(0,i)*mxinv(l,m))
     &           * (D4m00(i,j,k)- D4m00(0,j,k))
 48   continue
        En4(j,k,l,m) = En4(j,k,l,m)/4d0   
        if (sym.eq.1) then
          En4(l,j,k,m) = En4(j,k,l,m)
          En4(k,l,j,m) = En4(j,k,l,m)
          En4(j,l,k,m) = En4(j,k,l,m)
          En4(l,k,j,m) = En4(j,k,l,m)
          En4(k,j,l,m) = En4(j,k,l,m)
          En4(j,k,m,l) = En4(j,k,l,m)
          En4(l,j,m,k) = En4(j,k,l,m)
          En4(k,l,m,j) = En4(j,k,l,m)
          En4(j,l,m,k) = En4(j,k,l,m)
          En4(l,k,m,j) = En4(j,k,l,m)
          En4(k,j,m,l) = En4(j,k,l,m)
          En4(j,m,k,l) = En4(j,k,l,m)
          En4(l,m,j,k) = En4(j,k,l,m)
          En4(k,m,l,j) = En4(j,k,l,m)
          En4(j,m,l,k) = En4(j,k,l,m)
          En4(l,m,k,j) = En4(j,k,l,m)
          En4(k,m,j,l) = En4(j,k,l,m)
          En4(m,j,k,l) = En4(j,k,l,m)
          En4(m,l,j,k) = En4(j,k,l,m)
          En4(m,k,l,j) = En4(j,k,l,m)
          En4(m,j,l,k) = En4(j,k,l,m)
          En4(m,l,k,j) = En4(j,k,l,m)
          En4(m,k,j,l) = En4(j,k,l,m)
        end if
 47   continue

      fct(41)=En4(0,0,0,0)
      fct(42)=En4(0,0,1,1)
      fct(43)=En4(0,0,1,2)
      fct(44)=En4(0,0,1,3)
      fct(45)=En4(0,0,1,4)
      fct(46)=En4(0,0,2,2)
      fct(47)=En4(0,0,2,3)
      fct(48)=En4(0,0,2,4)
      fct(49)=En4(0,0,3,3)
      fct(50)=En4(0,0,3,4)
      fct(51)=En4(0,0,4,4)
      fct(52)=En4(1,1,1,1)
      fct(53)=En4(1,1,1,2)
      fct(54)=En4(1,1,1,3)
      fct(55)=En4(1,1,1,4)
      fct(56)=En4(1,1,2,2)
      fct(57)=En4(1,1,2,3)
      fct(58)=En4(1,1,2,4)
      fct(59)=En4(1,1,3,3)
      fct(60)=En4(1,1,3,4)
      fct(61)=En4(1,1,4,4)
      fct(62)=En4(1,2,2,2)
      fct(63)=En4(1,2,2,3)
      fct(64)=En4(1,2,2,4)
      fct(65)=En4(1,2,3,3)
      fct(66)=En4(1,2,3,4)
      fct(67)=En4(1,2,4,4)
      fct(68)=En4(1,3,3,3)
      fct(69)=En4(1,3,3,4)
      fct(70)=En4(1,3,4,4)
      fct(71)=En4(1,4,4,4)
      fct(72)=En4(2,2,2,2)
      fct(73)=En4(2,2,2,3)
      fct(74)=En4(2,2,2,4)
      fct(75)=En4(2,2,3,3)
      fct(76)=En4(2,2,3,4)
      fct(77)=En4(2,2,4,4)
      fct(78)=En4(2,3,3,3)
      fct(79)=En4(2,3,3,4)
      fct(80)=En4(2,3,4,4)
      fct(81)=En4(2,4,4,4)
      fct(82)=En4(3,3,3,3)
      fct(83)=En4(3,3,3,4)
      fct(84)=En4(3,3,4,4)
      fct(85)=En4(3,4,4,4)
      fct(86)=En4(4,4,4,4)

      if (rank.gt.4) then
        write(*,*) 'rank > 4 not implemented in cEp1234'
        write(*,*) 'rank = ',rank
        stop

      do 41 i1=1,3
         D5m0000(0,i1+1) = D5w0(0,0,0,0,i1) 
      do 41 i2=1,3
      do 41 i3=1,3
         D5m00(0,i1+1,i2+1,i3+1) = D5w0(0,0,i1,i2,i3)
      do 41 i4=1,3
         D4m(0,i1+1,i2+1,i3+1,i4+1) = D4w0(i1,i2,i3,i4)
 41   continue             

      D4m(0,1,2,2,2) = -D3m(0,2,2,2)
     &                 -D4m(0,2,2,2,2)-D4m(0,2,2,2,3)-D4m(0,2,2,2,4) 
      D4m(0,1,2,2,3) = -D3m(0,2,2,3)
     &                 -D4m(0,2,2,2,3)-D4m(0,2,2,3,3)-D4m(0,2,2,3,4) 
      D4m(0,1,2,2,4) = -D3m(0,2,2,4)
     &                 -D4m(0,2,2,2,4)-D4m(0,2,2,3,4)-D4m(0,2,2,4,4) 
      D4m(0,1,2,3,3) = -D3m(0,2,3,3)
     &                 -D4m(0,2,2,3,3)-D4m(0,2,3,3,3)-D4m(0,2,3,3,4) 
      D4m(0,1,2,3,4) = -D3m(0,2,3,4)
     &                 -D4m(0,2,2,3,4)-D4m(0,2,3,3,4)-D4m(0,2,3,4,4) 
      D4m(0,1,2,4,4) = -D3m(0,2,4,4)
     &                 -D4m(0,2,2,4,4)-D4m(0,2,3,4,4)-D4m(0,2,4,4,4) 
      D4m(0,1,3,3,3) = -D3m(0,3,3,3)
     &                 -D4m(0,2,3,3,3)-D4m(0,3,3,3,3)-D4m(0,3,3,3,4) 
      D4m(0,1,3,3,4) = -D3m(0,3,3,4)
     &                 -D4m(0,2,3,3,4)-D4m(0,3,3,3,4)-D4m(0,3,3,4,4) 
      D4m(0,1,3,4,4) = -D3m(0,3,4,4)
     &                 -D4m(0,2,3,4,4)-D4m(0,3,3,4,4)-D4m(0,3,4,4,4) 
      D4m(0,1,4,4,4) = -D3m(0,4,4,4)
     &                 -D4m(0,2,4,4,4)-D4m(0,3,4,4,4)-D4m(0,4,4,4,4) 
      D4m(0,1,1,2,2) = -D3m(0,1,2,2)
     &                 -D4m(0,1,2,2,2)-D4m(0,1,2,2,3)-D4m(0,1,2,2,4) 
      D4m(0,1,1,2,3) = -D3m(0,1,2,3)
     &                 -D4m(0,1,2,2,3)-D4m(0,1,2,3,3)-D4m(0,1,2,3,4) 
      D4m(0,1,1,2,4) = -D3m(0,1,2,4)
     &                 -D4m(0,1,2,2,4)-D4m(0,1,2,3,4)-D4m(0,1,2,4,4) 
      D4m(0,1,1,3,3) = -D3m(0,1,3,3)
     &                 -D4m(0,1,2,3,3)-D4m(0,1,3,3,3)-D4m(0,1,3,3,4) 
      D4m(0,1,1,3,4) = -D3m(0,1,3,4)
     &                 -D4m(0,1,2,3,4)-D4m(0,1,3,3,4)-D4m(0,1,3,4,4) 
      D4m(0,1,1,4,4) = -D3m(0,1,4,4)
     &                 -D4m(0,1,2,4,4)-D4m(0,1,3,4,4)-D4m(0,1,4,4,4) 
      D4m(0,1,1,1,2) = -D3m(0,1,1,2)
     &                 -D4m(0,1,1,2,2)-D4m(0,1,1,2,3)-D4m(0,1,1,2,4) 
      D4m(0,1,1,1,3) = -D3m(0,1,1,3)
     &                 -D4m(0,1,1,2,3)-D4m(0,1,1,3,3)-D4m(0,1,1,3,4) 
      D4m(0,1,1,1,4) = -D3m(0,1,1,4)
     &                 -D4m(0,1,1,2,4)-D4m(0,1,1,3,4)-D4m(0,1,1,4,4) 
      D4m(0,1,1,1,1) = -D3m(0,1,1,1)
     &                 -D4m(0,1,1,1,2)-D4m(0,1,1,1,3)-D4m(0,1,1,1,4) 

      D5m00(0,1,2,2) = -D4m00(0,2,2)
     &                 -D5m00(0,2,2,2)-D5m00(0,2,2,3)-D5m00(0,2,2,4)
      D5m00(0,1,2,3) = -D4m00(0,2,3)
     &                 -D5m00(0,2,2,3)-D5m00(0,2,3,3)-D5m00(0,2,3,4)
      D5m00(0,1,2,4) = -D4m00(0,2,4)
     &                 -D5m00(0,2,2,4)-D5m00(0,2,3,4)-D5m00(0,2,4,4)
      D5m00(0,1,3,3) = -D4m00(0,3,3)
     &                 -D5m00(0,2,3,3)-D5m00(0,3,3,3)-D5m00(0,3,3,4)
      D5m00(0,1,3,4) = -D4m00(0,3,4)
     &                 -D5m00(0,2,3,4)-D5m00(0,3,3,4)-D5m00(0,3,4,4)
      D5m00(0,1,4,4) = -D4m00(0,4,4)
     &                 -D5m00(0,2,4,4)-D5m00(0,3,4,4)-D5m00(0,4,4,4)
      D5m00(0,1,1,2) = -D4m00(0,1,2)
     &                 -D5m00(0,1,2,2)-D5m00(0,1,2,3)-D5m00(0,1,2,4)
      D5m00(0,1,1,3) = -D4m00(0,1,3)
     &                 -D5m00(0,1,2,3)-D5m00(0,1,3,3)-D5m00(0,1,3,4)
      D5m00(0,1,1,4) = -D4m00(0,1,4)
     &                 -D5m00(0,1,2,4)-D5m00(0,1,3,4)-D5m00(0,1,4,4)
      D5m00(0,1,1,1) = -D4m00(0,1,1)
     &                 -D5m00(0,1,1,2)-D5m00(0,1,1,3)-D5m00(0,1,1,4)

      D5m0000(0,1)   = -D4w0(0,0,0,0) 
     &                 -D5w0(0,0,0,0,1)-D5w0(0,0,0,0,2)-D5w0(0,0,0,0,3)

      do 43 i1=1,3
          D5m0000(1,i1+1)=D5w1(0,0,0,0,i1)
          D5m0000(4,i1)=D5w4(0,0,0,0,i1)
      do 43 i2=i1,3
      do 43 i3=i2,3
          D5m00(1,i1+1,i2+1,i3+1) = D5w1(0,0,i1,i2,i3) 
          D5m00(4,i1,i2,i3) = D5w4(0,0,i1,i2,i3)
      do 43 i4=i2,3
          D4m(1,i1+1,i2+1,i3+1,i4+1) = D4w1(i1,i2,i3,i4)
          D4m(4,i1,i2,i3,i4) = D4w4(i1,i2,i3,i4)
 43   continue    

      D4m(2,1,1,1,1) = D4w2(1,1,1,1)
      D4m(2,1,1,1,3) = D4w2(1,1,1,2)
      D4m(2,1,1,1,4) = D4w2(1,1,1,3)
      D4m(2,1,1,3,3) = D4w2(1,1,2,2)
      D4m(2,1,1,3,4) = D4w2(1,1,2,3)
      D4m(2,1,1,4,4) = D4w2(1,1,3,3)
      D4m(2,1,3,3,3) = D4w2(1,2,2,2)
      D4m(2,1,3,3,4) = D4w2(1,2,2,3)
      D4m(2,1,3,4,4) = D4w2(1,2,3,3)
      D4m(2,1,4,4,4) = D4w2(1,3,3,3)
      D4m(2,3,3,3,3) = D4w2(2,2,2,2)
      D4m(2,3,3,3,4) = D4w2(2,2,2,3)
      D4m(2,3,3,4,4) = D4w2(2,2,3,3)
      D4m(2,3,4,4,4) = D4w2(2,3,3,3)
      D4m(2,4,4,4,4) = D4w2(3,3,3,3)
      D4m(3,1,1,1,1) = D4w3(1,1,1,1)
      D4m(3,1,1,1,2) = D4w3(1,1,1,2)
      D4m(3,1,1,1,4) = D4w3(1,1,1,3)
      D4m(3,1,1,2,2) = D4w3(1,1,2,2)
      D4m(3,1,1,2,4) = D4w3(1,1,2,3)
      D4m(3,1,1,4,4) = D4w3(1,1,3,3)
      D4m(3,1,2,2,2) = D4w3(1,2,2,2)
      D4m(3,1,2,2,4) = D4w3(1,2,2,3)
      D4m(3,1,2,4,4) = D4w3(1,2,3,3)
      D4m(3,1,4,4,4) = D4w3(1,3,3,3)
      D4m(3,2,2,2,2) = D4w3(2,2,2,2)
      D4m(3,2,2,2,4) = D4w3(2,2,2,3)
      D4m(3,2,2,4,4) = D4w3(2,2,3,3)
      D4m(3,2,4,4,4) = D4w3(2,3,3,3)
      D4m(3,4,4,4,4) = D4w3(3,3,3,3)

      D5m00(2,1,1,1) = D5w2(0,0,1,1,1) 
      D5m00(2,1,1,3) = D5w2(0,0,1,1,2) 
      D5m00(2,1,1,4) = D5w2(0,0,1,1,3) 
      D5m00(2,1,3,3) = D5w2(0,0,1,2,2) 
      D5m00(2,1,3,4) = D5w2(0,0,1,2,3) 
      D5m00(2,1,4,4) = D5w2(0,0,1,3,3) 
      D5m00(2,3,3,3) = D5w2(0,0,2,2,2) 
      D5m00(2,3,3,4) = D5w2(0,0,2,2,3) 
      D5m00(2,3,4,4) = D5w2(0,0,2,3,3) 
      D5m00(2,4,4,4) = D5w2(0,0,3,3,3) 

      D5m00(3,1,1,1) = D5w3(0,0,1,1,1) 
      D5m00(3,1,1,2) = D5w3(0,0,1,1,2) 
      D5m00(3,1,1,4) = D5w3(0,0,1,1,3) 
      D5m00(3,1,2,2) = D5w3(0,0,1,2,2) 
      D5m00(3,1,2,4) = D5w3(0,0,1,2,3) 
      D5m00(3,1,4,4) = D5w3(0,0,1,3,3) 
      D5m00(3,2,2,2) = D5w3(0,0,2,2,2) 
      D5m00(3,2,2,4) = D5w3(0,0,2,2,3) 
      D5m00(3,2,4,4) = D5w3(0,0,2,3,3) 
      D5m00(3,4,4,4) = D5w3(0,0,3,3,3) 

      D5m0000(2,1)=D5w2(0,0,0,0,1)
      D5m0000(2,3)=D5w2(0,0,0,0,2)
      D5m0000(2,4)=D5w2(0,0,0,0,3)
      D5m0000(3,1)=D5w3(0,0,0,0,1)
      D5m0000(3,2)=D5w3(0,0,0,0,2)
      D5m0000(3,4)=D5w3(0,0,0,0,3)

      end if

      end if
      end if
      end if
      end if

 521  continue

      call cachewrite(fct,86,type)

 522  continue

c>      if (n.eq.4.and.usecache.eq.1)  then
c>        write(87,*) 'nto = ',ncall(1),countE,sym
c>        write(87,*) 'E0 = ',En0
c>        write(87,*) 'E1 = ',En1
c>        write(87,*) 'E2 = ',En2
c>        write(87,*) 'E3 = ',En3
c>c        write(87,*) 'E4 = ',En4(0,0,0,0)
c>c        write(87,*) 'E4 = ',En4(0,0,1,1)
c>        do i1=1,4
c>        do i2=1,4
c>        do i3=1,4
c>c        write(87,133) i1,i2,i3,En3(i1,i2,i3)
c>        end do
c>        end do
c>        end do
c>      else if (n.eq.4.and.usecache.eq.0) then
c>        write(88,*) 'nto = ',countE,countE,sym
c>        write(88,*) 'E0 = ',En0
c>        write(88,*) 'E1 = ',En1
c>        write(88,*) 'E2 = ',En2
c>        write(88,*) 'E3 = ',En3
c>c        write(88,*) 'E4 = ',En4(0,0,0,0)
c>c        write(88,*) 'E4 = ',En4(0,0,1,1)
c>        do i1=1,4
c>        do i2=1,4
c>        do i3=1,4
c>c        write(88,133) i1,i2,i3,En3(i1,i2,i3)
c>        end do
c>        end do
c>        end do
c>      end if

      end 


***********************************************************************
      subroutine cFp1234(p10,p21,p32,p43,p54,p50
     &                 ,p20,p31,p42,p53,p40,p51,p30,p41,p52
     &                 ,m02,m12,m22,m32,m42,m52
     &                 ,Fn0,Fn1,Fn2,Fn3,Fn4
     &                 ,rank)
c      entry Fp1234E(p10,p21,p32,p43,p54,p50
c     &                 ,p20,p31,p42,p53,p40,p51,p30,p41,p52
c     &                 ,m02,m12,m22,m32,m42,m52
c     &                 ,Fn0,Fn1,Fn2,Fn3,Fn4,En0w0,En1w0,En2w0,En3w0
c     &                 ,rank)
***********************************************************************
*     6-point tensor coefficient functions F0,F1,F2,F3,F4             *
*     up to 4 integration momenta in numerator                        *
*     g(mu,nu) terms in decomposition of Fn...                        *
*     version that needs only E's of rank-1                           *
*                        and selects between different possibilities  *
*     coefficients calculated up to *rank* momenta in numerator       *
*     definition of integrals:                                        *
*     {F0,F1,F2,F3,F4}= 1/(i pi^2)*\int d^4 q {1,q,q.q,q.q.q,q.q.q.q} *
*        1/[(q**2 - m02)*[(q+k1)**2-m12]*[(q+k2)**2 - m22]*           *
*        [(q+k3)**2-m32]*[(q+k4)**2-m42]*[(q+k5)**2 - m52]]           *
*     arguments: p10=k1**2, p21=(k2-k1)**2, p32=(k3-k2)**2,           *
*          p43=(k4-k3)**2, p54=(k5-k4)**2, p50=k5**2, p20=k2**2,      * 
*          p31=(k3-k1)**2, p42=(k4-k2)**2, p53=(k5-k3)**2, p40=k4**2, *
*          p51=(k5-k1)**2, p30=k3**2, p41=(k4-k1)**2, p52=(k5-k2)**2  *
*     definition of tensor coefficients:                              *
*          F1 = \sum_{i=1}^5 Fn1(i) *ki                               *
*          F2 = \sum_{i,j=1}^5 Fn2(i,j) *ki*kj  + Fn2(0,0) * g        *
*          F3 = \sum_{i,j,k=1}^5 Fn3(i,j,k) *ki*kj*kk                 *
*              + \sum_{i=1}^5 Fn3(0,0,i) * (ki*g + 2*sym)             *
*          F4 = \sum_{i,j,k,l=1}^5 Fn4(i,j,k) *ki*kj*kk*kl            *
*              + \sum_{i,j=1}^5 Fn4(0,0,i,j) * (ki*kj*g + 5*sym)      *
*              +   Fn4(0,0,0,0) * (g*g + 2*sym)                       *
*          En1w0 =  \sum_{i=1}^4 En1w0(i) *( k(i+1)-k1 )              *
*     5-point functions resulting from 6-point functions with         *
*     i-th denominator cancelled: EnXwI                               *
***********************************************************************
*    08.06.05 Ansgar Denner    last changed  12.12.05                 *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p32,p43,p54,p50
     &          ,p20,p31,p42,p53,p40,p51,p30,p41,p52
      complex*16 m02,m12,m22,m32,m42,m52
      complex*16 q10,q21,q32,q43,q54,q50
     &          ,q20,q31,q42,q53,q40,q51,q30,q41,q52
      complex*16 mm02,mm12,mm22,mm32,mm42,mm52
      complex*16 mx(0:5,0:5),mxinv(0:5,0:5)
      complex*16 mx0k(1:5,1:5),mx0kinv(1:5,1:5)
      complex*16 chdet,det,newdet
      complex*16 En0w5,En1w5(4),En2w5(0:4,0:4),En3w5(0:4,0:4,0:4)
      complex*16 En0w4,En1w4(4),En2w4(0:4,0:4),En3w4(0:4,0:4,0:4)
      complex*16 En0w3,En1w3(4),En2w3(0:4,0:4),En3w3(0:4,0:4,0:4)
      complex*16 En0w2,En1w2(4),En2w2(0:4,0:4),En3w2(0:4,0:4,0:4)
      complex*16 En0w1,En1w1(4),En2w1(0:4,0:4),En3w1(0:4,0:4,0:4)
      complex*16 En0w0,En1w0(4),En2w0(0:4,0:4),En3w0(0:4,0:4,0:4)
      complex*16 En4w5(0:4,0:4,0:4,0:4)
      complex*16 En4w4(0:4,0:4,0:4,0:4)
      complex*16 En4w3(0:4,0:4,0:4,0:4)
      complex*16 En4w2(0:4,0:4,0:4,0:4)
      complex*16 En4w1(0:4,0:4,0:4,0:4)
      complex*16 En4w0(0:4,0:4,0:4,0:4)
      complex*16 E0m(0:5),E1m(0:5,5),E2m(0:5,1:5,1:5)
     &          ,E3m(0:5,1:5,1:5,1:5),E4m(0:5,1:5,1:5,1:5,1:5)
      complex*16 E2m00(0:5),E3m00(0:5,5),E4m00(0:5,1:5,1:5)
     &          ,E4m0000(0:5)
      complex*16 Fn0,Fn1(5),Fn2(0:5,0:5),Fn3(0:5,0:5,0:5)
      complex*16 Fn4(0:5,0:5,0:5,0:5)
      complex*16 elimcminf2
      integer    rank,order,ordm1,ordm2,rankE
      integer    i,j,i1,i2,i3,i4,k,l,m,kbest
      integer    ltest,sym

      common /ltest/  ltest                  
      common /sym/    sym

      data       order,ordm1,ordm2 /6,5,4/ 
      data       E1m /30*0d0/,E2m/150*0d0/,E3m/750*0d0/,E4m/3750*0d0/
      data       E2m00/6*0d0/,E3m00/30*0d0/,E4m00/150*0d0/
      data       E4m0000/6*0d0/

c      write(*,*) 'cFp in',     p10,p21,p32,p43,p54,p50
c     &                 ,p20,p31,p42,p53,p40,p51,p30,p41,p52
c     &                 ,m02,m12,m22,m32,m42,m52

      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cFp1234'
        write(*,*) 'rank = ',rank
      end if

      order = 6
      ordm1 = 5
      ordm2 = 4

      mm02 = elimcminf2(m02)
      mm12 = elimcminf2(m12)
      mm22 = elimcminf2(m22)
      mm32 = elimcminf2(m32)
      mm42 = elimcminf2(m42)
      mm52 = elimcminf2(m52)
      q10  = elimcminf2(p10)
      q21  = elimcminf2(p21)
      q32  = elimcminf2(p32)
      q43  = elimcminf2(p43)
      q54  = elimcminf2(p54)
      q50  = elimcminf2(p50)
      q20  = elimcminf2(p20)
      q31  = elimcminf2(p31)
      q42  = elimcminf2(p42)
      q53  = elimcminf2(p53)
      q40  = elimcminf2(p40)
      q51  = elimcminf2(p51)
      q30  = elimcminf2(p30)
      q41  = elimcminf2(p41)
      q52  = elimcminf2(p52)
 
      mx(0,0) = 2d0*mm02
      mx(1,0) = q10 - mm12 + mm02
      mx(2,0) = q20 - mm22 + mm02
      mx(3,0) = q30 - mm32 + mm02
      mx(4,0) = q40 - mm42 + mm02
      mx(5,0) = q50 - mm52 + mm02
      mx(0,1) = mx(1,0)
      mx(0,2) = mx(2,0)
      mx(0,3) = mx(3,0)
      mx(0,4) = mx(4,0)
      mx(0,5) = mx(5,0)

      mx(1,1) = 2d0*q10
      mx(2,2) = 2d0*q20
      mx(3,3) = 2d0*q30
      mx(4,4) = 2d0*q40
      mx(5,5) = 2d0*q50
      mx(1,2) = q10+q20-q21
      mx(1,3) = q10+q30-q31
      mx(1,4) = q10+q40-q41
      mx(1,5) = q10+q50-q51
      mx(2,3) = q20+q30-q32
      mx(2,4) = q20+q40-q42
      mx(2,5) = q20+q50-q52
      mx(3,4) = q30+q40-q43
      mx(3,5) = q30+q50-q53
      mx(4,5) = q40+q50-q54
      mx(2,1) = mx(1,2)
      mx(3,1) = mx(1,3)
      mx(4,1) = mx(1,4)
      mx(5,1) = mx(1,5)
      mx(3,2) = mx(2,3)
      mx(4,2) = mx(2,4)
      mx(5,2) = mx(2,5)
      mx(4,3) = mx(3,4)
      mx(5,3) = mx(3,5)
      mx(5,4) = mx(4,5)

c      order = 2
c      ordm1 = 1
c      ordm2 = 0

c>      write(*,*) 'mx'
c>      do i=0,ordm1
c>      do j=0,ordm1
c>         write(*,*) i,j, mx(i,j)
c>      end do
c>      end do

      call chinv(order,mx,mxinv)

c>      write(*,*) 'mxinv'
c>      do i=0,ordm1
c>      do j=0,ordm1
c>         write(*,*) i,j, mxinv(i,j)
c>      end do
c>      end do

      do i=1,ordm1
      do j=1,ordm1
        mx0k(i,j) = mx(i,j-1)
      end do
      end do

c>      write(*,*) 'mx0k 5'
c>      do i=1,ordm1
c>      do j=1,ordm1
c>         write(*,*) i,j, mx0k(i,j)
c>      end do
c>      end do

      det = chdet(ordm1,mx0k)
      kbest = ordm1
c      write(*,*) 'det5',det

      do l=ordm1,2,-1
        do i=1,ordm1
          mx0k(i,l) = mx(i,l)
        end do
c>        write(*,*) 'mx0k',l
c>        do i=1,ordm1
c>          do j=1,ordm1
c>            write(*,*) i,j, mx0k(i,j)
c>          end do
c>        end do
        newdet =  chdet(order-1,mx0k)
c        write(*,*) 'detk',l,newdet,newdet/det
        if (cdabs(newdet).gt.cdabs(det)) then          
          kbest = l-1
          det = newdet
        end if
c        write(*,*) 'detk',l,k,newdet/det
      end do         

      do i=1,ordm1
        mx0k(i,1) = mx(i,1)
        mx0k(i,kbest) = mx(i,0)
      end do

c>      write(*,*) 'mx0k',kbest
c>      do i=1,ordm1
c>        do j=1,ordm1
c>          write(*,*) i,j, mx0k(i,j)
c>        end do
c>      end do
c>      write(*,*) 'det',det,chdet(order-1,mx0k)/det

      call chinv(ordm1,mx0k,mx0kinv)
      do i=1,ordm1
        mx0kinv(kbest,i) = 0d0
      end do

c  alternative calculation of mx0kinv from mxinv
c>      write(*,*) 'mx0kinv'
c>      do i=1,ordm1
c>        do j=1,ordm1
c>          write(*,*) i,j, mx0kinv(i,j)
c>          write(*,*) i,j,  mxinv(i,j) 
c>     &        - mxinv(j,kbest)* mxinv(i,0)/ mxinv(kbest,0)
c>        end do
c>      end do
c>
c>      do i=1,ordm1
c>        do j=1,ordm1
c>          mx0kinv(i,j) =  mxinv(i,j) 
c>     &        - mxinv(j,kbest)* mxinv(i,0)/ mxinv(kbest,0)
c>        end do
c>      end do

      rankE = max(rank-1,0)
      call cEp1234(p21,p32,p43,p54,p51,p31,p42,p53,p41,p52,
     &            m12,m22,m32,m42,m52,
     &            En0w0,En1w0,En2w0,En3w0,En4w0,rankE)
      call cEp1234(p20,p32,p43,p54,p50,p30,p42,p53,p40,p52,
     &            m02,m22,m32,m42,m52,
     &            En0w1,En1w1,En2w1,En3w1,En4w1,rankE)
      call cEp1234(p10,p31,p43,p54,p50,p30,p41,p53,p40,p51,
     &            m02,m12,m32,m42,m52,
     &            En0w2,En1w2,En2w2,En3w2,En4w2,rankE)
      call cEp1234(p10,p21,p42,p54,p50,p20,p41,p52,p40,p51,
     &            m02,m12,m22,m42,m52,
     &            En0w3,En1w3,En2w3,En3w3,En4w3,rankE)
      call cEp1234(p10,p21,p32,p53,p50,p20,p31,p52,p30,p51,
     &            m02,m12,m22,m32,m52,
     &            En0w4,En1w4,En2w4,En3w4,En4w4,rankE)
      call cEp1234(p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
     &            m02,m12,m22,m32,m42,
     &            En0w5,En1w5,En2w5,En3w5,En4w5,rankE)
 
      E0m(0) = En0w0
      E0m(1) = En0w1
      E0m(2) = En0w2
      E0m(3) = En0w3
      E0m(4) = En0w4
      E0m(5) = En0w5

      Fn0 = - mxinv(0,0)*E0m(0)
      do 5 i=1,ordm1
         Fn0 = Fn0 + mxinv(i,0)*(E0m(i)-E0m(0))
  5   continue

c>      write(*,*) 'cFp ',eta
c>      write(*,*) 'cFp ',En0w0
c>      write(*,*) 'cFp ',En0w1
c>      write(*,*) 'cFp ',En0w2
c>      write(*,*) 'cFp ',En0w3
c>      write(*,*) 'cFp ',En0w4
c>      write(*,*) 'cFp ',En0w5
c>      write(*,*) 'cFp ',Fn0

      if(rank.eq.0) return
 
      do 15 j=1,ordm1
        Fn1(j) = 0d0
      do 15 i=1,ordm1
        Fn1(j) = Fn1(j) + mx0kinv(j,i)*(E0m(i)- E0m(0))
 15   continue

      if(rank.eq.1) return

      E1m(0,1) = -En0w0 - En1w0(1) - En1w0(2) - En1w0(3) - En1w0(4)
      E1m(0,2) = En1w0(1)
      E1m(0,3) = En1w0(2)
      E1m(0,4) = En1w0(3)
      E1m(0,5) = En1w0(4)
      E1m(1,2) = En1w1(1)
      E1m(1,3) = En1w1(2)
      E1m(1,4) = En1w1(3)
      E1m(1,5) = En1w1(4)
      E1m(2,1) = En1w2(1)
      E1m(2,3) = En1w2(2)
      E1m(2,4) = En1w2(3)
      E1m(2,5) = En1w2(4)
      E1m(3,1) = En1w3(1)
      E1m(3,2) = En1w3(2)
      E1m(3,4) = En1w3(3)
      E1m(3,5) = En1w3(4)
      E1m(4,1) = En1w4(1)
      E1m(4,2) = En1w4(2)
      E1m(4,3) = En1w4(3)
      E1m(4,5) = En1w4(4)
      E1m(5,1) = En1w5(1)
      E1m(5,2) = En1w5(2)
      E1m(5,3) = En1w5(3)
      E1m(5,4) = En1w5(4)

      Fn2(0,0) = 0d0

      do 27 j=1,ordm1
      do 27 k=1,ordm1
        Fn2(j,k) = 0d0
      do 26 i=1,ordm1
        Fn2(j,k) = Fn2(j,k) + mx0kinv(j,i)*(E1m(i,k)- E1m(0,k))
     &                      + mx0kinv(k,i)*(E1m(i,j)- E1m(0,j))
 26   continue  
        Fn2(j,k) = Fn2(j,k)/2d0
        Fn2(k,j) = Fn2(j,k)
 27   continue

      if(rank.eq.2) return

      E2m00(0) = En2w0(0,0)
      E2m00(1) = En2w1(0,0)
      E2m00(2) = En2w2(0,0)
      E2m00(3) = En2w3(0,0)
      E2m00(4) = En2w4(0,0)
      E2m00(5) = En2w5(0,0)

      E2m(0,2,2) = En2w0(1,1) 
      E2m(0,2,3) = En2w0(1,2) 
      E2m(0,2,4) = En2w0(1,3) 
      E2m(0,2,5) = En2w0(1,4) 
      E2m(0,3,3) = En2w0(2,2) 
      E2m(0,3,4) = En2w0(2,3) 
      E2m(0,3,5) = En2w0(2,4) 
      E2m(0,4,4) = En2w0(3,3) 
      E2m(0,4,5) = En2w0(3,4) 
      E2m(0,5,5) = En2w0(4,4) 
      E2m(0,1,2) = -En1w0(1) - En2w0(1,1) - En2w0(1,2) - En2w0(1,3) 
     &            - En2w0(1,4)  
      E2m(0,1,3) = -En1w0(2) - En2w0(1,2) - En2w0(2,2) - En2w0(2,3) 
     &            - En2w0(2,4)  
      E2m(0,1,4) = -En1w0(3) - En2w0(1,3) - En2w0(2,3) - En2w0(3,3) 
     &            - En2w0(3,4)  
      E2m(0,1,5) = -En1w0(4) - En2w0(1,4) - En2w0(2,4) - En2w0(3,4) 
     &            - En2w0(4,4)  
      E2m(0,1,1) = -E1m(0,1) - E2m(0,1,2) - E2m(0,1,3) - E2m(0,1,4)
     &            - E2m(0,1,5)   

      E2m(1,2,2) = En2w1(1,1) 
      E2m(1,2,3) = En2w1(1,2) 
      E2m(1,2,4) = En2w1(1,3) 
      E2m(1,2,5) = En2w1(1,4) 
      E2m(1,3,3) = En2w1(2,2) 
      E2m(1,3,4) = En2w1(2,3) 
      E2m(1,3,5) = En2w1(2,4) 
      E2m(1,4,4) = En2w1(3,3) 
      E2m(1,4,5) = En2w1(3,4) 
      E2m(1,5,5) = En2w1(4,4) 
      E2m(2,1,1) = En2w2(1,1) 
      E2m(2,1,3) = En2w2(1,2) 
      E2m(2,1,4) = En2w2(1,3) 
      E2m(2,1,5) = En2w2(1,4) 
      E2m(2,3,3) = En2w2(2,2) 
      E2m(2,3,4) = En2w2(2,3) 
      E2m(2,3,5) = En2w2(2,4) 
      E2m(2,4,4) = En2w2(3,3) 
      E2m(2,4,5) = En2w2(3,4) 
      E2m(2,5,5) = En2w2(4,4) 
      E2m(3,1,1) = En2w3(1,1) 
      E2m(3,1,2) = En2w3(1,2) 
      E2m(3,1,4) = En2w3(1,3) 
      E2m(3,1,5) = En2w3(1,4) 
      E2m(3,2,2) = En2w3(2,2) 
      E2m(3,2,4) = En2w3(2,3) 
      E2m(3,2,5) = En2w3(2,4) 
      E2m(3,4,4) = En2w3(3,3) 
      E2m(3,4,5) = En2w3(3,4) 
      E2m(3,5,5) = En2w3(4,4) 
      E2m(4,1,1) = En2w4(1,1) 
      E2m(4,1,2) = En2w4(1,2) 
      E2m(4,1,3) = En2w4(1,3) 
      E2m(4,1,5) = En2w4(1,4) 
      E2m(4,2,2) = En2w4(2,2) 
      E2m(4,2,3) = En2w4(2,3) 
      E2m(4,2,5) = En2w4(2,4) 
      E2m(4,3,3) = En2w4(3,3) 
      E2m(4,3,5) = En2w4(3,4) 
      E2m(4,5,5) = En2w4(4,4) 
      E2m(5,1,1) = En2w5(1,1) 
      E2m(5,1,2) = En2w5(1,2) 
      E2m(5,1,3) = En2w5(1,3) 
      E2m(5,1,4) = En2w5(1,4) 
      E2m(5,2,2) = En2w5(2,2) 
      E2m(5,2,3) = En2w5(2,3) 
      E2m(5,2,4) = En2w5(2,4) 
      E2m(5,3,3) = En2w5(3,3) 
      E2m(5,3,4) = En2w5(3,4) 
      E2m(5,4,4) = En2w5(4,4)  

      do 36 j=1,ordm1
        Fn3(0,0,j) = 0d0
      do 35 i=1,ordm1
        Fn3(0,0,j) = Fn3(0,0,j) + mx0kinv(j,i)*(E2m00(i)-E2m00(0))
 35   continue
        Fn3(0,0,j) = Fn3(0,0,j)/3d0
 36   continue

c>      do i=0,ordm1
c>      do k=1,ordm1
c>      do l=k,ordm1
c>        E2m(i,l,k) = E2m(i,k,l)
c>      end do
c>      end do
c>      end do

      do 37 j=1,ordm1
      do 37 k=j,ordm1
      do 37 l=k,ordm1
        Fn3(j,k,l) = 0d0
      do 38 i=1,ordm1
        Fn3(j,k,l) = Fn3(j,k,l) + mx0kinv(j,i)*(E2m(i,k,l)- E2m(0,k,l))
     &                          + mx0kinv(k,i)*(E2m(i,j,l)- E2m(0,j,l))
     &                          + mx0kinv(l,i)*(E2m(i,j,k)- E2m(0,j,k))
 38   continue
        Fn3(j,k,l) =  Fn3(j,k,l)/3d0
        if (sym.eq.1) then
          Fn3(l,j,k) = Fn3(j,k,l)
          Fn3(k,l,j) = Fn3(j,k,l)
          Fn3(j,l,k) = Fn3(j,k,l)
          Fn3(l,k,j) = Fn3(j,k,l)
          Fn3(k,j,l) = Fn3(j,k,l)
        end if
 37    continue
        
      if(rank.eq.3) return

      E3m00(1,2)   = En3w1(0,0,1)
      E3m00(1,3)   = En3w1(0,0,2)
      E3m00(1,4)   = En3w1(0,0,3)
      E3m00(1,5)   = En3w1(0,0,4)
      E3m00(2,1)   = En3w2(0,0,1) 
      E3m00(2,3)   = En3w2(0,0,2)
      E3m00(2,4)   = En3w2(0,0,3)
      E3m00(2,5)   = En3w2(0,0,4)
      E3m00(3,1)   = En3w3(0,0,1) 
      E3m00(3,2)   = En3w3(0,0,2)
      E3m00(3,4)   = En3w3(0,0,3)
      E3m00(3,5)   = En3w3(0,0,4)
      E3m00(4,1)   = En3w4(0,0,1) 
      E3m00(4,2)   = En3w4(0,0,2)
      E3m00(4,3)   = En3w4(0,0,3)
      E3m00(4,5)   = En3w4(0,0,4)
      E3m00(5,1)   = En3w5(0,0,1) 
      E3m00(5,2)   = En3w5(0,0,2)
      E3m00(5,3)   = En3w5(0,0,3)
      E3m00(5,4)   = En3w5(0,0,4)
      E3m00(0,2)   = En3w0(0,0,1)
      E3m00(0,3)   = En3w0(0,0,2)
      E3m00(0,4)   = En3w0(0,0,3)
      E3m00(0,5)   = En3w0(0,0,4)
      E3m00(0,1)   = -E2m00(0) - E3m00(0,2) - E3m00(0,3) - E3m00(0,4) 
     &              - E3m00(0,5)

      do 31 i1=1,ordm2
      do 31 i2=i1,ordm2
      do 31 i3=i2,ordm2
          E3m(0,i1+1,i2+1,i3+1) = En3w0(i1,i2,i3)
          E3m(1,i1+1,i2+1,i3+1) = En3w1(i1,i2,i3)
          E3m(5,i1,i2,i3) = En3w5(i1,i2,i3)
 31   continue    

      E3m(0,1,2,2) = -E2m(0,2,2)-E3m(0,2,2,2)-E3m(0,2,2,3)-E3m(0,2,2,4) 
     &              - E3m(0,2,2,5)  
      E3m(0,1,2,3) = -E2m(0,2,3)-E3m(0,2,2,3)-E3m(0,2,3,3)-E3m(0,2,3,4) 
     &              - E3m(0,2,3,5)  
      E3m(0,1,2,4) = -E2m(0,2,4)-E3m(0,2,2,4)-E3m(0,2,3,4)-E3m(0,2,4,4) 
     &              - E3m(0,2,4,5)  
      E3m(0,1,2,5) = -E2m(0,2,5)-E3m(0,2,2,5)-E3m(0,2,3,5)-E3m(0,2,4,5) 
     &              - E3m(0,2,5,5)  
      E3m(0,1,3,3) = -E2m(0,3,3)-E3m(0,2,3,3)-E3m(0,3,3,3)-E3m(0,3,3,4) 
     &              - E3m(0,3,3,5)  
      E3m(0,1,3,4) = -E2m(0,3,4)-E3m(0,2,3,4)-E3m(0,3,3,4)-E3m(0,3,4,4) 
     &              - E3m(0,3,4,5)  
      E3m(0,1,3,5) = -E2m(0,3,5)-E3m(0,2,3,5)-E3m(0,3,3,5)-E3m(0,3,4,5) 
     &              - E3m(0,3,5,5)  
      E3m(0,1,4,4) = -E2m(0,4,4)-E3m(0,2,4,4)-E3m(0,3,4,4)-E3m(0,4,4,4) 
     &              - E3m(0,4,4,5)  
      E3m(0,1,4,5) = -E2m(0,4,5)-E3m(0,2,4,5)-E3m(0,3,4,5)-E3m(0,4,4,5) 
     &              - E3m(0,4,5,5)  
      E3m(0,1,5,5) = -E2m(0,5,5)-E3m(0,2,5,5)-E3m(0,3,5,5)-E3m(0,4,5,5) 
     &              - E3m(0,5,5,5)  
      E3m(0,1,1,2) = -E2m(0,1,2)-E3m(0,1,2,2)-E3m(0,1,2,3)-E3m(0,1,2,4) 
     &              - E3m(0,1,2,5)  
      E3m(0,1,1,3) = -E2m(0,1,3)-E3m(0,1,2,3)-E3m(0,1,3,3)-E3m(0,1,3,4) 
     &              - E3m(0,1,3,5)  
      E3m(0,1,1,4) = -E2m(0,1,4)-E3m(0,1,2,4)-E3m(0,1,3,4)-E3m(0,1,4,4) 
     &              - E3m(0,1,4,5)  
      E3m(0,1,1,5) = -E2m(0,1,5)-E3m(0,1,2,5)-E3m(0,1,3,5)-E3m(0,1,4,5) 
     &              - E3m(0,1,5,5)  
      E3m(0,1,1,1) = -E2m(0,1,1)-E3m(0,1,1,2)-E3m(0,1,1,3)-E3m(0,1,1,4) 
     &              - E3m(0,1,1,5)  

      E3m(2,1,1,1) = En3w2(1,1,1)
      E3m(2,1,1,3) = En3w2(1,1,2)
      E3m(2,1,1,4) = En3w2(1,1,3)
      E3m(2,1,1,5) = En3w2(1,1,4)
      E3m(2,1,3,3) = En3w2(1,2,2)
      E3m(2,1,3,4) = En3w2(1,2,3)
      E3m(2,1,3,5) = En3w2(1,2,4)
      E3m(2,1,4,4) = En3w2(1,3,3)
      E3m(2,1,4,5) = En3w2(1,3,4)
      E3m(2,1,5,5) = En3w2(1,4,4)
      E3m(2,3,3,3) = En3w2(2,2,2)
      E3m(2,3,3,4) = En3w2(2,2,3)
      E3m(2,3,3,5) = En3w2(2,2,4)
      E3m(2,3,4,4) = En3w2(2,3,3)
      E3m(2,3,4,5) = En3w2(2,3,4)
      E3m(2,3,5,5) = En3w2(2,4,4)
      E3m(2,4,4,4) = En3w2(3,3,3)
      E3m(2,4,4,5) = En3w2(3,3,4)
      E3m(2,4,5,5) = En3w2(3,4,4)
      E3m(2,5,5,5) = En3w2(4,4,4)

      E3m(3,1,1,1) = En3w3(1,1,1)
      E3m(3,1,1,2) = En3w3(1,1,2)
      E3m(3,1,1,4) = En3w3(1,1,3)
      E3m(3,1,1,5) = En3w3(1,1,4)
      E3m(3,1,2,2) = En3w3(1,2,2)
      E3m(3,1,2,4) = En3w3(1,2,3)
      E3m(3,1,2,5) = En3w3(1,2,4)
      E3m(3,1,4,4) = En3w3(1,3,3)
      E3m(3,1,4,5) = En3w3(1,3,4)
      E3m(3,1,5,5) = En3w3(1,4,4)
      E3m(3,2,2,2) = En3w3(2,2,2)
      E3m(3,2,2,4) = En3w3(2,2,3)
      E3m(3,2,2,5) = En3w3(2,2,4)
      E3m(3,2,4,4) = En3w3(2,3,3)
      E3m(3,2,4,5) = En3w3(2,3,4)
      E3m(3,2,5,5) = En3w3(2,4,4)
      E3m(3,4,4,4) = En3w3(3,3,3)
      E3m(3,4,4,5) = En3w3(3,3,4)
      E3m(3,4,5,5) = En3w3(3,4,4)
      E3m(3,5,5,5) = En3w3(4,4,4)

      E3m(4,1,1,1) = En3w4(1,1,1)
      E3m(4,1,1,2) = En3w4(1,1,2)
      E3m(4,1,1,3) = En3w4(1,1,3)
      E3m(4,1,1,5) = En3w4(1,1,4)
      E3m(4,1,2,2) = En3w4(1,2,2)
      E3m(4,1,2,3) = En3w4(1,2,3)
      E3m(4,1,2,5) = En3w4(1,2,4)
      E3m(4,1,3,3) = En3w4(1,3,3)
      E3m(4,1,3,5) = En3w4(1,3,4)
      E3m(4,1,5,5) = En3w4(1,4,4)
      E3m(4,2,2,2) = En3w4(2,2,2)
      E3m(4,2,2,3) = En3w4(2,2,3)
      E3m(4,2,2,5) = En3w4(2,2,4)
      E3m(4,2,3,3) = En3w4(2,3,3)
      E3m(4,2,3,5) = En3w4(2,3,4)
      E3m(4,2,5,5) = En3w4(2,4,4)
      E3m(4,3,3,3) = En3w4(3,3,3)
      E3m(4,3,3,5) = En3w4(3,3,4)
      E3m(4,3,5,5) = En3w4(3,4,4)
      E3m(4,5,5,5) = En3w4(4,4,4)

      Fn4(0,0,0,0) = 0d0

      do 46 j=1,ordm1
      do 46 k=j,ordm1
        Fn4(0,0,j,k) = 0d0
      do 45 i=1,ordm1
        Fn4(0,0,j,k) = Fn4(0,0,j,k)
     &      + mx0kinv(j,i)*(E3m00(i,k)-E3m00(0,k))
     &      + mx0kinv(k,i)*(E3m00(i,j)-E3m00(0,j))
 45   continue
        Fn4(0,0,j,k) = Fn4(0,0,j,k)/4d0
        Fn4(0,0,k,j) = Fn4(0,0,j,k)
 46   continue

c>      do i=0,ordm1
c>      do k=1,ordm1
c>      do l=k,ordm1
c>      do m=l,ordm1
c>        E3m(i,l,k,m) = E3m(i,k,l,m)
c>        E3m(i,k,m,l) = E3m(i,k,l,m)
c>        E3m(i,m,k,l) = E3m(i,k,l,m)
c>        E3m(i,l,m,k) = E3m(i,k,l,m)
c>        E3m(i,m,l,k) = E3m(i,k,l,m)
c>      end do
c>      end do
c>      end do
c>      end do

      do 47 j=1,ordm1
      do 47 k=j,ordm1
      do 47 l=k,ordm1
      do 47 m=l,ordm1
        Fn4(j,k,l,m) = 0d0
      do 48 i=1,ordm1
        Fn4(j,k,l,m) = Fn4(j,k,l,m) 
     &      + mx0kinv(j,i)*(E3m(i,k,l,m)- E3m(0,k,l,m))
     &      + mx0kinv(k,i)*(E3m(i,j,l,m)- E3m(0,j,l,m))
     &      + mx0kinv(l,i)*(E3m(i,j,k,m)- E3m(0,j,k,m))
     &      + mx0kinv(m,i)*(E3m(i,j,k,l)- E3m(0,j,k,l))
 48   continue
        Fn4(j,k,l,m) =  Fn4(j,k,l,m)/4d0
        if (sym.eq.1) then
          Fn4(l,j,k,m) = Fn4(j,k,l,m)
          Fn4(k,l,j,m) = Fn4(j,k,l,m)
          Fn4(j,l,k,m) = Fn4(j,k,l,m)
          Fn4(l,k,j,m) = Fn4(j,k,l,m)
          Fn4(k,j,l,m) = Fn4(j,k,l,m)
          Fn4(j,k,m,l) = Fn4(j,k,l,m)
          Fn4(l,j,m,k) = Fn4(j,k,l,m)
          Fn4(k,l,m,j) = Fn4(j,k,l,m)
          Fn4(j,l,m,k) = Fn4(j,k,l,m)
          Fn4(l,k,m,j) = Fn4(j,k,l,m)
          Fn4(k,j,m,l) = Fn4(j,k,l,m)
          Fn4(j,m,k,l) = Fn4(j,k,l,m)
          Fn4(l,m,j,k) = Fn4(j,k,l,m)
          Fn4(k,m,l,j) = Fn4(j,k,l,m)
          Fn4(j,m,l,k) = Fn4(j,k,l,m)
          Fn4(l,m,k,j) = Fn4(j,k,l,m)
          Fn4(k,m,j,l) = Fn4(j,k,l,m)
          Fn4(m,j,k,l) = Fn4(j,k,l,m)
          Fn4(m,l,j,k) = Fn4(j,k,l,m)
          Fn4(m,k,l,j) = Fn4(j,k,l,m)
          Fn4(m,j,l,k) = Fn4(j,k,l,m)
          Fn4(m,l,k,j) = Fn4(j,k,l,m)
          Fn4(m,k,j,l) = Fn4(j,k,l,m)
        end if
 47   continue

      if (rank.gt.4) then
        write(*,*) 'rank > 4 not implemented in cFp1234'
        write(*,*) 'rank = ',rank
        stop

      E4m0000(0) = En4w0(0,0,0,0) 
      E4m0000(1) = En4w1(0,0,0,0) 
      E4m0000(2) = En4w2(0,0,0,0) 
      E4m0000(3) = En4w3(0,0,0,0) 
      E4m0000(4) = En4w4(0,0,0,0) 
      E4m0000(5) = En4w5(0,0,0,0) 

      do 41 i1=1,ordm2
      do 41 i2=1,ordm2
         E4m00(0,i1+1,i2+1) = En4w0(0,0,i1,i2) 
         E4m00(1,i1+1,i2+1) = En4w1(0,0,i1,i2) 
         E4m00(5,i1,i2) = En4w5(0,0,i1,i2)
      do 41 i3=1,ordm2
      do 41 i4=1,ordm2
         E4m(0,i1+1,i2+1,i3+1,i4+1) = En4w0(i1,i2,i3,i4)
         E4m(1,i1+1,i2+1,i3+1,i4+1) = En4w1(i1,i2,i3,i4)
         E4m(5,i1,i2,i3,i4) = En4w5(i1,i2,i3,i4)
 41   continue    

      E4m00(0,1,2) = -E3m00(0,2)-E4m00(0,2,2)-E4m00(0,2,3)-E4m00(0,2,4)
     &              - E4m00(0,2,5)
      E4m00(0,1,3) = -E3m00(0,3)-E4m00(0,2,3)-E4m00(0,3,3)-E4m00(0,3,4)
     &              - E4m00(0,3,5)
      E4m00(0,1,4) = -E3m00(0,4)-E4m00(0,2,4)-E4m00(0,3,4)-E4m00(0,4,4)
     &              - E4m00(0,4,5)
      E4m00(0,1,5) = -E3m00(0,5)-E4m00(0,2,5)-E4m00(0,3,5)-E4m00(0,4,5)
     &              - E4m00(0,5,5)
      E4m00(0,1,1) = -E3m00(0,1)-E4m00(0,1,2)-E4m00(0,1,3)-E4m00(0,1,4)
     &              - E4m00(0,1,5)

      E4m00(2,1,1) = En4w2(0,0,1,1) 
      E4m00(2,1,3) = En4w2(0,0,1,2) 
      E4m00(2,1,4) = En4w2(0,0,1,3) 
      E4m00(2,1,5) = En4w2(0,0,1,4) 
      E4m00(2,3,3) = En4w2(0,0,2,2) 
      E4m00(2,3,4) = En4w2(0,0,2,3) 
      E4m00(2,3,5) = En4w2(0,0,2,4) 
      E4m00(2,4,4) = En4w2(0,0,3,3) 
      E4m00(2,4,5) = En4w2(0,0,3,4) 
      E4m00(2,5,5) = En4w2(0,0,4,4) 

      E4m00(3,1,1) = En4w3(0,0,1,1) 
      E4m00(3,1,2) = En4w3(0,0,1,2) 
      E4m00(3,1,4) = En4w3(0,0,1,3) 
      E4m00(3,1,5) = En4w3(0,0,1,4) 
      E4m00(3,2,2) = En4w3(0,0,2,2) 
      E4m00(3,2,4) = En4w3(0,0,2,3) 
      E4m00(3,2,5) = En4w3(0,0,2,4) 
      E4m00(3,4,4) = En4w3(0,0,3,3) 
      E4m00(3,4,5) = En4w3(0,0,3,4) 
      E4m00(3,5,5) = En4w3(0,0,4,4) 

      E4m00(4,1,1) = En4w4(0,0,1,1) 
      E4m00(4,1,2) = En4w4(0,0,1,2) 
      E4m00(4,1,3) = En4w4(0,0,1,3) 
      E4m00(4,1,5) = En4w4(0,0,1,4) 
      E4m00(4,2,2) = En4w4(0,0,2,2) 
      E4m00(4,2,3) = En4w4(0,0,2,3) 
      E4m00(4,2,5) = En4w4(0,0,2,4) 
      E4m00(4,3,3) = En4w4(0,0,3,3) 
      E4m00(4,3,5) = En4w4(0,0,3,4) 
      E4m00(4,5,5) = En4w4(0,0,4,4) 

      E4m(0,1,2,2,2) = -E3m(0,2,2,2)   - E4m(0,2,2,2,5) 
     &                - E4m(0,2,2,2,2) - E4m(0,2,2,2,3) - E4m(0,2,2,2,4)
      E4m(0,1,2,2,3) = -E3m(0,2,2,3)   - E4m(0,2,2,3,5) 
     &                - E4m(0,2,2,2,3) - E4m(0,2,2,3,3) - E4m(0,2,2,3,4) 
      E4m(0,1,2,2,4) = -E3m(0,2,2,4)   - E4m(0,2,2,4,5) 
     &                - E4m(0,2,2,2,4) - E4m(0,2,2,3,4) - E4m(0,2,2,4,4) 
      E4m(0,1,2,2,5) = -E3m(0,2,2,5)   - E4m(0,2,2,5,5) 
     &                - E4m(0,2,2,2,5) - E4m(0,2,2,3,5) - E4m(0,2,2,4,5) 
      E4m(0,1,2,3,3) = -E3m(0,2,3,3)   - E4m(0,2,3,3,5) 
     &                - E4m(0,2,2,3,3) - E4m(0,2,3,3,3) - E4m(0,2,3,3,4) 
      E4m(0,1,2,3,4) = -E3m(0,2,3,4)   - E4m(0,2,3,4,5) 
     &                - E4m(0,2,2,3,4) - E4m(0,2,3,3,4) - E4m(0,2,3,4,4) 
      E4m(0,1,2,3,5) = -E3m(0,2,3,5)   - E4m(0,2,3,5,5) 
     &                - E4m(0,2,2,3,5) - E4m(0,2,3,3,5) - E4m(0,2,3,4,5) 
      E4m(0,1,2,4,4) = -E3m(0,2,4,4)   - E4m(0,2,4,4,5) 
     &                - E4m(0,2,2,4,4) - E4m(0,2,3,4,4) - E4m(0,2,4,4,4) 
      E4m(0,1,2,4,5) = -E3m(0,2,4,5)   - E4m(0,2,4,5,5) 
     &                - E4m(0,2,2,4,5) - E4m(0,2,3,4,5) - E4m(0,2,4,4,5) 
      E4m(0,1,2,5,5) = -E3m(0,2,5,5)   - E4m(0,2,5,5,5) 
     &                - E4m(0,2,2,5,5) - E4m(0,2,3,5,5) - E4m(0,2,4,5,5) 
      E4m(0,1,3,3,3) = -E3m(0,3,3,3)   - E4m(0,3,3,3,5) 
     &                - E4m(0,2,3,3,3) - E4m(0,3,3,3,3) - E4m(0,3,3,3,4) 
      E4m(0,1,3,3,4) = -E3m(0,3,3,4)   - E4m(0,3,3,4,5) 
     &                - E4m(0,2,3,3,4) - E4m(0,3,3,3,4) - E4m(0,3,3,4,4) 
      E4m(0,1,3,3,5) = -E3m(0,3,3,5)   - E4m(0,3,3,5,5) 
     &                - E4m(0,2,3,3,5) - E4m(0,3,3,3,5) - E4m(0,3,3,4,5) 
      E4m(0,1,3,4,4) = -E3m(0,3,4,4)   - E4m(0,3,4,4,5) 
     &                - E4m(0,2,3,4,4) - E4m(0,3,3,4,4) - E4m(0,3,4,4,4) 
      E4m(0,1,3,4,5) = -E3m(0,3,4,5)   - E4m(0,3,4,5,5) 
     &                - E4m(0,2,3,4,5) - E4m(0,3,3,4,5) - E4m(0,3,4,4,5) 
      E4m(0,1,3,5,5) = -E3m(0,3,5,5)   - E4m(0,3,5,5,5) 
     &                - E4m(0,2,3,5,5) - E4m(0,3,3,5,5) - E4m(0,3,4,5,5) 
      E4m(0,1,4,4,4) = -E3m(0,4,4,4)   - E4m(0,4,4,4,5) 
     &                - E4m(0,2,4,4,4) - E4m(0,3,4,4,4) - E4m(0,4,4,4,4) 
      E4m(0,1,4,4,5) = -E3m(0,4,4,5)   - E4m(0,4,4,5,5) 
     &                - E4m(0,2,4,4,5) - E4m(0,3,4,4,5) - E4m(0,4,4,4,5) 
      E4m(0,1,4,5,5) = -E3m(0,4,5,5)   - E4m(0,4,5,5,5) 
     &                - E4m(0,2,4,5,5) - E4m(0,3,4,5,5) - E4m(0,4,4,5,5) 
      E4m(0,1,5,5,5) = -E3m(0,5,5,5)   - E4m(0,5,5,5,5) 
     &                - E4m(0,2,5,5,5) - E4m(0,3,5,5,5) - E4m(0,4,5,5,5) 
      E4m(0,1,1,2,2) = -E3m(0,1,2,2)   - E4m(0,1,2,2,5) 
     &                - E4m(0,1,2,2,2) - E4m(0,1,2,2,3) - E4m(0,1,2,2,4) 
      E4m(0,1,1,2,3) = -E3m(0,1,2,3)   - E4m(0,1,2,3,5) 
     &                - E4m(0,1,2,2,3) - E4m(0,1,2,3,3) - E4m(0,1,2,3,4) 
      E4m(0,1,1,2,4) = -E3m(0,1,2,4)   - E4m(0,1,2,4,5) 
     &                - E4m(0,1,2,2,4) - E4m(0,1,2,3,4) - E4m(0,1,2,4,4) 
      E4m(0,1,1,2,5) = -E3m(0,1,2,5)   - E4m(0,1,2,5,5) 
     &                - E4m(0,1,2,2,5) - E4m(0,1,2,3,5) - E4m(0,1,2,4,5) 
      E4m(0,1,1,3,3) = -E3m(0,1,3,3)   - E4m(0,1,3,3,5) 
     &                - E4m(0,1,2,3,3) - E4m(0,1,3,3,3) - E4m(0,1,3,3,4) 
      E4m(0,1,1,3,4) = -E3m(0,1,3,4)   - E4m(0,1,3,4,5) 
     &                - E4m(0,1,2,3,4) - E4m(0,1,3,3,4) - E4m(0,1,3,4,4) 
      E4m(0,1,1,3,5) = -E3m(0,1,3,5)   - E4m(0,1,3,5,5) 
     &                - E4m(0,1,2,3,5) - E4m(0,1,3,3,5) - E4m(0,1,3,4,5) 
      E4m(0,1,1,4,4) = -E3m(0,1,4,4)   - E4m(0,1,4,4,5) 
     &                - E4m(0,1,2,4,4) - E4m(0,1,3,4,4) - E4m(0,1,4,4,4) 
      E4m(0,1,1,4,5) = -E3m(0,1,4,5)   - E4m(0,1,4,5,5) 
     &                - E4m(0,1,2,4,5) - E4m(0,1,3,4,5) - E4m(0,1,4,4,5) 
      E4m(0,1,1,5,5) = -E3m(0,1,5,5)   - E4m(0,1,5,5,5) 
     &                - E4m(0,1,2,5,5) - E4m(0,1,3,5,5) - E4m(0,1,4,5,5) 
      E4m(0,1,1,1,2) = -E3m(0,1,1,2)   - E4m(0,1,1,2,5) 
     &                - E4m(0,1,1,2,2) - E4m(0,1,1,2,3) - E4m(0,1,1,2,4) 
      E4m(0,1,1,1,3) = -E3m(0,1,1,3)   - E4m(0,1,1,3,5) 
     &                - E4m(0,1,1,2,3) - E4m(0,1,1,3,3) - E4m(0,1,1,3,4) 
      E4m(0,1,1,1,4) = -E3m(0,1,1,4)   - E4m(0,1,1,4,5) 
     &                - E4m(0,1,1,2,4) - E4m(0,1,1,3,4) - E4m(0,1,1,4,4) 
      E4m(0,1,1,1,5) = -E3m(0,1,1,5)   - E4m(0,1,1,5,5) 
     &                - E4m(0,1,1,2,5) - E4m(0,1,1,3,5) - E4m(0,1,1,4,5) 
      E4m(0,1,1,1,1) = -E3m(0,1,1,1)   - E4m(0,1,1,1,5) 
     &                - E4m(0,1,1,1,2) - E4m(0,1,1,1,3) - E4m(0,1,1,1,4) 

      E4m(2,1,1,1,1) = En4w2(1,1,1,1)
      E4m(2,1,1,1,3) = En4w2(1,1,1,2)
      E4m(2,1,1,1,4) = En4w2(1,1,1,3)
      E4m(2,1,1,1,5) = En4w2(1,1,1,4)
      E4m(2,1,1,3,3) = En4w2(1,1,2,2)
      E4m(2,1,1,3,4) = En4w2(1,1,2,3)
      E4m(2,1,1,3,5) = En4w2(1,1,2,4)
      E4m(2,1,1,4,4) = En4w2(1,1,3,3)
      E4m(2,1,1,4,5) = En4w2(1,1,3,4)
      E4m(2,1,1,5,5) = En4w2(1,1,4,4)
      E4m(2,1,3,3,3) = En4w2(1,2,2,2)
      E4m(2,1,3,3,4) = En4w2(1,2,2,3)
      E4m(2,1,3,3,5) = En4w2(1,2,2,4)
      E4m(2,1,3,4,4) = En4w2(1,2,3,3)
      E4m(2,1,3,4,5) = En4w2(1,2,3,4)
      E4m(2,1,3,5,5) = En4w2(1,2,4,4)
      E4m(2,1,4,4,4) = En4w2(1,3,3,3)
      E4m(2,1,4,4,5) = En4w2(1,3,3,4)
      E4m(2,1,4,5,5) = En4w2(1,3,4,4)
      E4m(2,1,5,5,5) = En4w2(1,4,4,4)
      E4m(2,3,3,3,3) = En4w2(2,2,2,2)
      E4m(2,3,3,3,4) = En4w2(2,2,2,3)
      E4m(2,3,3,3,5) = En4w2(2,2,2,4)
      E4m(2,3,3,4,4) = En4w2(2,2,3,3)
      E4m(2,3,3,4,5) = En4w2(2,2,3,4)
      E4m(2,3,3,5,5) = En4w2(2,2,4,4)
      E4m(2,3,4,4,4) = En4w2(2,3,3,3)
      E4m(2,3,4,4,5) = En4w2(2,3,3,4)
      E4m(2,3,4,5,5) = En4w2(2,3,4,4)
      E4m(2,3,5,5,5) = En4w2(2,4,4,4)
      E4m(2,4,4,4,4) = En4w2(3,3,3,3)
      E4m(2,4,4,4,5) = En4w2(3,3,3,4)
      E4m(2,4,4,5,5) = En4w2(3,3,4,4)
      E4m(2,4,5,5,5) = En4w2(3,4,4,4)
      E4m(2,5,5,5,5) = En4w2(4,4,4,4)

      E4m(3,1,1,1,1) = En4w3(1,1,1,1)
      E4m(3,1,1,1,2) = En4w3(1,1,1,2)
      E4m(3,1,1,1,4) = En4w3(1,1,1,3)
      E4m(3,1,1,1,5) = En4w3(1,1,1,4)
      E4m(3,1,1,2,2) = En4w3(1,1,2,2)
      E4m(3,1,1,2,4) = En4w3(1,1,2,3)
      E4m(3,1,1,2,5) = En4w3(1,1,2,4)
      E4m(3,1,1,4,4) = En4w3(1,1,3,3)
      E4m(3,1,1,4,5) = En4w3(1,1,3,4)
      E4m(3,1,1,5,5) = En4w3(1,1,4,4)
      E4m(3,1,2,2,2) = En4w3(1,2,2,2)
      E4m(3,1,2,2,4) = En4w3(1,2,2,3)
      E4m(3,1,2,2,5) = En4w3(1,2,2,4)
      E4m(3,1,2,4,4) = En4w3(1,2,3,3)
      E4m(3,1,2,4,5) = En4w3(1,2,3,4)
      E4m(3,1,2,5,5) = En4w3(1,2,4,4)
      E4m(3,1,4,4,4) = En4w3(1,3,3,3)
      E4m(3,1,4,4,5) = En4w3(1,3,3,4)
      E4m(3,1,4,5,5) = En4w3(1,3,4,4)
      E4m(3,1,5,5,5) = En4w3(1,4,4,4)
      E4m(3,2,2,2,2) = En4w3(2,2,2,2)
      E4m(3,2,2,2,4) = En4w3(2,2,2,3)
      E4m(3,2,2,2,5) = En4w3(2,2,2,4)
      E4m(3,2,2,4,4) = En4w3(2,2,3,3)
      E4m(3,2,2,4,5) = En4w3(2,2,3,4)
      E4m(3,2,2,5,5) = En4w3(2,2,4,4)
      E4m(3,2,4,4,4) = En4w3(2,3,3,3)
      E4m(3,2,4,4,5) = En4w3(2,3,3,4)
      E4m(3,2,4,5,5) = En4w3(2,3,4,4)
      E4m(3,2,5,5,5) = En4w3(2,4,4,4)
      E4m(3,4,4,4,4) = En4w3(3,3,3,3)
      E4m(3,4,4,4,5) = En4w3(3,3,3,4)
      E4m(3,4,4,5,5) = En4w3(3,3,4,4)
      E4m(3,4,5,5,5) = En4w3(3,4,4,4)
      E4m(3,5,5,5,5) = En4w3(4,4,4,4)

      E4m(4,1,1,1,1) = En4w4(1,1,1,1)
      E4m(4,1,1,1,2) = En4w4(1,1,1,2)
      E4m(4,1,1,1,3) = En4w4(1,1,1,3)
      E4m(4,1,1,1,5) = En4w4(1,1,1,4)
      E4m(4,1,1,2,2) = En4w4(1,1,2,2)
      E4m(4,1,1,2,3) = En4w4(1,1,2,3)
      E4m(4,1,1,2,5) = En4w4(1,1,2,4)
      E4m(4,1,1,3,3) = En4w4(1,1,3,3)
      E4m(4,1,1,3,5) = En4w4(1,1,3,4)
      E4m(4,1,1,5,5) = En4w4(1,1,4,4)
      E4m(4,1,2,2,2) = En4w4(1,2,2,2)
      E4m(4,1,2,2,3) = En4w4(1,2,2,3)
      E4m(4,1,2,2,5) = En4w4(1,2,2,4)
      E4m(4,1,2,3,3) = En4w4(1,2,3,3)
      E4m(4,1,2,3,5) = En4w4(1,2,3,4)
      E4m(4,1,2,5,5) = En4w4(1,2,4,4)
      E4m(4,1,3,3,3) = En4w4(1,3,3,3)
      E4m(4,1,3,3,5) = En4w4(1,3,3,4)
      E4m(4,1,3,5,5) = En4w4(1,3,4,4)
      E4m(4,1,5,5,5) = En4w4(1,4,4,4)
      E4m(4,2,2,2,2) = En4w4(2,2,2,2)
      E4m(4,2,2,2,3) = En4w4(2,2,2,3)
      E4m(4,2,2,2,5) = En4w4(2,2,2,4)
      E4m(4,2,2,3,3) = En4w4(2,2,3,3)
      E4m(4,2,2,3,5) = En4w4(2,2,3,4)
      E4m(4,2,2,5,5) = En4w4(2,2,4,4)
      E4m(4,2,3,3,3) = En4w4(2,3,3,3)
      E4m(4,2,3,3,5) = En4w4(2,3,3,4)
      E4m(4,2,3,5,5) = En4w4(2,3,4,4)
      E4m(4,2,5,5,5) = En4w4(2,4,4,4)
      E4m(4,3,3,3,3) = En4w4(3,3,3,3)
      E4m(4,3,3,3,5) = En4w4(3,3,3,4)
      E4m(4,3,3,5,5) = En4w4(3,3,4,4)
      E4m(4,3,5,5,5) = En4w4(3,4,4,4)
      E4m(4,5,5,5,5) = En4w4(4,4,4,4)
      end if
      end

************************************************************************
c      function chdet(n,cai,ca,cd,cq)
      function chdet(n,cai)
************************************************************************
*     Calculation of the determinant of a complex n x n matrix         *
*     The matrix is first reduced via Householder transformations into *
*     the form A = QR                                                  *
*     where  Q is a orthogonal and R a upper triangular matrix.        *
*     See Stoer, Numerische Mathematik, chapter 4.7                    *
*     The determinant is then obtained as                              *
*     Det(a) = prod_i (-r_ii)                                          *
*     ca(i,j) j.ge.i contains the elements of u_i that form the        *
*             Householder matrix P_i = 1 - beta_i u x u^H              *
*     ca(i,j) j.lt.i contains the nondiagonal entries of R             *
*     cd contains the diagonal elements of R                           *
*----------------------------------------------------------------------*
*     20.04.04  Ansgar Denner     last changed  08.06.04               *
************************************************************************
      implicit   none
      integer    n
      complex*16 chdet
c      dimension  ai(n,n),a(n,n),d(n),q(n,n)
      real*8     aabs,sigma,rs,beta,norm
      complex*16 ca(n,n),cd(n),cq(n,n),cr(n,n),catest(n,n)
      complex*16 cai(n,n)
      complex*16 cdet,csum,cs
      integer    i,j,k
c
c saving of input
c
      do 10 j=1,n
        do 20 i=1,n
          ca(i,j) = cai(i,j)
 20     continue
 10   continue

c     write(*,*) 
c>      write(*,*) 'n= ',n
c>      write(*,*) 'ca'
c>      do  i=1,n
c>         write(*,2) (ca(i,j),j=1,n)
c>      end do
c> 2    format(5('(',g11.4,',',g11.4,') ':))

c
c  Householder transformation
c
      do 100 j=1,n
c
c  calculation of transformation matrix
c 
        aabs = cdabs(ca(j,j))
        sigma = aabs**2
        do 200 i=j+1,n
           sigma = sigma +  cdabs(ca(i,j))**2
 200    continue
        if (sigma.eq.0) then
          chdet = 0d0
          write(*,*) 'sigma = chdet =',chdet
          return
        end if
        rs = sqrt(sigma)
        if (aabs.gt.0d0) then
          cs = rs* ca(j,j)/aabs
        else
          cs = rs
        end if
        beta = 1d0/(sigma + rs * aabs)
        ca(j,j) = ca(j,j) + cs

c
c multiplication of a(i,j) with the transformation matrix
c
        cd(j) = -cs
        do 300 k=j+1,n
          csum = 0d0
          do 400 i=j,n
            csum = csum + dconjg(ca(i,j))*ca(i,k)
 400      continue
          csum = csum*beta
          do 500 i=j,n
            ca(i,k) = ca(i,k) - ca(i,j) * csum
 500      continue
 300    continue
c      write(*,*) 'j= ',j
c      write(*,*) 'caj'
c      do  l=1,n
c         write(*,2) (ca(l,m),m=1,n)
c      end do

 100  continue  

c     write(*,*) 'catrans'
c     do   i=1,n
c       write(*,2) (ca(i,k),k=1,n)
c     end do
c
c calculation of determinant
c
      cdet = 1d0
      do 600 i=1,n
        cdet = -cdet * cd(i)  
 600  continue
      chdet = cdet

c     write(*,*) 'chdet ',chdet

      return

c
c calculation of matrix Q
c
      do 700 j=n,1,-1
        norm = 0d0         
        do 800 i=j,n
          norm = norm + cdabs(ca(i,j))**2
 800    continue
        cq(j,j) = 1 - 2d0 * cdabs(ca(j,j))**2/norm
        do 900 i=j+1,n
          cq(i,j) = - 2d0 * ca(i,j)*dconjg(ca(j,j))/norm
 900    continue  
        do 1000 k=j+1,n
          csum = 0d0
          do 1100 i=j+1,n
            csum = csum + dconjg(ca(i,j))*cq(i,k)
 1100     continue
          csum = 2d0*csum/norm
          cq(j,k) =  - ca(j,j) * csum
          do 1200 i=j+1,n
            cq(i,k) = cq(i,k) - ca(i,j) * csum
 1200     continue             
 1000   continue
 700  continue
c     write(*,*) 'cq'
c     do   i=1,n
c       write(*,2) (cq(i,k),k=1,n)
c     end do

c
c determination of matrix R
c
      do 1300 i=1,n
      do 1300 k=1,n
        cr(i,k) = 0d0
 1300 continue
      do 1400 i=1,n
        cr(i,i) = cd(i)
      do 1400 k=i+1,n
        cr(i,k) = ca(i,k)
 1400 continue
c     write(*,*) 'cr'
c     do   i=1,n
c       write(*,2) (cr(i,k),k=1,n)
c     end do

c
c check of unitarty of Q
c
      do 1500 i=1,n
      do 1500 k=1,n
        catest(i,k) = 0d0
      do 1500 j=1,n
        catest(i,k) = catest(i,k) + cq(i,j)*dconjg(cq(k,j)) 
 1500 continue  
c     write(*,*) 'cqtest'
c     do   i=1,n
c        write(*,2) (catest(i,k),k=1,n)
c     end do

c
c check of decomposition
c
      do 1600 i=1,n
      do 1600 k=1,n
        catest(i,k) = 0d0
      do 1600 j=1,k
        catest(i,k) = catest(i,k) + cq(i,j)*cr(j,k) 
 1600 continue  
c     write(*,*) 'catest'
c     do   i=1,n
c        write(*,2) (catest(i,k),k=1,n)
c     end do

      end

************************************************************************
c      function chdet(n,cai,ca,cd,cq)
      subroutine chinv(n,cai,cainv)
************************************************************************
*     Calculation of the inverse of a complex n x n matrix             *
*     The matrix is first reduced via Householder transformations into *
*     the form A = QR                                                  *
*     where  Q is a hermitean and R a upper triangular matrix.         *
*     See Stoer, Numerische Mathematik, chapter 4.7                    *
*     The inverse is then obtained as                                  *
*     Ainv = Rinv * Qadj                                               *
*     ca(i,j) j.ge.i contains the elements of u_i that form the        *
*             Householder matrix P_i = 1 - beta_i u x u^H              *
*     ca(i,j) j.lt.i contains the nondiagonal entries of R             *
*     cd contains the diagonal elements of R                           *
*----------------------------------------------------------------------*
*     07.06.05  Ansgar Denner     last changed  16.06.05               *
************************************************************************
      implicit   none
      integer    n
      complex*16 chdet
c      dimension  ai(n,n),a(n,n),d(n),q(n,n)
      real*8     aabs,sigma,rs,beta,norm
      complex*16 ca(n,n),cd(n),cq(n,n),cr(n,n),catest(n,n),crtest(n,n)
      complex*16 cai(n,n),crinv(n,n),cainv(n,n)
c      complex*16 cdet
      complex*16 csum,cs
      integer    i,j,k
c
c saving of input
c
      do 10 j=1,n
        do 20 i=1,n
          ca(i,j) = cai(i,j)
 20     continue
 10   continue

c>      write(*,*) 
c>      write(*,*) 'n= ',n
c>      write(*,*) 'ca'
c>      do  i=1,n
c>         write(*,2) (ca(i,j),j=1,n)
c>      end do
 2    format(5('(',g11.4,',',g11.4,') ':))

c
c  Householder transformation
c
      do 100 j=1,n
c
c  calculation of transformation matrix
c 
        aabs = cdabs(ca(j,j))
        sigma = aabs**2
        do 200 i=j+1,n
           sigma = sigma +  cdabs(ca(i,j))**2
 200    continue
        if (sigma.eq.0) then
          chdet = 0d0
          write(*,*) 'sigma = chdet =',chdet
          return
        end if
        rs = sqrt(sigma)
        if (aabs.gt.0d0) then
          cs = rs* ca(j,j)/aabs
        else
          cs = rs
        end if
        beta = 1d0/(sigma + rs * aabs)
        ca(j,j) = ca(j,j) + cs

c
c multiplication of a(i,j) with the transformation matrix
c
        cd(j) = -cs
        do 300 k=j+1,n
          csum = 0d0
          do 400 i=j,n
            csum = csum + dconjg(ca(i,j))*ca(i,k)
 400      continue
          csum = csum*beta
          do 500 i=j,n
            ca(i,k) = ca(i,k) - ca(i,j) * csum
 500      continue
 300    continue
c      write(*,*) 'j= ',j
c      write(*,*) 'caj'
c      do  l=1,n
c         write(*,2) (ca(l,m),m=1,n)
c      end do

 100  continue  

c     write(*,*) 'catrans'
c     do   i=1,n
c       write(*,2) (ca(i,k),k=1,n)
c     end do
c
c calculation of determinant
c
c      cdet = 1d0
c      do 600 i=1,n
c        cdet = -cdet * cd(i)  
c 600  continue
c      chdet = cdet

cc     write(*,*) 'chdet ',chdet

c      return

c
c calculation of matrix Q
c
      do 700 j=n,1,-1
        norm = 0d0         
        do 800 i=j,n
          norm = norm + cdabs(ca(i,j))**2
 800    continue
        cq(j,j) = 1 - 2d0 * cdabs(ca(j,j))**2/norm
        do 900 i=j+1,n
          cq(i,j) = - 2d0 * ca(i,j)*dconjg(ca(j,j))/norm
 900    continue  
        do 1000 k=j+1,n
          csum = 0d0
          do 1100 i=j+1,n
            csum = csum + dconjg(ca(i,j))*cq(i,k)
 1100     continue
          csum = 2d0*csum/norm
          cq(j,k) =  - ca(j,j) * csum
          do 1200 i=j+1,n
            cq(i,k) = cq(i,k) - ca(i,j) * csum
 1200     continue             
 1000   continue
 700  continue
c>      write(*,*) 'cq'
c>      do   i=1,n
c>        write(*,2) (cq(i,k),k=1,n)
c>      end do

c
c determination of matrix R
c
      do 1300 i=1,n
      do 1300 k=1,n
        cr(i,k) = 0d0
 1300 continue
      do 1400 i=1,n
        cr(i,i) = cd(i)
      do 1400 k=i+1,n
        cr(i,k) = ca(i,k)
 1400 continue
c>      write(*,*) 'cr'
c>      do   i=1,n
c>        write(*,2) (cr(i,k),k=1,n)
c>      end do

c
c determination of inverse of R
c
      do 2000 i=1,n
      do 2000 k=1,n
        crinv(i,k) = 0d0
 2000 continue
      do 2100 k=1,n
        crinv(k,k) = 1d0/cr(k,k)
      do 2100 i=1,k-1
      do 2200 j=i,k-1
        crinv(i,k) = crinv(i,k) - crinv(i,j)* cr(j,k)/cr(k,k)
 2200 continue
 2100 continue
c>      write(*,*) 'crinv'
c>      do   i=1,n
c>        write(*,2) (crinv(i,k),k=1,n)
c>      end do

c
c determination of inverse of A
c
      do 2300 i=1,n
      do 2300 k=1,n
        cainv(i,k) = 0d0
      do 2300 j=1,n
        cainv(i,k) = cainv(i,k) + crinv(i,j)*dconjg(cq(k,j))
 2300 continue
c>      write(*,*) 'cainv'
c>      do   i=1,n
c>        write(*,2) (cainv(i,k),k=1,n)
c>      end do

      return

c
c check of unitarty of Q
c
      do 1500 i=1,n
      do 1500 k=1,n
        catest(i,k) = 0d0
      do 1500 j=1,n
        catest(i,k) = catest(i,k) + cq(i,j)*dconjg(cq(k,j)) 
 1500 continue  
      write(*,*) 'cqtest'
      do   i=1,n
         write(*,2) (catest(i,k),k=1,n)
      end do

c
c check of decomposition
c
      do 1600 i=1,n
      do 1600 k=1,n
        catest(i,k) = 0d0
      do 1600 j=1,k
        catest(i,k) = catest(i,k) + cq(i,j)*cr(j,k) 
 1600 continue  
      write(*,*) 'catest'
      do   i=1,n
        write(*,2) (catest(i,k),k=1,n)
      end do

c
c check of inverse of Q
c
      do 2400 i=1,n
      do 2400 k=1,n
        crtest(i,k) = 0d0
      do 2400 j=1,n
        crtest(i,k) = crtest(i,k) + cr(i,j)*crinv(j,k) 
 2400 continue  
c>      write(*,*) 'crinvtest'
c>      do   i=1,n
c>        write(*,2) (crtest(i,k),k=1,n)
c>      end do

c
c check of inverse of A
c
      do 2500 i=1,n
      do 2500 k=1,n
        catest(i,k) = 0d0
      do 2500 j=1,n
        catest(i,k) = catest(i,k) + cai(i,j)*cainv(j,k) 
 2500 continue  
c>      write(*,*) 'cainvtest'
c>      do   i=1,n
c>        write(*,2) (catest(i,k),k=1,n)
c>      end do


      end
***********************************************************************
      function cE0f(p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
     &                                    m02,m12,m22,m32,m42)
***********************************************************************
*     5-point scalar function                                         *
***********************************************************************
*    11.03.04 Ansgar Denner    last changed  03.05.04                 *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p32,p43,p40,p20,p31,p42,p30,p41
      complex*16 m02,m12,m22,m32,m42
      complex*16 cE0f
      complex*16 En0,En1(4),En2(0:4,0:4),En3(0:4,0:4,0:4)
      complex*16 En4(0:4,0:4,0:4,0:4)

c      write(*,*) 'cE0f in',p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
c     &                                    m02,m12,m22,m32,m42

      call cEp1234(p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
     &                                    m02,m12,m22,m32,m42,
     &                   En0,En1,En2,En3,En4,0)

      cE0f = En0

      end 

***********************************************************************
      function cF0f(p10,p21,p32,p43,p54,p50
     &                 ,p20,p31,p42,p53,p40,p51,p30,p41,p52
     &                 ,m02,m12,m22,m32,m42,m52)
***********************************************************************
*     6-point scalar function                                         *
***********************************************************************
*    11.03.04 Ansgar Denner    last changed  03.05.04                 *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p32,p43,p54,p50
     &          ,p20,p31,p42,p53,p40,p51,p30,p41,p52
      complex*16 m02,m12,m22,m32,m42,m52
      complex*16 cF0f
      complex*16 Fn0,Fn1(5),Fn2(0:5,0:5),Fn3(0:5,0:5,0:5)
      complex*16 Fn4(0:5,0:5,0:5,0:5)

C>      write(*,*) 'cF0f in ',p10,p21,p32,p43,p54,p50
C>     &                 ,p20,p31,p42,p53,p40,p51,p30,p41,p52
C>     &                 ,m0,m1,m2,m3,m4,m5

      call  cFp1234(p10,p21,p32,p43,p54,p50
     &                 ,p20,p31,p42,p53,p40,p51,p30,p41,p52
     &                 ,m02,m12,m22,m32,m42,m52
     &                 ,Fn0,Fn1,Fn2,Fn3,Fn4
     &                 ,0)

      cF0f = Fn0

      end 
