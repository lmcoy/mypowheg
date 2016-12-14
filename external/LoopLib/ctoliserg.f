*
*    accgy2 estimate is in general wrong for critical cases!
*
*
***********************************************************************
*                                                                     *
*     reduction of 3-, 4- and 5-point tensor functions                *
*     for small or vanishing Gram determinants                        *
*                                                                     *
***********************************************************************
*                                                                     *
*     last changed  24.01.06  Ansgar Denner                           *
*                                                                     *
***********************************************************************
* subroutines:                                                        *
* cBp12345                                                            *
* cCp12345,cCy12345,cCg12345,cCgy12345,cCgp12345                      *
* cDp12345,cDy12345,cDg12345,cDgy12345,cDgp12345                      *
* cacheoff,cacheon,setcachelevel                                      *
* cacheinit,cachewrite,cacheread                                      *
* writeinit,writecount                                                *
* functions:                                                          *
***********************************************************************
      subroutine cBp12345(p10,m02,m12,B0,B1,B2,B3,B4,B5,B6,
     &    rankin,switchin)
***********************************************************************
*     2-point tensor coefficient functions B0,B1,B2,B3,B4,B5          *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     definition of integrals:                                        *
*     {B0,B1,B2,B3}= 1/(i pi^2)*\int d^4 q {1,q,q.q,q.q.q} *          *
*        1/[(q**2 - m02)*[(q+k1)**2-m12] ]                            *
*     arguments:   p10=k1**2                                          *
*     definition of tensor coefficients:                              *
*          B1 = B1 *k1                                                *
*          B2 = B2(1,1) *k1*k1  + B2(0,0) * g                         *
*          B3 = B3(1,1,1) *k1*k1*k1 + B3(0,0,i) * (k1*g + 2*sym)      *
*          etc.                                                       *
*---------------------------------------------------------------------*
*     29.10.04  Ansgar Denner     last changed  25.02.05              *
***********************************************************************
      implicit   none
      complex*16 p10
      complex*16 m02,m12
      complex*16 B0,B1,B2(0:1,0:1),B3(0:1,0:1,0:1)
      complex*16 B4(0:1,0:1,0:1,0:1),B5(0:1,0:1,0:1,0:1,0:1)
      complex*16 B6(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 cBnf,cB00nf,cB0000nf
      real*8     mudim2
      integer    rankin,rank,switch,switchin
      integer    ltest

      common /uv/ mudim2

      common /ltest/ ltest
 
      integer flag
      data flag /0/
      save flag

      complex*16 fct(86),x(15)
      logical    nocalc
      integer    type
      character  name*11

c      write(*,*) 'cB12345 in',p10,m02,m12,rankin

      if (rankin.lt.0.or.(rankin.gt.5.and.flag.eq.0.or.
     &                  rankin.eq.5.and.switchin.gt.0)) then
        write(*,*) 'rank not implemented in cB5'
        write(*,*) 'rank   = ',rankin
        write(*,*) 'switch = ',switchin
        flag = 1
        stop
      end if

      x(1)=p10
      x(2)=m02
      x(3)=m12
      type = 4
      rank = rankin
      switch = switchin
      name='cBp12345'

      call cacheread(fct,x,15,3,type,rank,switch,name,nocalc)
      if(nocalc)then
        B0=fct(1)
        B1=fct(2)
        B2(0,0)=fct(3)
        B2(1,1)=fct(4)
        B3(0,0,1)=fct(5)
        B3(1,1,1)=fct(6)
        B4(0,0,0,0)=fct(7)
        B4(0,0,1,1)=fct(8)
        B4(1,1,1,1)=fct(9)
        B5(0,0,0,0,1)=fct(10)
        B5(0,0,1,1,1)=fct(11)
        B5(1,1,1,1,1)=fct(12)
        B6(0,0,0,0,1,1)=fct(13)
        B6(0,0,1,1,1,1)=fct(14)
        B6(1,1,1,1,1,1)=fct(15)
        return
      endif
      
      

      B0   = cBnf(0,p10,m02,m12)

      if (rank.eq.0) then
        if (switch.eq.-1) then 
          B1   = cBnf(1,p10,m02,m12)
        end if
        goto 521
      end if
        
      B1   = cBnf(1,p10,m02,m12)

      if (rank.eq.1) then
        if (switch.eq.-1) then 
          B2(1,1)   = cBnf(2,p10,m02,m12)
        else if (switch.eq.1) then 
          B2(0,0)   = cB00nf(0,p10,m02,m12)
        end if
        goto 521
      end if

      B2(1,1) = cBnf(2,p10,m02,m12)
      B2(0,0) = cB00nf(0,p10,m02,m12)

      if (rank.eq.2) then
        if (switch.eq.-1) then 
          B3(0,0,1)   = cB00nf(1,p10,m02,m12)
          B3(1,1,1)   = cBnf(3,p10,m02,m12)
        else if (switch.eq.1) then 
          B3(0,0,1)   = cB00nf(1,p10,m02,m12)
        end if
        goto 521
      end if

      B3(0,0,1) = cB00nf(1,p10,m02,m12)
      B3(1,1,1) = cBnf(3,p10,m02,m12)
     
      if (rank.eq.3) then
        if (switch.eq.-1) then 
          B4(0,0,1,1)   = cB00nf(2,p10,m02,m12)
          B4(1,1,1,1)   = cBnf(4,p10,m02,m12)
        else if (switch.eq.1) then 
          B4(0,0,0,0) = cB0000nf(0,p10,m02,m12)
          B4(0,0,1,1)   = cB00nf(2,p10,m02,m12)
        end if
        goto 521
      end if

      B4(0,0,0,0) = cB0000nf(0,p10,m02,m12)
      B4(0,0,1,1) = cB00nf(2,p10,m02,m12)
      B4(1,1,1,1) = cBnf(4,p10,m02,m12)
     
      if (rank.eq.4) then
        if (switch.eq.-1) then 
          B5(0,0,0,0,1) = cB0000nf(1,p10,m02,m12)
          B5(0,0,1,1,1)   = cB00nf(3,p10,m02,m12)
          B5(1,1,1,1,1)   = cBnf(5,p10,m02,m12)
        else if (switch.eq.1) then 
          B5(0,0,0,0,1) = cB0000nf(1,p10,m02,m12)
          B5(0,0,1,1,1)   = cB00nf(3,p10,m02,m12)
        end if
        goto 521
      end if

      B5(0,0,0,0,1) = cB0000nf(1,p10,m02,m12)
      B5(0,0,1,1,1) = cB00nf(3,p10,m02,m12)
      B5(1,1,1,1,1) = cBnf(5,p10,m02,m12)
     
      if (rank.eq.5) then
        if (switch.eq.-1) then 
          B6(0,0,0,0,1,1) = cB0000nf(2,p10,m02,m12)
          B6(0,0,1,1,1,1)   = cB00nf(4,p10,m02,m12)
          B6(1,1,1,1,1,1)   = cBnf(6,p10,m02,m12)
        else if (switch.eq.1) then 
          B6(0,0,0,0,1,1) = cB0000nf(2,p10,m02,m12)
          B6(0,0,1,1,1,1)   = cB00nf(4,p10,m02,m12)
        end if
        goto 521
      end if

      B6(0,0,0,0,1,1) = cB0000nf(2,p10,m02,m12)
      B6(0,0,1,1,1,1) = cB00nf(4,p10,m02,m12)
      B6(1,1,1,1,1,1) = cBnf(6,p10,m02,m12)


 521  fct(1)=B0
      fct(2)=B1
      fct(3)=B2(0,0)
      fct(4)=B2(1,1)
      fct(5)=B3(0,0,1)
      fct(6)=B3(1,1,1)
      fct(7)=B4(0,0,0,0)
      fct(8)=B4(0,0,1,1)
      fct(9)=B4(1,1,1,1)
      fct(10)=B5(0,0,0,0,1)
      fct(11)=B5(0,0,1,1,1)
      fct(12)=B5(1,1,1,1,1)
      fct(13)=B6(0,0,0,0,1,1)
      fct(14)=B6(0,0,1,1,1,1)
      fct(15)=B6(1,1,1,1,1,1)
      call cachewrite(fct,15,type)

      end
***********************************************************************

      subroutine cCp12345(p10,p21,p20,m02,m12,m22,
     &                C0,C1,C2,C3,C4,C5,C6,rankin)
***********************************************************************
*     3-point tensor coefficient functions C0,C1,C2,C3,C4,C5          *
*     up to 5 integration momenta in numerator                        *
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
*     20.10.04  Ansgar Denner     last changed  18.05.05              *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p20
      complex*16 q10,q21,q20
      complex*16 m02,m12,m22
      complex*16 mm02,mm12,mm22
      complex*16 C0,C1(2),C2(0:2,0:2),C3(0:2,0:2,0:2)
      complex*16 C4(0:2,0:2,0:2,0:2),C5(0:2,0:2,0:2,0:2,0:2)
      complex*16 C6(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 Ct0,Ct1(2),Ct2(0:2,0:2),Ct3(0:2,0:2,0:2)
      complex*16 Ct4(0:2,0:2,0:2,0:2),Ct5(0:2,0:2,0:2,0:2,0:2)
      complex*16 Ct6(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 f1,f2,zadjf1,zadjf2
      real*8     detz
      real*8     mf,mp2,mm2,mpm2,norm,rm02
      real*8     adetz,mzadjf,mzadjfd,mxadj,adety
      real*8     accstop,acctest,accinf,accB,accworst
      real*8     accpv,accg,accgy,accy,accgp,accgy2
      integer    rankin,switch,rank,rankacc,ordg3acc,ordgy3acc,ordgp3acc
      integer    ltest,sym
      integer    testout
      integer    errout,i1,i2,i3,i4,i5,atest,ordg,ordgy,lc
      integer    flag,counterr,counttest
      integer    Ccount(0:40)
      complex*16 elimcminf2
      logical    errwrite

      common /sym/    sym
      common /ltest/  ltest
      integer    version3,version4
      common /version/ version3,version4
      integer    ordg3,ordgy3,ordgp3
      integer    ordg4,ordgy4,ordgp4
      common /ctord3/ ordg3,ordgy3,ordgp3
      common /ctord4/ ordg4,ordgy4,ordgp4

      common /Ccount/ Ccount

      real*8     accbad1,accbad2,accbad3
      common /accbad/ accbad1,accbad2,accbad3


      real*8     acc
c 16.02.05: accstop changed from 1d-2 to 1d0
c 21.02.05: accstop changed from 1d0  to 1d40
c 23.01.08: accstop changed from 1d40 to 1d0
      data acc /1d-4/, acctest /1d-4/, accstop /1d0/, accB /1d-15/
c      data acc /1d-4/, acctest /1d-4/, accstop /1d40/, accB /1d-15/
c      data acc /1d-4/, acctest /1d-4/, accstop /1d-2/, accB /1d-15/
      data accworst /1d-2/, accinf /1d40/
      data testout /13/ 
      data errout /93/
c      data ordg /2/,ordgy /2/
      data ordg /2/,ordgy /1/
      data version3 /0/, version4 /0/
      data flag /0/, counterr /0/, counttest /0/
      save flag,accworst

c      call cCp1234sd(p10,p21,p20,m02,m12,m22,
c     &    C0,C1,C2,C3,C4,rank)
c      return


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
c>      common/cacheh/nf2,nf3
c>      integer nf2(maxt),nf3(maxt)
c>
c>      integer countC
c>      data countC /0/
c>      common /countC/ countC
c>
c>      if (cachelevel.ge.3) then
c>      countC = countC+1
c>      end if

c      sym = 1

c      write(*,*) 'cCp12345 in ',p10,p21,p20,m02,m12,m22,rankin

      x(1)=p10
      x(2)=p21
      x(3)=p20
      x(4)=m02
      x(5)=m12
      x(6)=m22
      type = 3
      rank = rankin
      switch = 0
      name='cCp12345'

c      rank = 5

      call cacheread(fct,x,50,6,type,rank,switch,name,nocalc)
      if(nocalc)then
        C0=fct(1)
        if (rank.eq.0) goto 511
        C1(1)=fct(2)
        C1(2)=fct(3)
        if (rank.eq.1) goto 511
c        if (rank.eq.1.and.switch.eq.0) goto 511
        C2(0,0)=fct(4)
c        if (rank.eq.1) goto 511
        C2(1,1)=fct(5)
        C2(1,2)=fct(6)
        C2(2,2)=fct(7)
        if (rank.eq.2) goto 511
c        if (rank.eq.2.and.switch.eq.0) goto 511
        C3(0,0,1)=fct(8)
        C3(0,0,2)=fct(9)
c        if (rank.eq.2) goto 511
        C3(1,1,1)=fct(10)
        C3(1,1,2)=fct(11)
        C3(1,2,2)=fct(12)
        C3(2,2,2)=fct(13)
        if (rank.eq.3) goto 511
c        if (rank.eq.3.and.switch.eq.0) goto 511
        C4(0,0,0,0)=fct(14)
        C4(0,0,1,1)=fct(15)
        C4(0,0,1,2)=fct(16)
        C4(0,0,2,2)=fct(17)
c        if (rank.eq.3) goto 511
        C4(1,1,1,1)=fct(18)
        C4(1,1,1,2)=fct(19)
        C4(1,1,2,2)=fct(20)
        C4(1,2,2,2)=fct(21)
        C4(2,2,2,2)=fct(22)
        if (rank.eq.4) goto 511
c        if (rank.eq.4.and.switch.eq.0) goto 511
        C5(0,0,0,0,1)=fct(23)
        C5(0,0,0,0,2)=fct(24)
        C5(0,0,1,1,1)=fct(25)
        C5(0,0,1,1,2)=fct(26)
        C5(0,0,1,2,2)=fct(27)
        C5(0,0,2,2,2)=fct(28)
c        if (rank.eq.4) goto 511
        C5(1,1,1,1,1)=fct(29)
        C5(1,1,1,1,2)=fct(30)
        C5(1,1,1,2,2)=fct(31)
        C5(1,1,2,2,2)=fct(32)
        C5(1,2,2,2,2)=fct(33)
        C5(2,2,2,2,2)=fct(34)
        if (rank.eq.5) goto 511
c        if (rank.eq.5.and.switch.eq.0) goto 511
        C6(0,0,0,0,0,0)=fct(35)
        C6(0,0,0,0,1,1)=fct(36)
        C6(0,0,0,0,1,2)=fct(37)
        C6(0,0,0,0,2,2)=fct(38)
        C6(0,0,1,1,1,1)=fct(39)
        C6(0,0,1,1,1,2)=fct(40)
        C6(0,0,1,1,2,2)=fct(41)
        C6(0,0,1,2,2,2)=fct(42)
        C6(0,0,2,2,2,2)=fct(43)
c        if (rank.eq.5) goto 511
        C6(1,1,1,1,1,1)=fct(44)
        C6(1,1,1,1,1,2)=fct(45)
        C6(1,1,1,1,2,2)=fct(46)
        C6(1,1,1,2,2,2)=fct(47)
        C6(1,1,2,2,2,2)=fct(48)
        C6(1,2,2,2,2,2)=fct(49)
        C6(2,2,2,2,2,2)=fct(50)

 511    continue

        if(sym.eq.1) then
      
          if(rank.le.1) goto 522

          C2(2,1)=C2(1,2) 
          
          if (rank.eq.2) goto 522
          
          C3(2,2,1) = C3(1,2,2)
          C3(2,1,2) = C3(1,2,2)
          C3(2,1,1) = C3(1,1,2)
          C3(1,2,1) = C3(1,1,2)
          
          if (rank.eq.3.and.switch.eq.0) goto 522
          
          C4(0,0,2,1)=C4(0,0,1,2) 
          
          if (rank.eq.3) goto 522
          
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
          
          if (rank.eq.4.and.switch.eq.0) goto 522
          
          C5(0,0,1,2,1) = C5(0,0,1,1,2)
          C5(0,0,2,1,1) = C5(0,0,1,1,2)
          C5(0,0,2,2,1) = C5(0,0,1,2,2)
          C5(0,0,2,1,2) = C5(0,0,1,2,2)
          
          if (rank.eq.4) goto 522
          
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

          if (rank.eq.5.and.switch.eq.0) goto 522

          C6(0,0,0,0,2,1) = C6(0,0,0,0,1,2)
          C6(0,0,1,1,2,1) = C6(0,0,1,1,1,2)
          C6(0,0,1,2,1,1) = C6(0,0,1,1,1,2)
          C6(0,0,2,1,1,1) = C6(0,0,1,1,1,2)
          C6(0,0,1,2,1,2) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,1,2) = C6(0,0,1,1,2,2)
          C6(0,0,1,2,2,1) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,2,1) = C6(0,0,1,1,2,2)
          C6(0,0,2,2,1,1) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,2,2) = C6(0,0,1,2,2,2)
          C6(0,0,2,2,1,2) = C6(0,0,1,2,2,2)
          C6(0,0,2,2,2,1) = C6(0,0,1,2,2,2)

        end if
        goto 522  
      endif

c      rankacc = 3
c      ordg3acc = ordg
c      ordgy3acc = ordgy

      Ccount(0) = Ccount(0) + 1

      ordg3  = min(5-rank,ordg)
      ordgy3 = min((6-rank)/2,ordgy)
      ordgp3 = min(5-rank,ordg)
      lc = 0 

      rankacc = rank
      ordg3acc = ordg3
c changed 18.05.05
c      ordgy3acc = ordgy3
      ordgy3acc = min((5-rank)/2,ordgy)
      ordgp3acc = ordgp3

      mm02 = elimcminf2(m02)
      mm12 = elimcminf2(m12)
      mm22 = elimcminf2(m22)
      q10  = elimcminf2(p10)
      q21  = elimcminf2(p21)
      q20  = elimcminf2(p20)

      f1 = q10-mm12+mm02
      f2 = q20-mm22+mm02
      detz = -(q10*q10 + q20*q20 + q21*q21
     &     -2d0*(q10*q20 + q10*q21 + q20*q21))
      detz  = 4d0*q10*q20 - (q21-q10-q20)*(q21-q10-q20)
      if (abs(detz/(4d0*q10*q20 + (q21-q10-q20)*(q21-q10-q20)))
     &    .lt.1d-4) then
        if (abs(q10-q20).lt.abs(q10-q21).and.
     &      abs(q10-q20).lt.abs(q20-q21)) then
          detz  =  4d0*q10*q21 - (q10-q20+q21)*(q10-q20+q21)
        end if
      end if
      adetz = abs(detz)
      mp2  = (dmax1(abs(q10),abs(q20),abs(q21)))
      mm2  = (dmax1(abs(mm12),abs(mm22),abs(mm02)))
      mpm2 = dmax1(mp2,mm2)
      mf   =  dmax1(abs(f1),abs(f2))
      norm = mpm2
      rm02 = abs(mm02)
      zadjf1 = 2*q20*f1 + (q21-q10-q20)*f2
      zadjf2 = 2*q10*f2 + (q21-q10-q20)*f1
      mzadjf = dmax1(abs(zadjf1),abs(zadjf2))
      mzadjfd= dmax1(mzadjf,adetz)
      mxadj  = dmax1(abs(4d0*mm02*q10 - f1*f1),
     &            abs(2d0*mm02*(q10+q20-q21) - f1*f2),
     &            abs(4d0*mm02*q20 - f2*f2))
      adety  = abs(2d0*mm02*detz  - zadjf1*f1 - zadjf2*f2)

c        write(*,*) 'cCp12345:  test  = ',2d0*mm02*(q21-q10-q20) - f1*f2
c     &             ,mm02,q21-q20-q10,f1,f2
c        write(*,*) 'cCp12345:  mxadj  = ',abs(4d0*mm02*q10 - f1*f1),
c     &            abs(2d0*mm02*(q10+q20-q21) - f1*f2),
c     &            abs(4d0*mm02*q20 - f2*f2)

c      write(*,*) 2*q20*f1 + (q21-q10-q20)*f2,2*q10*f2 + (q21-q10-q20)*f1
c      write(*,*) 2*q20*f1, (q21-q10-q20)*f2,2*q10*f2 , (q21-q10-q20)*f1

c      write(*,*) 'Cp test ',adetz/mp2,mzadjf/mp2,mxadj/mp2


      if (mzadjf.ne.0d0) then
        accg  = adetz/mzadjf*(mpm2**2/mzadjf)**rankacc
     &      *(adetz*mpm2**2/mzadjf**2)**ordg3acc
c      write(*,*) 'cCp12345 g  ', adetz/mzadjf,
      else 
        accg = accinf
      end if
      if (mxadj*mp2.ne.0d0) then
        accgy = (mzadjfd*mpm2/(mxadj*mp2))**(ordgy3acc+1)
c     write(*,*) 'cCp12345 gy ',(mzadjfd*mpm2/(mxadj*mp2))**(ordgy3acc+1)
      else
        accgy = accinf
      end if
      if (mp2*adety.ne.0d0) then
        accgy2 = (mzadjfd*dmax1(mxadj,mzadjf)/(mp2*adety))
     &             **(ordgy3acc+1)
      else
        accgy2 = accinf
      end if
      if (mod(rankacc,2).eq.1) then
        if (adety.ne.0d0.and.adety.ne.0) then
          accy  = accB*dmax1((mzadjf/adetz)**rankacc,
     &                 mzadjf/adetz*(mxadj/adetz)**((rankacc-1)/2),
     &                 norm*mzadjf/adety,
     &                 norm*mzadjf/adety*(mzadjf/adetz)**rankacc,
     &                 norm*mxadj/adety,
     &                 norm*mxadj/adety*(mzadjf/adetz)**(rankacc-1),
     &                 norm*mxadj/adety*(mxadj/adetz)**((rankacc-1)/2)
     &      )
        else
          accy = accinf
        end if
        if (adetz.ne.0) then
          accpv =  accB*(norm*mp2/adetz)
     &        *dmax1((mzadjf/adetz)**(rankacc-1),
     &        (rm02*mp2/adetz)**((rankacc-1)/2))
        else
          accpv = accinf
        end if
      else 
        if (adety.ne.0d0.and.adetz.ne.0) then
          accy  = accB*dmax1((mzadjf/adetz)**rankacc,
     &                  (mxadj/adetz)**(rankacc/2),
     &                  norm*mzadjf/adety,
     &                  norm*mzadjf/adety*(mzadjf/adetz)**rankacc,
     &                  norm*mxadj/adety,
     &                  norm*mxadj/adety*(mzadjf/adetz)**(rankacc-1),
     &                  norm*adetz/adety*(mxadj/adetz)**(rankacc/2),
     &                  norm*mzadjf/adety*(mxadj/adetz)**(rankacc/2)
     &      )
        else
          accy = accinf
        end if
        if (rankacc.eq.0) then
          accpv = accB
        else  if (adetz.ne.0) then
          accpv =  accB*(norm*mp2/adetz)*(mzadjf/adetz)
     &        *dmax1((mzadjf/adetz)**(rankacc-2),
     &        (rm02*mp2/adetz)**((rankacc-2)/2))
        else 
          accpv = accinf
        end if
      end if
      if (mf.ne.0) then
        accgp = (mp2/mf)**(1+ordgp3acc)
      else
        accgp = accinf
      end if
c@

c      if (.true.) then
      if (.false.) then
        write(*,*)
        write(*,*) 'cCp12345: vers3 = ',version3
        write(*,*) 'cCp12345: rank3 = ',rank,ordg3,ordgy3,ordgp3
        write(*,*) 'cCp12345: ranka = ',rankacc,ordg3acc,
     &      ordgy3acc,ordgp3acc
        write(*,*) 'cCp12345: mp2   = ',mp2,mpm2
        write(*,*) 'cCp12345: norm = ',norm
        write(*,*) 'cCp12345: detz  = ',adetz,adetz/norm**2,
     &      (norm*mp2)/adetz,(rm02*mp2)/adetz,mzadjf/adetz
        write(*,*) 'cCp12345: mzadjf= ',mzadjf,mzadjf/norm**2,
     &      mpm2**2/mzadjf,adetz/mzadjf
        write(*,*) 'cCp12345: mxadj = ',mxadj,mxadj/norm**2,
     &      mpm2**2/mxadj,mzadjfd/mxadj
        if (adety.ne.0d0) then
          write(*,*) 'cCp12345:  dety = ',adety,adety/norm**3,
     &      mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &        norm*mxadj/adety    
        else
          write(*,*) 'cCp12345:  dety = ',adety
        end if
        write(*,*) 'cCp12345:  mf   = ',mf,mf/norm,mp2/mf
        write(*,*) 'cCp12345:  acc  = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
        write(*,*) 'cCp12345:  p10  = ',p10,q10
        write(*,*) 'cCp12345:  p21  = ',p21,q21
        write(*,*) 'cCp12345:  p20  = ',p20,q20
        write(*,*) 'cCp12345:  m02  = ',m02,mm02
        write(*,*) 'cCp12345:  m12  = ',m12,mm12
        write(*,*) 'cCp12345:  m22  = ',m22,mm22
      end if

      if (accgy2.le.accg.and.accgy2.le.accgy.and.accgy2.le.accy
     &      .and.accgy2.le.accgy2.and.accgy2.lt.accpv) then
        Ccount(6) = Ccount(6) + 1
      end if

      if (version3.eq.0.or.version3.eq.5) then 

      if(accpv.lt.1d-10.or.
     &  accpv.le.accg.and.accpv.le.accgy.and.accpv.le.accy
     &      .and.accpv.le.accgp.and.accpv.lt.accstop) then
c        write(*,*) 'call Cpv'   
        call cCpv12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,rank)
        lc = 1
        Ccount(1) = Ccount(1) + 1
        if (accpv.gt.accworst) then
          errwrite = .true.
          accworst = accpv
        else
          errwrite = .false.       
        end if
        if (accpv.gt.accbad1) then
          Ccount(11) = Ccount(11) + 1
          if (accpv.gt.accbad2) then
            Ccount(21) = Ccount(21) + 1
            if (accpv.gt.accbad3) then
              Ccount(31) = Ccount(31) + 1
            end if  
          end if  
        end if  
      elseif (accg.le.accpv.and.accg.le.accgy.and.accg.le.accy
     &      .and.accg.le.accgp.and.accg.lt.accstop) then
c            write(*,*) 'call Cg'
        call cCg12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,C6,rank)
        lc = 2
        Ccount(2) = Ccount(2) + 1
        if (accg.gt.accworst) then
          errwrite = .true.
          accworst = accg
        else
          errwrite = .false.       
        end if
        if (accg.gt.accbad1) then
          Ccount(12) = Ccount(12) + 1
          if (accg.gt.accbad2) then
            Ccount(22) = Ccount(22) + 1
            if (accg.gt.accbad3) then
              Ccount(32) = Ccount(32) + 1
            end if  
          end if  
        end if  
      elseif (accgy.le.accpv.and.accgy.le.accg.and.accgy.le.accy
     &      .and.accgy.le.accgp.and.accgy.lt.accstop)then
c            write(*,*) 'call Cgy'
        call cCgy12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,C6,rank)
        lc = 3
        Ccount(3) = Ccount(3) + 1
        if (accgy.gt.accworst) then
          errwrite = .true.
          accworst = accgy
        else
          errwrite = .false.       
        end if
        if (accgy.gt.accbad1) then
          Ccount(13) = Ccount(13) + 1
          if (accgy.gt.accbad2) then
            Ccount(23) = Ccount(23) + 1
            if (accgy.gt.accbad3) then
              Ccount(33) = Ccount(33) + 1
            end if  
          end if  
        end if  
      elseif (accy.le.accpv.and.accy.le.accg.and.accy.le.accgy
     &      .and.accy.le.accgp.and.accy.lt.accstop)then
c            write(*,*) 'call Cy'
        call cCy12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,rank)
        lc = 4
        Ccount(4) = Ccount(4) + 1
        if (accy.gt.accworst) then
          errwrite = .true.
          accworst = accy
        else
          errwrite = .false.       
        end if
        if (accy.gt.accbad1) then
          Ccount(14) = Ccount(14) + 1
          if (accy.gt.accbad2) then
            Ccount(24) = Ccount(24) + 1
            if (accy.gt.accbad3) then
              Ccount(34) = Ccount(34) + 1
            end if  
          end if  
        end if  
      elseif (accgp.le.accpv.and.accgp.le.accg.and.accgp.le.accgy
     &      .and.accgp.le.accy.and.accgp.lt.accstop) then
c            write(*,*) 'call Cy'
        call cCgp12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,C6,rank)
        lc = 5
        Ccount(5) = Ccount(5) + 1
        if (accgp.gt.accworst) then
          errwrite = .true.
          accworst = accy
        else
          errwrite = .false.       
        end if
        if (accgp.gt.accbad1) then
          Ccount(15) = Ccount(15) + 1
          if (accgp.gt.accbad2) then
            Ccount(25) = Ccount(25) + 1
            if (accgp.gt.accbad3) then
              Ccount(35) = Ccount(35) + 1
            end if  
          end if  
        end if  
      else
        if (accpv.lt.accy) then
          call cCpv12345(p10,p21,p20,m02,m12,m22,
     &        C0,C1,C2,C3,C4,C5,rank)
          lc = 8
          Ccount(8) = Ccount(8) + 1
          if (accpv.gt.accworst) then
            errwrite = .true.
            accworst = accy
          else
            errwrite = .false.       
          end if
          if (accpv.gt.accbad1) then
            Ccount(18) = Ccount(18) + 1
            if (accpv.gt.accbad2) then
              Ccount(28) = Ccount(28) + 1
              if (accpv.gt.accbad3) then
                Ccount(38) = Ccount(38) + 1
              end if  
            end if  
          end if  
c        write(*,*) 'call Cpv 2'   
        else
          call cCy12345(p10,p21,p20,m02,m12,m22,
     &        C0,C1,C2,C3,C4,C5,rank)
          lc = 9
          Ccount(9) = Ccount(9) + 1
          if (accy.gt.accworst) then
            errwrite = .true.
            accworst = accy
          else
            errwrite = .false.       
          end if
c            write(*,*) 'call Cy 2'
          if (accy.gt.accbad1) then
            Ccount(19) = Ccount(19) + 1
            if (accy.gt.accbad2) then
              Ccount(29) = Ccount(29) + 1
              if (accy.gt.accbad3) then
                Ccount(39) = Ccount(39) + 1
              end if  
            end if  
          end if  
        end if
        if (counterr.le.20) then
          if (counterr.eq.0) then
            write(*,*)
            write(*,*) 'cCp12345:  Cijk calculation unstable'
            write(*,*) 'cCp12345:  wrong branch'
            write(*,*) 'cCp12345:  mf   = ',mf,mp2/mf
            write(*,*) 'cCp12345:  acc  = ',accpv,accg,accgy,accy,
     &          accgp,accgy2,accstop
          end if

          write(errout,*)
          write(errout,*) 'cCp12345:  Cijk calculation unstable'
          write(errout,*) 'cCp12345:  wrong branch'
          write(errout,*) 'cCp12345: vers3 = ',version3
          write(errout,*) 'cCp12345: rank3 = ',rank,ordg3,ordgy3,ordgp3
          write(errout,*) 'cCp12345: ranka = ',rankacc,ordg3acc
     &        ,ordgy3acc,ordgp3acc
          write(errout,*) 'cCp12345: mp2   = ',mp2,mpm2
          write(errout,*) 'cCp12345: norm  = ',norm
          write(errout,*) 'cCp12345: detz  = ',adetz,adetz/norm**2,
     &        (norm*mp2)/adetz,(rm02*mp2)/adetz,mzadjf/adetz
          write(errout,*) 'cCp12345: mzadjf= ',mzadjf,mzadjf/norm**2,
     &        mpm2**2/mzadjf,adetz/mzadjf
          write(errout,*) 'cCp12345: mxadj = ',mxadj,mxadj/norm**2,
     &        mpm2**2/mxadj,mzadjfd/mxadj
          if (adety.ne.0d0) then
            write(errout,*) 'cCp12345:  dety = ',adety,adety/norm**3,
     &          mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &          norm*mxadj/adety    
          else
            write(errout,*) 'cCp12345:  dety = ',adety
          end if
          write(errout,*) 'cCp12345:  mf   = ',mf,mf/norm,mp2/mf
          write(errout,*) 'cCp12345:  acc  = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
          write(errout,*) 'cCp12345:accstop= ',accstop
          write(errout,*) 'cCp12345:  p10  = ',p10,q10
          write(errout,*) 'cCp12345:  p21  = ',p21,q21
          write(errout,*) 'cCp12345:  p20  = ',p20,q20
          write(errout,*) 'cCp12345:  m02  = ',m02,mm02
          write(errout,*) 'cCp12345:  m12  = ',m12,mm12
          write(errout,*) 'cCp12345:  m22  = ',m22,mm22
          counterr = counterr + 1
        end if
c        stop
      end if      

      if (errwrite) then
        if (flag.eq.0) then
          write(*,*)
          write(*,*) 'cCp12345:  Cijk calculation unstable'
          flag = 1
        end if

        write(errout,*)
        write(errout,*) 'cCp12345:  Cijk calculation unstable'
        write(errout,*) 'cCp12345: vers3 = ',version3
        write(errout,*) 'cCp12345: rank3 = ',rank,ordg3,ordgy3,ordgp3
        write(errout,*) 'cCp12345: ranka = ',rankacc,ordg3acc
     &      ,ordgy3acc,ordgp3acc
        write(errout,*) 'cCp12345: mp2   = ',mp2,mpm2
        write(errout,*) 'cCp12345: norm  = ',norm
        write(errout,*) 'cCp12345: detz  = ',adetz,adetz/norm**2,
     &      (norm*mp2)/adetz,(rm02*mp2)/adetz,mzadjf/adetz
        write(errout,*) 'cCp12345: mzadjf= ',mzadjf,mzadjf/norm**2,
     &      mpm2**2/mzadjf,adetz/mzadjf
        write(errout,*) 'cCp12345: mxadj = ',mxadj,mxadj/norm**2,
     &      mpm2**2/mxadj,mzadjfd/mxadj
        if (adety.ne.0d0) then
          write(errout,*) 'cCp12345:  dety = ',adety,adety/norm**3,
     &        mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &        norm*mxadj/adety    
        else
          write(errout,*) 'cCp12345:  dety = ',adety
        end if
        write(errout,*) 'cCp12345:  mf   = ',mf,mf/norm,mp2/mf
        write(errout,*) 'cCp12345:  acc  = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
        write(errout,*) 'cCp12345:accwo  = ',accworst
        write(errout,*) 'cCp12345:  p10  = ',p10,q10
        write(errout,*) 'cCp12345:  p21  = ',p21,q21
        write(errout,*) 'cCp12345:  p20  = ',p20,q20
        write(errout,*) 'cCp12345:  m02  = ',m02,mm02
        write(errout,*) 'cCp12345:  m12  = ',m12,mm12
        write(errout,*) 'cCp12345:  m22  = ',m22,mm22
      end if

      else if (version3.eq.4) then 

      if(accpv.lt.1d-10.or.
     &  accpv.le.accg.and.accpv.le.accgy.and.accpv.le.accy
     &      .and.accpv.lt.accstop) then
c        write(*,*) 'call Cpv'   
        call cCpv12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,rank)
        lc = 1
        Ccount(1) = Ccount(1) + 1
        if (accpv.gt.accworst) then
          errwrite = .true.
          accworst = accpv
        else
          errwrite = .false.       
        end if
      elseif (accg.le.accpv.and.accg.le.accgy.and.accg.le.accy
     &      .and.accg.lt.accstop) then
c            write(*,*) 'call Cg'
        call cCg12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,C6,rank)
        lc = 2
        Ccount(2) = Ccount(2) + 1
        if (accg.gt.accworst) then
          errwrite = .true.
          accworst = accg
        else
          errwrite = .false.       
        end if
      elseif (accgy.le.accpv.and.accgy.le.accg.and.accgy.le.accy
     &      .and.accgy.lt.accstop)then
c            write(*,*) 'call Cgy'
        call cCgy12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,C6,rank)
        lc = 3
        Ccount(3) = Ccount(3) + 1
        if (accgy.gt.accworst) then
          errwrite = .true.
          accworst = accgy
        else
          errwrite = .false.       
        end if
      elseif (accy.le.accpv.and.accy.le.accg.and.accy.le.accgy
     &      .and.accy.lt.accstop)then
c            write(*,*) 'call Cy'
        call cCy12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,rank)
        lc = 4
        Ccount(4) = Ccount(4) + 1
        if (accy.gt.accworst) then
          errwrite = .true.
          accworst = accy
        else
          errwrite = .false.       
        end if
      else
        if (accpv.lt.accy) then
          call cCpv12345(p10,p21,p20,m02,m12,m22,
     &        C0,C1,C2,C3,C4,C5,rank)
          lc = 11
          Ccount(11) = Ccount(11) + 1
          if (accpv.gt.accworst) then
            errwrite = .true.
            accworst = accy
          else
            errwrite = .false.       
          end if
c        write(*,*) 'call Cpv 2'   
        else
          call cCy12345(p10,p21,p20,m02,m12,m22,
     &        C0,C1,C2,C3,C4,C5,rank)
          lc = 12
          Ccount(12) = Ccount(12) + 1
          if (accy.gt.accworst) then
            errwrite = .true.
            accworst = accy
          else
            errwrite = .false.       
          end if
c            write(*,*) 'call Cy 2'
        end if
        if (counterr.le.20) then
          if (counterr.eq.0) then
            write(*,*)
            write(*,*) 'cCp12345:  Cijk calculation unstable'
            write(*,*) 'cCp12345:  wrong branch'
            write(*,*) 'cCp12345:  acc  = ',
     &          accpv,accg,accgy,accy,accstop
          end if

          write(errout,*)
          write(errout,*) 'cCp12345:  Cijk calculation unstable'
          write(errout,*) 'cCp12345:  wrong branch'
          write(errout,*) 'cCp12345: vers3 = ',version3
          write(errout,*) 'cCp12345: rank3 = ',rank,ordg3,ordgy3,ordgp3
          write(errout,*) 'cCp12345: ranka = ',rankacc,ordg3acc
     &        ,ordgy3acc,ordgp3acc
          write(errout,*) 'cCp12345: mp2   = ',mp2,mpm2
          write(errout,*) 'cCp12345: norm  = ',norm
          write(errout,*) 'cCp12345: detz  = ',adetz,adetz/norm**2,
     &        (norm*mp2)/adetz,(rm02*mp2)/adetz,mzadjf/adetz
          write(errout,*) 'cCp12345: mzadjf= ',mzadjf,mzadjf/norm**2,
     &        mpm2**2/mzadjf,adetz/mzadjf
          write(errout,*) 'cCp12345: mxadj = ',mxadj,mxadj/norm**2,
     &        mpm2**2/mxadj,mzadjfd/mxadj
          if (adety.ne.0d0) then
            write(errout,*) 'cCp12345:  dety = ',adety,adety/norm**3,
     &          mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &          norm*mxadj/adety    
          else
            write(errout,*) 'cCp12345:  dety = ',adety
          end if
          write(errout,*) 'cCp12345:  mf   = ',mf,mf/norm,mp2/mf
          write(errout,*) 'cCp12345:  acc  = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
          write(errout,*) 'cCp12345:accstop= ',accstop
          write(errout,*) 'cCp12345:  p10  = ',p10,q10
          write(errout,*) 'cCp12345:  p21  = ',p21,q21
          write(errout,*) 'cCp12345:  p20  = ',p20,q20
          write(errout,*) 'cCp12345:  m02  = ',m02,mm02
          write(errout,*) 'cCp12345:  m12  = ',m12,mm12
          write(errout,*) 'cCp12345:  m22  = ',m22,mm22
          counterr = counterr + 1
        end if
c        stop
      end if      

      if (errwrite) then
        if (flag.eq.0) then
          write(*,*)
          write(*,*) 'cCp12345:  Cijk calculation unstable'
          flag = 1
        end if

        write(errout,*)
        write(errout,*) 'cCp12345:  Cijk calculation unstable'
        write(errout,*) 'cCp12345:  vers3= ',version3
        write(errout,*) 'cCp12345: vers3 = ',version3
        write(errout,*) 'cCp12345: rank3 = ',rank,ordg3,ordgy3,ordgp3
        write(errout,*) 'cCp12345: ranka = ',rankacc,ordg3acc
     &      ,ordgy3acc,ordgp3acc
        write(errout,*) 'cCp12345: mp2   = ',mp2,mpm2
        write(errout,*) 'cCp12345: norm  = ',norm
        write(errout,*) 'cCp12345: detz  = ',adetz,adetz/norm**2,
     &      (norm*mp2)/adetz,(rm02*mp2)/adetz,mzadjf/adetz
        write(errout,*) 'cCp12345: mzadjf= ',mzadjf,mzadjf/norm**2,
     &      mpm2**2/mzadjf,adetz/mzadjf
        write(errout,*) 'cCp12345: mxadj = ',mxadj,mxadj/norm**2,
     &      mpm2**2/mxadj,mzadjfd/mxadj
        if (adety.ne.0d0) then
          write(errout,*) 'cCp12345:  dety = ',adety,adety/norm**3,
     &        mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &        norm*mxadj/adety    
        else
          write(errout,*) 'cCp12345:  dety = ',adety
        end if
        write(errout,*) 'cCp12345:  mf   = ',mf,mf/norm,mp2/mf
        write(errout,*) 'cCp12345:  acc  = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
        write(errout,*) 'cCp12345:accwo  = ',accworst
        write(errout,*) 'cCp12345:  p10  = ',p10,q10
        write(errout,*) 'cCp12345:  p21  = ',p21,q21
        write(errout,*) 'cCp12345:  p20  = ',p20,q20
        write(errout,*) 'cCp12345:  m02  = ',m02,mm02
        write(errout,*) 'cCp12345:  m12  = ',m12,mm12
        write(errout,*) 'cCp12345:  m22  = ',m22,mm22
      end if

      else if  (version3.eq.1) then 
      
        call cCpv12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,rank)
        lc = 1
        Ccount(1) = Ccount(1) + 1


      else if (version3.eq.2) then

      if(accpv.lt.1d-12.or.
     &  accpv.lt.accg.and.accpv.lt.accstop) then
c        write(*,*) 'call Cpv'   
        call cCpv12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,rank)
        lc = 1
        Ccount(1) = Ccount(1) + 1
      elseif (accg.lt.accpv.and.accg.lt.accstop) then
c        write(*,*) 'call Cg'
        call cCg12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,C6,rank)
        lc = 2
        Ccount(2) = Ccount(2) + 1
      else
        call cCpv12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,rank)
        lc = 11
        Ccount(11) = Ccount(11) + 1

        if (counterr.le.20) then
          if (counterr.eq.0) then
            write(*,*) 'cCp12345:  Cijk calculation unstable'
          end if

          write(errout,*)
          write(errout,*) 'cCp12345: Cijk calculation unstable'
          write(errout,*) 'cCp12345: vers3 = ',version3
          write(errout,*) 'cCp12345: rank  = ',rank
          write(errout,*) 'cCp12345: mp2   = ',mp2,mpm2
          write(errout,*) 'cCp12345: detz  = ',adetz,(norm*mp2)/adetz,
     &        (rm02*mp2)/adetz,mzadjf/adetz
          write(errout,*) 'cCp12345: mzadjf= ',mzadjf,mzadjf/norm**2,
     &        mpm2**2/mzadjf,adetz/mzadjf
          write(errout,*) 'cCp12345: mxadj = ',mxadj,mxadj/norm**2,
     &        mpm2**2/mxadj,mzadjfd/mxadj
          write(errout,*) 'cCp12345: acc  = ',accpv,accg,accgy
          write(errout,*) 'cCp12345: p10  = ',p10,q10
          write(errout,*) 'cCp12345: p21  = ',p21,q21
          write(errout,*) 'cCp12345: p20  = ',p20,q20
          write(errout,*) 'cCp12345: m02  = ',m02,mm02
          write(errout,*) 'cCp12345: m12  = ',m12,mm12
          write(errout,*) 'cCp12345: m22  = ',m22,mm22
c          stop
          counterr = counterr + 1
        end if

      end if

      else if (version3.eq.3) then 

      if(accpv.lt.1d-12.or.
     &  accpv.lt.accg.and.accpv.lt.accgy.and.accpv.lt.accstop) then
c        write(*,*) 'call Cpv'   
        call cCpv12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,rank)
        lc = 1
      elseif (accg.lt.accpv.and.accg.lt.accgy.and.accg.lt.accstop) then
c            write(*,*) 'call Cg'
        call cCg12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,C6,rank)
        lc = 2
      elseif (accgy.lt.accpv.and.accgy.lt.accg.and.accgy.lt.accstop)then
c            write(*,*) 'call Cy'
        call cCgy12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,C6,rank)
        lc = 3
      else
        call cCpv12345(p10,p21,p20,m02,m12,m22,
     &      C0,C1,C2,C3,C4,C5,rank)
        lc = 1
        if (counterr.le.20) then
          if (counterr.eq.0) then
            write(*,*)
            write(*,*) 'cCp12345:  Cijk calculation unstable'
          end if

          write(errout,*)
          write(errout,*) 'cCp12345: Cijk calculation unstable'
          write(errout,*) 'cCp12345: vers3= ',version3
          write(errout,*) 'cCp12345: rank  = ',rank
          write(errout,*) 'cCp12345: mp2   = ',mp2,mpm2
          write(errout,*) 'cCp12345: detz  = ',adetz,(norm*mp2)/adetz,
     &        (rm02*mp2)/adetz,mzadjf/adetz
          write(errout,*) 'cCp12345: mzadjf= ',mzadjf,mzadjf/norm**2,
     &        mpm2**2/mzadjf,adetz/mzadjf
          write(errout,*) 'cCp12345: mxadj = ',mxadj,mxadj/norm**2,
     &        mpm2**2/mxadj,mzadjfd/mxadj
          write(errout,*) 'cCp12345: acc  = ',accpv,accg,accgy
          write(errout,*) 'cCp12345: p10  = ',p10,q10
          write(errout,*) 'cCp12345: p21  = ',p21,q21
          write(errout,*) 'cCp12345: p20  = ',p20,q20
          write(errout,*) 'cCp12345: m02  = ',m02,mm02
          write(errout,*) 'cCp12345: m12  = ',m12,mm12
          write(errout,*) 'cCp12345: m22  = ',m22,mm22
          counterr = counterr + 1
        end if
c        stop
      end if      

      else
        
        write(*,*) 'cCp12345:  version3 not defined'        
        write(*,*) 'cCp12345:  version3 = ',version3        

      end if      


c write to cash

      fct(1)=C0
      
      if (rank.eq.0) goto 521
      
      fct(2)=C1(1)
      fct(3)=C1(2)
      
      if (rank.eq.1) goto 521
c      if (rank.eq.1.and.switch.eq.0) goto 521
      
      fct(4)=C2(0,0)
      
c      if (rank.eq.1) goto 521
      
      fct(5) =C2(1,1)
      fct(6) =C2(1,2)
      fct(7) =C2(2,2)
      
      if (rank.eq.2) goto 521
c      if (rank.eq.2.and.switch.eq.0) goto 521
      
      fct(8)=C3(0,0,1)
      fct(9)=C3(0,0,2)
      
c      if (rank.eq.2) goto 521
      
      fct(10)=C3(1,1,1)
      fct(11)=C3(1,1,2)
      fct(12)=C3(1,2,2)
      fct(13)=C3(2,2,2)
      
      if (rank.eq.3) goto 521
c      if (rank.eq.3.and.switch.eq.0) goto 521
      
      fct(14)=C4(0,0,0,0)
      fct(15)=C4(0,0,1,1)
      fct(16)=C4(0,0,1,2)
      fct(17)=C4(0,0,2,2)
      
c      if (rank.eq.3) goto 521
      
      fct(18)=C4(1,1,1,1)
      fct(19)=C4(1,1,1,2)
      fct(20)=C4(1,1,2,2)
      fct(21)=C4(1,2,2,2)
      fct(22)=C4(2,2,2,2)
      
      if (rank.eq.4) goto 521
c      if (rank.eq.4.and.switch.eq.0) goto 521
      
      fct(23)=C5(0,0,0,0,1)
      fct(24)=C5(0,0,0,0,2)
      fct(25)=C5(0,0,1,1,1)
      fct(26)=C5(0,0,1,1,2)
      fct(27)=C5(0,0,1,2,2)
      fct(28)=C5(0,0,2,2,2)
      
c      if (rank.eq.4) goto 521
      
      fct(29)=C5(1,1,1,1,1)
      fct(30)=C5(1,1,1,1,2)
      fct(31)=C5(1,1,1,2,2)
      fct(32)=C5(1,1,2,2,2)
      fct(33)=C5(1,2,2,2,2)
      fct(34)=C5(2,2,2,2,2)

      if (rank.eq.5) goto 521

      fct(44)= C6(1,1,1,1,1,1)
      fct(45)= C6(1,1,1,1,1,2)
      fct(46)= C6(1,1,1,1,2,2)
      fct(47)= C6(1,1,1,2,2,2)
      fct(48)= C6(1,1,2,2,2,2)
      fct(49)= C6(1,2,2,2,2,2)
      fct(50)= C6(2,2,2,2,2,2)

 521  continue

      call cachewrite(fct,50,type)


c test part

c>      if (n.eq.4.and.usecache.eq.1)  then
c>        write(83,*) 'nto = ',ncall(3),countC,sym
c>        write(83,*) 'C0 = ',C0
c>        write(83,*) 'C1 = ',C1
c>        write(83,*) 'C2 = ',C2
c>        write(83,*) 'C3 = ',C3
c>c        write(83,*) 'C4 = ',C4(0,0,0,0)
c>c        write(83,*) 'C4 = ',C4(0,0,1,1)
c>        do i1=1,2
c>        do i2=1,2
c>        do i3=1,2
c>c        write(83,133) i1,i2,i3,D3(i1,i2,i3)
c>        end do
c>        end do
c>        end do
c>      else if (n.eq.4.and.usecache.eq.0) then
c>        write(84,*) 'nto = ',countC,countC,sym
c>        write(84,*) 'C0 = ',C0
c>        write(84,*) 'C1 = ',C1
c>        write(84,*) 'C2 = ',C2
c>        write(84,*) 'C3 = ',C3
c>c        write(84,*) 'C4 = ',C4(0,0,0,0)
c>c        write(84,*) 'C4 = ',C4(0,0,1,1)
c>        do i1=1,2
c>        do i2=1,2
c>        do i3=1,2
c>c        write(84,133) i1,i2,i3,C3(i1,i2,i3)
c>        end do
c>        end do
c>        end do
c>      end if

c      write(*,*) 'cCp12345 lc = ',lc
 
      atest = 0
      if (ltest.gt.1.and.counttest.lt.100) then
        call cachetempoff
        atest = 0
        counttest = counttest + 1
        if (lc.eq.1.and.accpv.gt.acctest) then
          if (accg.lt.accgy.and.accg.lt.accy) then
            call cCg12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,Ct6,rank)
          else if (accgy.lt.accg.and.accgy.lt.accy) then
            call cCgy12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,Ct6,rank)
          else 
            call cCy12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,rank)
          end if
          atest = 1 
        else if (lc.eq.2.and.accg.gt.acctest) then
          if (accpv.lt.accgy.and.accpv.lt.accy) then
            call cCpv12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,rank)
          else if (accgy.lt.accpv.and.accgy.lt.accy) then
            call cCgy12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,Ct6,rank)
          else 
            call cCy12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,rank)
          end if
          atest = 1 
        else if (lc.eq.3.and.accgy.gt.acctest) then
          if (accg.lt.accpv.and.accg.lt.accy) then
            call cCg12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,Ct6,rank)
          else if (accpv.lt.accg.and.accpv.lt.accy) then
            call cCpv12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,rank)
          else 
            call cCy12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,rank)
          end if
          atest = 1 
        else if (lc.eq.4.and.accy.gt.acctest) then
          if (accg.lt.accgy.and.accg.lt.accpv) then
            call cCg12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,Ct6,rank)
          else if (accgy.lt.accg.and.accgy.lt.accpv) then
            call cCgy12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,Ct6,rank)
          else 
            call cCpv12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,rank)
          end if
          atest = 1 
        end if        
    
c        write(*,*) 'atest ',atest
          
        if (ltest.eq.20) then 
          if(accpv.lt.1d-10.or.
     &        accpv.lt.accg.and.accpv.lt.accgy.and.accpv.lt.accy
     &        .and.accpv.lt.accstop) then
c     write(*,*) 'call Cpv'   
            call cCpv12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,rank)
          elseif (accg.lt.accpv.and.accg.lt.accgy.and.accg.lt.accy
     &          .and.accg.lt.accstop) then
c     write(*,*) 'call Cg'
            call cCg12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,Ct6,rank)
          elseif (accgy.lt.accpv.and.accgy.lt.accg.and.accgy.lt.accy
     &          .and.accgy.lt.accstop)then
c     write(*,*) 'call Cgy'
            call cCgy12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,Ct6,rank)
          elseif (accy.lt.accpv.and.accy.lt.accg.and.accy.lt.accgy
     &          .and.accy.lt.accstop)then
c     write(*,*) 'call Cy'
            call cCy12345(p10,p21,p20,m02,m12,m22,
     &          Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,rank)
          else
            if (accpv.lt.accy) then
              call cCpv12345(p10,p21,p20,m02,m12,m22,
     &            Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,rank)
c     write(*,*) 'call Cpv 2'   
            else
              call cCy12345(p10,p21,p20,m02,m12,m22,
     &            Ct0,Ct1,Ct2,Ct3,Ct4,Ct5,rank)
c     write(*,*) 'call Cy 2'
            end if
          end if
          atest = 1
        end if      
 

        if (atest.eq.1) then
        write(testout,*) 
        write(testout,*) 'cCp12345:  p10  = ',p10,q10
        write(testout,*) 'cCp12345:  p21  = ',p21,q21
        write(testout,*) 'cCp12345:  p20  = ',p20,q20
        write(testout,*) 'cCp12345:  m02  = ',m02,mm02
        write(testout,*) 'cCp12345:  m12  = ',m12,mm12
        write(testout,*) 'cCp12345:  m22  = ',m22,mm22
        write(testout,*) 'cCp12345:  rank = ',rank
        write(testout,*) 'cCp12345:  mp2  = ',mp2,mpm2
        write(testout,*) 'cCp12345:  norm = ',norm
        write(testout,*) 'cCp12345: detz  = ',adetz,adetz/norm**2,
     &      (norm*mp2)/adetz,(rm02*mp2)/adetz,mzadjf/adetz
        write(testout,*) 'cCp12345: mzadjf= ',mzadjf,mzadjf/norm**2,
     &      mpm2**2/mzadjf,adetz/mzadjf
        write(testout,*) 'cCp12345: mxadj = ',mxadj,mxadj/norm**2,
     &      mpm2**2/mxadj,mzadjfd/mxadj
        if (adety.ne.0d0) then
          write(testout,*) 'cCp12345:  dety = ',adety,adety/norm**3,
     &        mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &        norm*mxadj/adety    
        else
          write(testout,*) 'cCp12345:  dety = ',adety
        end if
        write(testout,*) 'cCp12345:  mf   = ',mf,mf/norm,mp2/mf
        write(testout,*) 'cCp12345:  acc  = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
        if(cdabs(Ct0/C0-1D0).gt.acc) then
          write(testout,441) C0
          write(testout,541) Ct0
          write(testout,99) cdabs(Ct0/C0-1D0)
        end if
 99    format(' deviation     = ',G20.14)
 441    format(' C0         = ',G20.14,' + i* ',G20.14)
 541    format(' C0         = ',G20.14,' + i* ',G20.14)
        if (rank.gt.1) then
        if(cdabs(Ct2(0,0)/C2(0,0)-1D0).gt.acc) then
          write(testout,446) C2(0,0)
          write(testout,546) Ct2(0,0)
          write(testout,99) cdabs(Ct2(0,0)/C2(0,0)-1D0)
        end if
 446    format(' C2(0,0)    = ',G20.14,' + i* ',G20.14)
 546    format(' C2(0,0)    = ',G20.14,' + i* ',G20.14)
        if (rank.gt.3) then      
        if(cdabs(Ct4(0,0,0,0)/C4(0,0,0,0)-1D0).gt.acc) then
          write(testout,448) C4(0,0,0,0)
          write(testout,548) Ct4(0,0,0,0)
          write(testout,99) cdabs(Ct4(0,0,0,0)/C4(0,0,0,0)-1D0)
        end if
 448    format(' C4(0,0,0,0)= ',G20.14,' + i* ',G20.14)
 548    format(' C4(0,0,0,0)= ',G20.14,' + i* ',G20.14)
      end if 
      end if 
      if (rank.gt.0) then
        do 101 i1=1,2
          if(cdabs(Ct1(i1)/C1(i1)-1D0).gt.acc) then
            write(testout,141) i1,C1(i1)
            write(testout,241) i1,Ct1(i1)
            write(testout,99) cdabs(Ct1(i1)/C1(i1)-1D0)
          end if
 141      format(' C1(',i1,')       = ',G20.14,' + i* ',G20.14)
 241      format(' C1(',i1,')       = ',G20.14,' + i* ',G20.14)
          if (rank.gt.2) then      
            if(cdabs(Ct3(0,0,i1)/C3(0,0,i1)-1D0).gt.acc) then
              write(testout,131) i1,C3(0,0,i1)
              write(testout,231) i1,Ct3(0,0,i1)
              write(testout,99) cdabs(Ct3(0,0,i1)/C3(0,0,i1)-1D0)
            end if
 131        format(' C3(0,0,',i1,')   = ',G20.14,' + i* ',G20.14)
 231        format(' C3(0,0,',i1,')   = ',G20.14,' + i* ',G20.14)
            if (rank.gt.4) then      
              if(cdabs(Ct5(0,0,0,0,i1)/C5(0,0,0,0,i1)-1D0).gt.acc) then
                write(testout,151) i1,C5(0,0,0,0,i1)
                write(testout,251) i1,Ct5(0,0,0,0,i1)
                write(testout,99) 
     &              cdabs(Ct5(0,0,0,0,i1)/C5(0,0,0,0,i1)-1D0)
              end if
 151        format(' C5(0,0,0,0,',i1,')   = ',G20.14,' + i* ',G20.14)
 251        format(' C5(0,0,0,0,',i1,')   = ',G20.14,' + i* ',G20.14)
            end if
          end if
          if (rank.gt.1) then
          do 102 i2=i1,2
            if(cdabs(Ct2(i1,i2)/C2(i1,i2)-1D0).gt.acc) then
              write(testout,142) i1,i2,C2(i1,i2)
              write(testout,242) i1,i2,Ct2(i1,i2)
              write(testout,99) cdabs(Ct2(i1,i2)/C2(i1,i2)-1D0)
            end if
 142        format(' C2(',i1,',',i1,')     = ',G20.14,' + i* ',G20.14)
 242        format(' C2(',i1,',',i1,')     = ',G20.14,' + i* ',G20.14)
            if (rank.gt.3) then
            if(cdabs(Ct4(0,0,i1,i2)/C4(0,0,i1,i2)-1D0).gt.acc) then
              write(testout,342) i1,i2,C4(0,0,i1,i2)
              write(testout,442) i1,i2,Ct4(0,0,i1,i2)
              write(testout,99) cdabs(Ct4(0,0,i1,i2)/C4(0,0,i1,i2)-1D0)
            end if
 342        format(' C4(0,0,',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 442        format(' C4(0,0,',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
          end if
          if (rank.gt.2) then
          do 103 i3=i2,2
            if(cdabs(Ct3(i1,i2,i3)/C3(i1,i2,i3)-1D0).gt.acc) then
              write(testout,143) i1,i2,i3,C3(i1,i2,i3)
              write(testout,243) i1,i2,i3,Ct3(i1,i2,i3)
              write(testout,99) cdabs(Ct3(i1,i2,i3)/C3(i1,i2,i3)-1D0)
            end if
 143    format(' C3(',i1,',',i1,',',i1,')   = ',G20.14,' + i* ',G20.14)
 243    format(' C3(',i1,',',i1,',',i1,')   = ',G20.14,' + i* ',G20.14)
            if (rank.gt.4) then 
              if(cdabs(Ct3(i1,i2,i3)/C3(i1,i2,i3)-1D0).gt.acc) then
                write(testout,153) i1,i2,i3,C5(0,0,i1,i2,i3)
                write(testout,253) i1,i2,i3,Ct5(0,0,i1,i2,i3)
                write(testout,99) 
     &              cdabs(Ct5(0,0,i1,i2,i3)/C5(0,0,i1,i2,i3)-1D0)
              end if
            end if
 153    format(' C5(0,0,',i1,',',i1,',',i1,')   = ',
     &   G20.14,' + i* ',G20.14)
 253    format(' C5(0,0,',i1,',',i1,',',i1,')   = ',
     &   G20.14,' + i* ',G20.14)
            if (rank.gt.3) then
            do 104 i4=i3,2
              if(cdabs(Ct4(i1,i2,i3,i4)/C4(i1,i2,i3,i4)-1D0).gt.acc)then
                write(testout,144) i1,i2,i3,i4,C4(i1,i2,i3,i4)
                write(testout,244) i1,i2,i3,i4,Ct4(i1,i2,i3,i4)
                write(testout,99)cdabs(Ct4(i1,i2,i3,i4)/C4(i1,i2,i3,i4)
     &              -1D0)
              end if
 144          format(' C4(',i1,',',i1,',',i1,',',i1,') = ',G20.14,
     &            ' + i* ',G20.14)
 244          format(' C4(',i1,',',i1,',',i1,',',i1,') = ',G20.14,
     &            ' + i* ',G20.14)
              if (rank.gt.4) then
              do 105 i5=i3,2
                if(cdabs(Ct5(i1,i2,i3,i4,i5)/C5(i1,i2,i3,i4,i5)-1D0)
     &              .gt.acc)then
                write(testout,155) i1,i2,i3,i4,i5,C5(i1,i2,i3,i4,i5)
                write(testout,255) i1,i2,i3,i4,i5,Ct5(i1,i2,i3,i4,i5)
                write(testout,99)
     &              cdabs(Ct5(i1,i2,i3,i4,i5)/C5(i1,i2,i3,i4,i5) -1D0)
              end if
 155          format(' C5(',i1,',',i1,',',i1,',',i1,',',i1,') = ',
     &            G20.14,' + i* ',G20.14)
 255          format(' C5(',i1,',',i1,',',i1,',',i1,',',i1,') = ',
     &            G20.14,' + i* ',G20.14)
 105          continue
              end if
 104        continue
            end if
 103      continue
          end if
 102    continue
        end if
 101  continue
      end if
      end if
      call cachereon
      end if

 522  continue

      end
                        
***********************************************************************
      subroutine cCy12345(p10,p21,p20,m02,m12,m22,
     &    C0,C1,C2,C3,C4,C5,rank)
***********************************************************************
*     3-point tensor coefficient functions C0,C1,C2,C3,C4,C5          *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     alternative reduction via Cayley matrix                         *
*     definition of integrals:                                        *
*     {C0,C1,C2,C3,C4,C5}=                                            *
*        1/(i pi^2)*\int d^4 q {1,q,q.q,q.q.q,q.q.q.q,q.q.q.q.q}      *
*        1/[q**2 - m02]*[(q+k1)**2-m12]*[(q+k2)**2 - m22]*            *
*        [(q+k3)**2-m32] ]                                            *
*     arguments:   p10=k1**2, p21=(k2-k1)**2, p32=(k3-k2)**2          *
*                      p30=k3**2, p31=(k3-k1)**2, p20=k2**2           *
*     definition of tensor coefficients:                              *
*          C1 = \sum_{i=1}^3 C1(i) *ki                                *
*          C2 = \sum_{i,j=1}^3 C2(i,j) *ki*kj  + C2(0,0) * g          *
*          C3 = \sum_{i,j,k=1}^3 C3(i,j,k) *ki*kj*kk                  *
*              + \sum_{i=1}^3 C3(0,0,i) * (ki*g + 2*sym)              *
*          C4 = \sum_{i,j,k,l=1}^3 C4(i,j,k) *ki*kj*kk*kl             *
*              + \sum_{i,j=1}^3 C4(0,0,i,j) * (ki*kj*g + 5*sym)       *
*              +   C4(0,0,0,0) * (g*g + 2*sym)                        *
*          C5 = \sum_{i,j,k,l,m=1}^3 C4(i,j,k) *ki*kj*kk*kl*km        *
*              + \sum_{i,j,k=1}^3 C4(0,0,i,j,k) * (ki*kj*kk*g +10*sym)*
*              +   \sum_{i=1}^3 C4(0,0,0,0,i) * (g*g*ki + 14*sym)     *
*          B1w0 =  \sum_{i=1}^2 B1w0(i) *( k(i+1)-k1 )                *
*     3-point functions resulting from 4-point functions with         *
*     i-th denominator cancelled: Bxwi                                *
***********************************************************************
*     30.12.04 Ansgar Cenner    last changed  24.01.06                *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p20
      complex*16 q10,q21,q20
      complex*16 m02,m12,m22
      complex*16 mm02,mm12,mm22
      complex*16 mx(0:2,0:2),mxinv(0:2,0:2)
c      complex*16 test(0:2,0:2)
      complex*16 B0w2,B1w2,B2w2(0:1,0:1),B3w2(0:1,0:1,0:1)
     &    ,B4w2(0:1,0:1,0:1,0:1),B5w2(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w2(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0w1,B1w1,B2w1(0:1,0:1),B3w1(0:1,0:1,0:1)
     &    ,B4w1(0:1,0:1,0:1,0:1),B5w1(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w1(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0w0,B1w0,B2w0(0:1,0:1),B3w0(0:1,0:1,0:1)
     &    ,B4w0(0:1,0:1,0:1,0:1),B5w0(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w0(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0m(0:2),B1m(0:2,2),B2m(0:2,0:2,0:2)
     &    ,B3m(0:2,0:2,0:2,0:2),B4m(0:2,0:2,0:2,0:2,0:2)
c     &    ,B5m(0:2,0:2,0:2,0:2,0:2,0:2),B6m(0:2,0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 S1hat(2),S2hat(2,2),S3hat(2,0:2,0:2),
     &    S4hat(2,0:2,0:2,2),S5hat(2,0:2,0:2,0:2,0:2)
      complex*16 S2mod(2,2),S3mod(2,0:2,0:2),
     &    S4mod(2,0:2,0:2,2),S5mod(2,0:2,0:2,0:2,0:2)
c      complex*16 S00,S001(2),S002(0:2,0:2),S003(0:2,0:2,0:2)
      complex*16 C200mod,C300mod(2),C400mod(2,2),C500mod(2,2,2)
      complex*16 C600mod(2,2,2,2)
      complex*16 C0,C1(2),C2(0:2,0:2),C3(0:2,0:2,0:2)
      complex*16 C4(0:2,0:2,0:2,0:2),C5(0:2,0:2,0:2,0:2,0:2)
      complex*16 cC0f
      complex*16 elimcminf2
c      real*8     maxzadjf,maxzadj
      integer    rank,rankB
      integer    i1,i2,i3,i4,i5,j,k
      integer    ltest,sym
c      integer    testout

      common /ltest/  ltest
      common /sym/ sym

      data       B1m /6*0d0/, B2m /27*0d0/, B3m/81*0d0/, B4m/243*0d0/
c     &         , B5m/729*0d0/
c      data       testout /35/

c      write(testout,*)'Cy12345 in',p10,p21,p20,m02,m12,m22
c     &    ,rank

      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cCy12345'
        write(*,*) 'rank = ',rank
        stop
      else if (rank.gt.5) then
        write(*,*) 'rank > 5 not implemented in cCy12345'
        write(*,*) 'rank = ',rank
        stop
      end if

c      write(*,*) 'cCy12345 not yet finished!'
c      stop

      C0 =cC0f(p10,p21,p20,m02,m12,m22)
c      write(testout,01) C0

      rankB = max(rank-1,0)

c  cBp must be called also for rank = 0 because of cache!

c     write(*,*) 'cCp1',p21,p32,p31,m12,m22,m32,rankB
      call cBp12345(p21,m12,m22,B0w0,B1w0,B2w0,B3w0,B4w0,
     &    B5w0,B6w0,rankB,0)
c     write(*,*) 'cCp2',p21,p32,p31,m12,m22,m32,rankB
      call cBp12345(p20,m02,m22,B0w1,B1w1,B2w1,B3w1,B4w1,
     &    B5w1,B6w1,rankB,0)
c     write(*,*) 'cCp3',p21,p32,p31,m12,m22,m32,rankB
      call cBp12345(p10,m02,m12,B0w2,B1w2,B2w2,B3w2,B4w2,
     &    B5w2,B6w2,rankB,0)
      
      if (rank.eq.0) goto 100

      mm02 = elimcminf2(m02)
      mm12 = elimcminf2(m12)
      mm22 = elimcminf2(m22)
      q10  = elimcminf2(p10)
      q21  = elimcminf2(p21)
      q20  = elimcminf2(p20)

c>      do i1=1,2
c>        C1(i1) = 0d0
c>        C3(0,0,i1) = 0d0
c>        C5(0,0,0,0,i1) = 0d0
c>        do i2=1,2
c>          C2(i1,i2) = 0d0
c>          C4(0,0,i1,i2) = 0d0
c>          do i3=1,2
c>            C3(i1,i2,i3) = 0d0
c>            C5(0,0,i1,i2,i3) = 0d0
c>            do i4=1,2
c>              C4(i1,i2,i3,i4) = 0d0
c>              do i5=1,2
c>                C5(i1,i2,i3,i4,i5) = 0d0
c>              end do
c>            end do
c>          end do
c>        end do
c>      end do
      
      mx(0,0) = 2d0*mm02
      mx(1,0) = q10 - mm12 + mm02
      mx(2,0) = q20 - mm22 + mm02
      mx(0,1) = mx(1,0)
      mx(0,2) = mx(2,0)

      mx(1,1) = 2d0*q10
      mx(2,2) = 2d0*q20
      mx(1,2) = q10+q20-q21
      mx(2,1) = mx(1,2)

c      do k=0,2
c      do j=0,2
c      write(*,*) 'mx =',k,j,mx(k,j)
c      end do
c      end do

      call chinv(3,mx,mxinv)

c>      write(*,*) 'x00 = ',mxinv(0,0)
c>      do k=1,2
c>      write(*,*) 'x0i = ',k,mxinv(0,k)
c>      write(*,*) 'x0i = ',k,mxinv(k,0)
c>      do j=1,2
c>      write(*,*) 'x0i = ',k,j,mxinv(k,j)
c>      end do
c>      end do

      B0m(0) = B0w0
      B0m(1) = B0w1
      B0m(2) = B0w2
      
      do j =1,2
        S1hat(j) = B0m(j) - B0m(0)
      end do

      if (rank.eq.0) goto 99
  
      B1m(0,2) = B1w0
      B1m(0,1) = -B0m(0) - B1m(0,2)
      B1m(1,2) = B1w1
      B1m(2,1) = B1w2
      
      do i1=1,2
        do j =1,2
          S2hat(j,i1) = B1m(j,i1) - B1m(0,i1)
        end do
      end do

      if (rank.eq.1) goto 99

      B2m(0,0,0) = B2w0(0,0)
      B2m(1,0,0) = B2w1(0,0)
      B2m(2,0,0) = B2w2(0,0)
        
      B2m(0,2,2) = B2w0(1,1)
      B2m(0,1,2) = -B1m(0,2) - B2m(0,2,2)
      B2m(0,1,1) = -B1m(0,1) - B2m(0,1,2)
      B2m(1,2,2) = B2w1(1,1)
      B2m(2,1,1) = B2w2(1,1)
      
c      B2m(0,2,1) = B2m(0,1,2) 

      do i1=1,2
        S3hat(i1,0,0) = B2m(i1,0,0) - B2m(0,0,0)
        do i2=i1,2
          do j =1,2
            S3hat(j,i1,i2) = B2m(j,i1,i2) - B2m(0,i1,i2)
          end do
        end do
      end do
      
      if (rank.eq.2) goto 99

      B3m(0,0,0,1) = -B2w0(0,0) - B3w0(0,0,1)
      B3m(0,0,0,2) =  B3w0(0,0,1)
      B3m(1,0,0,2) =  B3w1(0,0,1)
      B3m(2,0,0,1) =  B3w2(0,0,1)
      B3m(0,2,2,2) =  B3w0(1,1,1)
      
      B3m(0,1,2,2) = -B2m(0,2,2) - B3m(0,2,2,2)
      B3m(0,1,1,2) = -B2m(0,1,2) - B3m(0,1,2,2)
      B3m(0,1,1,1) = -B2m(0,1,1) - B3m(0,1,1,2)
      
      B3m(1,2,2,2) =  B3w1(1,1,1)
      B3m(2,1,1,1) =  B3w2(1,1,1)
      
c      B3m(0,1,2,1) =  B3m(0,1,1,2)
c      B3m(0,2,1,1) =  B3m(0,1,1,2)
c      B3m(0,2,2,1) =  B3m(0,1,2,2)
c      B3m(0,2,1,2) =  B3m(0,1,2,2)
      
      do j=1,2
        do i1=1,2
          S4hat(j,0,0,i1) = B3m(j,0,0,i1) - B3m(0,0,0,i1)
          do i2=i1,2
            do i3 =i2,2
              S4hat(j,i1,i2,i3) = B3m(j,i1,i2,i3) - B3m(0,i1,i2,i3)
            end do
          end do
        end do
      end do        

      if (rank.eq.3) goto 99

      B4m(0,0,0,0,0) = B4w0(0,0,0,0)
      B4m(1,0,0,0,0) = B4w1(0,0,0,0)
      B4m(2,0,0,0,0) = B4w2(0,0,0,0)

      B4m(0,0,0,2,2) = B4w0(0,0,1,1)
      B4m(0,0,0,1,2) = -B3m(0,0,0,2)-B4m(0,0,0,2,2)
      B4m(0,0,0,1,1) = -B3m(0,0,0,1)-B4m(0,0,0,1,2)
      
      B4m(1,0,0,2,2) = B4w1(0,0,1,1)
      B4m(2,0,0,1,1) = B4w2(0,0,1,1)

      B4m(0,2,2,2,2) = B4w0(1,1,1,1)
      B4m(0,1,2,2,2) = -B3m(0,2,2,2) - B4m(0,2,2,2,2)
      B4m(0,1,1,2,2) = -B3m(0,1,2,2) - B4m(0,1,2,2,2)
      B4m(0,1,1,1,2) = -B3m(0,1,1,2) - B4m(0,1,1,2,2)
      B4m(0,1,1,1,1) = -B3m(0,1,1,1) - B4m(0,1,1,1,2)

      B4m(1,2,2,2,2) = B4w1(1,1,1,1)
      B4m(2,1,1,1,1) = B4w2(1,1,1,1)
      
c>      do i1=1,2
c>        do i2=i1,2
c>          B4m(0,0,0,i2,i1) = B4m(0,0,0,i1,i2)
c>          B4m(1,0,0,i2,i1) = B4m(1,0,0,i1,i2)
c>          B4m(2,0,0,i2,i1) = B4m(2,0,0,i1,i2)
c>          do i3=i2,2
c>            do i4=i3,2
c>              do l=0,2
c>                B4m(l,i1,i3,i2,i4) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i2,i1,i3,i4) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i2,i3,i1,i4) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i3,i1,i2,i4) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i3,i2,i1,i4) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i1,i2,i4,i3) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i1,i4,i2,i3) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i2,i1,i4,i3) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i2,i4,i1,i3) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i4,i1,i2,i3) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i4,i2,i1,i3) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i1,i3,i4,i2) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i1,i4,i3,i2) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i3,i1,i4,i2) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i3,i4,i1,i2) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i4,i1,i3,i2) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i4,i3,i1,i2) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i3,i2,i4,i1) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i3,i4,i2,i1) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i2,i3,i4,i1) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i2,i4,i3,i1) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i4,i3,i2,i1) = B4m(l,i1,i2,i3,i4)
c>                B4m(l,i4,i2,i3,i1) = B4m(l,i1,i2,i3,i4)        
c>              end do
c>            end do
c>          end do
c>        end do
c>      end do

      do k =1,2
        S5hat(k,0,0,0,0) = B4m(k,0,0,0,0) - B4m(0,0,0,0,0)
        do i1=1,2
        do i2=i1,2
        do i3=i2,2
          S5hat(k,0,0,i1,i2) = B4m(k,0,0,i1,i2) - B4m(0,0,0,i1,i2)
        do i4=i3,2
          S5hat(k,i1,i2,i3,i4) = B4m(k,i1,i2,i3,i4) -
     &        B4m(0,i1,i2,i3,i4)
        end do
        end do
        end do
        end do
      end do
      
 99   continue

 01   format(' C0y = ',G20.14,' + i* ',G20.14)
 11   format(' C1y(',i1,') = ',G20.14,' + i* ',G20.14)
 20   format(' C2y(0,0) = ',G20.14,' + i* ',G20.14)
 22   format(' C2y(',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 31   format(' C3y(0,0,',i1,') = ',G20.14,' + i* ',G20.14)
 33   format(' C3y(',i1,',',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 40   format(' C4y(0,0,0,0) = ',G20.14,' + i* ',G20.14)
 42   format(' C4y(0,0,',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 44   format(' C4y(',i1,',',i1,',',i1,',',i1,') = ',
     &  G20.14,' + i* ',G20.14)
 51   format(' C5y(0,0,0,0,',i1,') = ',G20.14,' + i* ',G20.14)
 53   format(' C5y(0,0,',i1,',',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 55   format(' C5y(',i1,',',i1,',',i1,',',i1,',',i1,') = ',
     &  G20.14,' + i* ',G20.14)

      C200mod =  (C0 
     &         - mxinv(0,1)*S1hat(1) 
     &         - mxinv(0,2)*S1hat(2) )
     &           /mxinv(0,0)

      do i1 = 1,2
        C1(i1) = mxinv(i1,0)*C200mod 
     &         + mxinv(i1,1)*S1hat(1)
     &         + mxinv(i1,2)*S1hat(2)
c      write(testout,11) i1,C1(i1)
      end do
        
      if (rank.eq.1) goto 100

      C2(0,0) = (C200mod + B0m(0) + 1d0)/4d0
c      write(testout,20) C2(0,0)

      do i1=1,2
        do j =1,2
          S2mod(j,i1) = S2hat(j,i1)
        end do
        S2mod(i1,i1) = S2mod(i1,i1) - 2d0*C2(0,0)
      end do

      do i1 = 1,2
        C300mod(i1) =  (C1(i1) 
     &         - mxinv(0,1)*S2mod(1,i1) 
     &         - mxinv(0,2)*S2mod(2,i1) )
     &           /mxinv(0,0)
      end do

      do i1 = 1,2
      do i2 = i1,2
        C2(i1,i2) = mxinv(i1,0)*C300mod(i2) 
     &         + mxinv(i1,1)*S2mod(1,i2)
     &         + mxinv(i1,2)*S2mod(2,i2)
c        write(testout,22) i1,i2,C2(i1,i2)
      end do
      end do
      if (sym.eq.1) then
        C2(2,1) = C2(1,2)
      end if
        
      if (rank.eq.2) goto 100

      do i1=1,2
        C3(0,0,i1) = (C300mod(i1) + B1m(0,i1) - 1d0/3d0)/6d0
c        write(testout,31) i1,C3(0,0,i1)
      end do
      
      do i1=1,2
        do i2=i1,2
          do j =1,2
            S3mod(j,i1,i2) = S3hat(j,i1,i2)
          end do
          S3mod(i1,i1,i2) = S3mod(i1,i1,i2) - 2d0*C3(0,0,i2)
          S3mod(i2,i1,i2) = S3mod(i2,i1,i2) - 2d0*C3(0,0,i1)
        end do
      end do

      do i1 = 1,2
      do i2 = i1,2
        C400mod(i1,i2) =  (C2(i1,i2) 
     &         - mxinv(0,1)*S3mod(1,i1,i2) 
     &         - mxinv(0,2)*S3mod(2,i1,i2) )
     &           /mxinv(0,0)
      end do
      end do

      do i1 = 1,2
      do i2 = i1,2
      do i3 = i2,2
        C3(i1,i2,i3) = mxinv(i1,0)*C400mod(i2,i3) 
     &         + mxinv(i1,1)*S3mod(1,i2,i3)
     &         + mxinv(i1,2)*S3mod(2,i2,i3)
c        write(testout,33) i1,i2,i3,C3(i1,i2,i3)
      end do
      end do
      end do
      if (sym.eq.1) then
        C3(2,2,1) = C3(1,2,2)
        C3(2,1,2) = C3(1,2,2)
        C3(2,1,1) = C3(1,1,2)
        C3(1,2,1) = C3(1,1,2)
      end if

      if (rank.eq.3) goto 100

      C4(0,0,0,0) = (
     &    (C2(0,0) 
     &         - mxinv(0,1)*S3hat(1,0,0)  
     &         - mxinv(0,2)*S3hat(2,0,0) 
     &           )/mxinv(0,0)
     &    + B2m(0,0,0) 
     &    + 1d0/24d0*(4d0*(mm02+mm12+mm22)-(q10+q20+q21)))/8d0
c      write(testout,40) C4(0,0,0,0)

      do i1=1,2
      do i2=i1,2
        C4(0,0,i1,i2) = (C400mod(i1,i2) + B2m(0,i1,i2) + 1d0/12d0)/8d0
      end do
        C4(0,0,i1,i1) = C4(0,0,i1,i1) + 1d0/96d0
      end do
      if (sym.eq.1) then
        C4(0,0,2,1) = C4(0,0,1,2)
      end if

c      write(testout,42) ((i1,i2,C4(0,0,i1,i2),i1=1,2),i2=1,2)

      do i1=1,2
      do i2=i1,2
      do i3=i2,2
        do j =1,2
          S4mod(j,i1,i2,i3) = S4hat(j,i1,i2,i3)
        end do
        S4mod(i1,i1,i2,i3) = S4mod(i1,i1,i2,i3) -2d0*C4(0,0,i2,i3)
        S4mod(i2,i1,i2,i3) = S4mod(i2,i1,i2,i3) -2d0*C4(0,0,i1,i3)
        S4mod(i3,i1,i2,i3) = S4mod(i3,i1,i2,i3) -2d0*C4(0,0,i1,i2)
      end do
      end do
      end do
        
      do i1 = 1,2
      do i2 = i1,2
      do i3 = i2,2
        C500mod(i1,i2,i3) =  (C3(i1,i2,i3) 
     &         - mxinv(0,1)*S4mod(1,i1,i2,i3) 
     &         - mxinv(0,2)*S4mod(2,i1,i2,i3) )
     &           /mxinv(0,0)
      end do
      end do
      end do

      do i1 = 1,2
      do i2 = i1,2
      do i3 = i2,2
      do i4 = i3,2
        C4(i1,i2,i3,i4) = mxinv(i1,0)*C500mod(i2,i3,i4) 
     &         + mxinv(i1,1)*S4mod(1,i2,i3,i4)
     &         + mxinv(i1,2)*S4mod(2,i2,i3,i4)
c        write(testout,44) i1,i2,i3,i4,C4(i1,i2,i3,i4)
      end do
      end do
      end do
      end do
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

      if (rank.eq.4) goto 100

      do i1=1,2
        do j =1,2
          S4mod(j,0,0,i1) = S4hat(j,0,0,i1)
        end do
        S4mod(i1,0,0,i1) = S4mod(i1,0,0,i1) - 2d0*C4(0,0,0,0)
      end do

      do i1 = 1,2
        C5(0,0,0,0,i1) = (
     &    (C3(0,0,i1) 
     &         - mxinv(0,1)*S4mod(1,0,0,i1)  
     &         - mxinv(0,2)*S4mod(2,0,0,i1) 
     &           )/mxinv(0,0)
     &    - 1d0/120d0*(5d0*(mm02+mm12+mm22)-(q10+q20+2d0*q21))
     &    + B3m(0,0,0,i1) 
     &    )/10d0
      end do
      C5(0,0,0,0,1) = C5(0,0,0,0,1) - (5d0*mm12-q10)/1200d0
      C5(0,0,0,0,2) = C5(0,0,0,0,2) - (5d0*mm22-q20)/1200d0
 
c      write(testout,51) (i1,C5(0,0,0,0,i1),i1=1,2)
      
      do i1=1,2
      do i2=i1,2
      do i3=i2,2
        C5(0,0,i1,i2,i3) = (C500mod(i1,i2,i3) + B3m(0,i1,i2,i3)
     &                       - 1d0/30d0)/10d0
      end do
      end do
        C5(0,0,i1,i1,i1) = C5(0,0,i1,i1,i1) - 1d0/150D0
      end do
      if (sym.eq.1) then
        C5(0,0,1,2,1) = C5(0,0,1,1,2)
        C5(0,0,2,1,1) = C5(0,0,1,1,2)
        C5(0,0,2,2,1) = C5(0,0,1,2,2)
        C5(0,0,2,1,2) = C5(0,0,1,2,2)
      end if

c      write(testout,53) (((i1,i2,i3,C5(0,0,i1,i2,i3),
c     &    i1=1,2),i2=1,2),i3=1,2)

      do i1=1,2
      do i2=i1,2
      do i3=i2,2
      do i4=i3,2
        do j =1,2
          S5mod(j,i1,i2,i3,i4) = S5hat(j,i1,i2,i3,i4)
        end do
        S5mod(i1,i1,i2,i3,i4) = S5mod(i1,i1,i2,i3,i4) 
     &      - 2d0*C5(0,0,i2,i3,i4)
        S5mod(i2,i1,i2,i3,i4) = S5mod(i2,i1,i2,i3,i4) 
     &      - 2d0*C5(0,0,i1,i3,i4)
        S5mod(i3,i1,i2,i3,i4) = S5mod(i3,i1,i2,i3,i4) 
     &      - 2d0*C5(0,0,i1,i2,i4)
        S5mod(i4,i1,i2,i3,i4) = S5mod(i4,i1,i2,i3,i4) 
     &      - 2d0*C5(0,0,i1,i2,i3)
      end do
      end do
      end do
      end do

      do i1 = 1,2
      do i2 = i1,2
      do i3 = i2,2
      do i4 = i3,2
        C600mod(i1,i2,i3,i4) =  (C4(i1,i2,i3,i4) 
     &         - mxinv(0,1)*S5mod(1,i1,i2,i3,i4) 
     &         - mxinv(0,2)*S5mod(2,i1,i2,i3,i4) )
     &           /mxinv(0,0)
      end do
      end do
      end do
      end do

      do i1 = 1,2
      do i2 = i1,2
      do i3 = i2,2
      do i4 = i3,2
      do i5 = i4,2
        C5(i1,i2,i3,i4,i5) = mxinv(i1,0)*C600mod(i2,i3,i4,i5) 
     &         + mxinv(i1,1)*S5mod(1,i2,i3,i4,i5)
     &         + mxinv(i1,2)*S5mod(2,i2,i3,i4,i5)
c        write(testout,55) i1,i2,i3,i4,i5,C5(i1,i2,i3,i4,i5)
      end do
      end do
      end do
      end do
      end do
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


 100  continue

      end

***********************************************************************
      subroutine cCg12345(p10,p21,p20,m02,m12,m22,
     &                C0,C1,C2,C3,C4,C5,C6,rank)
***********************************************************************
*     3-point tensor coefficient functions C0,C1,C2,C3,C4             *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     expansion about Gram determinant = 0 to order ordg3             *
*     rank+ordg3 < 5                                                  *
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
*     i-th denominator cancelled: Bxwi                                *
*---------------------------------------------------------------------*
*     19.10.04  Ansgar Denner     last changed  23.02.05              *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p20
      complex*16 q10,q21,q20
      complex*16 m02,m12,m22
      complex*16 mm02,mm12,mm22
      real*8     k1k2
      complex*16 f(2),zadjf(2)
      complex*16 z(2,2),zadj(2,2),detz
      complex*16 B0w2,B1w2,B2w2(0:1,0:1),B3w2(0:1,0:1,0:1)
      complex*16 B4w2(0:1,0:1,0:1,0:1),B5w2(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w2(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0w1,B1w1,B2w1(0:1,0:1),B3w1(0:1,0:1,0:1)
      complex*16 B4w1(0:1,0:1,0:1,0:1),B5w1(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w1(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0w0,B1w0,B2w0(0:1,0:1),B3w0(0:1,0:1,0:1)
      complex*16 B4w0(0:1,0:1,0:1,0:1),B5w0(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w0(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0m(0:2),B1m(0:2,2),B2m(0:2,0:2,0:2)
     &          ,B3m(0:2,0:2,0:2,0:2),B4m(0:2,0:2,0:2,0:2,0:2)
     &          ,B5m(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 S1hat(2),S2hat(2,2),S3hat(2,0:2,0:2),S4hat(2,0:2,0:2,2)
     &          ,S5hat(2,0:2,0:2,0:2,0:2),S6hat(2,0:2,0:2,0:2,0:2,2)
      complex*16 S2mod(2,2),S3mod(2,0:2,0:2),S4mod(2,0:2,0:2,2)
     &          ,S5mod(2,0:2,0:2,0:2,0:2),S6mod(2,0:2,0:2,0:2,0:2,2)
      complex*16 S00,S001(2),S002(0:2,0:2),S003(0:2,0:2,0:2),
     &           S004(0:2,0:2,0:2,0:2)
      complex*16 C0,C1(2),C2(0:2,0:2),C3(0:2,0:2,0:2)
      complex*16 C4(0:2,0:2,0:2,0:2),C5(0:2,0:2,0:2,0:2,0:2)
      complex*16 C6(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 elimcminf2
c      complex*16 cC0f
      real*8     r1,r2
      integer    rank
      integer    i1,i2,i3,i4,i5,j,k,l,kt,lt,m,r,sgn
      integer    ltest,sym,flag
c      integer    testout
      save       flag

      integer    ordg3,ordgy3,ordgp3
      integer    ordg4,ordgy4,ordgp4
      common /ctord3/ ordg3,ordgy3,ordgp3
      common /ctord4/ ordg4,ordgy4,ordgp4
      common /ltest/  ltest
      common /sym/    sym
                                                                               
      data       B1m /6*0d0/, B2m /27*0d0/, B3m/81*0d0/, B4m/243*0d0/
     &         , B5m/729*0d0/
      data  flag /0/
c      data  testout /37/

c      write(testout,*) 'cCg12345 in',p10,p21,p20,m02,m12,m22,rank,ordg3
c      write(testout,*) 'cCg12345 in',C0,C1,C2
c      write(testout,*) 'cCg12345 in',C3
c      write(testout,*) 'cCg12345 in',C4


      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cCg12345'
        write(*,*) 'rank = ',rank
        stop
      else if (rank+ordg3.gt.5.and.flag.eq.0) then
        write(*,*) 'rank+ordg3 > 5 not implemented in cCg12345'
        write(*,*) 'rank  = ',rank
        write(*,*) 'ordg3 = ',ordg3
        flag = 1
        if (rank.gt.5) then
          write(*,*) 'rank > 5 not implemented in cCg12345'
          write(*,*) 'rank = ',rank
          stop
        end if
      end if

      call cBp12345(p21,m12,m22,B0w0,B1w0,B2w0,B3w0,B4w0,B5w0,B6w0
     &    ,rank+ordg3,0)
      call cBp12345(p20,m02,m22,B0w1,B1w1,B2w1,B3w1,B4w1,B5w1,B6w1
     &    ,rank+ordg3,0)
      call cBp12345(p10,m02,m12,B0w2,B1w2,B2w2,B3w2,B4w2,B5w2,B6w2
     &    ,rank+ordg3,0)

      mm02 = elimcminf2(m02)
      mm12 = elimcminf2(m12)
      mm22 = elimcminf2(m22)
      q10  = elimcminf2(p10)
      q21  = elimcminf2(p21)
      q20  = elimcminf2(p20)
 
      k1k2 = (q10+q20-q21)
      detz = 4d0*q10*q20-k1k2*k1k2

      if (abs(detz/( 4d0*q10*q20 + k1k2*k1k2)).lt.1d-4) then
        if (abs(q10-q20).lt.abs(q10-q21).and.
     &      abs(q10-q20).lt.abs(q20-q21)) then
          detz  =  4d0*q10*q21 - (q10-q20+q21)*(q10-q20+q21)
        end if
      end if

c      write(testout,*) 'detz = ',detz/(4d0*q10*q20+k1k2*k1k2)
      z(1,1) = 2d0*q10
      z(1,2) = k1k2
      z(2,1) = k1k2
      z(2,2) = 2d0*q20
      zadj(1,1) = 2d0*q20
      zadj(1,2) = -k1k2
      zadj(2,1) = -k1k2
      zadj(2,2) = 2d0*q10
      f(1) = q10+mm02-mm12
      f(2) = q20+mm02-mm22

      zadjf(1) = zadj(1,1)*f(1)+zadj(1,2)*f(2)
      zadjf(2) = zadj(2,1)*f(1)+zadj(2,2)*f(2)

c      write(testout,*) zadj(1,1), zadj(1,2), zadj(2,1), zadj(2,2)
c      write(testout,*) f(1),f(2)
c      write(testout,*) zadjf(1),zadjf(2)
      
c>      C1(1) = 0d0
c>      C1(2) = 0d0
c>
c>      C2(1,1) = 0d0
c>      C2(1,2) = 0d0
c>      C2(2,1) = 0d0
c>      C2(2,2) = 0d0
c>
c>      C3(1,1,1) = 0d0
c>      C3(1,1,2) = 0d0
c>      C3(1,2,1) = 0d0
c>      C3(2,1,1) = 0d0
c>      C3(1,2,2) = 0d0
c>      C3(2,1,2) = 0d0
c>      C3(2,2,1) = 0d0
c>      C3(2,2,2) = 0d0
c>
c>      C4(0,0,1,1) = 0d0
c>      C4(0,0,1,2) = 0d0
c>      C4(0,0,2,1) = 0d0
c>      C4(0,0,2,2) = 0d0
c>      C4(1,1,1,1) = 0d0
c>      C4(1,1,1,2) = 0d0
c>      C4(1,1,2,1) = 0d0
c>      C4(1,2,1,1) = 0d0
c>      C4(2,1,1,1) = 0d0
c>      C4(1,1,2,2) = 0d0
c>      C4(1,2,1,2) = 0d0
c>      C4(2,1,1,2) = 0d0
c>      C4(1,2,2,1) = 0d0
c>      C4(2,1,2,1) = 0d0
c>      C4(2,2,1,1) = 0d0
c>      C4(1,2,2,2) = 0d0
c>      C4(2,1,2,2) = 0d0
c>      C4(2,2,1,2) = 0d0
c>      C4(2,2,2,1) = 0d0
c>      C4(2,2,2,2) = 0d0

      do i1=1,2
        C1(i1) = 0d0
        C3(0,0,i1) = 0d0
        C5(0,0,0,0,i1) = 0d0
      do i2=1,2
        C2(i1,i2) = 0d0
        C4(0,0,i1,i2) = 0d0
        C6(0,0,0,0,i1,i2) = 0d0
      do i3=1,2
        C3(i1,i2,i3) = 0d0
        C5(0,0,i1,i2,i3) = 0d0
      do i4=1,2
        C4(i1,i2,i3,i4) = 0d0
        C6(0,0,i1,i2,i3,i4) = 0d0
      do i5=1,2
        C5(i1,i2,i3,i4,i5) = 0d0
c      do i6=1,2
c        C6(i1,i2,i3,i4,i5,i6) = 0d0
c      end do
      end do
      end do
      end do
      end do
      end do

      r1 = abs(z(1,1)/z(1,2))
      r2 = abs(z(2,2)/z(2,1))

      B0m(0) = B0w0
      B0m(1) = B0w1
      B0m(2) = B0w2

      S1hat(1) = B0m(1) - B0m(0)
      S1hat(2) = B0m(2) - B0m(0)

      if (rank+ordg3.eq.0) goto 99

      B1m(0,2) = B1w0
      B1m(0,1) = -B0m(0) - B1m(0,2)
      B1m(1,2) = B1w1
      B1m(2,1) = B1w2
        
      S2hat(1,1) =          - B1m(0,1)
      S2hat(1,2) = B1m(1,2) - B1m(0,2)
      S2hat(2,1) = B1m(2,1) - B1m(0,1)
      S2hat(2,2) =          - B1m(0,2)
                
      if (rank+ordg3.eq.1) goto 99

      B2m(0,0,0) = B2w0(0,0)
      B2m(1,0,0) = B2w1(0,0)
      B2m(2,0,0) = B2w2(0,0)
      B2m(0,2,2) = B2w0(1,1)
      B2m(0,1,2) = -B1m(0,2) - B2m(0,2,2)
      B2m(0,1,1) = -B1m(0,1) - B2m(0,1,2)
      B2m(1,2,2) = B2w1(1,1)
      B2m(2,1,1) = B2w2(1,1)

      B2m(0,2,1) =  B2m(0,1,2)
 
      S3hat(1,0,0) = B2m(1,0,0) - B2m(0,0,0)
      S3hat(2,0,0) = B2m(2,0,0) - B2m(0,0,0)
      S3hat(1,1,1) =            - B2m(0,1,1)
      S3hat(1,1,2) =            - B2m(0,1,2)
      S3hat(1,2,1) =            - B2m(0,1,2)
      S3hat(1,2,2) = B2m(1,2,2) - B2m(0,2,2)
      S3hat(2,1,1) = B2m(2,1,1) - B2m(0,1,1)
      S3hat(2,1,2) =            - B2m(0,1,2)
      S3hat(2,2,1) =            - B2m(0,1,2)
      S3hat(2,2,2) =            - B2m(0,2,2)

      if (rank+ordg3.eq.2) goto 99

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
      
      B3m(0,2,1,2) = B3m(0,1,2,2) 
cs      B3m(0,2,2,1) = B3m(0,1,2,2) 
cs      B3m(0,1,2,1) = B3m(0,1,1,2) 
      B3m(0,2,1,1) = B3m(0,1,1,2) 

      S4hat(1,0,0,1) =              - B3m(0,0,0,1)
      S4hat(2,0,0,1) = B3m(2,0,0,1) - B3m(0,0,0,1)
      S4hat(1,0,0,2) = B3m(1,0,0,2) - B3m(0,0,0,2)
      S4hat(2,0,0,2) =              - B3m(0,0,0,2)
        
      S4hat(1,1,1,1) =              - B3m(0,1,1,1)
cs      S4hat(1,1,2,1) =              - B3m(0,1,1,2)
      S4hat(1,2,1,1) =              - B3m(0,1,1,2)
cs      S4hat(1,2,2,1) =              - B3m(0,1,2,2)
      S4hat(2,1,1,1) = B3m(2,1,1,1) - B3m(0,1,1,1)
cs      S4hat(2,1,2,1) =              - B3m(0,1,1,2)
      S4hat(2,2,1,1) =              - B3m(0,1,1,2)
cs      S4hat(2,2,2,1) =              - B3m(0,1,2,2)
      S4hat(1,1,1,2) =              - B3m(0,1,1,2)
      S4hat(1,1,2,2) =              - B3m(0,1,2,2)
      S4hat(1,2,1,2) =              - B3m(0,1,2,2)
      S4hat(1,2,2,2) = B3m(1,2,2,2) - B3m(0,2,2,2)
      S4hat(2,1,1,2) =              - B3m(0,1,1,2)
      S4hat(2,1,2,2) =              - B3m(0,1,2,2)
      S4hat(2,2,1,2) =              - B3m(0,1,2,2)
      S4hat(2,2,2,2) =              - B3m(0,2,2,2)

      if (rank+ordg3.eq.3) goto 99

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
        

      S5hat(2,0,0,0,0) = B4m(2,0,0,0,0) - B4m(0,0,0,0,0)
      S5hat(1,0,0,0,0) = B4m(1,0,0,0,0) - B4m(0,0,0,0,0)
      do j=1,2
        do i1=1,2
        do i2=1,2
          S5hat(j,0,0,i1,i2) = B4m(j,0,0,i1,i2) - B4m(0,0,0,i1,i2)
          do i3=i2,2
          do i4=i3,2     
            S5hat(j,i1,i2,i3,i4) =             - B4m(0,i1,i2,i3,i4)
          end do
          end do
        end do
        end do
      end do
      S5hat(2,1,1,1,1) = B4m(2,1,1,1,1) - B4m(0,1,1,1,1)
      S5hat(1,2,2,2,2) = B4m(1,2,2,2,2) - B4m(0,2,2,2,2)

c      write(testout,*) 'S5hat ',S5hat(1,0,0,1,2),S5hat(2,0,0,1,2)
   
      if (rank+ordg3.eq.4) goto 99

      B5m(0,0,0,0,0,2) =  B5w0(0,0,0,0,1)
      B5m(0,0,0,0,0,1) = -B4m(0,0,0,0,0) - B5m(0,0,0,0,0,2)
      B5m(1,0,0,0,0,2) =  B5w1(0,0,0,0,1)
      B5m(2,0,0,0,0,1) =  B5w2(0,0,0,0,1)

      B5m(0,0,0,2,2,2) =  B5w0(0,0,1,1,1)
      B5m(0,0,0,1,2,2) = -B4m(0,0,0,2,2)-B5m(0,0,0,2,2,2)
      B5m(0,0,0,1,1,2) = -B4m(0,0,0,1,2)-B5m(0,0,0,1,2,2)
      B5m(0,0,0,1,1,1) = -B4m(0,0,0,1,1)-B5m(0,0,0,1,1,2)

      B5m(0,2,2,2,2,2) =  B5w0(1,1,1,1,1)
      B5m(0,1,2,2,2,2) = -B4m(0,2,2,2,2)-B5m(0,2,2,2,2,2)
      B5m(0,1,1,2,2,2) = -B4m(0,1,2,2,2)-B5m(0,1,2,2,2,2)
      B5m(0,1,1,1,2,2) = -B4m(0,1,1,2,2)-B5m(0,1,1,2,2,2)
      B5m(0,1,1,1,1,2) = -B4m(0,1,1,1,2)-B5m(0,1,1,1,2,2)
      B5m(0,1,1,1,1,1) = -B4m(0,1,1,1,1)-B5m(0,1,1,1,1,2)

      B5m(1,0,0,2,2,2) =  B5w1(0,0,1,1,1)
      B5m(2,0,0,1,1,1) =  B5w2(0,0,1,1,1)      
      
      B5m(1,2,2,2,2,2) =  B5w1(1,1,1,1,1)
      B5m(2,1,1,1,1,1) =  B5w2(1,1,1,1,1)
      
      B5m(0,0,0,2,1,1) =  B5m(0,0,0,1,1,2) 
cs      B5m(0,0,0,1,2,1) =  B5m(0,0,0,1,1,2) 
      B5m(0,0,0,2,1,2) =  B5m(0,0,0,1,2,2) 
cs      B5m(0,0,0,2,2,1) =  B5m(0,0,0,1,2,2) 

      B5m(0,2,1,1,1,1) =  B5m(0,1,1,1,1,2) 
cs      B5m(0,1,2,1,1,1) =  B5m(0,1,1,1,1,2) 
cs      B5m(0,1,1,2,1,1) =  B5m(0,1,1,1,1,2) 
cs      B5m(0,1,1,1,2,1) =  B5m(0,1,1,1,1,2) 
cs      B5m(0,1,1,2,1,2) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,1,1,2,2,1) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,1,2,1,1,2) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,1,2,1,2,1) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,1,2,2,1,1) =  B5m(0,1,1,1,2,2) 
      B5m(0,2,1,1,1,2) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,2,1,1,2,1) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,2,1,2,1,1) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,2,2,1,1,1) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,1,2,1,2,2) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,1,1,2,2) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,1,2,2,1,2) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,2,1,2,1,2) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,2,2,1,1,2) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,1,2,2,2,1) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,2,1,2,2,1) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,2,2,1,2,1) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,2,2,2,1,1) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,1,2,2,2) =  B5m(0,1,2,2,2,2) 
cs      B5m(0,2,2,1,2,2) =  B5m(0,1,2,2,2,2) 
cs      B5m(0,2,2,2,1,2) =  B5m(0,1,2,2,2,2) 
      B5m(0,2,2,2,2,1) =  B5m(0,1,2,2,2,2) 
        
      do i1=1,2
        do j=1,2
          S6hat(j,0,0,0,0,i1) = B5m(j,0,0,0,0,i1) - B5m(0,0,0,0,0,i1)
          do i2=1,2
          do i3=i2,2
            S6hat(j,0,0,i1,i2,i3) = 
     &          B5m(j,0,0,i1,i2,i3) - B5m(0,0,0,i1,i2,i3)
          do i4=i3,2
          do i5=i4,2     
            S6hat(j,i1,i2,i3,i4,i5) = -B5m(0,i1,i2,i3,i4,i5)
          end do
          end do
          end do
          end do
        end do
      end do
      S6hat(2,1,1,1,1,1) = B5m(2,1,1,1,1,1) - B5m(0,1,1,1,1,1)
      S6hat(1,2,2,2,2,2) = B5m(1,2,2,2,2,2) - B5m(0,2,2,2,2,2)


   99   continue

c choose reduction formulas with biggest denominators
      if (abs(zadjf(1)).ge.abs(zadjf(2))) then
        m = 1
      else 
        m = 2
      end if
      if (abs(z(1,1)).ge.abs(z(2,2)) .and.
     &    abs(z(1,1)).ge.abs(z(1,2)) ) then
        k = 1
        kt = 2
        l = 1
        lt = 2
        sgn = 1
c       sgn = (-1)**(k-l)
      else if (abs(z(2,2)).ge.abs(z(1,2))) then
        k = 2
        kt = 1
        l = 2
        lt = 1
        sgn = 1
      else 
        k = 1
        kt = 2
        l = 2
        lt = 1
        sgn = -1
      end if

c>      do 500 m = 1,2
c>      do 500 k = 1,2
c>      do 500 l = 1,2
c>        kt = 3-k
c>        lt = 3-l
c>        sgn = (-1)**(k-l)
c>        testout = testout+1
c>        write(testout,*) ' m = ',m
c>        write(testout,*) ' k = ',k
c>        write(testout,*) ' l = ',l

c      C0 = cC0f(p10,p21,p20,m02,m12,m22)

c      write(testout,*) 'ordg3 = ',ordg3
      do 100 r=0,ordg3+rank
c       do 100 r=0,0
c        write(testout,*) 'ordg3 = ',ordg3,r,z(k,l)

        C0 = (zadj(m,1)*S1hat(1) + zadj(m,2)*S1hat(2) 
     &      - detz*C1(m))/zadjf(m)
c        write(testout,9) 'C0 = ',C0
 9      format(1x,a10,2g24.16,2i2)  

        if (rank+ordg3-r.eq.0) goto 100
        
        S00 = 2d0*mm02*C0 + 2d0*B0m(0)
        
        C2(0,0) = (z(k,l)*(S00+1d0)
     &         - z(k,1)*S2hat(l,1) - z(k,2)*S2hat(l,2)
     &         + f(l)*S1hat(k)
     &         - f(k)*f(l)*C0
     &         - detz*C2(kt,lt)*sgn
     &           ) /(6d0*z(k,l))
c        write(testout,9) 'C00 = ',C2(0,0),k,l

        S2mod(1,1) =  S2hat(1,1) - 2d0*C2(0,0)
        S2mod(1,2) =  S2hat(1,2)
        S2mod(2,1) =  S2hat(2,1)
        S2mod(2,2) =  S2hat(2,2) - 2d0*C2(0,0)
      
        C1(1) = (zadj(m,1)*S2mod(1,1) + zadj(m,2)*S2mod(2,1) 
     &      - detz*C2(m,1))/zadjf(m)
        C1(2) = (zadj(m,1)*S2mod(1,2) + zadj(m,2)*S2mod(2,2) 
     &      - detz*C2(m,2))/zadjf(m)
c        write(testout,9) 'C1 = ',C1(1)
c        write(testout,9) 'C2 = ',C1(2)
        
        if (rank+ordg3-r.eq.1) goto 100
 
        S001(1) = 2d0*mm02*C1(1) + 2d0*B1m(0,1)
        S001(2) = 2d0*mm02*C1(2) + 2d0*B1m(0,2)

        C3(0,0,1) = z(k,l)*(S001(1)-1d0/3d0)
     &         - z(k,1)*S3hat(l,1,1) - z(k,2)*S3hat(l,2,1)
     &         + f(l)*S2hat(k,1)
     &         - f(k)*f(l)*C1(1)
     &         - detz*C3(kt,lt,1)*sgn
        if (k.eq.1) C3(0,0,1) = C3(0,0,1) - 2d0*f(l)*C2(0,0)
        if (l.eq.1) C3(0,0,1) = C3(0,0,1) - 2d0*f(k)*C2(0,0)
     &         + 2d0*S3hat(k,0,0)

        C3(0,0,1) = C3(0,0,1) /(10d0*z(k,l))
c        write(testout,9) 'C001 = ',C3(0,0,1)

        C3(0,0,2) = z(k,l)*(S001(2)-1d0/3d0)
     &         - z(k,1)*S3hat(l,1,2) - z(k,2)*S3hat(l,2,2)
     &         + f(l)*S2hat(k,2)
     &         - f(k)*f(l)*C1(2)
     &         - detz*C3(kt,lt,2)*sgn
        if (k.eq.2) C3(0,0,2) = C3(0,0,2) - 2d0*f(l)*C2(0,0)
        if (l.eq.2) C3(0,0,2) = C3(0,0,2) - 2d0*f(k)*C2(0,0)
     &         + 2d0*S3hat(k,0,0)
        C3(0,0,2) = C3(0,0,2) /(10d0*z(k,l))
c        write(testout,9) 'C002 = ',C3(0,0,2)
 
        S3mod(1,1,1) =  S3hat(1,1,1) - 4d0*C3(0,0,1)
        S3mod(1,1,2) =  S3hat(1,1,2) - 2d0*C3(0,0,2)
cs        S3mod(1,2,1) =  S3hat(1,2,1) - 2d0*C3(0,0,2)
        S3mod(1,2,2) =  S3hat(1,2,2)
        S3mod(2,1,1) =  S3hat(2,1,1)
        S3mod(2,1,2) =  S3hat(2,1,2) - 2d0*C3(0,0,1)
cs        S3mod(2,2,1) =  S3hat(2,2,1) - 2d0*C3(0,0,1)
        S3mod(2,2,2) =  S3hat(2,2,2) - 4d0*C3(0,0,2)

        C2(1,1) = (zadj(m,1)*S3mod(1,1,1) + zadj(m,2)*S3mod(2,1,1) 
     &      - detz*C3(m,1,1))/zadjf(m)
c        write(testout,9) 'C11 = ',C2(1,1)
        C2(1,2) = (zadj(m,1)*S3mod(1,1,2) + zadj(m,2)*S3mod(2,1,2) 
     &      - detz*C3(m,1,2))/zadjf(m)
c        write(testout,9) 'C12 = ',C2(1,2)
cs        C2(2,1) = (zadj(m,1)*S3mod(1,2,1) + zadj(m,2)*S3mod(2,2,1) 
cs     &      - detz*C3(m,2,1))/zadjf(m)
cs        write(testout,9) 'C21 = ',C2(2,1)
        C2(2,2) = (zadj(m,1)*S3mod(1,2,2) + zadj(m,2)*S3mod(2,2,2) 
     &      - detz*C3(m,2,2))/zadjf(m)
c        write(testout,9) 'C22 = ',C2(2,2)
        C2(2,1) = C2(1,2)

        if (rank+ordg3-r.eq.2) goto 100
        
c>        S002(0,0) = 2d0*mm02*C2(0,0) + 2d0*B2m(0,0,0)
c>        S002(1,1) = 2d0*mm02*C2(1,1) + 2d0*B2m(0,1,1)
c>        S002(1,2) = 2d0*mm02*C2(1,2) + 2d0*B2m(0,1,2)
c>        S002(2,1) = 2d0*mm02*C2(2,1) + 2d0*B2m(0,1,2)
c>        S002(2,2) = 2d0*mm02*C2(2,2) + 2d0*B2m(0,2,2)
c>
c>        C4(0,0,0,0) = (z(k,l)*(S002(0,0)+(mm02+mm12+mm22)/6d0
c>     &                        -(q10+q20+q21)/24d0)
c>     &         - z(k,1)*S4hat(l,0,0,1) - z(k,2)*S4hat(l,0,0,2)
c>     &         + f(l)*S3hat(k,0,0)
c>     &         - f(k)*f(l)*C2(0,0)
c>     &         - detz*C4(0,0,kt,lt)*sgn
c>     &           ) /(10d0*z(k,l))
c>c        write(testout,9) 'C0000 = ',C4(0,0,0,0) 
c>
c>        C4(0,0,1,1) = z(k,l)*(S002(1,1)+1d0/6d0)
c>     &         - z(k,1)*S4hat(l,1,1,1) - z(k,2)*S4hat(l,2,1,1)
c>     &         + f(l)*S3hat(k,1,1)
c>     &         - f(k)*f(l)*C2(1,1)
c>     &         - detz*C4(kt,lt,1,1)*sgn
c>        if (k.eq.1) then
c>          C4(0,0,1,1) = C4(0,0,1,1) - 4d0*f(l)*C3(0,0,1)
c>          if (l.eq.1) then
c>             C4(0,0,1,1) = C4(0,0,1,1) - 8d0*C4(0,0,0,0)
c>          end if
c>        end if
c>        if (l.eq.1) then
c>          C4(0,0,1,1) = C4(0,0,1,1) - 4d0*f(k)*C3(0,0,1)
c>     &         + 4d0*S4hat(k,0,0,1)
c>        end if
c>        C4(0,0,1,1) = C4(0,0,1,1) /(14d0*z(k,l))
c>c        write(testout,9) 'C0011 = ',C4(0,0,1,1)
c>
c>        C4(0,0,1,2) = z(k,l)*(S002(1,2)+1d0/12d0)
c>     &         - z(k,1)*S4hat(l,1,1,2) - z(k,2)*S4hat(l,2,1,2)
c>     &         + f(l)*S3hat(k,1,2)
c>     &         - f(k)*f(l)*C2(1,2)
c>     &         - detz*C4(kt,lt,1,2)*sgn
c>        if (k.eq.1) C4(0,0,1,2) = C4(0,0,1,2) - 2d0*f(l)*C3(0,0,2)
c>        if (k.eq.2) C4(0,0,1,2) = C4(0,0,1,2) - 2d0*f(l)*C3(0,0,1)
c>        if (l.eq.1) then 
c>          C4(0,0,1,2) = C4(0,0,1,2) - 2d0*f(k)*C3(0,0,2)
c>     &         + 2d0*S4hat(k,0,0,2)
c>          if (k.eq.2) C4(0,0,1,2) = C4(0,0,1,2) - 4d0*C4(0,0,0,0)
c>        end if
c>        if (l.eq.2) then
c>          C4(0,0,1,2) = C4(0,0,1,2) - 2d0*f(k)*C3(0,0,1)
c>     &         + 2d0*S4hat(k,0,0,1)
c>          if (k.eq.1) C4(0,0,1,2) = C4(0,0,1,2) - 4d0*C4(0,0,0,0)
c>        end if
c>        C4(0,0,1,2) = C4(0,0,1,2) /(14d0*z(k,l))
c>c        write(testout,9) 'C0012 = ',C4(0,0,1,2)
c>
c>        C4(0,0,2,1) = z(k,l)*(S002(2,1)+1d0/12d0)
c>     &         - z(k,1)*S4hat(l,1,2,1) - z(k,2)*S4hat(l,2,2,1)
c>     &         + f(l)*S3hat(k,2,1)
c>     &         - f(k)*f(l)*C2(2,1)
c>     &         - detz*C4(kt,lt,2,1)*sgn
c>        if (k.eq.1) C4(0,0,2,1) = C4(0,0,2,1) - 2d0*f(l)*C3(0,0,2)
c>        if (k.eq.2) C4(0,0,2,1) = C4(0,0,2,1) - 2d0*f(l)*C3(0,0,1)
c>        if (l.eq.1) then 
c>          C4(0,0,2,1) = C4(0,0,2,1) - 2d0*f(k)*C3(0,0,2)
c>     &         + 2d0*S4hat(k,0,0,2)
c>          if (k.eq.2) C4(0,0,2,1) = C4(0,0,2,1) - 4d0*C4(0,0,0,0)
c>        end if
c>        if (l.eq.2) then
c>          C4(0,0,2,1) = C4(0,0,2,1) - 2d0*f(k)*C3(0,0,1)
c>     &         + 2d0*S4hat(k,0,0,1)
c>          if (k.eq.1) C4(0,0,2,1) = C4(0,0,2,1) - 4d0*C4(0,0,0,0)
c>        end if
c>        C4(0,0,2,1) = C4(0,0,2,1) /(14d0*z(k,l))
c>c        write(testout,9) 'C0012 = ',C4(0,0,2,1)
c>
c>        C4(0,0,2,2) = z(k,l)*(S002(2,2)+1d0/6d0)
c>     &         - z(k,1)*S4hat(l,1,2,2) - z(k,2)*S4hat(l,2,2,2)
c>     &         + f(l)*S3hat(k,2,2)
c>     &         - f(k)*f(l)*C2(2,2)
c>     &         - detz*C4(kt,lt,2,2)*sgn
c>        if (k.eq.2) C4(0,0,2,2) = C4(0,0,2,2) - 4d0*f(l)*C3(0,0,2)
c>        if (l.eq.2) then
c>          C4(0,0,2,2) = C4(0,0,2,2) - 4d0*f(k)*C3(0,0,2)
c>     &         + 4d0*S4hat(k,0,0,2)
c>          if (k.eq.2)  C4(0,0,2,2) = C4(0,0,2,2) - 8d0*C4(0,0,0,0)
c>        end if
c>        C4(0,0,2,2) = C4(0,0,2,2) /(14d0*z(k,l))
c>c        write(testout,9) 'C0022 = ',C4(0,0,2,2)
c>
c>c        if (sym.eq.1) then
c>c          C4(0,0,1,2) = C4(0,0,2,1)
c>c        end if
c>
c>        S4mod(1,1,1,1) =  S4hat(1,1,1,1) - 6d0*C4(0,0,1,1)
c>        S4mod(1,1,2,1) =  S4hat(1,1,2,1) - 4d0*C4(0,0,1,2)
c>        S4mod(1,2,1,1) =  S4hat(1,2,1,1) - 4d0*C4(0,0,2,1)
c>        S4mod(1,2,2,1) =  S4hat(1,2,2,1) - 2d0*C4(0,0,2,2)
c>        S4mod(2,1,1,1) =  S4hat(2,1,1,1)
c>        S4mod(2,1,2,1) =  S4hat(2,1,2,1) - 2d0*C4(0,0,1,1)
c>        S4mod(2,2,1,1) =  S4hat(2,2,1,1) - 2d0*C4(0,0,1,1)
c>        S4mod(2,2,2,1) =  S4hat(2,2,2,1) - 4d0*C4(0,0,2,1)
c>        S4mod(1,1,1,2) =  S4hat(1,1,1,2) - 4d0*C4(0,0,1,2)
c>        S4mod(1,1,2,2) =  S4hat(1,1,2,2) - 2d0*C4(0,0,2,2)
c>        S4mod(1,2,1,2) =  S4hat(1,2,1,2) - 2d0*C4(0,0,2,2)
c>        S4mod(1,2,2,2) =  S4hat(1,2,2,2)
c>        S4mod(2,1,1,2) =  S4hat(2,1,1,2) - 2d0*C4(0,0,1,1)
c>        S4mod(2,1,2,2) =  S4hat(2,2,1,2) - 4d0*C4(0,0,1,2)
c>        S4mod(2,2,1,2) =  S4hat(2,2,1,2) - 4d0*C4(0,0,1,2)
c>        S4mod(2,2,2,2) =  S4hat(2,2,2,2) - 6d0*C4(0,0,2,2)
c>     
c>        C3(1,1,1) = (zadj(m,1)*S4mod(1,1,1,1) + zadj(m,2)*S4mod(2,1,1,1) 
c>     &    - detz*C4(m,1,1,1))/zadjf(1)
c>c        write(testout,9) 'C111 = ',C3(1,1,1)
c>        C3(1,1,2) = (zadj(m,1)*S4mod(1,1,1,2) + zadj(m,2)*S4mod(2,1,1,2) 
c>     &    - detz*C4(m,1,1,2))/zadjf(1)
c>c        write(testout,9) 'C112 = ',C3(1,1,2)
c>        C3(1,2,1) = (zadj(m,1)*S4mod(1,1,2,1) + zadj(m,2)*S4mod(2,1,2,1) 
c>     &    - detz*C4(m,1,2,1))/zadjf(1)
c>c        write(testout,9) 'C121 = ',C3(1,2,1)
c>        C3(1,2,2) = (zadj(m,1)*S4mod(1,1,2,2) + zadj(m,2)*S4mod(2,1,2,2) 
c>     &      - detz*C4(m,1,2,2))/zadjf(1)
c>c        write(testout,9) 'C122 = ',C3(1,2,2)
c>        C3(2,1,1) = (zadj(m,1)*S4mod(1,2,1,1) + zadj(m,2)*S4mod(2,2,1,1) 
c>     &      - detz*C4(m,2,1,1))/zadjf(1)
c>c        write(testout,9) 'C211 = ',C3(2,1,1)
c>        C3(2,1,2) = (zadj(m,1)*S4mod(1,2,1,2) + zadj(m,2)*S4mod(2,2,1,2) 
c>     &      - detz*C4(m,2,1,1))/zadjf(1)
c>c        write(testout,9) 'C212 = ',C3(2,1,2)
c>        C3(2,2,1) = (zadj(m,1)*S4mod(1,2,2,1) + zadj(m,2)*S4mod(2,2,2,1) 
c>     &      - detz*C4(m,2,2,1))/zadjf(1)
c>c        write(testout,9) 'C221 = ',C3(2,2,1)
c>        C3(2,2,2) = (zadj(m,1)*S4mod(1,2,2,2) + zadj(m,2)*S4mod(2,2,2,2) 
c>     &      - detz*C4(m,2,2,2))/zadjf(1)
c>c        write(testout,9) 'C222 = ',C3(2,2,2)
c>
        
        S002(0,0) = 2d0*mm02*C2(0,0) + 2d0*B2m(0,0,0)
        
        do i1=1,2
        do i2=i1,2
          S002(i1,i2) = 2d0*mm02*C2(i1,i2) + 2d0*B2m(0,i1,i2)
        end do
        end do
   
        C4(0,0,0,0) = (z(k,l)*(S002(0,0)+(mm02+mm12+mm22)/6d0
     &                        -(q10+q20+q21)/24d0)
     &         - z(k,1)*S4hat(l,0,0,1) - z(k,2)*S4hat(l,0,0,2)
     &         + f(l)*S3hat(k,0,0)
     &         - f(k)*f(l)*C2(0,0)
     &         - detz*C4(0,0,kt,lt)*sgn
     &           ) /(10d0*z(k,l))
c        write(testout,9) 'C0000 = ',C4(0,0,0,0) 

        do i1=1,2
          do i2=i1,2     
            C4(0,0,i1,i2) = z(k,l)*(S002(i1,i2)+1d0/12d0)
     &          - z(k,1)*S4hat(l,1,i1,i2) - z(k,2)*S4hat(l,2,i1,i2)
     &          + f(l)*S3hat(k,i1,i2)
     &         - f(k)*f(l)*C2(i1,i2)
     &          - detz*C4(kt,lt,i1,i2)*sgn
            if (i1.eq.i2)  C4(0,0,i1,i2) = C4(0,0,i1,i2) + z(k,l)/12d0
            if (i1.eq.l) then
              C4(0,0,i1,i2) = C4(0,0,i1,i2) + 2d0*S4hat(k,0,0,i2)
     &            - 2d0*f(k)*C3(0,0,i2)
              if (i2.eq.k) then
                C4(0,0,i1,i2) = C4(0,0,i1,i2) - 4d0*C4(0,0,0,0)
              end if
            end if
            if (i1.eq.k) then
              C4(0,0,i1,i2) = C4(0,0,i1,i2) - 2d0*f(l)*C3(0,0,i2)
              if (i2.eq.l) then
                C4(0,0,i1,i2) = C4(0,0,i1,i2) - 4d0*C4(0,0,0,0)
              end if
            end if
            if (i2.eq.l) then
              C4(0,0,i1,i2) = C4(0,0,i1,i2) + 2d0*S4hat(k,0,0,i1)
     &            - 2d0*f(k)*C3(0,0,i1)
            end if
            if (i2.eq.k) then
              C4(0,0,i1,i2) = C4(0,0,i1,i2) - 2d0*f(l)*C3(0,0,i1)
            end if
            C4(0,0,i1,i2) = C4(0,0,i1,i2) / (14d0*z(k,l))
c            write(testout,20) 'C00ij = ',i1,i2,C4(0,0,i1,i2)             
          end do
        end do
        C4(0,0,2,1) = C4(0,0,1,2)

        do j=1,2
        do i1=1,2
        do i2=i1,2     
        do i3=i2,2     
          S4mod(j,i1,i2,i3) =  S4hat(j,i1,i2,i3)
          if (j.eq.i1) S4mod(j,i1,i2,i3) =  S4mod(j,i1,i2,i3)
     &        -2d0*C4(0,0,i2,i3)
          if (j.eq.i2) S4mod(j,i1,i2,i3) =  S4mod(j,i1,i2,i3)
     &        -2d0*C4(0,0,i1,i3)
          if (j.eq.i3) S4mod(j,i1,i2,i3) =  S4mod(j,i1,i2,i3)
     &        -2d0*C4(0,0,i1,i2)
        end do
        end do
        end do
        end do
 
        do i1=1,2
        do i2=i1,2     
        do i3=i2,2     
          C3(i1,i2,i3) = (zadj(m,1)*S4mod(1,i1,i2,i3) 
     &                  + zadj(m,2)*S4mod(2,i1,i2,i3) 
     &      - detz*C4(m,i1,i2,i3)
     &        )/zadjf(m)
c          write(testout,30) 'Cijk = ',i1,i2,i3,C3(i1,i2,i3)
        end do 
        end do 
        end do 
c        write(testout,*) 'C122   ',r,m,
c     &   zadj(m,1),S4mod(1,1,2,2), + zadj(m,2),S4mod(2,1,2,2), 
c     &      - detz,C4(m,1,2,2),zadjf(m)
c        write(testout,9) 'C122 = ',C3(1,2,2)
 

        C3(2,1,1) = C3(1,1,2)
        C3(1,2,1) = C3(1,1,2)
        C3(2,1,2) = C3(1,2,2)
        C3(2,2,1) = C3(1,2,2)
       
        if (rank+ordg3-r.eq.3) goto 100

        C5(0,0,0,0,1) = z(k,l)*(-5d0*mm02 - 10d0*mm12 - 5d0*mm22 
     &      + q20 + 2d0*q10 + 2d0*q21)
     &      /(120d0)
        C5(0,0,0,0,2) = z(k,l)*(-5d0*mm02 - 10d0*mm22 - 5d0*mm12 
     &      + q10 + 2d0*q20 + 2d0*q21)
     &      /(120d0)

        do i1=1,2
          S003(0,0,i1) = 2d0*mm02*C3(0,0,i1) + 2d0*B3m(0,0,0,i1)          
        do i2=i1,2
        do i3=i2,2
          S003(i1,i2,i3) = 2d0*mm02*C3(i1,i2,i3) + 2d0*B3m(0,i1,i2,i3)
        end do
        end do

        C5(0,0,0,0,i1) =  C5(0,0,0,0,i1) +
     &      z(k,l)*S003(0,0,i1)
     &         - z(k,1)*S5hat(l,0,0,1,i1) - z(k,2)*S5hat(l,0,0,2,i1)
     &         + f(l)*S4hat(k,0,0,i1)
     &         - f(k)*f(l)*C3(0,0,i1)
     &         - detz*C5(0,0,kt,lt,i1)*sgn

        if (k.eq.i1) C5(0,0,0,0,i1) = C5(0,0,0,0,i1) 
     &      - 2d0*f(l)*C4(0,0,0,0)
        if (l.eq.i1) C5(0,0,0,0,i1) = C5(0,0,0,0,i1) 
     &      - 2d0*f(k)*C4(0,0,0,0)
     &      + 2d0*S5hat(k,0,0,0,0)

        C5(0,0,0,0,i1) =  C5(0,0,0,0,i1) /(14d0*z(k,l))

c        write(testout,10) 'C0000i = ',i1,C5(0,0,0,0,i1),k,l 
 10       format(1x,a10,i2,2g24.16,2i2)  
        end do

        do i1=1,2
        do i2=i1,2     
        do i3=i2,2     
          C5(0,0,i1,i2,i3) = z(k,l)*(S003(i1,i2,i3)-1d0/30d0)
     &        - z(k,1)*S5hat(l,1,i1,i2,i3)-z(k,2)*S5hat(l,2,i1,i2,i3)
     &        + f(l)*S4hat(k,i1,i2,i3)
     &        - f(k)*f(l)*C3(i1,i2,i3)
     &        - detz*C5(kt,lt,i1,i2,i3)*sgn
c          write(testout,*) 'C00ijk  ',
c     &         z(k,l)*(S003(i1,i2,i3)-1d0/30d0),
c     &        - z(k,1)*S5hat(l,1,i1,i2,i3),-z(k,2)*S5hat(l,2,i1,i2,i3),
c     &        + f(l)*S4hat(k,i1,i2,i3),
c     &        - f(k)*f(l)*C3(i1,i2,i3),
c     &        - detz*C5(kt,lt,i1,i2,i3)*sgn
          if (i1.eq.i2.and.i1.eq.i3)  
     &        C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) - z(k,l)/15d0
          if (i1.eq.l) then
            C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3)
     &          + 2d0*S5hat(k,0,0,i2,i3) - 2d0*f(k)*C4(0,0,i2,i3)
            if (i2.eq.k) then
              C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) - 4d0*C5(0,0,0,0,i3)
            end if
            if (i3.eq.k) then
              C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) - 4d0*C5(0,0,0,0,i2)
            end if
          end if
          if (i1.eq.k) then
            C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) 
     &          - 2d0*f(l)*C4(0,0,i2,i3)
          end if
          if (i2.eq.l) then
            C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3)
     &          + 2d0*S5hat(k,0,0,i1,i3) - 2d0*f(k)*C4(0,0,i1,i3)
            if (i1.eq.k) then
              C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) - 4d0*C5(0,0,0,0,i3)
            end if
            if (i3.eq.k) then
              C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) - 4d0*C5(0,0,0,0,i1)
            end if
          end if
          if (i2.eq.k) then
            C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) 
     &          - 2d0*f(l)*C4(0,0,i1,i3)
          end if
          if (i3.eq.l) then
            C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3)
     &          + 2d0*S5hat(k,0,0,i1,i2) - 2d0*f(k)*C4(0,0,i1,i2)
            if (i1.eq.k) then
              C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) - 4d0*C5(0,0,0,0,i2)
            end if
            if (i2.eq.k) then
              C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) - 4d0*C5(0,0,0,0,i1)
            end if
          end if
          if (i3.eq.k) then
            C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) 
     &          - 2d0*f(l)*C4(0,0,i1,i2)
          end if
          C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) / (18d0*z(k,l))
c          write(testout,30) 'C00ijk = ',i1,i2,i3,C5(0,0,i1,i2,i3)             
 30       format(1x,a10,3i2,2g24.16)  
        end do
        end do
        end do
        C5(0,0,2,1,1) = C5(0,0,1,1,2)
        C5(0,0,1,2,1) = C5(0,0,1,1,2)
        C5(0,0,2,2,1) = C5(0,0,1,2,2)
        C5(0,0,2,1,2) = C5(0,0,1,2,2)

        do j=1,2
        do i1=1,2
        do i2=i1,2     
        do i3=i2,2     
        do i4=i3,2     
          S5mod(j,i1,i2,i3,i4) =  S5hat(j,i1,i2,i3,i4)
          if (j.eq.i1) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i2,i3,i4)
          if (j.eq.i2) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i1,i3,i4)
          if (j.eq.i3) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i1,i2,i4)
          if (j.eq.i4) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i1,i2,i3)
        end do
        end do
        end do
        end do
        end do
c        write(testout,*) 'C1122   ',r,m,
c     &   S5hat(1,1,1,2,2),C5(0,0,1,2,2),S5mod(1,1,1,2,2) 
 
        do i1=1,2
        do i2=i1,2     
        do i3=i2,2
        do i4=i3,2
          C4(i1,i2,i3,i4) = (zadj(m,1)*S5mod(1,i1,i2,i3,i4) 
     &                  + zadj(m,2)*S5mod(2,i1,i2,i3,i4) 
     &      - detz*C5(m,i1,i2,i3,i4)
     &        )/zadjf(m)
c          write(testout,40) 'Cijkl = ',i1,i2,i3,i4,C4(i1,i2,i3,i4)
        end do 
        end do 
        end do 
        end do 
c        write(testout,*) 'C1122   ',r,m,
c     &   zadj(m,1),S5mod(1,1,1,2,2), + zadj(m,2),S5mod(2,1,1,2,2), 
c     &      - detz,C5(m,1,1,2,2),zadjf(m)
c        write(testout,9) 'C1122 = ',C4(1,1,2,2)
 


      C4(1,2,1,1) = C4(1,1,1,2)
      C4(2,1,1,1) = C4(1,1,1,2)
      C4(2,1,1,2) = C4(1,1,2,2)
      C4(1,2,1,2) = C4(1,1,2,2)
      C4(2,2,1,1) = C4(1,1,2,2)
      C4(2,1,2,2) = C4(1,2,2,2)
      C4(2,2,1,2) = C4(1,2,2,2)

        if (rank+ordg3-r.eq.4) goto 100

        S004(0,0,0,0) = 2d0*mm02*C4(0,0,0,0) + 2d0*B4m(0,0,0,0,0)
        
        C6(0,0,0,0,0,0) = (z(k,l)*(S004(0,0,0,0)
     &      +(15d0*(mm02*(mm02+mm12)+mm22*(mm22+mm02)+mm12*(mm12+mm22))
     &       - 3d0*(mm02*q21+mm12*q20+mm22*q10)
     &       - 6d0*(mm02*(q10+q20)+mm12*(q21+q10)+mm22*(q20+q21))
     &      +  q10*(q10+q20)+q21*(q21+q10)+q20*(q20+q21))/720d0 )
     &         - z(k,1)*S6hat(l,0,0,0,0,1) - z(k,2)*S6hat(l,0,0,0,0,2)
     &         + f(l)*S5hat(k,0,0,0,0)
     &         - f(k)*f(l)*C4(0,0,0,0)
     &         - detz*C6(0,0,0,0,kt,lt)*sgn
     &           ) /(14d0*z(k,l))

c        write(testout,*) 'C000000 = ',C6(0,0,0,0,0,0) 

        C6(0,0,0,0,1,1) = z(k,l)*(6d0*(mm02+mm22) + 18d0*mm12
     &      - q20 - 3d0*(q10+q21) )
     &      /(360d0)
        C6(0,0,0,0,1,2) = z(k,l)*(3d0*mm02 + 6d0*(mm22+mm12) 
     &      - q10 - q20 - 2d0*q21)
     &      /(360d0)
        C6(0,0,0,0,2,1) = C6(0,0,0,0,1,2)
        C6(0,0,0,0,2,2) = z(k,l)*(6d0*(mm02+mm12) + 18d0*mm22 
     &      - q10 - 3d0*(q20+q21) )
     &      /(360d0)

        do i1=1,2
        do i2=i1,2
          S004(0,0,i1,i2) = 2d0*mm02*C4(0,0,i1,i2) 
     &        + 2d0*B4m(0,0,0,i1,i2)

        do i3=i2,2
        do i4=i3,2
          S004(i1,i2,i3,i4) = 2d0*mm02*C4(i1,i2,i3,i4) 
     &        + 2d0*B4m(0,i1,i2,i3,i4)
        end do
        end do

        C6(0,0,0,0,i1,i2) =  C6(0,0,0,0,i1,i2) +
     &      z(k,l)*S004(0,0,i1,i2)
     &         - z(k,1)*S6hat(l,0,0,1,i1,i2) 
     &         - z(k,2)*S6hat(l,0,0,2,i1,i2)
     &         + f(l)*S5hat(k,0,0,i1,i2)
     &         - f(k)*f(l)*C4(0,0,i1,i2)
     &         - detz*C6(0,0,kt,lt,i1,i2)*sgn

        if (k.eq.i1) C6(0,0,0,0,i1,i2) = C6(0,0,0,0,i1,i2) 
     &      - 2d0*f(l)*C5(0,0,0,0,i2)
        if (l.eq.i1) then
          C6(0,0,0,0,i1,i2) = C6(0,0,0,0,i1,i2) 
     &      - 2d0*f(k)*C5(0,0,0,0,i2)
     &      + 2d0*S6hat(k,0,0,0,0,i2)
          if (k.eq.i2) then
            C6(0,0,0,0,i1,i2) = C6(0,0,0,0,i1,i2) 
     &          - 4d0*C6(0,0,0,0,0,0)
          end if
        end if
        if (k.eq.i2) C6(0,0,0,0,i1,i2) = C6(0,0,0,0,i1,i2) 
     &      - 2d0*f(l)*C5(0,0,0,0,i1)
        if (l.eq.i2) then
          C6(0,0,0,0,i1,i2) = C6(0,0,0,0,i1,i2) 
     &      - 2d0*f(k)*C5(0,0,0,0,i1)
     &      + 2d0*S6hat(k,0,0,0,0,i1)
          if (k.eq.i1) then
            C6(0,0,0,0,i1,i2) = C6(0,0,0,0,i1,i2) 
     &          - 4d0*C6(0,0,0,0,0,0)
          end if
        end if

        C6(0,0,0,0,i1,i2) =  C6(0,0,0,0,i1,i2) /(18d0*z(k,l))

c        write(testout,20) 'C0000ij = ',i1,i2,C6(0,0,0,0,i1,i2),k,l 
 20     format(1x,a10,2i2,2g24.16,2i2)  

        end do
        end do

        do i1=1,2
        do i2=i1,2     
        do i3=i2,2     
        do i4=i3,2     
          if (i1+i2+i3+i4.eq.4.or.i1+i2+i3+i4.eq.8) then  
            C6(0,0,i1,i2,i3,i4) = z(k,l)/15d0
          else if (i1+i2+i3+i4.eq.5.or.i1+i2+i3+i4.eq.7) then  
            C6(0,0,i1,i2,i3,i4) = z(k,l)/60d0
          else if (i1+i2+i3+i4.eq.6) then  
            C6(0,0,i1,i2,i3,i4) = z(k,l)/90d0
          end if
          C6(0,0,i1,i2,i3,i4) =  C6(0,0,i1,i2,i3,i4) 
     &        + z(k,l)*S004(i1,i2,i3,i4)
     &        - z(k,1)*S6hat(l,1,i1,i2,i3,i4)
     &        - z(k,2)*S6hat(l,2,i1,i2,i3,i4)
     &        + f(l)*S5hat(k,i1,i2,i3,i4)
     &        - f(k)*f(l)*C4(i1,i2,i3,i4)
c     &        - detz*C6(kt,lt,i1,i2,i3,i4)*sgn
          if (i1.eq.l) then
            C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4)
     &          + 2d0*S6hat(k,0,0,i2,i3,i4) - 2d0*f(k)*C5(0,0,i2,i3,i4)
            if (i2.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i3,i4)
            end if
            if (i3.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i2,i4)
            end if
            if (i4.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i2,i3)
            end if
          end if
          if (i1.eq.k) then
            C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &          - 2d0*f(l)*C5(0,0,i2,i3,i4)
          end if
          if (i2.eq.l) then
            C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4)
     &          + 2d0*S6hat(k,0,0,i1,i3,i4) - 2d0*f(k)*C5(0,0,i1,i3,i4)
            if (i1.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &         - 4d0*C6(0,0,0,0,i3,i4)
            end if
            if (i3.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i1,i4)
            end if
            if (i4.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i1,i3)
            end if
          end if
          if (i2.eq.k) then
            C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &          - 2d0*f(l)*C5(0,0,i1,i3,i4)
          end if
          if (i3.eq.l) then
            C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4)
     &          + 2d0*S6hat(k,0,0,i1,i2,i4) - 2d0*f(k)*C5(0,0,i1,i2,i4)
            if (i1.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i2,i4)
            end if
            if (i2.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i1,i4)
            end if
            if (i4.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i1,i2)
            end if
          end if
          if (i3.eq.k) then
            C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &          - 2d0*f(l)*C5(0,0,i1,i2,i4)
          end if
          if (i4.eq.l) then
            C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4)
     &          + 2d0*S6hat(k,0,0,i1,i2,i3) - 2d0*f(k)*C5(0,0,i1,i2,i3)
            if (i1.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i2,i3)
            end if
            if (i2.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i1,i3)
            end if
            if (i3.eq.k) then
              C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &            - 4d0*C6(0,0,0,0,i1,i2)
            end if
          end if
          if (i4.eq.k) then
            C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) 
     &          - 2d0*f(l)*C5(0,0,i1,i2,i3)
          end if
          C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) / (22d0*z(k,l))
c          write(testout,40) 'C00ijkl = ',i1,i2,i3,i4,C6(0,0,i1,i2,i3,i4)
 40       format(1x,a10,4i2,2g24.16)  
      end do
      end do
      end do
      end do

      do j=1,2
      do i1=1,2
      do i2=i1,2     
      do i3=i2,2     
      do i4=i3,2     
      do i5=i4,2     
        S6mod(j,i1,i2,i3,i4,i5) =  S6hat(j,i1,i2,i3,i4,i5)
        if (j.eq.i1) S6mod(j,i1,i2,i3,i4,i5) =  S6mod(j,i1,i2,i3,i4,i5)
     &      -2d0*C6(0,0,i2,i3,i4,i5)
        if (j.eq.i2) S6mod(j,i1,i2,i3,i4,i5) =  S6mod(j,i1,i2,i3,i4,i5)
     &      -2d0*C6(0,0,i1,i3,i4,i5)
        if (j.eq.i3) S6mod(j,i1,i2,i3,i4,i5) =  S6mod(j,i1,i2,i3,i4,i5)
     &      -2d0*C6(0,0,i1,i2,i4,i5) 
        if (j.eq.i4) S6mod(j,i1,i2,i3,i4,i5) =  S6mod(j,i1,i2,i3,i4,i5)
     &      -2d0*C6(0,0,i1,i2,i3,i5)
        if (j.eq.i5) S6mod(j,i1,i2,i3,i4,i5) =  S6mod(j,i1,i2,i3,i4,i5)
     &      -2d0*C6(0,0,i1,i2,i3,i4)
      end do
      end do
      end do
      end do
      end do
      end do
 
      do i1=1,2 
      do i2=i1,2     
      do i3=i2,2
      do i4=i3,2
      do i5=i4,2
        C5(i1,i2,i3,i4,i5) = (zadj(m,1)*S6mod(1,i1,i2,i3,i4,i5) 
     &      + zadj(m,2)*S6mod(2,i1,i2,i3,i4,i5) 
c     &      - detz*C6(m,i1,i2,i3,i4,i5)
     &      )/zadjf(m)
c        write(testout,50) 'Cijklm = ',i1,i2,i3,i4,i5,C5(i1,i2,i3,i4,i5)
 50       format(1x,a10,5i2,2g24.16)  
      end do 
      end do 
      end do 
      end do 
      end do 
      C5(2,1,1,1,1) = C5(1,1,1,1,2)
      C5(1,2,1,1,1) = C5(1,1,1,1,2)
      C5(2,1,1,1,2) = C5(1,1,1,2,2)
      C5(1,2,1,1,2) = C5(1,1,1,2,2)
      C5(2,2,1,1,1) = C5(1,1,1,2,2)
      C5(2,1,1,2,2) = C5(1,1,2,2,2)
      C5(1,2,1,2,2) = C5(1,1,2,2,2)
      C5(2,2,1,1,2) = C5(1,1,2,2,2)
      C5(2,1,2,2,2) = C5(1,2,2,2,2)
      C5(2,2,1,2,2) = C5(1,2,2,2,2)

 100   continue  

      if(sym.eq.1) then
c        if (rank.gt.1) then
c          C2(2,1) = C2(1,2)
c        if (rank.gt.2) then
c          C3(1,2,1) = C3(1,1,2)
c          C3(2,1,1) = C3(1,1,2)
c          C3(2,1,2) = C3(1,2,2)
c          C3(2,2,1) = C3(1,2,2)
        if (rank.gt.3) then
c          C4(0,0,2,1) = C4(0,0,1,2)
          C4(1,1,2,1) = C4(1,1,1,2)
c          C4(1,2,1,1) = C4(1,1,1,2)
c          C4(2,1,1,1) = C4(1,1,1,2)
c          C4(1,2,1,2) = C4(1,1,2,2)
c          C4(2,1,1,2) = C4(1,1,2,2)
          C4(1,2,2,1) = C4(1,1,2,2)
          C4(2,1,2,1) = C4(1,1,2,2)
c          C4(2,2,1,1) = C4(1,1,2,2)
c          C4(2,1,2,2) = C4(1,2,2,2)
c          C4(2,2,1,2) = C4(1,2,2,2)
          C4(2,2,2,1) = C4(1,2,2,2)
        if (rank.gt.4) then
c          C5(0,0,1,2,1) = C5(0,0,1,1,2)
c          C5(0,0,2,1,1) = C5(0,0,1,1,2)
c          C5(0,0,2,1,2) = C5(0,0,1,2,2)
c          C5(0,0,2,2,1) = C5(0,0,1,2,2)
          C5(1,1,1,2,1) = C5(1,1,1,1,2)
          C5(1,1,2,1,1) = C5(1,1,1,1,2)
c          C5(1,2,1,1,1) = C5(1,1,1,1,2)
c          C5(2,1,1,1,1) = C5(1,1,1,1,2)
c          C5(2,1,1,1,2) = C5(2,1,1,1,2)
          C5(1,1,2,1,2) = C5(1,1,1,2,2)
c          C5(1,2,1,1,2) = C5(1,1,1,2,2)
          C5(1,1,2,2,1) = C5(1,1,1,2,2)
          C5(1,2,1,2,1) = C5(1,1,1,2,2)
          C5(2,1,1,2,1) = C5(1,1,1,2,2)
          C5(1,2,2,1,1) = C5(1,1,1,2,2)
          C5(2,1,2,1,1) = C5(1,1,1,2,2)
c          C5(2,2,1,1,1) = C5(1,1,1,2,2)
c          C5(1,2,1,2,2) = C5(1,1,2,2,2)
c          C5(2,1,1,2,2) = C5(1,1,2,2,2)
          C5(1,2,2,1,2) = C5(1,1,2,2,2)
          C5(2,1,2,1,2) = C5(1,1,2,2,2)
c          C5(2,2,1,1,2) = C5(1,1,2,2,2)
          C5(1,2,2,2,1) = C5(1,1,2,2,2)
          C5(2,1,2,2,1) = C5(1,1,2,2,2)
          C5(2,2,1,2,1) = C5(1,1,2,2,2)
          C5(2,2,2,1,1) = C5(1,1,2,2,2)
c          C5(2,1,2,2,2) = C5(1,2,2,2,2)
c          C5(2,2,1,2,2) = C5(1,2,2,2,2)
          C5(2,2,2,1,2) = C5(1,2,2,2,2)
          C5(2,2,2,2,1) = C5(1,2,2,2,2)
        if (rank.gt.5) then
c          C6(0,0,0,0,2,1) = C6(0,0,0,0,1,2)
          C6(0,0,1,1,2,1) = C6(0,0,1,1,1,2)
          C6(0,0,1,2,1,1) = C6(0,0,1,1,1,2)
          C6(0,0,2,1,1,1) = C6(0,0,1,1,1,2)
          C6(0,0,1,2,1,2) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,1,2) = C6(0,0,1,1,2,2)
          C6(0,0,1,2,2,1) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,2,1) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,2,2) = C6(0,0,1,2,2,2)
          C6(0,0,2,2,1,2) = C6(0,0,1,2,2,2)
          C6(0,0,2,2,2,1) = C6(0,0,1,2,2,2)
        end if
        end if
        end if
c        end if
c        end if
      end if      

c 500  continue

      end

***********************************************************************

      subroutine cCgy12345(p10,p21,p20,m02,m12,m22,
     &                C0,C1,C2,C3,C4,C5,C6,rank)
***********************************************************************
*     3-point tensor coefficient functions C0,C1,C2,C3,C4,C5          *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     expansion about Gram determinant = 0                            *
*               and Cayley determinant = 0 to order ordgy3            *
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
*     i-th denominator cancelled: Bxwi                                *
*---------------------------------------------------------------------*
*     04.11.04  Ansgar Denner     last changed  25.02.05              *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p20
      complex*16 q10,q21,q20
      complex*16 m02,m12,m22
      complex*16 mm02,mm12,mm22
c      real*8     k1k2,mm0,mm1,mm2
      real*8     k1k2
      complex*16 f(2),zadjf(2),xadj(2,2)
      complex*16 z(2,2),zadj(2,2),detz
      complex*16 B0w2,B1w2,B2w2(0:1,0:1),B3w2(0:1,0:1,0:1)
      complex*16 B0w1,B1w1,B2w1(0:1,0:1),B3w1(0:1,0:1,0:1)
      complex*16 B0w0,B1w0,B2w0(0:1,0:1),B3w0(0:1,0:1,0:1)
      complex*16 B4w2(0:1,0:1,0:1,0:1),B5w2(0:1,0:1,0:1,0:1,0:1)
      complex*16 B6w2(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B4w1(0:1,0:1,0:1,0:1),B5w1(0:1,0:1,0:1,0:1,0:1)
      complex*16 B6w1(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B4w0(0:1,0:1,0:1,0:1),B5w0(0:1,0:1,0:1,0:1,0:1)
      complex*16 B6w0(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0m(0:2),B1m(0:2,2),B2m(0:2,0:2,0:2)
     &          ,B3m(0:2,0:2,0:2,0:2) ,B4m(0:2,0:2,0:2,0:2,0:2)
     &          ,B5m(0:2,0:2,0:2,0:2,0:2,0:2)
     &          ,B6m(0:2,0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 S1hat(2),S2hat(2,2),S3hat(2,0:2,0:2)
     &          ,S4hat(2,0:2,0:2,2) ,S5hat(2,0:2,0:2,0:2,0:2)
     &          ,S6hat(2,0:2,0:2,0:2,0:2,0:2)
     &          ,S7hat(2,0:2,0:2,0:2,0:2,0:2,0:2)
c      complex*16 S1mod(2),S2mod(2,2),S3mod(2,0:2,0:2)
      complex*16 S2mod(2,2),S3mod(2,0:2,0:2)
     &          ,S4mod(2,0:2,0:2,2),S5mod(2,0:2,0:2,0:2,0:2)
     &          ,S6mod(2,0:2,0:2,0:2,0:2,0:2)
c      complex*16 S00,S001(2)
      complex*16 S002(0:2,0:2),S003(0:2,0:2,0:2)
      complex*16 R4(2,2,2,2),R5(2,2,2,2,2),R6(2,2,2,2,2,2)
      complex*16 R7(2,2,2,2,2,2,2)
      complex*16 C0,C1(2),C2(0:2,0:2),C3(0:2,0:2,0:2)
      complex*16 C4(0:2,0:2,0:2,0:2),C5(0:2,0:2,0:2,0:2,0:2)
      complex*16 C6(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C7(0:2,0:2,0:2,0:2,0:2,0:2,0:2)
c      complex*16 cC0f 
      complex*16 elimcminf2
c      real*8     maxmp
      integer    rank,rankord,rankB
      integer    i1,i2,i3,i4,i5,i6,i7,r,j
      integer    a,b,c,d,at,bt,ct,dt,sgnab,sgncd
      integer    ltest,sym,flag
c      integer    testout

      integer    ordg3,ordgy3,ordgp3
      integer    ordg4,ordgy4,ordgp4
      common /ctord3/ ordg3,ordgy3,ordgp3
      common /ctord4/ ordg4,ordgy4,ordgp4
      common /ltest/  ltest
      common /sym/    sym
                                                                               
      data       B1m /6*0d0/, B2m /27*0d0/, B3m/81*0d0/, B4m/243*0d0/
     &          ,B5m /729*0d0/,B6m /2187*0d0/

c      data testout /39/
      data  flag /0/
      save  flag


c      write(testout,*) 'cCgy12345 in',p10,p21,p20,m02,m12,m22,rank
c      write(*,*) 'cCgy12345 in',C0,C1,C2
c      write(*,*) 'cCgy12345 in',C3
c      write(*,*) 'cCgy12345 in',C4

      rankord = rank + 2d0*ordgy3

      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cCgy12345'
        write(*,*) 'rank = ',rank
        stop
      else if (rankord.gt.5.and.flag.eq.0) then
        write(*,*) 'rank+2*ordgy3 > 5 not implemented in cCgy12345'
        write(*,*) 'rank   = ',rank
        write(*,*) 'ordgy3 = ',ordgy3
        flag = 1
        if (rank.gt.5) then
          write(*,*) 'rank > 5 not implemented in cCgy12345'
          write(*,*) 'rank = ',rank
          stop
        end if
      end if

      rankB = min(rankord+1,5)
      call cBp12345(p21,m12,m22,B0w0,B1w0,B2w0,B3w0,B4w0,B5w0,B6w0
     &    ,rankB,-1)
      call cBp12345(p20,m02,m22,B0w1,B1w1,B2w1,B3w1,B4w1,B5w1,B6w1
     &    ,rankB,-1)
      call cBp12345(p10,m02,m12,B0w2,B1w2,B2w2,B3w2,B4w2,B5w2,B6w2
     &    ,rankB,-1)

      mm02 = elimcminf2(m02)
      mm12 = elimcminf2(m12)
      mm22 = elimcminf2(m22)
      q10  = elimcminf2(p10)
      q21  = elimcminf2(p21)
      q20  = elimcminf2(p20)
 
      k1k2 = (q10+q20-q21)
      detz = 4d0*q10*q20-k1k2*k1k2

      if (abs(detz/( 4d0*q10*q20 + k1k2*k1k2)).lt.1d-4) then
        if (abs(q10-q20).lt.abs(q10-q21).and.
     &      abs(q10-q20).lt.abs(q20-q21)) then
          detz  =  4d0*q10*q21 - (q10-q20+q21)*(q10-q20+q21)
        end if
      end if

c      write(testout,*)'rank,ordgy3 = ',rank,ordgy3

      z(1,1) = 2d0*q10
      z(1,2) = k1k2
      z(2,1) = k1k2
      z(2,2) = 2d0*q20
      zadj(1,1) = 2d0*q20
      zadj(1,2) = -k1k2
      zadj(2,1) = -k1k2
      zadj(2,2) = 2d0*q10
      f(1) = q10+mm02-mm12
      f(2) = q20+mm02-mm22
c      write(testout,*)'zadj ',zadj(1,1),zadj(1,2),zadj(2,2)
c      write(testout,*)'f ',f(1),f(2)
 


      zadjf(1) = zadj(1,1)*f(1)+zadj(1,2)*f(2)
      zadjf(2) = zadj(2,1)*f(1)+zadj(2,2)*f(2)
c      write(testout,*)'dety = ',
c     & zadj(1,1)*f(1)+zadj(1,2)*f(2),zadj(2,1)*f(1)+zadj(2,2)*f(2),
c     & zadj(1,1)*f(1),zadj(1,2)*f(2),zadj(2,1)*f(1),zadj(2,2)*f(2)
      xadj(1,1) = 2d0*mm02*z(1,1) - f(1)*f(1)
      xadj(1,2) = 2d0*mm02*z(1,2) - f(1)*f(2)
      xadj(2,1) = 2d0*mm02*z(2,1) - f(1)*f(2)
      xadj(2,2) = 2d0*mm02*z(2,2) - f(2)*f(2)
c      write(testout,*) 'test = ', xadj(1,2),xadj(2,1),
c     &     mm02,z(1,2), f(1),f(2)
c      write(testout,*) 'detz = ',detz
c      write(testout,*) 'dety = ',abs(zadjf(1)),abs(zadjf(2))
c      write(testout,*) 'xadj   = ',abs(xadj(1,1)),abs(xadj(1,2)),
C     &      abs(xadj(2,2))
c      write(testout,*) 'detz/xadj = ',
c     &     detz/xadj(2,2), detz/xadj(1,2),
c     &     detz/xadj(1,1), detz/xadj(1,2)
c      write(testout,*)'dety/p2 = ',
c     &     zadjf(1)/dmax1(abs(p10),abs(p20),abs(p21))**2
c     &    ,zadjf(2)/dmax1(abs(p10),abs(p20),abs(p21))**2
c      write(testout,*) 'dety/xadj = ',
c     &     zadjf(1)/xadj(2,2), zadjf(1)/xadj(1,2),
c     &     zadjf(2)/xadj(1,1), zadjf(2)/xadj(1,2)
c      write(testout,*) 'xadj/p2 = ',
c     &    xadj(2,2)/dmax1(abs(p10),abs(p20),abs(p21))**2, 
c     &    xadj(1,2)/dmax1(abs(p10),abs(p20),abs(p21))**2,
c     &    xadj(1,1)/dmax1(abs(p10),abs(p20),abs(p21))**2
      
      do i1=1,2
        C1(i1) = 0d0
      do i2=1,2
        C2(i1,i2) = 0d0
        C4(0,0,i1,i2) = 0d0
      do i3=1,2
        C3(i1,i2,i3) = 0d0
        C5(0,0,i1,i2,i3) = 0d0
      do i4=1,2
        C4(i1,i2,i3,i4) = 0d0
      do i5=1,2
        C5(i1,i2,i3,i4,i5) = 0d0
      do i6=1,2
        C6(i1,i2,i3,i4,i5,i6) = 0d0
      do i7=1,2
        C7(i1,i2,i3,i4,i5,i6,i7) = 0d0
      end do
      end do
      end do
      end do
      end do
      end do
      end do

      B0m(0) = B0w0
      B0m(1) = B0w1
      B0m(2) = B0w2

      S1hat(1) = B0m(1) - B0m(0)
      S1hat(2) = B0m(2) - B0m(0)

      B1m(0,2) = B1w0
      B1m(0,1) = -B0m(0) - B1m(0,2)
      B1m(1,2) = B1w1
      B1m(2,1) = B1w2
 
c      S00 = 2d0*mm02*C0 + 2d0*B0m(0)
      S2hat(1,1) =          - B1m(0,1)
      S2hat(1,2) = B1m(1,2) - B1m(0,2)
      S2hat(2,1) = B1m(2,1) - B1m(0,1)
      S2hat(2,2) =          - B1m(0,2)

      if (rankord.gt.0) then

      B2m(0,0,0) = B2w0(0,0)
      B2m(1,0,0) = B2w1(0,0)
      B2m(2,0,0) = B2w2(0,0)
      B2m(0,2,2) = B2w0(1,1)
      B2m(0,1,2) = -B1m(0,2) - B2m(0,2,2)
      B2m(0,2,1) = B2m(0,1,2)
      B2m(0,1,1) = -B1m(0,1) - B2m(0,1,2)
      B2m(1,2,2) = B2w1(1,1)
      B2m(2,1,1) = B2w2(1,1)
 
c      S001(1) = 2d0*mm02*C1(1) + 2d0*B1m(0,1)
c      S001(2) = 2d0*mm02*C1(2) + 2d0*B1m(0,2)
      S3hat(1,1,1) =            - B2m(0,1,1)
      S3hat(1,1,2) =            - B2m(0,1,2)
      S3hat(1,2,1) =            - B2m(0,1,2)
      S3hat(1,2,2) = B2m(1,2,2) - B2m(0,2,2)
      S3hat(2,1,1) = B2m(2,1,1) - B2m(0,1,1)
      S3hat(2,1,2) =            - B2m(0,1,2)
      S3hat(2,2,1) =            - B2m(0,1,2)
      S3hat(2,2,2) =            - B2m(0,2,2)

      if (rankord.gt.1) then

      B3m(0,0,0,1) = -B2w0(0,0)-B3w0(0,0,1)
      B3m(0,0,0,2) = B3w0(0,0,1)
      B3m(1,0,0,2) = B3w1(0,0,1)
      B3m(2,0,0,1) = B3w2(0,0,1)
      B3m(0,2,2,2) = B3w0(1,1,1)

      B3m(0,1,2,2) = -B2m(0,2,2) - B3m(0,2,2,2)
      B3m(0,1,1,2) = -B2m(0,1,2) - B3m(0,1,2,2)
      B3m(0,1,1,1) = -B2m(0,1,1) - B3m(0,1,1,2)
      B3m(0,2,1,2) = B3m(0,1,2,2)
      B3m(0,2,2,1) = B3m(0,1,2,2)
      B3m(0,2,1,1) = B3m(0,1,1,2)
      B3m(0,1,2,1) = B3m(0,1,1,2)

      B3m(1,2,2,2) = B3w1(1,1,1)
      B3m(2,1,1,1) = B3w2(1,1,1)

      S4hat(1,1,1,1) =              - B3m(0,1,1,1)
      S4hat(1,1,1,2) =              - B3m(0,1,1,2)
      S4hat(1,1,2,1) =              - B3m(0,1,1,2)
      S4hat(1,2,1,1) =              - B3m(0,1,1,2)
      S4hat(1,1,2,2) =              - B3m(0,1,2,2)
      S4hat(1,2,1,2) =              - B3m(0,1,2,2)
      S4hat(1,2,2,1) =              - B3m(0,1,2,2)
      S4hat(1,2,2,2) = B3m(1,2,2,2) - B3m(0,2,2,2)
      S4hat(2,1,1,1) = B3m(2,1,1,1) - B3m(0,1,1,1)
      S4hat(2,1,1,2) =              - B3m(0,1,1,2)
      S4hat(2,1,2,1) =              - B3m(0,1,1,2)
      S4hat(2,2,1,1) =              - B3m(0,1,1,2)
      S4hat(2,1,2,2) =              - B3m(0,1,2,2)
      S4hat(2,2,1,2) =              - B3m(0,1,2,2)
      S4hat(2,2,2,1) =              - B3m(0,1,2,2)
      S4hat(2,2,2,2) =              - B3m(0,2,2,2)

      if (rankord.gt.2) then

      S3hat(1,0,0) = B2m(1,0,0) - B2m(0,0,0)
      S3hat(2,0,0) = B2m(2,0,0) - B2m(0,0,0)

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
      B4m(0,1,2,1,1) = B4m(0,1,1,1,2) 
      B4m(0,1,1,2,1) = B4m(0,1,1,1,2) 
      B4m(0,1,2,1,2) = B4m(0,1,1,2,2) 
      B4m(0,1,2,2,1) = B4m(0,1,1,2,2) 
      B4m(0,2,1,1,2) = B4m(0,1,1,2,2) 
      B4m(0,2,1,2,1) = B4m(0,1,1,2,2) 
      B4m(0,2,2,1,1) = B4m(0,1,1,2,2) 
      B4m(0,2,1,2,2) = B4m(0,1,2,2,2) 
      B4m(0,2,2,1,2) = B4m(0,1,2,2,2) 
      B4m(0,2,2,2,1) = B4m(0,1,2,2,2) 
        

c      S5hat(2,0,0,0,0) = B4m(2,0,0,0,0) - B4m(0,0,0,0,0)
c      S5hat(1,0,0,0,0) = B4m(1,0,0,0,0) - B4m(0,0,0,0,0)
      do j=1,2
        do i1=1,2
        do i2=1,2
c          S5hat(j,0,0,i1,i2) =  -B4m(0,0,0,i1,i2)
          do i3=1,2
          do i4=1,2     
            S5hat(j,i1,i2,i3,i4) = -B4m(0,i1,i2,i3,i4)
          end do
          end do
        end do
        end do
      end do
      S5hat(2,1,1,1,1) = B4m(2,1,1,1,1) - B4m(0,1,1,1,1)
      S5hat(1,2,2,2,2) = B4m(1,2,2,2,2) - B4m(0,2,2,2,2)
   
      if (rankord.gt.3) then

      S4hat(1,0,0,1) =              - B3m(0,0,0,1)
      S4hat(2,0,0,1) = B3m(2,0,0,1) - B3m(0,0,0,1)
      S4hat(1,0,0,2) = B3m(1,0,0,2) - B3m(0,0,0,2)
      S4hat(2,0,0,2) =              - B3m(0,0,0,2)

      B5m(0,0,0,0,0,2) =  B5w0(0,0,0,0,1)
      B5m(0,0,0,0,0,1) = -B4m(0,0,0,0,0) - B5m(0,0,0,0,0,2)
      B5m(1,0,0,0,0,2) =  B5w1(0,0,0,0,1)
      B5m(2,0,0,0,0,1) =  B5w2(0,0,0,0,1)

      B5m(0,0,0,2,2,2) =  B5w0(0,0,1,1,1)
      B5m(0,0,0,1,2,2) = -B4m(0,0,0,2,2)-B5m(0,0,0,2,2,2)
      B5m(0,0,0,1,1,2) = -B4m(0,0,0,1,2)-B5m(0,0,0,1,2,2)
      B5m(0,0,0,1,1,1) = -B4m(0,0,0,1,1)-B5m(0,0,0,1,1,2)

      B5m(0,2,2,2,2,2) =  B5w0(1,1,1,1,1)
      B5m(0,1,2,2,2,2) = -B4m(0,2,2,2,2)-B5m(0,2,2,2,2,2)
      B5m(0,1,1,2,2,2) = -B4m(0,1,2,2,2)-B5m(0,1,2,2,2,2)
      B5m(0,1,1,1,2,2) = -B4m(0,1,1,2,2)-B5m(0,1,1,2,2,2)
      B5m(0,1,1,1,1,2) = -B4m(0,1,1,1,2)-B5m(0,1,1,1,2,2)
      B5m(0,1,1,1,1,1) = -B4m(0,1,1,1,1)-B5m(0,1,1,1,1,2)

      B5m(1,0,0,2,2,2) =  B5w1(0,0,1,1,1)
      B5m(2,0,0,1,1,1) =  B5w2(0,0,1,1,1)      
      
      B5m(1,2,2,2,2,2) =  B5w1(1,1,1,1,1)
      B5m(2,1,1,1,1,1) =  B5w2(1,1,1,1,1)
      
      B5m(0,0,0,2,1,1) =  B5m(0,0,0,1,1,2) 
      B5m(0,0,0,1,2,1) =  B5m(0,0,0,1,1,2) 
      B5m(0,0,0,2,1,2) =  B5m(0,0,0,1,2,2) 
      B5m(0,0,0,2,2,1) =  B5m(0,0,0,1,2,2) 

      B5m(0,2,1,1,1,1) =  B5m(0,1,1,1,1,2) 
      B5m(0,1,2,1,1,1) =  B5m(0,1,1,1,1,2) 
      B5m(0,1,1,2,1,1) =  B5m(0,1,1,1,1,2) 
      B5m(0,1,1,1,2,1) =  B5m(0,1,1,1,1,2) 
      B5m(0,1,1,2,1,2) =  B5m(0,1,1,1,2,2) 
      B5m(0,1,1,2,2,1) =  B5m(0,1,1,1,2,2) 
      B5m(0,1,2,1,1,2) =  B5m(0,1,1,1,2,2) 
      B5m(0,1,2,1,2,1) =  B5m(0,1,1,1,2,2) 
      B5m(0,1,2,2,1,1) =  B5m(0,1,1,1,2,2) 
      B5m(0,2,1,1,1,2) =  B5m(0,1,1,1,2,2) 
      B5m(0,2,1,1,2,1) =  B5m(0,1,1,1,2,2) 
      B5m(0,2,1,2,1,1) =  B5m(0,1,1,1,2,2) 
      B5m(0,2,2,1,1,1) =  B5m(0,1,1,1,2,2) 
      B5m(0,1,2,1,2,2) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,1,1,2,2) =  B5m(0,1,1,2,2,2) 
      B5m(0,1,2,2,1,2) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,1,2,1,2) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,2,1,1,2) =  B5m(0,1,1,2,2,2) 
      B5m(0,1,2,2,2,1) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,1,2,2,1) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,2,1,2,1) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,2,2,1,1) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,1,2,2,2) =  B5m(0,1,2,2,2,2) 
      B5m(0,2,2,1,2,2) =  B5m(0,1,2,2,2,2) 
      B5m(0,2,2,2,1,2) =  B5m(0,1,2,2,2,2) 
      B5m(0,2,2,2,2,1) =  B5m(0,1,2,2,2,2) 
        
      do i1=1,2
        do j=1,2
c          S6hat(j,0,0,0,0,i1) = B5m(j,0,0,0,0,i1) - B5m(0,0,0,0,0,i1)
          do i2=1,2
          do i3=1,2
c            S6hat(j,0,0,i1,i2,i3) = 
c     &          B5m(j,0,0,i1,i2,i3) - B5m(0,0,0,i1,i2,i3)
          do i4=1,2
          do i5=1,2     
            S6hat(j,i1,i2,i3,i4,i5) = -B5m(0,i1,i2,i3,i4,i5)
          end do
          end do
          end do
          end do
        end do
      end do
      S6hat(2,1,1,1,1,1) = B5m(2,1,1,1,1,1) - B5m(0,1,1,1,1,1)
      S6hat(1,2,2,2,2,2) = B5m(1,2,2,2,2,2) - B5m(0,2,2,2,2,2)

      if (rankord.gt.4) then

      S5hat(2,0,0,0,0) = B4m(2,0,0,0,0) - B4m(0,0,0,0,0)
      S5hat(1,0,0,0,0) = B4m(1,0,0,0,0) - B4m(0,0,0,0,0)
      do j=1,2
        do i1=1,2
        do i2=1,2
          S5hat(j,0,0,i1,i2) =  -B4m(0,0,0,i1,i2)
        end do
        end do
      end do
      S5hat(2,0,0,1,1) = B4m(2,0,0,1,1) - B4m(0,0,0,1,1)
      S5hat(1,0,0,2,2) = B4m(1,0,0,2,2) - B4m(0,0,0,2,2)

      B6m(0,0,0,0,0,2,2) =  B6w0(0,0,0,0,1,1)
      B6m(0,0,0,0,0,1,2) = -B5m(0,0,0,0,0,2) - B6m(0,0,0,0,0,2,2)
      B6m(0,0,0,0,0,1,1) = -B5m(0,0,0,0,0,1) - B6m(0,0,0,0,0,1,2)
      B6m(1,0,0,0,0,2,2) =  B6w1(0,0,0,0,1,1)
      B6m(2,0,0,0,0,1,1) =  B6w2(0,0,0,0,1,1)

      B6m(0,0,0,2,2,2,2) =  B6w0(0,0,1,1,1,1)
      B6m(0,0,0,1,2,2,2) = -B5m(0,0,0,2,2,2)-B6m(0,0,0,2,2,2,2)
      B6m(0,0,0,1,1,2,2) = -B5m(0,0,0,1,2,2)-B6m(0,0,0,1,2,2,2)
      B6m(0,0,0,1,1,1,2) = -B5m(0,0,0,1,1,2)-B6m(0,0,0,1,1,2,2)
      B6m(0,0,0,1,1,1,1) = -B5m(0,0,0,1,1,1)-B6m(0,0,0,1,1,1,2)

      B6m(0,2,2,2,2,2,2) =  B6w0(1,1,1,1,1,1)
      B6m(0,1,2,2,2,2,2) = -B5m(0,2,2,2,2,2)-B6m(0,2,2,2,2,2,2)
      B6m(0,1,1,2,2,2,2) = -B5m(0,1,2,2,2,2)-B6m(0,1,2,2,2,2,2)
      B6m(0,1,1,1,2,2,2) = -B5m(0,1,1,2,2,2)-B6m(0,1,1,2,2,2,2)
      B6m(0,1,1,1,1,2,2) = -B5m(0,1,1,1,2,2)-B6m(0,1,1,1,2,2,2)
      B6m(0,1,1,1,1,1,2) = -B5m(0,1,1,1,1,2)-B6m(0,1,1,1,1,2,2)
      B6m(0,1,1,1,1,1,1) = -B5m(0,1,1,1,1,1)-B6m(0,1,1,1,1,1,2)

      B6m(1,0,0,2,2,2,2) =  B6w1(0,0,1,1,1,1)
      B6m(2,0,0,1,1,1,1) =  B6w2(0,0,1,1,1,1)      
      
      B6m(1,2,2,2,2,2,2) =  B6w1(1,1,1,1,1,1)
      B6m(2,1,1,1,1,1,1) =  B6w2(1,1,1,1,1,1)
      
      B6m(0,0,0,2,1,1,1) =  B6m(0,0,0,1,1,1,2) 
      B6m(0,0,0,1,2,1,1) =  B6m(0,0,0,1,1,1,2) 
      B6m(0,0,0,1,1,2,1) =  B6m(0,0,0,1,1,1,2) 
      B6m(0,0,0,2,2,1,1) =  B6m(0,0,0,1,1,2,2) 
      B6m(0,0,0,2,1,2,1) =  B6m(0,0,0,1,1,2,2) 
      B6m(0,0,0,1,2,2,1) =  B6m(0,0,0,1,1,2,2) 
      B6m(0,0,0,2,1,1,2) =  B6m(0,0,0,1,1,2,2) 
      B6m(0,0,0,1,2,1,2) =  B6m(0,0,0,1,1,2,2) 
      B6m(0,0,0,2,1,2,2) =  B6m(0,0,0,1,2,2,2) 
      B6m(0,0,0,2,2,1,2) =  B6m(0,0,0,1,2,2,2) 
      B6m(0,0,0,2,2,1,1) =  B6m(0,0,0,1,2,2,2) 

      B6m(0,2,1,1,1,1,1) =  B6m(0,1,1,1,1,1,2) 
      B6m(0,1,2,1,1,1,1) =  B6m(0,1,1,1,1,1,2) 
      B6m(0,1,1,2,1,1,1) =  B6m(0,1,1,1,1,1,2) 
      B6m(0,1,1,1,2,1,1) =  B6m(0,1,1,1,1,1,2) 
      B6m(0,1,1,1,1,2,1) =  B6m(0,1,1,1,1,1,2) 

      B6m(0,2,1,1,1,1,2) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,1,2,1,1,1,2) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,1,1,2,1,1,2) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,1,1,1,2,1,2) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,2,1,1,1,2,1) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,1,2,1,1,2,1) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,1,1,2,1,2,1) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,1,1,1,2,2,1) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,2,1,1,2,1,1) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,1,2,1,2,1,1) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,1,1,2,2,1,1) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,2,1,2,1,1,1) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,1,2,2,1,1,1) =  B6m(0,1,1,1,1,2,2) 
      B6m(0,2,2,1,1,1,1) =  B6m(0,1,1,1,1,2,2) 

      B6m(0,2,1,1,1,2,2) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,1,2,1,1,2,2) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,1,1,2,1,2,2) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,2,1,1,2,1,2) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,1,2,1,2,1,2) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,1,1,2,2,1,2) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,2,1,2,1,1,2) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,1,2,2,1,1,2) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,2,2,1,1,1,2) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,2,1,1,2,2,1) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,1,2,1,2,2,1) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,1,1,2,2,2,1) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,2,1,2,1,2,1) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,1,2,2,1,2,1) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,2,2,1,1,2,1) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,2,1,2,2,1,1) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,1,2,2,2,1,1) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,2,2,1,2,1,1) =  B6m(0,1,1,1,2,2,2) 
      B6m(0,2,2,2,1,1,1) =  B6m(0,1,1,1,2,2,2) 

      B6m(0,2,1,1,2,2,2) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,1,2,1,2,2,2) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,2,1,2,1,2,2) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,1,2,2,1,2,2) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,2,2,1,1,2,2) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,2,1,2,2,1,2) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,1,2,2,2,1,2) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,2,2,1,2,1,2) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,2,2,2,1,1,2) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,1,2,2,2,2,1) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,2,1,2,2,2,1) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,2,2,1,2,2,1) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,2,2,2,1,2,1) =  B6m(0,1,1,2,2,2,2) 
      B6m(0,2,2,2,2,1,1) =  B6m(0,1,1,2,2,2,2) 
        
      B6m(0,2,1,2,2,2,2) =  B6m(0,1,2,2,2,2,2) 
      B6m(0,2,2,1,2,2,2) =  B6m(0,1,2,2,2,2,2) 
      B6m(0,2,2,2,1,2,2) =  B6m(0,1,2,2,2,2,2) 
      B6m(0,2,2,2,2,1,2) =  B6m(0,1,2,2,2,2,2) 
      B6m(0,2,2,2,2,2,1) =  B6m(0,1,2,2,2,2,2) 

      do i1=1,2
      do i2=1,2
        do j=1,2
c          S7hat(j,0,0,0,0,i1,i2) = - B6m(0,0,0,0,0,i1,i2)
          do i3=1,2
          do i4=1,2
c            S7hat(j,0,0,i1,i2,i3,i4) =  - B6m(0,0,0,i1,i2,i3,i4)
          do i5=1,2
          do i6=1,2     
            S7hat(j,i1,i2,i3,i4,i5,i6) = -B6m(0,i1,i2,i3,i4,i5,i6)
          end do
          end do
          end do
          end do
        end do
      end do
      end do
c      S7hat(2,0,0,0,0,1,1) = B6m(2,0,0,0,0,1,1) - B6m(0,0,0,0,0,1,1)
c      S7hat(1,0,0,0,0,2,2) = B6m(1,0,0,0,0,2,2) - B6m(0,0,0,0,0,2,2)
c      S7hat(2,0,0,1,1,1,1) = B6m(2,0,0,1,1,1,1) - B6m(0,0,0,1,1,1,1)
c      S7hat(1,0,0,2,2,2,2) = B6m(1,0,0,2,2,2,2) - B6m(0,0,0,2,2,2,2)
      S7hat(2,1,1,1,1,1,1) = B6m(2,1,1,1,1,1,1) - B6m(0,1,1,1,1,1,1)
      S7hat(1,2,2,2,2,2,2) = B6m(1,2,2,2,2,2,2) - B6m(0,2,2,2,2,2,2)

      end if
      end if
      end if
      end if
      end if
 
c choose reduction formulas with biggest denominators
      if (abs(xadj(1,1)).ge.abs(xadj(1,2)).and.
     &    abs(xadj(1,1)).ge.abs(xadj(2,2))) then
        a=1
        b=1
        at=2
        bt=2
        sgnab = 1
      else if (abs(xadj(2,2)).ge.abs(xadj(1,2))) then
        a=2
        b=2
        at=1
        bt=1
        sgnab = 1
      else 
        a=1
        b=2
        at=2
        bt=1
        sgnab = -1
      end if
      if (abs(zadj(1,1)).ge.abs(zadj(2,2)).and.
     &     abs(zadj(1,1)).ge.abs(zadj(1,2))) then
        c=1
        d=1
        ct=2
        dt=2
        sgncd = 1
      else  if (abs(zadj(2,2)).ge.abs(zadj(1,2))) then
        c=2
        d=2
        ct=1
        dt=1
        sgncd = 1
      else
        c=1
        d=2
        ct=2
        dt=1
        sgncd = -1
      end if
c      write(testout,*) 'abc = ',a,b,c
c      write(testout,*) 'xadj  = ',maxxadj,
c     &    abs(xadj(1,1)),abs(xadj(2,2)),abs(xadj(1,2))
c      write(testout,*) 'zadj= ', abs(zadj(1,1)),abs(zadj(2,2))


c>      do 500 a = 1,2
c>      do 500 b = 1,2
c>      do 500 c = 1,2
c>        at = 3-a
c>        bt = 3-b
c>        ct = 3-c
c>        dt = 3-d
c>        sgnab = (-1)**(a-b)
c>        sgncd = (-1)**(c-d)
c>        testout = testout+1
c>        write(testout,*) ' a = ',a
c>        write(testout,*) ' b = ',b
c>        write(testout,*) ' c = ',c
c>        write(testout,*) ' d = ',d


      do 100 r=0,rankord/2
c      write(testout,*) 'ordgy3 = ',ordgy3,rankord,r,rankord-2*r
c      maxmp = dmax1(abs(q10),abs(q20),abs(q21),
c     &    abs(m02),abs(m12),abs(m22))
      C2(0,0) = ( zadj(d,1)*S2hat(1,c) + zadj(d,2)*S2hat(2,c)
     &          -zadjf(d)*C1(c) - detz * C2(d,c) )/(2d0*zadj(d,c))
c      write(testout,20) C2(0,0)
 20   format(' C2gy(0,0)           = ',G20.14,' + i* ',G20.14)

      C0 = ( z(a,b) * ( 4d0* C2(0,0) - 1d0 - B0m(0) )
     &      - f(b) * S1hat(a) 
     &      - sgnab* zadjf(at) * C1(bt) )
     &     /xadj(a,b)
c      write(testout,10) C0
 10   format(' C0gy                = ',G20.14,' + i* ',G20.14)

      if (rankord-2*r.gt.0) then

      C3(0,0,c) = ( zadj(d,1)*S3hat(1,c,c) + zadj(d,2)*S3hat(2,c,c)
     &          -zadjf(d)*C2(c,c) - detz*C3(d,c,c) )/(4d0*zadj(d,c))
      C3(0,0,ct) = ( zadj(d,1)*S3hat(1,c,ct) 
     &              + zadj(d,2)*S3hat(2,c,ct)
     &          -zadjf(d)*C2(c,ct) - detz*C3(d,c,ct)
     &            - 2d0*zadj(d,ct)*C3(0,0,c)  )/(2d0*zadj(d,c))

c      write(testout,31) (i1,C3(0,0,i1),i1=1,2)
 31   format(' C3gy(0,0,',i1,')         = ',G20.14,' + i* ',G20.14)

c      write(testout,*) 'C3gy(0,0,1) = ',C3(0,0,1)
c      write(testout,*) 'C3gy(0,0,2) = ',C3(0,0,2)

      S2mod(1,1) =  S2hat(1,1) - 2d0*C2(0,0)
      S2mod(1,2) =  S2hat(1,2)
      S2mod(2,1) =  S2hat(2,1)
      S2mod(2,2) =  S2hat(2,2) - 2d0*C2(0,0)
      
      C1(1) = ( z(a,b) * ( 6d0* C3(0,0,1) + 1d0/3d0 - B1m(0,1) )
     &      - f(b) * S2mod(a,1)
     &      - sgnab*zadjf(at) * C2(bt,1) )
     &     /xadj(a,b)
      C1(2) = ( z(a,b) * ( 6d0* C3(0,0,2) + 1d0/3d0 - B1m(0,2) )
     &      - f(b) * S2mod(a,2)
     &      - sgnab*zadjf(at) * C2(bt,2) )
     &     /xadj(a,b)

c      write(testout,11) (i1,C1(i1),i1=1,2)
 11   format(' C1gy(',i1,')             = ',G20.14,' + i* ',G20.14)

      if (rankord-2*r.gt.1) then

      do i1=1,2
      do i2=1,2
      do i3=1,2
        R4(d,i1,i2,i3) = .5d0*( zadj(d,1)*S4hat(1,i1,i2,i3) 
     &                  + zadj(d,2)*S4hat(2,i1,i2,i3)
     &                  -zadjf(d)*C3(i1,i2,i3) - detz*C4(d,i1,i2,i3) )
      end do
      end do
      end do
c      write(testout,802) (((i1,i2,i3,R4(c,i1,i2,i3),i1=1,2),i2=1,2),
c     &    i3=1,2)
c 802  format(' R4(c,',i1,',',i1,',',i1,') = ',G20.14,' + i* ',G20.14)

      C4(0,0,c,c) = R4(d,c,c,c)/(3d0*zadj(d,c))
      C4(0,0,c,3-c) = (R4(d,c,c,3-c)-C4(0,0,c,c)*zadj(d,3-c))
     &               /(2d0*zadj(d,c))
      C4(0,0,3-c,3-c) = (R4(d,c,3-c,3-c)-2d0*C4(0,0,c,3-c)*zadj(d,3-c))
     &               /(zadj(d,c))
      C4(0,0,3-c,c) = C4(0,0,c,3-c) 

c      write(testout,42) ((i1,i2,C4(0,0,i1,i2),i1=1,2),i2=1,2)
 42   format(' C4gy(0,0,',i1,',',i1,')       = ',G20.14,' + i* ',G20.14)

      S3mod(1,1,1) =  S3hat(1,1,1) - 4d0*C3(0,0,1)
      S3mod(1,1,2) =  S3hat(1,1,2) - 2d0*C3(0,0,2)
      S3mod(1,2,1) =  S3hat(1,2,1) - 2d0*C3(0,0,2)
      S3mod(1,2,2) =  S3hat(1,2,2)
      S3mod(2,1,1) =  S3hat(2,1,1)
      S3mod(2,1,2) =  S3hat(2,1,2) - 2d0*C3(0,0,1)
      S3mod(2,2,1) =  S3hat(2,2,1) - 2d0*C3(0,0,1)
      S3mod(2,2,2) =  S3hat(2,2,2) - 4d0*C3(0,0,2)
      
      do i1=1,2
      do i2=1,2
        C2(i1,i2) = ( z(a,b) * ( 8d0* C4(0,0,i1,i2) - 
     &                           (2d0-(i1-i2)**2)/12d0 
     &                          - B2m(0,i1,i2) )
     &      - f(b) * S3mod(a,i1,i2)
     &      - sgnab*zadjf(at) * C3(bt,i1,i2) )
     &     /xadj(a,b)
      end do
      end do

c      write(testout,22) ((i1,i2,C2(i1,i2),i1=1,2),i2=1,2)
 22   format(' C2gy(',i1,',',i1,')           = ',G20.14,' + i* ',G20.14)
    
      S002(0,0) = 2d0*mm02*C2(0,0) + 2d0*B2m(0,0,0)
      S002(1,1) = 2d0*mm02*C2(1,1) + 2d0*B2m(0,1,1)
      S002(1,2) = 2d0*mm02*C2(1,2) + 2d0*B2m(0,1,2)
      S002(2,1) = 2d0*mm02*C2(2,1) + 2d0*B2m(0,1,2)
      S002(2,2) = 2d0*mm02*C2(2,2) + 2d0*B2m(0,2,2)

      C4(0,0,0,0) = (z(dt,ct)*(S002(0,0)+(mm02+mm12+mm22)/6d0
     &                        -(q10+q20+q21)/24d0)
     &         - z(dt,1)*S4hat(ct,0,0,1) - z(dt,2)*S4hat(ct,0,0,2)
     &         + f(ct)*S3hat(dt,0,0)
     &         - f(dt)*f(ct)*C2(0,0)
     &         - sgncd*detz*C4(0,0,d,c)
     &           ) /(10d0*z(dt,ct))

c      write(testout,40) C4(0,0,0,0) 
 40   format(' C4gy(0,0,0,0)       = ',G20.14,' + i* ',G20.14)

      if (rankord-2*r.gt.2) then

      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
        R5(d,i1,i2,i3,i4) = .5d0*( zadj(d,1)*S5hat(1,i1,i2,i3,i4) 
     &              + zadj(d,2)*S5hat(2,i1,i2,i3,i4)
     &              - zadjf(d)*C4(i1,i2,i3,i4) - detz*C5(d,i1,i2,i3,i4)
     &      )
      end do
      end do
      end do
      end do

      C5(0,0,c,c,c) = R5(d,c,c,c,c)/(4d0*zadj(d,c))
      C5(0,0,c,c,3-c) = (R5(d,c,c,c,3-c)-C5(0,0,c,c,c)
     &                *zadj(d,3-c))
     &               /(3d0*zadj(d,c))
      C5(0,0,c,3-c,3-c) = (R5(d,c,c,3-c,3-c)-2d0*C5(0,0,c,c,3-c)
     &                *zadj(d,3-c))
     &               /(2d0*zadj(d,c))
      C5(0,0,3-c,3-c,3-c) = (R5(d,c,3-c,3-c,3-c)-3d0*C5(0,0,c,3-c,3-c)
     &                *zadj(d,3-c))
     &               /(zadj(d,c))
      C5(0,0,3-c,c,3-c) = C5(0,0,c,3-c,3-c)
      C5(0,0,3-c,3-c,c) = C5(0,0,c,3-c,3-c)
      C5(0,0,c,3-c,c)   = C5(0,0,c,c,3-c)
      C5(0,0,3-c,c,c)   = C5(0,0,c,c,3-c)

c      write(testout,53) (((i1,i2,i3,C5(0,0,i1,i2,i3),i1=1,2),
c     &                                        i2=1,2),i3=1,2)
 53   format(' C5gy(0,0,',i1,',',i1,',',i1,')     = ',
     &        G20.14,' + i* ',G20.14)

      S4mod(1,1,1,1) =  S4hat(1,1,1,1) - 6d0*C4(0,0,1,1)
      S4mod(1,1,2,1) =  S4hat(1,1,2,1) - 4d0*C4(0,0,1,2)
      S4mod(1,2,1,1) =  S4hat(1,2,1,1) - 4d0*C4(0,0,2,1)
      S4mod(1,2,2,1) =  S4hat(1,2,2,1) - 2d0*C4(0,0,2,2)
      S4mod(2,1,1,1) =  S4hat(2,1,1,1)
      S4mod(2,1,2,1) =  S4hat(2,1,2,1) - 2d0*C4(0,0,1,1)
      S4mod(2,2,1,1) =  S4hat(2,2,1,1) - 2d0*C4(0,0,1,1)
      S4mod(2,2,2,1) =  S4hat(2,2,2,1) - 4d0*C4(0,0,2,1)
      S4mod(1,1,1,2) =  S4hat(1,1,1,2) - 4d0*C4(0,0,1,2)
      S4mod(1,1,2,2) =  S4hat(1,1,2,2) - 2d0*C4(0,0,2,2)
      S4mod(1,2,1,2) =  S4hat(1,2,1,2) - 2d0*C4(0,0,2,2)
      S4mod(1,2,2,2) =  S4hat(1,2,2,2)
      S4mod(2,1,1,2) =  S4hat(2,1,1,2) - 2d0*C4(0,0,1,1)
      S4mod(2,1,2,2) =  S4hat(2,2,1,2) - 4d0*C4(0,0,1,2)
      S4mod(2,2,1,2) =  S4hat(2,2,1,2) - 4d0*C4(0,0,1,2)
      S4mod(2,2,2,2) =  S4hat(2,2,2,2) - 6d0*C4(0,0,2,2)
     
      do i1=1,2
      do i2=1,2
      do i3=1,2
        C3(i1,i2,i3) = ( z(a,b) * ( 10d0* C5(0,0,i1,i2,i3) + 
     &      (3d0-(i1-i2)**2-(i1-i3)**2-(i3-i2)**2)/30d0 
     &                          - B3m(0,i1,i2,i3) )
     &      - f(b) * S4mod(a,i1,i2,i3)
     &      - sgnab*zadjf(at) * C4(bt,i1,i2,i3) )
     &     /xadj(a,b)
      end do
      end do
      end do
c      write(testout,33) (((i1,i2,i3,C3(i1,i2,i3),i1=1,2),i2=1,2),i3=1,2)
 33   format(' C3gy(',i1,',',i1,',',i1,')         = ',
     &                   G20.14,' + i* ',G20.14)
    
      C5(0,0,0,0,1) = z(dt,ct)*(-5d0*mm02 - 10d0*mm12 - 5d0*mm22 
     &    + q20 + 2d0*q10 + 2d0*q21)
     &    /(120d0)
      C5(0,0,0,0,2) = z(dt,ct)*(-5d0*mm02 - 10d0*mm22 - 5d0*mm12 
     &    + q10 + 2d0*q20 + 2d0*q21)
     &    /(120d0)
      
      do i1=1,2
        S003(0,0,i1) = 2d0*mm02*C3(0,0,i1) + 2d0*B3m(0,0,0,i1)          
        do i2=i1,2
          do i3=i2,2
            S003(i1,i2,i3) = 2d0*mm02*C3(i1,i2,i3) + 2d0*B3m(0,i1,i2,i3)
          end do
        end do
        
        C5(0,0,0,0,i1) =  C5(0,0,0,0,i1) +
     &      z(dt,ct)*S003(0,0,i1)
     &      - z(dt,1)*S5hat(ct,0,0,1,i1) - z(dt,2)*S5hat(ct,0,0,2,i1)
     &      + f(ct)*S4hat(dt,0,0,i1)
     &      - f(dt)*f(ct)*C3(0,0,i1)
     &      - detz*C5(0,0,d,c,i1)*sgncd
        
        if (dt.eq.i1) C5(0,0,0,0,i1) = C5(0,0,0,0,i1) 
     &      - 2d0*f(ct)*C4(0,0,0,0)
        if (ct.eq.i1) C5(0,0,0,0,i1) = C5(0,0,0,0,i1) 
     &      - 2d0*f(dt)*C4(0,0,0,0)
     &      + 2d0*S5hat(dt,0,0,0,0)
        
        C5(0,0,0,0,i1) =  C5(0,0,0,0,i1) /(14d0*z(dt,ct))
      end do
c      write(testout,51) (i1,C5(0,0,0,0,i1),i1=1,2)
 51   format(' C5gy(0,0,0,0,',i1,')     = ',
     &    G20.14,' + i* ',G20.14)
      
      if (rankord-2*r.gt.3) then
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
        R6(d,i1,i2,i3,i4,i5) = .5d0*( zadj(d,1)*S6hat(1,i1,i2,i3,i4,i5) 
     &      + zadj(d,2)*S6hat(2,i1,i2,i3,i4,i5)
     &      - zadjf(d)*C5(i1,i2,i3,i4,i5) 
c     &      - detz*C6(d,i1,i2,i3,i4,i5) 
     &      )
      end do
      end do
      end do
      end do
      end do

      C6(0,0,c,c,c,c) = R6(d,c,c,c,c,c)/(5d0*zadj(d,c))
      C6(0,0,c,c,c,3-c) = (R6(d,c,c,c,c,3-c)-C6(0,0,c,c,c,c)
     &                *zadj(d,3-c))
     &               /(4d0*zadj(d,c))
      C6(0,0,c,c,3-c,3-c) = (R6(d,c,c,c,3-c,3-c)-2d0*C6(0,0,c,c,c,3-c)
     &                *zadj(d,3-c))
     &               /(3d0*zadj(d,c))
      C6(0,0,c,3-c,3-c,3-c) = (R6(d,c,c,3-c,3-c,3-c)
     &                -3d0*C6(0,0,c,c,3-c,3-c)*zadj(d,3-c))
     &               /(2d0*zadj(d,c))
      C6(0,0,3-c,3-c,3-c,3-c) = (R6(d,c,3-c,3-c,3-c,3-c)
     &                -4d0*C6(0,0,c,3-c,3-c,3-c)*zadj(d,3-c))
     &               /(zadj(d,c))
      C6(0,0,c,c,3-c,c)   = C6(0,0,c,c,c,3-c)
      C6(0,0,c,3-c,c,c)   = C6(0,0,c,c,c,3-c)
      C6(0,0,3-c,c,c,c)   = C6(0,0,c,c,c,3-c)
      C6(0,0,c,3-c,c,3-c)   = C6(0,0,c,c,3-c,3-c)
      C6(0,0,3-c,c,c,3-c)   = C6(0,0,c,c,3-c,3-c)
      C6(0,0,c,3-c,3-c,c)   = C6(0,0,c,c,3-c,3-c)
      C6(0,0,3-c,c,3-c,c)   = C6(0,0,c,c,3-c,3-c)
      C6(0,0,3-c,3-c,c,c)   = C6(0,0,c,c,3-c,3-c)
      C6(0,0,3-c,c,3-c,3-c) = C6(0,0,c,3-c,3-c,3-c)
      C6(0,0,3-c,3-c,c,3-c) = C6(0,0,c,3-c,3-c,3-c)
      C6(0,0,3-c,3-c,3-c,c) = C6(0,0,c,3-c,3-c,3-c)

c      write(testout,64) ((((i1,i2,i3,i4,C6(0,0,i1,i2,i3,i4),i1=1,2),
c     &                                      i2=1,2),i3=1,2),i4=1,2)
 64   format(' C6gy(0,0,',i1,',',i1,',',i1,',',i1,')   = ',
     &        G20.14,' + i* ',G20.14)
     
        do j=1,2
        do i1=1,2
        do i2=1,2     
        do i3=1,2     
        do i4=1,2     
          S5mod(j,i1,i2,i3,i4) =  S5hat(j,i1,i2,i3,i4)
          if (j.eq.i1) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i2,i3,i4)
          if (j.eq.i2) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i1,i3,i4)
          if (j.eq.i3) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i1,i2,i4)
          if (j.eq.i4) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i1,i2,i3)
        end do
        end do
        end do
        end do
        end do

      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
        if(i1+i2+i3+i4.eq.6) then
           C4(i1,i2,i3,i4) =  -z(a,b)/90d0
        else if(i1+i2+i3+i4.eq.5.or.i1+i2+i3+i4.eq.7) then
           C4(i1,i2,i3,i4) =  -z(a,b)/60d0
        else if(i1+i2+i3+i4.eq.4.or.i1+i2+i3+i4.eq.8) then
           C4(i1,i2,i3,i4) =  -z(a,b)/15d0
        end if
        C4(i1,i2,i3,i4) =  C4(i1,i2,i3,i4)
     &      + ( z(a,b) * ( 12d0* C6(0,0,i1,i2,i3,i4) 
     &                          - B4m(0,i1,i2,i3,i4) )
     &      - f(b) * S5mod(a,i1,i2,i3,i4)
     &      - sgnab*zadjf(at) * C5(bt,i1,i2,i3,i4) )
        C4(i1,i2,i3,i4) =  C4(i1,i2,i3,i4)  /xadj(a,b)
      end do
      end do
      end do
      end do
c      write(testout,44) ((((i1,i2,i3,i4,C4(i1,i2,i3,i4),i1=1,2),
c     &    i2=1,2),i3=1,2),i4=1,2)
 44   format(' C4gy(',i1,',',i1,',',i1,',',i1,')       = ',
     &                   G20.14,' + i* ',G20.14)

      if (rankord-2*r.gt.4) then


      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
      do i6=1,2
        R7(d,i1,i2,i3,i4,i5,i6) = .5d0*( 
     &        zadj(d,1)*S7hat(1,i1,i2,i3,i4,i5,i6) 
     &      + zadj(d,2)*S7hat(2,i1,i2,i3,i4,i5,i6)
c     &      - zadjf(d)*C6(i1,i2,i3,i4,i5,i6) 
c     &      - detz*C7(d,i1,i2,i3,i4,i5,i6) 
     &      )
      end do
      end do
      end do
      end do
      end do
      end do

      C7(0,0,c,c,c,c,c) = R7(d,c,c,c,c,c,c)/(6d0*zadj(d,c))
      C7(0,0,c,c,c,c,3-c) = (R7(d,c,c,c,c,c,3-c)
     &                 - C7(0,0,c,c,c,c,c)*zadj(d,3-c))
     &               /(5d0*zadj(d,c))
      C7(0,0,c,c,c,3-c,3-c) = (R7(d,c,c,c,c,3-c,3-c)
     &                 - 2d0*C7(0,0,c,c,c,c,3-c)*zadj(d,3-c))
     &               /(4d0*zadj(d,c))
      C7(0,0,c,c,3-c,3-c,3-c) = (R7(d,c,c,c,3-c,3-c,3-c)
     &                - 3d0*C7(0,0,c,c,c,3-c,3-c)*zadj(d,3-c))
     &               /(3d0*zadj(d,c))
      C7(0,0,c,3-c,3-c,3-c,3-c) = (R7(d,c,c,3-c,3-c,3-c,3-c)
     &                - 4d0*C7(0,0,c,c,3-c,3-c,3-c)*zadj(d,3-c))
     &               /(2d0*zadj(d,c))
      C7(0,0,3-c,3-c,3-c,3-c,3-c) = (R7(d,c,3-c,3-c,3-c,3-c,3-c)
     &                - 5d0*C7(0,0,c,3-c,3-c,3-c,3-c)*zadj(d,3-c))
     &               /(zadj(d,c))
      C7(0,0,c,c,c,3-c,c)       = C7(0,0,c,c,c,c,3-c)
      C7(0,0,c,c,3-c,c,c)       = C7(0,0,c,c,c,c,3-c)
      C7(0,0,c,3-c,c,c,c)       = C7(0,0,c,c,c,c,3-c)
      C7(0,0,3-c,c,c,c,c)       = C7(0,0,c,c,c,c,3-c)
      C7(0,0,c,c,3-c,c,3-c)     = C7(0,0,c,c,c,3-c,3-c)
      C7(0,0,c,3-c,c,c,3-c)     = C7(0,0,c,c,c,3-c,3-c)
      C7(0,0,3-c,c,c,c,3-c)     = C7(0,0,c,c,c,3-c,3-c)
      C7(0,0,c,c,3-c,3-c,c)     = C7(0,0,c,c,c,3-c,3-c)
      C7(0,0,c,3-c,c,3-c,c)     = C7(0,0,c,c,c,3-c,3-c)
      C7(0,0,3-c,c,c,3-c,c)     = C7(0,0,c,c,c,3-c,3-c)
      C7(0,0,c,3-c,3-c,c,c)     = C7(0,0,c,c,c,3-c,3-c)
      C7(0,0,3-c,c,3-c,c,c)     = C7(0,0,c,c,c,3-c,3-c)
      C7(0,0,3-c,3-c,c,c,c)     = C7(0,0,c,c,c,3-c,3-c)
      C7(0,0,c,3-c,c,3-c,3-c)   = C7(0,0,c,c,3-c,3-c,3-c)
      C7(0,0,3-c,c,c,3-c,3-c)   = C7(0,0,c,c,3-c,3-c,3-c)
      C7(0,0,c,3-c,3-c,c,3-c)   = C7(0,0,c,c,3-c,3-c,3-c)
      C7(0,0,3-c,c,3-c,c,3-c)   = C7(0,0,c,c,3-c,3-c,3-c)
      C7(0,0,3-c,3-c,c,c,3-c)   = C7(0,0,c,c,3-c,3-c,3-c)
      C7(0,0,c,3-c,3-c,3-c,c)   = C7(0,0,c,c,3-c,3-c,3-c)
      C7(0,0,3-c,c,3-c,3-c,c)   = C7(0,0,c,c,3-c,3-c,3-c)
      C7(0,0,3-c,3-c,c,3-c,c)   = C7(0,0,c,c,3-c,3-c,3-c)
      C7(0,0,3-c,3-c,3-c,c,c)   = C7(0,0,c,c,3-c,3-c,3-c)
      C7(0,0,3-c,c,3-c,3-c,3-c) = C7(0,0,c,3-c,3-c,3-c,3-c)
      C7(0,0,3-c,3-c,c,3-c,3-c) = C7(0,0,c,3-c,3-c,3-c,3-c)
      C7(0,0,3-c,3-c,3-c,c,3-c) = C7(0,0,c,3-c,3-c,3-c,3-c)
      C7(0,0,3-c,3-c,3-c,3-c,c) = C7(0,0,c,3-c,3-c,3-c,3-c)

c      write(testout,75) (((((i1,i2,i3,i4,i5,C7(0,0,i1,i2,i3,i4,i5),
c     &                        i1=1,2),i2=1,2),i3=1,2),i4=1,2),i5=1,2)
 75   format(' C7gy(0,0,',i1,',',i1,',',i1,',',i1,',',i1,') = ',
     &        G20.14,' + i* ',G20.14)
     
        do j=1,2
        do i1=1,2
        do i2=1,2     
        do i3=1,2     
        do i4=1,2     
        do i5=1,2     
          S6mod(j,i1,i2,i3,i4,i5) =  S6hat(j,i1,i2,i3,i4,i5)
          if (j.eq.i1) S6mod(j,i1,i2,i3,i4,i5) = S6mod(j,i1,i2,i3,i4,i5)
     &        -2d0*C6(0,0,i2,i3,i4,i5)
          if (j.eq.i2) S6mod(j,i1,i2,i3,i4,i5) = S6mod(j,i1,i2,i3,i4,i5)
     &        -2d0*C6(0,0,i1,i3,i4,i5)
          if (j.eq.i3) S6mod(j,i1,i2,i3,i4,i5) = S6mod(j,i1,i2,i3,i4,i5)
     &        -2d0*C6(0,0,i1,i2,i4,i5)
          if (j.eq.i4) S6mod(j,i1,i2,i3,i4,i5) = S6mod(j,i1,i2,i3,i4,i5)
     &        -2d0*C6(0,0,i1,i2,i3,i5)
          if (j.eq.i5) S6mod(j,i1,i2,i3,i4,i5) = S6mod(j,i1,i2,i3,i4,i5)
     &        -2d0*C6(0,0,i1,i2,i3,i4)
        end do
        end do
        end do
        end do
        end do
        end do

        do i1=1,2
        do i2=1,2     
        do i3=1,2     
        do i4=1,2     
        do i5=1,2     
c>      do i1=1,2
c>      do i2=i1,2
c>      do i3=i2,2
c>      do i4=i3,2
c>      do i5=i4,2
        if(i1+i2+i3+i4+i5.eq.7.or.i1+i2+i3+i4+i5.eq.8) then
           C5(i1,i2,i3,i4,i5) =  z(a,b)/210d0
        else if(i1+i2+i3+i4+i5.eq.6.or.i1+i2+i3+i4+i5.eq.9) then
           C5(i1,i2,i3,i4,i5) =  z(a,b)/105d0
        else if(i1+i2+i3+i4+i5.eq.5.or.i1+i2+i3+i4+i5.eq.10) then
           C5(i1,i2,i3,i4,i5) =  z(a,b)/21d0
        end if
        C5(i1,i2,i3,i4,i5) =  C5(i1,i2,i3,i4,i5)
     &      + ( z(a,b) * ( 14d0* C7(0,0,i1,i2,i3,i4,i5) 
     &                          - B5m(0,i1,i2,i3,i4,i5) )
     &      - f(b) * S6mod(a,i1,i2,i3,i4,i5)
c     &      - sgnab*zadjf(at) * C6(bt,i1,i2,i3,i4,i5) 
     &        )
        C5(i1,i2,i3,i4,i5) =  C5(i1,i2,i3,i4,i5)  /xadj(a,b)
      end do
      end do
      end do
      end do
      end do
      C5(2,1,1,1,1) = C5(1,1,1,1,2)
      C5(2,1,1,1,2) = C5(1,1,1,2,2)
      C5(2,1,1,2,2) = C5(1,1,2,2,2)
      C5(2,1,2,2,2) = C5(1,2,2,2,2)

c      write(testout,55) (((((i1,i2,i3,i4,i5,C5(i1,i2,i3,i4,i5),i1=1,2),
c     &    i2=1,2),i3=1,2),i4=1,2),i5=1,2)
 55   format(' C5gy(',i1,',',i1,',',i1,',',i1,',',i1,')     = ',
     &                   G20.14,' + i* ',G20.14)

      end if
      end if 
      end if
      end if
      end if

 100  continue

 500  continue

      if(sym.eq.1) then
        if (rank.gt.1) then
          C2(2,1) = C2(1,2)
        if (rank.gt.2) then
          C3(1,2,1) = C3(1,1,2)
          C3(2,1,1) = C3(1,1,2)
          C3(2,1,2) = C3(1,2,2)
          C3(2,2,1) = C3(1,2,2)
        if (rank.gt.3) then
          C4(0,0,2,1) = C4(0,0,1,2)
          C4(1,1,2,1) = C4(1,1,1,2)
          C4(1,2,1,1) = C4(1,1,1,2)
          C4(2,1,1,1) = C4(1,1,1,2)
          C4(1,2,1,2) = C4(1,1,2,2)
          C4(2,1,1,2) = C4(1,1,2,2)
          C4(1,2,2,1) = C4(1,1,2,2)
          C4(2,1,2,1) = C4(1,1,2,2)
          C4(2,2,1,1) = C4(1,1,2,2)
          C4(2,1,2,2) = C4(1,2,2,2)
          C4(2,2,1,2) = C4(1,2,2,2)
          C4(2,2,2,1) = C4(1,2,2,2)
        if (rank.gt.4) then
          C5(0,0,1,2,1) = C5(0,0,1,1,2)
          C5(0,0,2,1,1) = C5(0,0,1,1,2)
          C5(0,0,2,1,2) = C5(0,0,1,2,2)
          C5(0,0,2,2,1) = C5(0,0,1,2,2)
          C5(1,1,1,2,1) = C5(1,1,1,1,2)
          C5(1,1,2,1,1) = C5(1,1,1,1,2)
          C5(1,2,1,1,1) = C5(1,1,1,1,2)
          C5(2,1,1,1,1) = C5(1,1,1,1,2)
          C5(1,1,2,1,2) = C5(1,1,1,2,2)
          C5(1,2,1,1,2) = C5(1,1,1,2,2)
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
          C5(2,1,2,2,2) = C5(1,2,2,2,2)
          C5(2,2,1,2,2) = C5(1,2,2,2,2)
          C5(2,2,2,1,2) = C5(1,2,2,2,2)
          C5(2,2,2,2,1) = C5(1,2,2,2,2)
        if (rank.gt.5) then
          C6(0,0,0,0,2,1) = C6(0,0,0,0,1,2)
          C6(0,0,1,1,2,1) = C6(0,0,1,1,1,2)
          C6(0,0,1,2,1,1) = C6(0,0,1,1,1,2)
          C6(0,0,2,1,1,1) = C6(0,0,1,1,1,2)
          C6(0,0,1,2,1,2) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,1,2) = C6(0,0,1,1,2,2)
          C6(0,0,1,2,2,1) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,2,1) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,2,2) = C6(0,0,1,2,2,2)
          C6(0,0,2,2,1,2) = C6(0,0,1,2,2,2)
          C6(0,0,2,2,2,1) = C6(0,0,1,2,2,2)
        end if
        end if
        end if
        end if
        end if
      end if      

      end

***********************************************************************
      subroutine cCgp12345(p10,p21,p20,m02,m12,m22,
     &                C0,C1,C2,C3,C4,C5,C6,rank)
***********************************************************************
*     3-point tensor coefficient functions C0,C1,C2,C3,C4             *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     expansion about (pi pj) = 0 to order ordgp3                     *
*     rank+ordgp3 < 5                                                 *
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
*     i-th denominator cancelled: Bxwi                                *
*---------------------------------------------------------------------*
*     22.03.05  Ansgar Denner     last changed  23.03.05              *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p20
      complex*16 q10,q21,q20
      complex*16 m02,m12,m22
      complex*16 mm02,mm12,mm22
      real*8     k1k2
      complex*16 f(2),zadjf(2)
      complex*16 z(2,2),zadj(2,2),detz
      complex*16 B0w2,B1w2,B2w2(0:1,0:1),B3w2(0:1,0:1,0:1)
      complex*16 B4w2(0:1,0:1,0:1,0:1),B5w2(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w2(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0w1,B1w1,B2w1(0:1,0:1),B3w1(0:1,0:1,0:1)
      complex*16 B4w1(0:1,0:1,0:1,0:1),B5w1(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w1(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0w0,B1w0,B2w0(0:1,0:1),B3w0(0:1,0:1,0:1)
      complex*16 B4w0(0:1,0:1,0:1,0:1),B5w0(0:1,0:1,0:1,0:1,0:1)
     &    ,B6w0(0:1,0:1,0:1,0:1,0:1,0:1)
      complex*16 B0m(0:2),B1m(0:2,2),B2m(0:2,0:2,0:2)
     &          ,B3m(0:2,0:2,0:2,0:2),B4m(0:2,0:2,0:2,0:2,0:2)
     &          ,B5m(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 S1hat(2),S2hat(2,2),S3hat(2,0:2,0:2),S4hat(2,0:2,0:2,2)
     &          ,S5hat(2,0:2,0:2,0:2,0:2),S6hat(2,0:2,0:2,0:2,0:2,2)
      complex*16 S2mod(2,2),S3mod(2,0:2,0:2),S4mod(2,0:2,0:2,2)
     &          ,S5mod(2,0:2,0:2,0:2,0:2),S6mod(2,0:2,0:2,0:2,0:2,2)
      complex*16 S00,S001(2),S002(0:2,0:2),S003(0:2,0:2,0:2),
     &           S004(0:2,0:2,0:2,0:2)
      complex*16 C0,C1(2),C2(0:2,0:2),C3(0:2,0:2,0:2)
      complex*16 C4(0:2,0:2,0:2,0:2),C5(0:2,0:2,0:2,0:2,0:2)
      complex*16 C6(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 elimcminf2
c      complex*16 cC0f
      real*8     r1,r2
      integer    rank
      integer    i1,i2,i3,i4,i5,j,m,r
      integer    ltest,sym,flag
c      integer    testout
      save       flag

      integer    ordg3,ordgy3,ordgp3
      integer    ordg4,ordgy4,ordgp4
      common /ctord3/ ordg3,ordgy3,ordgp3
      common /ctord4/ ordg4,ordgy4,ordgp4
      common /ltest/  ltest
      common /sym/    sym
                                                                               
      data       B1m /6*0d0/, B2m /27*0d0/, B3m/81*0d0/, B4m/243*0d0/
     &         , B5m/729*0d0/
      data  flag /0/
c      data  testout /33/

c      write(testout,*) 'cCgp12345 in',p10,p21,p20,m02,m12,m22,rank,ordgp3
c      write(testout,*) 'cCgp12345 in',C0,C1,C2
c      write(testout,*) 'cCgp12345 in',C3
c      write(testout,*) 'cCgp12345 in',C4

      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cCgp12345'
        write(*,*) 'rank = ',rank
        stop
      else if (rank+ordgp3.gt.5.and.flag.eq.0) then
        write(*,*) 'rank+ordgp3 > 5 not implemented in cCgp12345'
        write(*,*) 'rank  = ',rank
        write(*,*) 'ordgp3= ',ordgp3
        flag = 1
        if (rank.gt.5) then
          write(*,*) 'rank > 5 not implemented in cCgp12345'
          write(*,*) 'rank = ',rank
          stop
        end if
      end if

      call cBp12345(p21,m12,m22,B0w0,B1w0,B2w0,B3w0,B4w0,B5w0,B6w0
     &    ,rank+ordgp3,0)
      call cBp12345(p20,m02,m22,B0w1,B1w1,B2w1,B3w1,B4w1,B5w1,B6w1
     &    ,rank+ordgp3,0)
      call cBp12345(p10,m02,m12,B0w2,B1w2,B2w2,B3w2,B4w2,B5w2,B6w2
     &    ,rank+ordgp3,0)

      mm02 = elimcminf2(m02)
      mm12 = elimcminf2(m12)
      mm22 = elimcminf2(m22)
      q10  = elimcminf2(p10)
      q21  = elimcminf2(p21)
      q20  = elimcminf2(p20)
 
      k1k2 = (q10+q20-q21)
      detz = 4d0*q10*q20-k1k2*k1k2

      if (abs(detz/( 4d0*q10*q20 + k1k2*k1k2)).lt.1d-4) then
        if (abs(q10-q20).lt.abs(q10-q21).and.
     &      abs(q10-q20).lt.abs(q20-q21)) then
          detz  =  4d0*q10*q21 - (q10-q20+q21)*(q10-q20+q21)
        end if
      end if

c      write(testout,*) 'detz = ',detz/(4d0*q10*q20+k1k2*k1k2)
      z(1,1) = 2d0*q10
      z(1,2) = k1k2
      z(2,1) = k1k2
      z(2,2) = 2d0*q20
      zadj(1,1) = 2d0*q20
      zadj(1,2) = -k1k2
      zadj(2,1) = -k1k2
      zadj(2,2) = 2d0*q10
      f(1) = q10+mm02-mm12
      f(2) = q20+mm02-mm22

      zadjf(1) = zadj(1,1)*f(1)+zadj(1,2)*f(2)
      zadjf(2) = zadj(2,1)*f(1)+zadj(2,2)*f(2)

c      write(testout,*) zadj(1,1), zadj(1,2), zadj(2,1), zadj(2,2)
c      write(testout,*) f(1),f(2)
c      write(testout,*) zadjf(1),zadjf(2)
      
c>      C1(1) = 0d0
c>      C1(2) = 0d0
c>
c>      C2(1,1) = 0d0
c>      C2(1,2) = 0d0
c>      C2(2,1) = 0d0
c>      C2(2,2) = 0d0
c>
c>      C3(1,1,1) = 0d0
c>      C3(1,1,2) = 0d0
c>      C3(1,2,1) = 0d0
c>      C3(2,1,1) = 0d0
c>      C3(1,2,2) = 0d0
c>      C3(2,1,2) = 0d0
c>      C3(2,2,1) = 0d0
c>      C3(2,2,2) = 0d0
c>
c>      C4(0,0,1,1) = 0d0
c>      C4(0,0,1,2) = 0d0
c>      C4(0,0,2,1) = 0d0
c>      C4(0,0,2,2) = 0d0
c>      C4(1,1,1,1) = 0d0
c>      C4(1,1,1,2) = 0d0
c>      C4(1,1,2,1) = 0d0
c>      C4(1,2,1,1) = 0d0
c>      C4(2,1,1,1) = 0d0
c>      C4(1,1,2,2) = 0d0
c>      C4(1,2,1,2) = 0d0
c>      C4(2,1,1,2) = 0d0
c>      C4(1,2,2,1) = 0d0
c>      C4(2,1,2,1) = 0d0
c>      C4(2,2,1,1) = 0d0
c>      C4(1,2,2,2) = 0d0
c>      C4(2,1,2,2) = 0d0
c>      C4(2,2,1,2) = 0d0
c>      C4(2,2,2,1) = 0d0
c>      C4(2,2,2,2) = 0d0

      do i1=1,2
        C1(i1) = 0d0
        C3(0,0,i1) = 0d0
        C5(0,0,0,0,i1) = 0d0
      do i2=1,2
        C2(i1,i2) = 0d0
        C4(0,0,i1,i2) = 0d0
        C6(0,0,0,0,i1,i2) = 0d0
      do i3=1,2
        C3(i1,i2,i3) = 0d0
        C5(0,0,i1,i2,i3) = 0d0
      do i4=1,2
        C4(i1,i2,i3,i4) = 0d0
        C6(0,0,i1,i2,i3,i4) = 0d0
      do i5=1,2
        C5(i1,i2,i3,i4,i5) = 0d0
c      do i6=1,2
c        C6(i1,i2,i3,i4,i5,i6) = 0d0
c      end do
      end do
      end do
      end do
      end do
      end do

      r1 = abs(z(1,1)/z(1,2))
      r2 = abs(z(2,2)/z(2,1))

      B0m(0) = B0w0
      B0m(1) = B0w1
      B0m(2) = B0w2

      S1hat(1) = B0m(1) - B0m(0)
      S1hat(2) = B0m(2) - B0m(0)

      if (rank+ordgp3.eq.0) goto 99

      B1m(0,2) = B1w0
      B1m(0,1) = -B0m(0) - B1m(0,2)
      B1m(1,2) = B1w1
      B1m(2,1) = B1w2
        
      S2hat(1,1) =          - B1m(0,1)
      S2hat(1,2) = B1m(1,2) - B1m(0,2)
      S2hat(2,1) = B1m(2,1) - B1m(0,1)
      S2hat(2,2) =          - B1m(0,2)
                
      if (rank+ordgp3.eq.1) goto 99

      B2m(0,0,0) = B2w0(0,0)
      B2m(1,0,0) = B2w1(0,0)
      B2m(2,0,0) = B2w2(0,0)
      B2m(0,2,2) = B2w0(1,1)
      B2m(0,1,2) = -B1m(0,2) - B2m(0,2,2)
      B2m(0,1,1) = -B1m(0,1) - B2m(0,1,2)
      B2m(1,2,2) = B2w1(1,1)
      B2m(2,1,1) = B2w2(1,1)

      B2m(0,2,1) =  B2m(0,1,2)
 
      S3hat(1,0,0) = B2m(1,0,0) - B2m(0,0,0)
      S3hat(2,0,0) = B2m(2,0,0) - B2m(0,0,0)
      S3hat(1,1,1) =            - B2m(0,1,1)
      S3hat(1,1,2) =            - B2m(0,1,2)
      S3hat(1,2,1) =            - B2m(0,1,2)
      S3hat(1,2,2) = B2m(1,2,2) - B2m(0,2,2)
      S3hat(2,1,1) = B2m(2,1,1) - B2m(0,1,1)
      S3hat(2,1,2) =            - B2m(0,1,2)
      S3hat(2,2,1) =            - B2m(0,1,2)
      S3hat(2,2,2) =            - B2m(0,2,2)

      if (rank+ordgp3.eq.2) goto 99

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
      
      B3m(0,2,1,2) = B3m(0,1,2,2) 
cs      B3m(0,2,2,1) = B3m(0,1,2,2) 
cs      B3m(0,1,2,1) = B3m(0,1,1,2) 
      B3m(0,2,1,1) = B3m(0,1,1,2) 

      S4hat(1,0,0,1) =              - B3m(0,0,0,1)
      S4hat(2,0,0,1) = B3m(2,0,0,1) - B3m(0,0,0,1)
      S4hat(1,0,0,2) = B3m(1,0,0,2) - B3m(0,0,0,2)
      S4hat(2,0,0,2) =              - B3m(0,0,0,2)
        
      S4hat(1,1,1,1) =              - B3m(0,1,1,1)
cs      S4hat(1,1,2,1) =              - B3m(0,1,1,2)
      S4hat(1,2,1,1) =              - B3m(0,1,1,2)
cs      S4hat(1,2,2,1) =              - B3m(0,1,2,2)
      S4hat(2,1,1,1) = B3m(2,1,1,1) - B3m(0,1,1,1)
cs      S4hat(2,1,2,1) =              - B3m(0,1,1,2)
      S4hat(2,2,1,1) =              - B3m(0,1,1,2)
cs      S4hat(2,2,2,1) =              - B3m(0,1,2,2)
      S4hat(1,1,1,2) =              - B3m(0,1,1,2)
      S4hat(1,1,2,2) =              - B3m(0,1,2,2)
      S4hat(1,2,1,2) =              - B3m(0,1,2,2)
      S4hat(1,2,2,2) = B3m(1,2,2,2) - B3m(0,2,2,2)
      S4hat(2,1,1,2) =              - B3m(0,1,1,2)
      S4hat(2,1,2,2) =              - B3m(0,1,2,2)
      S4hat(2,2,1,2) =              - B3m(0,1,2,2)
      S4hat(2,2,2,2) =              - B3m(0,2,2,2)

      if (rank+ordgp3.eq.3) goto 99

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
        

      S5hat(2,0,0,0,0) = B4m(2,0,0,0,0) - B4m(0,0,0,0,0)
      S5hat(1,0,0,0,0) = B4m(1,0,0,0,0) - B4m(0,0,0,0,0)
      do j=1,2
        do i1=1,2
        do i2=1,2
          S5hat(j,0,0,i1,i2) = B4m(j,0,0,i1,i2) - B4m(0,0,0,i1,i2)
          do i3=i2,2
          do i4=i3,2     
            S5hat(j,i1,i2,i3,i4) =             - B4m(0,i1,i2,i3,i4)
          end do
          end do
        end do
        end do
      end do
      S5hat(2,1,1,1,1) = B4m(2,1,1,1,1) - B4m(0,1,1,1,1)
      S5hat(1,2,2,2,2) = B4m(1,2,2,2,2) - B4m(0,2,2,2,2)

c      write(testout,*) 'S5hat ',S5hat(1,0,0,1,2),S5hat(2,0,0,1,2)
   
      if (rank+ordgp3.eq.4) goto 99

      B5m(0,0,0,0,0,2) =  B5w0(0,0,0,0,1)
      B5m(0,0,0,0,0,1) = -B4m(0,0,0,0,0) - B5m(0,0,0,0,0,2)
      B5m(1,0,0,0,0,2) =  B5w1(0,0,0,0,1)
      B5m(2,0,0,0,0,1) =  B5w2(0,0,0,0,1)

      B5m(0,0,0,2,2,2) =  B5w0(0,0,1,1,1)
      B5m(0,0,0,1,2,2) = -B4m(0,0,0,2,2)-B5m(0,0,0,2,2,2)
      B5m(0,0,0,1,1,2) = -B4m(0,0,0,1,2)-B5m(0,0,0,1,2,2)
      B5m(0,0,0,1,1,1) = -B4m(0,0,0,1,1)-B5m(0,0,0,1,1,2)

      B5m(0,2,2,2,2,2) =  B5w0(1,1,1,1,1)
      B5m(0,1,2,2,2,2) = -B4m(0,2,2,2,2)-B5m(0,2,2,2,2,2)
      B5m(0,1,1,2,2,2) = -B4m(0,1,2,2,2)-B5m(0,1,2,2,2,2)
      B5m(0,1,1,1,2,2) = -B4m(0,1,1,2,2)-B5m(0,1,1,2,2,2)
      B5m(0,1,1,1,1,2) = -B4m(0,1,1,1,2)-B5m(0,1,1,1,2,2)
      B5m(0,1,1,1,1,1) = -B4m(0,1,1,1,1)-B5m(0,1,1,1,1,2)

      B5m(1,0,0,2,2,2) =  B5w1(0,0,1,1,1)
      B5m(2,0,0,1,1,1) =  B5w2(0,0,1,1,1)      
      
      B5m(1,2,2,2,2,2) =  B5w1(1,1,1,1,1)
      B5m(2,1,1,1,1,1) =  B5w2(1,1,1,1,1)
      
      B5m(0,0,0,2,1,1) =  B5m(0,0,0,1,1,2) 
cs      B5m(0,0,0,1,2,1) =  B5m(0,0,0,1,1,2) 
      B5m(0,0,0,2,1,2) =  B5m(0,0,0,1,2,2) 
cs      B5m(0,0,0,2,2,1) =  B5m(0,0,0,1,2,2) 

      B5m(0,2,1,1,1,1) =  B5m(0,1,1,1,1,2) 
cs      B5m(0,1,2,1,1,1) =  B5m(0,1,1,1,1,2) 
cs      B5m(0,1,1,2,1,1) =  B5m(0,1,1,1,1,2) 
cs      B5m(0,1,1,1,2,1) =  B5m(0,1,1,1,1,2) 
cs      B5m(0,1,1,2,1,2) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,1,1,2,2,1) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,1,2,1,1,2) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,1,2,1,2,1) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,1,2,2,1,1) =  B5m(0,1,1,1,2,2) 
      B5m(0,2,1,1,1,2) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,2,1,1,2,1) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,2,1,2,1,1) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,2,2,1,1,1) =  B5m(0,1,1,1,2,2) 
cs      B5m(0,1,2,1,2,2) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,1,1,2,2) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,1,2,2,1,2) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,2,1,2,1,2) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,2,2,1,1,2) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,1,2,2,2,1) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,2,1,2,2,1) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,2,2,1,2,1) =  B5m(0,1,1,2,2,2) 
cs      B5m(0,2,2,2,1,1) =  B5m(0,1,1,2,2,2) 
      B5m(0,2,1,2,2,2) =  B5m(0,1,2,2,2,2) 
cs      B5m(0,2,2,1,2,2) =  B5m(0,1,2,2,2,2) 
cs      B5m(0,2,2,2,1,2) =  B5m(0,1,2,2,2,2) 
      B5m(0,2,2,2,2,1) =  B5m(0,1,2,2,2,2) 
        
      do i1=1,2
        do j=1,2
          S6hat(j,0,0,0,0,i1) = B5m(j,0,0,0,0,i1) - B5m(0,0,0,0,0,i1)
          do i2=1,2
          do i3=i2,2
            S6hat(j,0,0,i1,i2,i3) = 
     &          B5m(j,0,0,i1,i2,i3) - B5m(0,0,0,i1,i2,i3)
          do i4=i3,2
          do i5=i4,2     
            S6hat(j,i1,i2,i3,i4,i5) = -B5m(0,i1,i2,i3,i4,i5)
          end do
          end do
          end do
          end do
        end do
      end do
      S6hat(2,1,1,1,1,1) = B5m(2,1,1,1,1,1) - B5m(0,1,1,1,1,1)
      S6hat(1,2,2,2,2,2) = B5m(1,2,2,2,2,2) - B5m(0,2,2,2,2,2)


   99   continue

c choose reduction formulas with biggest denominators
      if (abs(f(1)).ge.abs(f(2))) then
        m = 1
      else 
        m = 2
      end if

c>      do 500 m = 1,2
c>      do 500 k = 1,2
c>      do 500 l = 1,2
c>        kt = 3-k
c>        lt = 3-l
c>        sgn = (-1)**(k-l)
c>        testout = testout+1
c>        write(testout,*) ' m = ',m
c>        write(testout,*) ' k = ',k
c>        write(testout,*) ' l = ',l

c      C0 = cC0f(p10,p21,p20,m02,m12,m22)
c      write(testout,9) 'C0 = ',C0
    

c      write(testout,*) 'ordgp3 = ',ordgp3
      do 100 r=0,ordgp3+rank
c       do 100 r=0,0
c        write(testout,*) 'ordgp3 = ',ordgp3,r,m,f(m)

        C0 = ( S1hat(m) - z(m,1)*C1(1)  - z(m,2)*C1(2))/f(m) 
c        write(testout,9) 'C0 = ',C0
 9      format(1x,a10,2g24.16,2i2)  

        if (rank+ordgp3-r.eq.0) goto 100
        
        S00 = 2d0*mm02*C0 + 2d0*B0m(0)
        
        C2(0,0) = ( S00+1d0
     &         - z(1,1)*C2(1,1) - 2d0*z(1,2)*C2(1,2)
     &         - z(2,2)*C2(2,2)
     &           ) /8d0
c        write(testout,9) 'C00 = ',C2(0,0)

        S2mod(1,1) =  S2hat(1,1) - 2d0*C2(0,0)
        S2mod(1,2) =  S2hat(1,2)
        S2mod(2,1) =  S2hat(2,1)
        S2mod(2,2) =  S2hat(2,2) - 2d0*C2(0,0)
      
        C1(1) = ( S2mod(m,1) - z(m,1)*C2(1,1)  - z(m,2)*C2(2,1))/f(m) 
        C1(2) = ( S2mod(m,2) - z(m,1)*C2(1,2)  - z(m,2)*C2(2,2))/f(m) 
c        write(testout,9) 'C1 = ',C1(1)
c        write(testout,9) 'C2 = ',C1(2)
        
        if (rank+ordgp3-r.eq.1) goto 100
 
        S001(1) = 2d0*mm02*C1(1) + 2d0*B1m(0,1)
        S001(2) = 2d0*mm02*C1(2) + 2d0*B1m(0,2)

        C3(0,0,1) = ( S001(1)-1d0/3d0
     &         - z(1,1)*C3(1,1,1) - 2d0*z(1,2)*C3(1,2,1)
     &         - z(2,2)*C3(2,2,1)
     &           ) /12d0
        C3(0,0,2) = ( S001(2)-1d0/3d0
     &         - z(1,1)*C3(1,1,2) - 2d0*z(1,2)*C3(1,2,2)
     &         - z(2,2)*C3(2,2,2)
     &           ) /12d0

c        write(testout,9) 'C001 = ',C3(0,0,1)
c        write(testout,9) 'C002 = ',C3(0,0,2)
 
        S3mod(1,1,1) =  S3hat(1,1,1) - 4d0*C3(0,0,1)
        S3mod(1,1,2) =  S3hat(1,1,2) - 2d0*C3(0,0,2)
cs        S3mod(1,2,1) =  S3hat(1,2,1) - 2d0*C3(0,0,2)
        S3mod(1,2,2) =  S3hat(1,2,2)
        S3mod(2,1,1) =  S3hat(2,1,1)
        S3mod(2,1,2) =  S3hat(2,1,2) - 2d0*C3(0,0,1)
cs        S3mod(2,2,1) =  S3hat(2,2,1) - 2d0*C3(0,0,1)
        S3mod(2,2,2) =  S3hat(2,2,2) - 4d0*C3(0,0,2)

        C2(1,1) = ( S3mod(m,1,1) 
     &      - z(m,1)*C3(1,1,1)  - z(m,2)*C3(2,1,1))/f(m) 
        C2(1,2) = ( S3mod(m,1,2) 
     &      - z(m,1)*C3(1,1,2)  - z(m,2)*C3(2,1,2))/f(m) 
cs        C2(2,1) = ( S3mod(m,2,1) 
cs     &      - z(m,1)*C3(1,2,1)  - z(m,2)*C3(2,2,1))/f(m) 
        C2(2,2) = ( S3mod(m,2,2) 
     &      - z(m,1)*C3(1,2,2)  - z(m,2)*C3(2,2,2))/f(m) 
        C2(2,1) = C2(1,2)

c        do i1=1,2
c          do i2=i1,2     
c            write(testout,20) 'Cij = ',i1,i2,C2(i1,i2)             
c          end do
c        end do

        if (rank+ordgp3-r.eq.2) goto 100
        
        do i1=1,2
        do i2=i1,2
          S002(i1,i2) = 2d0*mm02*C2(i1,i2) + 2d0*B2m(0,i1,i2)
        end do
        end do
   
        C4(0,0,1,1) = ( S002(1,1)+1d0/6d0
     &         - z(1,1)*C4(1,1,1,1) - 2d0*z(1,2)*C4(1,2,1,1)
     &         - z(2,2)*C4(2,2,1,1)
     &           ) /16d0
        C4(0,0,1,2) = ( S002(1,2)+1d0/12d0
     &         - z(1,1)*C4(1,1,1,2) - 2d0*z(1,2)*C4(1,2,1,2)
     &         - z(2,2)*C4(2,2,1,2)
     &           ) /16d0
        C4(0,0,2,1) = ( S002(2,1)+1d0/12d0
     &         - z(1,1)*C4(1,1,2,1) - 2d0*z(1,2)*C4(1,2,2,1)
     &         - z(2,2)*C4(2,2,2,1)
     &           ) /16d0
        C4(0,0,2,2) = ( S002(2,2)+1d0/6d0
     &         - z(1,1)*C4(1,1,2,2) - 2d0*z(1,2)*C4(1,2,2,2)
     &         - z(2,2)*C4(2,2,2,2)
     &           ) /16d0

        do i1=1,2
          do i2=i1,2     
            C4(0,0,i1,i2) = S002(i1,i2)+1d0/12d0
     &         - z(1,1)*C4(1,1,i1,i2) - 2d0*z(1,2)*C4(1,2,i1,i2)
     &         - z(2,2)*C4(2,2,i1,i2)
            if (i1.eq.i2)  C4(0,0,i1,i2) = C4(0,0,i1,i2) + 1d0/12d0
            C4(0,0,i1,i2) = C4(0,0,i1,i2) / 16d0
c           write(testout,20) 'C00ij = ',i1,i2,C4(0,0,i1,i2)             
          end do
        end do
        C4(0,0,2,1) = C4(0,0,1,2)

        do j=1,2
        do i1=1,2
        do i2=i1,2     
        do i3=i2,2     
          S4mod(j,i1,i2,i3) =  S4hat(j,i1,i2,i3)
          if (j.eq.i1) S4mod(j,i1,i2,i3) =  S4mod(j,i1,i2,i3)
     &        -2d0*C4(0,0,i2,i3)
          if (j.eq.i2) S4mod(j,i1,i2,i3) =  S4mod(j,i1,i2,i3)
     &        -2d0*C4(0,0,i1,i3)
          if (j.eq.i3) S4mod(j,i1,i2,i3) =  S4mod(j,i1,i2,i3)
     &        -2d0*C4(0,0,i1,i2)
        end do
        end do
        end do
        end do
 
        do i1=1,2
        do i2=i1,2     
        do i3=i2,2     
          C3(i1,i2,i3) = (S4mod(m,i1,i2,i3) 
     &      - z(m,1)*C4(1,i1,i2,i3)  - z(m,2)*C4(2,i1,i2,i3) 
     &        )/f(m)
c          write(testout,30) 'Cijk = ',i1,i2,i3,C3(i1,i2,i3)
        end do 
        end do 
        end do 

        C3(2,1,1) = C3(1,1,2)
        C3(1,2,1) = C3(1,1,2)
        C3(2,1,2) = C3(1,2,2)
        C3(2,2,1) = C3(1,2,2)
       
        if (rank+ordgp3-r.eq.3) goto 100

        S002(0,0) = 2d0*mm02*C2(0,0) + 2d0*B2m(0,0,0)
        
        C4(0,0,0,0) = (S002(0,0)+(mm02+mm12+mm22)/6d0
     &                        -(q10+q20+q21)/24d0
     &         - z(1,1)*C4(0,0,1,1) - 2d0*z(1,2)*C4(0,0,1,2)
     &         - z(2,2)*C4(0,0,2,2)
     &           ) /12d0
c        write(testout,9) 'C0000 = ',C4(0,0,0,0) 

        do i1=1,2
        do i2=i1,2     
        do i3=i2,2     
            C5(0,0,i1,i2,i3) = S003(i1,i2,i3)-1d0/30d0
     &         - z(1,1)*C5(1,1,i1,i2,i3) - 2d0*z(1,2)*C5(1,2,i1,i2,i3)
     &         - z(2,2)*C5(2,2,i1,i2,i3)
            if (i1.eq.i2.and.i1.eq.i3)  
     &          C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) - 1d0/15d0
            C5(0,0,i1,i2,i3) = C5(0,0,i1,i2,i3) / 20d0
c          write(testout,30) 'C00ijk = ',i1,i2,i3,C5(0,0,i1,i2,i3)             
 30       format(1x,a10,3i2,2g24.16)  
        end do
        end do
        end do
        C5(0,0,2,1,1) = C5(0,0,1,1,2)
        C5(0,0,1,2,1) = C5(0,0,1,1,2)
        C5(0,0,2,2,1) = C5(0,0,1,2,2)
        C5(0,0,2,1,2) = C5(0,0,1,2,2)

        do j=1,2
        do i1=1,2
        do i2=i1,2     
        do i3=i2,2     
        do i4=i3,2     
          S5mod(j,i1,i2,i3,i4) =  S5hat(j,i1,i2,i3,i4)
          if (j.eq.i1) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i2,i3,i4)
          if (j.eq.i2) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i1,i3,i4)
          if (j.eq.i3) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i1,i2,i4)
          if (j.eq.i4) S5mod(j,i1,i2,i3,i4) =  S5mod(j,i1,i2,i3,i4)
     &        -2d0*C5(0,0,i1,i2,i3)
        end do
        end do
        end do
        end do
        end do
 
        do i1=1,2
        do i2=i1,2     
        do i3=i2,2
        do i4=i3,2
          C4(i1,i2,i3,i4) = (S5mod(m,i1,i2,i3,i4) 
     &      - z(m,1)*C5(1,i1,i2,i3,i4)  - z(m,2)*C5(2,i1,i2,i3,i4) 
     &        )/f(m)
c          write(testout,40) 'Cijkl = ',i1,i2,i3,i4,C4(i1,i2,i3,i4)
        end do 
        end do 
        end do 
        end do 
 
        C4(1,2,1,1) = C4(1,1,1,2)
        C4(2,1,1,1) = C4(1,1,1,2)
        C4(2,1,1,2) = C4(1,1,2,2)
        C4(1,2,1,2) = C4(1,1,2,2)
        C4(2,2,1,1) = C4(1,1,2,2)
        C4(2,1,2,2) = C4(1,2,2,2)
        C4(2,2,1,2) = C4(1,2,2,2)

        if (rank+ordgp3-r.eq.4) goto 100


        C5(0,0,0,0,1) = (-5d0*mm02 - 10d0*mm12 - 5d0*mm22 
     &      + q20 + 2d0*q10 + 2d0*q21)
     &      /(120d0)
        C5(0,0,0,0,2) = (-5d0*mm02 - 10d0*mm22 - 5d0*mm12 
     &      + q10 + 2d0*q20 + 2d0*q21)
     &      /(120d0)

        do i1=1,2
          S003(0,0,i1) = 2d0*mm02*C3(0,0,i1) + 2d0*B3m(0,0,0,i1)          
        do i2=i1,2
        do i3=i2,2
          S003(i1,i2,i3) = 2d0*mm02*C3(i1,i2,i3) + 2d0*B3m(0,i1,i2,i3)
        end do
        end do

        C5(0,0,0,0,i1) =  C5(0,0,0,0,i1) +
     &      S003(0,0,i1)
     &         - z(1,1)*C5(0,0,1,1,i1) - 2d0*z(1,2)*C5(0,0,1,2,i1)
     &         - z(2,2)*C5(0,0,2,2,i1)
        C5(0,0,0,0,i1) =  C5(0,0,0,0,i1) / 16d0

c        write(testout,10) 'C0000i = ',i1,C5(0,0,0,0,i1)
 10       format(1x,a10,i2,2g24.16,2i2)  
        end do

        do i1=1,2
        do i2=i1,2     
        do i3=i2,2
        do i4=i3,2
          S004(i1,i2,i3,i4) = 2d0*mm02*C4(i1,i2,i3,i4) 
     &        + 2d0*B4m(0,i1,i2,i3,i4)
        end do
        end do
        end do
        end do

        do i1=1,2
        do i2=i1,2     
        do i3=i2,2     
        do i4=i3,2     
          if (i1+i2+i3+i4.eq.4.or.i1+i2+i3+i4.eq.8) then  
            C6(0,0,i1,i2,i3,i4) = 1d0/15d0
          else if (i1+i2+i3+i4.eq.5.or.i1+i2+i3+i4.eq.7) then  
            C6(0,0,i1,i2,i3,i4) = 1d0/60d0
          else if (i1+i2+i3+i4.eq.6) then  
            C6(0,0,i1,i2,i3,i4) = 1d0/90d0
          end if
            C6(0,0,i1,i2,i3,i4) =  C6(0,0,i1,i2,i3,i4) 
     &        + S004(i1,i2,i3,i4)-1d0/30d0
c     &        - z(1,1)*C6(1,1,i1,i2,i3,i4) 
c     &        - 2d0*z(1,2)*C6(1,2,i1,i2,i3,i4)
c     &        - z(2,2)*C6(2,2,i1,i2,i3,i4)
            C6(0,0,i1,i2,i3,i4) = C6(0,0,i1,i2,i3,i4) / 24d0
c          write(testout,40) 'C00ijkl = ',i1,i2,i3,i4,C6(0,0,i1,i2,i3,i4)
 40       format(1x,a10,4i2,2g24.16)  
      end do
      end do
      end do
      end do

      do j=1,2
      do i1=1,2
      do i2=i1,2     
      do i3=i2,2     
      do i4=i3,2     
      do i5=i4,2     
        S6mod(j,i1,i2,i3,i4,i5) =  S6hat(j,i1,i2,i3,i4,i5)
        if (j.eq.i1) S6mod(j,i1,i2,i3,i4,i5) =  S6mod(j,i1,i2,i3,i4,i5)
     &      -2d0*C6(0,0,i2,i3,i4,i5)
        if (j.eq.i2) S6mod(j,i1,i2,i3,i4,i5) =  S6mod(j,i1,i2,i3,i4,i5)
     &      -2d0*C6(0,0,i1,i3,i4,i5)
        if (j.eq.i3) S6mod(j,i1,i2,i3,i4,i5) =  S6mod(j,i1,i2,i3,i4,i5)
     &      -2d0*C6(0,0,i1,i2,i4,i5) 
        if (j.eq.i4) S6mod(j,i1,i2,i3,i4,i5) =  S6mod(j,i1,i2,i3,i4,i5)
     &      -2d0*C6(0,0,i1,i2,i3,i5)
        if (j.eq.i5) S6mod(j,i1,i2,i3,i4,i5) =  S6mod(j,i1,i2,i3,i4,i5)
     &      -2d0*C6(0,0,i1,i2,i3,i4)
      end do
      end do
      end do
      end do
      end do
      end do
 
      do i1=1,2 
      do i2=i1,2     
      do i3=i2,2
      do i4=i3,2
      do i5=i4,2
          C5(i1,i2,i3,i4,i5) = (S6mod(m,i1,i2,i3,i4,i5) 
c     &      - z(m,1)*C6(1,i1,i2,i3,i4,i5)  - z(m,2)*C6(2,i1,i2,i3,i4,i5) 
     &        )/f(m)
c        write(testout,50) 'Cijklm = ',i1,i2,i3,i4,i5,C5(i1,i2,i3,i4,i5)
 50       format(1x,a10,5i2,2g24.16)  
      end do 
      end do 
      end do 
      end do 
      end do 
      C5(2,1,1,1,1) = C5(1,1,1,1,2)
      C5(1,2,1,1,1) = C5(1,1,1,1,2)
      C5(2,1,1,1,2) = C5(1,1,1,2,2)
      C5(1,2,1,1,2) = C5(1,1,1,2,2)
      C5(2,2,1,1,1) = C5(1,1,1,2,2)
      C5(2,1,1,2,2) = C5(1,1,2,2,2)
      C5(1,2,1,2,2) = C5(1,1,2,2,2)
      C5(2,2,1,1,2) = C5(1,1,2,2,2)
      C5(2,1,2,2,2) = C5(1,2,2,2,2)
      C5(2,2,1,2,2) = C5(1,2,2,2,2)

        if (rank+ordgp3-r.eq.5) goto 100

c C6 incomplete!!!

        S004(0,0,0,0) = 2d0*mm02*C4(0,0,0,0) + 2d0*B4m(0,0,0,0,0)
        
        C6(0,0,0,0,0,0) = (S004(0,0,0,0)
     &      +(15d0*(mm02*(mm02+mm12)+mm22*(mm22+mm02)+mm12*(mm12+mm22))
     &       - 3d0*(mm02*q21+mm12*q20+mm22*q10)
     &       - 6d0*(mm02*(q10+q20)+mm12*(q21+q10)+mm22*(q20+q21))
     &      +  q10*(q10+q20)+q21*(q21+q10)+q20*(q20+q21))/720d0 
     &         - z(1,1)*C6(0,0,0,0,1,1) - 2d0*z(1,2)*C6(0,0,0,0,1,2)
     &         - z(2,2)*C6(0,0,0,0,2,2)
     &           ) / 16d0

c        write(testout,*) 'C000000 = ',C6(0,0,0,0,0,0) 

        C6(0,0,0,0,1,1) = (6d0*(mm02+mm22) + 18d0*mm12
     &      - q20 - 3d0*(q10+q21) )
     &      /(360d0)
        C6(0,0,0,0,1,2) = (3d0*mm02 + 6d0*(mm22+mm12) 
     &      - q10 - q20 - 2d0*q21)
     &      /(360d0)
        C6(0,0,0,0,2,1) = C6(0,0,0,0,1,2)
        C6(0,0,0,0,2,2) = (6d0*(mm02+mm12) + 18d0*mm22 
     &      - q10 - 3d0*(q20+q21) )
     &      /(360d0)

        do i1=1,2
        do i2=i1,2
          S004(0,0,i1,i2) = 2d0*mm02*C4(0,0,i1,i2) 
     &        + 2d0*B4m(0,0,0,i1,i2)

        C6(0,0,0,0,i1,i2) =  C6(0,0,0,0,i1,i2) +
     &           S004(0,0,i1,i2)
     &         - z(1,1)*C6(0,0,1,1,i1,i2) 
     &         - 2d0*z(1,2)*C6(0,0,1,2,i1,i2)
     &         - z(2,2)*C6(0,0,2,2,i1,i2)

        C6(0,0,0,0,i1,i2) =  C6(0,0,0,0,i1,i2) / 20d0

c        write(testout,20) 'C0000ij = ',i1,i2,C6(0,0,0,0,i1,i2)
 20     format(1x,a10,2i2,2g24.16,2i2)  

        end do
        end do

 100   continue  

      if(sym.eq.1) then
c        if (rank.gt.1) then
c          C2(2,1) = C2(1,2)
c        if (rank.gt.2) then
c          C3(1,2,1) = C3(1,1,2)
c          C3(2,1,1) = C3(1,1,2)
c          C3(2,1,2) = C3(1,2,2)
c          C3(2,2,1) = C3(1,2,2)
        if (rank.gt.3) then
c          C4(0,0,2,1) = C4(0,0,1,2)
          C4(1,1,2,1) = C4(1,1,1,2)
c          C4(1,2,1,1) = C4(1,1,1,2)
c          C4(2,1,1,1) = C4(1,1,1,2)
c          C4(1,2,1,2) = C4(1,1,2,2)
c          C4(2,1,1,2) = C4(1,1,2,2)
          C4(1,2,2,1) = C4(1,1,2,2)
          C4(2,1,2,1) = C4(1,1,2,2)
c          C4(2,2,1,1) = C4(1,1,2,2)
c          C4(2,1,2,2) = C4(1,2,2,2)
c          C4(2,2,1,2) = C4(1,2,2,2)
          C4(2,2,2,1) = C4(1,2,2,2)
        if (rank.gt.4) then
c          C5(0,0,1,2,1) = C5(0,0,1,1,2)
c          C5(0,0,2,1,1) = C5(0,0,1,1,2)
c          C5(0,0,2,1,2) = C5(0,0,1,2,2)
c          C5(0,0,2,2,1) = C5(0,0,1,2,2)
          C5(1,1,1,2,1) = C5(1,1,1,1,2)
          C5(1,1,2,1,1) = C5(1,1,1,1,2)
c          C5(1,2,1,1,1) = C5(1,1,1,1,2)
c          C5(2,1,1,1,1) = C5(1,1,1,1,2)
c          C5(2,1,1,1,2) = C5(2,1,1,1,2)
          C5(1,1,2,1,2) = C5(1,1,1,2,2)
c          C5(1,2,1,1,2) = C5(1,1,1,2,2)
          C5(1,1,2,2,1) = C5(1,1,1,2,2)
          C5(1,2,1,2,1) = C5(1,1,1,2,2)
          C5(2,1,1,2,1) = C5(1,1,1,2,2)
          C5(1,2,2,1,1) = C5(1,1,1,2,2)
          C5(2,1,2,1,1) = C5(1,1,1,2,2)
c          C5(2,2,1,1,1) = C5(1,1,1,2,2)
c          C5(1,2,1,2,2) = C5(1,1,2,2,2)
c          C5(2,1,1,2,2) = C5(1,1,2,2,2)
          C5(1,2,2,1,2) = C5(1,1,2,2,2)
          C5(2,1,2,1,2) = C5(1,1,2,2,2)
c          C5(2,2,1,1,2) = C5(1,1,2,2,2)
          C5(1,2,2,2,1) = C5(1,1,2,2,2)
          C5(2,1,2,2,1) = C5(1,1,2,2,2)
          C5(2,2,1,2,1) = C5(1,1,2,2,2)
          C5(2,2,2,1,1) = C5(1,1,2,2,2)
c          C5(2,1,2,2,2) = C5(1,2,2,2,2)
c          C5(2,2,1,2,2) = C5(1,2,2,2,2)
          C5(2,2,2,1,2) = C5(1,2,2,2,2)
          C5(2,2,2,2,1) = C5(1,2,2,2,2)
        if (rank.gt.5) then
c          C6(0,0,0,0,2,1) = C6(0,0,0,0,1,2)
          C6(0,0,1,1,2,1) = C6(0,0,1,1,1,2)
          C6(0,0,1,2,1,1) = C6(0,0,1,1,1,2)
          C6(0,0,2,1,1,1) = C6(0,0,1,1,1,2)
          C6(0,0,1,2,1,2) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,1,2) = C6(0,0,1,1,2,2)
          C6(0,0,1,2,2,1) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,2,1) = C6(0,0,1,1,2,2)
          C6(0,0,2,1,2,2) = C6(0,0,1,2,2,2)
          C6(0,0,2,2,1,2) = C6(0,0,1,2,2,2)
          C6(0,0,2,2,2,1) = C6(0,0,1,2,2,2)
        end if
        end if
        end if
c        end if
c        end if
      end if      

c 500  continue

      end

***********************************************************************
      subroutine cDp12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rankin,switchin)
***********************************************************************
*     4-point tensor coefficient functions D0,D1,D2,D3,D4,D5          *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
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
*     16.11.04 Ansgar Denner    last changed  20.05.05                *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p32,p30,p20,p31
      complex*16 m02,m12,m22,m32
      complex*16 D0,D1(3),D2(0:3,0:3),D3(0:3,0:3,0:3)
      complex*16 D4(0:3,0:3,0:3,0:3),D5(0:3,0:3,0:3,0:3,0:3)
      complex*16 Dt0,Dt1(3),Dt2(0:3,0:3),Dt3(0:3,0:3,0:3)
      complex*16 Dt4(0:3,0:3,0:3,0:3),Dt5(0:3,0:3,0:3,0:3,0:3)
c     complex*16 cD0f
      complex*16 zadj(3,3),zadjadjff(3,3),xadj(3,3),z(3,3)
      complex*16 k1k2,k1k3,k2k3
      complex*16 f1,f2,f3,zadjf1,zadjf2,zadjf3
      real*8     detz,rm02,h,mzadjadjff
      real*8     mf,mp2,mm2,mpm2,norm
      real*8     adetz,mzadjf,mzadjfd,mzadj,mxadj,adety
      real*8     accstop,acctest,accinf,accC,accworst
      real*8     accpv,accg,accgy,accy,accgp,accgy2 
      integer    rank,switch,rankin,switchin
      integer    rankacc,ordg4acc,ordgy4acc,ordgp4acc
      integer    ltest,sym
      integer    testout,errout,i1,i2,i3,i4,i5,k,atest,ordg,ordgy,lc
      integer    flag,flag1,count,counttest
      integer    Dcount(0:40)
      complex*16 q10,q21,q32,q30,q20,q31
      complex*16 mm02,mm12,mm22,mm32
      complex*16 elimcminf2
      logical    errwrite

      common /sym/    sym
      common /ltest/  ltest
      integer    version3,version4
      common /version/ version3,version4
      integer    ordg3,ordgy3,ordgp3
      integer    ordg4,ordgy4,ordgp4
      common /ctord3/ ordg3,ordgy3,ordgp3
      common /ctord4/ ordg4,ordgy4,ordgp4

      real*8     acc
c 16.02.05: accstop changed from 1d-2 to 1d-1
c 21.02.05: accstop changed from 1d0  to 1d40
c 23.01.08: accstop changed from 1d40  to 1d0
      data acc /1d-4/, acctest /1d-4/, accstop /1d0/, accC /1d-13/
c      data acc /1d-4/, acctest /1d-4/, accstop /1d40/, accC /1d-13/
c      data acc /1d-4/, acctest /1d-4/, accstop /1d-2/, accC /1d-13/
c      data acc /1d-14/, acctest /1d-4/, accstop /1d-2/, accC /1d-13/
c      data ordg /2/,ordgy /0/
      data accworst /1d-2/, accinf /1d40/
      data testout /14/, errout /94/
      data ordg /1/,ordgy /1/
c      data version4 /0/
      data flag /0/, flag1 /0/, count /0/, counttest /0/
      save flag,flag1,accworst

      common /Dcount/ Dcount

      real*8     accbad1,accbad2,accbad3
      common /accbad/ accbad1,accbad2,accbad3

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
c>      common/cacheh/nf2,nf3
c>      integer nf2(maxt),nf3(maxt)
c>
c>      integer countD
c>      data countD /0/
c>      common /countD/ countD
c>
c>      if (cachelevel.ge.2)then
c>      countD = countD+1
c>      end if

c      sym = 1

      x(1)=p10
      x(2)=p21
      x(3)=p32
      x(4)=p30
      x(5)=p20
      x(6)=p31
      x(7)=m02
      x(8)=m12
      x(9)=m22
      x(10)=m32
      type = 2
      rank = rankin
      switch = switchin
      name='cDp12345'

c      rank = 3

      call cacheread(fct,x,80,10,type,rank,switch,name,nocalc)
      if(nocalc)then
        D0=fct(1)
        if (rank.eq.0) goto 511
        D1(1)=fct(2)
        D1(2)=fct(3)
        D1(3)=fct(4)
        if (rank.eq.1.and.switch.eq.0) goto 511
        D2(0,0)=fct(5)
        if (rank.eq.1) goto 511
        D2(1,1)=fct(6)
        D2(1,2)=fct(7)
        D2(1,3)=fct(8)
        D2(2,2)=fct(9)
        D2(2,3)=fct(10)
        D2(3,3)=fct(11)
        if (rank.eq.2.and.switch.eq.0) goto 511
        D3(0,0,1)=fct(12)
        D3(0,0,2)=fct(13)
        D3(0,0,3)=fct(14)
        if (rank.eq.2) goto 511
        D3(1,1,1)=fct(15)
        D3(1,1,2)=fct(16)
        D3(1,1,3)=fct(17)
        D3(1,2,2)=fct(18)
        D3(1,2,3)=fct(19)
        D3(1,3,3)=fct(20)
        D3(2,2,2)=fct(21)
        D3(2,2,3)=fct(22)
        D3(2,3,3)=fct(23)
        D3(3,3,3)=fct(24)
        if (rank.eq.3.and.switch.eq.0) goto 511
        D4(0,0,0,0)=fct(25)
        D4(0,0,1,1)=fct(26)
        D4(0,0,1,2)=fct(27)
        D4(0,0,1,3)=fct(28)
        D4(0,0,2,2)=fct(29)
        D4(0,0,2,3)=fct(30)
        D4(0,0,3,3)=fct(31)
        if (rank.eq.3) goto 511
        D4(1,1,1,1)=fct(32)
        D4(1,1,1,2)=fct(33)
        D4(1,1,1,3)=fct(34)
        D4(1,1,2,2)=fct(35)
        D4(1,1,2,3)=fct(36)
        D4(1,1,3,3)=fct(37)
        D4(1,2,2,2)=fct(38)
        D4(1,2,2,3)=fct(39)
        D4(1,2,3,3)=fct(40)
        D4(1,3,3,3)=fct(41)
        D4(2,2,2,2)=fct(42)
        D4(2,2,2,3)=fct(43)
        D4(2,2,3,3)=fct(44)
        D4(2,3,3,3)=fct(45)
        D4(3,3,3,3)=fct(46)
        if (rank.eq.4.and.switch.eq.0) goto 511
        D5(0,0,0,0,1)=fct(47)
        D5(0,0,0,0,2)=fct(48)
        D5(0,0,0,0,3)=fct(49)
        D5(0,0,1,1,1)=fct(50)
        D5(0,0,1,1,2)=fct(51)
        D5(0,0,1,1,3)=fct(52)
        D5(0,0,1,2,2)=fct(53)
        D5(0,0,1,2,3)=fct(54)
        D5(0,0,1,3,3)=fct(55)
        D5(0,0,2,2,2)=fct(56)
        D5(0,0,2,2,3)=fct(57)
        D5(0,0,2,3,3)=fct(58)
        D5(0,0,3,3,3)=fct(59)
        if (rank.eq.4) goto 511
        D5(1,1,1,1,1)=fct(60)
        D5(1,1,1,1,2)=fct(61)
        D5(1,1,1,1,3)=fct(62)
        D5(1,1,1,2,2)=fct(63)
        D5(1,1,1,2,3)=fct(64)
        D5(1,1,1,3,3)=fct(65)
        D5(1,1,2,2,2)=fct(66)
        D5(1,1,2,2,3)=fct(67)
        D5(1,1,2,3,3)=fct(68)
        D5(1,1,3,3,3)=fct(69)
        D5(1,2,2,2,2)=fct(70)
        D5(1,2,2,2,3)=fct(71)
        D5(1,2,2,3,3)=fct(72)
        D5(1,2,3,3,3)=fct(73)
        D5(1,3,3,3,3)=fct(74)
        D5(2,2,2,2,2)=fct(75)
        D5(2,2,2,2,3)=fct(76)
        D5(2,2,2,3,3)=fct(77)
        D5(2,2,3,3,3)=fct(78)
        D5(2,3,3,3,3)=fct(79)
        D5(3,3,3,3,3)=fct(80)

 511    continue

        if(sym.eq.1) then
      
          if(rank.le.1) goto 522

          do i1=1,3
            do i2=i1+1,3
              D2(i2,i1)=D2(i1,i2) 
            enddo
          enddo
          
          if (rank.eq.2) goto 522
          
          do i1=1,3
            do i2=i1,3
              do i3=i2,3
                D3(i1,i3,i2)=D3(i1,i2,i3) 
                D3(i2,i1,i3)=D3(i1,i2,i3) 
                D3(i2,i3,i1)=D3(i1,i2,i3) 
                D3(i3,i1,i2)=D3(i1,i2,i3) 
                D3(i3,i2,i1)=D3(i1,i2,i3) 
              enddo
            enddo
          enddo
          
          if (rank.eq.3.and.switch.eq.0) goto 522
          
          do i1=1,3
            do i2=i1+1,3
              D4(0,0,i2,i1)=D4(0,0,i1,i2) 
            enddo
          enddo
          
          if (rank.eq.3) goto 522
          
          do i1=1,3
            do i2=i1,3
              do i3=i2,3
                do i4=i3,3
                  D4(i1,i2,i4,i3)=D4(i1,i2,i3,i4) 
                  D4(i1,i3,i2,i4)=D4(i1,i2,i3,i4) 
                  D4(i1,i3,i4,i2)=D4(i1,i2,i3,i4) 
                  D4(i1,i4,i2,i3)=D4(i1,i2,i3,i4) 
                  D4(i1,i4,i3,i2)=D4(i1,i2,i3,i4) 
                  D4(i2,i1,i3,i4)=D4(i1,i2,i3,i4) 
                  D4(i2,i1,i4,i3)=D4(i1,i2,i3,i4) 
                  D4(i2,i3,i1,i4)=D4(i1,i2,i3,i4) 
                  D4(i2,i3,i4,i1)=D4(i1,i2,i3,i4) 
                  D4(i2,i4,i1,i3)=D4(i1,i2,i3,i4) 
                  D4(i2,i4,i3,i1)=D4(i1,i2,i3,i4) 
                  D4(i3,i1,i2,i4)=D4(i1,i2,i3,i4) 
                  D4(i3,i1,i4,i2)=D4(i1,i2,i3,i4) 
                  D4(i3,i2,i1,i4)=D4(i1,i2,i3,i4) 
                  D4(i3,i2,i4,i1)=D4(i1,i2,i3,i4) 
                  D4(i3,i4,i1,i2)=D4(i1,i2,i3,i4) 
                  D4(i3,i4,i2,i1)=D4(i1,i2,i3,i4) 
                  D4(i4,i1,i2,i3)=D4(i1,i2,i3,i4) 
                  D4(i4,i1,i3,i2)=D4(i1,i2,i3,i4) 
                  D4(i4,i2,i1,i3)=D4(i1,i2,i3,i4) 
                  D4(i4,i2,i3,i1)=D4(i1,i2,i3,i4) 
                  D4(i4,i3,i1,i2)=D4(i1,i2,i3,i4) 
                  D4(i4,i3,i2,i1)=D4(i1,i2,i3,i4) 
                enddo
              enddo
            enddo
          enddo
          
          if (rank.eq.4.and.switch.eq.0) goto 522
          
          do i1=1,3
            do i2=i1,3
              do i3=i2,3
                D5(0,0,i1,i3,i2)=D5(0,0,i1,i2,i3) 
                D5(0,0,i2,i1,i3)=D5(0,0,i1,i2,i3) 
                D5(0,0,i2,i3,i1)=D5(0,0,i1,i2,i3) 
                D5(0,0,i3,i1,i2)=D5(0,0,i1,i2,i3) 
                D5(0,0,i3,i2,i1)=D5(0,0,i1,i2,i3) 
              enddo
            enddo
          enddo
          
          if (rank.eq.4) goto 522
          
          do i1=1,3
            do i2=i1,3
              do i3=i2,3
                do i4=i3,3
                  do k =1,i1
c                   D5(i1,i2,i3,i4,k) = D5(k,i1,i2,i3,i4)
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
                  end do
                end do
              end do
            end do
          end do
        end if
        goto 522  
      endif

      Dcount(0) = Dcount(0) + 1

c      rankacc = 3
c      ordg4acc = ordg
c      ordgy4acc = ordgy

      ordg4  = min(4-rank,ordg)
      ordgy4 = min((5-rank)/2,ordgy)
      ordgp4  = min(4-rank,ordg)
      lc = 0

      rankacc = rank
      ordg4acc = ordg4
      ordgy4acc = min(ordgy4,(4-rank)/2)
      ordgp4acc = ordgp4

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

      f1 = q10-mm12+mm02
      f2 = q20-mm22+mm02
      f3 = q30-mm32+mm02
      k1k2 = (q10+q20-q21)
      k1k3 = (q10+q30-q31)
      k2k3 = (q20+q30-q32)
      detz = (8d0*q10*q30*q20+2d0*k1k2*k1k3*k2k3
     &    -2d0*(q10*k2k3*k2k3+q20*k1k3*k1k3+q30*k1k2*k1k2))
      adetz = abs(detz)
      mp2  = (dmax1(abs(q10),abs(q20),abs(q21),
     &    abs(q30),abs(q31),abs(q32)))
      mm2  = (dmax1(abs(mm12),abs(mm22),abs(mm02),abs(mm32)))
      mpm2 = dmax1(mp2,mm2)
      mf   =  dmax1(abs(f1),abs(f2),abs(f3))
      norm = mpm2
      rm02 = abs(mm02)

      z(1,1) = 2d0*q10
      z(2,2) = 2d0*q20
      z(3,3) = 2d0*q30
      z(1,2) = q10+q20-q21
      z(1,3) = q10+q30-q31
      z(2,3) = q20+q30-q32

      zadj(1,1) = (4d0*q30*q20-k2k3*k2k3)
      zadj(1,2) = (k1k3*k2k3-2d0*q30*k1k2)
      zadj(1,3) = (k1k2*k2k3-2d0*q20*k1k3)
      zadj(2,2) = (4d0*q10*q30-k1k3*k1k3)
      zadj(2,3) = (k1k2*k1k3-2d0*q10*k2k3)
      zadj(3,3) = (4d0*q10*q20-k1k2*k1k2)
            
      mzadj =  dmax1(abs(zadj(1,1)),abs(zadj(2,2)),abs(zadj(3,3)),
     &               abs(zadj(1,2)),abs(zadj(1,3)),abs(zadj(2,3)))

      zadjf1 = zadj(1,1)*f1+zadj(1,2)*f2+zadj(1,3)*f3
      zadjf2 = zadj(1,2)*f1+zadj(2,2)*f2+zadj(2,3)*f3
      zadjf3 = zadj(1,3)*f1+zadj(2,3)*f2+zadj(3,3)*f3

c      write(*,*) 'cDp12345 ',adetz,zadjf1,zadjf2,zadjf3

      mzadjf  = dmax1(abs(zadjf1),abs(zadjf2),abs(zadjf3))
      mzadjfd = dmax1(mzadjf,adetz)

c --> changed 18.05.05
c>      xadj(1,1) = 2d0*mm02*zadj(1,1) - f2*f2*z(3,3)
c>     &         -f3*f3*z(2,2) + 2*f2*f3*z(2,3)
c>      xadj(2,2) = 2d0*mm02*zadj(2,2) - f1*f1*z(3,3)
c>     &         -f3*f3*z(1,1) + 2*f1*f3*z(1,3)
c>      xadj(3,3) = 2d0*mm02*zadj(3,3) - f1*f1*z(2,2)
c>     &         -f2*f2*z(1,1) + 2*f1*f2*z(1,2)
c>      xadj(1,2) = 2d0*mm02*zadj(1,2) + f2*f1*z(3,3)
c>     &         -f3*f1*z(2,3) - f2*f3*z(1,3)
c>     &         +f3*f3*z(1,2) 
c>      xadj(1,3) = 2d0*mm02*zadj(1,3) - f2*f1*z(2,3)
c>     &         +f3*f1*z(2,2) + f2*f2*z(1,3)
c>     &         -f3*f2*z(1,2) 
c>      xadj(2,3) = 2d0*mm02*zadj(2,3) + f1*f1*z(2,3)
c>     &         -f3*f1*z(1,2) - f1*f2*z(1,3)
c>     &         +f3*f2*z(1,1) 

      zadjadjff(1,1) = - f2*f2*z(3,3)
     &         -f3*f3*z(2,2) + 2*f2*f3*z(2,3)
      zadjadjff(2,2) = - f1*f1*z(3,3)
     &         -f3*f3*z(1,1) + 2*f1*f3*z(1,3)
      zadjadjff(3,3) = - f1*f1*z(2,2)
     &         -f2*f2*z(1,1) + 2*f1*f2*z(1,2)
      zadjadjff(1,2) = + f2*f1*z(3,3)
     &         -f3*f1*z(2,3) - f2*f3*z(1,3)
     &         +f3*f3*z(1,2) 
      zadjadjff(1,3) = - f2*f1*z(2,3)
     &         +f3*f1*z(2,2) + f2*f2*z(1,3)
     &         -f3*f2*z(1,2) 
      zadjadjff(2,3) = + f1*f1*z(2,3)
     &         -f3*f1*z(1,2) - f1*f2*z(1,3)
     &         +f3*f2*z(1,1) 

      mzadjadjff = dmax1(abs(zadjadjff(1,1)),abs(zadjadjff(1,2)),
     &    abs(zadjadjff(1,3)),abs(zadjadjff(2,2)),
     &    abs(zadjadjff(2,3)),abs(zadjadjff(3,3)))

      xadj(1,1) = 2d0*mm02*zadj(1,1) + zadjadjff(1,1)
      xadj(2,2) = 2d0*mm02*zadj(2,2) + zadjadjff(2,2)
      xadj(3,3) = 2d0*mm02*zadj(3,3) + zadjadjff(3,3)
      xadj(1,3) = 2d0*mm02*zadj(1,3) + zadjadjff(1,3)
      xadj(1,2) = 2d0*mm02*zadj(1,2) + zadjadjff(1,2)
      xadj(2,3) = 2d0*mm02*zadj(2,3) + zadjadjff(2,3)

c<--- changed 18.05.05

      mxadj = dmax1(abs(xadj(1,1)),abs(xadj(1,2)),abs(xadj(1,3)),
     &            abs(xadj(2,2)),abs(xadj(2,3)),abs(xadj(3,3)))
      adety = abs(2d0*mm02*detz  - zadjf1*f1 - zadjf2*f2 - zadjf3*f3)

c>      if (adetz.ne.0d0) then
c>        accpv = accC*(mpm2*mp2**2/adetz)**rankacc
c>      else
c>        accpv = accinf
c>      end if
      if (mzadjf.ne.0d0) then
c---> changed 18.05.05 
c        accg  = adetz/mzadjf*(mpm2**3/mzadjf)**rankacc
c     &      *(adetz*mpm2**3/mzadjf**2)**ordg4acc
c        accg  = adetz/mzadjf*(mf**2*mp2/mzadjf)**rankacc
c     &      *(adetz*mf**2*mp2/mzadjf**2)**ordg4acc
        h = dmax1(mzadjadjff/mzadjf,1d0)
        accg  = adetz/mzadjf*h**rankacc
     &      *(h*adetz/mzadjf)**ordg4acc
c<--- changed 18.05.05 
      else 
        accg = accinf
      end if
      if (mxadj*mp2.ne.0d0) then
        accgy = (mzadjfd*mpm2/(mxadj*mp2))**(ordgy4acc+1)
      else
        accgy = accinf
      end if
      if (mzadj*adety.ne.0d0) then
        accgy2 = (mzadjfd*dmax1(mxadj,mzadjf)/(mzadj*adety))
     &      **(ordgy4acc+1)
      else
        accgy2 = accinf
      end if
      if (mod(rankacc,2).eq.1) then
       if (adety.ne.0.and.adetz.ne.0) then
         accy  = accC*dmax1((mzadjf/adetz)**rankacc,
     &                 mzadjf/adetz*(mxadj/adetz)**((rankacc-1)/2),
     &                 norm*mzadjf/adety,
     &                 norm*mzadjf/adety*(mzadjf/adetz)**rankacc,
     &                 norm*mxadj/adety,
     &                 norm*mxadj/adety*(mzadjf/adetz)**(rankacc-1),
     &                 norm*mxadj/adety*(mxadj/adetz)**((rankacc-1)/2)
     &      )
c>         write(*,*) 'accy= ',
c>     &       accC*dmax1((mzadjf/adetz)**rankacc,
c>     &                 mzadjf/adetz*(mxadj/adetz)**((rankacc-1)/2),
c>     &                 norm*mzadjf/adety,
c>     &                 norm*mzadjf/adety*(mzadjf/adetz)**rankacc,
c>     &                 norm*mxadj/adety,
c>     &                 norm*mxadj/adety*(mzadjf/adetz)**(rankacc-1),
c>     &                norm*mxadj/adety*(mxadj/adetz)**((rankacc-1)/2)),
c>     &       accC*(mzadjf/adetz)**rankacc,
c>     &       accC*mzadjf/adetz*(mxadj/adetz)**((rankacc-1)/2),
c>     &       accC*norm*mzadjf/adety,
c>     &       accC*norm*mzadjf/adety*(mzadjf/adetz)**rankacc,
c>     &       accC*norm*mxadj/adety,
c>     &       accC*norm*mxadj/adety*(mzadjf/adetz)**(rankacc-1),
c>     &       accC*norm*mxadj/adety*(mxadj/adetz)**((rankacc-1)/2)
        else
          accy = accinf
        end if
        if (adetz.ne.0) then
          accpv =  accC*(norm*mzadj/adetz)
     &        *dmax1((mzadjf/adetz)**(rankacc-1),
     &        (rm02*mzadj/adetz)**((rankacc-1)/2))
        else 
          accpv = accinf
        end if
      else 
        accy  = accC*dmax1((mzadjf/adetz)**rankacc,
     &                 (mxadj/adetz)**(rankacc/2),
     &                 norm*mzadjf/adety,
     &                 norm*mzadjf/adety*(mzadjf/adetz)**rankacc,
     &                 norm*mxadj/adety,
     &                 norm*mxadj/adety*(mzadjf/adetz)**(rankacc-1),
     &                 norm*adetz/adety*(mxadj/adetz)**(rankacc/2),
     &                 norm*mzadjf/adety*(mxadj/adetz)**(rankacc/2)
     &      )
        if (rankacc.eq.0) then
          accpv = accC
        else if (adetz.ne.0) then
          accpv =  accC*(norm*mzadj/adetz)*(mzadjf/adetz)
     &        *dmax1((mzadjf/adetz)**(rankacc-2),
     &        (rm02*mzadj/adetz)**((rankacc-2)/2))
        else
          accpv = accinf
        end if
      end if
      if (mf.ne.0) then
        accgp = (mp2/mf)**(1+ordgp4acc)
      else
        accgp = accinf
      end if
c      accy = 0d0



c      accstop = 1d-2
c      acctest = 1d-4
c@@
c      if (.true.) then
      if (.false.) then
        write(*,*) 
        write(*,*) 'cDp12345: vers4 = ',version4
        write(*,*) 'cDp12345: rank4 = ',rank,ordg4,ordgy4,ordgp4
        write(*,*) 'cDp12345: ranka = ',rankacc,ordg4acc
     &      ,ordgy4acc,ordgp4acc
        write(*,*) 'cDp12345: mp2   = ',mp2,mpm2
        write(*,*) 'cDp12345: mzadj = ',norm,mzadj,mzadjadjff
        write(*,*) 'cDp12345: detz  = ',adetz,adetz/norm**3,
     &      (norm*mzadj)/adetz,(rm02*mzadj)/adetz,mzadjf/adetz
        write(*,*) 'cDp12345: mzadjf   = ',mzadjf,mzadjf/norm**3,
     &      mzadjadjff/mzadjf,adetz/mzadjf
        write(*,*) 'cDp12345: mxadj   = ',mxadj,mxadj/norm**3,
     &      mpm2/mp2,mzadjfd/mxadj
        write(*,*) 'cDp12345: dety  = ',adety,adety/norm**4,
     &      mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &      norm*mxadj/adety    
        write(*,*) 'cDp12345: mf    = ',mf,mf/norm,mp2/mf
        write(*,*) 'cDp12345: acc   = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
        write(*,*) 'cDp12345: p10   = ',p10,q10
        write(*,*) 'cDp12345: p21   = ',p21,q21
        write(*,*) 'cDp12345: p32   = ',p32,q32
        write(*,*) 'cDp12345: p30   = ',p30,q30
        write(*,*) 'cDp12345: p20   = ',p20,q20
        write(*,*) 'cDp12345: p31   = ',p31,q31
        write(*,*) 'cDp12345: m02   = ',m02,mm02
        write(*,*) 'cDp12345: m12   = ',m12,mm12
        write(*,*) 'cDp12345: m22   = ',m22,mm22
        write(*,*) 'cDp12345: m32   = ',m32,mm32
      end if

c      write(*,*) 'cDp12345 version',version4

      if (accgy2.le.accg.and.accgy2.le.accgy.and.accgy2.le.accy
     &      .and.accgy2.le.accgy2.and.accgy2.lt.accpv) then
        Dcount(6) = Dcount(6) + 1
        errwrite = .true.
        write(errout,*) 'cDp12345: accgy2 best!!!'
      end if

      if (version4.eq.0.or.version4.eq.5.or.version4.eq.6) then 

      if(accpv.lt.1d-12.or.
     &  accpv.le.accg.and.accpv.le.accgy.and.accpv.le.accy
     &       .and.accpv.le.accgp.and.accpv.lt.accstop) then
c        write(*,*) 'call Dpv'   
        call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank,switch)
        lc = 1
        Dcount(1) = Dcount(1) + 1
        if (accpv.gt.accworst) then
          errwrite = .true.
          accworst = accpv
        else
          errwrite = .false.       
        end if
        if (accpv.gt.accbad1) then
          Dcount(11) = Dcount(11) + 1
          if (accpv.gt.accbad2) then
            Dcount(21) = Dcount(21) + 1
            if (accpv.gt.accbad3) then
              Dcount(31) = Dcount(31) + 1
            end if  
          end if  
        end if  
      elseif (accg.le.accpv.and.accg.le.accgy.and.accg.le.accy
     &       .and.accg.le.accgp.and.accg.lt.accstop) then
c        write(*,*) 'call Dg'
        call cDg12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      D0,D1,D2,D3,D4,D5,rank)
        lc = 2
        Dcount(2) = Dcount(2) + 1
        if (accg.gt.accworst) then
          errwrite = .true.
          accworst = accg
        else
          errwrite = .false.       
        end if
        if (accg.gt.accbad1) then
          Dcount(12) = Dcount(12) + 1
          if (accg.gt.accbad2) then
            Dcount(22) = Dcount(22) + 1
            if (accg.gt.accbad3) then
              Dcount(32) = Dcount(32) + 1
            end if  
          end if  
        end if  
      elseif (accgy.le.accpv.and.accgy.le.accg.and.accgy.le.accy
     &       .and.accgy.le.accgp.and.accgy.lt.accstop) then
c            write(*,*) 'call Dgy'
        call cDgy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      D0,D1,D2,D3,D4,D5,rank)
        lc = 3
        Dcount(3) = Dcount(3) + 1
        if (accgy.gt.accworst) then
          errwrite = .true.
          accworst = accgy
        else
          errwrite = .false.       
        end if
        if (accgy.gt.accbad1) then
          Dcount(13) = Dcount(13) + 1
          if (accgy.gt.accbad2) then
            Dcount(23) = Dcount(23) + 1
            if (accgy.gt.accbad3) then
              Dcount(33) = Dcount(33) + 1
            end if  
          end if  
        end if  
      elseif (accy.le.accpv.and.accy.le.accg.and.accy.le.accgy
     &      .and.accy.le.accgp.and.accy.lt.accstop) then
c            write(*,*) 'call Dgy'
        call cDy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      D0,D1,D2,D3,D4,D5,rank,switch)
        lc = 4
        Dcount(4) = Dcount(4) + 1
        if (accy.gt.accworst) then
          errwrite = .true.
          accworst = accy
        else
          errwrite = .false.       
        end if
        if (accy.gt.accbad1) then
          Dcount(14) = Dcount(14) + 1
          if (accy.gt.accbad2) then
            Dcount(24) = Dcount(24) + 1
            if (accy.gt.accbad3) then
              Dcount(34) = Dcount(34) + 1
            end if  
          end if  
        end if  
      elseif (accgp.le.accpv.and.accgp.le.accg.and.accgp.le.accgy
     &      .and.accgp.le.accy.and.accgp.lt.accstop) then
c            write(*,*) 'call Dgy'
        call cDgp12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      D0,D1,D2,D3,D4,D5,rank)
        lc = 5
        Dcount(5) = Dcount(5) + 1
        if (accgp.gt.accworst) then
          errwrite = .true.
          accworst = accgp
        else
          errwrite = .false.       
        end if
        if (accgp.gt.accbad1) then
          Dcount(15) = Dcount(15) + 1
          if (accgp.gt.accbad2) then
            Dcount(25) = Dcount(25) + 1
            if (accgp.gt.accbad3) then
              Dcount(35) = Dcount(35) + 1
            end if  
          end if  
        end if  
      else
        if (count.le.20) then
          if (count.eq.0) then
            write(*,*) 
            write(*,*) 'cDp12345:  Dijk calculation unstable'
          end if

          write(errout,*)
          write(errout,*) 'cDp12345:  Dijk: calculation unstable'
          write(errout,*) 'cDp12345: vers4 = ',version4
          write(errout,*) 'cDp12345: rank4 = ',rank,ordg4,ordgy4,ordgp4
          write(errout,*) 'cDp12345: ranka = ',rankacc,ordg4acc
     &        ,ordgy4acc,ordgp4acc
          write(errout,*) 'cDp12345: mp2   = ',mp2,mpm2
          write(errout,*) 'cDp12345: mzadj = ',norm,mzadj,mzadjadjff
          write(errout,*) 'cDp12345: detz  = ',adetz,adetz/norm**3,
     &        (norm*mzadj)/adetz,(rm02*mzadj)/adetz,mzadjf/adetz
          write(errout,*) 'cDp12345: mzadjf   = ',mzadjf,mzadjf/norm**3,
     &        mzadjadjff/mzadjf,adetz/mzadjf
          write(errout,*) 'cDp12345: mxadj   = ',mxadj,mxadj/norm**3,
     &        mpm2/mp2,mzadjfd/mxadj
          write(errout,*) 'cDp12345: dety  = ',adety,adety/norm**4,
     &        mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &        norm*mxadj/adety    
          write(errout,*) 'cDp12345: mf    = ',mf,mf/norm,mp2/mf
          write(errout,*) 'cDp12345: acc   = ',
     &        accpv,accg,accgy,accy,accgp,accgy2

          write(errout,*) 'cDp12345:accwor = ',accworst
          write(errout,*) 'cDp12345:  p10  = ',p10,q10
          write(errout,*) 'cDp12345:  p21  = ',p21,q21
          write(errout,*) 'cDp12345:  p32  = ',p32,q32
          write(errout,*) 'cDp12345:  p30  = ',p30,q30
          write(errout,*) 'cDp12345:  p20  = ',p20,q20
          write(errout,*) 'cDp12345:  p31  = ',p31,q31
          write(errout,*) 'cDp12345:  m02  = ',m02,mm02
          write(errout,*) 'cDp12345:  m12  = ',m12,mm12
          write(errout,*) 'cDp12345:  m22  = ',m22,mm22
          write(errout,*) 'cDp12345:  m32  = ',m32,mm32
          count = count + 1
        end if
c        stop
        if (accpv.lt.accy) then
          call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &        D0,D1,D2,D3,D4,D5,rank,switch)
          lc = 8
          Dcount(8) = Dcount(8) + 1
          if (accpv.gt.accworst) then
            errwrite = .true.
            accworst = accpv
          else
            errwrite = .false.       
          end if
          if (accpv.gt.accbad1) then
            Dcount(18) = Dcount(18) + 1
            if (accpv.gt.accbad2) then
              Dcount(28) = Dcount(28) + 1
              if (accpv.gt.accbad3) then
                Dcount(38) = Dcount(38) + 1
              end if  
            end if  
          end if  
        else
          call cDy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &        D0,D1,D2,D3,D4,D5,rank,switch)
          lc = 9
          Dcount(9) = Dcount(9) + 1
          if (accy.gt.accworst) then
            errwrite = .true.
            accworst = accpv
          else
            errwrite = .false.       
          end if
          if (accy.gt.accbad1) then
            Dcount(19) = Dcount(19) + 1
            if (accy.gt.accbad2) then
              Dcount(29) = Dcount(29) + 1
              if (accy.gt.accbad3) then
                Dcount(39) = Dcount(39) + 1
              end if  
            end if  
          end if  
        end if
      end if      

c     write(*,*) 'errwrite = ',errwrite 
      if (errwrite) then
        if (flag.eq.0) then
          write(*,*)
          write(*,*) 'cDp12345:  Dijk calculation unstable'
          flag = 1
        end if

        write(errout,*)
        write(errout,*) 'cDp12345:  Dijk: calculation unstable'
        write(errout,*) 'cDp12345: vers4 = ',version4
        write(errout,*) 'cDp12345: rank4 = ',rank,ordg4,ordgy4,ordgp4
        write(errout,*) 'cDp12345: ranka = ',rankacc,ordg4acc
     &      ,ordgy4acc,ordgp4acc
        write(errout,*) 'cDp12345: mp2   = ',mp2,mpm2
        write(errout,*) 'cDp12345: mzadj = ',norm,mzadj,mzadjadjff
        write(errout,*) 'cDp12345: detz  = ',adetz,adetz/norm**3,
     &      (norm*mzadj)/adetz,(rm02*mzadj)/adetz,mzadjf/adetz
        write(errout,*) 'cDp12345: mzadjf   = ',mzadjf,mzadjf/norm**3,
     &      mzadjadjff/mzadjf,adetz/mzadjf
        write(errout,*) 'cDp12345: mxadj   = ',mxadj,mxadj/norm**3,
     &      mpm2/mp2,mzadjfd/mxadj
        write(errout,*) 'cDp12345: dety  = ',adety,adety/norm**4,
     &      mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &      norm*mxadj/adety    
        write(errout,*) 'cDp12345: mf    = ',mf,mf/norm,mp2/mf
        write(errout,*) 'cDp12345: acc   = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
        write(errout,*) 'cDp12345:accwor = ',accworst
        write(errout,*) 'cDp12345:  p10  = ',p10,q10
        write(errout,*) 'cDp12345:  p21  = ',p21,q21
        write(errout,*) 'cDp12345:  p32  = ',p32,q32
        write(errout,*) 'cDp12345:  p30  = ',p30,q30
        write(errout,*) 'cDp12345:  p20  = ',p20,q20
        write(errout,*) 'cDp12345:  p31  = ',p31,q31
        write(errout,*) 'cDp12345:  m02  = ',m02,mm02
        write(errout,*) 'cDp12345:  m12  = ',m12,mm12
        write(errout,*) 'cDp12345:  m22  = ',m22,mm22
        write(errout,*) 'cDp12345:  m32  = ',m32,mm32

      end if

      else if (version4.eq.4) then 

      if(accpv.lt.1d-12.or.
     &  accpv.le.accg.and.accpv.le.accgy.and.accpv.le.accy
     &      .and.accpv.lt.accstop) then
c        write(*,*) 'call Dpv'   
        call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank,switch)
        lc = 1
        Dcount(1) = Dcount(1) + 1
        if (accpv.gt.accworst) then
          errwrite = .true.
          accworst = accpv
        else
          errwrite = .false.       
        end if
        if (accpv.gt.accbad1) then
          Dcount(11) = Dcount(11) + 1
          if (accpv.gt.accbad2) then
            Dcount(21) = Dcount(21) + 1
            if (accpv.gt.accbad3) then
              Dcount(31) = Dcount(31) + 1
            end if  
          end if  
        end if  
      elseif (accg.le.accpv.and.accg.le.accgy.and.accg.le.accy
     &      .and.accg.lt.accstop) then
c        write(*,*) 'call Dg'
        call cDg12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      D0,D1,D2,D3,D4,D5,rank)
        lc = 2
        Dcount(2) = Dcount(2) + 1
        if (accg.gt.accworst) then
          errwrite = .true.
          accworst = accg
        else
          errwrite = .false.       
        end if
        if (accg.gt.accbad1) then
          Dcount(12) = Dcount(12) + 1
          if (accg.gt.accbad2) then
            Dcount(22) = Dcount(22) + 1
            if (accg.gt.accbad3) then
              Dcount(32) = Dcount(32) + 1
            end if  
          end if  
        end if  
      elseif (accgy.le.accpv.and.accgy.le.accg.and.accgy.le.accy
     &      .and.accgy.lt.accstop) then
c            write(*,*) 'call Dgy'
        call cDgy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      D0,D1,D2,D3,D4,D5,rank)
        lc = 3
        Dcount(3) = Dcount(3) + 1
        if (accgy.gt.accworst) then
          errwrite = .true.
          accworst = accgy
        else
          errwrite = .false.       
        end if
        if (accgy.gt.accbad1) then
          Dcount(13) = Dcount(13) + 1
          if (accgy.gt.accbad2) then
            Dcount(23) = Dcount(23) + 1
            if (accgy.gt.accbad3) then
              Dcount(33) = Dcount(33) + 1
            end if  
          end if  
        end if  
      elseif (accy.le.accpv.and.accy.le.accg.and.accy.le.accgy
     &      .and.accy.lt.accstop) then
c            write(*,*) 'call Dgy'
        call cDy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      D0,D1,D2,D3,D4,D5,rank,switch)
        lc = 4
        Dcount(4) = Dcount(4) + 1
        if (accy.gt.accworst) then
          errwrite = .true.
          accworst = accy
        else
          errwrite = .false.       
        end if
        if (accy.gt.accbad1) then
          Dcount(14) = Dcount(14) + 1
          if (accy.gt.accbad2) then
            Dcount(24) = Dcount(24) + 1
            if (accy.gt.accbad3) then
              Dcount(34) = Dcount(34) + 1
            end if  
          end if  
        end if  
      else
        if (count.le.20) then
          if (count.eq.0) then
            write(*,*) 
            write(*,*) 'cDp12345:  Dijk calculation unstable'
          end if

          write(errout,*)
          write(errout,*) 'cDp12345:  Dijk: calculation unstable'
          write(errout,*) 'cDp12345: vers4 = ',version4
          write(errout,*) 'cDp12345: rank4 = ',rank,ordg4,ordgy4,ordgp4
          write(errout,*) 'cDp12345: ranka = ',rankacc,ordg4acc
     &        ,ordgy4acc,ordgp4acc
          write(errout,*) 'cDp12345: mp2   = ',mp2,mpm2
          write(errout,*) 'cDp12345: mzadj = ',norm,mzadj,mzadjadjff
          write(errout,*) 'cDp12345: detz  = ',adetz,adetz/norm**3,
     &        (norm*mzadj)/adetz,(rm02*mzadj)/adetz,mzadjf/adetz
          write(errout,*) 'cDp12345: mzadjf   = ',mzadjf,mzadjf/norm**3,
     &        mzadjadjff/mzadjf,adetz/mzadjf
          write(errout,*) 'cDp12345: mxadj   = ',mxadj,mxadj/norm**3,
     &        mpm2/mp2,mzadjfd/mxadj
          write(errout,*) 'cDp12345: dety  = ',adety,adety/norm**4,
     &        mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &        norm*mxadj/adety    
          write(errout,*) 'cDp12345: mf    = ',mf,mf/norm,mp2/mf
          write(errout,*) 'cDp12345: acc   = ',
     &        accpv,accg,accgy,accy,accgp,accgy2
          write(errout,*) 'cDp12345:accwor = ',accworst
          write(errout,*) 'cDp12345:  p10  = ',p10,q10
          write(errout,*) 'cDp12345:  p21  = ',p21,q21
          write(errout,*) 'cDp12345:  p32  = ',p32,q32
          write(errout,*) 'cDp12345:  p30  = ',p30,q30
          write(errout,*) 'cDp12345:  p20  = ',p20,q20
          write(errout,*) 'cDp12345:  p31  = ',p31,q31
          write(errout,*) 'cDp12345:  m02  = ',m02,mm02
          write(errout,*) 'cDp12345:  m12  = ',m12,mm12
          write(errout,*) 'cDp12345:  m22  = ',m22,mm22
          write(errout,*) 'cDp12345:  m32  = ',m32,mm32
          count = count + 1
        end if
c        stop
        if (accpv.lt.accy) then
          call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &        D0,D1,D2,D3,D4,D5,rank,switch)
          lc = 8
          Dcount(8) = Dcount(8) + 1
          if (accpv.gt.accworst) then
            errwrite = .true.
            accworst = accpv
          else
            errwrite = .false.       
          end if
          if (accpv.gt.accbad1) then
            Dcount(18) = Dcount(18) + 1
            if (accpv.gt.accbad2) then
              Dcount(28) = Dcount(28) + 1
              if (accpv.gt.accbad3) then
                Dcount(38) = Dcount(38) + 1
              end if  
            end if  
          end if  
        else
          call cDy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &        D0,D1,D2,D3,D4,D5,rank,switch)
          lc = 9
          Dcount(9) = Dcount(9) + 1
          if (accy.gt.accworst) then
            errwrite = .true.
            accworst = accpv
          else
            errwrite = .false.       
          end if
          if (accy.gt.accbad1) then
            Dcount(19) = Dcount(19) + 1
            if (accy.gt.accbad2) then
              Dcount(29) = Dcount(29) + 1
              if (accy.gt.accbad3) then
                Dcount(39) = Dcount(39) + 1
              end if  
            end if  
          end if  
        end if
      end if      

c     write(*,*) 'errwrite = ',errwrite 
      if (errwrite) then
        if (flag.eq.0) then
          write(*,*)
          write(*,*) 'cDp12345:  Dijk calculation unstable'
          flag = 1
        end if

        write(errout,*)
        write(errout,*) 'cDp12345:  Dijk: calculation unstable'
        write(errout,*) 'cDp12345: vers4 = ',version4
        write(errout,*) 'cDp12345: rank4 = ',rank,ordg4,ordgy4,ordgp4
        write(errout,*) 'cDp12345: ranka = ',rankacc,ordg4acc
     &      ,ordgy4acc,ordgp4acc
        write(errout,*) 'cDp12345: mp2   = ',mp2,mpm2
        write(errout,*) 'cDp12345: mzadj = ',norm,mzadj,mzadjadjff
        write(errout,*) 'cDp12345: detz  = ',adetz,adetz/norm**3,
     &      (norm*mzadj)/adetz,(rm02*mzadj)/adetz,mzadjf/adetz
        write(errout,*) 'cDp12345: mzadjf   = ',mzadjf,mzadjf/norm**3,
     &      mzadjadjff/mzadjf,adetz/mzadjf
        write(errout,*) 'cDp12345: mxadj   = ',mxadj,mxadj/norm**3,
     &      mpm2/mp2,mzadjfd/mxadj
        write(errout,*) 'cDp12345: dety  = ',adety,adety/norm**4,
     &      mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &      norm*mxadj/adety    
        write(errout,*) 'cDp12345: mf    = ',mf,mf/norm,mp2/mf
        write(errout,*) 'cDp12345: acc   = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
        write(errout,*) 'cDp12345:accwor = ',accworst
        write(errout,*) 'cDp12345:  p10  = ',p10,q10
        write(errout,*) 'cDp12345:  p21  = ',p21,q21
        write(errout,*) 'cDp12345:  p32  = ',p32,q32
        write(errout,*) 'cDp12345:  p30  = ',p30,q30
        write(errout,*) 'cDp12345:  p20  = ',p20,q20
        write(errout,*) 'cDp12345:  p31  = ',p31,q31
        write(errout,*) 'cDp12345:  m02  = ',m02,mm02
        write(errout,*) 'cDp12345:  m12  = ',m12,mm12
        write(errout,*) 'cDp12345:  m22  = ',m22,mm22
        write(errout,*) 'cDp12345:  m32  = ',m32,mm32

      end if

      else if (version4.eq.1) then

        call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank,switch)
        lc = 1
        Dcount(1) = Dcount(1) + 1

      else if (version4.eq.2) then 

      if(accpv.lt.1d-12.or.
     &  accpv.lt.accg.and.accpv.lt.accstop) then
c        write(*,*) 'call Dpv'   
        call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank,switch)
        lc = 1
        Dcount(1) = Dcount(1) + 1
      elseif (accg.lt.accpv.and.accg.lt.accstop) then
c        write(*,*) 'call Dg 2'
        call cDg12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      D0,D1,D2,D3,D4,D5,rank)
        lc = 2
        Dcount(2) = Dcount(2) + 1
      else
        if (count.le.20) then
          if (count.eq.0) then
            write(*,*) 
            write(*,*) 'cDp12345:  Dijk calculation unstable'
          end if

          write(errout,*) 
          write(errout,*) 'cDp12345:  Dijk: calculation unstable'
          write(errout,*) 'cDp12345: vers4 = ',version4
          write(errout,*) 'cDp12345: rank4 = ',rank,ordg4,ordgy4,ordgp4
          write(errout,*) 'cDp12345: ranka = ',rankacc,ordg4acc
     &        ,ordgy4acc,ordgp4acc
          write(errout,*) 'cDp12345: mp2   = ',mp2,mpm2
          write(errout,*) 'cDp12345: mzadj = ',norm,mzadj,mzadjadjff
          write(errout,*) 'cDp12345: detz  = ',adetz,adetz/norm**3,
     &        (norm*mzadj)/adetz,(rm02*mzadj)/adetz,mzadjf/adetz
          write(errout,*) 'cDp12345: mzadjf   = ',mzadjf,mzadjf/norm**3,
     &        mzadjadjff/mzadjf,adetz/mzadjf
          write(errout,*) 'cDp12345: mxadj   = ',mxadj,mxadj/norm**3,
     &        mpm2/mp2,mzadjfd/mxadj
          write(errout,*) 'cDp12345: dety  = ',adety,adety/norm**4,
     &        mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &        norm*mxadj/adety    
          write(errout,*) 'cDp12345: mf    = ',mf,mf/norm,mp2/mf
          write(errout,*) 'cDp12345: acc   = ',
     &        accpv,accg,accgy,accy,accgp,accgy2
          write(errout,*) 'cDp12345:  p10  = ',p10
          write(errout,*) 'cDp12345:  p21  = ',p21
          write(errout,*) 'cDp12345:  p32  = ',p32
          write(errout,*) 'cDp12345:  p30  = ',p30
          write(errout,*) 'cDp12345:  p20  = ',p20
          write(errout,*) 'cDp12345:  p31  = ',p31
          write(errout,*) 'cDp12345:  m02  = ',m02
          write(errout,*) 'cDp12345:  m12  = ',m12
          write(errout,*) 'cDp12345:  m22  = ',m22
          write(errout,*) 'cDp12345:  m32  = ',m32
          count = count + 1
        end if
        call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank,switch)
c        stop
      end if      

      else if (version4.eq.3) then 

      if(accpv.lt.1d-12.or.
     &  accpv.lt.accg.and.accpv.lt.accgy.and.accpv.lt.accstop) then
c        write(*,*) 'call Dpv'   
        call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank,switch)
        lc = 1
      elseif (accg.lt.accpv.and.accg.lt.accgy.and.accg.lt.accstop) then
c        write(*,*) 'call Dg 2'
        call cDg12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      D0,D1,D2,D3,D4,D5,rank)
        lc = 2
      elseif (accgy.lt.accpv.and.accgy.lt.accg.and.accgy.lt.accstop)then
c            write(*,*) 'call Dgy'
        call cDgy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      D0,D1,D2,D3,D4,D5,rank)
        lc = 3
      else
        if (count.le.20) then
          if (count.eq.0) then
            write(*,*) 'cDp12345:  Dijk calculation unstable'
          end if

          write(errout,*)
          write(errout,*) 'cDp12345:  Dijk: calculation unstable'
          write(errout,*) 'cDp12345: vers4 = ',version4
          write(errout,*) 'cDp12345: rank4 = ',rank,ordg4,ordgy4,ordgp4
          write(errout,*) 'cDp12345: ranka = ',rankacc,ordg4acc
     &        ,ordgy4acc,ordgp4acc
          write(errout,*) 'cDp12345: mp2   = ',mp2,mpm2
          write(errout,*) 'cDp12345: mzadj = ',norm,mzadj,mzadjadjff
          write(errout,*) 'cDp12345: detz  = ',adetz,adetz/norm**3,
     &        (norm*mzadj)/adetz,(rm02*mzadj)/adetz,mzadjf/adetz
          write(errout,*) 'cDp12345: mzadjf   = ',mzadjf,mzadjf/norm**3,
     &        mzadjadjff/mzadjf,adetz/mzadjf
          write(errout,*) 'cDp12345: mxadj   = ',mxadj,mxadj/norm**3,
     &        mpm2/mp2,mzadjfd/mxadj
          write(errout,*) 'cDp12345: dety  = ',adety,adety/norm**4,
     &        mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &        norm*mxadj/adety    
          write(errout,*) 'cDp12345: mf    = ',mf,mf/norm,mp2/mf
          write(errout,*) 'cDp12345: acc   = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
          write(errout,*) 'cDp12345:  p10  = ',p10,q10
          write(errout,*) 'cDp12345:  p21  = ',p21,q21
          write(errout,*) 'cDp12345:  p32  = ',p32,q32
          write(errout,*) 'cDp12345:  p30  = ',p30,q30
          write(errout,*) 'cDp12345:  p20  = ',p20,q20
          write(errout,*) 'cDp12345:  p31  = ',p31,q31
          write(errout,*) 'cDp12345:  m02  = ',m02,mm02
          write(errout,*) 'cDp12345:  m12  = ',m12,mm12
          write(errout,*) 'cDp12345:  m22  = ',m22,mm22
          write(errout,*) 'cDp12345:  m32  = ',m32,mm32
          count = count + 1
        end if
        call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank,switch)
c        stop
      end if      

      else

        write(*,*) 'cDp12345:  version4 not defined'        
        write(*,*) 'cDp12345:  version4 = ',version4        

      end if

      if (Dcount(1).lt.2) then
c        write(*,*) 'cDp12345 lc = ',lc
      end if


c write to cash

      fct(1)=D0
      
      if (rank.eq.0) goto 521
      
      fct(2)=D1(1)
      fct(3)=D1(2)
      fct(4)=D1(3)
      
      if (rank.eq.1.and.switch.eq.0) goto 521
      
      fct(5)=D2(0,0)
      
      if (rank.eq.1) goto 521
      
      fct(6) =D2(1,1)
      fct(7) =D2(1,2)
      fct(8) =D2(1,3)
      fct(9) =D2(2,2)
      fct(10)=D2(2,3)
      fct(11)=D2(3,3)
      
      if (rank.eq.2.and.switch.eq.0) goto 521
      
      fct(12)=D3(0,0,1)
      fct(13)=D3(0,0,2)
      fct(14)=D3(0,0,3)
      
      if (rank.eq.2) goto 521
      
      fct(15)=D3(1,1,1)
      fct(16)=D3(1,1,2)
      fct(17)=D3(1,1,3)
      fct(18)=D3(1,2,2)
      fct(19)=D3(1,2,3)
      fct(20)=D3(1,3,3)
      fct(21)=D3(2,2,2)
      fct(22)=D3(2,2,3)
      fct(23)=D3(2,3,3)
      fct(24)=D3(3,3,3)
      
      if (rank.eq.3.and.switch.eq.0) goto 521
      
      fct(25)=D4(0,0,0,0)
      fct(26)=D4(0,0,1,1)
      fct(27)=D4(0,0,1,2)
      fct(28)=D4(0,0,1,3)
      fct(29)=D4(0,0,2,2)
      fct(30)=D4(0,0,2,3)
      fct(31)=D4(0,0,3,3)
      
      if (rank.eq.3) goto 521
      
      fct(32)=D4(1,1,1,1)
      fct(33)=D4(1,1,1,2)
      fct(34)=D4(1,1,1,3)
      fct(35)=D4(1,1,2,2)
      fct(36)=D4(1,1,2,3)
      fct(37)=D4(1,1,3,3)
      fct(38)=D4(1,2,2,2)
      fct(39)=D4(1,2,2,3)
      fct(40)=D4(1,2,3,3)
      fct(41)=D4(1,3,3,3)
      fct(42)=D4(2,2,2,2)
      fct(43)=D4(2,2,2,3)
      fct(44)=D4(2,2,3,3)
      fct(45)=D4(2,3,3,3)
      fct(46)=D4(3,3,3,3)
      
      if (rank.eq.4.and.switch.eq.0) goto 521
      
      fct(47)=D5(0,0,0,0,1)
      fct(48)=D5(0,0,0,0,2)
      fct(49)=D5(0,0,0,0,3)
      fct(50)=D5(0,0,1,1,1)
      fct(51)=D5(0,0,1,1,2)
      fct(52)=D5(0,0,1,1,3)
      fct(53)=D5(0,0,1,2,2)
      fct(54)=D5(0,0,1,2,3)
      fct(55)=D5(0,0,1,3,3)
      fct(56)=D5(0,0,2,2,2)
      fct(57)=D5(0,0,2,2,3)
      fct(58)=D5(0,0,2,3,3)
      fct(59)=D5(0,0,3,3,3)
      
      if (rank.eq.4) goto 521
      
      fct(60)=D5(1,1,1,1,1)
      fct(61)=D5(1,1,1,1,2)
      fct(62)=D5(1,1,1,1,3)
      fct(63)=D5(1,1,1,2,2)
      fct(64)=D5(1,1,1,2,3)
      fct(65)=D5(1,1,1,3,3)
      fct(66)=D5(1,1,2,2,2)
      fct(67)=D5(1,1,2,2,3)
      fct(68)=D5(1,1,2,3,3)
      fct(69)=D5(1,1,3,3,3)
      fct(70)=D5(1,2,2,2,2)
      fct(71)=D5(1,2,2,2,3)
      fct(72)=D5(1,2,2,3,3)
      fct(73)=D5(1,2,3,3,3)
      fct(74)=D5(1,3,3,3,3)
      fct(75)=D5(2,2,2,2,2)
      fct(76)=D5(2,2,2,2,3)
      fct(77)=D5(2,2,2,3,3)
      fct(78)=D5(2,2,3,3,3)
      fct(79)=D5(2,3,3,3,3)
      fct(80)=D5(3,3,3,3,3)

 521  continue

    
      call cachewrite(fct,80,type)

c>      if (n.eq.4.and.usecache.eq.1)  then
c>        write(85,*) 'nto = ',ncall(2),countD,sym
c>        write(85,*) 'D0 = ',D0
c>        write(85,*) 'D1 = ',D1
c>        write(85,*) 'D2 = ',D2
c>        write(85,*) 'D3 = ',D3
c>c        write(85,*) 'D4 = ',D4(0,0,0,0)
c>c        write(85,*) 'D4 = ',D4(0,0,1,1)
c>        do i1=1,3
c>        do i2=3,3
c>        do i3=2,2
c>c        write(85,133) i1,i2,i3,D3(i1,i2,i3)
c>        end do
c>        end do
c>        end do
c>      else if (n.eq.4.and.usecache.eq.0) then
c>        write(86,*) 'nto = ',countD,countD,sym
c>        write(86,*) 'D0 = ',D0
c>        write(86,*) 'D1 = ',D1
c>        write(86,*) 'D2 = ',D2
c>        write(86,*) 'D3 = ',D3
c>c        write(86,*) 'D4 = ',D4(0,0,0,0)
c>c        write(86,*) 'D4 = ',D4(0,0,1,1)
c>        do i1=1,3
c>        do i2=3,3
c>        do i3=2,2
c>c        write(86,133) i1,i2,i3,D3(i1,i2,i3)
c>        end do
c>        end do
c>        end do
c>      end if

c test part

c      ltest = 20

      
      if(ltest.gt.1)then
      atest = 0
      call cachetempoff
      if (ltest.gt.1.and.counttest.lt.100) then
        counttest = counttest + 1
        if (lc.eq.1.and.accpv.gt.acctest) then
c          if (accg.lt.accgy) then
            call cDg12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &          Dt0,Dt1,Dt2,Dt3,Dt4,Dt5,rank)
c         else
c          end if
          atest = 1 
        else if (lc.eq.2.and.accg.gt.acctest) then
c          if (accgy.lt.accpv) then
            call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &          Dt0,Dt1,Dt2,Dt3,Dt4,Dt5,rank,switch)
c          else
c          end if
          atest = 1 
c        else if (lc.eq.3.and.accgy.gt.acctest) then
c          if (accg.lt.accpv) then
c            call cDg12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
c     &          D0,D1,D2,D3,D4,D5,rank)
c          else
c            call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
c     &          D0,D1,D2,D3,D4,D5,rank,switch)
c          end if
c          atest = 1 
        end if        
      end if

      if (ltest.eq.20) then
      if(accpv.lt.1d-12.or.
     &  accpv.lt.accg.and.accpv.lt.accgy.and.accpv.lt.accy
     &      .and.accpv.lt.accstop) then
c        write(*,*) 'call Dpv'   
        call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    Dt0,Dt1,Dt2,Dt3,Dt4,Dt5,rank,switch)
        lc = 1
      elseif (accg.lt.accpv.and.accg.lt.accgy.and.accg.lt.accy
     &      .and.accg.lt.accstop) then
c        write(*,*) 'call Dg'
        call cDg12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      Dt0,Dt1,Dt2,Dt3,Dt4,Dt5,rank)
        lc = 2
      elseif (accgy.lt.accpv.and.accgy.lt.accg.and.accgy.lt.accy
     &      .and.accgy.lt.accstop) then
c            write(*,*) 'call Dgy'
        call cDgy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      Dt0,Dt1,Dt2,Dt3,Dt4,Dt5,rank)
        lc = 3
      elseif (accy.lt.accpv.and.accy.lt.accg.and.accy.lt.accgy
     &      .and.accy.lt.accstop) then
c            write(*,*) 'call Dgy'
        call cDy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &      Dt0,Dt1,Dt2,Dt3,Dt4,Dt5,rank,switch)
        lc = 4
      else
        if (accpv.lt.accy) then
          call cDpv12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &        Dt0,Dt1,Dt2,Dt3,Dt4,Dt5,rank,switch)
          lc = 5
        else
          call cDy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &        Dt0,Dt1,Dt2,Dt3,Dt4,Dt5,rank,switch)
          lc = 6
        end if
      end if      
      atest = 1
      end if

c        write(*,*) 'atest ',atest,ltest
          
        if (atest.eq.1) then

 99     format(' deviation     = ',G20.14)
 100    format(' D0         = ',G20.14,' + i* ',G20.14)
 200    format(' D0         = ',G20.14,' + i* ',G20.14)
 220    format(' D2(0,0)    = ',G20.14,' + i* ',G20.14)
 120    format(' D2(0,0)    = ',G20.14,' + i* ',G20.14)
 140    format(' D4(0,0,0,0)= ',G20.14,' + i* ',G20.14)
 240    format(' D4(0,0,0,0)= ',G20.14,' + i* ',G20.14)
 111    format(' D1(',i1,')       = ',G20.14,' + i* ',G20.14)
 211    format(' D1(',i1,')       = ',G20.14,' + i* ',G20.14)
 131    format(' D3(0,0,',i1,')   = ',G20.14,' + i* ',G20.14)
 231    format(' D3(0,0,',i1,')   = ',G20.14,' + i* ',G20.14)
 151    format(' D5(0,0,0,0,',i1,')   = ',G20.14,' + i* ',G20.14)
 251    format(' D5(0,0,0,0,',i1,')   = ',G20.14,' + i* ',G20.14)
 142    format(' D4(0,0,',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 242    format(' D4(0,0,',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 133    format(' D3(',i1,',',i1,',',i1,')   = ',G20.14,' + i* ',G20.14)
 233    format(' D3(',i1,',',i1,',',i1,')   = ',G20.14,' + i* ',G20.14)
 153    format(' D5(0,0,',i1,',',i1,',',i1,') = ',
     &      G20.14,' + i* ',G20.14)
 253    format(' D5(0,0,',i1,',',i1,',',i1,') = ',
     &      G20.14,' + i* ',G20.14)
 144    format(' D4(',i1,',',i1,',',i1,',',i1,') = ',G20.14,
     &            ' + i* ',G20.14)
 244    format(' D4(',i1,',',i1,',',i1,',',i1,') = ',G20.14,
     &            ' + i* ',G20.14)
 155    format(' D5(',i1,',',i1,',',i1,',',i1,',',i1,') = ',G20.14,
     &            ' + i* ',G20.14)
 255    format(' D5(',i1,',',i1,',',i1,',',i1,',',i1,') = ',G20.14,
     &            ' + i* ',G20.14)

        flag1 = 0
        if(cdabs(Dt0/D0-1D0).gt.acc) then
          write(testout,100) D0
          write(testout,200) Dt0
          write(testout,99) cdabs(Dt0/D0-1D0)
          flag1 = 1          
        end if
        if (rank.gt.1) then
          if(cdabs(Dt2(0,0)/D2(0,0)-1D0).gt.acc) then
            write(testout,120) D2(0,0)
            write(testout,220) Dt2(0,0)
            write(testout,99) cdabs(Dt2(0,0)/D2(0,0)-1D0)
            flag1 = 1          
          end if
          if (rank.gt.3) then      
            if(cdabs(Dt4(0,0,0,0)/D4(0,0,0,0)-1D0).gt.acc) then
              write(testout,140) D4(0,0,0,0)
              write(testout,240) Dt4(0,0,0,0)
              write(testout,99) cdabs(Dt4(0,0,0,0)/D4(0,0,0,0)-1D0)
              flag1 = 1          
            end if
          end if 
        end if 
        if (rank.gt.0) then
          do 101 i1=1,3
            if(cdabs(Dt1(i1)/D1(i1)-1D0).gt.acc) then
              write(testout,111) i1,D1(i1)
              write(testout,211) i1,Dt1(i1)
              write(testout,99) cdabs(Dt1(i1)/D1(i1)-1D0)
              flag1 = 1          
            end if
            if (rank.gt.2) then      
              if(cdabs(Dt3(0,0,i1)/D3(0,0,i1)-1D0).gt.acc) then
                write(testout,131) i1,D3(0,0,i1)
                write(testout,231) i1,Dt3(0,0,i1)
                write(testout,99) cdabs(Dt3(0,0,i1)/D3(0,0,i1)-1D0)
                flag1 = 1          
              end if
              if (rank.gt.4) then      
                if(cdabs(Dt5(0,0,0,0,i1)/D5(0,0,0,0,i1)-1D0).gt.acc)then
                  write(testout,151) i1,D5(0,0,0,0,i1)
                  write(testout,251) i1,Dt5(0,0,0,0,i1)
                  write(testout,99) 
     &                cdabs(Dt5(0,0,0,0,i1)/D5(0,0,0,0,i1)-1D0)
                  flag1 = 1          
                end if
              end if
            end if
            if (rank.gt.1) then
              do 102 i2=i1,3
                if(cdabs(Dt2(i1,i2)/D2(i1,i2)-1D0).gt.acc) then
                  write(testout,142) i1,i2,D2(i1,i2)
                  write(testout,242) i1,i2,Dt2(i1,i2)
                  write(testout,99) cdabs(Dt2(i1,i2)/D2(i1,i2)-1D0)
                  flag1 = 1          
                end if
                if (rank.gt.3) then
                  if(cdabs(Dt4(0,0,i1,i2)/D4(0,0,i1,i2)-1D0).gt.acc)then
                    write(testout,142) i1,i2,D4(0,0,i1,i2)
                    write(testout,242) i1,i2,Dt4(0,0,i1,i2)
                    write(testout,99) 
     &                  cdabs(Dt4(0,0,i1,i2)/D4(0,0,i1,i2)-1D0)
                    flag1 = 1          
                  end if
                end if
                if (rank.gt.2) then
                  do 103 i3=i2,3
                    if(cdabs(Dt3(i1,i2,i3)/D3(i1,i2,i3)-1D0).gt.acc)then
                      write(testout,133) i1,i2,i3,D3(i1,i2,i3)
                      write(testout,233) i1,i2,i3,Dt3(i1,i2,i3)
                      write(testout,99) 
     &                    cdabs(Dt3(i1,i2,i3)/D3(i1,i2,i3)-1D0)
                      flag1 = 1          
                    end if
                    if (rank.gt.4) then
                      if(cdabs(Dt5(0,0,i1,i2,i3)/D5(0,0,i1,i2,i3)-1D0)
     &                 .gt.acc) then
                        write(testout,153) i1,i2,i3,D5(0,0,i1,i2,i3)
                        write(testout,253) i1,i2,i3,Dt5(0,0,i1,i2,i3)
                        write(testout,99) 
     &                    cdabs(Dt5(0,0,i1,i2,i3)/D5(0,0,i1,i2,i3)-1D0)
                      flag1 = 1          
                      end if
                    end if
                    if (rank.gt.3) then
                      do 104 i4=i3,3
                        if(cdabs(Dt4(i1,i2,i3,i4)/D4(i1,i2,i3,i4)
     &                   -1D0).gt.acc)then
                          write(testout,144) i1,i2,i3,i4,D4(i1,i2,i3,i4)
                          write(testout,244)i1,i2,i3,i4,Dt4(i1,i2,i3,i4)
                          write(testout,99)
     &                        cdabs(Dt4(i1,i2,i3,i4)/D4(i1,i2,i3,i4)
     &                        -1D0)
                          flag1 = 1          
                        end if
                        if (rank.gt.4) then
                          do 105 i5=i4,3
                            if(cdabs(Dt5(i1,i2,i3,i4,i5)
     &                          /D5(i1,i2,i3,i4,i5)-1D0).gt.acc) then
                              write(testout,144) i1,i2,i3,i4,i5,
     &                            D5(i1,i2,i3,i4,i5)
                              write(testout,244) i1,i2,i3,i4,i5,
     &                            Dt5(i1,i2,i3,i4,i5)
                              write(testout,99)
     &                            cdabs(Dt5(i1,i2,i3,i4,i5)
     &                            /D5(i1,i2,i3,i4,i5)-1D0)
                              flag1 = 1          
                            end if
 105                      continue
                        end if
 104                  continue
                    end if
 103              continue
                end if
 102          continue
            end if
 101      continue
        end if
c       write(*,*) 'cDp flag1= ',flag1,nocalc
  
        if (flag1.eq.1) then
        write(testout,*) 'cDp12345:  p10  = ',p10,q10
        write(testout,*) 'cDp12345:  p21  = ',p21,q21
        write(testout,*) 'cDp12345:  p32  = ',p32,q32
        write(testout,*) 'cDp12345:  p30  = ',p30,q30
        write(testout,*) 'cDp12345:  p20  = ',p20,q20
        write(testout,*) 'cDp12345:  p31  = ',p31,q31
        write(testout,*) 'cDp12345:  m02  = ',m02,mm02
        write(testout,*) 'cDp12345:  m12  = ',m12,mm12
        write(testout,*) 'cDp12345:  m22  = ',m22,mm22
        write(testout,*) 'cDp12345:  m32  = ',m32,mm32
        write(testout,*) 'cDp12345: vers4 = ',version4
        write(testout,*) 'cDp12345: rank4 = ',rank,ordg4,ordgy4,ordgp4
        write(testout,*) 'cDp12345: ranka = ',rankacc,ordg4acc
     &      ,ordgy4acc,ordgp4acc
        write(testout,*) 'cDp12345: mp2   = ',mp2,mpm2
        write(testout,*) 'cDp12345: mzadj = ',norm,mzadj,mzadjadjff
        write(testout,*) 'cDp12345: detz  = ',adetz,adetz/norm**3,
     &      (norm*mzadj)/adetz,(rm02*mzadj)/adetz,mzadjf/adetz
        write(testout,*) 'cDp12345: mzadjf   = ',mzadjf,mzadjf/norm**3,
     &      mzadjadjff/mzadjf,adetz/mzadjf
        write(testout,*) 'cDp12345: mxadj   = ',mxadj,mxadj/norm**3,
     &      mpm2/mp2,mzadjfd/mxadj
        write(testout,*) 'cDp12345: dety  = ',adety,adety/norm**4,
     &      mzadjf/adetz,mxadj/adetz,norm*mzadjf/adety,
     &      norm*mxadj/adety    
        write(testout,*) 'cDp12345: mf    = ',mf,mf/norm,mp2/mf
        write(testout,*) 'cDp12345: acc   = ',
     &      accpv,accg,accgy,accy,accgp,accgy2
        write(testout,*)
        end if

      end if
      call cachereon
      end if

 522  continue

      end
                        

***********************************************************************
      subroutine cDy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank,switch)
***********************************************************************
*     4-point tensor coefficient functions D0,D1,D2,D3,D4,D5          *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     alternative reduction via Cayley matrix                         *
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
*     i-th denominator cancelled: Cxwi                                *
***********************************************************************
*     11.11.04 Ansgar Denner    last changed  24.01.06                *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p32,p30,p20,p31
      complex*16 q10,q21,q32,q30,q20,q31
      complex*16 m02,m12,m22,m32
      complex*16 mm02,mm12,mm22,mm32
      complex*16 mx(0:3,0:3),mxinv(0:3,0:3)
c      complex*16 test(0:3,0:3)
      complex*16 C0w3,C1w3(2),C2w3(0:2,0:2),C3w3(0:2,0:2,0:2)
     &    ,C4w3(0:2,0:2,0:2,0:2),C5w3(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w3(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w2,C1w2(2),C2w2(0:2,0:2),C3w2(0:2,0:2,0:2)
     &    ,C4w2(0:2,0:2,0:2,0:2),C5w2(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w2(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w1,C1w1(2),C2w1(0:2,0:2),C3w1(0:2,0:2,0:2)
     &    ,C4w1(0:2,0:2,0:2,0:2),C5w1(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w1(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w0,C1w0(2),C2w0(0:2,0:2),C3w0(0:2,0:2,0:2)
     &    ,C4w0(0:2,0:2,0:2,0:2),C5w0(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w0(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0m(0:3),C1m(0:3,3),C2m(0:3,0:3,0:3)
     &    ,C3m(0:3,0:3,0:3,0:3),C4m(0:3,0:3,0:3,0:3,0:3)
c     &    ,C5m(0:3,0:3,0:3,0:3,0:3,0:3),C6m(0:3,0:3,0:3,0:3,0:3,0:3,0:3)
      complex*16 S1hat(3),S2hat(3,3),S3hat(3,0:3,0:3),
     &    S4hat(3,0:3,0:3,3),S5hat(0:3,0:3,0:3,0:3,0:3)
      complex*16 S2mod(3,3),S3mod(3,0:3,0:3),
     &    S4mod(3,0:3,0:3,3),S5mod(0:3,0:3,0:3,0:3,0:3)
c      complex*16 S00,S001(3),S002(0:3,0:3),S003(0:3,0:3,0:3)
      complex*16 D200mod,D300mod(3),D400mod(3,3),D500mod(3,3,3)
      complex*16 D600mod(3,3,3,3)
      complex*16 D0,D1(3),D2(0:3,0:3),D3(0:3,0:3,0:3)
      complex*16 D4(0:3,0:3,0:3,0:3),D5(0:3,0:3,0:3,0:3,0:3)
      complex*16 cD0f
      complex*16 elimcminf2
c      real*8     maxzadjf,maxzadj
      integer    rank,switch,rankC
      integer    i1,i2,i3,i4,i5,j,k
      integer    ltest,sym
c      integer    testout

      common /ltest/  ltest
      common /sym/ sym

      data       C1m /12*0d0/,C2m /64*0d0/,C3m/256*0d0/,C4m/1024*0d0/
c      data       testout /45/

c      write(testout,*)'Dy12345 in',p10,p21,p32,p30,p20,p31,
c     &    m02,m12,m22,m32
c     &    ,rank

      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cDy12345'
        write(*,*) 'rank = ',rank
        stop
      else if (rank.gt.5) then
        write(*,*) 'rank > 5 not implemented in cDy12345'
        write(*,*) 'rank = ',rank
        stop
      else if (rank.eq.5.and.switch.ne.0) then
        write(*,*) 'rank = 5 and switch =/= 0',
     &      'not implemented in cDy12345'
        write(*,*) 'switch = ',switch
        stop
      end if

      D0 =cD0f(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32)
c      write(testout,01) D0

      rankC = max(rank-1,0)

c  cCp must be called also for rank = 0 because of cache!

c     write(*,*) 'cDp1',p21,p32,p31,m12,m22,m32,rankC
      call cCp12345(p21,p32,p31,m12,m22,m32,C0w0,C1w0,C2w0,C3w0,C4w0,
     &    C5w0,C6w0,rankC)
c     write(*,*) 'cDp2',p21,p32,p31,m12,m22,m32,rankC
      call cCp12345(p20,p32,p30,m02,m22,m32,C0w1,C1w1,C2w1,C3w1,C4w1,
     &    C5w1,C6w1,rankC)
c     write(*,*) 'cDp3',p21,p32,p31,m12,m22,m32,rankC
      call cCp12345(p10,p31,p30,m02,m12,m32,C0w2,C1w2,C2w2,C3w2,C4w2,
     &    C5w2,C6w2,rankC)
c     write(*,*) 'cDp4',p21,p32,p31,m12,m22,m32,rankC
      call cCp12345(p10,p21,p20,m02,m12,m22,C0w3,C1w3,C2w3,C3w3,C4w3,
     &    C5w3,C6w3,rankC)
      
      if (rank.eq.0) goto 100

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


c>      do i1=1,3
c>        D1(i1) = 0d0
c>        D3(0,0,i1) = 0d0
c>        D5(0,0,0,0,i1) = 0d0
c>        do i2=1,3
c>          D2(i1,i2) = 0d0
c>          D4(0,0,i1,i2) = 0d0
c>          do i3=1,3
c>            D3(i1,i2,i3) = 0d0
c>            D5(0,0,i1,i2,i3) = 0d0
c>            do i4=1,3
c>              D4(i1,i2,i3,i4) = 0d0
c>              do i5=1,3
c>                D5(i1,i2,i3,i4,i5) = 0d0
c>              end do
c>            end do
c>          end do
c>        end do
c>      end do

      mx(0,0) = 2d0*mm02
      mx(1,0) = q10 - mm12 + mm02
      mx(2,0) = q20 - mm22 + mm02
      mx(3,0) = q30 - mm32 + mm02
      mx(0,1) = mx(1,0)
      mx(0,2) = mx(2,0)
      mx(0,3) = mx(3,0)

      mx(1,1) = 2d0*q10
      mx(2,2) = 2d0*q20
      mx(3,3) = 2d0*q30
      mx(1,2) = q10+q20-q21
      mx(1,3) = q10+q30-q31
      mx(2,3) = q20+q30-q32
      mx(2,1) = mx(1,2)
      mx(3,1) = mx(1,3)
      mx(3,2) = mx(2,3)

c      do k=0,3
c      do j=0,3
c      write(*,*) 'mx =',k,j,mx(k,j)
c      end do
c      end do

      call chinv(4,mx,mxinv)

c>      write(*,*) 'x00 = ',mxinv(0,0)
c>      do k=1,3
c>      write(*,*) 'x0i = ',k,mxinv(0,k)
c>      write(*,*) 'x0i = ',k,mxinv(k,0)
c>      do j=1,3
c>      write(*,*) 'x0i = ',k,j,mxinv(k,j)
c>      end do
c>      end do


      C0m(0) = C0w0
      C0m(1) = C0w1
      C0m(2) = C0w2
      C0m(3) = C0w3
      
      do j =1,3
        S1hat(j) = C0m(j) - C0m(0)
      end do

      if (rank.eq.0) goto 99
  
      C1m(0,2) = C1w0(1)
      C1m(0,3) = C1w0(2)
      C1m(0,1) = -C0m(0) - C1m(0,2) - C1m(0,3)
      C1m(1,2) = C1w1(1)
      C1m(1,3) = C1w1(2)
      C1m(2,1) = C1w2(1)
      C1m(2,3) = C1w2(2)
      C1m(3,1) = C1w3(1)
      C1m(3,2) = C1w3(2)
      
      do i1=1,3
        do j =1,3
          S2hat(j,i1) = C1m(j,i1) - C1m(0,i1)
        end do
      end do

        if (rank.eq.1) goto 99

        C2m(0,0,0) = C2w0(0,0)
        C2m(1,0,0) = C2w1(0,0)
        C2m(2,0,0) = C2w2(0,0)
        C2m(3,0,0) = C2w3(0,0)
        
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
        
c        C2m(0,3,2) = C2m(0,2,3) 
c        C2m(0,3,1) = C2m(0,1,3) 
c        C2m(0,2,1) = C2m(0,1,2) 
c        C2m(1,3,2) = C2m(1,2,3) 
c        C2m(2,3,1) = C2m(2,1,3) 
c        C2m(3,2,1) = C2m(3,1,2) 

        do i1=1,3
          S3hat(i1,0,0) = C2m(i1,0,0) - C2m(0,0,0)
          do i2=i1,3
            do j =1,3
              S3hat(j,i1,i2) = C2m(j,i1,i2) - C2m(0,i1,i2)
            end do
          end do
        end do


        if (rank.eq.2) goto 99

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
        
c        C3m(0,1,2,1) = C3m(0,1,1,2)
c        C3m(0,2,1,1) = C3m(0,1,1,2)
c        C3m(0,1,3,1) = C3m(0,1,1,3)
c        C3m(0,3,1,1) = C3m(0,1,1,3)
c        C3m(0,2,2,1) = C3m(0,1,2,2)
c        C3m(0,2,1,2) = C3m(0,1,2,2)
c        C3m(0,3,3,1) = C3m(0,1,3,3)
c        C3m(0,3,1,3) = C3m(0,1,3,3)
c        C3m(0,2,3,2) = C3m(0,2,2,3)
c        C3m(0,3,2,2) = C3m(0,2,2,3)
c        C3m(0,3,3,2) = C3m(0,2,3,3)
c        C3m(0,3,2,3) = C3m(0,2,3,3)
        
c        C3m(0,1,3,2) = C3m(0,1,2,3)
c        C3m(0,2,1,3) = C3m(0,1,2,3)
c        C3m(0,2,3,1) = C3m(0,1,2,3)
c        C3m(0,3,1,2) = C3m(0,1,2,3)
c        C3m(0,3,2,1) = C3m(0,1,2,3)
        
c        C3m(3,1,2,1) = C3m(3,1,1,2)
c        C3m(3,2,1,1) = C3m(3,1,1,2)
c        C3m(2,1,3,1) = C3m(2,1,1,3)
c        C3m(2,3,1,1) = C3m(2,1,1,3)
c        C3m(3,2,2,1) = C3m(3,1,2,2)
c        C3m(3,2,1,2) = C3m(3,1,2,2)
c        C3m(2,3,3,1) = C3m(2,1,3,3)
c        C3m(2,3,1,3) = C3m(2,1,3,3)
c        C3m(1,2,3,2) = C3m(1,2,2,3)
c        C3m(1,3,2,2) = C3m(1,2,2,3)
c        C3m(1,3,3,2) = C3m(1,2,3,3)
c        C3m(1,3,2,3) = C3m(1,2,3,3)
        
        do j = 1,3
          do i1 = 1,3
            S4hat(j,0,0,i1) = C3m(j,0,0,i1) - C3m(0,0,0,i1)
          do i2 = i1,3
          do i3 = i2,3
            S4hat(j,i1,i2,i3) = C3m(j,i1,i2,i3) - C3m(0,i1,i2,i3)
          end do
          end do
          end do
        end do

        if (rank.eq.3) goto 99

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

c     >      do i1=1,2
c     >      do i2=i1,2
c     >        C4m(0,0,0,i1+1,i2+1) = C4w0(0,0,i1,i2)
c     >        C4m(1,0,0,i1+1,i2+1) = C4w0(0,0,i1,i2)
c     >        C4m(2,0,0,2*i1-1,2*i2-1) = C4w0(0,0,i1,i2)
c     >        C4m(3,0,0,i1+1,i2) = C4w0(0,0,i1,i2)
c     >      do i3=1,2
c     >      do i4=1,2
c     >        C4m(0,i1+1,i2+1,i3+1,i4+1) = C4w0(i1,i2,i3,i4)
c     >        C4m(1,i1+1,i2+1,i3+1,i4+1) = C4w0(i1,i2,i3,i4)
c     >        C4m(2,2*i1-1,2*i2-1,2*i3-1,2*i4-1) = C4w0(i1,i2,i3,i4)
c     >        C4m(3,i1+1,i2,i3,i4) = C4w0(i1,i2,i3,i4)
c     >      end do
c     >        C4m(0,1,i1+1,i2+1,i3+1) = -C3m(0,i1+1,i2+1,i3+1) 
c     >     &      - C4m(0,2,i1+1,i2+1,i3+1) 
c     >     &      - C4m(0,i1+1,i2+1,i3+1,3)
c     >      end do
c     >        C4m(0,1,1,i1+1,i2+1) = -C3m(0,1,i1+1,i2+1) 
c     >     &      - C4m(0,1,2,i1+1,i2+1) 
c     >     &      - C4m(0,1,i1+1,i2+1,3)
c     >      end do
c     >        C4m(0,1,1,1,i1+1) = -C3m(0,1,1,i1+1) 
c     >     &      - C4m(0,1,1,2,i1+1) 
c     >     &      - C4m(0,1,1,i1+1,3)
c     >      end do
c     >      C4m(0,1,1,1,1) = -C3m(0,1,1,1) 
c     >     &      - C4m(0,1,1,1,2) 
c     >     &      - C4m(0,1,1,1,3)

        do i1=1,3
          do i2=i1,3
            C4m(0,0,0,i2,i1) = C4m(0,0,0,i1,i2)
            C4m(1,0,0,i2,i1) = C4m(1,0,0,i1,i2)
            C4m(2,0,0,i2,i1) = C4m(2,0,0,i1,i2)
            C4m(3,0,0,i2,i1) = C4m(3,0,0,i1,i2)
            do i3=i2,3
              do i4=i3,3
                C4m(0,i1,i3,i2,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i1,i3,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i3,i1,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i1,i2,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i2,i1,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i2,i4,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i4,i2,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i1,i4,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i4,i1,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i1,i2,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i2,i1,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i3,i4,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i4,i3,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i1,i4,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i4,i1,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i1,i3,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i3,i1,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i2,i4,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i4,i2,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i3,i4,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i4,i3,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i3,i2,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i2,i3,i1) = C4m(0,i1,i2,i3,i4)

                C4m(1,i1,i3,i2,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i1,i3,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i3,i1,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i1,i2,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i2,i1,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i2,i4,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i4,i2,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i1,i4,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i4,i1,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i1,i2,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i2,i1,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i3,i4,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i4,i3,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i1,i4,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i4,i1,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i1,i3,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i3,i1,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i2,i4,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i4,i2,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i3,i4,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i4,i3,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i3,i2,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i2,i3,i1) = C4m(1,i1,i2,i3,i4)

                C4m(2,i1,i3,i2,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i1,i3,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i3,i1,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i1,i2,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i2,i1,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i2,i4,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i4,i2,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i1,i4,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i4,i1,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i1,i2,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i2,i1,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i3,i4,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i4,i3,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i1,i4,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i4,i1,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i1,i3,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i3,i1,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i2,i4,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i4,i2,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i3,i4,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i4,i3,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i3,i2,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i2,i3,i1) = C4m(2,i1,i2,i3,i4)

                C4m(3,i1,i3,i2,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i1,i3,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i3,i1,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i1,i2,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i2,i1,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i2,i4,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i4,i2,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i1,i4,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i4,i1,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i1,i2,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i2,i1,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i3,i4,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i4,i3,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i1,i4,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i4,i1,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i1,i3,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i3,i1,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i2,i4,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i4,i2,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i3,i4,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i4,i3,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i3,i2,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i2,i3,i1) = C4m(3,i1,i2,i3,i4)
              end do
            end do
          end do
        end do

        do k = 1,3
          S5hat(k,0,0,0,0) = C4m(k,0,0,0,0) - C4m(0,0,0,0,0)
          do i1 = 1,3
          do i2 = i1,3
            S5hat(k,0,0,i1,i2) = C4m(k,0,0,i1,i2) - C4m(0,0,0,i1,i2)
          do i3 = i2,3
          do i4 = i3,3
            S5hat(k,i1,i2,i3,i4) = C4m(k,i1,i2,i3,i4) -
     &          C4m(0,i1,i2,i3,i4)
          end do
          end do
          end do
          end do
        end do

 99     continue

 01   format(' D0y = ',G20.14,' + i* ',G20.14)
 11   format(' D1y(',i1,') = ',G20.14,' + i* ',G20.14)
 20   format(' D2y(0,0) = ',G20.14,' + i* ',G20.14)
 22   format(' D2y(',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 31   format(' D3y(0,0,',i1,') = ',G20.14,' + i* ',G20.14)
 33   format(' D3y(',i1,',',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 40   format(' D4y(0,0,0,0) = ',G20.14,' + i* ',G20.14)
 42   format(' D4y(0,0,',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 44   format(' D4y(',i1,',',i1,',',i1,',',i1,') = ',
     &  G20.14,' + i* ',G20.14)
 51   format(' D5y(0,0,0,0,',i1,') = ',G20.14,' + i* ',G20.14)
 53   format(' D5y(0,0,',i1,',',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 55   format(' D5y(',i1,',',i1,',',i1,',',i1,',',i1,') = ',
     &  G20.14,' + i* ',G20.14)

      D200mod =  (D0 
     &          - mxinv(0,1)*S1hat(1) 
     &          - mxinv(0,2)*S1hat(2) 
     &          - mxinv(0,3)*S1hat(3))
     &           /mxinv(0,0)

      do i1 = 1,3
        D1(i1) = mxinv(i1,0)*D200mod 
     &         + mxinv(i1,1)*S1hat(1)
     &         + mxinv(i1,2)*S1hat(2)
     &         + mxinv(i1,3)*S1hat(3)
c        write(testout,11) i1,D1(i1)
      end do
        
      if (rank.eq.1.and.switch.eq.0) goto 100

      D2(0,0) = (D200mod + C0m(0))/2d0
c      write(testout,20) D2(0,0)

      if (rank.eq.1) goto 100

      do i1=1,3
        do j =1,3
          S2mod(j,i1) = S2hat(j,i1)
        end do
        S2mod(i1,i1) = S2mod(i1,i1) - 2d0*D2(0,0)
      end do

      do i1 = 1,3
        D300mod(i1) =  (D1(i1) 
     &         - mxinv(0,1)*S2mod(1,i1) 
     &         - mxinv(0,2)*S2mod(2,i1) 
     &         - mxinv(0,3)*S2mod(3,i1))
     &           /mxinv(0,0)
      end do

      do i1 = 1,3
      do i2 = i1,3
        D2(i1,i2) = mxinv(i1,0)*D300mod(i2) 
     &         + mxinv(i1,1)*S2mod(1,i2)
     &         + mxinv(i1,2)*S2mod(2,i2)
     &         + mxinv(i1,3)*S2mod(3,i2)
c        write(testout,22) i1,i2,D2(i1,i2)
      end do
      end do

      if (sym.eq.1) then
        D2(2,1) = D2(1,2)
        D2(3,1) = D2(1,3)
        D2(3,2) = D2(2,3)
      end if
        
      if (rank.eq.2.and.switch.eq.0) goto 100

      do i1=1,3
        D3(0,0,i1) = (D300mod(i1) + C1m(0,i1))/4d0
c        write(testout,31) i1,D3(0,0,i1)
      end do
      
      if (rank.eq.2) goto 100

      do i1=1,3
        do i2=i1,3
          do j =1,3
            S3mod(j,i1,i2) = S3hat(j,i1,i2)
          end do
          S3mod(i1,i1,i2) = S3mod(i1,i1,i2) - 2d0*D3(0,0,i2)
          S3mod(i2,i1,i2) = S3mod(i2,i1,i2) - 2d0*D3(0,0,i1)
        end do
      end do

      do i1 = 1,3
      do i2 = i1,3
        D400mod(i1,i2) =  (D2(i1,i2) 
     &         - mxinv(0,1)*S3mod(1,i1,i2) 
     &         - mxinv(0,2)*S3mod(2,i1,i2) 
     &         - mxinv(0,3)*S3mod(3,i1,i2))
     &           /mxinv(0,0)
      end do
      end do

      do i1 = 1,3
      do i2 = i1,3
      do i3 = i2,3
        D3(i1,i2,i3) = mxinv(i1,0)*D400mod(i2,i3) 
     &         + mxinv(i1,1)*S3mod(1,i2,i3)
     &         + mxinv(i1,2)*S3mod(2,i2,i3)
     &         + mxinv(i1,3)*S3mod(3,i2,i3)
c        write(testout,33) i1,i2,i3,D3(i1,i2,i3)
      end do
      end do
      end do
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

      if (rank.eq.3.and.switch.eq.0) goto 100

      D4(0,0,0,0) = (
     &    (D2(0,0) 
     &         - mxinv(0,1)*S3hat(1,0,0)  
     &         - mxinv(0,2)*S3hat(2,0,0) 
     &         - mxinv(0,3)*S3hat(3,0,0)
     &           )/mxinv(0,0)
     &    + C2m(0,0,0) 
     &    + 1d0/6d0)/6d0
c      write(testout,40) D4(0,0,0,0)

      do i1=1,3
      do i2=i1,3
        D4(0,0,i1,i2) = (D400mod(i1,i2) + C2m(0,i1,i2))/6d0
      end do
      end do
      if (sym.eq.1) then
        D4(0,0,2,1) = D4(0,0,1,2)
        D4(0,0,3,1) = D4(0,0,1,3)
        D4(0,0,3,2) = D4(0,0,2,3)
      end if

      if (rank.eq.3) goto 100

      do i1=1,3
        do j =1,3
          S4mod(j,0,0,i1) = S4hat(j,0,0,i1)
        end do
        S4mod(i1,0,0,i1) = S4mod(i1,0,0,i1) - 2d0*D4(0,0,0,0)
      do i2=i1,3
      do i3=i2,3
        do j =1,3
          S4mod(j,i1,i2,i3) = S4hat(j,i1,i2,i3)
        end do
        S4mod(i1,i1,i2,i3) = S4mod(i1,i1,i2,i3) -2d0*D4(0,0,i2,i3)
        S4mod(i2,i1,i2,i3) = S4mod(i2,i1,i2,i3) -2d0*D4(0,0,i1,i3)
        S4mod(i3,i1,i2,i3) = S4mod(i3,i1,i2,i3) -2d0*D4(0,0,i1,i2)
      end do
      end do
      end do
        
      do i1 = 1,3
      do i2 = i1,3
      do i3 = i2,3
        D500mod(i1,i2,i3) =  (D3(i1,i2,i3) 
     &         - mxinv(0,1)*S4mod(1,i1,i2,i3) 
     &         - mxinv(0,2)*S4mod(2,i1,i2,i3) 
     &         - mxinv(0,3)*S4mod(3,i1,i2,i3))
     &           /mxinv(0,0)
      end do
      end do
      end do

      do i1 = 1,3
      do i2 = i1,3
      do i3 = i2,3
      do i4 = i3,3
        D4(i1,i2,i3,i4) = mxinv(i1,0)*D500mod(i2,i3,i4) 
     &         + mxinv(i1,1)*S4mod(1,i2,i3,i4)
     &         + mxinv(i1,2)*S4mod(2,i2,i3,i4)
     &         + mxinv(i1,3)*S4mod(3,i2,i3,i4)
c        write(testout,44) i1,i2,i3,i4,D4(i1,i2,i3,i4)
      end do
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
      end do
      end do
      end do

      if (rank.eq.4.and.switch.eq.0) goto 100

      do i1 = 1,3
        D5(0,0,0,0,i1) = (
     &    (D3(0,0,i1) 
     &         - mxinv(0,1)*S4mod(1,0,0,i1)  
     &         - mxinv(0,2)*S4mod(2,0,0,i1) 
     &         - mxinv(0,3)*S4mod(3,0,0,i1)
     &           )/mxinv(0,0)
     &    + C3m(0,0,0,i1) 
     &    - 1d0/24d0)/8d0
c          write(testout,51) i1,D5(0,0,0,0,i1)
      end do
      
      do i1=1,3
      do i2=i1,3
      do i3=i2,3
        D5(0,0,i1,i2,i3) = (D500mod(i1,i2,i3) + C3m(0,i1,i2,i3))/8d0
c          write(testout,53) i1,i2,i3,D5(0,0,i1,i2,i3)
      end do
      end do
      end do
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

      if (rank.eq.4) goto 100

      do i1=1,3
      do i2=i1,3
      do i3=i2,3
      do i4=i3,3
        do j =1,3
          S5mod(j,i1,i2,i3,i4) = S5hat(j,i1,i2,i3,i4)
        end do
        S5mod(i1,i1,i2,i3,i4) = S5mod(i1,i1,i2,i3,i4) 
     &      - 2d0*D5(0,0,i2,i3,i4)
        S5mod(i2,i1,i2,i3,i4) = S5mod(i2,i1,i2,i3,i4) 
     &      - 2d0*D5(0,0,i1,i3,i4)
        S5mod(i3,i1,i2,i3,i4) = S5mod(i3,i1,i2,i3,i4) 
     &      - 2d0*D5(0,0,i1,i2,i4)
        S5mod(i4,i1,i2,i3,i4) = S5mod(i4,i1,i2,i3,i4) 
     &      - 2d0*D5(0,0,i1,i2,i3)
      end do
      end do
      end do
      end do

      do i1 = 1,3
      do i2 = i1,3
      do i3 = i2,3
      do i4 = i3,3
        D600mod(i1,i2,i3,i4) =  (D4(i1,i2,i3,i4) 
     &         - mxinv(0,1)*S5mod(1,i1,i2,i3,i4) 
     &         - mxinv(0,2)*S5mod(2,i1,i2,i3,i4) 
     &         - mxinv(0,3)*S5mod(3,i1,i2,i3,i4))
     &           /mxinv(0,0)
      end do
      end do
      end do
      end do

      do i1 = 1,3
      do i2 = i1,3
      do i3 = i2,3
      do i4 = i3,3
      do i5 = i4,3
        D5(i1,i2,i3,i4,i5) = mxinv(i1,0)*D600mod(i2,i3,i4,i5) 
     &         + mxinv(i1,1)*S5mod(1,i2,i3,i4,i5)
     &         + mxinv(i1,2)*S5mod(2,i2,i3,i4,i5)
     &         + mxinv(i1,3)*S5mod(3,i2,i3,i4,i5)
c        write(testout,55) i1,i2,i3,i4,i5,D5(i1,i2,i3,i4,i5)
      end do
        if (sym.eq.1) then
          do k =1,i1
c           D5(i1,i2,i3,i4,k) = D5(k,i1,i2,i3,i4)
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
          end do
        end if
      end do
      end do
      end do
      end do

 100  continue

      end

***********************************************************************
      subroutine cDg12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank)
***********************************************************************
*     4-point tensor coefficient functions D0,D1,D2,D3,D4,D5          *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     expansion about Gram determinant = 0 to order ordg4             *
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
*     i-th denominator cancelled: Cxwi                                *
***********************************************************************
*     11.11.04 Ansgar Denner    last changed  13.05.05                *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p32,p30,p20,p31
      complex*16 q10,q21,q32,q30,q20,q31
      complex*16 m02,m12,m22,m32
      complex*16 mm02,mm12,mm22,mm32,f(3),zadjf(3)
      complex*16 k1k2,k1k3,k2k3,z(3,3),zadj(3,3),detz
      complex*16 zadj2(3,3,3,3)
      complex*16 C0w3,C1w3(2),C2w3(0:2,0:2),C3w3(0:2,0:2,0:2)
     &    ,C4w3(0:2,0:2,0:2,0:2),C5w3(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w3(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w2,C1w2(2),C2w2(0:2,0:2),C3w2(0:2,0:2,0:2)
     &    ,C4w2(0:2,0:2,0:2,0:2),C5w2(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w2(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w1,C1w1(2),C2w1(0:2,0:2),C3w1(0:2,0:2,0:2)
     &    ,C4w1(0:2,0:2,0:2,0:2),C5w1(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w1(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w0,C1w0(2),C2w0(0:2,0:2),C3w0(0:2,0:2,0:2)
     &    ,C4w0(0:2,0:2,0:2,0:2),C5w0(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w0(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0m(0:3),C1m(0:3,3),C2m(0:3,0:3,0:3)
     &    ,C3m(0:3,0:3,0:3,0:3),C4m(0:3,0:3,0:3,0:3,0:3)
c     &    ,C5m(0:3,0:3,0:3,0:3,0:3,0:3),C6m(0:3,0:3,0:3,0:3,0:3,0:3,0:3)
      complex*16 S1hat(3),S2hat(3,3),S3hat(3,0:3,0:3),
     &    S4hat(3,0:3,0:3,3),S5hat(0:3,0:3,0:3,0:3,0:3)
      complex*16 S2mod(3,3),S3mod(3,0:3,0:3),
     &    S4mod(3,0:3,0:3,3),S5mod(0:3,0:3,0:3,0:3,0:3)
      complex*16 S00,S001(3),S002(0:3,0:3),S003(0:3,0:3,0:3)
      complex*16 D0,D1(3),D2(0:3,0:3),D3(0:3,0:3,0:3)
      complex*16 D4(0:3,0:3,0:3,0:3),D5(0:3,0:3,0:3,0:3,0:3)
c     complex*16 cD0f
      complex*16 elimcminf2
      real*8     maxzadjf,maxzadj
      integer    rank
      integer    i1,i2,i3,i4,j,k,l,m,a,b,r
      integer    ltest,sym
c      integer    testout

      integer    ordg3,ordgy3,ordgp3
      integer    ordg4,ordgy4,ordgp4
      common /ctord3/ ordg3,ordgy3,ordgp3
      common /ctord4/ ordg4,ordgy4,ordgp4
      common /ltest/  ltest
      common /sym/ sym

      data       C1m /12*0d0/,C2m /64*0d0/,C3m/256*0d0/,C4m/1024*0d0/
      data       zadj2/81*0d0/
c      data       testout /47/

c      write(*,*)'Dg12345 in',p10,p21,p32,p30,p20,p31,
c      write(testout,*)'Dg12345 in',p10,p21,p32,p30,p20,p31,
c     &    m02,m12,m22,m32
c     &    ,rank,ordg4

      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cDg12345'
        write(*,*) 'rank = ',rank
        stop
      else if (rank+ordg4.gt.4) then
        write(*,*) 'rank+ordg4 > 4 not implemented in cDg12345'
        write(*,*) 'rank  = ',rank
        write(*,*) 'ordg4 = ',ordg4
        if (rank.gt.4) then
          write(*,*) 'rank > 4 not implemented in cDg12345'
          write(*,*) 'rank = ',rank
          stop
        end if
      end if

c      write(*,*) 'cDg12345: ordg4 = ',ordg4

c     calculate coefficients
      

c     write(*,*) 'cDp1',p21,p32,p31,m12,m22,m32,rank
      call cCp12345(p21,p32,p31,m12,m22,m32,C0w0,C1w0,C2w0,C3w0,C4w0,
     &    C5w0,C6w0,rank+ordg4)
c     write(*,*) 'cDp2',p21,p32,p31,m12,m22,m32,rank
      call cCp12345(p20,p32,p30,m02,m22,m32,C0w1,C1w1,C2w1,C3w1,C4w1,
     &    C5w1,C6w1,rank+ordg4)
c     write(*,*) 'cDp3',p21,p32,p31,m12,m22,m32,rank
      call cCp12345(p10,p31,p30,m02,m12,m32,C0w2,C1w2,C2w2,C3w2,C4w2,
     &    C5w2,C6w2,rank+ordg4)
c     write(*,*) 'cDp4',p21,p32,p31,m12,m22,m32,rank
      call cCp12345(p10,p21,p20,m02,m12,m22,C0w3,C1w3,C2w3,C3w3,C4w3,
     &    C5w3,C6w3,rank+ordg4)
      
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

      z(1,1) = 2d0*q10
      z(2,2) = 2d0*q20
      z(3,3) = 2d0*q30
      z(1,2) = q10+q20-q21
      z(1,3) = q10+q30-q31
      z(2,3) = q20+q30-q32
      z(2,1) = z(1,2)
      z(3,1) = z(1,3)
      z(3,2) = z(2,3)
      
      k1k2 = (q10+q20-q21)
      k1k3 = (q10+q30-q31)
      k2k3 = (q20+q30-q32)
      detz = 8d0*q10*q30*q20+2D0*k1k2*k1k3*k2k3
     &    -2d0*(q10*k2k3*k2k3+q20*k1k3*k1k3+q30*k1k2*k1k2)

c      write(testout,*) 'detz = ',detz,detz/
c     &    dmax1(abs(k1k2),abs(k1k3),abs(k2k3),
c     &    2*abs(q10),2*abs(q20),2*abs(q30))**3

      zadj(1,1) = (4d0*q30*q20-k2k3*k2k3)
      zadj(1,2) = (k1k3*k2k3-2d0*q30*k1k2)
      zadj(1,3) = (k1k2*k2k3-2d0*q20*k1k3)
      zadj(2,1) = zadj(1,2)
      zadj(2,2) = (4d0*q10*q30-k1k3*k1k3)
      zadj(2,3) = (k1k2*k1k3-2d0*q10*k2k3)
      zadj(3,1) = zadj(1,3)
      zadj(3,2) = zadj(2,3)
      zadj(3,3) = (4d0*q10*q20-k1k2*k1k2)

      zadj2(1,2,1,2) = -z(3,3)
      zadj2(1,2,2,1) = z(3,3)
      zadj2(2,1,1,2) = z(3,3)
      zadj2(2,1,2,1) = -z(3,3)
      zadj2(1,3,1,3) = -z(2,2)
      zadj2(1,3,3,1) = z(2,2)
      zadj2(3,1,1,3) = z(2,2)
      zadj2(3,1,3,1) = -z(2,2)
      zadj2(3,2,3,2) = -z(1,1)
      zadj2(3,2,2,3) = z(1,1)
      zadj2(2,3,3,2) = z(1,1)
      zadj2(2,3,2,3) = -z(1,1)
      zadj2(1,2,1,3) = z(3,2)
      zadj2(1,2,3,1) = -z(3,2)
      zadj2(2,1,1,3) = -z(3,2)
      zadj2(2,1,3,1) = z(3,2)
      zadj2(1,3,1,2) = z(2,3)
      zadj2(1,3,2,1) = -z(2,3)
      zadj2(3,1,1,2) = -z(2,3)
      zadj2(3,1,2,1) = z(2,3)
      zadj2(3,2,1,2) = z(1,3)
      zadj2(3,2,2,1) = -z(1,3)
      zadj2(2,3,1,2) = -z(1,3)
      zadj2(2,3,2,1) = z(1,3)
      zadj2(1,2,3,2) = z(3,1)
      zadj2(1,2,2,3) = -z(3,1)
      zadj2(2,1,3,2) = -z(3,1)
      zadj2(2,1,2,3) = z(3,1)
      zadj2(1,3,2,3) = z(2,1)
      zadj2(1,3,3,2) = -z(2,1)
      zadj2(3,1,2,3) = -z(2,1)
      zadj2(3,1,3,2) = z(2,1)
      zadj2(3,2,3,1) = z(1,2)
      zadj2(3,2,1,3) = -z(1,2)
      zadj2(2,3,3,1) = -z(1,2)
      zadj2(2,3,1,3) = z(1,2)

      f(1) = q10+mm02-mm12
      f(2) = q20+mm02-mm22
      f(3) = q30+mm02-mm32
      
      zadjf(1) = zadj(1,1)*f(1)+zadj(1,2)*f(2)+zadj(1,3)*f(3)
      zadjf(2) = zadj(2,1)*f(1)+zadj(2,2)*f(2)+zadj(2,3)*f(3)
      zadjf(3) = zadj(3,1)*f(1)+zadj(3,2)*f(2)+zadj(3,3)*f(3)

c      write(testout,*) 'zadjf ',zadjf
c      write(testout,*) 'zadjf ',detz/zadjf(1),detz/zadjf(2),detz/zadjf(3)

      do i1=1,3
        D1(i1) = 0d0
        D3(0,0,i1) = 0d0
        D5(0,0,0,0,i1) = 0d0
        do i2=1,3
          D2(i1,i2) = 0d0
          D4(0,0,i1,i2) = 0d0
          do i3=1,3
            D3(i1,i2,i3) = 0d0
            D5(0,0,i1,i2,i3) = 0d0
            do i4=1,3
              D4(i1,i2,i3,i4) = 0d0
c              do i5=1,3
c                D5(i1,i2,i3,i4,i5) = 0d0
c              end do
            end do
          end do
        end do
      end do

      C0m(0) = C0w0
      C0m(1) = C0w1
      C0m(2) = C0w2
      C0m(3) = C0w3
      
      do j =1,3
        S1hat(j) = C0m(j) - C0m(0)
      end do

      if (rank+ordg4.eq.0) goto 99
  
      C1m(0,2) = C1w0(1)
      C1m(0,3) = C1w0(2)
      C1m(0,1) = -C0m(0) - C1m(0,2) - C1m(0,3)
      C1m(1,2) = C1w1(1)
      C1m(1,3) = C1w1(2)
      C1m(2,1) = C1w2(1)
      C1m(2,3) = C1w2(2)
      C1m(3,1) = C1w3(1)
      C1m(3,2) = C1w3(2)
      
      do i1=1,3
        do j =1,3
          S2hat(j,i1) = C1m(j,i1) - C1m(0,i1)
        end do
      end do

c     >      do i1=1,3
c     >      do i2=1,3
c     >      do k=1,3
c     >        do l=1,3
c     >          test = zadj(k,l)*z(i1,i2)
c     >          do a=1,3
c     >            do b=1,3
c     >              if(b.ne.k.and.a.ne.l) then
c     >                test  = test + zadj2(k,b,l,a)*z(b,i1)*z(a,i2)
c     >              end if
c     >            end do
c     >          end do
c     >          write(testout,*) 'detz = ',i1,i2,k,l,test
c     >        end do
c     >      end do
c     >      end do
c     >      end do

        if (rank+ordg4.eq.1) goto 99

        C2m(0,0,0) = C2w0(0,0)
        C2m(1,0,0) = C2w1(0,0)
        C2m(2,0,0) = C2w2(0,0)
        C2m(3,0,0) = C2w3(0,0)
        
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
        
        C2m(0,3,2) = C2m(0,2,3) 
        C2m(0,3,1) = C2m(0,1,3) 
        C2m(0,2,1) = C2m(0,1,2) 
        C2m(1,3,2) = C2m(1,2,3) 
        C2m(2,3,1) = C2m(2,1,3) 
        C2m(3,2,1) = C2m(3,1,2) 

        do i1=1,3
          S3hat(i1,0,0) = C2m(i1,0,0) - C2m(0,0,0)
          do i2=1,3
            do j =1,3
              S3hat(j,i1,i2) = C2m(j,i1,i2) - C2m(0,i1,i2)
            end do
          end do
        end do


        if (rank+ordg4.eq.2) goto 99

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
        
        C3m(0,1,2,1) = C3m(0,1,1,2)
        C3m(0,2,1,1) = C3m(0,1,1,2)
        C3m(0,1,3,1) = C3m(0,1,1,3)
        C3m(0,3,1,1) = C3m(0,1,1,3)
        C3m(0,2,2,1) = C3m(0,1,2,2)
        C3m(0,2,1,2) = C3m(0,1,2,2)
        C3m(0,3,3,1) = C3m(0,1,3,3)
        C3m(0,3,1,3) = C3m(0,1,3,3)
        C3m(0,2,3,2) = C3m(0,2,2,3)
        C3m(0,3,2,2) = C3m(0,2,2,3)
        C3m(0,3,3,2) = C3m(0,2,3,3)
        C3m(0,3,2,3) = C3m(0,2,3,3)
        
        C3m(0,1,3,2) = C3m(0,1,2,3)
        C3m(0,2,1,3) = C3m(0,1,2,3)
        C3m(0,2,3,1) = C3m(0,1,2,3)
        C3m(0,3,1,2) = C3m(0,1,2,3)
        C3m(0,3,2,1) = C3m(0,1,2,3)
        
        C3m(3,1,2,1) = C3m(3,1,1,2)
        C3m(3,2,1,1) = C3m(3,1,1,2)
        C3m(2,1,3,1) = C3m(2,1,1,3)
        C3m(2,3,1,1) = C3m(2,1,1,3)
        C3m(3,2,2,1) = C3m(3,1,2,2)
        C3m(3,2,1,2) = C3m(3,1,2,2)
        C3m(2,3,3,1) = C3m(2,1,3,3)
        C3m(2,3,1,3) = C3m(2,1,3,3)
        C3m(1,2,3,2) = C3m(1,2,2,3)
        C3m(1,3,2,2) = C3m(1,2,2,3)
        C3m(1,3,3,2) = C3m(1,2,3,3)
        C3m(1,3,2,3) = C3m(1,2,3,3)
        
        do i1=1,3
          do i2=1,3
            S4hat(i1,0,0,i2) = C3m(i1,0,0,i2) - C3m(0,0,0,i2)
            do i3=1,3
              do j =1,3
                S4hat(j,i1,i2,i3) = C3m(j,i1,i2,i3) - C3m(0,i1,i2,i3)
              end do
            end do
          end do
        end do
        

        if (rank+ordg4.eq.3) goto 99

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

c     >      do i1=1,2
c     >      do i2=i1,2
c     >        C4m(0,0,0,i1+1,i2+1) = C4w0(0,0,i1,i2)
c     >        C4m(1,0,0,i1+1,i2+1) = C4w0(0,0,i1,i2)
c     >        C4m(2,0,0,2*i1-1,2*i2-1) = C4w0(0,0,i1,i2)
c     >        C4m(3,0,0,i1+1,i2) = C4w0(0,0,i1,i2)
c     >      do i3=1,2
c     >      do i4=1,2
c     >        C4m(0,i1+1,i2+1,i3+1,i4+1) = C4w0(i1,i2,i3,i4)
c     >        C4m(1,i1+1,i2+1,i3+1,i4+1) = C4w0(i1,i2,i3,i4)
c     >        C4m(2,2*i1-1,2*i2-1,2*i3-1,2*i4-1) = C4w0(i1,i2,i3,i4)
c     >        C4m(3,i1+1,i2,i3,i4) = C4w0(i1,i2,i3,i4)
c     >      end do
c     >        C4m(0,1,i1+1,i2+1,i3+1) = -C3m(0,i1+1,i2+1,i3+1) 
c     >     &      - C4m(0,2,i1+1,i2+1,i3+1) 
c     >     &      - C4m(0,i1+1,i2+1,i3+1,3)
c     >      end do
c     >        C4m(0,1,1,i1+1,i2+1) = -C3m(0,1,i1+1,i2+1) 
c     >     &      - C4m(0,1,2,i1+1,i2+1) 
c     >     &      - C4m(0,1,i1+1,i2+1,3)
c     >      end do
c     >        C4m(0,1,1,1,i1+1) = -C3m(0,1,1,i1+1) 
c     >     &      - C4m(0,1,1,2,i1+1) 
c     >     &      - C4m(0,1,1,i1+1,3)
c     >      end do
c     >      C4m(0,1,1,1,1) = -C3m(0,1,1,1) 
c     >     &      - C4m(0,1,1,1,2) 
c     >     &      - C4m(0,1,1,1,3)

        do i1=1,3
          do i2=i1,3
            C4m(0,0,0,i2,i1) = C4m(0,0,0,i1,i2)
            C4m(1,0,0,i2,i1) = C4m(1,0,0,i1,i2)
            C4m(2,0,0,i2,i1) = C4m(2,0,0,i1,i2)
            C4m(3,0,0,i2,i1) = C4m(3,0,0,i1,i2)
            do i3=i2,3
              do i4=i3,3
                C4m(0,i1,i3,i2,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i1,i3,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i3,i1,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i1,i2,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i2,i1,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i2,i4,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i4,i2,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i1,i4,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i4,i1,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i1,i2,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i2,i1,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i3,i4,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i4,i3,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i1,i4,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i4,i1,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i1,i3,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i3,i1,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i2,i4,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i4,i2,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i3,i4,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i4,i3,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i3,i2,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i2,i3,i1) = C4m(0,i1,i2,i3,i4)

                C4m(1,i1,i3,i2,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i1,i3,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i3,i1,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i1,i2,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i2,i1,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i2,i4,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i4,i2,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i1,i4,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i4,i1,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i1,i2,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i2,i1,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i3,i4,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i4,i3,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i1,i4,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i4,i1,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i1,i3,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i3,i1,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i2,i4,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i4,i2,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i3,i4,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i4,i3,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i3,i2,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i2,i3,i1) = C4m(1,i1,i2,i3,i4)

                C4m(2,i1,i3,i2,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i1,i3,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i3,i1,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i1,i2,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i2,i1,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i2,i4,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i4,i2,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i1,i4,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i4,i1,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i1,i2,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i2,i1,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i3,i4,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i4,i3,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i1,i4,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i4,i1,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i1,i3,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i3,i1,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i2,i4,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i4,i2,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i3,i4,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i4,i3,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i3,i2,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i2,i3,i1) = C4m(2,i1,i2,i3,i4)

                C4m(3,i1,i3,i2,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i1,i3,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i3,i1,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i1,i2,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i2,i1,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i2,i4,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i4,i2,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i1,i4,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i4,i1,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i1,i2,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i2,i1,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i3,i4,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i4,i3,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i1,i4,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i4,i1,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i1,i3,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i3,i1,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i2,i4,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i4,i2,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i3,i4,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i4,i3,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i3,i2,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i2,i3,i1) = C4m(3,i1,i2,i3,i4)
              end do
            end do
          end do
        end do

        do i1=1,3
          S5hat(i1,0,0,0,0) = C4m(i1,0,0,0,0) - C4m(0,0,0,0,0)
          do i2=1,3
            do i3=1,3
              S5hat(i1,0,0,i2,i3) = C4m(i1,0,0,i2,i3) - C4m(0,0,0,i2,i3)
              do i4=1,3
                do k =1,3
                  S5hat(k,i1,i2,i3,i4) = C4m(k,i1,i2,i3,i4) -
     &                C4m(0,i1,i2,i3,i4)
                end do
              end do
            end do
          end do
      end do

 99     continue

c choose reduction formulas with biggest denominators

      maxzadjf = abs(zadjf(1))
      m = 1
      if (abs(zadjf(2)).gt.maxzadjf) then
        maxzadjf = abs(zadjf(2))
        m = 2
      end if
      if (abs(zadjf(3)).gt.maxzadjf) then
        maxzadjf = abs(zadjf(3))
        m = 3
      end if
        
      maxzadj = abs(zadj(1,1))
      k = 1
      l = 1
      if (abs(zadj(2,2)).gt.maxzadj) then
        maxzadj = abs(zadj(2,2))
        k = 2
        l = 2
      end if
      if (abs(zadj(3,3)).gt.maxzadj) then
        maxzadj = abs(zadj(3,3))
        k = 3
        l = 3
      end if
      if (abs(zadj(1,2)).gt.maxzadj) then
        maxzadj = abs(zadj(1,2))
        k = 1
        l = 2
      end if
      if (abs(zadj(1,3)).gt.maxzadj) then
        maxzadj = abs(zadj(1,3))
        k = 1
        l = 3
      end if
      if (abs(zadj(2,3)).gt.maxzadj) then
        maxzadj = abs(zadj(2,3))
        k = 2
        l = 3
      end if

c      write(*,*) 'klm',k,l,m

c      write(testout,*) 'ordg4 = ',ordg4
      do 100 r=0,ordg4+rank
c        write(testout,*) 'ordg4 = ',ordg4,r,m
        
        D0 = (zadj(m,1)*S1hat(1) + zadj(m,2)*S1hat(2) 
     &      + zadj(m,3)*S1hat(3) 
     &      - detz*D1(m))/zadjf(m)
c        write(testout,9) 'D0  = ',D0
 9      format(1x,a10,2g24.16,2i2)  
 10     format(1x,a10,i2,2g24.16,2i2)  
 20     format(1x,a10,2i2,2g24.16,2i2)  
 30     format(1x,a10,3i2,2g24.16)  
 40     format(1x,a10,4i2,2g24.16)  
 50     format(1x,a10,5i2,2g24.16)  


        if (rank+ordg4-r.eq.0) goto 100

        S00 = 2d0*mm02*D0 + 2d0*C0m(0)

        D2(0,0) = zadj(k,l)*S00 - detz*D2(k,l)
        do a=1,3
          D2(0,0) = D2(0,0) - zadj(k,l)*S2hat(a,a) 
     &        + zadj(a,l)*S2hat(a,k)
          do b=1,3
            if(b.ne.k.and.a.ne.l) then
              D2(0,0) = D2(0,0) - zadj2(k,b,l,a)*( 
     &            + f(b)*S1hat(a)
     &            - f(a)*f(b)*D0   )
            end if
          end do
        end do
        D2(0,0) = D2(0,0) /(4d0*zadj(k,l))
c        write(testout,9) 'D00  = ',D2(0,0)

        do i1=1,3
          do j =1,3
            S2mod(j,i1) = S2hat(j,i1)
          end do
          S2mod(i1,i1) = S2mod(i1,i1) - 2d0*D2(0,0)
        end do

        do i1=1,3
          D1(i1) = (zadj(m,1)*S2mod(1,i1) + zadj(m,2)*S2mod(2,i1) 
     &        + zadj(m,3)*S2mod(3,i1) 
     &        - detz*D2(m,i1))/zadjf(m)
c          write(testout,10) 'D1  = ',i1,D1(i1)
        end do

        if (rank+ordg4-r.eq.1) goto 100

        do i1=1,3
          S001(i1) = 2d0*mm02*D1(i1) + 2d0*C1m(0,i1)
        end do

        do i1=1,3
          D3(0,0,i1) = zadj(k,l)*S001(i1) - detz*D3(k,l,i1)
          do a=1,3
            D3(0,0,i1) = D3(0,0,i1) - zadj(k,l)*S3hat(a,a,i1) 
     &          + zadj(a,l)*S3hat(a,k,i1)
            do b=1,3
              if(b.ne.k.and.a.ne.l) then
                D3(0,0,i1) = D3(0,0,i1) - zadj2(k,b,l,a)*( 
     &              + f(b)*S2hat(a,i1)
     &              - f(a)*f(b)*D1(i1)   )
              end if
            end do
            if(i1.ne.k.and.a.ne.l) then
              D3(0,0,i1) = D3(0,0,i1) - 2d0*zadj2(k,i1,l,a)*( 
     &            + S3hat(a,0,0)
     &            - f(a)*D2(0,0)  )
            end if
          end do
          do b=1,3
            if(b.ne.k.and.i1.ne.l) then
              D3(0,0,i1) = D3(0,0,i1) - 2d0*zadj2(k,b,l,i1)*( 
     &            - f(b)*D2(0,0)   )
            end if
          end do
          D3(0,0,i1) = D3(0,0,i1) /(8d0*zadj(k,l))
c          write(testout,20) 'D00i  = ',i1,D3(0,0,i1)
        end do
        
        do i1=1,3
          do i2=1,3
            do j =1,3
              S3mod(j,i1,i2) = S3hat(j,i1,i2)
            end do
            S3mod(i1,i1,i2) = S3mod(i1,i1,i2) - 2d0*D3(0,0,i2)
            S3mod(i2,i1,i2) = S3mod(i2,i1,i2) - 2d0*D3(0,0,i1)
          end do
        end do
        
        do i1=1,3
          do i2=1,3
            D2(i1,i2) = (zadj(m,1)*S3mod(1,i1,i2) 
     &          + zadj(m,2)*S3mod(2,i1,i2) 
     &          + zadj(m,3)*S3mod(3,i1,i2) 
     &          - detz*D3(m,i1,i2))/zadjf(m)
c            write(testout,20) 'D2  = ',i1,i2,D2(i1,i2)
          end do
        end do

c     >      if (sym.eq.1) then
c     >        D2(2,1) = D2(1,2)
c     >        D2(3,1) = D2(1,3)
c     >        D2(3,2) = D2(2,3)
c     >      end if

        if (rank+ordg4-r.eq.2) goto 100

        
        S002(0,0) = 2d0*mm02*D2(0,0) + 2d0*C2m(0,0,0)

        do i1=1,3
          do i2=1,3
            S002(i1,i2) = 2d0*mm02*D2(i1,i2) + 2d0*C2m(0,i1,i2)
          end do
        end do
        
        D4(0,0,0,0) = zadj(k,l)*(S002(0,0) + 1d0/6d0) 
     &      - detz*D4(0,0,k,l)
        do a=1,3
          D4(0,0,0,0) = D4(0,0,0,0) - zadj(k,l)*S4hat(a,0,0,a) 
     &        + zadj(a,l)*S4hat(a,0,0,k)
          do b=1,3
            if(b.ne.k.and.a.ne.l) then
              D4(0,0,0,0) = D4(0,0,0,0) - zadj2(k,b,l,a)*( 
     &            + f(b)*S3hat(a,0,0)
     &            - f(a)*f(b)*D2(0,0)   )
            end if
          end do
        end do
        D4(0,0,0,0) = D4(0,0,0,0) /(8d0*zadj(k,l))
c        write(testout,10) 'D0000  = ',D4(0,0,0,0)
        
        do i1=1,3
          do i2=1,3
            D4(0,0,i1,i2) = zadj(k,l)*S002(i1,i2) 
     &          - detz*D4(k,l,i1,i2)
            do a=1,3
              D4(0,0,i1,i2) = D4(0,0,i1,i2) 
     &            - zadj(k,l)*S4hat(a,a,i1,i2) 
     &            + zadj(a,l)*S4hat(a,k,i1,i2)
              do b=1,3
                if(b.ne.k.and.a.ne.l) then
                  D4(0,0,i1,i2) = D4(0,0,i1,i2) - zadj2(k,b,l,a)*( 
     &                + f(b)*S3hat(a,i1,i2)
     &                - f(a)*f(b)*D2(i1,i2)   )
                end if
              end do
              if(i1.ne.k.and.a.ne.l) then
                D4(0,0,i1,i2) = D4(0,0,i1,i2) -
     &              2d0*zadj2(k,i1,l,a)*( 
     &              + S4hat(a,0,0,i2)
     &              - f(a)*D3(0,0,i2)  )
              end if
              if(i2.ne.k.and.a.ne.l) then
                D4(0,0,i1,i2) = D4(0,0,i1,i2)-2d0*zadj2(k,i2,l,a)*( 
     &              + S4hat(a,0,0,i1)
     &              - f(a)*D3(0,0,i1)  )
              end if
            end do
            do b=1,3
              if(b.ne.k.and.i1.ne.l) then
                D4(0,0,i1,i2) = D4(0,0,i1,i2)-2d0*zadj2(k,b,l,i1)*( 
     &              - f(b)*D3(0,0,i2)   )
              end if
              if(b.ne.k.and.i2.ne.l) then
                D4(0,0,i1,i2) = D4(0,0,i1,i2)-2d0*zadj2(k,b,l,i2)*( 
     &              - f(b)*D3(0,0,i1)   )
              end if
            end do
            if(i1.ne.k.and.i2.ne.l) then
              D4(0,0,i1,i2) = D4(0,0,i1,i2) - 4d0*zadj2(k,i1,l,i2)*( 
     &            - D4(0,0,0,0)   )
            end if
            if(i2.ne.k.and.i1.ne.l) then
              D4(0,0,i1,i2) = D4(0,0,i1,i2) - 4d0*zadj2(k,i2,l,i1)*( 
     &            - D4(0,0,0,0)   )
            end if
            D4(0,0,i1,i2) = D4(0,0,i1,i2) /(12d0*zadj(k,l))
c            write(testout,20) 'D00ij  = ',i1,i2,D4(0,0,i1,i2)
          end do
        end do
        
        do i1=1,3
          do i2=1,3
            do i3=1,3
              do j =1,3
                S4mod(j,i1,i2,i3) = S4hat(j,i1,i2,i3)
              end do
              S4mod(i1,i1,i2,i3) = S4mod(i1,i1,i2,i3) -2d0*D4(0,0,i2,i3)
              S4mod(i2,i1,i2,i3) = S4mod(i2,i1,i2,i3) -2d0*D4(0,0,i1,i3)
              S4mod(i3,i1,i2,i3) = S4mod(i3,i1,i2,i3) -2d0*D4(0,0,i1,i2)
            end do
          end do
        end do
        
        do i1=1,3
          do i2=1,3
            do i3=1,3
              D3(i1,i2,i3) = (zadj(m,1)*S4mod(1,i1,i2,i3) 
     &            + zadj(m,2)*S4mod(2,i1,i2,i3) 
     &            + zadj(m,3)*S4mod(3,i1,i2,i3) 
     &            - detz*D4(m,i1,i2,i3))/zadjf(m)
c              write(testout,30) 'D3  = ',i1,i2,i3,D3(i1,i2,i3)
            end do
          end do
        end do

c     >      if (sym.eq.1) then
c     >        D3(2,1,1) = D3(1,1,2)
c     >        D3(1,2,1) = D3(1,1,2)
c     >        D3(3,1,1) = D3(1,1,3)
c     >        D3(1,3,1) = D3(1,1,3)
c     >        D3(2,1,2) = D3(1,2,2)
c     >        D3(2,2,1) = D3(1,2,2)
c     >        D3(3,1,3) = D3(1,3,3)
c     >        D3(3,3,1) = D3(1,3,3)
c     >        D3(3,3,2) = D3(2,3,3)
c     >        D3(3,2,3) = D3(2,3,3)
c     >        D3(3,2,2) = D3(2,2,3)
c     >        D3(2,3,2) = D3(2,2,3)
c     >        D3(2,1,3) = D3(1,2,3)
c     >        D3(1,3,2) = D3(1,2,3)
c     >        D3(3,1,2) = D3(1,2,3)
c     >        D3(2,3,1) = D3(1,2,3)
c     >        D3(3,2,1) = D3(1,2,3)
c     >      end if

        if (rank+ordg4-r.eq.3) goto 100

        do i1=1,3
          S003(0,0,i1) = 2d0*mm02*D3(0,0,i1) + 2d0*C3m(0,0,0,i1)
          do i2=1,3
            do i3=1,3
              S003(i1,i2,i3) = 2d0*mm02*D3(i1,i2,i3) 
     &            + 2d0*C3m(0,i1,i2,i3)
            end do
          end do
        end do

        do i1=1,3
          D5(0,0,0,0,i1) = zadj(k,l)*(S003(0,0,i1) - 1d0/24d0) 
     &        - detz*D5(0,0,k,l,i1)
          do a=1,3
            D5(0,0,0,0,i1) = D5(0,0,0,0,i1) 
     &          - zadj(k,l)*S5hat(a,0,0,a,i1) 
     &          + zadj(a,l)*S5hat(a,0,0,k,i1)
            do b=1,3
              if(b.ne.k.and.a.ne.l) then
                D5(0,0,0,0,i1) = D5(0,0,0,0,i1) - zadj2(k,b,l,a)*( 
     &              + f(b)*S4hat(a,0,0,i1)
     &              - f(a)*f(b)*D3(0,0,i1)   )
              end if
            end do
            if(i1.ne.k.and.a.ne.l) then
              D5(0,0,0,0,i1) = D5(0,0,0,0,i1) 
     &            - 2d0*zadj2(k,i1,l,a)*( 
     &            + S5hat(a,0,0,0,0)
     &            - f(a)*D4(0,0,0,0)  )
            end if
          end do
          do b=1,3
            if(b.ne.k.and.i1.ne.l) then
              D5(0,0,0,0,i1) = D5(0,0,0,0,i1) 
     &            - 2d0*zadj2(k,b,l,i1)*( 
     &            - f(b)*D4(0,0,0,0)   )
            end if
          end do
          D5(0,0,0,0,i1) = D5(0,0,0,0,i1) /(12d0*zadj(k,l))
c          write(testout,10) 'D0000i  = ',i1,D5(0,0,0,0,i1)
        end do
      
        do i1=1,3
          do i2=1,3
            do i3=1,3
              D5(0,0,i1,i2,i3) = zadj(k,l)*S003(i1,i2,i3) 
c     &            - detz*D5(k,l,i1,i2,i3)
              do a=1,3
                D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &              - zadj(k,l)*S5hat(a,a,i1,i2,i3) 
     &              + zadj(a,l)*S5hat(a,k,i1,i2,i3)
                do b=1,3
                  if(b.ne.k.and.a.ne.l) then
                    D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &                  - zadj2(k,b,l,a)*( 
     &                  + f(b)*S4hat(a,i1,i2,i3)
     &                  - f(a)*f(b)*D3(i1,i2,i3)   )
                  end if
                end do
                if(i1.ne.k.and.a.ne.l) then
                  D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &                - 2d0*zadj2(k,i1,l,a)*( 
     &                + S5hat(a,0,0,i2,i3)
     &                - f(a)*D4(0,0,i2,i3)  )
                end if
                if(i2.ne.k.and.a.ne.l) then
                  D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &                - 2d0*zadj2(k,i2,l,a)*( 
     &                + S5hat(a,0,0,i1,i3)
     &                - f(a)*D4(0,0,i1,i3)  )
                end if
                if(i3.ne.k.and.a.ne.l) then
                  D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &                - 2d0*zadj2(k,i3,l,a)*( 
     &                + S5hat(a,0,0,i1,i2)
     &                - f(a)*D4(0,0,i1,i2)  )
                end if
              end do
              do b=1,3
                if(b.ne.k.and.i1.ne.l) then
                  D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &                - 2d0*zadj2(k,b,l,i1)*( 
     &                - f(b)*D4(0,0,i2,i3)   )
                end if
                if(b.ne.k.and.i2.ne.l) then
                  D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &                - 2d0*zadj2(k,b,l,i2)*( 
     &                - f(b)*D4(0,0,i1,i3)   )
                end if
                if(b.ne.k.and.i3.ne.l) then
                  D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &                - 2d0*zadj2(k,b,l,i3)*( 
     &                - f(b)*D4(0,0,i1,i2)   )
                end if
              end do
              if(i1.ne.k.and.i2.ne.l) then
                D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &              + 4d0*zadj2(k,i1,l,i2) * D5(0,0,0,0,i3)   
c     write(*,*) 'D5 31 ', 
c     &       + 4d0*zadj2(k,i1,l,i2) * D5(0,0,0,0,i3)   
              end if
              if(i2.ne.k.and.i1.ne.l) then
                D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &              + 4d0*zadj2(k,i2,l,i1) * D5(0,0,0,0,i3)   
c     write(*,*) 'D5 32 ', 
c     &       + 4d0*zadj2(k,i2,l,i1) * D5(0,0,0,0,i3),   
c     &       + 4d0*zadj2(k,i2,l,i1) , D5(0,0,0,0,i3)   
              end if
              if(i1.ne.k.and.i3.ne.l) then
                D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &              + 4d0*zadj2(k,i1,l,i3) * D5(0,0,0,0,i2)   
c     write(*,*) 'D5 33 ', 
c     &       + 4d0*zadj2(k,i1,l,i3) * D5(0,0,0,0,i2)   
              end if
              if(i3.ne.k.and.i1.ne.l) then
                D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &              + 4d0*zadj2(k,i3,l,i1) * D5(0,0,0,0,i2)   
c     write(*,*) 'D5 34 ', 
c     &       + 4d0*zadj2(k,i3,l,i1) * D5(0,0,0,0,i2),   
c     &       + 4d0*zadj2(k,i3,l,i1) , D5(0,0,0,0,i2)   
              end if
              if(i3.ne.k.and.i2.ne.l) then
                D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &              + 4d0*zadj2(k,i3,l,i2) * D5(0,0,0,0,i1)   
c     write(*,*) 'D5 35 ', 
c     &       + 4d0*zadj2(k,i3,l,i2) * D5(0,0,0,0,i1)   
              end if
              if(i2.ne.k.and.i3.ne.l) then
                D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) 
     &              + 4d0*zadj2(k,i2,l,i3) * D5(0,0,0,0,i1)   
c     write(*,*) 'D5 36 ', 
c     &       + 4d0*zadj2(k,i2,l,i3) * D5(0,0,0,0,i1)   
              end if
c     write(*,*) 'D5 4 ',D5(0,0,i1,i2,i3) 
              D5(0,0,i1,i2,i3) = D5(0,0,i1,i2,i3) /(16d0*zadj(k,l))
c              write(testout,30) 
c     &            'D00ijk  = ',i1,i2,i3,D5(0,0,i1,i2,i3)
c     write(*,*) 'D5 5 ',D5(0,0,i1,i2,i3) 
            end do
          end do
        end do

c     >      if (sym.eq.1) then
c     >        D5(0,0,2,1,1) = D5(0,0,1,1,2)
c     >        D5(0,0,1,2,1) = D5(0,0,1,1,2)
c     >        D5(0,0,3,1,1) = D5(0,0,1,1,3)
c     >        D5(0,0,1,3,1) = D5(0,0,1,1,3)
c     >        D5(0,0,2,1,2) = D5(0,0,1,2,2)
c     >        D5(0,0,2,2,1) = D5(0,0,1,2,2)
c     >        D5(0,0,3,1,3) = D5(0,0,1,3,3)
c     >        D5(0,0,3,3,1) = D5(0,0,1,3,3)
c     >        D5(0,0,3,3,2) = D5(0,0,2,3,3)
c     >        D5(0,0,3,2,3) = D5(0,0,2,3,3)
c     >        D5(0,0,3,2,2) = D5(0,0,2,2,3)
c     >        D5(0,0,2,3,2) = D5(0,0,2,2,3)
c     >        D5(0,0,2,1,3) = D5(0,0,1,2,3)
c     >        D5(0,0,1,3,2) = D5(0,0,1,2,3)
c     >        D5(0,0,3,1,2) = D5(0,0,1,2,3)
c     >        D5(0,0,2,3,1) = D5(0,0,1,2,3)
c     >        D5(0,0,3,2,1) = D5(0,0,1,2,3)
c     >      end if

        do i1=1,3
          do i2=1,3
            do i3=1,3
              do i4=1,3
                do j =1,3
                  S5mod(j,i1,i2,i3,i4) = S5hat(j,i1,i2,i3,i4)
                end do
                S5mod(i1,i1,i2,i3,i4) = S5mod(i1,i1,i2,i3,i4) 
     &              - 2d0*D5(0,0,i2,i3,i4)
                S5mod(i2,i1,i2,i3,i4) = S5mod(i2,i1,i2,i3,i4) 
     &              - 2d0*D5(0,0,i1,i3,i4)
                S5mod(i3,i1,i2,i3,i4) = S5mod(i3,i1,i2,i3,i4) 
     &              - 2d0*D5(0,0,i1,i2,i4)
                S5mod(i4,i1,i2,i3,i4) = S5mod(i4,i1,i2,i3,i4) 
     &              - 2d0*D5(0,0,i1,i2,i3)
              end do
            end do
          end do
        end do

        do i1=1,3
          do i2=1,3
            do i3=1,3
              do i4=1,3
                  D4(i1,i2,i3,i4) = (zadj(m,1)*S5mod(1,i1,i2,i3,i4) 
     &                + zadj(m,2)*S5mod(2,i1,i2,i3,i4) 
     &                + zadj(m,3)*S5mod(3,i1,i2,i3,i4) 
c     &                - detz*D5(m,i1,i2,i3,i4)
     &              )/zadjf(m)
c                  write(testout,40) 'D4  = ',i1,i2,i3,i4,D4(i1,i2,i3,i4)
              end do
            end do
          end do
        end do

 100  continue

      end

***********************************************************************
      subroutine cDgy12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank)
***********************************************************************
*     4-point tensor coefficient functions D0,D1,D2,D3,D4,D5          *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     expansion about Gram determinant = 0                            *
*               and Cayley determinant = 0 to order ordgy4            *
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
*     i-th denominator cancelled: Cxwi                                *
***********************************************************************
*     22.12.04 Ansgar Denner    last changed  27.09.05                *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p32,p30,p20,p31
      complex*16 q10,q21,q32,q30,q20,q31
      complex*16 m02,m12,m22,m32
      complex*16 mm02,mm12,mm22,mm32,f(3),zadjf(3)
      complex*16 k1k2,k1k3,k2k3,z(3,3),zadj(3,3),detz
      complex*16 zadj2(3,3,3,3),z2f(3,3,3),z2ff(3,3)
      complex*16 C0w3,C1w3(2),C2w3(0:2,0:2),C3w3(0:2,0:2,0:2)
     &    ,C4w3(0:2,0:2,0:2,0:2),C5w3(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w3(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w2,C1w2(2),C2w2(0:2,0:2),C3w2(0:2,0:2,0:2)
     &    ,C4w2(0:2,0:2,0:2,0:2),C5w2(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w2(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w1,C1w1(2),C2w1(0:2,0:2),C3w1(0:2,0:2,0:2)
     &    ,C4w1(0:2,0:2,0:2,0:2),C5w1(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w1(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w0,C1w0(2),C2w0(0:2,0:2),C3w0(0:2,0:2,0:2)
     &    ,C4w0(0:2,0:2,0:2,0:2),C5w0(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w0(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0m(0:3),C1m(0:3,3),C2m(0:3,0:3,0:3)
     &    ,C3m(0:3,0:3,0:3,3),C4m(0:3,0:3,0:3,0:3,0:3)
     &    ,C5m(0:3,0:3,0:3,0:3,0:3,3)
c     &    ,C6m(0:3,0:3,0:3,0:3,0:3,0:3,0:3)
      complex*16 S1hat(3),S2hat(3,3),S3hat(3,0:3,0:3),
     &    S4hat(3,0:3,0:3,3),S5hat(0:3,0:3,0:3,0:3,0:3),
     &    S6hat(3,0:3,0:3,0:3,0:3,3)
c     &    S6hat(3,0:3,0:3,0:3,0:3,3),S7hat(0:3,0:3,0:3,0:3,0:3,0:3,0:3)
      complex*16 S2mod(3,3),S3mod(3,0:3,0:3),
     &    S4mod(3,0:3,0:3,3),S5mod(0:3,0:3,0:3,0:3,0:3)
      complex*16 S002(0:3,0:3)
c      complex*16 S001(3),S002(0:3,0:3),S003(0:3,0:3,0:3)
      complex*16 R3(3,3,3),R4(3,3,3,3),R5(3,3,3,3,3),R6(3,3,3,3,3,3)
      complex*16 D0,D1(3),D2(0:3,0:3),D3(0:3,0:3,0:3)
      complex*16 D4(0:3,0:3,0:3,0:3),D5(0:3,0:3,0:3,0:3,0:3)
      complex*16 D6(0:3,0:3,0:3,0:3,0:3,0:3)
c      complex*16 D7(0:3,0:3,0:3,0:3,0:3,0:3,0:3)
c     complex*16 cD0f
      complex*16 elimcminf2
      real*8     maxzadj,maxz2ff
      integer    rank,rankord,flag,rankC
      integer    i1,i2,i3,i4,i5,i,j,k,l,l2,l3,m,r
      integer    ltest,sym
      integer    a,b

      integer    ordg3,ordgy3,ordgp3
      integer    ordg4,ordgy4,ordgp4
      common /ctord3/ ordg3,ordgy3,ordgp3
      common /ctord4/ ordg4,ordgy4,ordgp4
      common /ltest/  ltest
      common /sym/ sym

      data       C1m /12*0d0/,C2m /64*0d0/,C3m/192*0d0/,C4m/1024*0d0/
      data       C5m /3072*0d0/
c      data       C6m /16384*0d0/
      data       zadj2/81*0d0/
      data       flag /0/
c      integer testout
c      data       testout /49/
      data       i,k,l,m /4*0/

c      write(*,*)'Dgy12345 in',p10,p21,p32,p30,p20,p31,
c     &    m02,m12,m22,m32
c     &    ,rank,ordgy4
c      write(testout,*)'Dgy12345 in',p10,p21,p32,p30,p20,p31,
c     &    m02,m12,m22,m32
c     &    ,rank,ordgy4

      rankord = rank + 2*ordgy4

      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cDgy12345'
        write(*,*) 'rank = ',rank
        stop
      else if (rankord.gt.4.and.flag.eq.0) then
        write(*,*) 'rank+2*ordgy4 > 4 not implemented in cDgy12345'
        write(*,*) 'rank   = ',rank
        write(*,*) 'ordgy4 = ',ordgy4
        flag = 1
        if (rank.gt.4) then
          write(*,*) 'rank > 4 not implemented in cDgy12345'
          write(*,*) 'rank = ',rank
          stop
        end if
      end if

c      write(*,*) ' cDgy12345: ordgy4 = ',ordgy4
c      write(*,*) ' cDgy12345: rankord = ',rankord

c--->   27.09.05    "rank+2*ordgy4" replaced by "rankC"
c                   in calls of  cCp12345
      rankC = min(rankord+1,5)
c     write(*,*) 'cDp1',p21,p32,p31,m12,m22,m32,rankC
      call cCp12345(p21,p32,p31,m12,m22,m32,C0w0,C1w0,C2w0,C3w0,C4w0,
     &    C5w0,C6w0,rankC)
c     write(*,*) 'cDp2',p21,p32,p31,m12,m22,m32,rankC
      call cCp12345(p20,p32,p30,m02,m22,m32,C0w1,C1w1,C2w1,C3w1,C4w1,
     &    C5w1,C6w1,rankC)
c     write(*,*) 'cDp3',p21,p32,p31,m12,m22,m32,rankC
      call cCp12345(p10,p31,p30,m02,m12,m32,C0w2,C1w2,C2w2,C3w2,C4w2,
     &    C5w2,C6w2,rankC)
c     write(*,*) 'cDp4',p21,p32,p31,m12,m22,m32,rankC
      call cCp12345(p10,p21,p20,m02,m12,m22,C0w3,C1w3,C2w3,C3w3,C4w3,
     &    C5w3,C6w3,rankC)
c <----      

c      write(*,*) 'cDgy C0w0',C0w0,C1w0
c      write(*,*) 'cDgy C0w0',C0w1,C1w1
c      write(*,*) 'cDgy C0w0',C0w2,C1w2
c      write(*,*) 'cDgy C0w0',C0w3,C1w3

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

      z(1,1) = 2d0*q10
      z(2,2) = 2d0*q20
      z(3,3) = 2d0*q30
      z(1,2) = q10+q20-q21
      z(1,3) = q10+q30-q31
      z(2,3) = q20+q30-q32
      z(2,1) = z(1,2)
      z(3,1) = z(1,3)
      z(3,2) = z(2,3)
      
      k1k2 = (q10+q20-q21)
      k1k3 = (q10+q30-q31)
      k2k3 = (q20+q30-q32)
      detz = 8d0*q10*q30*q20+2D0*k1k2*k1k3*k2k3
     &    -2d0*(q10*k2k3*k2k3+q20*k1k3*k1k3+q30*k1k2*k1k2)

c      write(testout,*) 'detz = ',detz,detz/
c     &    dmax1(abs(k1k2),abs(k1k3),abs(k2k3),
c     &    2*abs(q10),2*abs(q20),2*abs(q30))**3

      zadj(1,1) = (4d0*q30*q20-k2k3*k2k3)
      zadj(1,2) = (k1k3*k2k3-2d0*q30*k1k2)
      zadj(1,3) = (k1k2*k2k3-2d0*q20*k1k3)
      zadj(2,1) = zadj(1,2)
      zadj(2,2) = (4d0*q10*q30-k1k3*k1k3)
      zadj(2,3) = (k1k2*k1k3-2d0*q10*k2k3)
      zadj(3,1) = zadj(1,3)
      zadj(3,2) = zadj(2,3)
      zadj(3,3) = (4d0*q10*q20-k1k2*k1k2)

      zadj2(1,2,1,2) = -z(3,3)
      zadj2(1,2,2,1) = z(3,3)
      zadj2(2,1,1,2) = z(3,3)
      zadj2(2,1,2,1) = -z(3,3)
      zadj2(1,3,1,3) = -z(2,2)
      zadj2(1,3,3,1) = z(2,2)
      zadj2(3,1,1,3) = z(2,2)
      zadj2(3,1,3,1) = -z(2,2)
      zadj2(3,2,3,2) = -z(1,1)
      zadj2(3,2,2,3) = z(1,1)
      zadj2(2,3,3,2) = z(1,1)
      zadj2(2,3,2,3) = -z(1,1)
      zadj2(1,2,1,3) = z(3,2)
      zadj2(1,2,3,1) = -z(3,2)
      zadj2(2,1,1,3) = -z(3,2)
      zadj2(2,1,3,1) = z(3,2)
      zadj2(1,3,1,2) = z(2,3)
      zadj2(1,3,2,1) = -z(2,3)
      zadj2(3,1,1,2) = -z(2,3)
      zadj2(3,1,2,1) = z(2,3)
      zadj2(3,2,1,2) = z(1,3)
      zadj2(3,2,2,1) = -z(1,3)
      zadj2(2,3,1,2) = -z(1,3)
      zadj2(2,3,2,1) = z(1,3)
      zadj2(1,2,3,2) = z(3,1)
      zadj2(1,2,2,3) = -z(3,1)
      zadj2(2,1,3,2) = -z(3,1)
      zadj2(2,1,2,3) = z(3,1)
      zadj2(1,3,2,3) = z(2,1)
      zadj2(1,3,3,2) = -z(2,1)
      zadj2(3,1,2,3) = -z(2,1)
      zadj2(3,1,3,2) = z(2,1)
      zadj2(3,2,3,1) = z(1,2)
      zadj2(3,2,1,3) = -z(1,2)
      zadj2(2,3,3,1) = -z(1,2)
      zadj2(2,3,1,3) = z(1,2)

      f(1) = q10+mm02-mm12
      f(2) = q20+mm02-mm22
      f(3) = q30+mm02-mm32
      
      zadjf(1) = zadj(1,1)*f(1)+zadj(1,2)*f(2)+zadj(1,3)*f(3)
      zadjf(2) = zadj(2,1)*f(1)+zadj(2,2)*f(2)+zadj(2,3)*f(3)
      zadjf(3) = zadj(3,1)*f(1)+zadj(3,2)*f(2)+zadj(3,3)*f(3)

      do i1 = 1,3
      do i2 = 1,3
        z2ff(i1,i2) = 2d0*mm02*zadj(i1,i2)
      do i3 = 1,3
        z2f(i1,i2,i3) = 0d0
      do i4 = 1,3
        z2ff(i1,i2) = z2ff(i1,i2) + zadj2(i1,i3,i2,i4)*f(i3)*f(i4)
        z2f(i1,i2,i3) = z2f(i1,i2,i3) + zadj2(i1,i4,i2,i3)*f(i4)
      end do
      end do
      end do
      end do

c      write(testout,*) 'zadjf ',zadjf
c      write(testout,*) 'zadjf ',detz/zadjf(1),detz/zadjf(2),detz/zadjf(3)

      do i1=1,3
        D1(i1) = 0d0
        D3(0,0,i1) = 0d0
        D5(0,0,0,0,i1) = 0d0
        do i2=1,3
          D2(i1,i2) = 0d0
          D4(0,0,i1,i2) = 0d0
          D6(0,0,0,0,i1,i2) = 0d0
          do i3=1,3
            D3(i1,i2,i3) = 0d0
            D5(0,0,i1,i2,i3) = 0d0
            do i4=1,3
              D4(i1,i2,i3,i4) = 0d0
              D6(0,0,i1,i2,i3,i4) = 0d0
c              do i5=1,3
c                D5(i1,i2,i3,i4,i5) = 0d0
c                do i6=1,3
c                  D6(i1,i2,i3,i4,i5,i6) = 0d0
c                end do
c              end do
            end do
          end do
        end do
      end do

      C0m(0) = C0w0
      C0m(1) = C0w1
      C0m(2) = C0w2
      C0m(3) = C0w3
      
      do j =1,3
        S1hat(j) = C0m(j) - C0m(0)
      end do

      C1m(0,2) = C1w0(1)
      C1m(0,3) = C1w0(2)
      C1m(0,1) = -C0m(0) - C1m(0,2) - C1m(0,3)
      C1m(1,2) = C1w1(1)
      C1m(1,3) = C1w1(2)
      C1m(2,1) = C1w2(1)
      C1m(2,3) = C1w2(2)
      C1m(3,1) = C1w3(1)
      C1m(3,2) = C1w3(2)
      
      do i1=1,3
        do j =1,3
          S2hat(j,i1) = C1m(j,i1) - C1m(0,i1)
c          write(testout,*) 'C1m ',j,i1,C1m(j,i1), - C1m(0,i1)
        end do
      end do
      

      if (rankord.eq.0) goto 99
  
c     >      do i1=1,3
c     >      do i2=1,3
c     >      do k=1,3
c     >        do l=1,3
c     >          test = zadj(k,l)*z(i1,i2)
c     >          do a=1,3
c     >            do b=1,3
c     >              if(b.ne.k.and.a.ne.l) then
c     >                test  = test + zadj2(k,b,l,a)*z(b,i1)*z(a,i2)
c     >              end if
c     >            end do
c     >          end do
c     >          write(testout,*) 'detz = ',i1,i2,k,l,test
c     >        end do
c     >      end do
c     >      end do
c     >      end do

        C2m(0,0,0) = C2w0(0,0)
        C2m(1,0,0) = C2w1(0,0)
        C2m(2,0,0) = C2w2(0,0)
        C2m(3,0,0) = C2w3(0,0)
        
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
        
        C2m(0,3,2) = C2m(0,2,3) 
        C2m(0,3,1) = C2m(0,1,3) 
        C2m(0,2,1) = C2m(0,1,2) 
        C2m(1,3,2) = C2m(1,2,3) 
        C2m(2,3,1) = C2m(2,1,3) 
        C2m(3,2,1) = C2m(3,1,2) 

        do i1=1,3
          do i2=1,3
            do j =1,3
              S3hat(j,i1,i2) = C2m(j,i1,i2) - C2m(0,i1,i2)
            end do
          end do
        end do

        if (rankord.eq.1) goto 99

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
        
        C3m(0,1,2,1) = C3m(0,1,1,2)
        C3m(0,2,1,1) = C3m(0,1,1,2)
        C3m(0,1,3,1) = C3m(0,1,1,3)
        C3m(0,3,1,1) = C3m(0,1,1,3)
        C3m(0,2,2,1) = C3m(0,1,2,2)
        C3m(0,2,1,2) = C3m(0,1,2,2)
        C3m(0,3,3,1) = C3m(0,1,3,3)
        C3m(0,3,1,3) = C3m(0,1,3,3)
        C3m(0,2,3,2) = C3m(0,2,2,3)
        C3m(0,3,2,2) = C3m(0,2,2,3)
        C3m(0,3,3,2) = C3m(0,2,3,3)
        C3m(0,3,2,3) = C3m(0,2,3,3)
        
        C3m(0,1,3,2) = C3m(0,1,2,3)
        C3m(0,2,1,3) = C3m(0,1,2,3)
        C3m(0,2,3,1) = C3m(0,1,2,3)
        C3m(0,3,1,2) = C3m(0,1,2,3)
        C3m(0,3,2,1) = C3m(0,1,2,3)
        
        C3m(3,1,2,1) = C3m(3,1,1,2)
        C3m(3,2,1,1) = C3m(3,1,1,2)
        C3m(2,1,3,1) = C3m(2,1,1,3)
        C3m(2,3,1,1) = C3m(2,1,1,3)
        C3m(3,2,2,1) = C3m(3,1,2,2)
        C3m(3,2,1,2) = C3m(3,1,2,2)
        C3m(2,3,3,1) = C3m(2,1,3,3)
        C3m(2,3,1,3) = C3m(2,1,3,3)
        C3m(1,2,3,2) = C3m(1,2,2,3)
        C3m(1,3,2,2) = C3m(1,2,2,3)
        C3m(1,3,3,2) = C3m(1,2,3,3)
        C3m(1,3,2,3) = C3m(1,2,3,3)
        
        do i1=1,3
          do i2=1,3
            S4hat(i1,0,0,i2) = C3m(i1,0,0,i2) - C3m(0,0,0,i2)
            do i3=1,3
              do j =1,3
                S4hat(j,i1,i2,i3) = C3m(j,i1,i2,i3) - C3m(0,i1,i2,i3)
              end do
            end do
          end do
        end do
        
        if (rankord.eq.2) goto 99

        do i1=1,3
          S3hat(i1,0,0) = C2m(i1,0,0) - C2m(0,0,0)
        end do

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

c     >      do i1=1,2
c     >      do i2=i1,2
c     >        C4m(0,0,0,i1+1,i2+1) = C4w0(0,0,i1,i2)
c     >        C4m(1,0,0,i1+1,i2+1) = C4w0(0,0,i1,i2)
c     >        C4m(2,0,0,2*i1-1,2*i2-1) = C4w0(0,0,i1,i2)
c     >        C4m(3,0,0,i1+1,i2) = C4w0(0,0,i1,i2)
c     >      do i3=1,2
c     >      do i4=1,2
c     >        C4m(0,i1+1,i2+1,i3+1,i4+1) = C4w0(i1,i2,i3,i4)
c     >        C4m(1,i1+1,i2+1,i3+1,i4+1) = C4w0(i1,i2,i3,i4)
c     >        C4m(2,2*i1-1,2*i2-1,2*i3-1,2*i4-1) = C4w0(i1,i2,i3,i4)
c     >        C4m(3,i1+1,i2,i3,i4) = C4w0(i1,i2,i3,i4)
c     >      end do
c     >        C4m(0,1,i1+1,i2+1,i3+1) = -C3m(0,i1+1,i2+1,i3+1) 
c     >     &      - C4m(0,2,i1+1,i2+1,i3+1) 
c     >     &      - C4m(0,i1+1,i2+1,i3+1,3)
c     >      end do
c     >        C4m(0,1,1,i1+1,i2+1) = -C3m(0,1,i1+1,i2+1) 
c     >     &      - C4m(0,1,2,i1+1,i2+1) 
c     >     &      - C4m(0,1,i1+1,i2+1,3)
c     >      end do
c     >        C4m(0,1,1,1,i1+1) = -C3m(0,1,1,i1+1) 
c     >     &      - C4m(0,1,1,2,i1+1) 
c     >     &      - C4m(0,1,1,i1+1,3)
c     >      end do
c     >      C4m(0,1,1,1,1) = -C3m(0,1,1,1) 
c     >     &      - C4m(0,1,1,1,2) 
c     >     &      - C4m(0,1,1,1,3)

      do j=0,3
      do i1=1,3
      do i2=i1,3
        C4m(j,0,0,i2,i1) = C4m(j,0,0,i1,i2)
      do i3=i2,3
      do i4=i3,3
        C4m(j,i1,i3,i2,i4) = C4m(j,i1,i2,i3,i4)
        C4m(j,i2,i1,i3,i4) = C4m(j,i1,i2,i3,i4)
        C4m(j,i2,i3,i1,i4) = C4m(j,i1,i2,i3,i4)
        C4m(j,i3,i1,i2,i4) = C4m(j,i1,i2,i3,i4)
        C4m(j,i3,i2,i1,i4) = C4m(j,i1,i2,i3,i4)
        C4m(j,i1,i2,i4,i3) = C4m(j,i1,i2,i3,i4)
        C4m(j,i1,i4,i2,i3) = C4m(j,i1,i2,i3,i4)
        C4m(j,i2,i1,i4,i3) = C4m(j,i1,i2,i3,i4)
        C4m(j,i2,i4,i1,i3) = C4m(j,i1,i2,i3,i4)
        C4m(j,i4,i1,i2,i3) = C4m(j,i1,i2,i3,i4)
        C4m(j,i4,i2,i1,i3) = C4m(j,i1,i2,i3,i4)
        C4m(j,i1,i3,i4,i2) = C4m(j,i1,i2,i3,i4)
        C4m(j,i1,i4,i3,i2) = C4m(j,i1,i2,i3,i4)
        C4m(j,i3,i1,i4,i2) = C4m(j,i1,i2,i3,i4)
        C4m(j,i3,i4,i1,i2) = C4m(j,i1,i2,i3,i4)
        C4m(j,i4,i1,i3,i2) = C4m(j,i1,i2,i3,i4)
        C4m(j,i4,i3,i1,i2) = C4m(j,i1,i2,i3,i4)
        C4m(j,i3,i2,i4,i1) = C4m(j,i1,i2,i3,i4)
        C4m(j,i3,i4,i2,i1) = C4m(j,i1,i2,i3,i4)
        C4m(j,i2,i3,i4,i1) = C4m(j,i1,i2,i3,i4)
        C4m(j,i2,i4,i3,i1) = C4m(j,i1,i2,i3,i4)
        C4m(j,i4,i3,i2,i1) = C4m(j,i1,i2,i3,i4)
        C4m(j,i4,i2,i3,i1) = C4m(j,i1,i2,i3,i4)
      end do
      end do
      end do
      end do
      end do

      do i1=1,3
c          S5hat(i1,0,0,0,0) = C4m(i1,0,0,0,0) - C4m(0,0,0,0,0)
      do i2=1,3
      do i3=1,3
c              S5hat(i1,0,0,i2,i3) = C4m(i1,0,0,i2,i3) - C4m(0,0,0,i2,i3)
      do i4=1,3
        do k =1,3
          S5hat(k,i1,i2,i3,i4) = C4m(k,i1,i2,i3,i4) -
     &        C4m(0,i1,i2,i3,i4)
        end do
      end do
      end do
      end do
      end do

c      write(*,*) 'Dgy ',S5hat(1,2,2,2,2), C4m(1,2,2,2,2), -
c     &        C4m(0,2,2,2,2)
c      write(*,*) 'Dgy ',S5hat(1,2,2,2,2), C1m(1,2),
c     &        C1m(1,1)

      if (rankord.eq.3) goto 99

      do i1=1,3
      do i2=1,3
          S4hat(i1,0,0,i2) = C3m(i1,0,0,i2) - C3m(0,0,0,i2)
      end do
      end do

      C5m(0,0,0,0,0,2) = C5w0(0,0,0,0,1)
      C5m(0,0,0,0,0,3) = C5w0(0,0,0,0,2)
      C5m(0,0,0,0,0,1) = -C4m(0,0,0,0,0) - C5m(0,0,0,0,0,2) 
     &                                   - C5m(0,0,0,0,0,3)
      C5m(1,0,0,0,0,2) = C5w1(0,0,0,0,1)
      C5m(1,0,0,0,0,3) = C5w1(0,0,0,0,2)
      C5m(2,0,0,0,0,1) = C5w2(0,0,0,0,1)
      C5m(2,0,0,0,0,3) = C5w2(0,0,0,0,2)
      C5m(3,0,0,0,0,1) = C5w3(0,0,0,0,1)
      C5m(3,0,0,0,0,2) = C5w3(0,0,0,0,2)

      C5m(0,0,0,2,2,2) = C5w0(0,0,1,1,1)
      C5m(0,0,0,2,2,3) = C5w0(0,0,1,1,2)
      C5m(0,0,0,2,3,3) = C5w0(0,0,1,2,2)
      C5m(0,0,0,3,3,3) = C5w0(0,0,2,2,2)
      C5m(0,0,0,1,2,2) = -C4m(0,0,0,2,2) - C5m(0,0,0,2,2,2)
     &                                   - C5m(0,0,0,2,2,3)
      C5m(0,0,0,1,2,3) = -C4m(0,0,0,2,3) - C5m(0,0,0,2,2,3)
     &                                   - C5m(0,0,0,2,3,3)
      C5m(0,0,0,1,3,3) = -C4m(0,0,0,3,3) - C5m(0,0,0,2,3,3)
     &                                   - C5m(0,0,0,3,3,3)
      C5m(0,0,0,1,1,2) = -C4m(0,0,0,1,2) - C5m(0,0,0,1,2,2)
     &                                   - C5m(0,0,0,1,2,3)
      C5m(0,0,0,1,1,3) = -C4m(0,0,0,1,3) - C5m(0,0,0,1,2,3)
     &                                   - C5m(0,0,0,1,3,3)
      C5m(0,0,0,1,1,1) = -C4m(0,0,0,1,1) - C5m(0,0,0,1,1,2)
     &                                   - C5m(0,0,0,1,1,3)

      C5m(1,0,0,2,2,2) = C5w1(0,0,1,1,1)
      C5m(1,0,0,2,2,3) = C5w1(0,0,1,1,2)
      C5m(1,0,0,2,3,3) = C5w1(0,0,1,2,2)
      C5m(1,0,0,3,3,3) = C5w1(0,0,2,2,2)
      C5m(2,0,0,1,1,1) = C5w2(0,0,1,1,1)
      C5m(2,0,0,1,1,3) = C5w2(0,0,1,1,2)
      C5m(2,0,0,1,3,3) = C5w2(0,0,1,2,2)
      C5m(2,0,0,3,3,3) = C5w2(0,0,2,2,2)
      C5m(3,0,0,1,1,1) = C5w3(0,0,1,1,1)
      C5m(3,0,0,1,1,2) = C5w3(0,0,1,1,2)
      C5m(3,0,0,1,2,2) = C5w3(0,0,1,2,2)
      C5m(3,0,0,2,2,2) = C5w3(0,0,2,2,2)

      C5m(0,2,2,2,2,2) = C5w0(1,1,1,1,1)
      C5m(0,2,2,2,2,3) = C5w0(1,1,1,1,2)
      C5m(0,2,2,2,3,3) = C5w0(1,1,1,2,2)
      C5m(0,2,2,3,3,3) = C5w0(1,1,2,2,2)
      C5m(0,2,3,3,3,3) = C5w0(1,2,2,2,2)
      C5m(0,3,3,3,3,3) = C5w0(2,2,2,2,2)

      C5m(0,1,2,2,2,2) = -C4m(0,2,2,2,2) - C5m(0,2,2,2,2,2) 
     &                                   - C5m(0,2,2,2,2,3)
      C5m(0,1,2,2,2,3) = -C4m(0,2,2,2,3) - C5m(0,2,2,2,2,3) 
     &                                   - C5m(0,2,2,2,3,3)
      C5m(0,1,2,2,3,3) = -C4m(0,2,2,3,3) - C5m(0,2,2,2,3,3) 
     &                                   - C5m(0,2,2,3,3,3)
      C5m(0,1,2,3,3,3) = -C4m(0,2,3,3,3) - C5m(0,2,2,3,3,3) 
     &                                   - C5m(0,2,3,3,3,3)
      C5m(0,1,3,3,3,3) = -C4m(0,3,3,3,3) - C5m(0,2,3,3,3,3) 
     &                                   - C5m(0,3,3,3,3,3)
      C5m(0,1,1,2,2,2) = -C4m(0,1,2,2,2) - C5m(0,1,2,2,2,2) 
     &                                   - C5m(0,1,2,2,2,3)
      C5m(0,1,1,2,2,3) = -C4m(0,1,2,2,3) - C5m(0,1,2,2,2,3) 
     &                                   - C5m(0,1,2,2,3,3)
      C5m(0,1,1,2,3,3) = -C4m(0,1,2,3,3) - C5m(0,1,2,2,3,3) 
     &                                   - C5m(0,1,2,3,3,3)
      C5m(0,1,1,3,3,3) = -C4m(0,1,3,3,3) - C5m(0,1,2,3,3,3) 
     &                                   - C5m(0,1,3,3,3,3)
      C5m(0,1,1,1,2,2) = -C4m(0,1,1,2,2) - C5m(0,1,1,2,2,2) 
     &                                   - C5m(0,1,1,2,2,3)
      C5m(0,1,1,1,2,3) = -C4m(0,1,1,2,3) - C5m(0,1,1,2,2,3) 
     &                                   - C5m(0,1,1,2,3,3)
      C5m(0,1,1,1,3,3) = -C4m(0,1,1,3,3) - C5m(0,1,1,2,3,3) 
     &                                   - C5m(0,1,1,3,3,3)
      C5m(0,1,1,1,1,2) = -C4m(0,1,1,1,2) - C5m(0,1,1,1,2,2) 
     &                                   - C5m(0,1,1,1,2,3)
      C5m(0,1,1,1,1,3) = -C4m(0,1,1,1,3) - C5m(0,1,1,1,2,3)
     &                                   - C5m(0,1,1,1,3,3)
      C5m(0,1,1,1,1,1) = -C4m(0,1,1,1,1) - C5m(0,1,1,1,1,2) 
     &                                   - C5m(0,1,1,1,1,3)

      C5m(1,2,2,2,2,2) = C5w1(1,1,1,1,1)
      C5m(1,2,2,2,2,3) = C5w1(1,1,1,1,2)
      C5m(1,2,2,2,3,3) = C5w1(1,1,1,2,2)
      C5m(1,2,2,3,3,3) = C5w1(1,1,2,2,2)
      C5m(1,2,3,3,3,3) = C5w1(1,2,2,2,2)
      C5m(1,3,3,3,3,3) = C5w1(2,2,2,2,2)
      C5m(2,1,1,1,1,1) = C5w2(1,1,1,1,1)
      C5m(2,1,1,1,1,3) = C5w2(1,1,1,1,2)
      C5m(2,1,1,1,3,3) = C5w2(1,1,1,2,2)
      C5m(2,1,1,3,3,3) = C5w2(1,1,2,2,2)
      C5m(2,1,3,3,3,3) = C5w2(1,2,2,2,2)
      C5m(2,3,3,3,3,3) = C5w2(2,2,2,2,2)
      C5m(3,1,1,1,1,1) = C5w3(1,1,1,1,1)
      C5m(3,1,1,1,1,2) = C5w3(1,1,1,1,2)
      C5m(3,1,1,1,2,2) = C5w3(1,1,1,2,2)
      C5m(3,1,1,2,2,2) = C5w3(1,1,2,2,2)
      C5m(3,1,2,2,2,2) = C5w3(1,2,2,2,2)
      C5m(3,2,2,2,2,2) = C5w3(2,2,2,2,2)

c     >      do i1=1,2
c     >      do i2=i1,2
c     >        C5m(0,0,0,i1+1,i2+1) = C5w0(0,0,i1,i2)
c     >        C5m(1,0,0,i1+1,i2+1) = C5w0(0,0,i1,i2)
c     >        C5m(2,0,0,2*i1-1,2*i2-1) = C5w0(0,0,i1,i2)
c     >        C5m(3,0,0,i1+1,i2) = C5w0(0,0,i1,i2)
c     >      do i3=1,2
c     >      do i4=1,2
c     >        C5m(0,i1+1,i2+1,i3+1,i4+1) = C5w0(i1,i2,i3,i4)
c     >        C5m(1,i1+1,i2+1,i3+1,i4+1) = C5w0(i1,i2,i3,i4)
c     >        C5m(2,2*i1-1,2*i2-1,2*i3-1,2*i4-1) = C5w0(i1,i2,i3,i4)
c     >        C5m(3,i1+1,i2,i3,i4) = C5w0(i1,i2,i3,i4)
c     >      end do
c     >        C5m(0,1,i1+1,i2+1,i3+1) = -C3m(0,i1+1,i2+1,i3+1) 
c     >     &      - C5m(0,2,i1+1,i2+1,i3+1) 
c     >     &      - C5m(0,i1+1,i2+1,i3+1,3)
c     >      end do
c     >        C5m(0,1,1,i1+1,i2+1) = -C3m(0,1,i1+1,i2+1) 
c     >     &      - C5m(0,1,2,i1+1,i2+1) 
c     >     &      - C5m(0,1,i1+1,i2+1,3)
c     >      end do
c     >        C5m(0,1,1,1,i1+1) = -C3m(0,1,1,i1+1) 
c     >     &      - C5m(0,1,1,2,i1+1) 
c     >     &      - C5m(0,1,1,i1+1,3)
c     >      end do
c     >      C5m(0,1,1,1,1) = -C3m(0,1,1,1) 
c     >     &      - C5m(0,1,1,1,2) 
c     >     &      - C5m(0,1,1,1,3)

      do j=0,3
      do i1=1,3
      do i2=i1,3
      do i3=i2,3
        C5m(j,0,0,i1,i3,i2) = C5m(j,0,0,i1,i2,i3)
        C5m(j,0,0,i2,i1,i3) = C5m(j,0,0,i1,i2,i3)
        C5m(j,0,0,i2,i3,i1) = C5m(j,0,0,i1,i2,i3)
        C5m(j,0,0,i3,i1,i2) = C5m(j,0,0,i1,i2,i3)
        C5m(j,0,0,i3,i2,i1) = C5m(j,0,0,i1,i2,i3)
      do i4=i3,3
      do i5=i4,3
        C5m(j,i1,i3,i2,i4,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i1,i3,i4,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i3,i1,i4,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i1,i2,i4,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i2,i1,i4,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i2,i4,i3,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i4,i2,i3,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i1,i4,i3,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i4,i1,i3,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i1,i2,i3,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i2,i1,i3,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i3,i4,i2,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i4,i3,i2,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i1,i4,i2,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i4,i1,i2,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i1,i3,i2,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i3,i1,i2,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i2,i4,i1,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i4,i2,i1,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i3,i4,i1,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i4,i3,i1,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i3,i2,i1,i5) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i2,i3,i1,i5) = C5m(j,i1,i2,i3,i4,i5)

        C5m(j,i1,i2,i3,i5,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i3,i2,i5,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i1,i3,i5,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i3,i1,i5,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i1,i2,i5,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i2,i1,i5,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i2,i5,i3,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i5,i2,i3,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i1,i5,i3,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i5,i1,i3,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i1,i2,i3,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i2,i1,i3,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i3,i5,i2,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i5,i3,i2,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i1,i5,i2,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i5,i1,i2,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i1,i3,i2,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i3,i1,i2,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i2,i5,i1,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i5,i2,i1,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i3,i5,i1,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i5,i3,i1,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i3,i2,i1,i4) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i2,i3,i1,i4) = C5m(j,i1,i2,i3,i4,i5)

        C5m(j,i1,i2,i5,i4,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i5,i2,i4,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i1,i5,i4,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i5,i1,i4,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i1,i2,i4,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i2,i1,i4,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i2,i4,i5,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i4,i2,i5,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i1,i4,i5,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i4,i1,i5,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i1,i2,i5,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i2,i1,i5,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i5,i4,i2,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i4,i5,i2,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i1,i4,i2,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i4,i1,i2,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i1,i5,i2,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i5,i1,i2,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i2,i4,i1,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i4,i2,i1,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i5,i4,i1,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i4,i5,i1,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i5,i2,i1,i3) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i2,i5,i1,i3) = C5m(j,i1,i2,i3,i4,i5)

        C5m(j,i1,i5,i3,i4,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i3,i5,i4,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i1,i3,i4,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i3,i1,i4,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i1,i5,i4,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i5,i1,i4,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i5,i4,i3,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i4,i5,i3,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i1,i4,i3,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i4,i1,i3,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i1,i5,i3,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i5,i1,i3,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i3,i4,i5,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i1,i4,i3,i5,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i1,i4,i5,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i4,i1,i5,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i1,i3,i5,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i3,i1,i5,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i5,i4,i1,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i4,i5,i1,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i3,i4,i1,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i4,i3,i1,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i3,i5,i1,i2) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i5,i3,i1,i2) = C5m(j,i1,i2,i3,i4,i5)

        C5m(j,i5,i2,i3,i4,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i3,i2,i4,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i5,i3,i4,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i3,i5,i4,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i5,i2,i4,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i2,i5,i4,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i2,i4,i3,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i4,i2,i3,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i5,i4,i3,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i4,i5,i3,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i5,i2,i3,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i2,i5,i3,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i3,i4,i2,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i5,i4,i3,i2,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i5,i4,i2,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i4,i5,i2,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i5,i3,i2,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i3,i5,i2,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i2,i4,i5,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i3,i4,i2,i5,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i3,i4,i5,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i2,i4,i3,i5,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i3,i2,i5,i1) = C5m(j,i1,i2,i3,i4,i5)
        C5m(j,i4,i2,i3,i5,i1) = C5m(j,i1,i2,i3,i4,i5)
      end do
      end do
      end do
      end do
      end do
      end do

      do i1=1,3
      do i2=1,3
      do i3=1,3
      do i4=1,3
      do i5=1,3
        do k =1,3
            S6hat(k,i1,i2,i3,i4,i5) = C5m(k,i1,i2,i3,i4,i5) -
     &                C5m(0,i1,i2,i3,i4,i5)
        end do
      end do
      end do
      end do
      end do
      end do

      if (rankord.eq.4) goto 99

 99     continue

c choose reduction formulas with biggest denominators

      maxz2ff = 0d0
      maxzadj = 0d0
      do i1=1,3
      do i2=1,3
        if (abs(z2ff(i1,i2)).gt.maxz2ff) then
          i = i1
          k = i2
          maxz2ff = abs(z2ff(i1,i2))
        end if
        if (abs(zadj(i1,i2)).gt.maxzadj) then
          m = i1
          l = i2
          maxzadj = abs(zadj(i1,i2))
        end if
      end do
      end do
      l2 = mod(l,3)+1
      l3 = mod(l+1,3)+1
        
c      write(testout,*) 'iklm best',i,k,m,l,l2,l3

 01   format(' D0gy = ',G20.14,' + i* ',G20.14)
 11   format(' D1gy(',i1,') = ',G20.14,' + i* ',G20.14)
 20   format(' D2gy(0,0) = ',G20.14,' + i* ',G20.14)
 22   format(' D2gy(',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 31   format(' D3gy(0,0,',i1,') = ',G20.14,' + i* ',G20.14)
 33   format(' D3gy(',i1,',',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 40   format(' D4gy(0,0,0,0) = ',G20.14,' + i* ',G20.14)
 42   format(' D4gy(0,0,',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 44   format(' D4gy(',i1,',',i1,',',i1,',',i1,') = ',
     &    G20.14,' + i* ',G20.14)
 51   format(' D5gy(0,0,0,0,',i1,') = ',G20.14,' + i* ',G20.14)
 53   format(' D5gy(0,0,',i1,',',i1,',',i1,') = ',
     &    G20.14,' + i* ',G20.14)
 55   format(' D5gy(',i1,',',i1,',',i1,',',i1,',',i1,') = ',
     &    G20.14,' + i* ',G20.14)
 60   format(' D6gy(0,0,0,0,0,0) = ',G20.14,' + i* ',G20.14)
 62   format(' D6gy(0,0,0,0,',i1,',',i1,') = ',G20.14,' + i* ',G20.14)
 64   format(' D6gy(0,0,',i1,',',i1,',',i1,',',i1,') = ',
     &    G20.14,' + i* ',G20.14)
 66   format(' D6gy(',i1,',',i1,',',i1,',',i1,',',i1,',',i1,') = ',
     &    G20.14,' + i* ',G20.14)

c      do 101 i=2,3
c      do 101 k=2,3       
c      do 101 l=1,3
c      do 101 m=1,3       
      
c      write(testout,*) 'z2ff',z2ff
c      write(testout,*) 'zadj',zadj
c      write(testout,*) 'ik',i,k,z2ff(i,k)
c      write(testout,*) 'lm',m,l,zadj(m,l)

c      write(testout,*) 'ordgy4 = ',ordgy4
      do 100 r=0,rankord/2
c        write(testout,*) 'ordgy4 = ',ordgy4,r,m
        
      D2(0,0) = ( zadj(m,1)*S2hat(1,l) 
     &          + zadj(m,2)*S2hat(2,l)
     &          + zadj(m,3)*S2hat(3,l)
     &          -zadjf(m)*D1(l) 
     &      - detz * D2(m,l) )/(2d0*zadj(m,l))
c      write(testout,20) 
c     &    zadj(m,1)*S2hat(1,l), 
c     &          + zadj(m,2)*S2hat(2,l),
c     &          + zadj(m,3)*S2hat(3,l),
c     &          -zadjf(m)*D1(l) ,
c     &      - detz * D2(m,l) , (2d0*zadj(m,l))
c      write(testout,20) D2(0,0)

      D0 = ( zadj(i,k) * ( 2d0* D2(0,0) - C0m(0) )
     &      + z2f(i,k,1) * S1hat(1)
     &      + z2f(i,k,2) * S1hat(2)
     &      + z2f(i,k,3) * S1hat(3)
     &      - zadjf(k) * D1(i) )
     &     /z2ff(i,k)
c      write(testout,1) D0,
c     &   zadj(i,k) * ( 2d0* D2(0,0) - C0m(0) ),
c     &    z2f(i,k,1) * S1hat(1),z2f(i,k,2) * S1hat(2), 
c     &    z2f(i,k,3) * S1hat(3),- zadjf(k) * D1(i),z2ff(i,k)
        
      if (rankord-2*r.le.0) goto 100

      do i2=1,3
        R3(m,l,i2) = .5d0*( zadj(m,1)*S3hat(1,l,i2)
     &                  + zadj(m,2)*S3hat(2,l,i2)
     &                  + zadj(m,3)*S3hat(3,l,i2)
     &                  -zadjf(m)*D2(l,i2) - detz*D3(m,l,i2) )
      end do

c>      D3(0,0,l) = ( zadj(m,1)*S3hat(1,l,l) 
c>     &          + zadj(m,2)*S3hat(2,l,l)
c>     &          + zadj(m,3)*S3hat(3,l,l)
c>     &          -zadjf(m)*D2(l,l) 
c>     &      - detz * D3(m,l,l) )/(4d0*zadj(m,l))
c>
c>      D3(0,0,l2) = ( zadj(m,1)*S3hat(1,l,l2) 
c>     &          + zadj(m,2)*S3hat(2,l,l2)
c>     &          + zadj(m,3)*S3hat(3,l,l2)
c>     &          - 2d0*zadj(m,l2)*D3(0,0,l)
c>     &          -zadjf(m)*D2(l,l2) 
c>     &      - detz * D3(m,l,l2) )/(2d0*zadj(m,l))
c>
c>      D3(0,0,l3) = ( zadj(m,1)*S3hat(1,l,l3) 
c>     &          + zadj(m,2)*S3hat(2,l,l3)
c>     &          + zadj(m,3)*S3hat(3,l,l3)
c>     &          - 2d0*zadj(m,l3)*D3(0,0,l)
c>     &          -zadjf(m)*D2(l,l3) 
c>     &      - detz * D3(m,l,l3) )/(2d0*zadj(m,l))

      D3(0,0,l) =  R3(m,l,l)/(2d0*zadj(m,l))
      D3(0,0,l2) =  (R3(m,l,l2)-zadj(m,l2)*D3(0,0,l))/(zadj(m,l))
      D3(0,0,l3) =  (R3(m,l,l3)-zadj(m,l3)*D3(0,0,l))/(zadj(m,l))


c      do i1=1,3
c        write(testout,31) i1,D3(0,0,i1)
c      end do

      do i1=1,3
        do j =1,3
          S2mod(j,i1) = S2hat(j,i1)
        end do
        S2mod(i1,i1) = S2mod(i1,i1) - 2d0*D2(0,0)
      end do

      do i1=1,3
        D1(i1) = ( zadj(i,k) * ( 4d0* D3(0,0,i1) - C1m(0,i1) )
     &      + z2f(i,k,1) * S2mod(1,i1)
     &      + z2f(i,k,2) * S2mod(2,i1)
     &      + z2f(i,k,3) * S2mod(3,i1)
     &      - zadjf(k) * D2(i,i1) )
     &      /z2ff(i,k)
c        write(testout,*) i1,
c     &       zadj(i,k) * ( 4d0* D3(0,0,i1) - C1m(0,i1) ),
c     &      + z2f(i,k,1) , S2mod(1,i1),
c     &      + z2f(i,k,2) , S2mod(2,i1),
c     &      + z2f(i,k,3) , S2mod(3,i1),
c     &      - zadjf(k) * D2(i,i1)
c        write(testout,11) i1,D1(i1)
      end do

      if (rankord-2*r.le.1) goto 100

      do i2=1,3
      do i3=1,3
        R4(m,l,i2,i3) = .5d0*( zadj(m,1)*S4hat(1,l,i2,i3)
     &                  + zadj(m,2)*S4hat(2,l,i2,i3)
     &                  + zadj(m,3)*S4hat(3,l,i2,i3)
     &                  -zadjf(m)*D3(l,i2,i3) - detz*D4(m,l,i2,i3) )
      end do
      end do

      D4(0,0,l,l) = R4(m,l,l,l)/(3d0*zadj(m,l))
      D4(0,0,l,l2) = (R4(m,l,l,l2)-D4(0,0,l,l)*zadj(m,l2))
     &               /(2d0*zadj(m,l))
      D4(0,0,l2,l2) = (R4(m,l,l2,l2)-2d0*D4(0,0,l,l2)*zadj(m,l2))
     &               /(zadj(m,l))
      D4(0,0,l,l3) = (R4(m,l,l,l3)-D4(0,0,l,l)*zadj(m,l3))
     &               /(2d0*zadj(m,l))
      D4(0,0,l3,l3) = (R4(m,l,l3,l3)-2d0*D4(0,0,l,l3)*zadj(m,l3))
     &               /(zadj(m,l))
      D4(0,0,l2,l3) = (R4(m,l,l2,l3)-D4(0,0,l,l2)*zadj(m,l3)
     &                              -D4(0,0,l,l3)*zadj(m,l2))
     &               /(zadj(m,l))

      D4(0,0,l2,l) = D4(0,0,l,l2)
      D4(0,0,l3,l) = D4(0,0,l,l3)
      D4(0,0,l3,l2)= D4(0,0,l2,l3)

      do i1=1,3
      do i2=1,3
c        write(testout,42) i1,i2,D4(0,0,i1,i2)
      end do
      end do

      do i1=1,3
      do i2=1,3
        do j =1,3
          S3mod(j,i1,i2) = S3hat(j,i1,i2)
        end do
        S3mod(i1,i1,i2) = S3mod(i1,i1,i2) - 2d0*D3(0,0,i2)
        S3mod(i2,i1,i2) = S3mod(i2,i1,i2) - 2d0*D3(0,0,i1)
      end do
      end do

      do i1=1,3
      do i2=1,3
        D2(i1,i2) = ( zadj(i,k) * ( 6d0* D4(0,0,i1,i2) - C2m(0,i1,i2) )
     &      + z2f(i,k,1) * S3mod(1,i1,i2)
     &      + z2f(i,k,2) * S3mod(2,i1,i2)
     &      + z2f(i,k,3) * S3mod(3,i1,i2)
     &      - zadjf(k) * D3(i,i1,i2) )
     &      /z2ff(i,k)
c        write(testout,22) i1,i2,D2(i1,i2)
      end do
      end do


      if (rankord-2*r.le.2) goto 100

      do i2=1,3
      do i3=1,3
      do i4=1,3
        R5(m,l,i2,i3,i4) = .5d0*( zadj(m,1)*S5hat(1,l,i2,i3,i4)
     &                  + zadj(m,2)*S5hat(2,l,i2,i3,i4)
     &                  + zadj(m,3)*S5hat(3,l,i2,i3,i4)
     &                  -zadjf(m)*D4(l,i2,i3,i4) 
c     &                  - detz*D5(m,l,i2,i3,i4) 
     &                   )
      end do
      end do
      end do
c        write(*,*) 'TEST',
c     &   R5(2,2,2,2,2),zadj(2,1)*S5hat(1,2,2,2,2),
c     &                  + zadj(2,2)*S5hat(2,2,2,2,2),
c     &                  + zadj(2,3)*S5hat(3,2,2,2,2),
c     &                  -zadjf(2)*D4(2,2,2,2), - detz*D5(2,2,2,2,2) 

      D5(0,0,l,l,l) = R5(m,l,l,l,l)/(4d0*zadj(m,l))
      D5(0,0,l,l,l2) = (R5(m,l,l,l,l2)-D5(0,0,l,l,l)*zadj(m,l2))
     &               /(3d0*zadj(m,l))
      D5(0,0,l,l2,l2) = (R5(m,l,l,l2,l2)-2d0*D5(0,0,l,l,l2)*zadj(m,l2))
     &               /(2d0*zadj(m,l))
      D5(0,0,l2,l2,l2) = (R5(m,l,l2,l2,l2)
     &                   -3d0*D5(0,0,l,l2,l2)*zadj(m,l2))
     &               /(zadj(m,l))
      D5(0,0,l,l,l3) = (R5(m,l,l,l,l3)-D5(0,0,l,l,l)*zadj(m,l3))
     &               /(3d0*zadj(m,l))
      D5(0,0,l,l3,l3) = (R5(m,l,l,l3,l3)-2d0*D5(0,0,l,l,l3)*zadj(m,l3))
     &               /(2d0*zadj(m,l))
      D5(0,0,l3,l3,l3) = (R5(m,l,l3,l3,l3)
     &                   -3d0*D5(0,0,l,l3,l3)*zadj(m,l3))
     &               /(zadj(m,l))
      D5(0,0,l,l2,l3) = (R5(m,l,l,l2,l3)-D5(0,0,l,l,l2)*zadj(m,l3)
     &                              -D5(0,0,l,l,l3)*zadj(m,l2))
     &               /(2d0*zadj(m,l))
      D5(0,0,l2,l2,l3) = (R5(m,l,l2,l2,l3)-D5(0,0,l,l2,l2)*zadj(m,l3)
     &                              -2d0*D5(0,0,l,l2,l3)*zadj(m,l2))
     &               /(zadj(m,l))
      D5(0,0,l2,l3,l3) = (R5(m,l,l2,l3,l3)-D5(0,0,l,l3,l3)*zadj(m,l2)
     &                              -2d0*D5(0,0,l,l2,l3)*zadj(m,l3))
     &               /(zadj(m,l))

      D5(0,0,l,l2,l) = D5(0,0,l,l,l2)
      D5(0,0,l2,l,l) = D5(0,0,l,l,l2)
      D5(0,0,l,l3,l) = D5(0,0,l,l,l3)
      D5(0,0,l3,l,l) = D5(0,0,l,l,l3)
      D5(0,0,l3,l3,l)= D5(0,0,l,l3,l3)
      D5(0,0,l3,l,l3)= D5(0,0,l,l3,l3)
      D5(0,0,l2,l2,l)= D5(0,0,l,l2,l2)
      D5(0,0,l2,l,l2)= D5(0,0,l,l2,l2)
      D5(0,0,l,l3,l2)= D5(0,0,l,l2,l3)
      D5(0,0,l2,l,l3)= D5(0,0,l,l2,l3)
      D5(0,0,l2,l3,l)= D5(0,0,l,l2,l3)
      D5(0,0,l3,l,l2)= D5(0,0,l,l2,l3)
      D5(0,0,l3,l2,l)= D5(0,0,l,l2,l3)
      D5(0,0,l3,l2,l2)= D5(0,0,l2,l2,l3)
      D5(0,0,l2,l3,l2)= D5(0,0,l2,l2,l3)
      D5(0,0,l3,l2,l3)= D5(0,0,l2,l3,l3)
      D5(0,0,l3,l3,l2)= D5(0,0,l2,l3,l3)

c      do i1=1,3
c      do i2=1,3
c      do i3=1,3
c        write(testout,53) i1,i2,i3,D5(0,0,i1,i2,i3)
c      end do
c      end do
c      end do

      do i1=1,3
      do i2=1,3
      do i3=1,3
        do j =1,3
          S4mod(j,i1,i2,i3) = S4hat(j,i1,i2,i3)
        end do
        S4mod(i1,i1,i2,i3) = S4mod(i1,i1,i2,i3) -2d0*D4(0,0,i2,i3)
        S4mod(i2,i1,i2,i3) = S4mod(i2,i1,i2,i3) -2d0*D4(0,0,i1,i3)
        S4mod(i3,i1,i2,i3) = S4mod(i3,i1,i2,i3) -2d0*D4(0,0,i1,i2)
      end do
      end do
      end do

      do i1=1,3
      do i2=1,3
      do i3=1,3
        D3(i1,i2,i3) = ( zadj(i,k) * ( 8d0* D5(0,0,i1,i2,i3) 
     &                             - C3m(0,i1,i2,i3) )
     &      + z2f(i,k,1) * S4mod(1,i1,i2,i3)
     &      + z2f(i,k,2) * S4mod(2,i1,i2,i3)
     &      + z2f(i,k,3) * S4mod(3,i1,i2,i3)
     &      - zadjf(k) * D4(i,i1,i2,i3) )
     &      /z2ff(i,k)
c        write(testout,33) i1,i2,i3,D3(i1,i2,i3)
      end do
      end do
      end do

      if (rankord-2*r.le.3) goto 100

        S002(0,0) = 2d0*mm02*D2(0,0) + 2d0*C2m(0,0,0)

        do i1=1,3
          do i2=1,3
            S002(i1,i2) = 2d0*mm02*D2(i1,i2) + 2d0*C2m(0,i1,i2)
          end do
        end do
        
        D4(0,0,0,0) = zadj(m,l)*(S002(0,0) + 1d0/6d0) 
     &      - detz*D4(0,0,m,l)
        do a=1,3
          D4(0,0,0,0) = D4(0,0,0,0) - zadj(m,l)*S4hat(a,0,0,a) 
     &        + zadj(a,l)*S4hat(a,0,0,m)
          do b=1,3
            if(b.ne.k.and.a.ne.l) then
              D4(0,0,0,0) = D4(0,0,0,0) - zadj2(m,b,l,a)*( 
     &            + f(b)*S3hat(a,0,0)
     &            - f(a)*f(b)*D2(0,0)   )
            end if
          end do
        end do
        D4(0,0,0,0) = D4(0,0,0,0) /(8d0*zadj(m,l))
c        write(testout,40) D4(0,0,0,0)

      do i2=1,3
      do i3=1,3
      do i4=1,3
      do i5=1,3
        R6(m,l,i2,i3,i4,i5) = .5d0*( zadj(m,1)*S6hat(1,l,i2,i3,i4,i5)
     &                  + zadj(m,2)*S6hat(2,l,i2,i3,i4,i5)
     &                  + zadj(m,3)*S6hat(3,l,i2,i3,i4,i5)
c     &                  - zadjf(m)*D5(l,i2,i3,i4,i5) 
c     &                  - detz*D6(m,l,i2,i3,i4,i5) 
     &                   )
      end do
      end do
      end do
      end do

      D6(0,0,l,l,l,l) = R6(m,l,l,l,l,l)/(5d0*zadj(m,l))
      D6(0,0,l,l,l,l2) = (R6(m,l,l,l,l,l2)-D6(0,0,l,l,l,l)*zadj(m,l2))
     &               /(4d0*zadj(m,l))
      D6(0,0,l,l,l2,l2) = (R6(m,l,l,l,l2,l2)
     &                    - 2d0*D6(0,0,l,l,l,l2)*zadj(m,l2))
     &               /(3d0*zadj(m,l))
      D6(0,0,l,l2,l2,l2) = (R6(m,l,l,l2,l2,l2)
     &                    - 3d0*D6(0,0,l,l,l2,l2)*zadj(m,l2))
     &               /(2d0*zadj(m,l))
      D6(0,0,l2,l2,l2,l2) = (R6(m,l,l2,l2,l2,l2)
     &                    - 4d0*D6(0,0,l,l2,l2,l2)*zadj(m,l2))
     &               /(zadj(m,l))
      D6(0,0,l,l,l,l3) = (R6(m,l,l,l,l,l3)-D6(0,0,l,l,l,l)*zadj(m,l3))
     &               /(4d0*zadj(m,l))
      D6(0,0,l,l,l3,l3) = (R6(m,l,l,l,l3,l3)
     &                   - 2d0*D6(0,0,l,l,l,l3)*zadj(m,l3))
     &               /(3d0*zadj(m,l))
      D6(0,0,l,l3,l3,l3) = (R6(m,l,l,l3,l3,l3)
     &                   - 3d0*D6(0,0,l,l,l3,l3)*zadj(m,l3))
     &               /(2d0*zadj(m,l))
      D6(0,0,l3,l3,l3,l3) = (R6(m,l,l3,l3,l3,l3)
     &                   - 4d0*D6(0,0,l,l3,l3,l3)*zadj(m,l3))
     &               /(zadj(m,l))
      D6(0,0,l,l,l2,l3) = (R6(m,l,l,l,l2,l3)
     &                    - D6(0,0,l,l,l,l2)*zadj(m,l3)
     &                    - D6(0,0,l,l,l,l3)*zadj(m,l2))
     &               /(3d0*zadj(m,l))
      D6(0,0,l,l2,l2,l3) = (R6(m,l,l,l2,l2,l3)
     &                     - D6(0,0,l,l,l2,l2)*zadj(m,l3)
     &                     - 2d0*D6(0,0,l,l,l2,l3)*zadj(m,l2))
     &               /(2d0*zadj(m,l))
      D6(0,0,l,l2,l3,l3) = (R6(m,l,l,l2,l3,l3)
     &                     - D6(0,0,l,l,l3,l3)*zadj(m,l2)
     &                     - 2d0*D6(0,0,l,l,l2,l3)*zadj(m,l3))
     &               /(2d0*zadj(m,l))
      D6(0,0,l2,l2,l2,l3) = (R6(m,l,l2,l2,l2,l3)
     &                     - D6(0,0,l,l2,l2,l2)*zadj(m,l3)
     &                     - 3d0*D6(0,0,l,l2,l2,l3)*zadj(m,l2))
     &               /(zadj(m,l))
      D6(0,0,l2,l3,l3,l3) = (R6(m,l,l2,l3,l3,l3)
     &                     - D6(0,0,l,l3,l3,l3)*zadj(m,l2)
     &                     - 3d0*D6(0,0,l,l2,l3,l3)*zadj(m,l3))
     &               /(zadj(m,l))
      D6(0,0,l2,l2,l3,l3) = (R6(m,l,l2,l2,l3,l3)
     &                     - 2d0*D6(0,0,l,l2,l3,l3)*zadj(m,l2)
     &                     - 2d0*D6(0,0,l,l2,l2,l3)*zadj(m,l3))
     &               /(zadj(m,l))

      D6(0,0,l,l,l2,l)    = D6(0,0,l,l,l,l2)
      D6(0,0,l,l2,l,l)    = D6(0,0,l,l,l,l2)
      D6(0,0,l2,l,l,l)    = D6(0,0,l,l,l,l2)
      D6(0,0,l,l,l3,l)    = D6(0,0,l,l,l,l3)
      D6(0,0,l,l3,l,l)    = D6(0,0,l,l,l,l3)
      D6(0,0,l3,l,l,l)    = D6(0,0,l,l,l,l3)
      D6(0,0,l,l3,l3,l)   = D6(0,0,l,l,l3,l3)
      D6(0,0,l,l3,l,l3)   = D6(0,0,l,l,l3,l3)
      D6(0,0,l3,l,l,l3)   = D6(0,0,l,l,l3,l3)
      D6(0,0,l3,l,l3,l)   = D6(0,0,l,l,l3,l3)
      D6(0,0,l3,l3,l,l)   = D6(0,0,l,l,l3,l3)
      D6(0,0,l,l2,l2,l)   = D6(0,0,l,l,l2,l2)
      D6(0,0,l,l2,l,l2)   = D6(0,0,l,l,l2,l2)
      D6(0,0,l2,l,l,l2)   = D6(0,0,l,l,l2,l2)
      D6(0,0,l2,l,l2,l)   = D6(0,0,l,l,l2,l2)
      D6(0,0,l2,l2,l,l)   = D6(0,0,l,l,l2,l2)
      D6(0,0,l,l,l3,l2)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l,l2,l,l3)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l,l2,l3,l)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l,l3,l,l2)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l,l3,l2,l)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l2,l,l,l3)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l2,l,l3,l)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l2,l3,l,l)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l3,l,l,l2)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l3,l,l2,l)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l3,l2,l,l)   = D6(0,0,l,l,l2,l3)
      D6(0,0,l,l3,l2,l2)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l,l2,l3,l2)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l3,l,l2,l2)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l2,l,l3,l2)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l2,l,l2,l3)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l3,l2,l,l2)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l2,l3,l,l2)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l2,l2,l,l3)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l3,l2,l2,l)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l2,l3,l2,l)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l2,l2,l3,l)  = D6(0,0,l,l2,l2,l3)
      D6(0,0,l,l3,l2,l3)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l,l3,l3,l2)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l2,l,l3,l3)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l3,l,l2,l3)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l3,l,l3,l2)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l2,l3,l,l3)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l3,l2,l,l3)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l3,l3,l,l2)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l2,l3,l3,l)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l3,l2,l3,l)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l3,l3,l2,l)  = D6(0,0,l,l2,l3,l3)
      D6(0,0,l2,l2,l3,l2) = D6(0,0,l2,l2,l2,l3)
      D6(0,0,l2,l3,l2,l2) = D6(0,0,l2,l2,l2,l3)
      D6(0,0,l3,l2,l2,l2) = D6(0,0,l2,l2,l2,l3)
      D6(0,0,l3,l2,l3,l3) = D6(0,0,l2,l3,l3,l3)
      D6(0,0,l3,l3,l2,l3) = D6(0,0,l2,l3,l3,l3)
      D6(0,0,l3,l3,l3,l2) = D6(0,0,l2,l3,l3,l3)
      D6(0,0,l2,l3,l2,l3) = D6(0,0,l2,l2,l3,l3)
      D6(0,0,l2,l3,l3,l2) = D6(0,0,l2,l2,l3,l3)
      D6(0,0,l3,l2,l2,l3) = D6(0,0,l2,l2,l3,l3)
      D6(0,0,l3,l2,l3,l2) = D6(0,0,l2,l2,l3,l3)
      D6(0,0,l3,l3,l2,l2) = D6(0,0,l2,l2,l3,l3)

c      do i1=1,3
c      do i2=1,3
c      do i3=1,3
c      do i4=1,3
c        write(testout,64) i1,i2,i3,i4,D6(0,0,i1,i2,i3,i4)
c      end do
c      end do
c      end do
c      end do

      do i1=1,3
      do i2=1,3
      do i3=1,3
      do i4=1,3
        do j =1,3
          S5mod(j,i1,i2,i3,i4) = S5hat(j,i1,i2,i3,i4)
        end do
        S5mod(i1,i1,i2,i3,i4) = S5mod(i1,i1,i2,i3,i4) 
     &      - 2d0*D5(0,0,i2,i3,i4)
        S5mod(i2,i1,i2,i3,i4) = S5mod(i2,i1,i2,i3,i4) 
     &      - 2d0*D5(0,0,i1,i3,i4)
        S5mod(i3,i1,i2,i3,i4) = S5mod(i3,i1,i2,i3,i4) 
     &      - 2d0*D5(0,0,i1,i2,i4)
        S5mod(i4,i1,i2,i3,i4) = S5mod(i4,i1,i2,i3,i4) 
     &      - 2d0*D5(0,0,i1,i2,i3)
      end do
      end do
      end do
      end do

      do i1=1,3
      do i2=1,3
      do i3=1,3
      do i4=1,3
        D4(i1,i2,i3,i4) = ( zadj(i,k) * ( 8d0* D6(0,0,i1,i2,i3,i4) 
     &                             - C4m(0,i1,i2,i3,i4) )
     &      + z2f(i,k,1) * S5mod(1,i1,i2,i3,i4)
     &      + z2f(i,k,2) * S5mod(2,i1,i2,i3,i4)
     &      + z2f(i,k,3) * S5mod(3,i1,i2,i3,i4)
c     &      - zadjf(k) * D5(i,i1,i2,i3,i4) 
     &      )/z2ff(i,k)
c        write(testout,44) i1,i2,i3,i4,D4(i1,i2,i3,i4)
      end do
      end do
      end do
      end do

 100  continue

 101  continue

c      write(*,*) 'D0gy = ',D0

      end

***********************************************************************
      subroutine cDgp12345(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &    D0,D1,D2,D3,D4,D5,rank)
***********************************************************************
*     4-point tensor coefficient functions D0,D1,D2,D3,D4,D5          *
*     up to 5 integration momenta in numerator                        *
*     coefficients calculated up to *rank* momenta in numerator       *
*     expansion about (pi pj) = 0 to order ordgp4                     *
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
*     i-th denominator cancelled: Cxwi                                *
***********************************************************************
*     04.04.05 Ansgar Denner    last changed  05.04.05                *
***********************************************************************
      implicit   none
      complex*16 p10,p21,p32,p30,p20,p31
      complex*16 q10,q21,q32,q30,q20,q31
      complex*16 m02,m12,m22,m32
      complex*16 mm02,mm12,mm22,mm32,f(3),zadjf(3)
      complex*16 k1k2,k1k3,k2k3,z(3,3),zadj(3,3),detz
      complex*16 zadj2(3,3,3,3)
      complex*16 C0w3,C1w3(2),C2w3(0:2,0:2),C3w3(0:2,0:2,0:2)
     &    ,C4w3(0:2,0:2,0:2,0:2),C5w3(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w3(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w2,C1w2(2),C2w2(0:2,0:2),C3w2(0:2,0:2,0:2)
     &    ,C4w2(0:2,0:2,0:2,0:2),C5w2(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w2(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w1,C1w1(2),C2w1(0:2,0:2),C3w1(0:2,0:2,0:2)
     &    ,C4w1(0:2,0:2,0:2,0:2),C5w1(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w1(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0w0,C1w0(2),C2w0(0:2,0:2),C3w0(0:2,0:2,0:2)
     &    ,C4w0(0:2,0:2,0:2,0:2),C5w0(0:2,0:2,0:2,0:2,0:2)
     &    ,C6w0(0:2,0:2,0:2,0:2,0:2,0:2)
      complex*16 C0m(0:3),C1m(0:3,3),C2m(0:3,0:3,0:3)
     &    ,C3m(0:3,0:3,0:3,0:3),C4m(0:3,0:3,0:3,0:3,0:3)
c     &    ,C5m(0:3,0:3,0:3,0:3,0:3,0:3),C6m(0:3,0:3,0:3,0:3,0:3,0:3,0:3)
      complex*16 S1hat(3),S2hat(3,3),S3hat(3,0:3,0:3),
     &    S4hat(3,0:3,0:3,3),S5hat(0:3,0:3,0:3,0:3,0:3)
      complex*16 S2mod(3,3),S3mod(3,0:3,0:3),
     &    S4mod(3,0:3,0:3,3),S5mod(0:3,0:3,0:3,0:3,0:3)
      complex*16 S00,S001(3),S002(0:3,0:3),S003(0:3,0:3,0:3)
      complex*16 D0,D1(3),D2(0:3,0:3),D3(0:3,0:3,0:3)
      complex*16 D4(0:3,0:3,0:3,0:3),D5(0:3,0:3,0:3,0:3,0:3)
c      complex*16 cD0f
      complex*16 elimcminf2
      real*8     maxf
      integer    rank
      integer    i1,i2,i3,i4,j,k,m,r
      integer    ltest,sym
c      integer    testout

      integer    ordg3,ordgy3,ordgp3
      integer    ordg4,ordgy4,ordgp4
      common /ctord3/ ordg3,ordgy3,ordgp3
      common /ctord4/ ordg4,ordgy4,ordgp4
      common /ltest/  ltest
      common /sym/ sym

      data       C1m /12*0d0/,C2m /64*0d0/,C3m/256*0d0/,C4m/1024*0d0/
      data       zadj2/81*0d0/
c      data       testout /43/

c      write(testout,*)'Dg12345 in',p10,p21,p32,p30,p20,p31,
c     &    m02,m12,m22,m32
c     &    ,rank,ordgp4

      if (rank.lt.0) then
        write(*,*) 'rank < 0 not implemented in cDgp12345'
        write(*,*) 'rank = ',rank
        stop
      else if (rank+ordgp4.gt.4) then
        write(*,*) 'rank+ordgp4 > 4 not implemented in cDgp12345'
        write(*,*) 'rank  = ',rank
        write(*,*) 'ordgp4 = ',ordgp4
        if (rank.gt.4) then
          write(*,*) 'rank > 4 not implemented in cDgp12345'
          write(*,*) 'rank = ',rank
          stop
        end if
      end if

c     calculate coefficients
      

c     write(*,*) 'cDp1',p21,p32,p31,m12,m22,m32,rank
      call cCp12345(p21,p32,p31,m12,m22,m32,C0w0,C1w0,C2w0,C3w0,C4w0,
     &    C5w0,C6w0,rank+ordgp4)
c     write(*,*) 'cDp2',p21,p32,p31,m12,m22,m32,rank
      call cCp12345(p20,p32,p30,m02,m22,m32,C0w1,C1w1,C2w1,C3w1,C4w1,
     &    C5w1,C6w1,rank+ordgp4)
c     write(*,*) 'cDp3',p21,p32,p31,m12,m22,m32,rank
      call cCp12345(p10,p31,p30,m02,m12,m32,C0w2,C1w2,C2w2,C3w2,C4w2,
     &    C5w2,C6w2,rank+ordgp4)
c     write(*,*) 'cDp4',p21,p32,p31,m12,m22,m32,rank
      call cCp12345(p10,p21,p20,m02,m12,m22,C0w3,C1w3,C2w3,C3w3,C4w3,
     &    C5w3,C6w3,rank+ordgp4)
      
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

      z(1,1) = 2d0*q10
      z(2,2) = 2d0*q20
      z(3,3) = 2d0*q30
      z(1,2) = q10+q20-q21
      z(1,3) = q10+q30-q31
      z(2,3) = q20+q30-q32
      z(2,1) = z(1,2)
      z(3,1) = z(1,3)
      z(3,2) = z(2,3)
      
      k1k2 = (q10+q20-q21)
      k1k3 = (q10+q30-q31)
      k2k3 = (q20+q30-q32)
      detz = 8d0*q10*q30*q20+2D0*k1k2*k1k3*k2k3
     &    -2d0*(q10*k2k3*k2k3+q20*k1k3*k1k3+q30*k1k2*k1k2)

c      write(testout,*) 'detz = ',detz,detz/
c     &    dmax1(abs(k1k2),abs(k1k3),abs(k2k3),
c     &    2*abs(q10),2*abs(q20),2*abs(q30))**3

      zadj(1,1) = (4d0*q30*q20-k2k3*k2k3)
      zadj(1,2) = (k1k3*k2k3-2d0*q30*k1k2)
      zadj(1,3) = (k1k2*k2k3-2d0*q20*k1k3)
      zadj(2,1) = zadj(1,2)
      zadj(2,2) = (4d0*q10*q30-k1k3*k1k3)
      zadj(2,3) = (k1k2*k1k3-2d0*q10*k2k3)
      zadj(3,1) = zadj(1,3)
      zadj(3,2) = zadj(2,3)
      zadj(3,3) = (4d0*q10*q20-k1k2*k1k2)

      zadj2(1,2,1,2) = -z(3,3)
      zadj2(1,2,2,1) = z(3,3)
      zadj2(2,1,1,2) = z(3,3)
      zadj2(2,1,2,1) = -z(3,3)
      zadj2(1,3,1,3) = -z(2,2)
      zadj2(1,3,3,1) = z(2,2)
      zadj2(3,1,1,3) = z(2,2)
      zadj2(3,1,3,1) = -z(2,2)
      zadj2(3,2,3,2) = -z(1,1)
      zadj2(3,2,2,3) = z(1,1)
      zadj2(2,3,3,2) = z(1,1)
      zadj2(2,3,2,3) = -z(1,1)
      zadj2(1,2,1,3) = z(3,2)
      zadj2(1,2,3,1) = -z(3,2)
      zadj2(2,1,1,3) = -z(3,2)
      zadj2(2,1,3,1) = z(3,2)
      zadj2(1,3,1,2) = z(2,3)
      zadj2(1,3,2,1) = -z(2,3)
      zadj2(3,1,1,2) = -z(2,3)
      zadj2(3,1,2,1) = z(2,3)
      zadj2(3,2,1,2) = z(1,3)
      zadj2(3,2,2,1) = -z(1,3)
      zadj2(2,3,1,2) = -z(1,3)
      zadj2(2,3,2,1) = z(1,3)
      zadj2(1,2,3,2) = z(3,1)
      zadj2(1,2,2,3) = -z(3,1)
      zadj2(2,1,3,2) = -z(3,1)
      zadj2(2,1,2,3) = z(3,1)
      zadj2(1,3,2,3) = z(2,1)
      zadj2(1,3,3,2) = -z(2,1)
      zadj2(3,1,2,3) = -z(2,1)
      zadj2(3,1,3,2) = z(2,1)
      zadj2(3,2,3,1) = z(1,2)
      zadj2(3,2,1,3) = -z(1,2)
      zadj2(2,3,3,1) = -z(1,2)
      zadj2(2,3,1,3) = z(1,2)

      f(1) = q10+mm02-mm12
      f(2) = q20+mm02-mm22
      f(3) = q30+mm02-mm32
      
      zadjf(1) = zadj(1,1)*f(1)+zadj(1,2)*f(2)+zadj(1,3)*f(3)
      zadjf(2) = zadj(2,1)*f(1)+zadj(2,2)*f(2)+zadj(2,3)*f(3)
      zadjf(3) = zadj(3,1)*f(1)+zadj(3,2)*f(2)+zadj(3,3)*f(3)

c      write(testout,*) 'zadjf ',zadjf
c      write(testout,*) 'zadjf ',detz/zadjf(1),detz/zadjf(2),detz/zadjf(3)

      do i1=1,3
        D1(i1) = 0d0
c        D3(0,0,i1) = 0d0
c        D5(0,0,0,0,i1) = 0d0
        do i2=1,3
          D2(i1,i2) = 0d0
c          D4(0,0,i1,i2) = 0d0
          do i3=1,3
            D3(i1,i2,i3) = 0d0
c            D5(0,0,i1,i2,i3) = 0d0
            do i4=1,3
              D4(i1,i2,i3,i4) = 0d0
c              do i5=1,3
c                D5(i1,i2,i3,i4,i5) = 0d0
c              end do
            end do
          end do
        end do
      end do

      C0m(0) = C0w0
      C0m(1) = C0w1
      C0m(2) = C0w2
      C0m(3) = C0w3
      
      do j =1,3
        S1hat(j) = C0m(j) - C0m(0)
      end do

      if (rank+ordgp4.eq.0) goto 99
  
      C1m(0,2) = C1w0(1)
      C1m(0,3) = C1w0(2)
      C1m(0,1) = -C0m(0) - C1m(0,2) - C1m(0,3)
      C1m(1,2) = C1w1(1)
      C1m(1,3) = C1w1(2)
      C1m(2,1) = C1w2(1)
      C1m(2,3) = C1w2(2)
      C1m(3,1) = C1w3(1)
      C1m(3,2) = C1w3(2)
      
      do i1=1,3
        do j =1,3
          S2hat(j,i1) = C1m(j,i1) - C1m(0,i1)
        end do
      end do

c     >      do i1=1,3
c     >      do i2=1,3
c     >      do k=1,3
c     >        do l=1,3
c     >          test = zadj(k,l)*z(i1,i2)
c     >          do a=1,3
c     >            do b=1,3
c     >              if(b.ne.k.and.a.ne.l) then
c     >                test  = test + zadj2(k,b,l,a)*z(b,i1)*z(a,i2)
c     >              end if
c     >            end do
c     >          end do
c     >          write(testout,*) 'detz = ',i1,i2,k,l,test
c     >        end do
c     >      end do
c     >      end do
c     >      end do

        if (rank+ordgp4.eq.1) goto 99

        C2m(0,0,0) = C2w0(0,0)
        C2m(1,0,0) = C2w1(0,0)
        C2m(2,0,0) = C2w2(0,0)
        C2m(3,0,0) = C2w3(0,0)
        
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
        
        C2m(0,3,2) = C2m(0,2,3) 
        C2m(0,3,1) = C2m(0,1,3) 
        C2m(0,2,1) = C2m(0,1,2) 
        C2m(1,3,2) = C2m(1,2,3) 
        C2m(2,3,1) = C2m(2,1,3) 
        C2m(3,2,1) = C2m(3,1,2) 

        do i1=1,3
          S3hat(i1,0,0) = C2m(i1,0,0) - C2m(0,0,0)
          do i2=1,3
            do j =1,3
              S3hat(j,i1,i2) = C2m(j,i1,i2) - C2m(0,i1,i2)
            end do
          end do
        end do


        if (rank+ordgp4.eq.2) goto 99

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
        
        C3m(0,1,2,1) = C3m(0,1,1,2)
        C3m(0,2,1,1) = C3m(0,1,1,2)
        C3m(0,1,3,1) = C3m(0,1,1,3)
        C3m(0,3,1,1) = C3m(0,1,1,3)
        C3m(0,2,2,1) = C3m(0,1,2,2)
        C3m(0,2,1,2) = C3m(0,1,2,2)
        C3m(0,3,3,1) = C3m(0,1,3,3)
        C3m(0,3,1,3) = C3m(0,1,3,3)
        C3m(0,2,3,2) = C3m(0,2,2,3)
        C3m(0,3,2,2) = C3m(0,2,2,3)
        C3m(0,3,3,2) = C3m(0,2,3,3)
        C3m(0,3,2,3) = C3m(0,2,3,3)
        
        C3m(0,1,3,2) = C3m(0,1,2,3)
        C3m(0,2,1,3) = C3m(0,1,2,3)
        C3m(0,2,3,1) = C3m(0,1,2,3)
        C3m(0,3,1,2) = C3m(0,1,2,3)
        C3m(0,3,2,1) = C3m(0,1,2,3)
        
        C3m(3,1,2,1) = C3m(3,1,1,2)
        C3m(3,2,1,1) = C3m(3,1,1,2)
        C3m(2,1,3,1) = C3m(2,1,1,3)
        C3m(2,3,1,1) = C3m(2,1,1,3)
        C3m(3,2,2,1) = C3m(3,1,2,2)
        C3m(3,2,1,2) = C3m(3,1,2,2)
        C3m(2,3,3,1) = C3m(2,1,3,3)
        C3m(2,3,1,3) = C3m(2,1,3,3)
        C3m(1,2,3,2) = C3m(1,2,2,3)
        C3m(1,3,2,2) = C3m(1,2,2,3)
        C3m(1,3,3,2) = C3m(1,2,3,3)
        C3m(1,3,2,3) = C3m(1,2,3,3)
        
        do i1=1,3
          do i2=1,3
            S4hat(i1,0,0,i2) = C3m(i1,0,0,i2) - C3m(0,0,0,i2)
            do i3=1,3
              do j =1,3
                S4hat(j,i1,i2,i3) = C3m(j,i1,i2,i3) - C3m(0,i1,i2,i3)
              end do
            end do
          end do
        end do
        

        if (rank+ordgp4.eq.3) goto 99

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

c     >      do i1=1,2
c     >      do i2=i1,2
c     >        C4m(0,0,0,i1+1,i2+1) = C4w0(0,0,i1,i2)
c     >        C4m(1,0,0,i1+1,i2+1) = C4w0(0,0,i1,i2)
c     >        C4m(2,0,0,2*i1-1,2*i2-1) = C4w0(0,0,i1,i2)
c     >        C4m(3,0,0,i1+1,i2) = C4w0(0,0,i1,i2)
c     >      do i3=1,2
c     >      do i4=1,2
c     >        C4m(0,i1+1,i2+1,i3+1,i4+1) = C4w0(i1,i2,i3,i4)
c     >        C4m(1,i1+1,i2+1,i3+1,i4+1) = C4w0(i1,i2,i3,i4)
c     >        C4m(2,2*i1-1,2*i2-1,2*i3-1,2*i4-1) = C4w0(i1,i2,i3,i4)
c     >        C4m(3,i1+1,i2,i3,i4) = C4w0(i1,i2,i3,i4)
c     >      end do
c     >        C4m(0,1,i1+1,i2+1,i3+1) = -C3m(0,i1+1,i2+1,i3+1) 
c     >     &      - C4m(0,2,i1+1,i2+1,i3+1) 
c     >     &      - C4m(0,i1+1,i2+1,i3+1,3)
c     >      end do
c     >        C4m(0,1,1,i1+1,i2+1) = -C3m(0,1,i1+1,i2+1) 
c     >     &      - C4m(0,1,2,i1+1,i2+1) 
c     >     &      - C4m(0,1,i1+1,i2+1,3)
c     >      end do
c     >        C4m(0,1,1,1,i1+1) = -C3m(0,1,1,i1+1) 
c     >     &      - C4m(0,1,1,2,i1+1) 
c     >     &      - C4m(0,1,1,i1+1,3)
c     >      end do
c     >      C4m(0,1,1,1,1) = -C3m(0,1,1,1) 
c     >     &      - C4m(0,1,1,1,2) 
c     >     &      - C4m(0,1,1,1,3)

        do i1=1,3
          do i2=i1,3
            C4m(0,0,0,i2,i1) = C4m(0,0,0,i1,i2)
            C4m(1,0,0,i2,i1) = C4m(1,0,0,i1,i2)
            C4m(2,0,0,i2,i1) = C4m(2,0,0,i1,i2)
            C4m(3,0,0,i2,i1) = C4m(3,0,0,i1,i2)
            do i3=i2,3
              do i4=i3,3
                C4m(0,i1,i3,i2,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i1,i3,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i3,i1,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i1,i2,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i2,i1,i4) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i2,i4,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i4,i2,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i1,i4,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i4,i1,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i1,i2,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i2,i1,i3) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i3,i4,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i1,i4,i3,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i1,i4,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i4,i1,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i1,i3,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i3,i1,i2) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i2,i4,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i3,i4,i2,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i3,i4,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i2,i4,i3,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i3,i2,i1) = C4m(0,i1,i2,i3,i4)
                C4m(0,i4,i2,i3,i1) = C4m(0,i1,i2,i3,i4)

                C4m(1,i1,i3,i2,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i1,i3,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i3,i1,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i1,i2,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i2,i1,i4) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i2,i4,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i4,i2,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i1,i4,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i4,i1,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i1,i2,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i2,i1,i3) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i3,i4,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i1,i4,i3,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i1,i4,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i4,i1,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i1,i3,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i3,i1,i2) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i2,i4,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i3,i4,i2,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i3,i4,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i2,i4,i3,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i3,i2,i1) = C4m(1,i1,i2,i3,i4)
                C4m(1,i4,i2,i3,i1) = C4m(1,i1,i2,i3,i4)

                C4m(2,i1,i3,i2,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i1,i3,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i3,i1,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i1,i2,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i2,i1,i4) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i2,i4,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i4,i2,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i1,i4,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i4,i1,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i1,i2,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i2,i1,i3) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i3,i4,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i1,i4,i3,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i1,i4,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i4,i1,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i1,i3,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i3,i1,i2) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i2,i4,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i3,i4,i2,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i3,i4,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i2,i4,i3,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i3,i2,i1) = C4m(2,i1,i2,i3,i4)
                C4m(2,i4,i2,i3,i1) = C4m(2,i1,i2,i3,i4)

                C4m(3,i1,i3,i2,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i1,i3,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i3,i1,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i1,i2,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i2,i1,i4) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i2,i4,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i4,i2,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i1,i4,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i4,i1,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i1,i2,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i2,i1,i3) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i3,i4,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i1,i4,i3,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i1,i4,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i4,i1,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i1,i3,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i3,i1,i2) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i2,i4,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i3,i4,i2,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i3,i4,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i2,i4,i3,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i3,i2,i1) = C4m(3,i1,i2,i3,i4)
                C4m(3,i4,i2,i3,i1) = C4m(3,i1,i2,i3,i4)
              end do
            end do
          end do
        end do

        do i1=1,3
          S5hat(i1,0,0,0,0) = C4m(i1,0,0,0,0) - C4m(0,0,0,0,0)
          do i2=1,3
            do i3=1,3
              S5hat(i1,0,0,i2,i3) = C4m(i1,0,0,i2,i3) - C4m(0,0,0,i2,i3)
              do i4=1,3
                do k =1,3
                  S5hat(k,i1,i2,i3,i4) = C4m(k,i1,i2,i3,i4) -
     &                C4m(0,i1,i2,i3,i4)
                end do
              end do
            end do
          end do
      end do

 99     continue

c choose reduction formulas with biggest denominators

      maxf = abs(f(1))
      m = 1
      if (abs(f(2)).gt.maxf) then
        maxf = abs(f(2))
        m = 2
      end if
      if (abs(f(3)).gt.maxf) then
        maxf = abs(f(3))
        m = 3
      end if
        
c      write(testout,*) 'm',m
c      write(testout,*) 'rank   = ',rank
c      write(testout,*) 'ordgp4 = ',ordgp4
      do 100 r=0,ordgp4+rank
c        write(testout,*) 'r     = ',r
        
        D0 = ( S1hat(m) - z(m,1)*D1(1)  - z(m,2)*D1(2)
     &                  - z(m,3)*D1(3) )/f(m)
c        write(testout,9) 'D0  = ',
c     &  cD0f(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32)
c        write(testout,9) 'D0  = ',D0
 9      format(1x,a10,2g24.16,2i2)  
 10     format(1x,a10,i2,2g24.16,2i2)  
 20     format(1x,a10,2i2,2g24.16,2i2)  
 30     format(1x,a10,3i2,2g24.16)  
 40     format(1x,a10,4i2,2g24.16)  
 50     format(1x,a10,5i2,2g24.16)  


        if (rank+ordgp4-r.eq.0) goto 100

        S00 = 2d0*( mm02*D0 + C0m(0) )

        D2(0,0) = ( S00
     &         - z(1,1)*D2(1,1) - 2d0*z(1,2)*D2(1,2)
     &         - z(2,2)*D2(2,2) - 2d0*z(1,3)*D2(1,3)
     &         - z(3,3)*D2(3,3) - 2d0*z(2,3)*D2(2,3)
     &           ) /8d0

c        write(testout,9) 'D00  = ',D2(0,0)

        do i1=1,3
          do j =1,3
            S2mod(j,i1) = S2hat(j,i1)
          end do
          S2mod(i1,i1) = S2mod(i1,i1) - 2d0*D2(0,0)
        end do

        do i1=1,3
          D1(i1) = (S2mod(m,i1) 
     &        - z(m,1)*D2(1,i1)  - z(m,2)*D2(2,i1) - z(m,3)*D2(3,i1)  
     &        )/f(m)
c          write(testout,10) 'D1  = ',i1,D1(i1)
        end do

        if (rank+ordgp4-r.eq.1) goto 100

        do i1=1,3
          S001(i1) = 2d0*( mm02*D1(i1) + C1m(0,i1) )
        end do

        do i1=1,3
          D3(0,0,i1) =  ( S001(i1)
     &         - z(1,1)*D3(1,1,i1) - 2d0*z(1,2)*D3(1,2,i1)
     &         - z(2,2)*D3(2,2,i1) - 2d0*z(1,3)*D3(1,3,i1)
     &         - z(3,3)*D3(3,3,i1) - 2d0*z(2,3)*D3(2,3,i1)
     &           ) /12d0
c          write(testout,10) 'D00i  = ',i1,D3(0,0,i1)
        end do
        
        do i1=1,3
          do i2=1,3
            do j =1,3
              S3mod(j,i1,i2) = S3hat(j,i1,i2)
            end do
            S3mod(i1,i1,i2) = S3mod(i1,i1,i2) - 2d0*D3(0,0,i2)
            S3mod(i2,i1,i2) = S3mod(i2,i1,i2) - 2d0*D3(0,0,i1)
          end do
        end do
        
        do i1=1,3
          do i2=1,3
            D2(i1,i2) =  (S3mod(m,i1,i2) 
     &          - z(m,1)*D3(1,i1,i2)  - z(m,2)*D3(2,i1,i2) 
     &          - z(m,3)*D3(3,i1,i2)  
     &        )/f(m)
c            write(testout,20) 'D2  = ',i1,i2,D2(i1,i2)
          end do
        end do

c     >      if (sym.eq.1) then
c     >        D2(2,1) = D2(1,2)
c     >        D2(3,1) = D2(1,3)
c     >        D2(3,2) = D2(2,3)
c     >      end if

        if (rank+ordgp4-r.eq.2) goto 100
      
        do i1=1,3
          do i2=1,3
            S002(i1,i2) = 2d0*( mm02*D2(i1,i2) + C2m(0,i1,i2) )
          end do
        end do

        do i1=1,3
        do i2=1,3
          D4(0,0,i1,i2) =  ( S002(i1,i2)
     &         - z(1,1)*D4(1,1,i1,i2) - 2d0*z(1,2)*D4(1,2,i1,i2)
     &         - z(2,2)*D4(2,2,i1,i2) - 2d0*z(1,3)*D4(1,3,i1,i2)
     &         - z(3,3)*D4(3,3,i1,i2) - 2d0*z(2,3)*D4(2,3,i1,i2)
     &           ) /16d0
c          write(testout,20) 'D00ij  = ',i1,i2,D4(0,0,i1,i2)
        end do
        end do
                
        do i1=1,3
          do i2=1,3
            do i3=1,3
              do j =1,3
                S4mod(j,i1,i2,i3) = S4hat(j,i1,i2,i3)
              end do
              S4mod(i1,i1,i2,i3) = S4mod(i1,i1,i2,i3) -2d0*D4(0,0,i2,i3)
              S4mod(i2,i1,i2,i3) = S4mod(i2,i1,i2,i3) -2d0*D4(0,0,i1,i3)
              S4mod(i3,i1,i2,i3) = S4mod(i3,i1,i2,i3) -2d0*D4(0,0,i1,i2)
            end do
          end do
        end do
        
        do i1=1,3
          do i2=1,3
            do i3=1,3
              D3(i1,i2,i3) =  (S4mod(m,i1,i2,i3) 
     &          - z(m,1)*D4(1,i1,i2,i3)  - z(m,2)*D4(2,i1,i2,i3) 
     &          - z(m,3)*D4(3,i1,i2,i3)  
     &        )/f(m)
c              write(testout,30) 'D3  = ',i1,i2,i3,D3(i1,i2,i3)
            end do
          end do
        end do

c     >      if (sym.eq.1) then
c     >        D3(2,1,1) = D3(1,1,2)
c     >        D3(1,2,1) = D3(1,1,2)
c     >        D3(3,1,1) = D3(1,1,3)
c     >        D3(1,3,1) = D3(1,1,3)
c     >        D3(2,1,2) = D3(1,2,2)
c     >        D3(2,2,1) = D3(1,2,2)
c     >        D3(3,1,3) = D3(1,3,3)
c     >        D3(3,3,1) = D3(1,3,3)
c     >        D3(3,3,2) = D3(2,3,3)
c     >        D3(3,2,3) = D3(2,3,3)
c     >        D3(3,2,2) = D3(2,2,3)
c     >        D3(2,3,2) = D3(2,2,3)
c     >        D3(2,1,3) = D3(1,2,3)
c     >        D3(1,3,2) = D3(1,2,3)
c     >        D3(3,1,2) = D3(1,2,3)
c     >        D3(2,3,1) = D3(1,2,3)
c     >        D3(3,2,1) = D3(1,2,3)
c     >      end if

        if (rank+ordgp4-r.eq.3) goto 100

        S002(0,0) = 2d0*( mm02*D2(0,0) + C2m(0,0,0) )

        D4(0,0,0,0) = ( S002(0,0) + 1d0/6d0 
     &         - z(1,1)*D4(0,0,1,1) - 2d0*z(1,2)*D4(0,0,1,2)
     &         - z(2,2)*D4(0,0,2,2) - 2d0*z(1,3)*D4(0,0,1,3)
     &         - z(3,3)*D4(0,0,3,3) - 2d0*z(2,3)*D4(0,0,2,3)
     &           ) /12d0
c        write(testout,9) 'D0000  = ',D4(0,0,0,0)

      
        do i1=1,3
          do i2=1,3
            do i3=1,3
              S003(i1,i2,i3) =
     &             2d0*( mm02*D3(i1,i2,i3) + C3m(0,i1,i2,i3) )
            end do
          end do
        end do

        do i1=1,3
          do i2=1,3
            do i3=1,3
              D5(0,0,i1,i2,i3) =  ( S003(i1,i2,i3)
c     &         - z(1,1)*D5(1,1,i1,i2,i3) - 2d0*z(1,2)*D5(1,2,i1,i2,i3)
c     &         - z(2,2)*D5(2,2,i1,i2,i3) - 2d0*z(1,3)*D5(1,3,i1,i2,i3)
c     &         - z(3,3)*D5(3,3,i1,i2,i3) - 2d0*z(2,3)*D5(2,3,i1,i2,i3)
     &           ) /20d0
c          write(testout,30) 'D00ijk = ',i1,i2,i3,D5(0,0,i1,i2,i3)
            end do
          end do
        end do

c     >      if (sym.eq.1) then
c     >        D5(0,0,2,1,1) = D5(0,0,1,1,2)
c     >        D5(0,0,1,2,1) = D5(0,0,1,1,2)
c     >        D5(0,0,3,1,1) = D5(0,0,1,1,3)
c     >        D5(0,0,1,3,1) = D5(0,0,1,1,3)
c     >        D5(0,0,2,1,2) = D5(0,0,1,2,2)
c     >        D5(0,0,2,2,1) = D5(0,0,1,2,2)
c     >        D5(0,0,3,1,3) = D5(0,0,1,3,3)
c     >        D5(0,0,3,3,1) = D5(0,0,1,3,3)
c     >        D5(0,0,3,3,2) = D5(0,0,2,3,3)
c     >        D5(0,0,3,2,3) = D5(0,0,2,3,3)
c     >        D5(0,0,3,2,2) = D5(0,0,2,2,3)
c     >        D5(0,0,2,3,2) = D5(0,0,2,2,3)
c     >        D5(0,0,2,1,3) = D5(0,0,1,2,3)
c     >        D5(0,0,1,3,2) = D5(0,0,1,2,3)
c     >        D5(0,0,3,1,2) = D5(0,0,1,2,3)
c     >        D5(0,0,2,3,1) = D5(0,0,1,2,3)
c     >        D5(0,0,3,2,1) = D5(0,0,1,2,3)
c     >      end if

        do i1=1,3
          do i2=1,3
            do i3=1,3
              do i4=1,3
                do j =1,3
                  S5mod(j,i1,i2,i3,i4) = S5hat(j,i1,i2,i3,i4)
                end do
                S5mod(i1,i1,i2,i3,i4) = S5mod(i1,i1,i2,i3,i4) 
     &              - 2d0*D5(0,0,i2,i3,i4)
                S5mod(i2,i1,i2,i3,i4) = S5mod(i2,i1,i2,i3,i4) 
     &              - 2d0*D5(0,0,i1,i3,i4)
                S5mod(i3,i1,i2,i3,i4) = S5mod(i3,i1,i2,i3,i4) 
     &              - 2d0*D5(0,0,i1,i2,i4)
                S5mod(i4,i1,i2,i3,i4) = S5mod(i4,i1,i2,i3,i4) 
     &              - 2d0*D5(0,0,i1,i2,i3)
              end do
            end do
          end do
        end do

        do i1=1,3
          do i2=1,3
            do i3=1,3
              do i4=1,3
                  D4(i1,i2,i3,i4) =  (S5mod(m,i1,i2,i3,i4) 
c     &          - z(m,1)*D5(1,i1,i2,i3,i4)  - z(m,2)*D5(2,i1,i2,i3,i4) 
c     &          - z(m,3)*D5(3,i1,i2,i3,i4)  
     &        )/f(m)
c                 write(testout,40) 'D4  = ',i1,i2,i3,i4,D4(i1,i2,i3,i4)
              end do
            end do
          end do
        end do

c D5 incomplete!!!
       if (rank+ordgp3-r.eq.4) goto 100

        do i1=1,3
          S003(0,0,i1) = 2d0*( mm02*D3(0,0,i1) + C3m(0,0,0,i1) )
        end do

        do i1=1,3
          D5(0,0,0,0,i1) = ( S003(0,0,i1) - 1d0/24d0 
     &         - z(1,1)*D5(0,0,1,1,i1) - 2d0*z(1,2)*D5(0,0,1,2,i1)
     &         - z(2,2)*D5(0,0,2,2,i1) - 2d0*z(1,3)*D5(0,0,1,3,i1)
     &         - z(3,3)*D5(0,0,3,3,i1) - 2d0*z(2,3)*D5(0,0,2,3,i1)
     &           ) /12d0
c          write(testout,10) 'D0000i  = ',i1,D5(0,0,0,0,i1)
        end do

 100  continue

      end



c cache

***********************************************************************
*                                                                     *
*     fast evaluation of tensor integrals                             *
*                                                                     *
*     written by Markus Roth,          3.01.05                        *
*     adapted by Ansgar Denner,       18.01.05                        *
*     last modified by Ansgar Denner, 31.01.05                        *
*                                                                     *
* subroutines:                                                        *
* cacheon, cacheoff, setcachelevel, cachetempoff, cachereon           *
* cacheinit, cacheread, cachewrite                                    *
***********************************************************************
      subroutine cacheon
***********************************************************************
*     switch cache for loop integrals on                              *
*---------------------------------------------------------------------*
*     20.01.05  Ansgar Denner     last changed  31.01.05              *
***********************************************************************
      implicit none
c fast
      integer maxtype,maxcalc,maxcall,maxfct,maxx
c      parameter(maxtype=5,maxcalc=500,maxcall=2000,maxfct=86,maxx=15)
      parameter(maxtype=5,maxcalc=580,maxcall=2400,maxfct=86,maxx=15)
c      common  /cachelim/  maxtype,maxcalc,maxcall,maxfct,maxx
      complex*16 f(maxtype,maxcalc,maxfct),arg(maxtype,maxcalc,maxx)
      integer grank(maxtype,maxcalc),gswitch(maxtype,maxcalc)
      integer lrank(maxtype,maxcalc),lswitch(maxtype,maxcalc)
      integer ncalc(maxtype),calc(maxtype,maxcall)
      integer pointer(maxtype,maxcall),ncall(maxtype),n
      character*11 name(maxtype)
      integer usecache,usecachesave,cachelevel,cacheleveli
      common/cache/f,arg,grank,lrank,gswitch,lswitch,
     &             ncalc,calc,pointer,ncall,n,name
      common/cachem/usecache,usecachesave,cachelevel,cacheleveli

      usecache = 1
      usecachesave = 1
c      write(*,*) 'cacheon'
      end
***********************************************************************
      subroutine cacheoff
***********************************************************************
*     switch cache for loop integrals off                             *
*---------------------------------------------------------------------*
*     20.01.05  Ansgar Denner     last changed  31.01.05              *
***********************************************************************
      implicit none
c fast
      integer maxtype,maxcalc,maxcall,maxfct,maxx
c      parameter(maxtype=5,maxcalc=500,maxcall=2000,maxfct=86,maxx=15)
      parameter(maxtype=5,maxcalc=580,maxcall=2400,maxfct=86,maxx=15)
c      common  /cachelim/  maxtype,maxcalc,maxcall,maxfct,maxx
      complex*16 f(maxtype,maxcalc,maxfct),arg(maxtype,maxcalc,maxx)
      integer grank(maxtype,maxcalc),gswitch(maxtype,maxcalc)
      integer lrank(maxtype,maxcalc),lswitch(maxtype,maxcalc)
      integer ncalc(maxtype),calc(maxtype,maxcall)
      integer pointer(maxtype,maxcall),ncall(maxtype),n
      character*11 name(maxtype)
      integer usecache,usecachesave,cachelevel,cacheleveli
      common/cache/f,arg,grank,lrank,gswitch,lswitch,
     &             ncalc,calc,pointer,ncall,n,name
      common/cachem/usecache,usecachesave,cachelevel,cacheleveli

      usecache = 0
      usecachesave = 0
      end
***********************************************************************
      subroutine setcachelevel(cachelevelin)
***********************************************************************
*     switch level for loop integral cache                            *
*     E = 1, D = 2, C = 3, B = 4                                      *
*---------------------------------------------------------------------*
*     20.01.05  Ansgar Denner     last changed  31.01.05              *
***********************************************************************
      implicit none
      integer  cachelevelin
c fast
      integer maxtype,maxcalc,maxcall,maxfct,maxx
c      parameter(maxtype=5,maxcalc=500,maxcall=2000,maxfct=86,maxx=15)
      parameter(maxtype=5,maxcalc=580,maxcall=2400,maxfct=86,maxx=15)
c      common  /cachelim/  maxtype,maxcalc,maxcall,maxfct,maxx
      complex*16 f(maxtype,maxcalc,maxfct),arg(maxtype,maxcalc,maxx)
      integer grank(maxtype,maxcalc),gswitch(maxtype,maxcalc)
      integer lrank(maxtype,maxcalc),lswitch(maxtype,maxcalc)
      integer ncalc(maxtype),calc(maxtype,maxcall)
      integer pointer(maxtype,maxcall),ncall(maxtype),n
      character*11 name(maxtype)
      integer usecache,usecachesave,cachelevel,cacheleveli
      common/cache/f,arg,grank,lrank,gswitch,lswitch,
     &             ncalc,calc,pointer,ncall,n,name
      common/cachem/usecache,usecachesave,cachelevel,cacheleveli

      cachelevel = cachelevelin
      cacheleveli = cachelevelin
      end
***********************************************************************
      subroutine cachetempoff
***********************************************************************
*     switch cache for loop integrals temporarily off                 *
*---------------------------------------------------------------------*
*     20.01.05  Ansgar Denner     last changed  31.01.05              *
***********************************************************************
      implicit none
c fast
      integer maxtype,maxcalc,maxcall,maxfct,maxx
c      parameter(maxtype=5,maxcalc=500,maxcall=2000,maxfct=86,maxx=15)
      parameter(maxtype=5,maxcalc=580,maxcall=2400,maxfct=86,maxx=15)
c      common  /cachelim/  maxtype,maxcalc,maxcall,maxfct,maxx
      complex*16 f(maxtype,maxcalc,maxfct),arg(maxtype,maxcalc,maxx)
      integer grank(maxtype,maxcalc),gswitch(maxtype,maxcalc)
      integer lrank(maxtype,maxcalc),lswitch(maxtype,maxcalc)
      integer ncalc(maxtype),calc(maxtype,maxcall)
      integer pointer(maxtype,maxcall),ncall(maxtype),n
      character*11 name(maxtype)
      integer usecache,usecachesave,cachelevel,cacheleveli
      common/cache/f,arg,grank,lrank,gswitch,lswitch,
     &             ncalc,calc,pointer,ncall,n,name
      common/cachem/usecache,usecachesave,cachelevel,cacheleveli

      usecache = 0
c      write(*,*) 'cacheon'
      end
***********************************************************************
      subroutine cachereon
***********************************************************************
*     switch cache for loop integrals back on                         *
*---------------------------------------------------------------------*
*     20.01.05  Ansgar Denner     last changed  31.01.05              *
***********************************************************************
      implicit none
c fast
      integer maxtype,maxcalc,maxcall,maxfct,maxx
c      parameter(maxtype=5,maxcalc=500,maxcall=2000,maxfct=86,maxx=15)
      parameter(maxtype=5,maxcalc=580,maxcall=2400,maxfct=86,maxx=15)
c      common  /cachelim/  maxtype,maxcalc,maxcall,maxfct,maxx
      complex*16 f(maxtype,maxcalc,maxfct),arg(maxtype,maxcalc,maxx)
      integer grank(maxtype,maxcalc),gswitch(maxtype,maxcalc)
      integer lrank(maxtype,maxcalc),lswitch(maxtype,maxcalc)
      integer ncalc(maxtype),calc(maxtype,maxcall)
      integer pointer(maxtype,maxcall),ncall(maxtype),n
      character*11 name(maxtype)
      integer usecache,usecachesave,cachelevel,cacheleveli
      common/cache/f,arg,grank,lrank,gswitch,lswitch,
     &             ncalc,calc,pointer,ncall,n,name
      common/cachem/usecache,usecachesave,cachelevel,cacheleveli

      usecache = usecachesave
c      write(*,*) 'cacheon'
      end

***********************************************************************
      subroutine cacheinit(nout)
***********************************************************************
*     initialization of fast evaluation of tensor integrals           *
*                                                                     *
*     written by Markus Roth, 3.1.2005                                *
*     adapted by Ansgar Denner                                        *
*---------------------------------------------------------------------*
*     20.01.05  Ansgar Denner     last changed  31.01.05              *
***********************************************************************
      implicit none
c local variables
      integer i1,i2,nout
c fast
      integer maxtype,maxcalc,maxcall,maxfct,maxx
c      parameter(maxtype=5,maxcalc=500,maxcall=2000,maxfct=86,maxx=15)
      parameter(maxtype=5,maxcalc=580,maxcall=2400,maxfct=86,maxx=15)
      complex*16 f(maxtype,maxcalc,maxfct),arg(maxtype,maxcalc,maxx)
      integer grank(maxtype,maxcalc),gswitch(maxtype,maxcalc)
      integer lrank(maxtype,maxcalc),lswitch(maxtype,maxcalc)
      integer ncalc(maxtype),calc(maxtype,maxcall)
      integer pointer(maxtype,maxcall),ncall(maxtype),n
      character*11 name(maxtype)
      integer usecache,usecachesave,cachelevel,cacheleveli
      common/cache/f,arg,grank,lrank,gswitch,lswitch,
     &             ncalc,calc,pointer,ncall,n,name
      common/cachem/usecache,usecachesave,cachelevel,cacheleveli

      common/cacheh/ncalc2,ncalc3
      integer ncalc2(maxtype),ncalc3(maxtype)

      data usecache/0/,usecachesave/0/,
     &     cachelevel /20/, cacheleveli /20/, n /0/

      integer countC
      common /countC/ countC
      integer countD
      common /countD/ countD
      integer countE
      common /countE/ countE

c initialization for the first call
      if(n.eq.0)then
        do i1=1,maxtype
          ncalc(i1)=0
          do i2=1,maxcall
            calc(i1,i2)=0
          enddo
        enddo
      endif
c output of initialization
      if(nout.ne.0.and.n.le.maxtype.and.usecache.gt.0)then
        if(n.eq.0)then
          write(nout,'(a)')' '
          write(nout,'(a)')
     *      ' Tensor integral      Type        Calls     Calculations'
        elseif(ncalc(n).ne.0)then
          write(nout,'(1x,a11,3(4x,i9))')name(n),n,ncall(n),ncalc(n)
        endif 
      endif

c initialization for each event
      if(usecache.gt.0) then
        if (n.le.maxtype+10) then
          do i1=max(n-10,1),min(maxtype,n-1)
            do i2=1,maxcalc
              if (lrank(i1,i2).lt.grank(i1,i2)) then
                grank(i1,i2) = lrank(i1,i2)
                gswitch(i1,i2) = lswitch(i1,i2)
              else if (lrank(i1,i2).eq.grank(i1,i2)
     &                 .and.lswitch(i1,i2).eq.0) then
                gswitch(i1,i2) = 0
              end if   
            end do
          end do
c>          do i1=max(n-10,1),min(maxtype,n-1)
c>            do i2=1,maxcalc
c>              if (n.eq.1) then
c>                write(11,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              else if (n.eq.2) then
c>                write(12,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              else if (n.eq.3) then
c>                write(13,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              else if (n.eq.4) then
c>                write(14,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              else if (n.eq.5) then
c>                write(15,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              else if (n.eq.6) then
c>                write(16,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              else if (n.eq.7) then
c>                write(17,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              else if (n.eq.8) then
c>                write(18,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              else if (n.eq.9) then
c>                write(19,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              else if (n.eq.10) then
c>                write(20,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              else if (n.eq.11) then
c>                write(21,*) i1,i2, grank(i1,i2),gswitch(i1,i2)
c>              end if
c>            end do
c>          end do
        end if

        n=n+1

c      write(*,*) 'cacheinit:      n = ',n
c      write(*,*) 'cacheinit:     ncalc = ',ncalc(1),ncalc(2),ncalc(3),ncalc(4)
c      write(*,*) 'cacheinit:     ncall = ',ncall(1),ncall(2),ncall(3),ncall(4)

        do i1=1,maxtype 
          ncall(i1)=0
c        ncalc2(i1)=0
c        ncalc3(i1)=0
        end do

      end if

      countC = 0
      countD = 0
      countE = 0

      end

***********************************************************************
      subroutine cacheread(fct,x,nfct,nx,type,nrank,nswitch,
     &    namefct,nocalc)
***********************************************************************
*     fast evaluation of tensor integrals                             *
*                                                                     *
*     written by Markus Roth, 3.1.2005                                *
*     adapted by Ansgar Denner                                        *
*---------------------------------------------------------------------*
*     20.01.05  Ansgar Denner     last changed  31.01.05              *
***********************************************************************
      implicit none
c local variables
      integer maxtype,maxcalc,maxcall,maxfct,maxx
c      common  /cachelim/  maxtype,maxcalc,maxcall,maxfct,maxx
c      parameter(maxtype=5,maxcalc=500,maxcall=2000,maxfct=86,maxx=15)
      parameter(maxtype=5,maxcalc=580,maxcall=2400,maxfct=86,maxx=15)
      complex*16 fct(maxfct),x(maxx)
      integer i1,i2,nfct,nx,type,nrank,nswitch
      logical nocalc,same
      character*11 namefct
c fast
      complex*16 f(maxtype,maxcalc,maxfct),arg(maxtype,maxcalc,maxx)
      integer grank(maxtype,maxcalc),gswitch(maxtype,maxcalc)
      integer lrank(maxtype,maxcalc),lswitch(maxtype,maxcalc)
      integer ncalc(maxtype),calc(maxtype,maxcall)
      integer pointer(maxtype,maxcall),ncall(maxtype),n
      character*11 name(maxtype)
      integer usecache,usecachesave,cachelevel,cacheleveli
      common/cache/f,arg,grank,lrank,gswitch,lswitch,
     &             ncalc,calc,pointer,ncall,n,name
      common/cachem/usecache,usecachesave,cachelevel,cacheleveli

c      common/cacheh/ncalc2,ncalc3
c      integer ncalc2(maxtype),ncalc3(maxtype)

      nocalc=.false.
c      write(*,*) 'cacheread in ',
c     &    usecache,usecachesave,cachelevel,cacheleveli,
c     &    n,type,cachelevel.lt.type,n.lt.type
      if(cachelevel.lt.type.or.n.lt.type.or.usecache.eq.0)return
c      write(*,*)  'cacheread2 ',n,type,ncall(type)+1
      if(n.gt.type+10)then
c use full optimization
        ncall(type)=ncall(type)+1
c        if (pointer(type,ncall(type)).eq.0) then
c        write(*,*) 'cacheread1 ',type,n,ncall(type),calc(type,ncall(type))
c        write(*,*) 'cacheread2 ',nrank,pointer(type,ncall(type))
c        write(*,*) 'cacheread3 ',usecache,usecachesave,cachelevel,cacheleveli
c        write(*,*) 'cacheread4 ',grank(type,pointer(type,ncall(type)))
c        end if

        if(calc(type,ncall(type)).eq.1)then
c calculation of function
c          ncalc2(type) = ncalc2(type) + 1
c          ncalc3(type) = ncalc3(type) + 1
c          if (pointer(type,ncall(type)).ne.ncalc2(type)) then
c            write(*,*)  'cache error',
c     &          pointer(type,ncall(type)),ncalc2(type)
c          end if
          if (nrank.eq.grank(type,pointer(type,ncall(type)))) then
            if (nswitch.eq.0) then
              nswitch =  gswitch(type,pointer(type,ncall(type)))
            else if (nswitch.eq.
     &           -gswitch(type,pointer(type,ncall(type)))) then
              write(*,*) 'increase of nrank should not happen! a '
              nrank = nrank + 1
              nswitch = 0
            end if
          else if (nrank.lt.grank(type,pointer(type,ncall(type)))) then
            nrank = grank(type,pointer(type,ncall(type)))
            nswitch = gswitch(type,pointer(type,ncall(type)))
          end if 
          lrank(type,pointer(type,ncall(type))) = nrank
          lswitch(type,pointer(type,ncall(type))) = nswitch
c     check whether calculated rank is sufficient
        else if (nrank.lt.lrank(type,pointer(type,ncall(type))).or.
     &       nrank.eq.lrank(type,pointer(type,ncall(type)))
     &        .and.(nswitch.eq.lswitch(type,pointer(type,ncall(type)))
     &        .or.nswitch.eq.0)) then
c using already calculated result to fasten computation
          nocalc=.true.
          do i1=1,nfct
            fct(i1)=f(type,pointer(type,ncall(type)),i1)
          enddo
        else
c recalculate with cache off
          cacheleveli = cachelevel
          cachelevel = type
c          write(*,*) 'cachelevel changed to ',cachelevel
c          ncalc3(type) = ncalc3(type) + 1
          if (nrank.eq.lrank(type,pointer(type,ncall(type)))) then
            if (lswitch(type,pointer(type,ncall(type))).eq.0) then
              lswitch(type,pointer(type,ncall(type))) = nswitch
            else 
              write(*,*) 'increase of nrank should not happen! b '
              nrank = nrank + 1
              nswitch = 0
              lrank(type,pointer(type,ncall(type))) = nrank
              lswitch(type,pointer(type,ncall(type))) = 0
            end if
          endif
        end if
      else if(n.gt.type)then
c try to minimize grank and gswitch
c not yet full optimization
        ncall(type)=ncall(type)+1
c        if (pointer(type,ncall(type)).eq.0) then
c        write(*,*) 'cacheread1 ',type,n,ncall(type),calc(type,ncall(type))
c        write(*,*) 'cacheread2 ',nrank,pointer(type,ncall(type))
c        write(*,*) 'cacheread3 ',usecache,usecachesave,cachelevel,cacheleveli
c        write(*,*) 'cacheread4 ',grank(type,pointer(type,ncall(type)))
c        end if

        if(calc(type,ncall(type)).eq.1)then
c calculation of function
          lrank(type,pointer(type,ncall(type))) = nrank
          lswitch(type,pointer(type,ncall(type))) = nswitch
c          ncalc2(type) = ncalc2(type) + 1
c          ncalc3(type) = ncalc3(type) + 1
c          if (pointer(type,ncall(type)).ne.ncalc2(type)) then
c            write(*,*)  'cache error',
c     &          pointer(type,ncall(type)),ncalc2(type)
c          end if
          if (nrank.eq.grank(type,pointer(type,ncall(type)))) then
            if (nswitch.eq.0) then
              nswitch =  gswitch(type,pointer(type,ncall(type)))
            else if (nswitch.eq.
     &           -gswitch(type,pointer(type,ncall(type)))) then
              write(*,*) 'increase of nrank should not happen! c '
              nrank = nrank + 1
              nswitch = 0
            end if
          else if (nrank.lt.grank(type,pointer(type,ncall(type)))) then
            nrank = grank(type,pointer(type,ncall(type)))
            nswitch = gswitch(type,pointer(type,ncall(type)))
          end if 
        else 
c store maximal rank to lrank
          if (nrank.gt.lrank(type,pointer(type,ncall(type)))) then
            lrank(type,pointer(type,ncall(type))) = nrank
            lswitch(type,pointer(type,ncall(type))) = nswitch
          else if (nrank.eq.lrank(type,pointer(type,ncall(type)))) then
            if (nswitch.ne.lswitch(type,pointer(type,ncall(type)))
     &          .and.nswitch.ne.0) then
              if (lswitch(type,pointer(type,ncall(type))).eq.0) then
                lswitch(type,pointer(type,ncall(type))) = nswitch
              else 
                lrank(type,pointer(type,ncall(type))) = 
     &              lrank(type,pointer(type,ncall(type))) + 1
                lswitch(type,pointer(type,ncall(type))) = 0
              end if
            end if
          end if              
c     check whether calculated rank is sufficient
          if (nrank.lt.grank(type,pointer(type,ncall(type))).or.
     &       nrank.eq.grank(type,pointer(type,ncall(type)))
     &        .and.(nswitch.eq.gswitch(type,pointer(type,ncall(type)))
     &        .or.nswitch.eq.0)) then
c using already calculated result to fasten computation
            nocalc=.true.
            do i1=1,nfct
              fct(i1)=f(type,pointer(type,ncall(type)),i1)
            enddo
          else
c recalculate with cache off
            cacheleveli = cachelevel
            cachelevel = type
c     write(*,*) 'cachelevel changed to ',cachelevel
c            ncalc3(type) = ncalc3(type) + 1
            if (nrank.eq.grank(type,pointer(type,ncall(type)))) then
              if (gswitch(type,pointer(type,ncall(type))).ne.0) then
                write(*,*) 'increase of nrank should not happen! d '
                write(*,*) 'type,rank ',type,nrank,
     &              grank(type,pointer(type,ncall(type)))
                write(*,*) 'switch ',nswitch,
     &              gswitch(type,pointer(type,ncall(type)))
                nrank = nrank + 1
                nswitch = 0
              end if
            else
              lrank(type,pointer(type,ncall(type))) = nrank
              lswitch(type,pointer(type,ncall(type))) = nswitch
            endif
          end if
        end if
      elseif(n.eq.type)then
c      write(*,*)  'cacheread n=type'
c     error messages
        if(type.gt.maxtype)then
          write(*,'(a)')' fastin: maxtype too small'
          stop
        endif
        if(ncall(type)+1.gt.maxcall)then
          write(*,'(a)')' fastin: maxcall too small'
          write(*,*) type,ncall(type),maxcall
          stop
        endif
        if(ncalc(type)+1.gt.maxcalc)then
          write(*,'(a)')' fastin: maxcalc too small'
          stop
        endif
        if(nfct.gt.maxfct)then
          write(*,'(a)')' fastin: maxfct too small'
          stop
        endif
        if(nx.gt.maxx)then
          write(*,'(a)')' fastin: maxx too small'
          stop
        endif
c initialization for function type
        ncall(type)=ncall(type)+1
        do i1=1,ncalc(type)
          same=.true.
          do i2=1,nx
            if(arg(type,i1,i2).ne.x(i2))same=.false.
          enddo
          if(same)then
            pointer(type,ncall(type))=i1
            calc(type,ncall(type))=0
c            grank(type,i1)=max(grank(type,i1),nrank)
            if (nrank.gt.grank(type,i1)) then
              grank(type,i1) = nrank
              gswitch(type,i1) = nswitch
            else if (nrank.eq.grank(type,i1)) then
              if (nswitch.ne.gswitch(type,i1).and.nswitch.ne.0) then
                if (gswitch(type,i1).eq.0) then
                  gswitch(type,i1) = nswitch
                else 
                  write(*,*) 'increase of nrank should not happen! e '
                  nrank = nrank + 1
                  nswitch = 0
                  grank(type,i1) = nrank
                  gswitch(type,i1) = 0
                end if
              end if
            end if              
            return
          endif
        enddo
        ncalc(type)=ncalc(type)+1
        pointer(type,ncall(type))=ncalc(type)
        calc(type,ncall(type))=1
        grank(type,ncalc(type))=nrank
        gswitch(type,ncalc(type))=nswitch
        name(type)=namefct
        do i1=1,nx
          arg(type,ncalc(type),i1)=x(i1)
        enddo
      endif

      return

      write(*,10)
      write(*,12) 'cacheread: usecache =',usecache,usecachesave
     &                     ,cachelevel,cacheleveli
      write(*,10) 'cacheread:      n = ',n
      write(*,10) 'cacheread:   type = ',type
      write(*,10) 'cacheread:  ncalc = ',ncalc(type)
c      write(*,10) 'cacheread: ncalc3 = ',ncalc3(type)
c      write(*,10) 'cacheread: ncalc2 = ',ncalc2(type)
      write(*,10) 'cacheread:  ncall = ',ncall(type)
      write(*,10) 'cacheread:     pt = ',pointer(type,ncall(type))
      write(*,11) 'cacheread:   name = ',name(type)
      if (n.ge.type) then
      write(*,10) 'cacheread:  grank = ',grank(type,ncalc(type))
      write(*,10) 'cacheread: gswitch= ',gswitch(type,ncalc(type))
      write(*,10) 'cacheread:  nrank = ',nrank
      end if
 10   format(A21,i4)
 11   format(A21,a11)
 12   format(A22,3i4)
 
      end

***********************************************************************
      subroutine cachewrite(fct,nfct,type)
***********************************************************************
*     fast evaluation of tensor integrals                             *
*                                                                     *
*     written by Markus Roth, 3.1.2005                                *
*     adapted by Ansgar Denner                                        *
*---------------------------------------------------------------------*
*     20.01.05  Ansgar Denner     last changed  31.01.05              *
***********************************************************************
      implicit none
c local variables
      integer maxtype,maxcalc,maxcall,maxfct,maxx
c      parameter(maxtype=5,maxcalc=500,maxcall=2000,maxfct=86,maxx=15)
      parameter(maxtype=5,maxcalc=580,maxcall=2400,maxfct=86,maxx=15)
c      common  /cachelim/  maxtype,maxcalc,maxcall,maxfct,maxx
      complex*16 fct(maxfct)
      integer i1,nfct,type
c fast
      complex*16 f(maxtype,maxcalc,maxfct),arg(maxtype,maxcalc,maxx)
      integer grank(maxtype,maxcalc),gswitch(maxtype,maxcalc)
      integer lrank(maxtype,maxcalc),lswitch(maxtype,maxcalc)
      integer ncalc(maxtype),calc(maxtype,maxcall)
      integer pointer(maxtype,maxcall),ncall(maxtype),n
      character*11 name(maxtype)
      integer usecache,usecachesave,cachelevel,cacheleveli
      common/cache/f,arg,grank,lrank,gswitch,lswitch,
     &             ncalc,calc,pointer,ncall,n,name
      common/cachem/usecache,usecachesave,cachelevel,cacheleveli

c      common/cacheh/ncalc2,ncalc3
c      integer ncalc2(maxtype),ncalc3(maxtype)


      if(cachelevel.lt.type.or.n.lt.type.or.usecache.eq.0)return
      if(calc(type,ncall(type)).eq.1)then
c storing result of function type
        do i1=1,nfct
          f(type,pointer(type,ncall(type)),i1)=fct(i1)
        enddo
      endif
c      if (usecache.lt.usecachesavei) write(*,*) 'usecache reset to',
c     &      usecachesavei
 
      cachelevel = cacheleveli

      return

      write(*,10)
      write(*,12) 'cachewrite: usecache = ',usecache,usecachesave,
     &                   cachelevel,cacheleveli
      write(*,10) 'cachewrite:      n = ',n
      write(*,10) 'cachewrite:   type = ',type
      write(*,10) 'cachewrite:  ncalc = ',ncalc(type)
c      write(*,10) 'cachewrite: ncalc3 = ',ncalc3(type)
c      write(*,10) 'cachewrite: ncalc2 = ',ncalc2(type)
      write(*,10) 'cachewrite:  ncall = ',ncall(type)
      write(*,10) 'cachewrite:     pt = ',pointer(type,ncall(type))
c      write(*,11) 'cachewrite:   name = ',name(type)
      write(*,10) 'cachewrite:  grank = ',grank(type,ncalc(type))
      write(*,10) 'cachewrite: gswitch= ',gswitch(type,ncalc(type))
 10   format(A21,i4)
 11   format(A21,a11)
 12   format(A22,3i4)
      end



***********************************************************************
      subroutine countinit
***********************************************************************
*     initialize counts for different branches in C and D reduction   *
*---------------------------------------------------------------------*
*     23.03.05  Ansgar Denner     last changed  24.03.05              *
***********************************************************************
      implicit   none

      integer    Ccount(0:40)
      integer    Dcount(0:40)
      integer    i

      common /Ccount/ Ccount
      common /Dcount/ Dcount

      real*8     accbad1,accbad2,accbad3
      common /accbad/ accbad1,accbad2,accbad3

      accbad1  = 1d-4
      accbad2  = 1d-2
      accbad3  = 1d0
      do i=0,40
        Ccount(i) = 0
        Dcount(i) = 0
      end do 

      end 

***********************************************************************
      subroutine writecount(out)
***********************************************************************
*     output counts for different branches in C and D reduction       *
*---------------------------------------------------------------------*
*     23.03.05  Ansgar Denner     last changed  23.05.05              *
*                            cosmetic changes 01.06.06 Ansgar Denner  *
***********************************************************************
      implicit   none
      integer    out
      integer    Ccount(0:40)
      integer    Dcount(0:40)
      integer    Ccount1,Ccount2,Ccount3,Ccount4
      integer    Dcount1,Dcount2,Dcount3,Dcount4
      integer    i
      real*8     accbad1,accbad2,accbad3

      common /Ccount/ Ccount
      common /Dcount/ Dcount
      common /accbad/ accbad1,accbad2,accbad3

      Ccount1 = 0
      Ccount2 = 0
      Ccount3 = 0
      Ccount4 = 0
      Dcount1 = 0
      Dcount2 = 0
      Dcount3 = 0
      Dcount4 = 0
      do i=1,9
        Ccount1 = Ccount1 + Ccount(i)
        Ccount2 = Ccount2 + Ccount(i+10)
        Ccount3 = Ccount3 + Ccount(i+20)
        Ccount4 = Ccount4 + Ccount(i+30)
        Dcount1 = Dcount1 + Dcount(i)
        Dcount2 = Dcount2 + Dcount(i+10)
        Dcount3 = Dcount3 + Dcount(i+20)
        Dcount4 = Dcount4 + Dcount(i+30)
      end do 
      Ccount1 = Ccount1 - Ccount(6)
      Ccount2 = Ccount2 - Ccount(6+10)
      Ccount3 = Ccount3 - Ccount(6+20)
      Ccount4 = Ccount4 - Ccount(6+30)
      Dcount1 = Dcount1 - Dcount(6)
      Dcount2 = Dcount2 - Dcount(6+10)
      Dcount3 = Dcount3 - Dcount(6+20)
      Dcount4 = Dcount4 - Dcount(6+30)
      
      if (Dcount1.eq.0d0)  Dcount1=1d0
      if (Ccount1.eq.0d0)  Ccount1=1d0

      write(out,100)  
 100  format (/' Numbers for calls of different branches in C and D',
     &     ' reduction'/)
      
      write(out,300) (Ccount(i),dble(Ccount(i))/Ccount1*1d2,i=1,9),
     &         Ccount(0),dble(Ccount(0))/Ccount1*1d2
      write(out,400) (Dcount(i),dble(Dcount(i))/Dcount1*1d2,i=1,9),
     &         Dcount(0),dble(Dcount(0))/Dcount1*1d2

 300  format(' #calls Cpv = ',i20,' or ',F10.5,' %'/
     &       ' #calls Cg  = ',i20,' or ',F10.5,' %'/
     &       ' #calls Cgy = ',i20,' or ',F10.5,' %'/
     &       ' #calls Cy  = ',i20,' or ',F10.5,' %'/
     &       ' #calls Cgp = ',i20,' or ',F10.5,' %'/
c     &       ' #calls     = ',i20,' or ',F10.5,' %'/
     &       ' #calls     = ',i20,' or ',F10.5,' %'/
     &       ' #calls     = ',i20,' or ',F10.5,' %'/
     &       ' #calls Cpv = ',i20,' or ',F10.5,' %'/
     &       ' #calls Cy  = ',i20,' or ',F10.5,' %'/
     &       ' #calls C   = ',i20,' or ',F10.5,' %'/)
 400  format(' #calls Dpv = ',i20,' or ',F10.5,' %'/
     &       ' #calls Dg  = ',i20,' or ',F10.5,' %'/
     &       ' #calls Dgy = ',i20,' or ',F10.5,' %'/
     &       ' #calls Dy  = ',i20,' or ',F10.5,' %'/
     &       ' #calls Dgp = ',i20,' or ',F10.5,' %'/
     &       ' #calls     = ',i20,' or ',F10.5,' %'/
     &       ' #calls     = ',i20,' or ',F10.5,' %'/
c     &       ' #calls     = ',i20,' or ',F10.5,' %'/
c     &       ' #calls     = ',i20,' or ',F10.5,' %'/
     &       ' #calls Dpv = ',i20,' or ',F10.5,' %'/
     &       ' #calls Dy  = ',i20,' or ',F10.5,' %'/
     &       ' #calls D   = ',i20,' or ',F10.5,' %'/)

 101  format (/' Numbers for calls of different branches in C and D',
     &     ' reduction'/ ' accuracy worse than',G12.4)

      write(out,101) accbad1
      write(out,300) (Ccount(i),dble(Ccount(i))/Ccount1*1d2,i=11,19),
     &         Ccount2,dble(Ccount2)/Ccount1*1d2
      write(out,400) (Dcount(i),dble(Dcount(i))/Dcount1*1d2,i=11,19),
     &         Dcount2,dble(Dcount2)/Dcount1*1d2

      write(out,101) accbad2
      write(out,300) (Ccount(i),dble(Ccount(i))/Ccount1*1d2,i=21,29),
     &         Ccount3,dble(Ccount3)/Ccount1*1d2
      write(out,400) (Dcount(i),dble(Dcount(i))/Dcount1*1d2,i=21,29),
     &         Dcount3,dble(Dcount3)/Dcount1*1d2

      write(out,101) accbad3
      write(out,300) (Ccount(i),dble(Ccount(i))/Ccount1*1d2,i=31,39),
     &         Ccount4,dble(Ccount4)/Ccount1*1d2
      write(out,400) (Dcount(i),dble(Dcount(i))/Dcount1*1d2,i=31,39),
     &         Dcount4,dble(Dcount4)/Dcount1*1d2
      end 
