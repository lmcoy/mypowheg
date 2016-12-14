#include "global_coli.h"

#define QCONT

************************************************************************
*
*  BEWARE OF REAL ARGUMENTS OF  CDLN AND CSPCOE, ...
*
*  dependence on msshift not yet tested at all
*
************************************************************************
*                                                                      *
*     Scalar 3- point functions                                        *
*                                                                      *
************************************************************************
*                                                                      *
*     last changed  21.08.08  Ansgar Denner                            *
*                                                                      *
************************************************************************
* Subroutines:                                                         *
* Functions:                                                           *
* C0_coli                                                              *
* C0reg_coli,C0m0_coli,C0m0_colia,C02m0_coli,C03m0_coli                *
* C0ms0ir1_coli, C0ms1ir0_coli, C0ms1ir1_coli, C0ms2ir1_coli           *
************************************************************************

************************************************************************
      function C0_coli(p12,p23,p13,m12,m22,m32)
************************************************************************
*     scalar 3-point function                                          *
*                                                                      *
*                      p23                                             *
*                       |                                              *
*                       |                                              *
*                      / \                                             *
*                     /   \                                            *
*              m22   /     \   m32                                     *
*                   /       \                                          *
*                  /         \                                         *
*       p12  ---------------------  p13                                *
*                      m12                                             *
*                                                                      *
************************************************************************
*     20.08.08 Ansgar Denner        last changed  05.08.08             *
************************************************************************
      implicit   none
      complex*16 C0_coli
      complex*16 p12,p23,p13
      complex*16 m12,m22,m32
      complex*16 ps12,ps23,ps13
      complex*16 ms12,ms22,ms32
      complex*16 p2(3,3),m2(3)

      complex*16 elimminf2_coli
      complex*16 C0reg_coli 
      complex*16 C0m0_coli,C02m0_coli,C03m0_coli
      complex*16 C0ms0ir1_coli
      complex*16 C0ms1ir0_coli,C0ms1ir1_coli
      complex*16 C0ms2ir1_coli

      integer i,j,k
      integer nsm,nsoft,ncoll
      logical smallm2(3),smallp2(3,3),soft(3,3,3),coll(3,3),onsh(3,3)

      data coll /9*.false./, onsh /9*.false./,
     &    soft /27*.false./ ,smallm2 /3*.false./, smallp2 /9*.false./

#include "params_coli.h"
#include "checkparams_coli.h"

 100  format(((a)))
 101  format(a22,g24.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))

c      write(*,*) ' C0_coli in',p12,p23,p13,m12,m22,m32

      p2(1,2) = p12
      p2(1,3) = p13
      p2(2,3) = p23
      m2(1)   = m12
      m2(2)   = m22
      m2(3)   = m32

c determine infinitesimal parameters
      ms12 = elimminf2_coli(m12)
      ms22 = elimminf2_coli(m22)
      ms32 = elimminf2_coli(m32)
      ps12 = elimminf2_coli(p12)
      ps23 = elimminf2_coli(p23)
      ps13 = elimminf2_coli(p13)

c      write(*,*) ' C0_coli small',ps12,ps23,ps13,ms12,ms22,ms32

      smallm2(1)=ms12.eq.cd0
      smallm2(2)=ms22.eq.cd0
      smallm2(3)=ms32.eq.cd0
      smallp2(1,2)=ps12.eq.cd0
      smallp2(1,3)=ps13.eq.cd0
      smallp2(2,3)=ps23.eq.cd0

      nsm=0
      do i=1,3
        if(smallm2(i)) nsm=nsm+1
      enddo

c determine on-shell momenta
      onsh(1,2)=ps12.eq.ms22
      onsh(1,3)=ps13.eq.ms32
      onsh(2,1)=ps12.eq.ms12
      onsh(2,3)=ps23.eq.ms32
      onsh(3,1)=ps13.eq.ms12
      onsh(3,2)=ps23.eq.ms22

c determine collinear singularities
      ncoll=0
      do i=1,2
      do j=2,3
        coll(i,j)=smallm2(i).and.smallm2(j).and.smallp2(i,j)
        if(coll(i,j).and.(.not.
     &      ( (p2(i,j).eq.cd0.and.m2(i).eq.m2(j)).or.
     &        (m2(i).eq.cd0.and.onsh(i,j)).or.  
     &        (m2(j).eq.cd0.and.onsh(j,i)) ))) then
          write(*,*) ' C0_coli: structure of collinear singularity ',
     &        '    not supported'
          write(*,111) ' C0_coli: p12 = ',p12
          write(*,111) ' C0_coli: p23 = ',p23
          write(*,111) ' C0_coli: p13 = ',p13
          write(*,111) ' C0_coli: m12 = ',m12
          write(*,111) ' C0_coli: m22 = ',m22
          write(*,111) ' C0_coli: m32 = ',m32
          write(*,*)   ' C0_coli: i,j = ',i,j
          write(*,*)   ' C0_coli: t1  = ',p2(i,j).eq.cd0,m2(i).eq.m2(j)
          write(*,*)   ' C0_coli: t2  = ',m2(i).eq.cd0,onsh(i,j)
          write(*,*)   ' C0_coli: t3  = ',m2(j).eq.cd0,onsh(j,i)
          stop
        endif
        coll(j,i)=coll(i,j)
        if(coll(i,j)) ncoll=ncoll+1
      enddo
      enddo
c determine soft singularities
      nsoft=0 
      do i=1,3
        do j=1,2
          if(i.ne.j)then
            do k=j,3
              if(k.ne.i.and.k.ne.j)then
                soft(i,j,k)=smallm2(i).and.onsh(i,j).and.onsh(i,k)
                soft(i,k,j)=soft(i,j,k)
                if(soft(i,j,k)) nsoft=nsoft+1
              endif
            enddo
          endif
        enddo
      enddo

c      write(*,*) 'C0_coli: ncoll,nsoft,nsm= ',ncoll,nsoft,nsm

c regular cases
      if(ncoll.eq.0.and.nsoft.eq.0)then
c     regular case
        if(nsm.eq.0)then
            C0_coli = C0reg_coli(ps12,ps23,ps13,ms12,ms22,ms32)
        elseif (nsm.eq.1)then
          if (ms12.eq.cd0) then
            C0_coli = C0m0_coli(ps12,ps23,ps13,cd0,ms22,ms32)
          else if (ms22.eq.cd0) then
            C0_coli = C0m0_coli(ps23,ps13,ps12,cd0,ms32,ms12)
          else if (ms32.eq.cd0) then
            C0_coli = C0m0_coli(ps13,ps12,ps23,cd0,ms12,ms22)
#ifdef CHECK
          else
            write(*,*) ' C0_coli: inconsistency 0ir 0ms 1m0'
#endif
          endif
        elseif(nsm.eq.2)then
          if (ms22.ne.cd0) then
            C0_coli = C02m0_coli(ps12,ps23,ps13,cd0,ms22,cd0)
          else if (ms32.ne.cd0) then
            C0_coli = C02m0_coli(ps23,ps13,ps12,cd0,ms32,cd0)
          else if (ms12.ne.cd0) then
            C0_coli = C02m0_coli(ps13,ps12,ps23,cd0,ms12,cd0)
#ifdef CHECK
          else
            write(*,*) ' C0_coli: inconsistency 0ir 0ms 2m0'
#endif
          endif
        elseif (nsm.eq.3)then
          C0_coli = C03m0_coli(ps13,ps12,ps23,cd0,cd0,cd0)
#ifdef CHECK
        else
          write(*,*) ' C0_coli: inconsistency 0ir 0ms'
#endif
        endif

c IR singular cases without collinear singularities
      elseif(ncoll.eq.0.and.nsoft.ge.1)then
        if(soft(1,2,3))then
          C0_coli = C0ms0ir1_coli(ps12,ps23,ps13,m12,ms22,ms32)
        elseif(soft(2,3,1))then
          C0_coli = C0ms0ir1_coli(ps23,ps13,p12,m22,ms32,ms12)
        elseif(soft(3,1,2))then
          C0_coli = C0ms0ir1_coli(ps13,ps12,ps23,m32,ms12,ms22)
#ifdef CHECK
        else
          write(*,*) ' C0_coli: inconsistency 1ir 0ms'
#endif
        endif

c single collinear case
      elseif(ncoll.eq.1)then
        if(coll(1,2))then
          if(soft(1,2,3))then
            C0_coli = C0ms1ir1_coli(p12,ps23,ps13,m12,m22,ms32)
          elseif(soft(2,3,1))then
            C0_coli = C0ms1ir1_coli(p12,ps13,ps23,m22,m12,ms32)
          else
            C0_coli = C0ms1ir0_coli(p12,ps23,ps13,m12,m22,ms32)
          endif
        elseif(coll(2,3))then
          if(soft(2,3,1))then
            C0_coli = C0ms1ir1_coli(p23,ps13,ps12,m22,m32,ms12)
          elseif(soft(3,1,2))then
            C0_coli = C0ms1ir1_coli(p23,ps12,ps13,m32,m22,ms12)
          else
            C0_coli = C0ms1ir0_coli(p23,ps13,ps12,m22,m32,ms12)
          endif
        elseif(coll(1,3))then
          if(soft(3,1,2))then
            C0_coli = C0ms1ir1_coli(p13,ps12,ps23,m32,m12,ms22)
          elseif(soft(1,2,3))then
            C0_coli = C0ms1ir1_coli(p13,ps23,ps12,m12,m32,ms22)
          else
            C0_coli = C0ms1ir0_coli(p13,ps12,ps23,m32,m12,ms22)
          endif
#ifdef CHECK
        else
          write(*,*) ' C0_coli: inconsistency 1ms'
#endif
        endif

c double collinear case
      elseif(ncoll.eq.2)then
        if(.not.coll(2,3))then
            C0_coli = C0ms2ir1_coli(p12,ps23,p13,m12,m22,m32)
        elseif(.not.coll(1,3))then
            C0_coli = C0ms2ir1_coli(p23,ps13,p12,m22,m32,m12)
        elseif(.not.coll(1,2))then
            C0_coli = C0ms2ir1_coli(p13,ps12,p23,m32,m12,m22)
        else
          write(*,*) ' C0_coli: inconsistency 2ms'
        endif

#ifdef CHECK
      else
        write(*,*) ' C0_coli: case with more than 2 collinear ',
     &      'singularities not implemented'
        stop
#endif
      endif

c      write(*,*) ' C0_coli out ', C0_coli     

      end

************************************************************************
      function C0reg_coli(p12,p23,p13,m12,m22,m32)
************************************************************************
*  scalar 3-point function, regular case                               *
*                                                                      *
*                      p23                                             *
*                       |                                              *
*                       |                                              *
*                      / \                                             *
*                     /   \                                            *
*              m22   /     \   m32                                     *
*                   /       \                                          *
*                  /         \                                         *
*       p12  ---------------------  p13                                *
*                      m12                                             *
*                                                                      *
*     y(i,1) = y0(i), y(i,2) = y0(i)-1                                 *
*     y(i,3) = y(i-) = y0(i) - x(i-), y(i,4) = y(i+) = y0(i) - x(i+)   *
*     x0(i,1) = -x(i-), x0(i,2) = -x(i+)                               *
*     x1(i,1) = 1-x(i-), x1(i,2) = 1-x(i+)                             *
*----------------------------------------------------------------------*
*  10.05.04 Ansgar Denner      last changed 28.10.08                   *
************************************************************************
      implicit   none
      complex*16 p12,p23,p13
      complex*16 m12,m22,m32
      complex*16 mi2(5),mi(5),swap
      complex*16 q2(5)
      complex*16 a,b(3),y(3,4),yh,yh0,yh1,yh2,yh3,yh4,yh5
      complex*16 x0(3,2),x1(3,2),deleta
      real*8     iy(3,4),ix0(3,2),iy1y3,iy2y3,iy1y4,iy2y4
      complex*16 cspenc_coli,cspcon_coli,eta2s_coli,eta2_coli,cln_coli
      complex*16 C0reg_coli,C0reg_check
c      complex*16 C0m0_coli

      integer i,j,k
      
#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps13
      complex*16 ms12,ms22,ms32
      complex*16 eta
      logical flag(0:7)
      data    flag /8*.true./
      save    flag  
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) ' C0reg_coli in',p12,p23,p13,m12,m22,m32

#ifdef CHECK
 100  format(((a)))
 101  format(a22,g24.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        ms12 = elimminf2_coli(m12)
        ms22 = elimminf2_coli(m22)
        ms32 = elimminf2_coli(m32)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps13 = elimminf2_coli(p13)

c>        if(m12.eq.cd0.and.p12.eq.m22.and.p13.eq.m32
c>     &  .or.m22.eq.cd0.and.p12.eq.m12.and.p23.eq.m32
c>     &  .or.m32.eq.cd0.and.p13.eq.m12.and.p23.eq.m22
        if(m12.eq.cd0.or.m22.eq.cd0.or.m32.eq.cd0  
     &      .or.ps12.eq.cd0.and.ms12.eq.cd0.and.ms22.eq.cd0
     &      .or.ps13.eq.cd0.and.ms12.eq.cd0.and.ms32.eq.cd0
     &      .or.ps23.eq.cd0.and.ms22.eq.cd0.and.ms32.eq.cd0
     &      .or.ms12.eq.cd0.and.ps12.eq.ms22.and.ps13.eq.ms32
     &      .or.ms22.eq.cd0.and.ps12.eq.ms12.and.ps23.eq.ms32
     &      .or.ms32.eq.cd0.and.ps13.eq.ms12.and.ps23.eq.ms22
     &      .or.dimag(p12).ne.0d0.or.dimag(p13).ne.0d0
     &      .or.dimag(p23).ne.0d0
     &   ) then
          write(*,100) ' C0reg_coli called improperly:'
          write(*,100) '     or case not implemented    '
          write(*,111) ' C0reg_coli: p12 = ',p12
          write(*,111) ' C0reg_coli: p23 = ',p23
          write(*,111) ' C0reg_coli: p13 = ',p13
          write(*,111) ' C0reg_coli: m12 = ',m12
          write(*,111) ' C0reg_coli: m22 = ',m22
          write(*,111) ' C0reg_coli: m32 = ',m32
          stop
        end if
      end if
#endif

c C0reg_coli becomes unstable for zero masses and finite width
c in the would-be IR limit

c      if (m12.eq.cd0) then
c        C0reg_coli = C0m0_coli(p12,p23,p13,m12,m22,m32)
c      else if (m22.eq.cd0) then
c        C0reg_coli = C0m0_coli(p23,p13,p12,m22,m32,m12)
c      else if (m32.eq.cd0) then
c        C0reg_coli = C0m0_coli(p13,p12,p23,m32,m12,m22)
c      else

      C0reg_coli=cd0
      C0reg_check=cd0
      
      mi2(1)=m12
      mi2(2)=m22
      mi2(3)=m32
      mi2(4)=mi2(1)
      mi2(5)=mi2(2)
      mi(1)=cdsqrt(m12)
      mi(2)=cdsqrt(m22)
      mi(3)=cdsqrt(m32)
      mi(4)=mi(1)
      mi(5)=mi(2)
      q2(1)=p12
      q2(2)=p23
      q2(3)=p13
      q2(4)=q2(1)
      q2(5)=q2(2)
      
c      write(*,*) ' C0reg_coli '
c      write(*,*) ' C0reg_coli q ',p12,p23,p13
c      write(*,*) ' C0reg_coli m2',m12,m22,m32
      if (p12*p23*p13.ne.cd0) then

#ifdef ALLCHECK
c ppttbb
        if(flag(1))then
          write(*,100) ' C0reg_coli: p12,p13,p23=/=0 '
          write(*,100) '    case not yet tested in physical process'
          write(*,111) ' C0reg_coli: p12 = ',p12
          write(*,111) ' C0reg_coli: p23 = ',p23
          write(*,111) ' C0reg_coli: p13 = ',p13
          write(*,111) ' C0reg_coli: m12 = ',m12
          write(*,111) ' C0reg_coli: m22 = ',m22
          write(*,111) ' C0reg_coli: m32 = ',m32
          flag(1)=.false.
        endif
#endif

c ---> program not correct for imaginary alpha!!!
c ---> olieta(1/y4,1/y3) not correct for 1/y3/y4 real!
c        if(q2(1)*q2(1)+q2(2)*q2(2)+q2(3)*q2(3)
c     *       -2d0*(q2(1)*q2(2)+q2(2)*q2(3)+q2(3)*q2(1)).lt.0d0) then
c          write(*,*) 'cC0: imaginary alpha not implemented '
c          write(*,*) 'cC0: unphysical configuration'
c          write(99,*) 'cC0: imaginary alpha not implemented '
c          write(99,*) 'cC0: unphysical configuration'
c          stop
c        end if
        a = cdsqrt(p12*p12 + p23*p23 + p13*p13
     &       - 2d0*(p12*p23 + p23*p13 + p12*p13))
        
        do 10 i=1,3
          j = i+1
          k = i+2
          yh = q2(i)* (q2(i)-q2(j)-q2(k)
     *         +2d0*mi2(k)-mi2(i)-mi2(j))
     *         -(mi2(i)-mi2(j))*(q2(j)-q2(k))
          y(i,1)=(yh+a*( q2(i)-mi2(i)+mi2(j)))/(2d0*a*q2(i))
          y(i,2)=(yh+a*(-q2(i)-mi2(i)+mi2(j)))/(2d0*a*q2(i))
c          b(i)=cdsqrt(q2(i)*q2(i)+mi2(i)*mi2(i)+mi2(j)*mi2(j)
c     *         -2d0*(q2(i)*mi2(i)+mi2(i)*mi2(j)+mi2(j)*q2(i)))
          if (mi2(i)*mi2(j).ne.cd0) then
            b(i)=cdsqrt((q2(i)-(mi(i)+mi(j))**2)
     &           *(q2(i)-(mi(i)-mi(j))**2))
          else if (mi2(i).ne.cd0) then
            b(i) = q2(i) - mi2(i)
          else 
            b(i) = q2(i) - mi2(j)
          end if
          y(i,3)=(yh+a*b(i))/(2d0*a*q2(i))
          y(i,4)=(yh-a*b(i))/(2d0*a*q2(i))

c>          write(*,*) ' C0reg_coli q ',i,q2(i),q2(j),q2(k)
c>          write(*,*) ' C0reg_coli m ',i,mi2(i),mi2(j),mi2(k)
c>          write(*,*) ' C0reg_coli b ',i,yh,a,b(i),a*b(i)
c>          write(*,*) ' C0reg_coli y ',i,y(i,1),y(i,2),y(i,3),y(i,4)

          yh0 = q2(i)*q2(j)*q2(k) + q2(i)*(q2(i)-q2(j)-q2(k))*mi2(k)
     *         + q2(i)*(mi2(k)-mi2(i))*(mi2(k)-mi2(j))
     *         - (q2(j)-q2(k))*mi2(k)*(mi2(i)-mi2(j))
     *         + (mi2(i)-mi2(j))*(q2(j)*mi2(i)-q2(k)*mi2(j))
          yh3 = yh0 + q2(j)*(q2(j)-q2(k)-q2(i))*mi2(i)
     *         + q2(k)*(q2(k)-q2(i)-q2(j))*mi2(j)
          if(cdabs(y(i,3)).lt.cdabs(y(i,4))) then
            y(i,3) = yh3/(a*a*q2(i)*y(i,4))
          else
            y(i,4) = yh3/(a*a*q2(i)*y(i,3))
          end if
          if(dreal(b(i)).ne.0d0) then
            iy(i,3) = dreal(b(i))
            iy(i,4) = -iy(i,3)
          else
            iy(i,3) =  1d0
            iy(i,4) = -1d0
          end if
          yh1=(yh-a*( q2(i)-mi2(i)+mi2(j)))/(2d0*a*q2(i))
          if(cdabs(y(i,1)).lt.cdabs(yh1)) then
            yh4 = yh0 + q2(j)*(q2(j)-q2(k)-q2(i))*mi2(i)
     *        + (q2(k)*(q2(i)+q2(j))-(q2(i)-q2(j))*(q2(i)-q2(j)))*mi2(j)
            y(i,1) = yh4/(a*a*q2(i)*yh1)
          end if
          yh2=(yh-a*(-q2(i)-mi2(i)+mi2(j)))/(2d0*a*q2(i))
          if(cdabs(y(i,2)).lt.cdabs(yh2)) then
            yh5 = yh0 + q2(k)*(q2(k)-q2(i)-q2(j))*mi2(j)
     *        + (q2(j)*(q2(k)+q2(i))-(q2(k)-q2(i))*(q2(k)-q2(i)))*mi2(i)
            y(i,2) = yh5/(a*a*q2(i)*yh2)
          end if

          x0(i,1)=-(q2(i)-mi2(i)+mi2(j)-b(i))/(2d0*q2(i))
          x0(i,2)=-(q2(i)-mi2(i)+mi2(j)+b(i))/(2d0*q2(i))
          x1(i,1)=-(-q2(i)-mi2(i)+mi2(j)-b(i))/(2d0*q2(i))
          x1(i,2)=-(-q2(i)-mi2(i)+mi2(j)+b(i))/(2d0*q2(i))
c>        write(*,*) ' C0reg_coli test',y(i,1),y(i,2)+1,
c>     &        y(i,3),y(i,1)+x0(i,1),
c>     &        y(i,4),y(i,1)+x0(i,2),
c>     &        x1(i,1),1+x0(i,1),
c>     &        x1(i,2),1+x0(i,2)

          if (abs(x0(i,1)).lt.abs(x0(i,2))) then
            x0(i,1) =  mi2(j)/(q2(i)* x0(i,2))
          else
            x0(i,2) =  mi2(j)/(q2(i)* x0(i,1))
          end if
          if (abs(x1(i,1)).lt.abs(x1(i,2))) then
            x1(i,1) =  mi2(i)/(q2(i)* x1(i,2))
          else
            x1(i,2) =  mi2(i)/(q2(i)* x1(i,1))
          end if 
c-> 15.01.07:  iy2y3,iy2y4,iy1y3,iy1y4 introduced 
c (07.09.08 wrong signs?, irrelevant?)
          if(y(i,2).ne.0d0)then
            iy2y3=iy(i,3)/y(i,2)
          else
            iy2y3=0d0
          endif
          if(y(i,2).ne.0d0)then
            iy2y4=iy(i,4)/y(i,2)
          else
            iy2y3=0d0
          endif
          if(y(i,1).ne.0d0)then
            iy1y3=iy(i,3)/y(i,1)
          else
            iy1y3=0d0
          endif
          if(y(i,1).ne.0d0)then
            iy1y4=iy(i,4)/y(i,1)
          else
            iy1y4=0d0
          endif
          C0reg_coli = C0reg_coli
     &        +cspcon_coli(x1(i,1),1d0/y(i,3),x1(i,1)/y(i,3),
     &                          y(i,2)/y(i,3),iy(i,3),-iy(i,3))
     &         -cspcon_coli(x0(i,1),1d0/y(i,3),x0(i,1)/y(i,3),
     &                          y(i,1)/y(i,3),iy(i,3),-iy(i,3))
     &         +cspcon_coli(x1(i,2),1d0/y(i,4),x1(i,2)/y(i,4),
     &                          y(i,2)/y(i,4),iy(i,4),-iy(i,4))
     &         -cspcon_coli(x0(i,2),1d0/y(i,4),x0(i,2)/y(i,4),
     &                          y(i,1)/y(i,4),iy(i,4),-iy(i,4))
     &        -(eta2_coli(x0(i,1),x0(i,2),mi2(j)/q2(i)
     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i)))
     &         - eta2_coli(y(i,3),y(i,4), yh3/(a*a*q2(i))
     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i))))
     &         * cln_coli(y(i,2)/y(i,1),1d0)

c          write(*,*) 'p1',x0(i,1),x0(i,2),x0(i,1)*x0(i,2),mi2(j)/q2(i)
c          write(*,*) 'p2',y(i,3),y(i,4), y(i,3)*y(i,4),yh3/(a*a*q2(i))

c>          write(*,*) ' C0reg_coli x ',
c>     &       x1(i,1),x0(i,1),x1(i,2),x0(i,2)
c>          write(*,*) ' C0reg_coli y ',
c>     &       y(i,3),y(i,4)
c>          write(*,*) ' C0reg_coli x/y ',
c>     &       x1(i,1)/y(i,3),x0(i,1)/y(i,3),x1(i,2)/y(i,4),x0(i,2)/y(i,4)
c>          write(*,*) ' C0reg_coli y/y ',
c>     &       y(i,2)/y(i,3),y(i,1)/y(i,3),y(i,2)/y(i,4),y(i,1)/y(i,4)
c>          write(*,*) ' C0reg_coli s ',
c>     &         +cspcon_coli(x1(i,1),1d0/y(i,3),x1(i,1)/y(i,3),y(i,2)/y(i,3)
c>     &                          ,iy(i,3),-iy(i,3),-iy2y3),
c>     &         -cspcon_coli(x0(i,1),1d0/y(i,3),x0(i,1)/y(i,3),y(i,1)/y(i,3)
c>     &                          ,iy(i,3),-iy(i,3),-iy1y4),
c>     &         +cspcon_coli(x1(i,2),1d0/y(i,4),x1(i,2)/y(i,4),y(i,2)/y(i,4)
c>     &                          ,iy(i,4),-iy(i,4),-iy2y3),
c>     &         -cspcon_coli(x0(i,2),1d0/y(i,4),x0(i,2)/y(i,4),y(i,1)/y(i,4)
c>     &                          ,iy(i,4),-iy(i,4),-iy1y4),
c>     &        -(eta2_coli(x0(i,1),x0(i,2),mi2(j)/q2(i)
c>     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i)))
c>     &         - eta2_coli(y(i,3),y(i,4), yh3/(a*a*q2(i))
c>     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i))))
c>     &         * cln_coli(y(i,2)/y(i,1),1d0)
c>
c>          write(*,*) ' C0reg_coli eta ',
c>     &         eta2_coli(x0(i,1),x0(i,2),mi2(j)/q2(i)
c>     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i))),
c>     &         - eta2_coli(y(i,3),y(i,4), yh3/(a*a*q2(i))
c>     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i)))
c>          write(*,*) ' C0reg_coli etaargs5 ',
c>     &        x0(i,1),x0(i,2),mi2(j)/q2(i)
c>     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i))
c>          write(*,*) ' C0reg_coli etaargs6 ',
c>     &         y(i,3),y(i,4), yh3/(a*a*q2(i))
c>     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i))
    

          if(dimag(yh3/(a*a*q2(i))).lt.0d0.and.dreal(q2(i)).lt.0d0)
     &         C0reg_coli=C0reg_coli+2d0*pi*dcmplx(0d0,1d0)
     &         *cln_coli(y(i,2)/y(i,1),1d0)

c>          if(dimag(yh3/(a*a*q2(i))).lt.0d0.and.dreal(q2(i)).lt.0d0)
c>     &         write(*,*) 2d0*pi*dcmplx(0d0,1d0)
c>     &         *cln_coli(y(i,2)/y(i,1),1d0)     
c         write(*,*) ' C0reg_coli ',i,C0reg_coli

 10     continue
        C0reg_coli=C0reg_coli/a

c        write(*,*) ' C0reg_coli = ',C0reg_coli


#ifdef CHECK
        if(conscheck)then

          a = cdsqrt(p12*p12 + p23*p23 + p13*p13
     &        - 2d0*(p12*p23 + p23*p13 + p12*p13))
        
          do 11 i=1,3
            j = i+1
            k = i+2
            yh = q2(i)* (q2(i)-q2(j)-q2(k)
     *          +2d0*mi2(k)-mi2(i)-mi2(j))
     *          -(mi2(i)-mi2(j))*(q2(j)-q2(k))
            y(i,1)=(yh+a*( q2(i)-mi2(i)+mi2(j)))/(2d0*a*q2(i))
            y(i,2)=(yh+a*(-q2(i)-mi2(i)+mi2(j)))/(2d0*a*q2(i))
c     b(i)=cdsqrt(q2(i)*q2(i)+mi2(i)*mi2(i)+mi2(j)*mi2(j)
c     *         -2d0*(q2(i)*mi2(i)+mi2(i)*mi2(j)+mi2(j)*q2(i)))
            if (mi2(i)*mi2(j).ne.cd0) then
              b(i)=cdsqrt((q2(i)-(mi(i)+mi(j))**2)
     &            *(q2(i)-(mi(i)-mi(j))**2))
            else if (mi2(i).ne.cd0) then
              b(i) = q2(i) - mi2(i)
            else 
              b(i) = q2(i) - mi2(j)
            end if
            y(i,3)=(yh+a*b(i))/(2d0*a*q2(i))
            y(i,4)=(yh-a*b(i))/(2d0*a*q2(i))
            
c>          write(*,*) ' C0reg_check q ',i,q2(i),q2(j),q2(k)
c>          write(*,*) ' C0reg_check m ',i,mi2(i),mi2(j),mi2(k)
c>          write(*,*) ' C0reg_check b ',i,yh,a,b(i),a*b(i)
c>          write(*,*) ' C0reg_check y ',i,y(i,1),y(i,2),y(i,3),y(i,4)

            yh0 = q2(i)*q2(j)*q2(k) + q2(i)*(q2(i)-q2(j)-q2(k))*mi2(k)
     *          + q2(i)*(mi2(k)-mi2(i))*(mi2(k)-mi2(j))
     *          - (q2(j)-q2(k))*mi2(k)*(mi2(i)-mi2(j))
     *          + (mi2(i)-mi2(j))*(q2(j)*mi2(i)-q2(k)*mi2(j))
            yh3 = yh0 + q2(j)*(q2(j)-q2(k)-q2(i))*mi2(i)
     *          + q2(k)*(q2(k)-q2(i)-q2(j))*mi2(j)
            if(cdabs(y(i,3)).lt.cdabs(y(i,4))) then
              y(i,3) = yh3/(a*a*q2(i)*y(i,4))
            else
              y(i,4) = yh3/(a*a*q2(i)*y(i,3))
            end if
            if(dreal(b(i)).ne.0d0) then
              iy(i,3) = dreal(b(i))
              iy(i,4) = -iy(i,3)
            else
              iy(i,3) =  1d0
              iy(i,4) = -1d0
            end if
            yh1=(yh-a*( q2(i)-mi2(i)+mi2(j)))/(2d0*a*q2(i))
            if(cdabs(y(i,1)).lt.cdabs(yh1)) then
              yh4 = yh0 + q2(j)*(q2(j)-q2(k)-q2(i))*mi2(i)
     *            + (q2(k)*(q2(i)+q2(j))-(q2(i)-q2(j))*(q2(i)-q2(j)))
     &            *mi2(j)
              y(i,1) = yh4/(a*a*q2(i)*yh1)
            end if
            yh2=(yh-a*(-q2(i)-mi2(i)+mi2(j)))/(2d0*a*q2(i))
            if(cdabs(y(i,2)).lt.cdabs(yh2)) then
              yh5 = yh0 + q2(k)*(q2(k)-q2(i)-q2(j))*mi2(j)
     *            + (q2(j)*(q2(k)+q2(i))-(q2(k)-q2(i))*(q2(k)-q2(i)))
     &            *mi2(i)
              y(i,2) = yh5/(a*a*q2(i)*yh2)
            end if
c changed 29.09.2008
c            iy1y3 = -dreal(y(i,1)/(y(i,3)*y(i,3)*b(i))) 
c            iy1y4 =  dreal(y(i,1)/(y(i,4)*y(i,4)*b(i))) 
c            iy2y3 = -dreal(y(i,2)/(y(i,3)*y(i,3)*b(i))) 
c            iy2y4 =  dreal(y(i,2)/(y(i,4)*y(i,4)*b(i))) 
            iy1y3 = -dreal(y(i,1)*b(i)) 
            if(iy1y3.eq.0d0) iy1y3=1d0
            iy1y4 =  -iy1y3
            iy2y3 = -dreal(y(i,2)*b(i))
            if(iy2y3.eq.0d0) iy2y3=1d0
            iy2y4 =  -iy2y3
            
            C0reg_check = C0reg_check + cspenc_coli(y(i,2)/y(i,3),iy2y3)
     &          + cspenc_coli(y(i,2)/y(i,4),iy2y4)
     &          - cspenc_coli(y(i,1)/y(i,3),iy1y3)
     &          - cspenc_coli(y(i,1)/y(i,4),iy1y4)
            
c>        write(*,*) ' C0reg_coli yh ',i,yh1,yh2
c>        write(*,*) ' C0reg_coli y ',i,y(i,1),y(i,2),y(i,3),y(i,4)
c>        write(*,*) ' C0reg_coli y/y',i,
c>     *              y(i,2)/y(i,3),y(i,2)/y(i,4)
c>     *             ,y(i,1)/y(i,3),y(i,1)/y(i,4)
c>        write(*,*) ' C0reg_coli -iy*y',i,
c>     *              iy2y3,iy2y4,iy1y3,iy1y4
c>        write(*,*) ' C0reg_check so',i,
c>     &         + cspenc_coli(y(i,2)/y(i,3),iy2y3),
c>     &         + cspenc_coli(y(i,2)/y(i,4),iy2y4),
c>     &         - cspenc_coli(y(i,1)/y(i,3),iy1y3),
c>     &         - cspenc_coli(y(i,1)/y(i,4),iy1y4)
c>        write(*,*) ' C0reg_check ',i,C0reg_check

            x0(i,1)=-(q2(i)-mi2(i)+mi2(j)-b(i))/(2d0*q2(i))
            x0(i,2)=-(q2(i)-mi2(i)+mi2(j)+b(i))/(2d0*q2(i))
            x1(i,1)=-(-q2(i)-mi2(i)+mi2(j)-b(i))/(2d0*q2(i))
            x1(i,2)=-(-q2(i)-mi2(i)+mi2(j)+b(i))/(2d0*q2(i))

c>        write(*,*) ' C0reg_coli test',y(i,1),y(i,2)+1,
c>     &        y(i,3),y(i,1)+x0(i,1),
c>     &        y(i,4),y(i,1)+x0(i,2),
c>     &        x1(i,1),1+x0(i,1),
c>     &        x1(i,2),1+x0(i,2)

             if (abs(x0(i,1)).lt.abs(x0(i,2))) then
              x0(i,1) =  mi2(j)/(q2(i)* x0(i,2))
            else
              x0(i,2) =  mi2(j)/(q2(i)* x0(i,1))
            end if
            if (abs(x1(i,1)).lt.abs(x1(i,2))) then
              x1(i,1) =  mi2(i)/(q2(i)* x1(i,2))
            else
              x1(i,2) =  mi2(i)/(q2(i)* x1(i,1))
            end if 
            C0reg_check = C0reg_check 
     &          + eta2s_coli(x1(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy2y3)
     &          *cln_coli(y(i,2)/y(i,3),iy2y3)
     &          - eta2s_coli(x0(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy1y3)
     &          *cln_coli(y(i,1)/y(i,3),iy1y3)
     &          + eta2s_coli(x1(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy2y4)
     &          *cln_coli(y(i,2)/y(i,4),iy2y4)
     &          - eta2s_coli(x0(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy1y4)
     &          *cln_coli(y(i,1)/y(i,4),iy1y4)
     &          -(eta2_coli(x0(i,1),x0(i,2),mi2(j)/q2(i),iy(i,3),iy(i,4)
     &          ,-dreal(1d0/q2(i)))
     &          - eta2_coli(y(i,3),y(i,4),yh3/(a*a*q2(i)),
     &          iy(i,3),iy(i,4),-dreal(1d0/q2(i))))
     &          * cln_coli(y(i,2)/y(i,1),1d0)
c>        write(*,*) ' C0reg_coli x',i,
c>     &         x1(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy2y3,
c>     &         y(i,2)/y(i,3),
c>     &         x0(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy1y3,
c>     &         y(i,1)/y(i,3),
c>     &         x1(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy2y4,
c>     &         y(i,2)/y(i,4),
c>     &         x0(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy1y4,
c>     &         y(i,1)/y(i,4),
c>     &         x0(i,1),x0(i,2),iy(i,3),iy(i,4),-dreal(1d0/q2(i)),
c>     &         y(i,3),y(i,4),iy(i,3),iy(i,4),-dreal(1d0/q2(i))

c        write(*,*) ' C0reg_coli x2',i,
c     &         eta2s_coli(x0(i,1),x0(i,2),iy(i,3),iy(i,4),-1d0/q2(i)),
c     &         - eta2s_coli(y(i,3),y(i,4),iy(i,3),iy(i,4),-1d0/q2(i)), 
c     &       x0(i,1),x0(i,2),y(i,3),y(i,4)
c>            write(*,*) ' C0reg_check se',i,
c>     &         + eta2s_coli(x1(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy2y3)
c>     &         *cln_coli(y(i,2)/y(i,3),iy2y3),
c>     &         - eta2s_coli(x0(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy1y3)
c>     &         *cln_coli(y(i,1)/y(i,3),iy1y3),
c>     &         + eta2s_coli(x1(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy2y4)
c>     &         *cln_coli(y(i,2)/y(i,4),iy2y4),
c>     &         - eta2s_coli(x0(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy1y4)
c>     &         *cln_coli(y(i,1)/y(i,4),iy1y4),
c>     &          -(eta2_coli(x0(i,1),x0(i,2),mi2(j)/q2(i),iy(i,3),iy(i,4)
c>     &          ,-dreal(1d0/q2(i)))
c>     &          - eta2_coli(y(i,3),y(i,4),yh3/(a*a*q2(i)),
c>     &          iy(i,3),iy(i,4),-dreal(1d0/q2(i))))
c>     &          * cln_coli(y(i,2)/y(i,1),1d0)
c>         write(*,*) ' C0reg_check e',i,
c>     &         + eta2s_coli(x1(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy2y3),
c>     &         - eta2s_coli(x0(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy1y3),
c>     &         + eta2s_coli(x1(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy2y4),
c>     &         - eta2s_coli(x0(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy1y4),
c>     &          -eta2s_coli(x0(i,1),x0(i,2),mi2(j)/q2(i),iy(i,3),iy(i,4),
c>     &          dreal(-1d0/q2(i))),
c>     &          + eta2s_coli(y(i,3),y(i,4),iy(i,3),iy(i,4),
c>     &          -dreal(1d0/q2(i)))

c>        write(*,*) ' C0reg_check ea1',i,x1(i,1),1d0/y(i,3),x1(i,1)/y(i,3),
c>     &               iy(i,3),-iy(i,3),-iy2y3
c>        write(*,*) ' C0reg_check ea2',i,x0(i,1),1d0/y(i,3),x0(i,1)/y(i,3),
c>     &               iy(i,3),-iy(i,3),-iy1y3
c>        write(*,*) ' C0reg_check ea3',i,x1(i,2),1d0/y(i,4),x1(i,2)/y(i,4),
c>     &               iy(i,4),-iy(i,4),-iy2y4
c>        write(*,*) ' C0reg_check ea4',i,x0(i,2),1d0/y(i,4),x0(i,2)/y(i,4),
c>     &               iy(i,4),-iy(i,4),-iy1y4
c>         write(*,*) ' C0reg_check ea5',i,x0(i,1),x0(i,2),mi2(j)/q2(i),
c>     &               iy(i,3),iy(i,4),dreal(-1d0/q2(i))
c>         write(*,*) ' C0reg_check ea6',i,y(i,3),y(i,4),yh3/(a*a*q2(i)),
c>     &               iy(i,3),iy(i,4),dreal(-1d0/q2(i))
        
            if(dimag(yh3/(a*a*q2(i))).lt.0d0.and.dreal(q2(i)).lt.0d0)
     &          C0reg_check=C0reg_check+2d0*pi*dcmplx(0d0,1d0)
     &          *cln_coli(y(i,2)/y(i,1),1d0)

c>            if(dimag(yh3/(a*a*q2(i))).lt.0d0.and.dreal(q2(i)).lt.0d0)
c>     &          write(*,*) ' C0reg_check se2',i,2d0*pi*dcmplx(0d0,1d0)
c>     &         *cln_coli(y(i,2)/y(i,1),1d0)

c          write(*,*) ' C0reg_check ',i,C0reg_check

 11       continue
          C0reg_check=C0reg_check/a

          if(cdabs(C0reg_coli/C0reg_check-1d0).gt.testacc) then
            write(testout,100) ' C0reg_coli inconsistent'
            write(testout,111) ' C0reg_coli: p12 = ',p12
            write(testout,111) ' C0reg_coli: p23 = ',p23
            write(testout,111) ' C0reg_coli: p13 = ',p13
            write(testout,111) ' C0reg_coli: m12 = ',m12
            write(testout,111) ' C0reg_coli: m22 = ',m22
            write(testout,111) ' C0reg_coli: m32 = ',m32
            write(testout,111) ' C0reg_coli  = ',C0reg_coli
            write(testout,111) ' C0reg_check = ',C0reg_check
            write(testout,101) ' diff       = ',
     &          cdabs(C0reg_coli/C0reg_check-1d0)
            
            write(*,100) ' C0reg_coli inconsistent'
            write(*,111) ' C0reg_coli: p12 = ',p12
            write(*,111) ' C0reg_coli: p23 = ',p23
            write(*,111) ' C0reg_coli: p13 = ',p13
            write(*,111) ' C0reg_coli: m12 = ',m12
            write(*,111) ' C0reg_coli: m22 = ',m22
            write(*,111) ' C0reg_coli: m32 = ',m32
            write(*,111) ' C0reg_coli  = ',C0reg_coli
            write(*,111) ' C0reg_check = ',C0reg_check
            write(*,111) ' diff        = ',
     &          cdabs(C0reg_coli/C0reg_check-1d0)
c     stop
          endif
        endif
#endif
        
      else if (p23*p12+p13*p23+p12*p13.ne.0d0) then


        if (p23.eq.0d0) then
          mi2(1)=m22
          mi2(2)=m32
          mi2(3)=m12
          mi2(4)=mi2(1)
          mi2(5)=mi2(2)
          swap = mi(1)
          mi(1)=mi(2)
          mi(2)=mi(3)
          mi(3)=swap
          mi(4)=mi(1)
          mi(5)=mi(2)
          q2(1)=p23
          q2(2)=p13
          q2(3)=p12
          q2(4)=q2(1)
          q2(5)=q2(2)
        else if (p13.eq.0d0) then
          mi2(1)=m32
          mi2(2)=m12
          mi2(3)=m22
          mi2(4)=mi2(1)
          mi2(5)=mi2(2)
          swap = mi(1)
          mi(1)=mi(3)
          mi(3)=mi(2)
          mi(2)=swap
          mi(4)=mi(1)
          mi(5)=mi(2)
          q2(1)=p13
          q2(2)=p12
          q2(3)=p23
          q2(4)=q2(1)
          q2(5)=q2(2)
        end if
 
        a=q2(2)-q2(3)

        if(a.ne.0d0)then

#ifdef ALLCHECK
c ppttbb
          if(flag(2))then
            write(*,100) ' C0reg_coli: p12=0, p23=/=p13 or symm. '
            write(*,100) '    case not yet tested in physical process'
            write(*,111) ' C0reg_coli: p12 = ',p12
            write(*,111) ' C0reg_coli: p23 = ',p23
            write(*,111) ' C0reg_coli: p13 = ',p13
            write(*,111) ' C0reg_coli: m12 = ',m12
            write(*,111) ' C0reg_coli: m22 = ',m22
            write(*,111) ' C0reg_coli: m32 = ',m32
            flag(2)=.false.
          endif
#endif
 
          y(2,1) = ( q2(2)-q2(3)+mi2(1)-mi2(2)) / a
          y(2,2) = ( mi2(1)-mi2(2)) / a
          y(3,1) = ( mi2(2)-mi2(1)) / a
          y(3,2) = ( q2(3)-q2(2)+mi2(2)-mi2(1)) / a
 
          do 30 i=2,3
            j=i+1
            k=i+2
c          b(i) = cdsqrt(dcmplx(q2(i)*q2(i)+mi2(i)*mi2(i)+mi2(j)*mi2(j)
c     *         -2d0*(q2(i)*mi2(i)+mi2(i)*mi2(j)+mi2(j)*q2(i))))
            b(i)=cdsqrt((q2(i)-(mi(i)+mi(j))**2)*(q2(i)
     &          -(mi(i)-mi(j))**2))
            yh = q2(i)* (q2(i)-q2(j)-q2(k)
     *          +2d0*mi2(k)-mi2(i)-mi2(j))
     *          -(mi2(i)-mi2(j))*(q2(j)-q2(k))
            y(i,3) = (yh+a*b(i))/(2d0*a*q2(i))
            y(i,4) = (yh-a*b(i))/(2d0*a*q2(i))

c>        write(*,*) ' C0reg_coli q ',i,q2(i),q2(j),q2(k)
c>        write(*,*) ' C0reg_coli m ',i,mi(i),mi(j),mi(k)
c>        write(*,*) ' C0reg_coli m2',i,mi2(i),mi2(j),mi2(k)
c>        write(*,*) ' C0reg_coli b ',i,yh,a,b(i),a*b(i)
c>        write(*,*) ' C0reg_coli y ',i,y(i,1),y(i,2),y(i,3),y(i,4)

            yh0 = q2(i)*q2(j)*q2(k) + q2(i)*(q2(i)-q2(j)-q2(k))*mi2(k)
     *          + q2(i)*(mi2(k)-mi2(i))*(mi2(k)-mi2(j))
     *          - (q2(j)-q2(k))*mi2(k)*(mi2(i)-mi2(j))
     *          + (mi2(i)-mi2(j))*(q2(j)*mi2(i)-q2(k)*mi2(j))
            yh3 = yh0 + q2(j)*(q2(j)-q2(k)-q2(i))*mi2(i)
     *          + q2(k)*(q2(k)-q2(i)-q2(j))*mi2(j)
            if(cdabs(y(i,3)).lt.cdabs(y(i,4))) then
              y(i,3) = yh3/(a*a*q2(i)*y(i,4))
            else
              y(i,4) = yh3/(a*a*q2(i)*y(i,3))
            end if
            if(dreal(b(i)).ne.0d0) then
              iy(i,3) = dreal(b(i))
              iy(i,4) = -iy(i,3)
            else
              iy(i,3) =  1d0
              iy(i,4) = -1d0
            end if
            x0(i,1)=-(q2(i)-mi2(i)+mi2(j)-b(i))/(2d0*q2(i))
            x0(i,2)=-(q2(i)-mi2(i)+mi2(j)+b(i))/(2d0*q2(i))
            x1(i,1)=-(-q2(i)-mi2(i)+mi2(j)-b(i))/(2d0*q2(i))
            x1(i,2)=-(-q2(i)-mi2(i)+mi2(j)+b(i))/(2d0*q2(i))

c>        write(*,*) ' C0reg_coli test',y(i,1),y(i,2)+1,
c>     &        y(i,3),y(i,1)+x0(i,1),
c>     &        y(i,4),y(i,1)+x0(i,2),
c>     &        x1(i,1),1+x0(i,1),
c>     &        x1(i,2),1+x0(i,2)

            if (abs(x0(i,1)).lt.abs(x0(i,2))) then
              x0(i,1) =  mi2(j)/(q2(i)* x0(i,2))
            else
              x0(i,2) =  mi2(j)/(q2(i)* x0(i,1))
            end if
            if (abs(x1(i,1)).lt.abs(x1(i,2))) then
              x1(i,1) =  mi2(i)/(q2(i)* x1(i,2))
            else
              x1(i,2) =  mi2(i)/(q2(i)* x1(i,1))
            end if 

c-> 15.01.07:  iy2y3,iy2y4,iy1y3,iy1y4 introduced
c (07.09.08 wrong signs?, irrelevant?)
            if(y(i,2).ne.0d0)then
              iy2y3=iy(i,3)/y(i,2)
            else
              iy2y3=0d0
            endif
            if(y(i,2).ne.0d0)then
              iy2y4=iy(i,4)/y(i,2)
            else
              iy2y3=0d0
            endif
            if(y(i,1).ne.0d0)then
              iy1y3=iy(i,3)/y(i,1)
            else
              iy1y3=0d0
            endif
            if(y(i,1).ne.0d0)then
              iy1y4=iy(i,4)/y(i,1)
            else
              iy1y4=0d0
            endif

c>          write(*,*) ' C0reg_coli args1',
c>     &         x1(i,1),1d0/y(i,3),x1(i,1)/y(i,3),y(i,2)/y(i,3)
c>     &                          ,iy(i,3),-iy(i,3)
c>          write(*,*) ' C0reg_coli args2',
c>     &         x1(i,2),1d0/y(i,4),x1(i,2)/y(i,4),y(i,2)/y(i,4)
c>     &                          ,iy(i,4),-iy(i,4)
c>          write(*,*) ' C0reg_coli args3',
c>     &        x0(i,1),1d0/y(i,3),x0(i,1)/y(i,3),y(i,1)/y(i,3)
c>     &                          ,iy(i,3),-iy(i,3)
c>          write(*,*) ' C0reg_coli args4',
c>     &         x0(i,2),1d0/y(i,4),x0(i,2)/y(i,4),y(i,1)/y(i,4)
c>     &                          ,iy(i,4),-iy(i,4)
c>
c>          write(*,*) ' C0reg_coli s1',
c>     &         +cspcon_coli(x1(i,1),1d0/y(i,3),x1(i,1)/y(i,3),y(i,2)/y(i,3)
c>     &                          ,iy(i,3),-iy(i,3),-iy2y3)
c>          write(*,*) ' C0reg_coli s2',
c>     &         +cspcon_coli(x1(i,2),1d0/y(i,4),x1(i,2)/y(i,4),y(i,2)/y(i,4)
c>     &                          ,iy(i,4),-iy(i,4),-iy1y3)
c>          write(*,*) ' C0reg_coli s3',
c>     &         -cspcon_coli(x0(i,1),1d0/y(i,3),x0(i,1)/y(i,3),y(i,1)/y(i,3)
c>     &                          ,iy(i,3),-iy(i,3),-iy2y4)
c>          write(*,*) ' C0reg_coli s4',
c>     &         -cspcon_coli(x0(i,2),1d0/y(i,4),x0(i,2)/y(i,4),y(i,1)/y(i,4)
c>     &                          ,iy(i,4),-iy(i,4),-iy1y4)

            C0reg_coli = C0reg_coli
     &         +cspcon_coli(x1(i,1),1d0/y(i,3),x1(i,1)/y(i,3)
     &                          ,y(i,2)/y(i,3),iy(i,3),-iy(i,3))
     &         +cspcon_coli(x1(i,2),1d0/y(i,4),x1(i,2)/y(i,4)
     &                          ,y(i,2)/y(i,4),iy(i,4),-iy(i,4))
     &         -cspcon_coli(x0(i,1),1d0/y(i,3),x0(i,1)/y(i,3)
     &                          ,y(i,1)/y(i,3),iy(i,3),-iy(i,3))
     &         -cspcon_coli(x0(i,2),1d0/y(i,4),x0(i,2)/y(i,4)
     &                          ,y(i,1)/y(i,4),iy(i,4),-iy(i,4))

c--> changed 15.01.07
c     &        -(eta2_coli(x0(i,1),x0(i,2),mi2(j)/q2(i)
c     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i)))
c     &         - eta2_coli(y(i,3),y(i,4), yh3/(a*a*q2(i))
c     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i))))
c     &         * cln_coli(y(i,2)/y(i,1),1d0)

            deleta = eta2_coli(x0(i,1),x0(i,2),mi2(j)/q2(i)
     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i)))
     &         - eta2_coli(y(i,3),y(i,4), yh3/(a*a*q2(i))
     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i)))

c       write(*,*) ' C0reg_coli deleta',deleta

            if(deleta.ne.0d0)then
              C0reg_coli = C0reg_coli
     &            -deleta*cln_coli(y(i,2)/y(i,1),1d0)

c>          write(*,*) ' C0reg_coli s5',
c>     &            -deleta*cln_coli(y(i,2)/y(i,1),1d0)

            endif

c>       write(*,*) ' C0reg_coli eta',i,eta2_coli(x0(i,1),x0(i,2),mi2(j)/q2(i)
c>     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i))),
c>     &         - eta2_coli(y(i,3),y(i,4), yh3/(a*a*q2(i))
c>     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i)))
c>       write(*,*) ' C0reg_coli etaargs1',i,
c>     &     x0(i,1),x0(i,2),mi2(j)/q2(i)
c>     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i))
c>       write(*,*) ' C0reg_coli etaargs2',i,
c>     &         y(i,3),y(i,4), yh3/(a*a*q2(i))
c>     &           ,iy(i,3),iy(i,4),-dreal(1d0/q2(i))

            if(dimag(yh3/(a*a*q2(i))).lt.0d0.and.dreal(q2(i)).lt.0d0)
     &          C0reg_coli=C0reg_coli+2d0*pi*dcmplx(0d0,1d0)
     &          *cln_coli(y(i,2)/y(i,1),1d0)

c>            if(dimag(yh3/(a*a*q2(i))).lt.0d0.and.dreal(q2(i)).lt.0d0)
c>     &          write(*,*) ' C0reg_coli s5',
c>     &          2d0*pi*dcmplx(0d0,1d0)*cln_coli(y(i,2)/y(i,1),1d0)

c            write(*,*) ' C0reg_coli ',i,C0reg_coli

 30       continue
          C0reg_coli=C0reg_coli/a

c          write(*,*) ' C0reg_coli = ',C0reg_coli
         
#ifdef CHECK
          if(conscheck)then

            y(2,1) = 2d0*q2(2)*( q2(2)-q2(3)+mi2(1)-mi2(2))
            y(2,2) = 2d0*q2(2)*( mi2(1)-mi2(2))
            y(3,1) = 2d0*q2(3)*( mi2(2)-mi2(1))
            y(3,2) = 2d0*q2(3)*( q2(3)-q2(2)+mi2(2)-mi2(1))
            
            do 20 i=2,3
              j=i+1
              k=i+2
c          b(i) = cdsqrt(dcmplx(q2(i)*q2(i)+mi2(i)*mi2(i)+mi2(j)*mi2(j)
c     *         -2d0*(q2(i)*mi2(i)+mi2(i)*mi2(j)+mi2(j)*q2(i))))
              b(i)=cdsqrt((q2(i)-(mi(i)+mi(j))**2)*(q2(i)
     &            -(mi(i)-mi(j))**2))
              yh = q2(i)* (q2(i)-q2(j)-q2(k)
     *            +2d0*mi2(k)-mi2(i)-mi2(j))
     *            -(mi2(i)-mi2(j))*(q2(j)-q2(k))
              y(i,3) = yh+a*b(i)
              y(i,4) = yh-a*b(i)

c>              write(*,*) ' C0reg_check q ',i,q2(i),q2(j),q2(k)
c>            write(*,*) ' C0reg_check m ',i,mi(i),mi(j),mi(k)
c>            write(*,*) ' C0reg_check m2',i,mi2(i),mi2(j),mi2(k)
c>            write(*,*) ' C0reg_check b ',i,yh,a,b(i),a*b(i)
c>              write(*,*) ' C0reg_check y ',i,y(i,1),y(i,2),y(i,3),y(i,4)

              yh0 = q2(i)*q2(j)*q2(k) + q2(i)*(q2(i)-q2(j)-q2(k))*mi2(k)
     *            + q2(i)*(mi2(k)-mi2(i))*(mi2(k)-mi2(j))
     *            - (q2(j)-q2(k))*mi2(k)*(mi2(i)-mi2(j))
     *            + (mi2(i)-mi2(j))*(q2(j)*mi2(i)-q2(k)*mi2(j))
              yh3 = yh0 + q2(j)*(q2(j)-q2(k)-q2(i))*mi2(i)
     *            + q2(k)*(q2(k)-q2(i)-q2(j))*mi2(j)
              if(cdabs(y(i,3)).lt.cdabs(y(i,4))) then
                y(i,3) = 4d0*yh3*q2(i)/y(i,4)
              else
                y(i,4) = 4d0*yh3*q2(i)/y(i,3)
              end if
              if(dreal(b(i)).ne.0d0) then
                iy(i,3) = dreal(b(i)*a*q2(i))
                iy(i,4) = -iy(i,3)
              else
                iy(i,3) =  1d0
                iy(i,4) = -1d0
              end if

c       write(*,*)  ' C0reg_coli  6'
c>              write(*,*) ' C0reg_coli y ',i,y(i,1),y(i,2),y(i,3),y(i,4)
c>              write(*,*) ' C0reg_coli iy/y',i,
c>     *            -dreal(y(i,2))*iy(i,3),-dreal(y(i,2))*iy(i,4)
c>     *            ,-dreal(y(i,1))*iy(i,3),-dreal(y(i,1))*iy(i,4)
c>              write(*,*) ' C0reg_coli y/y',i,
c>     *            y(i,2)/y(i,3),y(i,2)/y(i,4)
c>     *            ,y(i,1)/y(i,3),y(i,1)/y(i,4)

              iy2y3=-dreal(y(i,2))*iy(i,3)
              iy2y4=-dreal(y(i,2))*iy(i,4)
              iy1y3=-dreal(y(i,1))*iy(i,3)
              iy1y4=-dreal(y(i,1))*iy(i,4)

c>              write(*,*) ' C0reg_check s',i,
c>     &            + cspenc_coli(y(i,2)/y(i,3),iy2y3),
c>     &            + cspenc_coli(y(i,2)/y(i,4),iy2y4),
c>     &            - cspenc_coli(y(i,1)/y(i,3),iy1y3),
c>     &            - cspenc_coli(y(i,1)/y(i,4),iy1y4)

              C0reg_check = C0reg_check 
     &            + cspenc_coli(y(i,2)/y(i,3),iy2y3)
     &            + cspenc_coli(y(i,2)/y(i,4),iy2y4)
     &            - cspenc_coli(y(i,1)/y(i,3),iy1y3)
     &            - cspenc_coli(y(i,1)/y(i,4),iy1y4)
c
c              write(*,*) ' C0reg_check ',i,C0reg_check

              x0(i,1)=-(q2(i)-mi2(i)+mi2(j)-b(i))*a
              x0(i,2)=-(q2(i)-mi2(i)+mi2(j)+b(i))*a
              x1(i,1)=-(-q2(i)-mi2(i)+mi2(j)-b(i))*a
              x1(i,2)=-(-q2(i)-mi2(i)+mi2(j)+b(i))*a

              if (abs(x0(i,1)).lt.abs(x0(i,2))) then
                x0(i,1) =  4d0*a*a*mi2(j)*q2(i)/x0(i,2)
              else
                x0(i,2) =  4d0*a*a*mi2(j)*q2(i)/x0(i,1)
              end if
              if (abs(x1(i,1)).lt.abs(x1(i,2))) then
                x1(i,1) =   4d0*a*a*mi2(i)*q2(i)/x1(i,2)
              else
                x1(i,2) =   4d0*a*a*mi2(i)*q2(i)/x1(i,1)
              end if 

c>              write(*,*) ' C0reg_check test',y(i,1),y(i,2)+1,
c>     &            y(i,3),y(i,1)+x0(i,1),
c>     &            y(i,4),y(i,1)+x0(i,2),
c>     &            x1(i,1),1+x0(i,1),
c>     &            x1(i,2),1+x0(i,2)

c>              write(*,*) ' C0reg_check ea1',i,x1(i,1),1d0/y(i,3),
c>     &            x1(i,1)/y(i,3),iy(i,3),-iy(i,3),-iy2y3
c>              write(*,*) ' C0reg_check ea2',i,x0(i,1),1d0/y(i,3),
c>     &            x0(i,1)/y(i,3),iy(i,3),-iy(i,3),-iy1y3
c>              write(*,*) ' C0reg_check ea3',i,x1(i,2),1d0/y(i,4),
c>     &            x1(i,2)/y(i,4),iy(i,4),-iy(i,4),-iy2y4
c>              write(*,*) ' C0reg_check ea4',i,x0(i,2),1d0/y(i,4),
c>     &            x0(i,2)/y(i,4),iy(i,4),-iy(i,4),-iy1y4
c>              write(*,*) ' C0reg_check ea5',i,x0(i,1),x0(i,2),
c>     &            mi2(j)/q2(i),iy(i,3),iy(i,4),dreal(-1d0/q2(i))
c>              write(*,*) ' C0reg_check ea6',i,y(i,3),y(i,4),
c>     &            4d0*yh3*q2(i),iy(i,3),iy(i,4),dreal(-1d0/q2(i))
c>
c>              write(*,*) ' C0reg_check e',i,
c>     &         + eta2s_coli(x1(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy2y3),
c>     &         - eta2s_coli(x0(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy1y3),
c>     &         + eta2s_coli(x1(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy2y4),
c>     &         - eta2s_coli(x0(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy1y4),
c>     &         -eta2_coli(x0(i,1),x0(i,2),4d0*a*a*mi2(j)*q2(i),
c>     &          iy(i,3),iy(i,4),dreal(-1d0/q2(i))),
c>     &          + eta2s_coli(y(i,3),y(i,4),iy(i,3),iy(i,4),
c>     &          -dreal(1d0/q2(i)))

              deleta=
     &          +eta2s_coli(x1(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy2y3)
              if(deleta.ne.cd0)then
                C0reg_check = C0reg_check 
     &              +deleta*cln_coli(y(i,2)/y(i,3),iy2y3)

c>                write(*,*) ' C0reg_check se1',i,
c>     &              +eta2s_coli(x1(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),
c>     &              -iy2y3)*cln_coli(y(i,2)/y(i,3),iy2y3)

              endif
              deleta=
     &            - eta2s_coli(x0(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),
     &            -iy1y3)
              if(deleta.ne.cd0)then
                C0reg_check = C0reg_check 
     &              +deleta *cln_coli(y(i,1)/y(i,3),iy1y3)

c>                write(*,*) ' C0reg_check se2',i,
c>     &         - eta2s_coli(x0(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy1y3)
c>     &         *cln_coli(y(i,1)/y(i,3),iy1y3)

              endif
              deleta=
     &            + eta2s_coli(x1(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),
     &            -iy2y4)
              if(deleta.ne.cd0)then
                C0reg_check = C0reg_check 
     &              +deleta*cln_coli(y(i,2)/y(i,4),iy2y4)

c>                write(*,*) ' C0reg_check se3',i,
c>     &         + eta2s_coli(x1(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy2y4)
c>     &         *cln_coli(y(i,2)/y(i,4),iy2y4)

              endif
              deleta=
     &            - eta2s_coli(x0(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),
     &            -iy1y4)
              if(deleta.ne.cd0)then
                C0reg_check = C0reg_check 
     &              +deleta*cln_coli(y(i,1)/y(i,4),iy1y4)

c>                write(*,*) ' C0reg_check se4',i,
c>     &         - eta2s_coli(x0(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy1y4)
c>     &         *cln_coli(y(i,1)/y(i,4),iy1y4)

              endif
              deleta=
     &            -(eta2_coli(x0(i,1),x0(i,2),4d0*a*a*mi2(j)*q2(i),
     &            iy(i,3),iy(i,4),-dreal(1d0/q2(i)))
     &            - eta2_coli(y(i,3),y(i,4),4d0*yh3*q2(i),
     &            iy(i,3),iy(i,4),-dreal(1d0/q2(i))))
              if(deleta.ne.cd0)then
                C0reg_check = C0reg_check 
     &              +deleta* cln_coli(y(i,2)/y(i,1),1d0)

c>                write(*,*) ' C0reg_check se5',i,
c>     &          -(eta2_coli(x0(i,1),x0(i,2),4d0*a*a*mi2(j)*q2(i),
c>     &              iy(i,3),iy(i,4),-dreal(1d0/q2(i)))
c>     &          - eta2_coli(y(i,3),y(i,4),4d0*yh3*q2(i),
c>     &          iy(i,3),iy(i,4),-dreal(1d0/q2(i))))
c>     &          * cln_coli(y(i,2)/y(i,1),1d0)

              endif

c changed 06.09.2008
c>              C0reg_check = C0reg_check 
c>     &            + eta2s_coli(x1(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy2y3)
c>     &            *cln_coli(y(i,2)/y(i,3),iy2y3)
c>     &            - eta2s_coli(x0(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy1y3)
c>     &            *cln_coli(y(i,1)/y(i,3),iy1y3)
c>     &            + eta2s_coli(x1(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy2y4)
c>     &            *cln_coli(y(i,2)/y(i,4),iy2y4)
c>     &            - eta2s_coli(x0(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy1y4)
c>     &            *cln_coli(y(i,1)/y(i,4),iy1y4)
c>     &            -(eta2_coli(x0(i,1),x0(i,2),mi2(j)/q2(i),
c>     &            iy(i,3),iy(i,4),-dreal(1d0/q2(i)))
c>     &            - eta2_coli(y(i,3),y(i,4),yh3/(a*a*q2(i)),
c>     &            iy(i,3),iy(i,4),-dreal(1d0/q2(i))))
c>     &            * cln_coli(y(i,2)/y(i,1),1d0)

c>        write(*,*) ' C0reg_coli x',i,
c>     &         x1(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy2y3,
c>     &         y(i,2)/y(i,3),
c>     &         x0(i,1),1d0/y(i,3),iy(i,3),-iy(i,3),-iy1y3,
c>     &         y(i,1)/y(i,3),
c>     &         x1(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy2y4,
c>     &         y(i,2)/y(i,4),
c>     &         x0(i,2),1d0/y(i,4),iy(i,4),-iy(i,4),-iy1y4,
c>     &         y(i,1)/y(i,4),
c>     &         x0(i,1),x0(i,2),iy(i,3),iy(i,4),-dreal(1d0/q2(i)),
c>     &         y(i,3),y(i,4),iy(i,3),iy(i,4),-dreal(1d0/q2(i))

c        write(*,*) ' C0reg_coli x2',i,
c     &         eta2s_coli(x0(i,1),x0(i,2),iy(i,3),iy(i,4),-1d0/q2(i)),
c     &         - eta2s_coli(y(i,3),y(i,4),iy(i,3),iy(i,4),-1d0/q2(i)), 
c     &       x0(i,1),x0(i,2),y(i,3),y(i,4)

        
            if(dimag(yh3/(a*a*q2(i))).lt.0d0.and.dreal(q2(i)).lt.0d0)
     &          C0reg_check=C0reg_check+2d0*pi*dcmplx(0d0,1d0)
     &          *cln_coli(y(i,2)/y(i,1),1d0)

c>            if(dimag(yh3/(a*a*q2(i))).lt.0d0.and.dreal(q2(i)).lt.0d0)
c>     &          write(*,*) ' C0reg_check se6',i,2d0*pi*dcmplx(0d0,1d0)
c>     &         *cln_coli(y(i,2)/y(i,1),1d0)
c>
c>            write(*,*) ' C0reg_check ',i,C0reg_check

 20         continue
            C0reg_check=C0reg_check/a

c            write(*,*) ' C0reg_check = ',C0reg_check
 


            if(cdabs(C0reg_coli/C0reg_check-1d0).gt.testacc) then
              write(testout,100) ' C0reg_coli inconsistent'
              write(testout,111) ' C0reg_coli: p12 = ',p12
              write(testout,111) ' C0reg_coli: p23 = ',p23
              write(testout,111) ' C0reg_coli: p13 = ',p13
              write(testout,111) ' C0reg_coli: m12 = ',m12
              write(testout,111) ' C0reg_coli: m22 = ',m22
              write(testout,111) ' C0reg_coli: m32 = ',m32
              write(testout,111) ' C0reg_coli  = ',C0reg_coli
              write(testout,111) ' C0reg_check = ',C0reg_check
              write(testout,101) ' diff       = ',
     &            cdabs(C0reg_coli/C0reg_check-1d0)
              
              write(*,100) ' C0reg_coli inconsistent'
              write(*,111) ' C0reg_coli: p12 = ',p12
              write(*,111) ' C0reg_coli: p23 = ',p23
              write(*,111) ' C0reg_coli: p13 = ',p13
              write(*,111) ' C0reg_coli: m12 = ',m12
              write(*,111) ' C0reg_coli: m22 = ',m22
              write(*,111) ' C0reg_coli: m32 = ',m32
              write(*,111) ' C0reg_coli  = ',C0reg_coli
              write(*,111) ' C0reg_check = ',C0reg_check
              write(*,111) ' diff        = ',
     &            cdabs(C0reg_coli/C0reg_check-1d0)
c     stop
            endif
          endif
#endif


        else

#ifdef CHECK
          if(flag(3))then
            write(*,100) ' C0reg_coli: p12=0, p23=p13 or symm. '
            write(*,100) '    case not yet tested at all'
            write(*,111) ' C0reg_coli: p12 = ',p12
            write(*,111) ' C0reg_coli: p23 = ',p23
            write(*,111) ' C0reg_coli: p13 = ',p13
            write(*,111) ' C0reg_coli: m12 = ',m12
            write(*,111) ' C0reg_coli: m22 = ',m22
            write(*,111) ' C0reg_coli: m32 = ',m32
            flag(3)=.false.
          endif
#endif

          do 31 i=2,3
            j=i+1
            k=i+2
            b(i)=cdsqrt((q2(i)-(mi(i)+mi(j))**2)
     &          *(q2(i)-(mi(i)-mi(j))**2))
            if(dreal(b(i)).ne.0d0) then
              ix0(i,1) = dreal(b(i))
              ix0(i,2) = -ix0(i,1)
            else
              ix0(i,1) =  1d0
              ix0(i,2) = -1d0
            end if
            x0(i,1)=-(q2(i)-mi2(i)+mi2(j)-b(i))/(2d0*q2(i))
            x0(i,2)=-(q2(i)-mi2(i)+mi2(j)+b(i))/(2d0*q2(i))
            x1(i,1)=-(-q2(i)-mi2(i)+mi2(j)-b(i))/(2d0*q2(i))
            x1(i,2)=-(-q2(i)-mi2(i)+mi2(j)+b(i))/(2d0*q2(i))

c>          write(*,*) ' C0reg_coli b(i) ',q2(i),mi2(i),mi2(j),b(i)
c>          write(*,*) ' C0reg_coli x0(i,1) ',i,x0(i,1),ix0(i,1)
c>          write(*,*) ' C0reg_coli x1(i,1) ',i,x1(i,1),ix0(i,1)
c>          write(*,*) ' C0reg_coli x0(i,2) ',i,x0(i,2),ix0(i,2)
c>          write(*,*) ' C0reg_coli x1(i,2) ',i,x1(i,2),ix0(i,2)

            if (abs(x0(i,1)).lt.abs(x0(i,2))) then
              x0(i,1) =  mi2(j)/(q2(i)* x0(i,2))
            else
              x0(i,2) =  mi2(j)/(q2(i)* x0(i,1))
            end if
            if (abs(x1(i,1)).lt.abs(x1(i,2))) then
              x1(i,1) =  mi2(i)/(q2(i)* x1(i,2))
            else
              x1(i,2) =  mi2(i)/(q2(i)* x1(i,1))
            end if 
 31       continue


c>          write(*,*) ' C0reg_coli x0(2,1) ',x0(2,1),ix0(2,1)
c>          write(*,*) ' C0reg_coli x1(2,1) ',x1(2,1),ix0(2,1)
c>          write(*,*) ' C0reg_coli x0(2,2) ',x0(2,2),ix0(2,2)
c>          write(*,*) ' C0reg_coli x1(2,2) ',x1(2,2),ix0(2,2)
c>          write(*,*) ' C0reg_coli x0(3,1) ',x0(3,1),ix0(3,1)
c>          write(*,*) ' C0reg_coli x1(3,1) ',x1(3,1),ix0(3,1)
c>          write(*,*) ' C0reg_coli x0(3,2) ',x0(3,2),ix0(3,2)
c>          write(*,*) ' C0reg_coli x1(3,2) ',x1(3,2),ix0(3,2)
c>
c>          write(*,*) ' C0reg_coli cln_coli ',cln_coli(x1(2,1),ix0(2,1))
c>          write(*,*) ' C0reg_coli cln_coli ',cln_coli(x0(2,1),ix0(2,1))
c>          write(*,*) ' C0reg_coli cln_coli ',cln_coli(x1(2,2),ix0(2,2))
c>          write(*,*) ' C0reg_coli cln_coli ',cln_coli(x0(2,2),ix0(2,2))
c>          write(*,*) ' C0reg_coli cln_coli ',cln_coli(x1(3,1),ix0(3,1))
c>          write(*,*) ' C0reg_coli cln_coli ',cln_coli(x0(3,1),ix0(3,1))
c>          write(*,*) ' C0reg_coli cln_coli ',cln_coli(x1(3,2),ix0(3,2))
c>          write(*,*) ' C0reg_coli cln_coli ',cln_coli(x0(3,2),ix0(3,2))
c>
c>          write(*,*) ' C0reg_coli con ',
c>     &         x1(2,1)*(cln_coli(x1(2,1),ix0(2,1))-cln_coli(x0(2,1),ix0(2,1)))
c>          write(*,*) ' C0reg_coli con ',
c>     &        +x1(2,2)*(cln_coli(x1(2,2),ix0(2,2))-cln_coli(x0(2,2),ix0(2,2)))
c>          write(*,*) ' C0reg_coli con ',
c>     &        x0(3,1)*(cln_coli(x0(3,1),ix0(3,1))-cln_coli(x1(3,1),ix0(3,1)))
c>          write(*,*) ' C0reg_coli con ',
c>     &        x0(3,2)*(cln_coli(x0(3,2),ix0(3,2))-cln_coli(x1(3,2),ix0(3,2)))

          C0reg_coli = 
     &         x1(2,1)
     &        *(cln_coli(x1(2,1),ix0(2,1))-cln_coli(x0(2,1),ix0(2,1)))
     &        +x1(2,2)
     &        *(cln_coli(x1(2,2),ix0(2,2))-cln_coli(x0(2,2),ix0(2,2)))
     &        +x0(3,1)
     &        *(cln_coli(x0(3,1),ix0(3,1))-cln_coli(x1(3,1),ix0(3,1)))
     &        +x0(3,2)
     &        *(cln_coli(x0(3,2),ix0(3,2))-cln_coli(x1(3,2),ix0(3,2)))

c          write(*,*) ' C0reg_coli sum  ',C0reg_coli

          C0reg_coli=C0reg_coli/(mi2(1)-mi2(2))
c          write(*,*) ' C0reg_coli res  ',C0reg_coli
        end if


c 29.07.08: eta-terms added
      else
        if (p23.ne.0d0) then
          mi2(1)=m32
          mi2(2)=m12
          mi2(3)=m22
          mi2(4)=mi2(1)
          mi2(5)=mi2(2)
          swap = mi(1)
          mi(1)=mi(3)
          mi(3)=mi(2)
          mi(2)=swap
          mi(4)=mi(1)
          mi(5)=mi(2)
          q2(1)=p13
          q2(2)=p12
          q2(3)=p23
          q2(4)=q2(1)
          q2(5)=q2(2)
        else if (p12.ne.0d0) then
          mi2(1)=m22
          mi2(2)=m32
          mi2(3)=m12
          mi2(4)=mi2(1)
          mi2(5)=mi2(2)
          swap = mi(1)
          mi(1)=mi(2)
          mi(2)=mi(3)
          mi(3)=swap
          mi(4)=mi(1)
          mi(5)=mi(2)
          q2(1)=p23
          q2(2)=p13
          q2(3)=p12
          q2(4)=q2(1)
          q2(5)=q2(2)
        end if
 
        a=-q2(3)

c        write(*,*) ' C0reg_coli q ',3,q2(3),q2(1),q2(2)
c        write(*,*) ' C0reg_coli m ',3,mi2(3),mi2(1),mi2(2)

        if(dimag(a).ne.0d0) then
          write(*,*) 'olietas required'
          write(99,*) 'olietas required'
        end if
 
c        b(3) = cdsqrt(q2(3)*q2(3)+mi2(3)*mi2(3)+mi2(1)*mi2(1)
c     *         -2d0*(q2(3)*mi2(3)+mi2(3)*mi2(1)+mi2(1)*q2(3)))
        b(3)=cdsqrt((q2(3)-(mi(3)+mi(1))**2)*(q2(3)-(mi(3)-mi(1))**2))
        yh = q2(3)* (q2(3) + 2d0*mi2(2)-mi2(3)-mi2(1))
        y(3,1)=2d0*q2(3)*( mi2(2)-mi2(1))
        y(3,2)=2d0*q2(3)*( q2(3)+mi2(2)-mi2(1))
        y(3,3)=yh+a*b(3)
        y(3,4)=yh-a*b(3)
        yh0 = q2(3)*q2(3)*mi2(2)
     *           + q2(3)*(mi2(2)-mi2(3))*(mi2(2)-mi2(1))
        if(cdabs(y(3,3)).lt.cdabs(y(3,4))) then
           y(3,3) = 4d0*yh0*q2(3)/y(3,4)
        else
           y(3,4) = 4d0*yh0*q2(3)/y(3,3)
        end if
c       write(*,*) ' C0reg_coli y ',3,y(3,1),y(3,2),y(3,3),y(3,4)

c>        write(*,*) ' C0reg_check b ',3,yh,a,b(3),a*b(3)
c>        write(*,*) ' C0reg_check y ',3,y(3,1),y(3,2),y(3,3),y(3,4)
c>        write(*,*) ' C0reg_check y/y',3,
c>     &      y(3,2)/y(3,3),y(3,2)/y(3,4)
c>     &      ,y(3,1)/y(3,3),y(3,1)/y(3,4)

        if (mi2(2).ne.mi2(3)) then

#ifdef CHECK
          if(flag(4))then
            write(*,100) ' C0reg_coli: p12=0=p23, m22=/=m32. '
            write(*,100) '    case not yet tested in physical process'
            write(*,111) ' C0reg_coli: p12 = ',p12
            write(*,111) ' C0reg_coli: p23 = ',p23
            write(*,111) ' C0reg_coli: p13 = ',p13
            write(*,111) ' C0reg_coli: m12 = ',m12
            write(*,111) ' C0reg_coli: m22 = ',m22
            write(*,111) ' C0reg_coli: m32 = ',m32
            flag(4)=.false.
          endif
#endif

          y(2,1) = (-q2(3)+mi2(1)-mi2(2)) / a
          y(2,2) = ( mi2(1)-mi2(2)) / a
          y(2,3) = (mi2(1)-mi2(2)-q2(3)*mi2(2)/(mi2(2)-mi2(3))) / a
          iy(2,3) = -1d0/(mi2(2)-mi2(3))
          x1(2,1) = y(2,3)-y(2,2)
          x1(2,1) = mi2(2)/(mi2(2)-mi2(3))
          x0(2,1) = y(2,3)-y(2,1)
          x0(2,1) = mi2(3)/(mi2(2)-mi2(3))

c>          write(*,*) ' C0reg_coli test',y(2,1),y(2,2)+1,
c>     &        y(2,3),y(2,1)+x0(2,1),
c>     &        x1(2,1),1+x0(2,1)

          if(y(2,2).ne.0d0)then
            iy2y3=iy(2,3)/y(2,2)
          else
            iy2y3=0d0
          endif
          if(y(2,2).ne.0d0)then
            iy2y4=iy(2,4)/y(2,2)
          else
            iy2y4=0d0
          endif

          C0reg_coli = 
     &       cspcon_coli(x1(2,1),1d0/y(2,3),x1(2,1)/y(2,3),y(2,2)/y(2,3)
     &      ,iy(2,3),-iy(2,3))
     &      -cspcon_coli(x0(2,1),1d0/y(2,3),x0(2,1)/y(2,3),y(2,1)/y(2,3)
     &      ,iy(2,3),-iy(2,3))

c          write(*,*) ' C0reg_coli s2',
c     &       cspcon_coli(x1(2,1),1d0/y(2,3),x1(2,1)/y(2,3),y(2,2)/y(2,3)
c     &      ,iy(2,3),-iy(2,3))
c     &      -cspcon_coli(x0(2,1),1d0/y(2,3),x0(2,1)/y(2,3),y(2,1)/y(2,3)
c     &      ,iy(2,3),-iy(2,3))

          deleta = eta2_coli(mi2(2)-mi2(3),x0(2,1),mi2(3)
     &        ,-1d0,iy(2,3),-1d0)
     &        - eta2_coli(mi2(2)-mi2(3),y(2,3),(mi2(2)-mi2(3))*y(2,3)
     &        ,-1d0,iy(2,3),-1d0)

c       write(*,*) ' C0reg_coli deleta',deleta

          if(deleta.ne.0d0)then
            C0reg_coli = C0reg_coli
     &          -deleta*cln_coli(y(2,2)/y(2,1),1d0)

c            write(*,*) ' C0reg_coli es 2',
c     &          -deleta*cln_coli(y(2,2)/y(2,1),1d0)

          endif

c          write(*,*) ' C0reg_coli s2',C0reg_coli

        else

#ifdef ALLCHECK
c ppttbb 
          if(flag(5))then
            write(*,100) ' C0reg_coli: p12=0=p23, m22=m32 '
            write(*,100) '    case not yet tested in physical process'
            write(*,111) ' C0reg_coli: p12 = ',p12
            write(*,111) ' C0reg_coli: p23 = ',p23
            write(*,111) ' C0reg_coli: p13 = ',p13
            write(*,111) ' C0reg_coli: m12 = ',m12
            write(*,111) ' C0reg_coli: m22 = ',m22
            write(*,111) ' C0reg_coli: m32 = ',m32
            flag(5)=.false.
          endif
#endif

        endif

        y(3,1) = ( mi2(2)-mi2(1)) / a
        y(3,2) = ( q2(3)-q2(2)+mi2(2)-mi2(1)) / a
 
        b(3)=cdsqrt((q2(3)-(mi(3)+mi(1))**2)*(q2(3)-(mi(3)-mi(1))**2))
        yh = q2(3)*(q2(3)+2d0*mi2(2)-mi2(3)-mi2(1))
        y(3,3) = (yh+a*b(3))/(2d0*a*q2(3))
        y(3,4) = (yh-a*b(3))/(2d0*a*q2(3))
c>        write(*,*) ' C0reg_coli q ',i,q2(i),q2(j),q2(k)
c>        write(*,*) ' C0reg_coli m ',i,mi(i),mi(j),mi(k)
c>        write(*,*) ' C0reg_coli m2',i,mi2(i),mi2(j),mi2(k)
c>        write(*,*) ' C0reg_coli b ',i,yh,a,b(i),a*b(i)
c>        write(*,*) ' C0reg_coli y ',i,y(i,1),y(i,2),y(i,3),y(i,4)
        yh0 = q2(3)*q2(3)*mi2(2)
     &           + q2(3)*(mi2(2)-mi2(3))*(mi2(2)-mi2(1))
        yh3 = yh0 
        if(cdabs(y(3,3)).lt.cdabs(y(3,4))) then
          y(3,3) = yh3/(a*a*q2(3)*y(3,4))
        else
          y(3,4) = yh3/(a*a*q2(3)*y(3,3))
        end if
        if(dreal(b(3)).ne.0d0) then
          iy(3,3) = dreal(b(3))
          iy(3,4) = -iy(3,3)
        else
          iy(3,3) =  1d0
          iy(3,4) = -1d0
        end if
        x0(3,1)=-(q2(3)-mi2(3)+mi2(1)-b(3))/(2d0*q2(3))
        x0(3,2)=-(q2(3)-mi2(3)+mi2(1)+b(3))/(2d0*q2(3))
        x1(3,1)=-(-q2(3)-mi2(3)+mi2(1)-b(3))/(2d0*q2(3))
        x1(3,2)=-(-q2(3)-mi2(3)+mi2(1)+b(3))/(2d0*q2(3))

c>        write(*,*) ' C0reg_coli test',y(3,1),y(3,2)+1,
c>     &        y(3,3),y(3,1)+x0(3,1),
c>     &        y(3,4),y(3,1)+x0(3,2),
c>     &        x1(3,1),1+x0(3,1),
c>     &        x1(3,2),1+x0(3,2)
c>        write(*,*) ' C0reg_coli iy',iy(3,3),iy(3,4),b(3)

        if (abs(x0(3,1)).lt.abs(x0(3,2))) then
          x0(3,1) =  mi2(1)/(q2(3)* x0(3,2))
        else
          x0(3,2) =  mi2(1)/(q2(3)* x0(3,1))
        end if
        if (abs(x1(3,1)).lt.abs(x1(3,2))) then
          x1(3,1) =  mi2(3)/(q2(3)* x1(3,2))
        else
          x1(3,2) =  mi2(3)/(q2(3)* x1(3,1))
        end if 

        if(y(3,2).ne.0d0)then
          iy2y3=iy(3,3)/y(3,2)
        else
          iy2y3=0d0
        endif
        if(y(3,2).ne.0d0)then
          iy2y4=iy(3,4)/y(3,2)
        else
          iy2y3=0d0
        endif
        if(y(3,1).ne.0d0)then
          iy1y3=iy(3,3)/y(3,1)
        else
          iy1y3=0d0
        endif
        if(y(3,1).ne.0d0)then
          iy1y4=iy(3,4)/y(3,1)
        else
          iy1y4=0d0
        endif
            
c>          write(*,*) ' C0reg_coli args1',
c>     &         x1(3,1),1d0/y(3,3),x1(3,1)/y(3,3),y(3,2)/y(3,3)
c>     &                          ,iy(3,3),-iy(3,3),-iy2y3
c>          write(*,*) ' C0reg_coli args2',
c>     &         x1(3,2),1d0/y(3,4),x1(3,2)/y(3,4),y(3,2)/y(3,4)
c>     &                          ,iy(3,4),-iy(3,4)
c>          write(*,*) ' C0reg_coli args3',
c>     &        x0(3,1),1d0/y(3,3),x0(3,1)/y(3,3),y(3,1)/y(3,3)
c>     &                          ,iy(3,3),-iy(3,3)
c>          write(*,*) ' C0reg_coli args4',
c>     &         x0(3,2),1d0/y(3,4),x0(3,2)/y(3,4),y(3,1)/y(3,4)
c>     &                          ,iy(3,4),-iy(3,4)
c>
c>        write(*,*) ' C0reg_coli s1',
c>     &      +cspcon_coli(x1(3,1),1d0/y(3,3),x1(3,1)/y(3,3),y(3,2)/y(3,3)
c>     &      ,iy(3,3),-iy(3,3),-iy2y3),cspenc_coli(y(3,2)/y(3,3), iy2y3)
c>        write(*,*) ' C0reg_coli s2',
c>     &      +cspcon_coli(x1(3,2),1d0/y(3,4),x1(3,2)/y(3,4),y(3,2)/y(3,4)
c>     &      ,iy(3,4),-iy(3,4),-iy1y3),+cspenc_coli(y(3,2)/y(3,4), iy2y4)
c>        write(*,*) ' C0reg_coli s3',
c>     &      -cspcon_coli(x0(3,1),1d0/y(3,3),x0(3,1)/y(3,3),y(3,1)/y(3,3)
c>     &      ,iy(3,3),-iy(3,3),-iy2y4),-cspenc_coli(y(3,1)/y(3,3), iy1y3)
c>        write(*,*) ' C0reg_coli s4',
c>     &      -cspcon_coli(x0(3,2),1d0/y(3,4),x0(3,2)/y(3,4),y(3,1)/y(3,4)
c>     &      ,iy(3,4),-iy(3,4),-iy1y4),-cspenc_coli(y(3,1)/y(3,4), iy1y4)
c>        write(*,*) ' C0reg_coli s5',
c>     &      -(eta2_coli(x0(3,1),x0(3,2),mi2(1)/q2(3)
c>     &      ,iy(3,3),iy(3,4),-dreal(1d0/q2(3)))
c>     &      - eta2_coli(y(3,3),y(3,4), yh3/(a*a*q2(3))
c>     &      ,iy(3,3),iy(3,4),-dreal(1d0/q2(3))))
c>     &      * cln_coli(y(3,2)/y(3,1),1d0)
        
        C0reg_coli = C0reg_coli
     &      +cspcon_coli(x1(3,1),1d0/y(3,3),x1(3,1)/y(3,3),y(3,2)/y(3,3)
     &      ,iy(3,3),-iy(3,3))
     &      +cspcon_coli(x1(3,2),1d0/y(3,4),x1(3,2)/y(3,4),y(3,2)/y(3,4)
     &      ,iy(3,4),-iy(3,4))
     &      -cspcon_coli(x0(3,1),1d0/y(3,3),x0(3,1)/y(3,3),y(3,1)/y(3,3)
     &      ,iy(3,3),-iy(3,3))
     &      -cspcon_coli(x0(3,2),1d0/y(3,4),x0(3,2)/y(3,4),y(3,1)/y(3,4)
     &      ,iy(3,4),-iy(3,4))

c>        write(*,*)   ' C0reg_coli s 3',      
c>     &      +cspcon_coli(x1(3,1),1d0/y(3,3),x1(3,1)/y(3,3),y(3,2)/y(3,3)
c>     &      ,iy(3,3),-iy(3,3))
c>     &      +cspcon_coli(x1(3,2),1d0/y(3,4),x1(3,2)/y(3,4),y(3,2)/y(3,4)
c>     &      ,iy(3,4),-iy(3,4))
c>     &      -cspcon_coli(x0(3,1),1d0/y(3,3),x0(3,1)/y(3,3),y(3,1)/y(3,3)
c>     &      ,iy(3,3),-iy(3,3))
c>     &      -cspcon_coli(x0(3,2),1d0/y(3,4),x0(3,2)/y(3,4),y(3,1)/y(3,4)
c>     &      ,iy(3,4),-iy(3,4))

        deleta = eta2_coli(x0(3,1),x0(3,2),mi2(1)/q2(3)
     &      ,iy(3,3),iy(3,4),-dreal(1d0/q2(3)))
     &      - eta2_coli(y(3,3),y(3,4), yh3/(a*a*q2(3))
     &      ,iy(3,3),iy(3,4),-dreal(1d0/q2(3)))

c       write(*,*) ' C0reg_coli deleta',deleta

        if(deleta.ne.0d0)then
          C0reg_coli = C0reg_coli
     &        -deleta*cln_coli(y(3,2)/y(3,1),1d0)

c          write(*,*)   ' C0reg_coli es 3',      
c     &        -deleta*cln_coli(y(3,2)/y(3,1),1d0)

        endif

        if(dimag(yh3/(a*a*q2(3))).lt.0d0.and.dreal(q2(3)).lt.0d0)
     &      C0reg_coli=C0reg_coli+2d0*pi*dcmplx(0d0,1d0)
     &      *cln_coli(y(3,2)/y(3,1),1d0)

c            write(*,*) ' C0reg_coli etaterm', y(3,3),y(3,4), 
c     &          y(3,3)*y(3,4),yh3/(a*a*q2(3)),y(3,2),y(3,1)

        
c        write(*,*) ' C0reg_coli s',C0reg_coli

        C0reg_coli=C0reg_coli/a

#ifdef CHECK
        if(conscheck)then

          if (mi2(2).ne.mi2(3)) then
            y(2,1) = (-q2(3)+mi2(1)-mi2(2)) / a
            y(2,2) = ( mi2(1)-mi2(2)) / a
            y(2,3) = (mi2(1)-mi2(2)-q2(3)*mi2(2)/(mi2(2)-mi2(3))) / a
            iy(2,3) = -1d0/(mi2(2)-mi2(3))
            x1(2,1) = y(2,3)-y(2,2)
            x1(2,1) = mi2(2)/(mi2(2)-mi2(3))
            x0(2,1) = y(2,3)-y(2,1)
            x0(2,1) = mi2(3)/(mi2(2)-mi2(3))
            
            if(y(2,2).ne.0d0)then
              iy2y3=iy(2,3)/y(2,2)
            else
              iy2y3=0d0
            endif
            if(y(2,1).ne.0d0)then
              iy1y3=iy(2,3)/y(2,1)
            else
              iy1y3=0d0
            endif
            if(y(2,2).ne.0d0)then
              iy2y4=iy(2,4)/y(2,2)
            else
              iy2y3=0d0
            endif

c>          write(*,*) ' C0reg_coli test',y(2,1),y(2,2)+1,
c>     &        y(2,3),y(2,1)+x0(2,1),
c>     &        x1(2,1),1+x0(2,1)

            C0reg_check=cspenc_coli(y(2,2)/y(2,3),iy2y3)
     &          -cspenc_coli(y(2,1)/y(2,3),iy1y3)

c            write(*,*) ' C0reg_check 2',C0reg_check

            eta = eta2s_coli(x1(2,1),1d0/y(2,3),iy(2,3),-iy(2,3),-iy2y3)
            if(eta.ne.cd0)then
              C0reg_check = C0reg_check 
     &            +eta*cln_coli(y(2,2)/y(2,3),iy2y3)

c              write(*,*) ' C0reg_check se 2a',
c     &            +eta*cln_coli(y(2,2)/y(2,3),iy2y3)

            endif            
            eta = - eta2s_coli(x0(2,1),1d0/y(2,3),iy(2,3),-iy(2,3),
     &          -iy1y3)
            if(eta.ne.cd0)then
              C0reg_check = C0reg_check 
     &            +eta*cln_coli(y(2,1)/y(2,3),iy1y3)

c              write(*,*) ' C0reg_check se 2b',
c     &            +eta*cln_coli(y(2,1)/y(2,3),iy1y3)

            endif            
            eta = -( eta2_coli(mi2(2)-mi2(3),x0(2,1),mi2(3)
     &          ,-1d0,iy(2,3),-1d0)
     &          - eta2_coli(mi2(2)-mi2(3),y(2,3),(mi2(2)-mi2(3))*y(2,3)
     &          ,-1d0,iy(2,3),-1d0))
            if(eta.ne.cd0)then
              C0reg_check = C0reg_check 
     &            +eta* cln_coli(y(2,2)/y(2,1),1d0)

c              write(*,*) ' C0reg_check se 2c',
c     &            +eta* cln_coli(y(2,2)/y(2,1),1d0)

            endif            

          endif

c          write(*,*) ' C0reg_check aft eta2',C0reg_check


          y(3,1) = ( mi2(2)-mi2(1)) / a
          y(3,2) = ( q2(3)-q2(2)+mi2(2)-mi2(1)) / a
          
          b(3)=cdsqrt((q2(3)-(mi(3)+mi(1))**2)*(q2(3)-(mi(3)-mi(1))**2))
          yh = q2(3)*(q2(3)+2d0*mi2(2)-mi2(3)-mi2(1))
          y(3,3) = (yh+a*b(3))/(2d0*a*q2(3))
          y(3,4) = (yh-a*b(3))/(2d0*a*q2(3))
c>        write(*,*) ' C0reg_coli q ',i,q2(i),q2(j),q2(k)
c>        write(*,*) ' C0reg_coli m ',i,mi(i),mi(j),mi(k)
c>        write(*,*) ' C0reg_coli m2',i,mi2(i),mi2(j),mi2(k)
c>        write(*,*) ' C0reg_coli b ',i,yh,a,b(i),a*b(i)
c>        write(*,*) ' C0reg_coli y ',i,y(i,1),y(i,2),y(i,3),y(i,4)
          yh0 = q2(3)*q2(3)*mi2(2)
     &        + q2(3)*(mi2(2)-mi2(3))*(mi2(2)-mi2(1))
          yh3 = yh0 
          if(cdabs(y(3,3)).lt.cdabs(y(3,4))) then
            y(3,3) = yh3/(a*a*q2(3)*y(3,4))
          else
            y(3,4) = yh3/(a*a*q2(3)*y(3,3))
          end if
          if(dreal(b(3)).ne.0d0) then
            iy(3,3) = dreal(b(3))
            iy(3,4) = -iy(3,3)
          else
            iy(3,3) =  1d0
            iy(3,4) = -1d0
          end if
          if(y(3,2).ne.0d0)then
            iy2y3=iy(3,3)/y(3,2)
          else
            iy2y3=0d0
          endif
          if(y(3,2).ne.0d0)then
            iy2y4=iy(3,4)/y(3,2)
          else
            iy2y3=0d0
          endif
          if(y(3,1).ne.0d0)then
            iy1y3=iy(3,3)/y(3,1)
          else
            iy1y3=0d0
          endif
          if(y(3,1).ne.0d0)then
            iy1y4=iy(3,4)/y(3,1)
          else
            iy1y4=0d0
          endif

c       write(*,*) ' C0reg_coli y ',3,y(3,1),y(3,2),y(3,3),y(3,4)
 
          C0reg_check=C0reg_check
     &        +cspenc_coli(y(3,2)/y(3,3), iy2y3)
     &        +cspenc_coli(y(3,2)/y(3,4), iy2y4)
     &        -cspenc_coli(y(3,1)/y(3,3), iy1y3)
     &        -cspenc_coli(y(3,1)/y(3,4), iy1y4)

c          write(*,*) ' C0reg_check ',3,C0reg_check

          x0(3,1)=-(q2(3)-mi2(3)+mi2(1)-b(3))/(2d0*q2(3))
          x0(3,2)=-(q2(3)-mi2(3)+mi2(1)+b(3))/(2d0*q2(3))
          x1(3,1)=-(-q2(3)-mi2(3)+mi2(1)-b(3))/(2d0*q2(3))
          x1(3,2)=-(-q2(3)-mi2(3)+mi2(1)+b(3))/(2d0*q2(3))

c>          write(*,*) ' C0reg_coli test',y(3,1),y(3,2)+1,
c>     &        y(3,3),y(3,1)+x0(3,1),
c>     &        y(3,4),y(3,1)+x0(3,2),
c>     &        x1(3,1),1+x0(3,1),
c>     &        x1(3,2),1+x0(3,2)

          if (abs(x0(3,1)).lt.abs(x0(3,2))) then
            x0(3,1) =  mi2(1)/(q2(3)* x0(3,2))
          else
            x0(3,2) =  mi2(1)/(q2(3)* x0(3,1))
          end if
          if (abs(x1(3,1)).lt.abs(x1(3,2))) then
            x1(3,1) =  mi2(3)/(q2(3)* x1(3,2))
          else
            x1(3,2) =  mi2(3)/(q2(3)* x1(3,1))
          end if 

          eta =  eta2s_coli(x1(3,1),1d0/y(3,3),iy(3,3),-iy(3,3),-iy2y3)
          if(eta.ne.cd0)then
            C0reg_check = C0reg_check 
     &          +eta*cln_coli(y(3,2)/y(3,3),iy2y3)
            
c            write(*,*) ' C0reg_check se 3a ',
c     &          +eta*cln_coli(y(3,2)/y(3,3),iy2y3)
            
          endif
          eta =  - eta2s_coli(x0(3,1),1d0/y(3,3),iy(3,3),-iy(3,3),
     &        -iy1y3)
          if(eta.ne.cd0)then
            C0reg_check = C0reg_check 
     &          +eta*cln_coli(y(3,1)/y(3,3),iy1y3)
            
c            write(*,*) ' C0reg_check se 3b ',
c     &          +eta*cln_coli(y(3,1)/y(3,3),iy1y3)
            
          endif
          eta = + eta2s_coli(x1(3,2),1d0/y(3,4),iy(3,4),-iy(3,4),-iy2y4)
          if(eta.ne.cd0)then
            C0reg_check = C0reg_check 
     &          +eta*cln_coli(y(3,2)/y(3,4),iy2y4)
            
c            write(*,*) ' C0reg_check se 3c ',
c     &          +eta*cln_coli(y(3,2)/y(3,4),iy2y4)
            
          endif
          eta = - eta2s_coli(x0(3,2),1d0/y(3,4),iy(3,4),-iy(3,4),-iy1y4)
          if(eta.ne.cd0)then
            C0reg_check = C0reg_check 
     &          +eta*cln_coli(y(3,1)/y(3,4),iy1y4)
            
c            write(*,*) ' C0reg_check se 3d ',
c     &          +eta*cln_coli(y(3,1)/y(3,4),iy1y4)
            
          endif
          
          eta=-(eta2s_coli(x0(3,1),x0(3,2),iy(3,3),iy(3,4),
     &        -dreal(1d0/q2(3)))
     &        - eta2s_coli(y(3,3),y(3,4),iy(3,3),iy(3,4),
     &        -dreal(1d0/q2(3))))
          if(eta.ne.cd0)then
            C0reg_check = C0reg_check +
     &          eta*cln_coli(y(3,2)/y(3,1),1d0)
            
c>            write(*,*) ' C0reg_check se 3e',
c>     &         -(eta2s_coli(x0(3,1),x0(3,2),iy(3,3),iy(3,4),
c>     &          -dreal(1d0/q2(3)))
c>     &          - eta2s_coli(y(3,3),y(3,4),iy(3,3),iy(3,4),
c>     &          -dreal(1d0/q2(3))))
c>     &          * cln_coli(y(3,2)/y(3,1),1d0)
            
          endif 
          
c            write(*,*) ' C0reg_check etaterm', y(3,3),y(3,4), 
c     &          y(3,3)*y(3,4),yh3/(a*a*q2(3)),y(3,2),y(3,1)

          if(dimag(yh3/(a*a*q2(3))).lt.0d0.and.dreal(q2(3)).lt.0d0) then
            C0reg_check=C0reg_check+2d0*pi*dcmplx(0d0,1d0)
     &          *cln_coli(y(3,2)/y(3,1),1d0)
            
c            write(*,*) ' C0reg_check se 3f', +2d0*pi*dcmplx(0d0,1d0)
c     &          *cln_coli(y(3,2)/y(3,1),1d0)
            
          endif
          
c          write(*,*) ' C0reg_check ',3,C0reg_check

          C0reg_check=C0reg_check/a

          if(cdabs(C0reg_coli/C0reg_check-1d0).gt.testacc) then
            write(testout,100) ' C0reg_coli inconsistent'
            write(testout,111) ' C0reg_coli: p12 = ',p12
            write(testout,111) ' C0reg_coli: p23 = ',p23
            write(testout,111) ' C0reg_coli: p13 = ',p13
            write(testout,111) ' C0reg_coli: m12 = ',m12
            write(testout,111) ' C0reg_coli: m22 = ',m22
            write(testout,111) ' C0reg_coli: m32 = ',m32
            write(testout,111) ' C0reg_coli  = ',C0reg_coli
            write(testout,111) ' C0reg_check = ',C0reg_check
            write(testout,101) ' diff       = ',
     &          cdabs(C0reg_coli/C0reg_check-1d0)
            
            write(*,100) ' C0reg_coli inconsistent'
            write(*,111) ' C0reg_coli: p12 = ',p12
            write(*,111) ' C0reg_coli: p23 = ',p23
            write(*,111) ' C0reg_coli: p13 = ',p13
            write(*,111) ' C0reg_coli: m12 = ',m12
            write(*,111) ' C0reg_coli: m22 = ',m22
            write(*,111) ' C0reg_coli: m32 = ',m32
            write(*,111) ' C0reg_coli  = ',C0reg_coli
            write(*,111) ' C0reg_check = ',C0reg_check
            write(*,111) ' diff        = ',
     &          cdabs(C0reg_coli/C0reg_check-1d0)
c            stop
          endif
        endif
#endif
      endif

c        write(*,*) ' C0reg_check= ',C0reg_check,C0reg_check*a
c        write(*,*) ' C0reg_coli=  ',C0reg_coli,C0reg_coli*a
       

      end

************************************************************************
      function C0m0_coli(p12,p23,p13,m12,m22,m32)
************************************************************************
*     regular scalar 3-point function with m12 = 0                     *
*     assumes no mass singularities and no ir singularities            *
*                                                                      *
*                            p23                                       *
*                             |                                        *
*                             |                                        *
*                            / \                                       *
*                           /   \                                      *
*                   m22    /     \   m32                               *
*                         /       \                                    *
*                        /         \                                   *
*              p12  ---------------------  p13                         *
*                          m12 = 0                                     *
*                                                                      *
*----------------------------------------------------------------------*
*     13.05.04 Ansgar Denner        last changed  26.09.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p13,m12,m22,m32
      complex*16 C0m0_coli

      complex*16 m2,m3
      complex*16 k
      complex*16 r,det,rp,detp
      real*8     ir,irp
      complex*16 mm22,mm32,q13,q12,sq12,sq13
      real*8     ix1,iz1,ix3
      complex*16 x1,x2,z1,z2,x3,invx3
      complex*16 cspenc_coli,cspcos_coli,cln_coli
      complex*16 C0m0_check

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps13,ms12,ms22,ms32
      logical flag(0:8)
      data    flag /9*.true./
      save    flag  
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) ' C0m0_coli in ',p12,p23,p13,m12,m22,m32
 
#ifdef CHECK
 100  format((a))
 101  format(a22,g24.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        ms12 = elimminf2_coli(m12)
        ms22 = elimminf2_coli(m22)
        ms32 = elimminf2_coli(m32)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps13 = elimminf2_coli(p13)

        if(m12.ne.cd0.or.ps12.eq.ms22.and.ps13.eq.ms32
     &        .or.ps23.eq.cd0.and.ms32.eq.cd0.and.ms22.eq.cd0) then
          write(*,100) ' C0m0_coli called improperly'
          write(*,100) '     or case not implemented'
          write(*,111) ' C0m0_coli: p12 = ',p12,ps12
          write(*,111) ' C0m0_coli: p23 = ',p23,ps23
          write(*,111) ' C0m0_coli: p13 = ',p13,ps13
          write(*,111) ' C0m0_coli: m12 = ',m12,ms12
          write(*,111) ' C0m0_coli: m22 = ',m22,ms22
          write(*,111) ' C0m0_coli: m32 = ',m32,ms32
          stop
        end if
      end if
#endif        

      if (p12.ne.cd0.and.p13.ne.cd0.and.p23.ne.cd0) then

#ifdef ALLCHECK
c ppttbb
        if(flag(1))then
          write(*,100) ' C0m0_coli: p12,p13,p23=/=0 '
          write(*,100) '    case not yet tested in physical process'
          write(*,111) ' C0m0_coli: p12 = ',p12
          write(*,111) ' C0m0_coli: p23 = ',p23
          write(*,111) ' C0m0_coli: p13 = ',p13
          write(*,111) ' C0m0_coli: m12 = ',m12
          write(*,111) ' C0m0_coli: m22 = ',m22
          write(*,111) ' C0m0_coli: m32 = ',m32
          flag(1)=.false.
        endif
#endif
          
        if (m32.ne.p13) then
          mm32 = m32
          mm22 = m22
          q12 = p12
          q13 = p13
        else
          mm32 = m22
          mm22 = m32
          q12 = p13
          q13 = p12
        end if          

        m2 = sqrt(mm22)
        m3 = sqrt(mm32)
        det = sqrt((p23-(m2+m3)*(m2+m3))*(p23-(m2-m3)*(m2-m3))) 

        k = (mm32+mm22-p23)/(m2*m3)
        det = sqrt((p23-(m2+m3)*(m2+m3))*
     &      (p23-(m2-m3)*(m2-m3))) * dsign(1d0,dreal(k))
        r = 0.5d0*k + det/(2d0*m2*m3)

c          write(*,*) ' C0m0_coli k  = ',k,r,dreal(k),dreal(k).le.-2d0

        if(dreal(k).le.-2d0+calacc) then
          ir = sign(1d0,1d0-abs(r))
          if(ir.eq.0d0) ir=1d0            ! check
        else
          ir = 0D0
        end if
        sq13 = sqrt(q13)
        sq12 = sqrt(q12)
        if (dreal(q12).ge.0d0.and.dreal(q13).ge.0d0) then
          detp = sqrt((p23-(sq12+sq13)*(sq12+sq13))*
     &        (p23-(sq12-sq13)*(sq12-sq13)))
        else 
          detp = sqrt((q13-p23)**2 + q12*q12 - 2d0*q12*(q13+p23))
        end if

c>        write(*,*) ' C0m0_coli: ps ',q13,q12,p23
c>        write(*,*) ' C0m0_coli: detp= ',detp,
c>     &    sqrt((p23-(sq12+sq13)*(sq12+sq13))*
c>     &        (p23-(sq12-sq13)*(sq12-sq13))),
c>     &      sqrt((q13-p23)**2 + q12*q12 - 2d0*q12*(q13+p23))


c>        kp = (q13+q12-p23)/(sq12*sq13)
c>        rp = 0.5d0*kp 
c>     &      + detp* dsign(1d0,dreal(q13+q12-p23))/(2d0*sq12*sq13)
c>
c>        if (kp.ne.0d0) then
c>          rp = kp/2d0*(1d0+sqrt(dcmplx(1d0-4d0/kp**2)))
c>        else
c>          rp = dcmplx(0d0,1d0)
c>        end if

        detp=detp* dsign(1d0,dreal(q13+q12-p23))
        rp=((q13+q12-p23)+detp)/(2d0*sq12*sq13)
        irp =-1d0

c>        write(*,*) 'C0m0 detp ',detp,detp* dsign(1d0,dreal(q13+q12-p23))
c>
c>        write(*,*) ' C0m0_coli ir  = ',ir,irp
c>        write(*,*) ' C0m0_coli r,rp= ',r,rp
c>        write(*,*) ' C0m0_coli r,rp= ',m2*r/m3,sq13*rp/sq12
c>        write(*,*) ' C0m0_coli arg5= ',
c>     &      m2*r/m3,sq13*rp/sq12,ir,irp
c>        write(*,*) ' C0m0_coli arg6= ',
c>     &      m2/(m3*r),sq13/(sq12*rp),-ir,-irp
c>        write(*,*) ' C0m0_coli arg7= ',
c>     &      m2*r/m3,sq13/(sq12*rp),ir,-irp
c>        write(*,*) ' C0m0_coli arg8= ',(mm22-q12)/(mm32-q13),
c>     &      sq13/(sq12*rp),dreal(mm22-q12-mm32+q13),-irp
c>        write(*,*) ' C0m0_coli arg9= ',(mm22-q12)/(mm32-q13)
c>     &      ,sq13*rp/sq12,dreal(mm22-q12-mm32+q13),irp
c>        
c>        write(*,*) 'C0m0_coli c1 ',cspenc_coli(1d0-sq13/(sq12*rp),irp) 
c>        write(*,*) 'C0m0_coli c2 ',
c>     &      - cspenc_coli(1d0-sq13*rp/sq12,-irp)
c>        write(*,*) 'C0m0_coli c3 ',
c>     &      +(cln_coli(sq13*rp/sq12,irp)
c>     &      -cln_coli(sq13/(sq12*rp),-irp))
c>     &      *cln_coli((mm32-q13)/(mm32),-1d0)
c>        write(*,*) 'C0m0_coli c4 ',
c>     &      + cspcos_coli(m2/(m3*r),sq13*rp/sq12,-ir,irp)
c>        write(*,*) 'C0m0_coli c5 ',
c>     &      + cspcos_coli(m2*r/m3,sq13*rp/sq12,ir,irp)
c>        write(*,*) 'C0m0_coli c6 ',
c>     &      - cspcos_coli(m2/(m3*r),sq13/(sq12*rp),-ir,-irp)
c>        write(*,*) 'C0m0_coli c7 ',
c>     &      - cspcos_coli(m2*r/m3,sq13/(sq12*rp),ir,-irp)
c>        write(*,*) 'C0m0_coli c8 ',
c>     &      + cspcos_coli((mm22-q12)/(mm32-q13),sq13/(sq12*rp)
c>     &      ,dreal(mm22-q12-mm32+q13),-irp)
c>        write(*,*) 'C0m0_coli c9 ',
c>     &      - cspcos_coli((mm22-q12)/(mm32-q13),sq13*rp/sq12
c>     &      ,dreal(mm22-q12-mm32+q13),irp)
c>
c>        write(*,*) 'C0m0_coli c ',cspenc_coli(1d0-sq13/(sq12*rp),irp), 
c>     &      - cspenc_coli(1d0-sq13*rp/sq12,-irp),
c>     &      +(cln_coli(sq13*rp/sq12,irp)
c>     &      -cln_coli(sq13/(sq12*rp),-irp))
c>     &      *cln_coli((mm32-q13)/(mm32),-1d0),
c>     &      + cspcos_coli(m2/(m3*r),sq13*rp/sq12,-ir,irp),
c>     &      + cspcos_coli(m2*r/m3,sq13*rp/sq12,ir,irp),
c>     &      - cspcos_coli(m2/(m3*r),sq13/(sq12*rp),-ir,-irp),
c>     &      - cspcos_coli(m2*r/m3,sq13/(sq12*rp),ir,-irp),
c>     &      + cspcos_coli((mm22-q12)/(mm32-q13),sq13/(sq12*rp)
c>     &      ,dreal(mm22-q12-mm32+q13),-irp),
c>     &      - cspcos_coli((mm22-q12)/(mm32-q13),sq13*rp/sq12
c>     &      ,dreal(mm22-q12-mm32+q13),irp)
 
        C0m0_coli = cspenc_coli(1d0-sq13/(sq12*rp),irp) 
     &      - cspenc_coli(1d0-sq13*rp/sq12,-irp)
     &      +(cln_coli(sq13*rp/sq12,irp)
     &      -cln_coli(sq13/(sq12*rp),-irp))
     &      *cln_coli((mm32-q13)/(mm32),-1d0)
     &      + cspcos_coli(m2/(m3*r),sq13*rp/sq12,-ir,irp)
     &      + cspcos_coli(m2*r/m3,sq13*rp/sq12,ir,irp)
     &      - cspcos_coli(m2/(m3*r),sq13/(sq12*rp),-ir,-irp)
     &      - cspcos_coli(m2*r/m3,sq13/(sq12*rp),ir,-irp)
     &      + cspcos_coli((mm22-q12)/(mm32-q13),sq13/(sq12*rp)
     &      ,dreal(mm22-q12-mm32+q13),-irp)
     &        - cspcos_coli((mm22-q12)/(mm32-q13),sq13*rp/sq12
     &      ,dreal(mm22-q12-mm32+q13),irp)

c>      write(*,*) ' C0m0_coli s',
c>     &      cspenc_coli(1d0-sq13/(sq12*rp),irp) 
c>     &      - cspenc_coli(1d0-sq13*rp/sq12,-irp)
c>     &        +(cln_coli(sq13*rp/sq12,irp)
c>     &        -cln_coli(sq13/(sq12*rp),irp))
c>     &        *cln_coli((mm32-q13)/(mm32),-1d0)
c>     &      + cspcos_coli(m2/(m3*r),sq13*rp/sq12,-ir,irp)
c>     &      + cspcos_coli(m2*r/m3,sq13*rp/sq12,ir,irp)
c>     &      - cspcos_coli(m2/(m3*r),sq13/(sq12*rp),-ir,-irp)
c>     &      - cspcos_coli(m2*r/m3,sq13/(sq12*rp),ir,-irp)
c>     &      + cspcos_coli((mm22-q12)/(mm32-q13),sq13/(sq12*rp)
c>     &      ,dreal(mm22-q12-mm32+q13),-irp)
c>     &      - cspcos_coli((mm22-q12)/(mm32-q13),sq13*rp/sq12
c>     &      ,dreal(mm22-q12-mm32+q13),irp)


#ifdef CHECK
        if(detp.eq.cd0) then
          write(*,100) ' C0m0_coli'
          write(*,100) '     case with zero determinant not implemented'
          write(*,111) ' C0m0_coli: p12 = ',p12,ps12
          write(*,111) ' C0m0_coli: p23 = ',p23,ps23
          write(*,111) ' C0m0_coli: p13 = ',p13,ps13
          write(*,111) ' C0m0_coli: m12 = ',m12,ms12
          write(*,111) ' C0m0_coli: m22 = ',m22,ms22
          write(*,111) ' C0m0_coli: m32 = ',m32,ms32
          stop
        end if
#endif        

        C0m0_coli = C0m0_coli/detp

c        write(*,*) ' C0m0_coli= ',C0m0_coli

#ifdef CHECK
        if(conscheck)then
          det = sqrt((p23-(m2+m3)*(m2+m3))*(p23-(m2-m3)*(m2-m3))) 
          x2 = (p23-mm32-mm22-det)/(2d0*mm32)
          x1 = (p23-mm32-mm22+det)/(2d0*mm32)
          if (abs(x2).lt.abs(x1)) then
            x2 = mm22/(mm32*x1)
          else
            x1 = mm22/(mm32*x2)
          end if          

          if(det.ne.cd0)then
            ix1 = dreal(x1*det)
          else
            ix1=1d0
          endif          

          if (dreal(q12).ge.0d0.and.dreal(q13).ge.0d0) then
            sq13 = sqrt(q13)
            sq12 = sqrt(q12) 
            detp = sqrt((p23-(sq12+sq13)*(sq12+sq13))*
     &          (p23-(sq12-sq13)*(sq12-sq13)))
          else 
            detp = sqrt((q13-p23)**2 + q12*q12 - 2d0*q12*(q13+p23))
          end if
          z2 = (p23-q13-q12-detp)/(2d0*q13)
          z1 = (p23-q13-q12+detp)/(2d0*q13)
          if (abs(z2).lt.abs(z1)) then
            z2 = q12/(q13*z1)
          else
            z1 = q12/(q13*z2)
          end if
          iz1 = dreal(z1/detp)
          
          C0m0_check = cspenc_coli(1d0+1d0/z2,iz1) 
     &        - cspenc_coli(1d0+1d0/z1,-iz1)
     &        + (cln_coli(-z2,iz1)-cln_coli(-z1,-iz1))
     &        *cln_coli((mm32-q13)/mm32,-1d0)
     &        + cspcos_coli(-x1,-1d0/z1,-ix1,iz1)
     &        + cspcos_coli(-x2,-1d0/z1,ix1,iz1)
     &        - cspcos_coli(-x1,-1d0/z2,-ix1,-iz1)
     &        - cspcos_coli(-x2,-1d0/z2,ix1,-iz1)
     &        + cspcos_coli((mm22-q12)/(mm32-q13),-1d0/z2
     &        ,dreal(q13-mm32+mm22-q12),-iz1)
     &        - cspcos_coli((mm22-q12)/(mm32-q13),-1d0/z1
     &        ,dreal(q13-mm32+mm22-q12),iz1)
          C0m0_check = C0m0_check/detp
          
c>      write(*,*) ' C0m0_coli x',x2,x1,z2,z1,detp
c>
c>      write(*,*) ' C0m0_check c',
c>     &        cspenc_coli(1d0+1d0/z2,iz1), 
c>     &       - cspenc_coli(1d0+1d0/z1,-iz1),
c>     &       + (cln_coli(-z2,iz1)-cln_coli(-z1,-iz1))
c>     &           *cln_coli((mm32-q13)/mm32,-1d0),
c>     &       + cspcos_coli(-x1,-1d0/z1,-ix1,iz1),
c>     &       + cspcos_coli(-x2,-1d0/z1,ix1,iz1),
c>     &       - cspcos_coli(-x1,-1d0/z2,-ix1,-iz1),
c>     &       - cspcos_coli(-x2,-1d0/z2,ix1,-iz1),
c>     &       + cspcos_coli((mm22-q12)/(mm32-q13),-1d0/z2
c>     &       ,dreal(q13-mm32+mm22-q12),-iz1),
c>     &       - cspcos_coli((mm22-q12)/(mm32-q13),-1d0/z1
c>     &       ,dreal(q13-mm32+mm22-q12),iz1)

c>      write(*,*) ' C0m0_check s',
c>     &        cspenc_coli(1d0+1d0/z2,iz1) 
c>     &       - cspenc_coli(1d0+1d0/z1,-iz1)
c>     &       + (cln_coli(-z2,iz1)-cln_coli(-z1,-iz1))
c>     &           *cln_coli((mm32-q13)/mm32,-1d0)
c>     &       + cspcos_coli(-x1,-1d0/z1,-ix1,iz1)
c>     &       + cspcos_coli(-x2,-1d0/z1,ix1,iz1)
c>     &       - cspcos_coli(-x1,-1d0/z2,-ix1,-iz1)
c>     &       - cspcos_coli(-x2,-1d0/z2,ix1,-iz1)
c>     &       + cspcos_coli((mm22-q12)/(mm32-q13),-1d0/z2
c>     &       ,dreal(q13-mm32+mm22-q12),-iz1)
c>     &       - cspcos_coli((mm22-q12)/(mm32-q13),-1d0/z1
c>     &       ,dreal(q13-mm32+mm22-q12),iz1),det

c      write(*,*) ' C0m0_check ',C0m0_check

          if(cdabs(C0m0_coli/C0m0_check-1d0).gt.testacc) then
            write(testout,100) ' C0m0_coli inconsistent'
            write(testout,111) ' C0m0_coli: p12 = ',p12
            write(testout,111) ' C0m0_coli: p23 = ',p23
            write(testout,111) ' C0m0_coli: p13 = ',p13
            write(testout,111) ' C0m0_coli: m12 = ',m12
            write(testout,111) ' C0m0_coli: m22 = ',m22
            write(testout,111) ' C0m0_coli: m32 = ',m32
            write(testout,111) ' C0m0_coli  = ',C0m0_coli
            write(testout,111) ' C0m0_check = ',C0m0_check
            write(testout,101) ' diff       = ',
     &          cdabs(C0m0_coli/C0m0_check-1d0)
            
            write(*,100) ' C0m0_coli inconsistent'
            write(*,111) ' C0m0_coli: p12 = ',p12
            write(*,111) ' C0m0_coli: p23 = ',p23
            write(*,111) ' C0m0_coli: p13 = ',p13
            write(*,111) ' C0m0_coli: m12 = ',m12
            write(*,111) ' C0m0_coli: m22 = ',m22
            write(*,111) ' C0m0_coli: m32 = ',m32
            write(*,111) ' C0m0_coli  = ',C0m0_coli
            write(*,111) ' C0m0_check = ',C0m0_check
            write(*,111) ' diff        = ',
     &          cdabs(C0m0_coli/C0m0_check-1d0)
c            stop
          endif
        endif
#endif 

      else if (p23.eq.cd0) then
        if(p12.ne.p13) then
          if (m32.ne.p13.and.p12.ne.cd0) then
            mm32 = m32
            mm22 = m22
            q12 = p12
            q13 = p13
          else
            mm32 = m22
            mm22 = m32
            q12 = p13
            q13 = p12
          end if          
          x1 = -mm22/mm32
          x3 = -(mm22-q12)/(mm32-q13)
          ix1 = dreal(mm32 - mm22)
          ix3 = dreal(mm32-q13 - mm22+q12) 
          
#ifdef ALLCHECK
c ppttbb dim
          if(flag(2))then
            write(*,100)  ' C0m0_coli: p23=0, p12=/=p13 '
            write(*,100) '    case not yet tested in physical process'
            write(*,111) ' C0m0_coli: p12 = ',p12
            write(*,111) ' C0m0_coli: p23 = ',p23
            write(*,111) ' C0m0_coli: p13 = ',p13
            write(*,111) ' C0m0_coli: m12 = ',m12
            write(*,111) ' C0m0_coli: m22 = ',m22
            write(*,111) ' C0m0_coli: m32 = ',m32
            flag(2)=.false.
          endif
#endif

          C0m0_coli = cspenc_coli(1d0+x3,ix3) - cspenc_coli(1+x1,ix1)
          if(q13.ne.cd0) then
            z1 = -q12/q13
            iz1 = dreal(q12 - q13)
            C0m0_coli = C0m0_coli 
     &           + cln_coli(-z1,-iz1)
     &          *cln_coli(mm32/(mm32-q13),dreal(p13))
     &           + cspcos_coli(-x1,-1d0/z1,-ix1,iz1)
     &           - cspcos_coli(-x3,-1d0/z1,-ix3,iz1)
          end if
          C0m0_coli = C0m0_coli/(q13-q12)

        elseif(p12.ne.cd0)then


          if (m32.ne.m22) then

#ifdef CHECK
            if(flag(3))then
              write(*,100) ' C0m0_coli: p23=0, p12=p13 '
              write(*,100) '    case not yet tested at all'
              write(*,111) ' C0m0_coli: p12 = ',p12
              write(*,111) ' C0m0_coli: p23 = ',p23
              write(*,111) ' C0m0_coli: p13 = ',p13
              write(*,111) ' C0m0_coli: m12 = ',m12
              write(*,111) ' C0m0_coli: m22 = ',m22
              write(*,111) ' C0m0_coli: m32 = ',m32
              flag(3)=.false.
            endif
#endif

            C0m0_coli = (1d0-m22/p12)*cln_coli(1d0-m22/p12,dreal(p12))
     &          -(1d0-m32/p12)*cln_coli(1d0-m32/p12,dreal(p12))
     &          +m22/p12*cln_coli(-m22/p12,dreal(p12))
     &          -m32/p12*cln_coli(-m32/p12,dreal(p12))
            C0m0_coli = C0m0_coli/(m32-m22)
          else

#ifdef CHECK
            if(flag(7))then
              write(*,100) ' C0m0_coli: p23=0, p12=p13, m22=m32 '
              write(*,100) '    case not yet tested at all'
              write(*,111) ' C0m0_coli: p12 = ',p12
              write(*,111) ' C0m0_coli: p23 = ',p23
              write(*,111) ' C0m0_coli: p13 = ',p13
              write(*,111) ' C0m0_coli: m12 = ',m12
              write(*,111) ' C0m0_coli: m22 = ',m22
              write(*,111) ' C0m0_coli: m32 = ',m32
              flag(7)=.false.
            endif
#endif

            C0m0_coli = 
     &          (cln_coli(1d0-m32/p12,dreal(p12))
     &          -cln_coli(-m32/p12,dreal(p12)))/p12
          endif           

        elseif(m22.ne.m32)then

#ifdef CHECK
          if(flag(4))then
            write(*,100) ' C0m0_coli: p23=0=p12=p13 '
            write(*,100) '    case not yet tested at all'
            write(*,111) ' C0m0_coli: p12 = ',p12
            write(*,111) ' C0m0_coli: p23 = ',p23
            write(*,111) ' C0m0_coli: p13 = ',p13
            write(*,111) ' C0m0_coli: m12 = ',m12
            write(*,111) ' C0m0_coli: m22 = ',m22
            write(*,111) ' C0m0_coli: m32 = ',m32
            flag(4)=.false.
          endif
#endif

          C0m0_coli = log(m32/m22)/(m22-m32)
        else

#ifdef CHECK
          if(flag(8))then
            write(*,100) ' C0m0_coli: p23=0=p12=p13 m32=m22'
            write(*,100) '    case not yet tested at all'
            write(*,111) ' C0m0_coli: p12 = ',p12
            write(*,111) ' C0m0_coli: p23 = ',p23
            write(*,111) ' C0m0_coli: p13 = ',p13
            write(*,111) ' C0m0_coli: m12 = ',m12
            write(*,111) ' C0m0_coli: m22 = ',m22
            write(*,111) ' C0m0_coli: m32 = ',m32
            flag(8)=.false.
          endif
#endif

          C0m0_coli = -1d0/m32
        end if

      else

        if (p13.ne.cd0) then
          mm32 = m32
          mm22 = m22
c          q12 = p12
          q13 = p13
        else
          mm32 = m22
          mm22 = m32
c          q12 = p13
          q13 = p12
        end if          
        m2 = sqrt(mm22)
        m3 = sqrt(mm32)
        det = sqrt((p23-(m2+m3)*(m2+m3))*(p23-(m2-m3)*(m2-m3))) 
        x2 = (p23-mm32-mm22-det)/(2d0*mm32)
        x1 = (p23-mm32-mm22+det)/(2d0*mm32)
        if (abs(x2).lt.abs(x1)) then
          x2 = mm22/(mm32*x1)
        else
          x1 = mm22/(mm32*x2)
        end if
        if(det.ne.cd0)then
          ix1 = dreal(x1/det)
        else
          ix1 = 1d0
        endif

        if (q13.ne.cd0) then  

#ifdef ALLCHECK 
c ppenj
          if(flag(5))then
            write(*,100) ' C0m0_coli: p12=0, p23,p13 =/=0 '
            write(*,100) '    case not yet tested in physical process'
            write(*,111) ' C0m0_coli: p12 = ',p12
            write(*,111) ' C0m0_coli: p23 = ',p23
            write(*,111) ' C0m0_coli: p13 = ',p13
            write(*,111) ' C0m0_coli: m12 = ',m12
            write(*,111) ' C0m0_coli: m22 = ',m22
            write(*,111) ' C0m0_coli: m32 = ',m32
            flag(5)=.false.
          endif
#endif

          invx3 = -(mm32-q13)/mm22
          ix3 = dreal(mm32-q13-mm22)  
          
          z1 = p23/q13-1d0
          iz1 = dreal(q13-p23)
          
          C0m0_coli = cspcos_coli(-1d0/x1,-z1,ix1,-iz1)
     &         +cspcos_coli(-1d0/x2,-z1,-ix1,-iz1)
     &         -cspenc_coli(1d0+z1,iz1)
     &         -cspcos_coli(-invx3,-z1,ix3,-iz1)
          C0m0_coli = C0m0_coli/(q13-p23)

c>          write(*,*) ' C0m0_coli x',x1,ix1,x2,x3,ix3,z1,iz1
c>
c>          write(*,*) ' C0m0_coli c',
c>     &         cspcos_coli(-1d0/x1,-z1,ix1,-iz1),
c>     &         +cspcos_coli(-1d0/x2,-z1,-ix1,-iz1),
c>     &         -cspenc_coli(1d0+z1,iz1),
c>     &         -cspcos_coli(-invx3,-z1,ix3,-iz1),
c>     &       (q13-p23)
c>
c>          write(*,*) ' C0m0_coli ',C0m0_coli

        else

#ifdef ALLCHECK
c ppenj
          if(flag(6))then
            write(*,100) ' C0m0_coli: p12=0=p13 , p23 =/=0 '
            write(*,100) '    case not yet tested in physical process'
            write(*,111) ' C0m0_coli: p12 = ',p12
            write(*,111) ' C0m0_coli: p23 = ',p23
            write(*,111) ' C0m0_coli: p13 = ',p13
            write(*,111) ' C0m0_coli: m12 = ',m12
            write(*,111) ' C0m0_coli: m22 = ',m22
            write(*,111) ' C0m0_coli: m32 = ',m32
            flag(6)=.false.
          endif
#endif

          C0m0_coli = cln_coli(-x1,-ix1)**2 +  cln_coli(-x2,ix1)**2
     &         - cdlog(mm22/mm32)**2
          C0m0_coli = C0m0_coli/(2d0*p23)       
        end if
      end if
c      write(*,*) ' C0m0_coli ',C0m0_coli
      end
              
              
************************************************************************
      function C02m0_coli(p12,p23,p13,m12,m22,m32)
************************************************************************
*     regular scalar 3-point function with m12 = m32 = 0               *
*     assumes no mass singularities and no ir singularities            *
*                                                                      *
*                            p23                                       *
*                             |                                        *
*                             |                                        *
*                            / \                                       *
*                           /   \                                      *
*                    m22   /     \   m32 = 0                           *
*                         /       \                                    *
*                        /         \                                   *
*              p12  ---------------------  p13                         *
*                          m12 = 0                                     *
*                                                                      *
*----------------------------------------------------------------------*
*     18.05.04 Ansgar Denner        last changed  25.09.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p13,m12,m22,m32
      complex*16 C02m0_coli

      complex*16 sp23,sp12
      complex*16 q23,q12
      real*8     ix1,iz1,ix3
      complex*16 x1,z1,z2,detp,x3
      complex*16 cspenc_coli,cspcos_coli,cln_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps13,ms12,ms22,ms32
      logical flag(0:5)
      data    flag /6*.true./
      save    flag  
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) ' C02m0_coli in ',p12,p23,p13,m12,m22,m32
 
#ifdef CHECK
 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        ms12 = elimminf2_coli(m12)
        ms22 = elimminf2_coli(m22)
        ms32 = elimminf2_coli(m32)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps13 = elimminf2_coli(p13)

        if(m12.ne.cd0.or.m32.ne.cd0
     &      .or.m22.eq.cd0.and.p23.eq.cd0
     &      .or.m22.eq.cd0.and.p12.eq.cd0
     &      .or.p13.eq.cd0) then
          write(*,100) ' C02m0_coli called improperly'
          write(*,100) '     or case not implemented '
          write(*,111) ' C02m0_coli: p12 = ',p12
          write(*,111) ' C02m0_coli: p23 = ',p23
          write(*,111) ' C02m0_coli: p13 = ',p13
          write(*,111) ' C02m0_coli: m12 = ',m12
          write(*,111) ' C02m0_coli: m22 = ',m22
          write(*,111) ' C02m0_coli: m32 = ',m32
          stop
        end if
      end if
#endif
        
      if (p12.ne.cd0.and.p13.ne.cd0.and.p23.ne.cd0) then
        if (p12.ne.m22) then
          q12 = p12
          q23 = p23
        else
          q12 = p23
          q23 = p12
        end if          
        
#ifdef ALLCHECK
c ppttbb
        if(flag(1))then
          write(*,100) ' C02m0_coli: p12,p23,p13=/=0'
          write(*,100) '    case not yet tested in physical process'
          write(*,111) ' C02m0_coli: p12 = ',p12
          write(*,111) ' C02m0_coli: p23 = ',p23
          write(*,111) ' C02m0_coli: p13 = ',p13
          write(*,111) ' C02m0_coli: m12 = ',m12
          write(*,111) ' C02m0_coli: m22 = ',m22
          write(*,111) ' C02m0_coli: m32 = ',m32
          flag(1)=.false.
        endif
#endif

        if (dreal(q12).ge.0d0.and.dreal(p13).ge.0d0) then
          sp12 = sqrt(q12)
          sp23 = sqrt(q23) 
          detp = sqrt((p13-(sp12+sp23)*(sp12+sp23))*
     &         (p13-(sp12-sp23)*(sp12-sp23)))
        else 
          detp = sqrt((p13-q23)**2 + q12*q12 - 2d0*q12*(p13+q23))
        end if
        z2 = (q12-p13-q23-detp)/(2d0*p13)
        z1 = (q12-p13-q23+detp)/(2d0*p13)
        if (abs(z2).lt.abs(z1)) then
          z2 = q23/(p13*z1)
        else
          z1 = q23/(p13*z2)
        end if
        iz1 = dreal(z1/detp)
        
        x1  = (m22-q23)/p13
        ix1 = -p13 - m22 + q23
        x3  = (q12-m22)/m22
        ix3 = 1d0

c>        C02m0_coli = cspenc_coli(1d0+z2,-iz1) 
c>     &       - cspenc_coli(1d0+z1,iz1)
c>     &       + (cln_coli(-z2,iz1)-cln_coli(-z1,-iz1))
c>     &           *cln_coli((q12-m22)/p13,p13-q12+m22)
c>     &       + cspcos_coli(-x1,-1d0/z1,-ix1,iz1)
c>     &       - cspcos_coli(-x1,-1d0/z2,-ix1,-iz1)
c>     &       + cspcos_coli(-1d0/x3,-1d0/z2,ix3,-iz1)
c>     &       - cspcos_coli(-1d0/x3,-1d0/z1,ix3,iz1)
c>        C02m0_coli = -C02m0_coli/detp
c>
c>        write(*,*) ' C02m0_coli 1',C02m0_coli
c>
        C02m0_coli = cspenc_coli(1d0+z2,-iz1) 
     &       - cspenc_coli(1d0+z1,iz1)
     &       - ((cln_coli(-z2,iz1))**2 - (cln_coli(-z1,-iz1))**2)/2d0 
     &       + (cln_coli(-z2,iz1)-cln_coli(-z1,-iz1))
     &           *cln_coli(-m22/p13,p13)
     &       + cspcos_coli(-x1,-1d0/z1,-ix1,iz1)
     &       - cspcos_coli(-x1,-1d0/z2,-ix1,-iz1)
     &       - cspcos_coli(-x3,-z2,-ix3,iz1)
     &       + cspcos_coli(-x3,-z1,-ix3,-iz1)
        C02m0_coli = -C02m0_coli/detp

c        write(*,*) ' C02m0_coli 2',C02m0_coli
        
      else if (p13.ne.cd0) then
        if (p23.ne.cd0) then
c          q12 = p12
          q23 = p23
        else
c          q12 = p23
          q23 = p12
        end if          

        if (abs(q23).lt.abs(q23-m22)) then

          z1  = -q23/p13
          iz1 = dreal(q23-p13)
          x1  = p13/(m22-q23)
          ix1 = dreal(p13 + m22 - q23)
          
          C02m0_coli = -cspenc_coli(1d0+x1,ix1)
          if (q23.ne.cd0) then

#ifdef ALLCHECK
c ppenj
            if(flag(2))then
              write(*,100) ' C02m0_coli: p12=0 a'
              write(*,100) '    case not yet tested in physical process'
              write(*,111) ' C02m0_coli: p12 = ',p12
              write(*,111) ' C02m0_coli: p23 = ',p23
              write(*,111) ' C02m0_coli: p13 = ',p13
              write(*,111) ' C02m0_coli: m12 = ',m12
              write(*,111) ' C02m0_coli: m22 = ',m22
              write(*,111) ' C02m0_coli: m32 = ',m32
              flag(2)=.false.
            endif
#endif

            C02m0_coli = C02m0_coli  
     &           + cln_coli(-z1,-iz1)
     &          *cln_coli(m22/(m22-q23),dreal(q23))
     &           + cspcos_coli(-x1,-z1,-ix1,-iz1)
          else

#ifdef ALLCHECK
c ppenj
            if(flag(3))then
              write(*,100) ' C02m0_coli: p12=0=p23'
              write(*,100) '    case not yet tested in physical process'
              write(*,111) ' C02m0_coli: p12 = ',p12
              write(*,111) ' C02m0_coli: p23 = ',p23
              write(*,111) ' C02m0_coli: p13 = ',p13
              write(*,111) ' C02m0_coli: m12 = ',m12
              write(*,111) ' C02m0_coli: m22 = ',m22
              write(*,111) ' C02m0_coli: m32 = ',m32
              flag(3)=.false.
            endif
#endif

            C02m0_coli = C02m0_coli + pi2_6
          end if
          C02m0_coli = C02m0_coli/(p13-q23)
        else

          z1  = -q23/p13
          iz1 = dreal(q23-p13)
          x1  = (m22-q23)/p13
          ix1 = dreal(q23-p13 - m22)
          
          C02m0_coli = -cln_coli(-z1,-iz1)**2/2d0 
     &         + cln_coli(-z1,-iz1)*cln_coli(-m22/p13,dreal(p13+m22))
          if (q23.ne.m22) then

#ifdef ALLCHECK
c ppenj
            if(flag(4))then
              write(*,100) ' C02m0_coli: p12=0 b'
              write(*,100) '    case not yet tested in physical process'
              write(*,111) ' C02m0_coli: p12 = ',p12
              write(*,111) ' C02m0_coli: p23 = ',p23
              write(*,111) ' C02m0_coli: p13 = ',p13
              write(*,111) ' C02m0_coli: m12 = ',m12
              write(*,111) ' C02m0_coli: m22 = ',m22
              write(*,111) ' C02m0_coli: m32 = ',m32
              flag(4)=.false.
            endif
#endif

            C02m0_coli = C02m0_coli  
     &           + cspenc_coli(1d0+x1,ix1)
     &           - cspcos_coli(-x1,-1d0/z1,-ix1,iz1)

          else

#ifdef CHECK
            if(flag(5))then
              write(*,100) ' C02m0_coli: p12=0, p23=m22'
              write(*,100) '    case not yet tested in physical process'
              write(*,111) ' C02m0_coli: p12 = ',p12
              write(*,111) ' C02m0_coli: p23 = ',p23
              write(*,111) ' C02m0_coli: p13 = ',p13
              write(*,111) ' C02m0_coli: m12 = ',m12
              write(*,111) ' C02m0_coli: m22 = ',m22
              write(*,111) ' C02m0_coli: m32 = ',m32
              flag(5)=.false.
            endif
#endif

            C02m0_coli = C02m0_coli 
          end if
          C02m0_coli = C02m0_coli/(p13-q23)
        end if 
#ifdef CHECK
      else
        write(*,100) ' C0m00 called with improper arguments'
        write(*,100) '     or case not implemented    '
        write(*,111) ' C02m0_coli: p12 = ',p12
        write(*,111) ' C02m0_coli: p23 = ',p23
        write(*,111) ' C02m0_coli: p13 = ',p13
        write(*,111) ' C02m0_coli: m12 = ',m12
        write(*,111) ' C02m0_coli: m22 = ',m22
        write(*,111) ' C02m0_coli: m32 = ',m32
        stop
#endif
      end if

c      write(*,*) ' C02m0_coli out ', C02m0_coli     

      end
              
************************************************************************
      function C03m0_coli(p12,p23,p13,m12,m22,m32)
************************************************************************
*     regular scalar 3-point function with m12 = m22 = m32 = 0         *
*     assumes no mass singularities and no ir singularities            *
*                                                                      *
*                            p23                                       *
*                             |                                        *
*                             |                                        *
*                            / \                                       *
*                           /   \                                      *
*                m22 = 0   /     \   m32 = 0                           *
*                         /       \                                    *
*                        /         \                                   *
*              p12  ---------------------  p13                         *
*                          m12 = 0                                     *
*                                                                      *
*----------------------------------------------------------------------*
*     18.05.04 Ansgar Denner        last changed  28.08.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p13,m12,m22,m32
      complex*16 C03m0_coli

      complex*16 sp12,sp13
      real*8     ix1,iz1
      complex*16 x1,z1,z2,detp
      complex*16 cspenc_coli,cspcos_coli,cln_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps13,ms12,ms22,ms32
      logical flag(0:1)
      data    flag /2*.true./
      save    flag  
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(30,*) ' C03m0_coli in ',p12,p23,p13,m12,m22,m32,ltest
 
#ifdef CHECK
 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        ms12 = elimminf2_coli(m12)
        ms22 = elimminf2_coli(m22)
        ms32 = elimminf2_coli(m32)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps13 = elimminf2_coli(p13)

        if(m12.ne.cd0.or.m22.ne.cd0.or.m32.ne.cd0
     &        .or.ps12.eq.cd0.or.ps13.eq.cd0.or.ps23.eq.cd0
     &        ) then
          write(*,100) ' C03m0_coli called improperly'
          write(*,100) '     or case not implemented '
          write(*,111) ' C03m0_coli: p12 = ',p12
          write(*,111) ' C03m0_coli: p23 = ',p23
          write(*,111) ' C03m0_coli: p13 = ',p13
          write(*,111) ' C03m0_coli: m12 = ',m12
          write(*,111) ' C03m0_coli: m22 = ',m22
          write(*,111) ' C03m0_coli: m32 = ',m32
          stop
        end if
      end if
#endif
        
#ifdef ALLCHECK
c ppttbb
        if(flag(1))then
          write(*,100) ' C03m0_coli: '
          write(*,100) '    case not yet tested in physical process'
          write(*,111) ' C03m0_coli: p12 = ',p12
          write(*,111) ' C03m0_coli: p23 = ',p23
          write(*,111) ' C03m0_coli: p13 = ',p13
          write(*,111) ' C03m0_coli: m12 = ',m12
          write(*,111) ' C03m0_coli: m22 = ',m22
          write(*,111) ' C03m0_coli: m32 = ',m32
          flag(1)=.false.
        endif
#endif

      if (dreal(p12).ge.0d0.and.dreal(p13).ge.0d0) then
        sp13 = sqrt(p13)
        sp12 = sqrt(p12) 
        detp = sqrt((p23-(sp12+sp13)*(sp12+sp13))*
     &       (p23-(sp12-sp13)*(sp12-sp13)))
      else 
        detp = sqrt((p13-p23)**2 + p12*p12 - 2d0*p12*(p13+p23))
      end if
      z2 = (p12-p13-p23-detp)/(2d0*p13)
      z1 = (p12-p13-p23+detp)/(2d0*p13)
      if (abs(z2).lt.abs(z1)) then
        z2 = p23/(p13*z1)
      else
        z1 = p23/(p13*z2)
      end if
      iz1 = dreal(z1/detp)
        
      x1  = -p23/p13
      ix1 = p23 - p13

      C03m0_coli = cspenc_coli(1d0+z2,-iz1) 
     &       - cspenc_coli(1d0+z1,iz1)
     &       + (cln_coli(-z2,iz1)-cln_coli(-z1,-iz1))
     &           *cln_coli(p12/p13,p13-p12)
     &       + cspcos_coli(-x1,-1d0/z1,-ix1,iz1)
     &       - cspcos_coli(-x1,-1d0/z2,-ix1,-iz1)

      C03m0_coli = -C03m0_coli/detp
        
      end
************************************************************************
      function C0ms0ir1_coli(p12,p23,p13,m12,m22,m32)
************************************************************************
*     ir-divergent scalar 3-point function                             *
*     assumes    m12 = 0, p12 = m22, p13 = m32                         * 
*     and no mass singularities                                        *
*                                                                      *
*                            p23                                       *
*                             |                                        *
*                             |                                        *
*                            / \                                       *
*                           /   \                                      *
*                    m22   /     \   m32                               *
*                         /       \                                    *
*                        /         \                                   *
*       m22 = p12  ---------------------  p13 = m32                    *
*                            la2                                       *
*                                                                      *
*----------------------------------------------------------------------*
*     15.03.04 Ansgar Denner        last changed  28.08.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p13
      complex*16 m12,m22,m32
      complex*16 C0ms0ir1_coli

      complex*16 a,h1,h2,h3
      complex*16 cln_coli,cspenc_coli,cspcon_coli

#ifdef CHECK
      complex*16 C0ms0ir1_check
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps13
      complex*16 ms12,ms22,ms32
      complex*16 m2,m3,k23,r23
      real*8     ir23
      logical flag(0:4)
      data    flag /5*.true./
      save    flag  
#endif

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

c      write(*,*) ' C0ms0ir1_coli in ',p12,p23,p13,m12,m22,m32

#ifdef CHECK
 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        ms12 = elimminf2_coli(m12)
        ms22 = elimminf2_coli(m22)
        ms32 = elimminf2_coli(m32)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps13 = elimminf2_coli(p13)

        if(ms12.ne.cd0.or.p12.ne.m22.or.p13.ne.m32
     &      .or.ms22.eq.cd0.or.ms32.eq.cd0
     &      .or.ps23.eq.cd0.and.ms32.eq.cd0.and.ms22.eq.cd0
     &      .or.dimag(p12).ne.0d0.or.dimag(p13).ne.0d0
     &      .or.dimag(p23).ne.0d0) then
          write(*,100) ' C0ms0ir1_coli called improperly'
          write(*,100) '     or case not implemented    '
          write(*,111) ' C0ms0ir1_coli: p12 = ',p12
          write(*,111) ' C0ms0ir1_coli: p23 = ',p23
          write(*,111) ' C0ms0ir1_coli: p13 = ',p13
          write(*,111) ' C0ms0ir1_coli: m12 = ',m12
          write(*,111) ' C0ms0ir1_coli: m22 = ',m22
          write(*,111) ' C0ms0ir1_coli: m32 = ',m32
          stop
        end if
      end if
#endif

c W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349  (iv)
c according to Dittmaier Nucl. Phys. B675 (2003) 447   (B.5,B.6)          

        if (p23.ne.cd0) then
          a=sqrt(p12*p12+p23*p23+p13*p13-2d0*(p12*p23+p23*p13+p13*p12))
          if (dreal(p23-p12-p13).le.0d0) then
            h1=(a-p23+p12+p13)/2d0
          else
            h1=-4d0*p12*p13/(a+p23-p12-p13)/2d0
          end if
          if (dreal(p23-p12+p13).le.0d0) then
            h2=(a-p23+p12-p13)/2d0
          else
            h2=-4d0*p23*p13/(a+p23-p12+p13)/2d0
          end if
          if (dreal(p23+p12-p13).le.0d0) then
            h3=(a-p23-p12+p13)/2d0
          else
            h3=-4d0*p12*p23/(a+p23+p12-p13)/2d0
          end if

          C0ms0ir1_coli=cspenc_coli(h2/a,-1d0)
     *         +cspenc_coli(h3/a,-1d0)-pi2_6
     *         -.5d0*cln_coli(-h2/p23,-1d0)**2
     *         -.5d0*cln_coli(-h3/p23,-1d0)**2
     *         +.25d0*cln_coli(-p12/p23,-1d0)**2
     *         +.25d0*cln_coli(-p13/p23,-1d0)**2
     *         -cln_coli(-a/p23,-1d0)
     *         *(cln_coli(-h1/p23,-1d0)-cln_coli(-a/p23,-1d0))
          if(m12.eq.cd0)then

#ifdef ALLCHECK
c ppttbb
            if(flag(1))then
              write(*,100) ' C0ms0ir1_coli: (B.6)'
              write(*,100) '    case not yet tested in physical process'
              write(*,111) ' C0ms0ir1_coli: p12 = ',p12
              write(*,111) ' C0ms0ir1_coli: p23 = ',p23
              write(*,111) ' C0ms0ir1_coli: p13 = ',p13
              write(*,111) ' C0ms0ir1_coli: m12 = ',m12
              write(*,111) ' C0ms0ir1_coli: m22 = ',m22
              write(*,111) ' C0ms0ir1_coli: m32 = ',m32
              flag(1)=.false.
            endif
#endif

            C0ms0ir1_coli=C0ms0ir1_coli
     *          +cln_coli(-muir2/p23,-1d0)
     *          *cln_coli(h1/sqrt(p12*p13),1d0)
#ifdef SING
     *          +delta1ir
     *          *cln_coli(h1/sqrt(p12*p13),1d0)
#endif
          else

#ifdef CHECK
            if(flag(2))then
              write(*,100) ' C0ms0ir1_coli: (B.5)',
     &          '    case not yet tested in physical process'
              write(*,111) ' C0ms0ir1_coli: p12 = ',p12
              write(*,111) ' C0ms0ir1_coli: p23 = ',p23
              write(*,111) ' C0ms0ir1_coli: p13 = ',p13
              write(*,111) ' C0ms0ir1_coli: m12 = ',m12
              write(*,111) ' C0ms0ir1_coli: m22 = ',m22
              write(*,111) ' C0ms0ir1_coli: m32 = ',m32
              flag(2)=.false.
            endif
#endif

            C0ms0ir1_coli=C0ms0ir1_coli
     *          +cln_coli(-m12/p23,-1d0)            
     *          *cln_coli(h1/sqrt(p12*p13),1d0)
          endif
          C0ms0ir1_coli=C0ms0ir1_coli/a
        else
           if(m12.eq.cd0)then

#ifdef CHECK
            if(flag(3))then
              write(*,100) ' C0ms0ir1_coli: (B.6) for s=0',
     &          '    case not yet tested in physical process'
              write(*,111) ' C0ms0ir1_coli: p12 = ',p12
              write(*,111) ' C0ms0ir1_coli: p23 = ',p23
              write(*,111) ' C0ms0ir1_coli: p13 = ',p13
              write(*,111) ' C0ms0ir1_coli: m12 = ',m12
              write(*,111) ' C0ms0ir1_coli: m22 = ',m22
              write(*,111) ' C0ms0ir1_coli: m32 = ',m32
              flag(3)=.false.
            endif
#endif

             C0ms0ir1_coli=-.25d0*log(p13*p12/muir2**2)*log(p12/p13)
#ifdef SING
     *           +.5d0*delta1ir*log(p12/p13)
     *           *cln_coli(h1/sqrt(p12*p13),1d0)
#endif
           else

#ifdef CHECK
            if(flag(4))then
              write(*,100) ' C0ms0ir1_coli: (B.5) for s=0',
     &          '    case not yet tested in physical process'
              write(*,111) ' C0ms0ir1_coli: p12 = ',p12
              write(*,111) ' C0ms0ir1_coli: p23 = ',p23
              write(*,111) ' C0ms0ir1_coli: p13 = ',p13
              write(*,111) ' C0ms0ir1_coli: m12 = ',m12
              write(*,111) ' C0ms0ir1_coli: m22 = ',m22
              write(*,111) ' C0ms0ir1_coli: m32 = ',m32
              flag(4)=.false.
            endif
#endif
            
            C0ms0ir1_coli=-.25d0*log(p13*p12/m12**2)*log(p12/p13)
          endif
          C0ms0ir1_coli=C0ms0ir1_coli/(p12-p13)
        end if

#ifdef CHECK
        if(conscheck) then
          
c W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349  (iv)
          m2 = sqrt(m22)
          m3 = sqrt(m32)
          k23=(m22+m32-p23)/(m2*m3)
          if (k23.ne.cd0) then
c     |r23| > 0 yields wrong results !!!!!!!  
c     analytic continuation requires |r23| < 1    ????
            r23 = 2d0/(k23*(1d0+sqrt(dcmplx(1d0-4d0/k23**2))))
          else
            r23 = dcmplx(0d0,1d0)
          end if
          if(dreal(k23).lt.-2d0) then
            ir23 = sign(1d1,1d0-abs(r23))
          else
            ir23 = 0d0
          end if
            
          if(p23.ne.(m2-m3)**2) then
            C0ms0ir1_check = 
     &          cln_coli(r23,ir23)*(
     &          2d0*cln_coli(1d0-r23*r23,-dreal(r23)*ir23)
     &          -.5d0*cln_coli(r23,ir23)) 
     &          -pi2_6
     &          + cspenc_coli(r23*r23,dreal(r23)*ir23)
     &          + .5d0*(log(m2/m3))**2 
     &          + cspcon_coli(r23,m2/m3,r23*m2/m3,1d0-r23*m2/m3,
     &          ir23,0d0)
     &          + cspcon_coli(r23,m3/m2,r23*m3/m2,1d0-r23*m3/m2,
     &          ir23,0d0)

            if(m12.eq.cd0)then
              C0ms0ir1_check = C0ms0ir1_check
     &            +cln_coli(r23,ir23)*log((m2*m3)/muir2) 
#ifdef SING
     &            - delta1ir* cln_coli(r23,ir23)
#endif
            else
              C0ms0ir1_check = C0ms0ir1_check
     &            +cln_coli(r23,ir23)*log((m2*m3)/m12) 
            endif            
            C0ms0ir1_check = C0ms0ir1_check
     &          *r23/(m2*m3*(1-r23*r23))
          elseif(m32.ne.m22) then
            if(m12.eq.cd0)then
              C0ms0ir1_check = 
     &            log(muir2/(m2*m3))-2d0-(m3+m2)/(m3-m2)*log(m2/m3)
#ifdef SING
     &            + delta1ir
#endif
            else
              C0ms0ir1_check = 
     &            log(m12/(m2*m3))-2d0-(m3+m2)/(m3-m2)*log(m2/m3)
            endif
            C0ms0ir1_check = C0ms0ir1_check/(2d0*m2*m3)
          else
            if(m12.eq.cd0)then
              C0ms0ir1_check = log(muir2/(m2*m3))
#ifdef SING
     &            + delta1ir
#endif
            else
              C0ms0ir1_check = log(m12/(m2*m3))
            endif
            C0ms0ir1_check = C0ms0ir1_check/(2d0*m2*m3)
          endif
            
          if(cdabs(C0ms0ir1_coli/C0ms0ir1_check-1d0).gt.testacc) then
            write(testout,100) ' C0ms0ir1_coli inconsistent'
            write(testout,111) ' C0ms0ir1_coli: p12 = ',p12
            write(testout,111) ' C0ms0ir1_coli: p23 = ',p23
            write(testout,111) ' C0ms0ir1_coli: p13 = ',p13
            write(testout,111) ' C0ms0ir1_coli: m12 = ',m12
            write(testout,111) ' C0ms0ir1_coli: m22 = ',m22
            write(testout,111) ' C0ms0ir1_coli: m32 = ',m32
            write(testout,111) ' C0ms0ir1_coli  = ',C0ms0ir1_coli
            write(testout,111) ' C0ms0ir1_check = ',C0ms0ir1_check
            write(testout,111) ' diff        = ',
     &          cdabs(C0ms0ir1_coli/C0ms0ir1_check-1d0)
            
            write(*,100) ' C0ms0ir1_coli inconsistent'
            write(*,111) ' C0ms0ir1_coli: p12 = ',p12
            write(*,111) ' C0ms0ir1_coli: p23 = ',p23
            write(*,111) ' C0ms0ir1_coli: p13 = ',p13
            write(*,111) ' C0ms0ir1_coli: m12 = ',m12
            write(*,111) ' C0ms0ir1_coli: m22 = ',m22
            write(*,111) ' C0ms0ir1_coli: m32 = ',m32
            write(*,111) ' C0ms0ir1_coli  = ',C0ms0ir1_coli
            write(*,111) ' C0ms0ir1_check = ',C0ms0ir1_check
            write(*,111) ' diff        = ',
     &          cdabs(C0ms0ir1_coli/C0ms0ir1_check-1d0)
            stop
          endif
          
        endif
#endif

        end

************************************************************************
      function C0ms1ir0_coli(p12,p23,p13,m12,m22,m32)
************************************************************************
*     single mass-singular 3-point function                            *
*                                                                      *
*     assumes single mass singularity                                  *
*     and no soft singularity with                                     *
*     p12, m12, m22 small                                              *
*     p23, p13, m32  finite                                            *
*                                                                      *
*                            p23                                       *
*                             |                                        *
*                             |                                        *
*                            / \                                       *
*                           /   \                                      *
*              m22 small   /     \   m32                               *
*                         /       \                                    *
*                        /         \                                   *
*             p12  ---------------------  p13                          *
*            small        m12 small                                    *
*                                                                      *
*----------------------------------------------------------------------*
*     20.08.08 Ansgar Denner        last changed  30.09.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p13
      complex*16 m12,m22,m32
      complex*16 C0ms1ir0_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps13
      complex*16 ms12,ms22,ms32
      complex*16 C0ms1ir0_check
      complex*16 x1,x2,z1,z2
      real*8     ix1,iz1
      logical flag(0:6)
      data    flag /7*.true./
      save    flag  
#endif

      complex*16 mm22,q13,q23
      complex*16 cln_coli,cspenc_coli,cspcon_coli
c      integer    nsoft
c      integer    i,j,k
c      logical    soft(3),onsh(3,3)

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

#ifdef CHECK
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      ms12 = elimminf2_coli(m12)
      ms22 = elimminf2_coli(m22)
      ms32 = elimminf2_coli(m32)
      ps12 = elimminf2_coli(p12)
      ps23 = elimminf2_coli(p23)
      ps13 = elimminf2_coli(p13)
      if(argcheck)then
        if(ms12.ne.cd0.or.ms22.ne.cd0.or.ps12.ne.cd0.or.
     &      ps13.eq.ms32.and.ps23.eq.ms32.or.
     &      ps13.eq.cd0.and.ms32.eq.cd0.or.
     &      ps23.eq.cd0.and.ms32.eq.cd0.or.
     &      ps23.eq.ps13
     &      .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &      .or.dimag(p13).ne.0d0) then
          write(*,100) ' C0ms1ir0_coli called improperly'
          write(*,100) '     or case not implemented    '
          write(*,111) ' C0ms1ir0_coli: p12 = ',p12,ps12
          write(*,111) ' C0ms1ir0_coli: p23 = ',p23,ps23
          write(*,111) ' C0ms1ir0_coli: p13 = ',p13,ps13
          write(*,111) ' C0ms1ir0_coli: m12 = ',m12,ms12
          write(*,111) ' C0ms1ir0_coli: m22 = ',m22,ms22
          write(*,111) ' C0ms1ir0_coli: m32 = ',m32,ms32
          stop
        end if
      end if
#endif

      if(m12.eq.cd0.and.m22.eq.cd0)then
* no soft singularity and 2 zero masses
*
*                   | p23             
*                   |   
*                  / \         
*              0  /2 3\ m32       
*                /  1  \         
*          0 ----------------  p13
*                   0             
*

#ifdef CHECK
        if(m12.ne.cd0.or.m22.ne.cd0.or.p12.ne.cd0
     &      .or.p13.eq.m32.or.p23.eq.m32)
     &      then
          write(*,100) ' C0ms1ir0_coli: case with 2 zero masses',
     &        '    wrong arguments'
          write(*,111) ' C0ms1ir0_coli: p12 = ',p12
          write(*,111) ' C0ms1ir0_coli: p23 = ',p23
          write(*,111) ' C0ms1ir0_coli: p13 = ',p13
          write(*,111) ' C0ms1ir0_coli: m12 = ',m12
          write(*,111) ' C0ms1ir0_coli: m22 = ',m22
          write(*,111) ' C0ms1ir0_coli: m32 = ',m32
          stop
        endif
#endif
          
c Dittmaier Nucl. Phys. B675 (2003) 447   (B.4)          
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002    triangle 3

        if(abs(p23).gt.abs(m32).and.abs(p13).gt.abs(m32)) then

#ifdef ALLCHECK
c ppenj
          if(flag(1))then
            write(*,100) ' C0ms1ir0_coli: (B.4) small m32',
     &          '    case with 2 zero masses',
     &          '    not yet tested in physical process'
            write(*,111) ' C0ms1ir0_coli: p12 = ',p12
            write(*,111) ' C0ms1ir0_coli: p23 = ',p23
            write(*,111) ' C0ms1ir0_coli: p13 = ',p13
            write(*,111) ' C0ms1ir0_coli: m12 = ',m12
            write(*,111) ' C0ms1ir0_coli: m22 = ',m22
            write(*,111) ' C0ms1ir0_coli: m32 = ',m32
            flag(1)=.false.
          endif
#endif
          
          C0ms1ir0_coli =
#ifdef SING
     &        delta1ir*(cln_coli((m32-p23),-1d0)
     &        -cln_coli((m32-p13),-1d0))
#endif
     &        +.5d0*(cln_coli((-p13)/muir2,-1d0))**2
     &        -.5d0*(cln_coli((-p23)/muir2,-1d0))**2
     &        -cspenc_coli(m32/p13,-1d0)
     &        +cspenc_coli(m32/p23,-1d0)
          if(m32.ne.cd0)then

#ifdef ALLCHECK
c ppttbb dim
            if(flag(6))then
              write(*,100) ' C0ms1ir0_coli: (B.4) small finite m32',
     &            '    case with 2 zero masses',
     &            '    not yet tested in physical process'
              write(*,111) ' C0ms1ir0_coli: p12 = ',p12
              write(*,111) ' C0ms1ir0_coli: p23 = ',p23
              write(*,111) ' C0ms1ir0_coli: p13 = ',p13
              write(*,111) ' C0ms1ir0_coli: m12 = ',m12
              write(*,111) ' C0ms1ir0_coli: m22 = ',m22
              write(*,111) ' C0ms1ir0_coli: m32 = ',m32
              flag(6)=.false.
            endif
#endif
            
            C0ms1ir0_coli =  C0ms1ir0_coli
     &          +cln_coli((m32-p13)/(-p13),1d0)*
     &          (cln_coli((m32-p13)/m32,-1d0)+cln_coli(-p13/muir2,-1d0))
     &          -cln_coli((m32-p23)/(-p23),1d0)*
     &          (cln_coli((m32-p23)/m32,-1d0)+cln_coli(-p23/muir2,-1d0))

         endif
        else

#ifdef CHECK
          if(flag(5))then
            write(*,100) ' C0ms1ir0_coli: (B.4) large m32',
     &          '    case with 2 zero masses',
     &          '    not yet tested in physical process'
            write(*,111) ' C0ms1ir0_coli: p12 = ',p12
            write(*,111) ' C0ms1ir0_coli: p23 = ',p23
            write(*,111) ' C0ms1ir0_coli: p13 = ',p13
            write(*,111) ' C0ms1ir0_coli: m12 = ',m12
            write(*,111) ' C0ms1ir0_coli: m22 = ',m22
            write(*,111) ' C0ms1ir0_coli: m32 = ',m32
            flag(5)=.false.
          endif
#endif
          
          C0ms1ir0_coli =
#ifdef SING
     &        delta1ir*(cln_coli((m32-p23)/m32,-1d0)
     &        -cln_coli((m32-p13)/m32,-1d0))
#endif
     &        -(cln_coli((m32-p23)/m32,-1d0))**2
     &        +(cln_coli((m32-p13)/m32,-1d0))**2
     &        -cln_coli(m32/muir2,-1d0)*(cln_coli((m32-p23)/m32,-1d0)
     &        -cln_coli((m32-p13)/m32,-1d0))
     &        +cspenc_coli(p13/m32,1d0)
     &        -cspenc_coli(p23/m32,1d0)
        endif 
        
        C0ms1ir0_coli = C0ms1ir0_coli/(p13-p23)
        
      elseif(m12.eq.cd0.or.m22.eq.cd0)then
* no soft singularity and 1 zero masses
*
*                   | p23             
*                   |   
*                  / \         
*       m22 small /2 3\ m32       
*                /  1  \         
*    p12=m22 ----------------  p13
*    small          0             
*

#ifdef CHECK
        if(m12.ne.cd0.and.m12.ne.p12.or.
     &      m22.ne.cd0.and.m22.ne.p12.or.
     &      m22*m12.ne.cd0.or.
     &      m12.eq.cd0.and.p13.eq.m32.or.
     &      m22.eq.cd0.and.p23.eq.m32)
     &      then
          write(*,100) ' C0ms1ir0_coli: case with 1 zero mass',
     &        ': wrong arguments'
          write(*,111) ' C0ms1ir0_coli: p12 = ',p12
          write(*,111) ' C0ms1ir0_coli: p23 = ',p23
          write(*,111) ' C0ms1ir0_coli: p13 = ',p13
          write(*,111) ' C0ms1ir0_coli: m12 = ',m12
          write(*,111) ' C0ms1ir0_coli: m22 = ',m22
          write(*,111) ' C0ms1ir0_coli: m32 = ',m32
          stop
        endif
#endif
        

        if(m12.eq.cd0)then        
          mm22=m22*coliminfscale2
          q23=p23   
          q13=p13
        else   
          mm22=m12*coliminfscale2
          q23=p13   
          q13=p23
        endif

c     Dittmaier Nucl. Phys. B675 (2003) 447   (B.2)          
        
        if (m32.ne.cd0) then

#ifdef ALLCHECK
c ppenj
          if(flag(2))then
            write(*,100) ' C0ms1ir0_coli: (B.2) with m32=/=0',
     &          '    case with 1 zero masses',
     &          '    not yet tested in physical process'
            write(*,111) ' C0ms1ir0_coli: p12 = ',p12
            write(*,111) ' C0ms1ir0_coli: p23 = ',p23
            write(*,111) ' C0ms1ir0_coli: p13 = ',p13
            write(*,111) ' C0ms1ir0_coli: m12 = ',m12
            write(*,111) ' C0ms1ir0_coli: m22 = ',m22
            write(*,111) ' C0ms1ir0_coli: m32 = ',m32
            flag(2)=.false.
          endif          
#endif

c>          C0ms1ir0_coli = cspcon_coli(mm22/(m32-q23),(q13-q23)/mm22,
c>     &        (q13-q23)/(m32-q23),(m32-q13)/(m32-q23),1d0,
c>     &        -dreal(q13-q23))
c>     &        + cspcon_coli(m32/(m32-q23),(q13-q23)/q13,
c>     &        (q13-q23)*m32/((m32-q23)*q13),
c>     &        (m32-q13)*q23/((m32-q23)*q13),
c>     &        1d0,-dreal(q13-q23))
c>     &        - cspenc_coli(q23/q13,dreal(q13-q23)) - pi2_6
c>     &        + (cln_coli(mm22/(q13-q23),dreal(q13-q23)) 
c>     &        - cln_coli((q13-q23)/q13,-dreal(q13-q23)))
c>     &        * cln_coli((m32-q23)/(m32-q13),dreal(q13-q23)) 


          C0ms1ir0_coli = cspenc_coli(q23/m32,1d0)
     &        - cspenc_coli(q13/m32,1d0)
     &        - 2d0*cspenc_coli((q23-q13)/(m32-q13),q23-q13)
     &        + cln_coli((m32-q23)/mm22,-1d0)
     &         *cln_coli((m32-q23)/m32,-1d0)
     &        - cln_coli((m32-q13)/mm22,-1d0)
     &         *cln_coli((m32-q13)/m32,-1d0)
          C0ms1ir0_coli = C0ms1ir0_coli/(q23-q13)

        else

#ifdef ALLCHECK
c ppenj
          if(flag(3))then
            write(*,100) ' C0ms1ir0_coli: (B.2) with m32=0',
     &          '    case with 1 zero masses',
     &          '    not yet tested in physical process'
            write(*,100) ' C0ms1ir0_coli: m32 = 0'
            write(*,111) ' C0ms1ir0_coli: p12 = ',p12
            write(*,111) ' C0ms1ir0_coli: p23 = ',p23
            write(*,111) ' C0ms1ir0_coli: p13 = ',p13
            write(*,111) ' C0ms1ir0_coli: m12 = ',m12
            write(*,111) ' C0ms1ir0_coli: m22 = ',m22
            write(*,111) ' C0ms1ir0_coli: m32 = ',m32
            flag(3)=.false.
          endif
#endif
          
          C0ms1ir0_coli = cspcon_coli(mm22/(-q23),(q13-q23)/mm22,
     &        (q13-q23)/(-q23),(-q13)/(-q23),1d0,
     &        -dreal(q13-q23))
     &        - cspenc_coli(q23/q13,dreal(q13-q23))
     &        + (cln_coli(mm22/(q13-q23),dreal(q13-q23)) 
     &        - cln_coli((q13-q23)/q13,-dreal(q13-q23)))
     &        * cln_coli(q23/q13,dreal(q13-q23))
          C0ms1ir0_coli = C0ms1ir0_coli/(q13-q23)

        end if
       
#ifdef CHECK
        if(conscheck) then
          
c     ocoli: C0ms0e
          if (m32.ne.cd0) then
            
            if(q13.eq.cd0) return

            x1  = (q23-m32)/m32
            ix1 = dreal(q23-m32)
            x2  = mm22/(q23-m32)
            z1  = (q23-q13)/q13
            iz1 = dreal(q13-q23)
            z2  = mm22/(q23-q13)
            
            C0ms1ir0_check = 
     &          cspcon_coli(-x2,-1d0/z2,x2/z2,1d0-x2/z2,ix1,-iz1)
     &          + cspcon_coli(-1d0/x1,-z1,z1/x1,1d0-z1/x1,ix1,-iz1)
     &          - cspenc_coli(1d0+z1,iz1) - pi2_6
     &          +(cln_coli(-z2,iz1) - cln_coli(-z1,-iz1))
     &          * (cln_coli(m32/(m32-q13),dreal(q13)) 
     &          + cln_coli(-x1,-ix1))
            C0ms1ir0_check = C0ms1ir0_check/(q13-q23)
            
          else
            ix1 = dreal(q23)
            x2  = mm22/q23
            z1  = (q23-q13)/q13
            iz1 = dreal(q13-q23)
            z2  = mm22/(q23-q13)
            
            C0ms1ir0_check = 
     &          cspcon_coli(-x2,-1d0/z2,x2/z2,1d0-x2/z2,ix1,-iz1)
     &          -cspenc_coli(1d0+z1,iz1)
     &          +(cln_coli(-z2,iz1) - cln_coli(-z1,-iz1))
     &          *cln_coli(q23/q13,dreal(q13-q23))

            C0ms1ir0_check = C0ms1ir0_check/(q13-q23)
            
          end if

          if(cdabs(C0ms1ir0_coli/C0ms1ir0_check-1d0).gt.testacc) then
            write(testout,100) ' C0ms1ir0_coli inconsistent'
            write(testout,111) ' C0ms1ir0_coli: p12 = ',p12
            write(testout,111) ' C0ms1ir0_coli: p23 = ',p23
            write(testout,111) ' C0ms1ir0_coli: p13 = ',p13
            write(testout,111) ' C0ms1ir0_coli: m12 = ',m12
            write(testout,111) ' C0ms1ir0_coli: m22 = ',m22
            write(testout,111) ' C0ms1ir0_coli: m32 = ',m32
            write(testout,111) ' C0ms1ir0_coli  = ',C0ms1ir0_coli
            write(testout,111) ' C0ms1ir0_check = ',C0ms1ir0_check
            write(testout,111) ' diff        = ',
     &          cdabs(C0ms1ir0_coli/C0ms1ir0_check-1d0)
            
            write(*,100) ' C0ms1ir0_coli inconsistent'
            write(*,111) ' C0ms1ir0_coli: p12 = ',p12
            write(*,111) ' C0ms1ir0_coli: p23 = ',p23
            write(*,111) ' C0ms1ir0_coli: p13 = ',p13
            write(*,111) ' C0ms1ir0_coli: m12 = ',m12
            write(*,111) ' C0ms1ir0_coli: m22 = ',m22
            write(*,111) ' C0ms1ir0_coli: m32 = ',m32
            write(*,111) ' C0ms1ir0_coli  = ',C0ms1ir0_coli
            write(*,111) ' C0ms1ir0_check = ',C0ms1ir0_check
            write(*,101) ' diff        = ',
     &          cdabs(C0ms1ir0_coli/C0ms1ir0_check-1d0)
            stop
          endif
          
c Dittmaier Nucl. Phys. B675 (2003) 447   (B.2)          
          if (m32.ne.cd0) then
            
            C0ms1ir0_check = 
     &          -2d0*cspcon_coli(m32-q23,1d0/(m32-q13),
     &          (m32-q23)/(m32-q13),(q23-q13)/(m32-q13),-1d0,1d0)
c     &          2d0*cspcon_coli(m32-q13,1d0/(m32-q23),
c     &          (m32-q13)/(m32-q23),(q13-q23)/(m32-q23),-1d0,1d0)
     &          + cspenc_coli(q23/m32,1d0)
     &          - cspenc_coli(q13/m32,1d0)
     &          + cln_coli((m32-q23)/mm22,-1d0)
     &          *cln_coli((m32-q23)/m32,-1d0)
     &          - cln_coli((m32-q13)/mm22,-1d0)
     &          *cln_coli((m32-q13)/m32,-1d0)
            C0ms1ir0_check = C0ms1ir0_check/(q23-q13)
            
            if(cdabs(C0ms1ir0_coli/C0ms1ir0_check-1d0).gt.testacc) then
              write(testout,100) ' C0ms1ir0_coli inconsistent'
              write(testout,111) ' C0ms1ir0_coli: p12 = ',p12
              write(testout,111) ' C0ms1ir0_coli: p23 = ',p23
              write(testout,111) ' C0ms1ir0_coli: p13 = ',p13
              write(testout,111) ' C0ms1ir0_coli: m12 = ',m12
              write(testout,111) ' C0ms1ir0_coli: m22 = ',m22
              write(testout,111) ' C0ms1ir0_coli: m32 = ',m32
              write(testout,111) ' C0ms1ir0_coli  = ',C0ms1ir0_coli
              write(testout,111) ' C0ms1ir0_check2= ',C0ms1ir0_check
              write(testout,101) ' diff        = ',
     &            cdabs(C0ms1ir0_coli/C0ms1ir0_check-1d0)
              
              write(*,100)   ' C0ms1ir0_coli inconsistent'
              write(*,111) ' C0ms1ir0_coli: p12 = ',p12
              write(*,111) ' C0ms1ir0_coli: p23 = ',p23
              write(*,111) ' C0ms1ir0_coli: p13 = ',p13
              write(*,111) ' C0ms1ir0_coli: m12 = ',m12
              write(*,111) ' C0ms1ir0_coli: m22 = ',m22
              write(*,111) ' C0ms1ir0_coli: m32 = ',m32
              write(*,111) ' C0ms1ir0_coli  = ',C0ms1ir0_coli
              write(*,111) ' C0ms1ir0_check2= ',C0ms1ir0_check
              write(*,101) ' diff        = ',
     &            cdabs(C0ms1ir0_coli/C0ms1ir0_check-1d0)
              stop
            endif
            
          endif
        endif
#endif
          
      else
*     no soft singularity and 0 zero masses
*
*                   | p23             
*                   |   
*                  / \         
*            m22  /2 3\ m32       
*                /  1  \         
*          0 ----------------  p13
*                m12=m22             
*

#ifdef CHECK
        if(m12.eq.cd0.or.m22.eq.cd0.or.m22.ne.m12
     &      .or.p12.ne.cd0.or.
     &      ps23.eq.cd0.and.m32.eq.cd0.or.
     &      ps13.eq.cd0.and.m32.eq.cd0)
     &      then
          write(*,100) ' C0ms1ir0_coli: case with no zero mass',
     &        ': wrong arguments'
          write(*,111) ' C0ms1ir0_coli: p12 = ',p12,ps12
          write(*,111) ' C0ms1ir0_coli: p23 = ',p23,ps23
          write(*,111) ' C0ms1ir0_coli: p13 = ',p13,ps13
          write(*,111) ' C0ms1ir0_coli: m12 = ',m12
          write(*,111) ' C0ms1ir0_coli: m22 = ',m22
          write(*,111) ' C0ms1ir0_coli: m32 = ',m32
	  write(*,*) ' C0ms1ir0_coli: test ',
     &	  m12.eq.cd0,m22.eq.cd0,m22.ne.m12
     &      ,p12.ne.cd0,
     &      ps23.eq.cd0.and.m32.eq.cd0,
     &      ps13.eq.cd0.and.m32.eq.cd0
          stop
        endif
#endif

#ifdef ALLCHECK
c ppenj        
        if(flag(4))then
          write(*,100) ' C0ms1ir0_coli: (B.3)',
     &          '    case with 0 zero masses',
     &          '    not yet tested in physical process'
          write(*,111) ' C0ms1ir0_coli: p12 = ',p12
          write(*,111) ' C0ms1ir0_coli: p23 = ',p23
          write(*,111) ' C0ms1ir0_coli: p13 = ',p13
          write(*,111) ' C0ms1ir0_coli: m12 = ',m12
          write(*,111) ' C0ms1ir0_coli: m22 = ',m22
          write(*,111) ' C0ms1ir0_coli: m32 = ',m32
          flag(4)=.false.
        endif       
#endif

        C0ms1ir0_coli=cspenc_coli(-p23/(m32-p23),-1d0)
     &      - cspenc_coli(-p13/(m32-p13),-1d0)
     &      - .5d0*cln_coli((m32-p13)/(m32-p23),dreal(p23-p13))**2
     &      + cln_coli((m32-p13)/(m32-p23),dreal(p23-p13))
     &      *cln_coli((m32-p13)/(sqrt(m12*m22)*coliminfscale2),-1d0)
        
        C0ms1ir0_coli = C0ms1ir0_coli/(p13-p23)
        
#ifdef CHECK
        if(conscheck.and.m32.ne.cd0) then


c Dittmaier Nucl. Phys. B675 (2003) 447   (B.3)          
          C0ms1ir0_check = 
     &        + cspenc_coli(p23/m32,1d0)
     &        - cspenc_coli(p13/m32,1d0)
     &        + (cln_coli((m32-p23)/sqrt(m22*coliminfscale2*m32),
     &        -1d0))**2
     &        - (cln_coli((m32-p13)/sqrt(m22*coliminfscale2*m32),
     &        -1d0))**2
          C0ms1ir0_check = C0ms1ir0_check/(p23-p13)
          
          if(cdabs(C0ms1ir0_coli/C0ms1ir0_check-1d0).gt.testacc) then
            write(testout,100) ' C0ms1ir0_coli inconsistent'
            write(testout,111) ' C0ms1ir0_coli: p12 = ',p12
            write(testout,111) ' C0ms1ir0_coli: p23 = ',p23
            write(testout,111) ' C0ms1ir0_coli: p13 = ',p13
            write(testout,111) ' C0ms1ir0_coli: m12 = ',m12
            write(testout,111) ' C0ms1ir0_coli: m22 = ',m22
            write(testout,111) ' C0ms1ir0_coli: m32 = ',m32
            write(testout,111) ' C0ms1ir0_coli  = ',C0ms1ir0_coli
            write(testout,111) ' C0ms1ir0_check = ',C0ms1ir0_check
            write(testout,101) ' diff        = ',
     &          cdabs(C0ms1ir0_coli/C0ms1ir0_check-1d0)
            
            write(*,100) ' C0ms1ir0_coli inconsistent'
            write(*,111) ' C0ms1ir0_coli: p12 = ',p12
            write(*,111) ' C0ms1ir0_coli: p23 = ',p23
            write(*,111) ' C0ms1ir0_coli: p13 = ',p13
            write(*,111) ' C0ms1ir0_coli: m12 = ',m12
            write(*,111) ' C0ms1ir0_coli: m22 = ',m22
            write(*,111) ' C0ms1ir0_coli: m32 = ',m32
            write(*,111) ' C0ms1ir0_coli  = ',C0ms1ir0_coli
            write(*,111) ' C0ms1ir0_check = ',C0ms1ir0_check
            write(*,101) ' diff        = ',
     &          cdabs(C0ms1ir0_coli/C0ms1ir0_check-1d0)
c            stop
          endif
          
        endif
#endif
        
      endif  
      
      end
************************************************************************
      function C0ms1ir1_coli(p12,p23,p13,m12,m22,m32)
************************************************************************
*     single mass-singular 3-point function                            *
*                                                                      *
*     assumes single mass singularity                                  *
*     and single soft singularity with                                 *
*     p12, m12, m22 small                                              *
*     p23, p13=m32 finite                                              *
*                                                                      *
*                            p23                                       *
*                             |                                        *
*                             |                                        *
*                            / \                                       *
*                           /   \                                      *
*              m22 small   /     \   m32                               *
*                         /       \                                    *
*                        /         \                                   *
*       m22 = p12  ---------------------  p13 = m32                    *
*       small             m12 small                                    *
*                                                                      *
*----------------------------------------------------------------------*
*     21.08.08 Ansgar Denner        last changed  03.10.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p13
      complex*16 m12,m22,m32
      complex*16 C0ms1ir1_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps13
      complex*16 ms12,ms22,ms32
      logical flag(0:4)
      data    flag /5*.true./
      save    flag  
#endif

c      complex*16 mm32,m2(3)
      complex*16 cln_coli,cspenc_coli
c      integer    nsoft
c      integer    i,j,k
c      logical    soft(3),onsh(3,3)

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

c      write(*,*) ' C0ms1ir1_coli in', p12,p23,p13,m12,m22,m32,flag

#ifdef CHECK
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if(argcheck)then
        ms12 = elimminf2_coli(m12)
        ms22 = elimminf2_coli(m22)
        ms32 = elimminf2_coli(m32)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps13 = elimminf2_coli(p13)

        if(ms12.ne.cd0.or.ms22.ne.cd0.or.ps12.ne.cd0
     &        .or.ms32.eq.cd0.or.p13.ne.m32.or.p23.eq.m32
     &        .or.ps23.eq.ms32
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p13).ne.0d0) then
          write(*,100) ' C0ms1ir1_coli called improperly'
          write(*,100) '     or case not implemented    '
          write(*,111) ' C0ms1ir1_coli: p12 = ',p12,ps12
          write(*,111) ' C0ms1ir1_coli: p23 = ',p23,ps23
          write(*,111) ' C0ms1ir1_coli: p13 = ',p13,ps13
          write(*,111) ' C0ms1ir1_coli: m12 = ',m12,ms12
          write(*,111) ' C0ms1ir1_coli: m22 = ',m22,ms22
          write(*,111) ' C0ms1ir1_coli: m32 = ',m32,ms32
          write(*,111) ' C0ms1ir1_coli: tests ',
     &     ms12.ne.cd0,ms22.ne.cd0,ps12.ne.cd0
     &        ,ms32.eq.cd0,p13.ne.m32,p23.eq.m32
     &        ,ps23.eq.ms32
     &        ,dimag(p12).ne.0d0,dimag(p23).ne.0d0
     &        ,dimag(p13).ne.0d0   
          stop
        end if
      end if
#endif

      if(p12.eq.m22)then
        if(m22.eq.cd0)then
* one soft singularity and 2 zero masses
*
*                   | p23             
*                   |   
*                  / \         
*              0  /2 3\ m32       
*                /  1  \         
*          0 ----------------  p13=m32
*                   0             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.cd0.or.m32.eq.cd0
     &        .or.p12.ne.cd0.or.p13.ne.m32.or.p23.eq.m32)
     &        then
            write(*,100) ' C0ms1ir1_coli:  (B.12)',
     &          '    case with 1 soft singularity and 2 zero masses',
     &          '    wrong arguments'
            write(*,111) ' C0ms1ir1_coli: p12 = ',p12
            write(*,111) ' C0ms1ir1_coli: p23 = ',p23
            write(*,111) ' C0ms1ir1_coli: p13 = ',p13
            write(*,111) ' C0ms1ir1_coli: m12 = ',m12
            write(*,111) ' C0ms1ir1_coli: m22 = ',m22
            write(*,111) ' C0ms1ir1_coli: m32 = ',m32
            stop
          endif
#endif
        
#ifdef ALLCHECK
c ppttbb dim
          if(flag(1))then
            write(*,100) ' C0ms1ir1_coli: (B.12)',
     &          '    case with 1 soft singularity and 2 zero masses',
     &          '    not yet tested in physical process'
            write(*,111) ' C0ms1ir1_coli: p12 = ',p12
            write(*,111) ' C0ms1ir1_coli: p23 = ',p23
            write(*,111) ' C0ms1ir1_coli: p13 = ',p13
            write(*,111) ' C0ms1ir1_coli: m12 = ',m12
            write(*,111) ' C0ms1ir1_coli: m22 = ',m22
            write(*,111) ' C0ms1ir1_coli: m32 = ',m32
            flag(1)=.false.
          endif
#endif

c Dittmaier Nucl. Phys. B675 (2003) 447   (B.12)          
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002    triangle 4
          C0ms1ir1_coli =
#ifdef SING
     &        .5d0*delta2ir
     &        -delta1ir*(cln_coli((m32-p23)/muir2,-1d0)
     &        -.5d0*cln_coli(m32/muir2,-1d0))
#endif
     &        +.5d0*(cln_coli((m32-p23)/muir2,-1d0))**2
     &        -.25d0*(cln_coli(m32/muir2,-1d0))**2
     &        -cspenc_coli(-p23/(m32-p23),-1d0)
          
          C0ms1ir1_coli = C0ms1ir1_coli/(p23-m32)
        
        else
* one soft singularity and 2 zero masses
*
*                   | p23             
*                   |   
*                  / \         
*            m22  /2 3\ m32       
*                /  1  \         
*    p12=m22 ----------------  p13=m32
*                   0             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.eq.cd0.or.m32.eq.cd0
     &        .or.p12.ne.m22.or.p13.ne.m32)
     &        then
            write(*,100) ' C0ms1ir1_coli:  (B.8)',
     &          '    case with 1 soft singularity and 1 zero mass',
     &          '    wrong arguments'
            write(*,111) ' C0ms1ir1_coli: p12 = ',p12
            write(*,111) ' C0ms1ir1_coli: p23 = ',p23
            write(*,111) ' C0ms1ir1_coli: p13 = ',p13
            write(*,111) ' C0ms1ir1_coli: m12 = ',m12
            write(*,111) ' C0ms1ir1_coli: m22 = ',m22
            write(*,111) ' C0ms1ir1_coli: m32 = ',m32
            stop
          endif
#endif
        
#ifdef CHECK
          if(flag(2))then
            write(*,100) ' C0ms1ir1_coli: (B.8)',
     &          '    case with 1 soft singularity and 1 zero mass',
     &          '    not yet tested in physical process'
            write(*,111) ' C0ms1ir1_coli: p12 = ',p12
            write(*,111) ' C0ms1ir1_coli: p23 = ',p23
            write(*,111) ' C0ms1ir1_coli: p13 = ',p13
            write(*,111) ' C0ms1ir1_coli: m12 = ',m12
            write(*,111) ' C0ms1ir1_coli: m22 = ',m22
            write(*,111) ' C0ms1ir1_coli: m32 = ',m32
            flag(2)=.false.
          endif
#endif        

c W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 iva
c according to Dittmaier Nucl. Phys. B675 (2003) 447   (B.8)          
c Ellis, Zanderighi, JHEP 0802 (2008) 002    triangle 6

c single mass-singular case  m2 = small
          if(p23.ne.cd0) then
            C0ms1ir1_coli = cspenc_coli(p23/m32,1d0)
     &        + cln_coli((m32-p23)/(sqrt(m22*m32)*coliminfscale),-1d0)
     &          * cln_coli((m32-p23)*sqrt(m22/m32)*coliminfscale/muir2
     &          ,-1d0)
#ifdef SING
     &        - cln_coli((m32-p23)/(sqrt(m22*m32)*coliminfscale),-1d0)
     &          *delta1ir
     &          - 0.25*colishiftms2  
#endif
            C0ms1ir1_coli=-C0ms1ir1_coli/(m32-p23)
          else
            C0ms1ir1_coli=0.5d0*log(muir2/(sqrt(m22*m32)*coliminfscale))
     &          *log(m32/(m22*coliminfscale2))
#ifdef SING
     &          + 0.25d0*colishiftms2
#endif
            C0ms1ir1_coli=C0ms1ir1_coli/m32
          end if
          
        endif
      
      elseif(m22.eq.cd0.and.p12.eq.m12)then

* 0 soft singularity and 1 zero masses
*
*                   | p23             
*                   |   
*                  / \         
*              0  /2 3\ m32       
*                /  1  \         
*    p12=m12 ----------------  p13=m32
*                  m12            
*

#ifdef CHECK
        if(m22.ne.cd0.or.m12.eq.cd0.or.m32.eq.cd0
     &      .or.p12.ne.m12.or.p13.ne.m32.or.p23.eq.m32)
     &      then
          write(*,100) ' C0ms1ir1_coli: (B.10)',
     &          '    case with 0 soft singularity and 1 zero mass',
     &          '    wrong arguments'
          write(*,111) ' C0ms1ir1_coli: p12 = ',p12
          write(*,111) ' C0ms1ir1_coli: p23 = ',p23
          write(*,111) ' C0ms1ir1_coli: p13 = ',p13
          write(*,111) ' C0ms1ir1_coli: m12 = ',m12
          write(*,111) ' C0ms1ir1_coli: m22 = ',m22
          write(*,111) ' C0ms1ir1_coli: m32 = ',m32
          stop
        endif
#endif
        
#ifdef CHECK
        if(flag(3))then
          write(*,100) ' C0ms1ir1_coli: (B.10)',
     &          '    case with 0 soft singularity and 1 zero mass',
     &          '    not yet tested in physical process'
          write(*,111) ' C0ms1ir1_coli: p12 = ',p12
          write(*,111) ' C0ms1ir1_coli: p23 = ',p23
          write(*,111) ' C0ms1ir1_coli: p13 = ',p13
          write(*,111) ' C0ms1ir1_coli: m12 = ',m12
          write(*,111) ' C0ms1ir1_coli: m22 = ',m22
          write(*,111) ' C0ms1ir1_coli: m32 = ',m32
          flag(3)=.false.
        endif
#endif        

c case: Dittmaier Nucl. Phys. B675 (2003) 447  (B.10)

        C0ms1ir1_coli = -cspenc_coli(p23/m32,1d0)
     &      - cln_coli((m32-p23)/(m12*coliminfscale2),-1d0)
     &      * cln_coli((m32-p23)/m32,-1d0)
     &      - 0.25d0*(log(m32/(m12*coliminfscale2)))**2 
     &      - 2.5d0*pi2_6
#ifdef SING
     &      - 0.25*colishiftms2
#endif
        C0ms1ir1_coli=C0ms1ir1_coli/(m32-p23)

      elseif(p12.eq.cd0.and.m22.eq.m12)then

* 0 soft singularity and 0 zero masses
*
*                   | p23             
*                   |   
*                  / \         
*         m22=m12 /2 3\ m32       
*                /  1  \         
*      p12=0 ----------------  p13=m32
*                  m12            
*

#ifdef CHECK
        if(m22.ne.m12.or.m12.eq.cd0.or.m32.eq.cd0
     &      .or.p12.ne.cd0.or.p13.ne.m32.or.p23.eq.m32)
     &      then
          write(*,100) ' C0ms1ir1_coli: (B.11)',
     &          '    case with 0 soft singularity and 0 zero mass',
     &          '    not yet tested in physical process'
          write(*,111) ' C0ms1ir1_coli: p12 = ',p12
          write(*,111) ' C0ms1ir1_coli: p23 = ',p23
          write(*,111) ' C0ms1ir1_coli: p13 = ',p13
          write(*,111) ' C0ms1ir1_coli: m12 = ',m12
          write(*,111) ' C0ms1ir1_coli: m22 = ',m22
          write(*,111) ' C0ms1ir1_coli: m32 = ',m32
          stop
        endif
#endif
        
#ifdef CHECK
        if(flag(4))then
          write(*,100) ' C0ms1ir1_coli: (B.11)',
     &          '    case with 0 soft singularity and 0 zero mass',
     &          '    not yet tested in physical process'
          write(*,111) ' C0ms1ir1_coli: p12 = ',p12
          write(*,111) ' C0ms1ir1_coli: p23 = ',p23
          write(*,111) ' C0ms1ir1_coli: p13 = ',p13
          write(*,111) ' C0ms1ir1_coli: m12 = ',m12
          write(*,111) ' C0ms1ir1_coli: m22 = ',m22
          write(*,111) ' C0ms1ir1_coli: m32 = ',m32
          flag(4)=.false.
        endif
#endif        

c case: Dittmaier Nucl. Phys. B675 (2003) 447  (B.10)

        C0ms1ir1_coli = -cspenc_coli(p23/m32,1d0)
     &      - (cln_coli((m32-p23)/sqrt(m32*m12*coliminfscale2),-1d0))**2
     &      - .5d0*pi2_6
#ifdef SING
     &      - 0.25*colishiftms2
#endif
        C0ms1ir1_coli=C0ms1ir1_coli/(m32-p23)

#ifdef CHECK
      else
        write(*,100) ' C0ms1ir1_coli called with improper arguments'
        write(*,100) '     or case not implemented    '
        write(*,111) ' C0ms1ir1_coli: p12 = ',p12
        write(*,111) ' C0ms1ir1_coli: p23 = ',p23
        write(*,111) ' C0ms1ir1_coli: p13 = ',p13
        write(*,111) ' C0ms1ir1_coli: m12 = ',m12
        write(*,111) ' C0ms1ir1_coli: m22 = ',m22
        write(*,111) ' C0ms1ir1_coli: m32 = ',m32
        stop
#endif

      endif
      end
************************************************************************
      function C0ms2ir1_coli(p12,p23,p13,m12,m22,m32)
************************************************************************
*     mass-singular 3-point function                                   *
*                                                                      *
*     assumes double mass singularity with                             *
*     p12, p13, m12, m22, m32 small                                    *
*     p23 finite                                                       *
*                                                                      *
*                            p23                                       *
*                             |                                        *
*                             |                                        *
*                            / \                                       *
*                           /   \                                      *
*              m22 small   /     \   m32 small                         *
*                         /       \                                    *
*                        /         \                                   *
*       m22 = p12  ---------------------  p13 = m32                    *
*       small             m12 small       small                        *
*                                                                      *
*----------------------------------------------------------------------*
*     20.08.08 Ansgar Denner        last changed  26.09.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p13
      complex*16 m12,m22,m32
      complex*16 C0ms2ir1_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps13
      complex*16 ms12,ms22,ms32
      logical flag(0:6)
      data    flag /7*.true./
      save    flag  
#endif

      complex*16 mm32,m2(3)
      complex*16 cln_coli
      integer    nsoft
      integer    i,j,k
      logical    soft(3),onsh(3,3)

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

#ifdef CHECK
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if(argcheck)then
        ms12 = elimminf2_coli(m12)
        ms22 = elimminf2_coli(m22)
        ms32 = elimminf2_coli(m32)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps13 = elimminf2_coli(p13)

        if(ms12.ne.cd0.or.ms22.ne.cd0.or.ms32.ne.cd0
     &        .or.ps12.ne.cd0.or.ps13.ne.cd0.or.ps23.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p13).ne.0d0) then
          write(*,100) ' C0ms2ir1_coli called improperly:'
          write(*,100) '     or case not implemented    '
          write(*,111) ' C0ms2ir1_coli: p12 = ',p12,ps12
          write(*,111) ' C0ms2ir1_coli: p23 = ',p23,ps23
          write(*,111) ' C0ms2ir1_coli: p13 = ',p13,ps13
          write(*,111) ' C0ms2ir1_coli: m12 = ',m12,ms12
          write(*,111) ' C0ms2ir1_coli: m22 = ',m22,ms22
          write(*,111) ' C0ms2ir1_coli: m32 = ',m32,ms32
          stop
        end if
      end if
#endif

      m2(1)=m12
      m2(2)=m22
      m2(3)=m32

c determine on-shell momenta
      onsh(1,2)=p12.eq.m22
      onsh(1,3)=p13.eq.m32
      onsh(2,1)=p12.eq.m12
      onsh(2,3)=p23.eq.m32
      onsh(3,1)=p13.eq.m12
      onsh(3,2)=p23.eq.m22


c count/determine soft singularities
      nsoft=0 
      do i=1,3
        do j=1,2
          if(i.ne.j)then
            do k=j,3
              if(k.ne.i.and.k.ne.j)then
                soft(i)=m2(i).eq.cd0.and.onsh(i,j).and.onsh(i,k)
c                soft(i,k,j)=soft(i,j,k)
                if(soft(i)) nsoft=nsoft+1
              endif
            enddo
          endif
        enddo
      enddo

      if(nsoft.eq.1)then

        if(m22.eq.cd0.and.m32.eq.cd0)then
* one soft singularity and 3 zero masses
*
*                   | p23             
*                   |   
*                  / \         
*              0  /2 3\ 0       
*                /  1  \         
*          0 ----------------  0
*                   0             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.cd0.or.m32.ne.cd0
     &        .or.p12.ne.cd0.or.p13.ne.cd0.or.p23.eq.cd0)
     &        then
            write(*,100) ' C0ms2ir1_coli:  (B.16)',
     &          '    case with 1 soft singularity and 3 zero masses ',
     &          '    wrong arguments'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            stop
          endif
#endif

#ifdef ALLCHECK
c ppenj          
          if(flag(1))then
            write(*,100) ' C0ms2ir1_coli: (B.16)',
     &          '    case with 1 soft singularity and 3 zero masses ',
     &          '    not yet tested in physical process'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            flag(1)=.false.
          endif
#endif

          
c Dittmaier Nucl. Phys. B675 (2003) 447   (B.16)
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002    triangle 1
          C0ms2ir1_coli =
#ifdef SING
     &        delta2ir
     &        -delta1ir*cln_coli(-p23/muir2,-1d0)
#endif
     &        +.5d0*(cln_coli(-p23/muir2,-1d0))**2
     &        - pi2_6
          
          C0ms2ir1_coli = C0ms2ir1_coli/p23

        elseif(m22.eq.cd0.or.m32.eq.cd0)then
* one soft singularity and 2 zero masses
*
*                   | p23             
*                   |   
*                  / \         
*              0  /2 3\ m32       
*                /  1  \         
*          0 ----------------  p13=m32
*                   0             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.p12.or.m32.ne.p13
     &        .or.m22*m32.ne.cd0.or.p23.eq.cd0)
     &        then
            write(*,100) ' C0ms2ir1_coli: (B.12) m->0',
     &          '    case with 1 soft singularity and 2 zero masses',
     &          '    wrong arguments'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(2))then
            write(*,100) ' C0ms2ir1_coli: (B.12) m->0',
     &          '    case with 1 soft singularity and 2 zero masses',
     &          '    not yet tested for physical process'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            flag(2)=.false.
          endif
#endif

          mm32=(m22+m32)*coliminfscale2
          
c Dittmaier Nucl. Phys. B675 (2003) 447   (B.12)   m -> 0
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002     triangle 4
          C0ms2ir1_coli =
#ifdef SING
     &        .5d0*delta2ir
     &        -delta1ir*(cln_coli(-p23/muir2,-1d0)
     &        -.5d0*cln_coli(mm32/muir2,-1d0))
     &        -.25d0*colishiftms2
#endif
     &        +.5d0*(cln_coli(-p23/muir2,-1d0))**2
     &        -.25d0*(cln_coli(mm32/muir2,-1d0))**2
     &        - pi2_6 
 
          C0ms2ir1_coli = C0ms2ir1_coli/p23

        else
* one soft singularity and 1 zero mass
*
*                   | p23             
*                   |   
*                  / \         
*            m22  /2 3\ m32       
*                /  1  \         
*     p12=m22 ----------------  p13=m32
*                   0             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.p12.or.m32.ne.p13
     &        .or.m22.eq.cd0.or.m32.eq.cd0.or.p23.eq.cd0)
     &        then
            write(*,100) ' C0ms2ir1_coli: (B.9)',
     &          '    case with 1 soft singularity and 1 zero mass:',
     &          '    wrong arguments'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            stop
          endif
#endif
          

#ifdef ALLCHECK
c ppenj
c ppttbb
          if(flag(3))then
            write(*,100) ' C0ms2ir1_coli: (B.9)',
     &          '    case with 1 soft singularity and 1 zero mass ',
     &          '    not yet tested in physical process'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            flag(3)=.false.
          endif
#endif

c  according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 ivb
c Dittmaier Nucl. Phys. B675 (2003) 447   (B.9)
          C0ms2ir1_coli = pi2_6
     *        +.25d0*cln_coli(-m22*coliminfscale2/p23,1d0)**2
     *        +.25d0*cln_coli(-m32*coliminfscale2/p23,1d0)**2
     *        + cln_coli(-p23/muir2,-1d0)
     &        *cln_coli(-sqrt(m22*m32)*coliminfscale2/p23,1d0)
#ifdef SING
     &        +.5d0*colishiftms2
     &        -delta1ir
     &        *cln_coli(-sqrt(m22*m32)*coliminfscale2/p23,1d0)
#endif
          C0ms2ir1_coli=-C0ms2ir1_coli/p23

        endif

      elseif(nsoft.eq.0)then

        if(m22.eq.cd0.and.m32.eq.cd0)then
* no soft singularity and 2 zero masses
*
*                   | p23             
*                   |   
*                  / \         
*              0  /2 3\ 0       
*                /  1  \         
*     p12=m12 ----------------  p13=m12
*                   m12             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.ne.cd0.or.m32.ne.cd0
     &        .or.p12.ne.m12.or.p13.ne.m12.or.p23.eq.cd0)
     &        then
            write(*,100) ' C0ms2ir1_coli: (B.13)',
     &          '    case with 0 soft singularity and 2 zero masses: ',
     &          '    wrong arguments'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            stop
          endif
#endif
          
#ifdef CHECK
c ppttbb
          if(flag(4))then
            write(*,100) ' C0ms2ir1_coli: (B.13)',
     &          '    case with 0 soft singularity and 2 zero masses ',
     &          '    not yet tested in physical process'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            flag(4)=.false.
          endif
#endif

c Dittmaier Nucl. Phys. B675 (2003) 447   (B.13)          
c ocoli: cC0ms00ee
          C0ms2ir1_coli =
#ifdef SING
     &        + colishiftms2 
#endif
     &        + cln_coli(-p23/(m12*coliminfscale2),-1d0)**2 
     &        + 8d0*pi2_6
 
          C0ms2ir1_coli = C0ms2ir1_coli/(2d0*p23)


          else if(m22.eq.cd0.or.m32.eq.cd0)then
* no soft singularity and 1 zero mass
*
*                   | p23             
*                   |   
*                  / \         
*              0  /2 3\ m32=m12       
*                /  1  \         
*     p12=m12 ----------------  0
*                   m12             
*

#ifdef CHECK
          if(m12.eq.cd0.or.p23.eq.cd0.or.
     &        .not.(p12.eq.m12.and.m32.eq.m12.or.
     &          p13.eq.m12.and.m22.eq.m12)) then
            write(*,100) ' C0ms2ir1_coli: (B.14)',
     &          '    case with 0 soft singularity and 1 zero masses: ',
     &          '    wrong arguments'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            stop
          endif
#endif

#ifdef ALLCHECK
c ppenj          
          if(flag(5))then
            write(*,100) ' C0ms2ir1_coli: (B.14)',
     &          '    case with 0 soft singularity and 1 zero masses ',
     &          '    not yet tested in physical process'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            flag(5)=.false.
          endif
#endif

c Dittmaier Nucl. Phys. B675 (2003) 447   (B.14)          
c ocoli: cC0ms0ee
          C0ms2ir1_coli = cln_coli(-p23/(m12*coliminfscale2),-1d0)**2 
     &        + 4d0*pi2_6
#ifdef SING
     &        + colishiftms2 
#endif
 
          C0ms2ir1_coli = C0ms2ir1_coli/(2d0*p23)


          else
* no soft singularity and no zero mass
*
*                   | p23             
*                   |   
*                  / \         
*        m22=m12  /2 3\ m32=m12       
*                /  1  \         
*         0  ----------------  0
*                   m12             
*

#ifdef CHECK
          if(m12.eq.cd0.or.p23.eq.cd0.or.
     &       p12.ne.cd0.or.p13.ne.cd0.or.
     &       m22.ne.m12.or.m32.ne.m12) then
            write(*,100) ' C0ms2ir1_coli: (B.15)',
     &          '    case with 0 soft singularity and 0 zero masses: ',
     &          '    wrong arguments'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(6))then
            write(*,100) ' C0ms2ir1_coli: (B.15)',
     &          '    case with 0 soft singularity and 0 zero masses: ',
     &          '    not yet tested in physical process'
            write(*,111) ' C0ms2ir1_coli: p12 = ',p12
            write(*,111) ' C0ms2ir1_coli: p23 = ',p23
            write(*,111) ' C0ms2ir1_coli: p13 = ',p13
            write(*,111) ' C0ms2ir1_coli: m12 = ',m12
            write(*,111) ' C0ms2ir1_coli: m22 = ',m22
            write(*,111) ' C0ms2ir1_coli: m32 = ',m32
            flag(6)=.false.
          endif
#endif

c Dittmaier Nucl. Phys. B675 (2003) 447   (B.15)          
c ocoli: cC0msee
          C0ms2ir1_coli = cln_coli(-p23/(m12*coliminfscale2),-1d0)**2 
#ifdef SING
     &        +colishiftms2
#endif
 
          C0ms2ir1_coli = C0ms2ir1_coli/(2d0*p23)

        endif

#ifdef CHECK
      else
        write(*,100) ' C0ms2ir1_coli called with improper arguments'
        write(*,100) '     or case not implemented    '
        write(*,111) ' C0ms2ir1_coli: p12 = ',p12
        write(*,111) ' C0ms2ir1_coli: p23 = ',p23
        write(*,111) ' C0ms2ir1_coli: p13 = ',p13
        write(*,111) ' C0ms2ir1_coli: m12 = ',m12
        write(*,111) ' C0ms2ir1_coli: m22 = ',m22
        write(*,111) ' C0ms2ir1_coli: m32 = ',m32
        stop
#endif
      endif
      end
