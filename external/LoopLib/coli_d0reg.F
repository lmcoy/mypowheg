#include "global_coli.h"

c#define WRITERP

************************************************************************
*                                                                      *
*     Regular scalar 4-point function                                  *
*                                                                      *
************************************************************************
*                                                                      *
*     last changed  29.12.08  Ansgar Denner                            *
*                                                                      *
************************************************************************
* Subroutines:                                                         *
* Functions:                                                           *
* D0regrp_coli,D0regrm_coli                                            *
************************************************************************
      function D0regrp_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*  scalar 4-point function  for r13 real positive                      *
*  regular case  based on general result of                            *
*        A.Denner, U.Nierste and R.Scharf, Nucl. Phys. B367 (1991) 637 *
*                                                                      *
*                     m22                                              *
*       q12  ---------------------  q23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*              m12| 1     3 | m32                                      *
*                 |         |                                          *
*                 |    4    |                                          *
*       q14  ---------------------  q34                                *
*                     m42                                              *
*                                                                      *
*----------------------------------------------------------------------*
*  29.12.08 Ansgar Denner       last changed 29.12.08 Ansgar Denner    *
************************************************************************
      implicit   none
      complex*16 D0regrp_coli
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42

      complex*16 l12,l13,l14,l23,l24,l34
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 r21,r31,r41,r32,r42,r43
      complex*16 a,b,c,d,det
      complex*16 x(2,4)
      complex*16 ch0,ch1,ch2,ch3,l1,l2
      complex*16 eta
      real*8     ir12,ir13,ir14,ir23,ir24,ir34
      real*8     ix(2,4),ipop(2)
      real*8     test0,test1,test2,test3
      real*8     test01,test23
      real*8     u,v
      complex*16 cspcos_coli,cln_coli,eta2s_coli
      integer    i

      logical    flag2,flag3
      save       flag2,flag3
      data       flag2 /.true./, flag3 /.true./
c      logical    errflag
c      save       errflag
c      data       errflag /.true./

#ifdef CHECK
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      complex*16 elimminf2_coli
      complex*16 D0regrp_check
c     complex*16 cD0reg
#endif
      integer    j
      logical flag(0:1)
      data    flag /2*.true./
      save    flag  

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) ' D0regrp_coli in ',
c     &    p12,p23,p34,p14,p13,p24,m12,m22,m32,m42

c      D0regrp_coli=cD0reg(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
c      write(*,*) ' D0regrp_coli: cD0reg= ',D0regrp_coli


 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
#ifdef CHECK
      if (argcheck) then
        ms12 = elimminf2_coli(m12)
        ms22 = elimminf2_coli(m22)
        ms32 = elimminf2_coli(m32)
        ms42 = elimminf2_coli(m42)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps34 = elimminf2_coli(p34)
        ps14 = elimminf2_coli(p14)
        ps24 = elimminf2_coli(p24)
        ps13 = elimminf2_coli(p13)

        if(ms32.eq.0d0.or.ms12.eq.cd0.or.ms22.eq.cd0.or.ms42.eq.cd0
     &      .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &      .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &      .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0regrp_coli called improperly:'
          write(*,111)' D0regrp_coli: p12 = ',p12
          write(*,111)' D0regrp_coli: p23 = ',p23
          write(*,111)' D0regrp_coli: p34 = ',p34
          write(*,111)' D0regrp_coli: p14 = ',p14
          write(*,111)' D0regrp_coli: p13 = ',p13
          write(*,111)' D0regrp_coli: p24 = ',p24
          write(*,111)' D0regrp_coli: m12 = ',m12
          write(*,111)' D0regrp_coli: m22 = ',m22
          write(*,111)' D0regrp_coli: m32 = ',m32
          write(*,111)' D0regrp_coli: m42 = ',m42
          write(*,*)' D0regrp_coli: test= ',
     &        m32.ne.0d0,m12.eq.cd0,m22.eq.cd0,m42.eq.cd0,
     &        p23.eq.m22.and.p13.eq.m12
     &        ,p23.eq.m22.and.p34.eq.m42
     &        ,p13.eq.m12.and.p34.eq.m42
     &        ,ps23.eq.cd0.and.ms22.eq.cd0
     &        ,ps34.eq.cd0.and.ms42.eq.cd0
     &        ,ps13.eq.cd0.and.ms12.eq.cd0
     &        ,ps12.eq.cd0.and.ms12.eq.cd0.and.ms22.eq.cd0
     &        ,ps14.eq.cd0.and.ms12.eq.cd0.and.ms42.eq.cd0
     &        ,ps24.eq.cd0.and.ms22.eq.cd0.and.ms22.eq.cd0
     &        ,dimag(p12).ne.0d0,dimag(p23).ne.0d0
     &        ,dimag(p34).ne.0d0,dimag(p14).ne.0d0
     &        ,dimag(p24).ne.0d0,dimag(p13).ne.0d0
          stop
        end if
      end if
#endif

      if(flag(1))then
        write(*,100)' D0regrp_coli: ',
c     &      '    not yet tested at all'
     &      '    not yet tested in physical process'
        write(*,111)' D0regrp_coli: p12 = ',p12
        write(*,111)' D0regrp_coli: p23 = ',p23
        write(*,111)' D0regrp_coli: p34 = ',p34
        write(*,111)' D0regrp_coli: p14 = ',p14
        write(*,111)' D0regrp_coli: p24 = ',p24
        write(*,111)' D0regrp_coli: p13 = ',p13
        write(*,111)' D0regrp_coli: m12 = ',m12
        write(*,111)' D0regrp_coli: m22 = ',m22
        write(*,111)' D0regrp_coli: m32 = ',m32
        write(*,111)' D0regrp_coli: m42 = ',m42
        flag(1)=.false.
      endif

      l12 = (m12+m22-p12)
      l13 = (m12+m32-p13)
      l14 = (m12+m42-p14)
      l23 = (m22+m32-p23)
      l24 = (m22+m42-p24)
      l34 = (m32+m42-p34)

      if(l12.ne.cd0)then
        r12 = l12/(2d0*m22)*(1d0+sqrt(1d0-4d0*m12*m22/l12**2))
      else
        r12 = dcmplx(0d0,1d0)*sqrt(m12/m22)
      endif
      r21 = r12*m22/m12
      if(l13.ne.cd0)then
        r13 = l13/(2d0*m32)*(1d0+sqrt(1d0-4d0*m12*m32/l13**2))
      else
        r13 = dcmplx(0d0,1d0)*sqrt(m12/m32)
      endif
      r31 = r13*m32/m12
      if(l14.ne.cd0)then
        r14 = l14/(2d0*m42)*(1d0+sqrt(1d0-4d0*m12*m42/l14**2))
      else
        r14 = dcmplx(0d0,1d0)*sqrt(m12/m42)
      endif
      r41 = r14*m42/m12
      if(l23.ne.cd0)then
        r23 = l23/(2d0*m32)*(1d0+sqrt(1d0-4d0*m22*m32/l23**2))
      else
        r23 = dcmplx(0d0,1d0)*sqrt(m22/m32)
      endif
      r32 = r23*m32/m22
      if(l24.ne.cd0)then
        r24 = l24/(2d0*m42)*(1d0+sqrt(1d0-4d0*m22*m42/l24**2))
      else
        r24 = dcmplx(0d0,1d0)*sqrt(m22/m42)
      endif
      r42 = r24*m42/m22
      if(l34.ne.cd0)then
        r34 = l34/(2d0*m42)*(1d0+sqrt(1d0-4d0*m32*m42/l34**2))
      else
        r34 = dcmplx(0d0,1d0)*sqrt(m32/m42)
      endif
      r43 = r34*m42/m32


c for testing
      r31 = 1d0/r13
      r13 = m12/m32/r13
      r42 = 1d0/r24
      r24 = m22/m42/r24

#ifdef CHECK
      if(conscheck)then
        if( dimag(r13).eq.cd0.and.dreal(r13).lt.0d0.or.
     &      dimag(m22*r31*r31).gt.0d0.or.dimag(m42*r31*r31).gt.0d0.or.
     &      dimag(r31*l23).gt.0d0.or.dimag(r31*l34).gt.0d0.or.
     &      dimag(r31*r31*l24).gt.0d0) then
          write(*,100)'  D0regrp_coli: analytic continuation doubtful'
          write(*,111)' D0regrp_coli: p12 = ',p12
          write(*,111)' D0regrp_coli: p23 = ',p23
          write(*,111)' D0regrp_coli: p34 = ',p34
          write(*,111)' D0regrp_coli: p14 = ',p14
          write(*,111)' D0regrp_coli: p13 = ',p13
          write(*,111)' D0regrp_coli: p24 = ',p24
          write(*,111)' D0regrp_coli: m12 = ',m12
          write(*,111)' D0regrp_coli: m22 = ',m22
          write(*,111)' D0regrp_coli: m32 = ',m32
          write(*,111)' D0regrp_coli: m42 = ',m42
          write(*,*)' D0regrp_coli: test= ',
     &        dimag(r13).eq.cd0.and.dreal(r13).lt.0d0,
     &        dimag(m22*r31*r31).gt.0d0,dimag(m42*r31*r31).gt.0d0,
     &        dimag(r31*l23).gt.0d0,dimag(r31*l34).gt.0d0,
     &        dimag(r31*r31*l24).gt.0d0
c          stop

        endif
      endif
#endif      

      ir12 = sign(1d0,1d0-abs(r12*r21))
      ir13 = sign(1d0,1d0-abs(r13*r31))
      ir14 = sign(1d0,1d0-abs(r14*r41))
      ir23 = sign(1d0,1d0-abs(r23*r32))
      ir24 = sign(1d0,1d0-abs(r24*r42))
      ir34 = sign(1d0,1d0-abs(r34*r43))

      a   =  m42*(l34/r42-l23 + (l12-l14/r42)*r31)
      b   =  m12*m42*(1d0/r13-r31)*(1d0/r42-r24)+l12*l34-l14*l23
      c   =  m12*(l34*r24-l23 + (l12-l14*r24)/r13)
      d   =  l23-r31*l12 - (l34-r31*l14)*r24
      det =  sqrt(b*b-4d0*a*c)
      x(1,4) = (-b+sqrt(b*b-4d0*a*c))/(2d0*a)
      x(2,4) = (-b-sqrt(b*b-4d0*a*c))/(2d0*a)
      if(abs(x(1,4)).gt.abs(x(2,4))) then
        x(2,4) = c/(a*x(1,4))
      else
        x(1,4) = c/(a*x(2,4))
      end if

      ix(1,4) = -sign(1d0,dreal(d))
      ix(2,4) = +sign(1d0,dreal(d))
      ix(1,1) =  sign(1d0,ix(1,4)*dreal(r24)) ! needed for correct
      ix(2,1) =  sign(1d0,ix(2,4)*dreal(r24)) ! continuation

#ifdef WRITERP
      write(*,*) 'D0regrp_coli l ',l12,l13,l14
      write(*,*) 'D0regrp_coli l ',l23,l24,l34
      write(*,*) 'D0regrp_coli r ',r12,r13,r14
      write(*,*) 'D0regrp_coli r ',r21,r31,r41
      write(*,*) 'D0regrp_coli r ',r23,r24,r34
      write(*,*) 'D0regrp_coli r ',r32,r42,r43
      write(*,*) 'D0regrp_coli x ',x(1,4),x(2,4)
      write(*,*) 'D0regrp_coli ir ',ir12,ir13,ir14
      write(*,*) 'D0regrp_coli ir ',ir23,ir24,ir34
      write(*,*) 'D0regrp_coli ix ',ix(1,4),ix(2,4)
#endif

      x(1,1) = x(1,4)/r24
      x(2,1) = x(2,4)/r24
c      x(1,2) = x(1,4)/r24*r13
c      x(2,2) = x(2,4)/r24*r13
c      x(1,3) = x(1,4)*r13
c      x(2,3) = x(2,4)*r13
 
      D0regrp_coli = dcmplx(0d0)
      do i=1,2
        eta=eta2s_coli(-x(i,4),1d0/r24,-ix(i,4),-ir24,-ix(i,1))

        if(eta.ne.0d0)then
          ch0   = m42*(r31*x(i,4))**2+l34*r31*x(i,4)+m32
          test0 = abs(ch0)/max(abs(m42*(r31*x(i,4))**2),
     &        abs(l34*r31*x(i,4)),abs(m32))
          ch1   = m42*(x(i,4))**2+l14*x(i,4)+m12
          test1 = abs(ch1)/
     &        max(abs(m42*(x(i,4))**2),abs(l14*x(i,4)),abs(m12))
          test01=min(test0,test1)
          ch2   = m22*(r31*x(i,1))**2+l23*r31*x(i,1)+m32
          test2 = abs(ch2)/max(abs(m22*(r31*x(i,1))**2),
     &        abs(l23*r31*x(i,1)),abs(m32))
          ch3   = m22*(x(i,1))**2+l12*x(i,1)+m12
          test3 = abs(ch3)/
     &        max(abs(m22*(x(i,1))**2),abs(l12*x(i,1)),abs(m12))
          test23 = min(test2,test3)
          if(dreal(x(i,4)).gt.0d0)then
            ipop(i)=dreal(ch0-ch1)
          elseif(dreal(x(i,1)).gt.0d0)then
            ipop(i)=dreal(ch2-ch3)
          endif
          if (test23.gt.test01) then
            l1 = cln_coli(ch2/ch3,ipop(i))
          else
            l1 = cln_coli(ch0/ch1,ipop(i))      
          end if


#ifdef CHECK
          if(conscheck)then
            if(dimag(r24).ne.0d0) then
              v=dimag(x(i,4))/dimag(r24)
              u=dimag(x(i,4)/r24)/dimag(1d0/r24)
              l2=cln_coli(m32+l23*r31*v+l34*r31*u+m22*(r31*v)**2
     &            +m42*(r31*u)**2+l24*r31*r31*u*v,-1d0)
     &            -cln_coli(m12+l12*v+l14*u+m22*v*v
     &            +m42*u*u+l24*u*v,-1d0)
              
#ifdef WRITERP
              write(*,*) 'D0regrp_coli comp1 l2= ',l2,l1,l2/l1
#endif              

            elseif(dimag(x(i,4)).ne.0d0.or.dreal(x(i,4)).gt.0d0) then
              if(test01.gt.test23)then
                l2=cln_coli(ch0,-1d0)-cln_coli(ch1,-1d0)
              elseif(test0.gt.test1)then
                l2=cln_coli(ch0,-1d0)-cln_coli(ch0*ch3/ch2,-1d0)
              else
                l2=cln_coli(ch1*ch2/ch3,-1d0)-cln_coli(ch1,-1d0)
              endif

#ifdef WRITERP
              write(*,*) 'D0regrp_coli comp2 l2= ',l2,l1,l2/l1
#endif

            elseif(dreal(x(i,1)).gt.0d0)then
              if(test23.gt.test01)then
                l2=cln_coli(ch2,-1d0)-cln_coli(ch3,-1d0)
              elseif(test2.gt.test3)then
                l2=cln_coli(ch2,-1d0)-cln_coli(ch2*ch1/ch0,-1d0)
              else
                l2=cln_coli(ch3*ch0/ch1,-1d0)-cln_coli(ch3,-1d0)
              endif

#ifdef WRITERP
              write(*,*) 'D0regrp_coli comp3 l2= ',l2,l1,l2/l1
#endif

            else
              if (flag3) then
                write(*,*) 
                write(*,*) 'D0regrp_coli: case for l2 not implemented'
                flag3 = .false.
              end if
              write(errout,*) 
              write(errout,*) 
     &            'D0regrp_coli: case for l2 not implemented '
              write(errout,*) 
     &            'D0regrp_coli: args ',p12,p23,p34,p14,p13,p24,
     &            m12,m22,m32,m42
              write(errout,*) 'r24 = ',r24
              write(errout,*) 'x(i,4) = ',x(i,4)
              
              write(errout,*) 'x(i,1) = ',x(i,1)
            endif 
            
c              write(*,*) 'D0regrp_coli l2= ',l2,l1,l2/l1

            if(abs(l2/l1-1d0).gt.testacc) then
              if (flag2) then
                write(*,*) 
                write(*,*) 'D0regrp_coli: eta terms possibly wrong'
                flag2 = .false.
              end if
              write(errout,*) 
              write(errout,*) 'D0regrp_coli: eta terms possibly wrong'
              write(errout,*) 
     &            'D0regrp_coli args ',p12,p23,p34,p14,p13,p24,
     &            m12,m22,m32,m42
              write(errout,*) 
     &            'D0regrp_coli chi ',ch0,ch1,ch2,ch3
              write(errout,*) 'l1 = ',l1
              write(errout,*) 'l2 = ',l2
              write(errout,*) 'diff   = ',abs((l2-l1)/l1)
            endif 
          endif
#endif


        if(eta.ne.0d0)then
          D0regrp_coli = D0regrp_coli + (2*i-3) * (
     &        + eta*(l1-log(m32/m12)) ) 
        endif
            
#ifdef WRITERP
          write(*,*) 'D0regrp_coli eta i ',
     &        eta*(l1-log(m32/m12))
#endif
          
        endif
        
        D0regrp_coli = D0regrp_coli + (2*i-3) * (
     &      - cspcos_coli(-x(i,4),r43*r31,-ix(i,4),ir34)
     &      - cspcos_coli(-x(i,4),r31/r34,-ix(i,4),-ir34)
     &      + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14)
     &      + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14)
     &      - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12)
     &      - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12)
     &      + cspcos_coli(-x(i,1),r31*r23,-ix(i,1),ir23)
     &      + cspcos_coli(-x(i,1),r31/r32,-ix(i,1),-ir23))
        
#ifdef WRITERP
        write(*,*) 'D0regrp_coli arg1 i ',
     &      -x(i,4),r43*r31,-ix(i,4),ir34
        write(*,*) 'D0regrp_coli arg2 i ',
     &      -x(i,4),r31/r34,-ix(i,4),-ir34
        write(*,*) 'D0regrp_coli arg3 i ',
     &      -x(i,4),r41,-ix(i,4),ir14
        write(*,*) 'D0regrp_coli arg4 i ',
     &      -x(i,4),1d0/r14,-ix(i,4),-ir14
        write(*,*) 'D0regrp_coli arg5 i ',
     &      -x(i,1),r21,-ix(i,1),ir12
        write(*,*) 'D0regrp_coli arg6 i ',
     &      -x(i,1),1d0/r12,-ix(i,1),-ir12
        write(*,*) 'D0regrp_coli arg7 i ',
     &      -x(i,1),r31*r23,-ix(i,1),ir23
        write(*,*) 'D0regrp_coli arg8 i ',
     &      -x(i,1),r31/r32,-ix(i,1),-ir23
        
        write(*,*) 'D0regrp_coli contr i ',
     &      - cspcos_coli(-x(i,4),r43*r31,-ix(i,4),ir34),
     &      - cspcos_coli(-x(i,4),r31/r34,-ix(i,4),-ir34),
     &      + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14),
     &      + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14),
     &      - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12),
     &      - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12),
     &      + cspcos_coli(-x(i,1),r31*r23,-ix(i,1),ir23),
     &      + cspcos_coli(-x(i,1),r31/r32,-ix(i,1),-ir23)
        write(*,*) 'D0regrp_coli sum i ',
     &      - cspcos_coli(-x(i,4),r43*r31,-ix(i,4),ir34)
     &      - cspcos_coli(-x(i,4),r31/r34,-ix(i,4),-ir34)
     &      + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14)
     &      + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14)
     &      - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12)
     &      - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12)
     &      + cspcos_coli(-x(i,1),r31*r23,-ix(i,1),ir23)
     &      + cspcos_coli(-x(i,1),r31/r32,-ix(i,1),-ir23)
#endif
        
      enddo
      
#ifdef WRITERP
      write(*,*) 'D0regrp_coli sum ',D0regrp_coli
#endif
      
      D0regrp_coli = D0regrp_coli/det
      
#ifdef WRITERP
      write(*,*) 'D0regrp_coli=   ',D0regrp_coli
#endif

#ifdef CHECK

      if (conscheck.and.(
     &    dimag(x(1,4)).eq.0d0.or.dimag(x(2,4)).eq.0d0)) then

        do j=1,3
          if(j.eq.1) then
            ix(1,4) = +sign(1d0,dreal(d))
            ix(2,4) = -sign(1d0,dreal(d))
          elseif(j.eq.2) then
            ix(1,4) = +sign(1d0,dreal(d))
            ix(2,4) = +sign(1d0,dreal(d))
          elseif(j.eq.3) then
            ix(1,4) = -sign(1d0,dreal(d))
            ix(2,4) = -sign(1d0,dreal(d))
          endif
          ix(1,1) =  sign(1d0,ix(1,4)*dreal(r24)) ! needed for correct
          ix(2,1) =  sign(1d0,ix(2,4)*dreal(r24)) ! continuation


          D0regrp_check = dcmplx(0d0)
          do i=1,2
            eta=eta2s_coli(-x(i,4),1d0/r24,-ix(i,4),-ir24,-ix(i,1))
            
            if(eta.ne.0d0)then
              ch0   = m42*(r31*x(i,4))**2+l34*r31*x(i,4)+m32
              test0 = abs(ch0)/max(abs(m42*(r31*x(i,4))**2),
     &            abs(l34*r31*x(i,4)),abs(m32))
              ch1   = m42*(x(i,4))**2+l14*x(i,4)+m12
              test1 = abs(ch1)/
     &            max(abs(m42*(x(i,4))**2),abs(l14*x(i,4)),abs(m12))
              test01=min(test0,test1)
              ch2   = m22*(r31*x(i,1))**2+l23*r31*x(i,1)+m32
              test2 = abs(ch2)/max(abs(m22*(r31*x(i,1))**2),
     &            abs(l23*r31*x(i,1)),abs(m32))
              ch3   = m22*(x(i,1))**2+l12*x(i,1)+m12
              test3 = abs(ch3)/
     &            max(abs(m22*(x(i,1))**2),abs(l12*x(i,1)),abs(m12))
              test23 = min(test2,test3)
              if(dreal(x(i,4)).gt.0d0)then
                ipop(i)=dreal(ch0-ch1)
              elseif(dreal(x(i,1)).gt.0d0)then
                ipop(i)=dreal(ch2-ch3)
              endif
              if (test23.gt.test01) then
                l1 = cln_coli(ch2/ch3,ipop(i))
              else
                l1 = cln_coli(ch0/ch1,ipop(i))      
              end if
              
              if(conscheck)then
                if(dimag(r24).ne.0d0) then
                  v=dimag(x(i,4))/dimag(r24)
                  u=dimag(x(i,4)/r24)/dimag(1d0/r24)
                  l2=cln_coli(m32+l23*r31*v+l34*r31*u+m22*(r31*v)**2
     &                +m42*(r31*u)**2+l24*r31*r31*u*v,-1d0)
     &                -cln_coli(m12+l12*v+l14*u+m22*v*v
     &                +m42*u*u+l24*u*v,-1d0)
                  
#ifdef WRITERP
                  write(*,*) 'D0regrp_check comp1 l2= ',l2,l1,l2/l1
#endif
                  
                elseif(dimag(x(i,4)).ne.0d0.or.dreal(x(i,4)).gt.0d0)then
                  if(test01.gt.test23)then
                    l2=cln_coli(ch0,-1d0)-cln_coli(ch1,-1d0)
                  elseif(test0.gt.test1)then
                    l2=cln_coli(ch0,-1d0)-cln_coli(ch0*ch3/ch2,-1d0)
                  else
                    l2=cln_coli(ch1*ch2/ch3,-1d0)-cln_coli(ch1,-1d0)
                  endif
                  
#ifdef WRITERP
                  write(*,*) 'D0regrp_check comp2 l2= ',l2,l1,l2/l1
#endif
                  
                elseif(dreal(x(i,1)).gt.0d0)then
                  if(test23.gt.test01)then
                    l2=cln_coli(ch2,-1d0)-cln_coli(ch3,-1d0)
                  elseif(test2.gt.test3)then
                    l2=cln_coli(ch2,-1d0)-cln_coli(ch2*ch1/ch0,-1d0)
                  else
                    l2=cln_coli(ch3*ch0/ch1,-1d0)-cln_coli(ch3,-1d0)
                  endif
                  
#ifdef WRITERP
                  write(*,*) 'D0regrp_check comp3 l2= ',l2,l1,l2/l1
#endif
                  
                else
                  if (flag3) then
                    write(*,*) 
                    write(*,*) '
     &                  D0regrp_check: case for l2 not implemented'
                    flag3 = .false.
                  end if
                  write(errout,*) 
                  write(errout,*) 
     &                'D0regrp_check: case for l2 not implemented '
                  write(errout,*) 
     &                'D0regrp_check: args ',p12,p23,p34,p14,p13,p24,
     &                m12,m22,m32,m42
                  write(errout,*) 'r24 = ',r24
                  write(errout,*) 'x(i,4) = ',x(i,4)
                  
                  write(errout,*) 'x(i,1) = ',x(i,1)
                endif 
                
c     write(*,*) 'D0regrp_check l2= ',l2,l1,l2/l1
                
                if(abs(l2/l1-1d0).gt.testacc) then
                  if (flag2) then
                    write(*,*) 
                    write(*,*) 'D0regrp_check: eta terms possibly wrong'
                    flag2 = .false.
                  end if
                  write(errout,*) 
                  write(errout,*) 
     &                'D0regrp_check: eta terms possibly wrong'
                  write(errout,*) 
     &                'D0regrp_check args ',p12,p23,p34,p14,p13,p24,
     &                m12,m22,m32,m42
                  write(errout,*) 
     &                'D0regrp_check chi ',ch0,ch1,ch2,ch3
                  write(errout,*) 'l1 = ',l1
                  write(errout,*) 'l2 = ',l2
                  write(errout,*) 'diff   = ',abs((l2-l1)/l1)
                endif 
              endif
              
              
              if(eta.ne.0d0)then
                D0regrp_check = D0regrp_check + (2*i-3) * (
     &              + eta*(l1-log(m32/m12)) ) 
              endif 
              
#ifdef WRITERP
              write(*,*) 'D0regrp_check eta i ',
     &            eta*(l1-log(m32/m12))
#endif
              
            endif
            
            D0regrp_check = D0regrp_check + (2*i-3) * (
     &          - cspcos_coli(-x(i,4),r43*r31,-ix(i,4),ir34)
     &          - cspcos_coli(-x(i,4),r31/r34,-ix(i,4),-ir34)
     &          + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14)
     &          + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14)
     &          - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12)
     &          - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12)
     &          + cspcos_coli(-x(i,1),r31*r23,-ix(i,1),ir23)
     &          + cspcos_coli(-x(i,1),r31/r32,-ix(i,1),-ir23))
            
#ifdef WRITERP
            write(*,*) 'D0regrp_check arg1 i ',
     &          -x(i,4),r43*r31,-ix(i,4),ir34
            write(*,*) 'D0regrp_check arg2 i ',
     &          -x(i,4),r31/r34,-ix(i,4),-ir34
            write(*,*) 'D0regrp_check arg3 i ',
     &          -x(i,4),r41,-ix(i,4),ir14
            write(*,*) 'D0regrp_check arg4 i ',
     &          -x(i,4),1d0/r14,-ix(i,4),-ir14
            write(*,*) 'D0regrp_check arg5 i ',
     &          -x(i,1),r21,-ix(i,1),ir12
            write(*,*) 'D0regrp_check arg6 i ',
     &          -x(i,1),1d0/r12,-ix(i,1),-ir12
            write(*,*) 'D0regrp_check arg7 i ',
     &          -x(i,1),r31*r23,-ix(i,1),ir23
            write(*,*) 'D0regrp_check arg8 i ',
     &          -x(i,1),r31/r32,-ix(i,1),-ir23
            
            write(*,*) 'D0regrp_check contr i ',
     &          - cspcos_coli(-x(i,4),r43*r31,-ix(i,4),ir34),
     &          - cspcos_coli(-x(i,4),r31/r34,-ix(i,4),-ir34),
     &          + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14),
     &          + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14),
     &          - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12),
     &          - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12),
     &          + cspcos_coli(-x(i,1),r31*r23,-ix(i,1),ir23),
     &          + cspcos_coli(-x(i,1),r31/r32,-ix(i,1),-ir23)
            write(*,*) 'D0regrp_check sum i ',
     &          - cspcos_coli(-x(i,4),r43*r31,-ix(i,4),ir34)
     &          - cspcos_coli(-x(i,4),r31/r34,-ix(i,4),-ir34)
     &          + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14)
     &          + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14)
     &          - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12)
     &          - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12)
     &          + cspcos_coli(-x(i,1),r31*r23,-ix(i,1),ir23)
     &          + cspcos_coli(-x(i,1),r31/r32,-ix(i,1),-ir23)
#endif
            
          enddo
          
#ifdef WRITERP
          write(*,*) 'D0regrp_check sum ',D0regrp_check
#endif
          
          D0regrp_check = D0regrp_check/det
          
#ifdef WRITERP
          write(*,*) 'D0regrp_check =  ',D0regrp_check
#endif

          if(cdabs(D0regrp_coli/D0regrp_check-1d0).gt.testacc) then
            write(testout,*) ' D0regrp_coli inconsistent ',j
            write(testout,111) ' D0regrp_coli: p12 = ',p12
            write(testout,111) ' D0regrp_coli: p23 = ',p23
            write(testout,111) ' D0regrp_coli: p34 = ',p34
            write(testout,111) ' D0regrp_coli: p14 = ',p14
            write(testout,111) ' D0regrp_coli: p24 = ',p24
            write(testout,111) ' D0regrp_coli: p13 = ',p13
            write(testout,111) ' D0regrp_coli: m12 = ',m12
            write(testout,111) ' D0regrp_coli: m22 = ',m22
            write(testout,111) ' D0regrp_coli: m32 = ',m32
            write(testout,111) ' D0regrp_coli: m42 = ',m42
            write(testout,111) ' D0regrp_coli  = ',D0regrp_coli
            write(testout,111) ' D0regrp_check = ',D0regrp_check
            write(testout,101) ' diff        = ',
     &          cdabs(D0regrp_coli/D0regrp_check-1d0)
            
            write(*,*)   ' D0regrp_coli inconsistent ',i
            write(*,111) ' D0regrp_coli: p12 = ',p12
            write(*,111) ' D0regrp_coli: p23 = ',p23
            write(*,111) ' D0regrp_coli: p34 = ',p34
            write(*,111) ' D0regrp_coli: p14 = ',p14
            write(*,111) ' D0regrp_coli: p24 = ',p24
            write(*,111) ' D0regrp_coli: p13 = ',p13
            write(*,111) ' D0regrp_coli: m12 = ',m12
            write(*,111) ' D0regrp_coli: m22 = ',m22
            write(*,111) ' D0regrp_coli: m32 = ',m32
            write(*,111) ' D0regrp_coli: m42 = ',m42
            write(*,111) ' D0regrp_coli  = ',D0regrp_coli
            write(*,111) ' D0regrp_check = ',D0regrp_check
            write(*,101) ' diff        = ',
     &          cdabs(D0regrp_coli/D0regrp_check-1d0)
            stop
          endif
          
        enddo
      endif
#endif

      end


#ifdef CHUNK

************************************************************************
      function D0regrm_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*  scalar 4-point function  for r13 real positive                      *
*  regular case  based on general result of                            *
*        A.Denner, U.Nierste and R.Scharf, Nucl. Phys. B367 (1991) 637 *
*                                                                      *
*                     m22                                              *
*       q12  ---------------------  q23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*              m12| 1     3 | m32                                      *
*                 |         |                                          *
*                 |    4    |                                          *
*       q14  ---------------------  q34                                *
*                     m42                                              *
*                                                                      *
*----------------------------------------------------------------------*
*  29.12.08 Ansgar Denner       last changed 29.12.08 Ansgar Denner    *
************************************************************************
      implicit   none
      complex*16 D0regrm_coli
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42

      complex*16 l12,l13,l14,l23,l24,l34
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 r21,r31,r41,r32,r42,r43
      complex*16 a,b,c,d,det
      complex*16 x(2,4)
      complex*16 ch0,ch1,ch2,ch3,l1,l2,l3,l4,l3t,l4t
      complex*16 eta1,eta2,eta3,eta4
      real*8     ir12,ir13,ir14,ir23,ir24,ir34
      real*8     ix(2,4),iqbod1,iqbod2
      real*8     pol12,pol14,pol23,pol34,q1,q2,q3,q4,qb1od,qb2od
      real*8     t12,t23,t34,t14,t1,t2,t3,t4,tb1od,tb2od,td,tqb1,tqb2
      real*8     test01,test23
      real*8     u,v
      complex*16 cspcos_coli,cln_coli,eta2s_coli
      integer    i

      logical    flag2,flag3
      save       flag2,flag3
      data       flag2 /.true./, flag3 /.true./
c      logical    errflag
c      save       errflag
c      data       errflag /.true./

#ifdef CHECK
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      complex*16 elimminf2_coli
      complex*16 D0regrm_check
c     complex*16 cD0reg
      integer    j
      logical flag(0:1)
      data    flag /2*.true./
      save    flag  
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) ' D0regrm_coli in ',
c     &    p12,p23,p34,p14,p13,p24,m12,m22,m32,m42

c      D0regrm_coli=cD0reg(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
c      write(*,*) ' D0regrm_coli: cD0reg= ',D0regrm_coli


#ifdef CHECK
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        ms12 = elimminf2_coli(m12)
        ms22 = elimminf2_coli(m22)
        ms32 = elimminf2_coli(m32)
        ms42 = elimminf2_coli(m42)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps34 = elimminf2_coli(p34)
        ps14 = elimminf2_coli(p14)
        ps24 = elimminf2_coli(p24)
        ps13 = elimminf2_coli(p13)

        if(ms32.eq.0d0.or.ms12.eq.cd0.or.ms22.eq.cd0.or.ms42.eq.cd0
     &      .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &      .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &      .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0regrm_coli called improperly:'
          write(*,111)' D0regrm_coli: p12 = ',p12
          write(*,111)' D0regrm_coli: p23 = ',p23
          write(*,111)' D0regrm_coli: p34 = ',p34
          write(*,111)' D0regrm_coli: p14 = ',p14
          write(*,111)' D0regrm_coli: p13 = ',p13
          write(*,111)' D0regrm_coli: p24 = ',p24
          write(*,111)' D0regrm_coli: m12 = ',m12
          write(*,111)' D0regrm_coli: m22 = ',m22
          write(*,111)' D0regrm_coli: m32 = ',m32
          write(*,111)' D0regrm_coli: m42 = ',m42
          write(*,*)' D0regrm_coli: test= ',
     &        m32.ne.0d0,m12.eq.cd0,m22.eq.cd0,m42.eq.cd0,
     &        p23.eq.m22.and.p13.eq.m12
     &        ,p23.eq.m22.and.p34.eq.m42
     &        ,p13.eq.m12.and.p34.eq.m42
     &        ,ps23.eq.cd0.and.ms22.eq.cd0
     &        ,ps34.eq.cd0.and.ms42.eq.cd0
     &        ,ps13.eq.cd0.and.ms12.eq.cd0
     &        ,ps12.eq.cd0.and.ms12.eq.cd0.and.ms22.eq.cd0
     &        ,ps14.eq.cd0.and.ms12.eq.cd0.and.ms42.eq.cd0
     &        ,ps24.eq.cd0.and.ms22.eq.cd0.and.ms22.eq.cd0
     &        ,dimag(p12).ne.0d0,dimag(p23).ne.0d0
     &        ,dimag(p34).ne.0d0,dimag(p14).ne.0d0
     &        ,dimag(p24).ne.0d0,dimag(p13).ne.0d0
          stop
        end if
      end if
#endif

      if(flag(1))then
        write(*,100)' D0regrm_coli: ',
     &      '    not yet tested at all'
c     &      '    not yet tested in physical process'
        write(*,111)' D0regrm_coli: p12 = ',p12
        write(*,111)' D0regrm_coli: p23 = ',p23
        write(*,111)' D0regrm_coli: p34 = ',p34
        write(*,111)' D0regrm_coli: p14 = ',p14
        write(*,111)' D0regrm_coli: p24 = ',p24
        write(*,111)' D0regrm_coli: p13 = ',p13
        write(*,111)' D0regrm_coli: m12 = ',m12
        write(*,111)' D0regrm_coli: m22 = ',m22
        write(*,111)' D0regrm_coli: m32 = ',m32
        write(*,111)' D0regrm_coli: m42 = ',m42
        flag(1)=.false.
      endif

      l12 = (m12+m22-p12)
      l13 = (m12+m32-p13)
      l14 = (m12+m42-p14)
      l23 = (m22+m32-p23)
      l24 = (m22+m42-p24)
      l34 = (m32+m42-p34)

      if(l12.ne.cd0)then
        r12 = l12/(2d0*m22)*(1d0+sqrt(1d0-4d0*m12*m22/l12**2))
      else
        r12 = dcmplx(0d0,1d0)*sqrt(m12/m22)
      endif
      r21 = r12*m22/m12
      if(l13.ne.cd0)then
        r13 = l13/(2d0*m32)*(1d0+sqrt(1d0-4d0*m12*m32/l13**2))
      else
        r13 = dcmplx(0d0,1d0)*sqrt(m12/m32)
      endif
      r31 = r13*m32/m12
      if(l14.ne.cd0)then
        r14 = l14/(2d0*m42)*(1d0+sqrt(1d0-4d0*m12*m42/l14**2))
      else
        r14 = dcmplx(0d0,1d0)*sqrt(m12/m42)
      endif
      r41 = r14*m42/m12
      if(l23.ne.cd0)then
        r23 = l23/(2d0*m32)*(1d0+sqrt(1d0-4d0*m22*m32/l23**2))
      else
        r23 = dcmplx(0d0,1d0)*sqrt(m22/m32)
      endif
      r32 = r23*m32/m22
      if(l24.ne.cd0)then
        r24 = l24/(2d0*m42)*(1d0+sqrt(1d0-4d0*m22*m42/l24**2))
      else
        r24 = dcmplx(0d0,1d0)*sqrt(m22/m42)
      endif
      r42 = r24*m42/m22
      if(l34.ne.cd0)then
        r34 = l34/(2d0*m42)*(1d0+sqrt(1d0-4d0*m32*m42/l34**2))
      else
        r34 = dcmplx(0d0,1d0)*sqrt(m32/m42)
      endif
      r43 = r34*m42/m32


c for testing
c      r31 = 1d0/r13
c      r13 = m12/m32/r13
c      r42 = 1d0/r24
c      r24 = m22/m42/r24

#ifdef CHECK
      if(conscheck)then
        if( dimag(r13).eq.cd0.and.dreal(r13).gt.0d0.or.
     &      dimag(l23-r31*l12).gt.0d0.or.dimag(l34-r31*l14).gt.0d0.or.
     &      dimag(m12*(1d0/r13-r31)).gt.0d0.or.
     &      dimag(m12*(1d0/r13-r31)*r31).gt.0d0.or.
     &      dimag(m32-m12*r31*r31).gt.0d0.or.
     &      dimag(2d0*m32+l13*r31).gt.0d0) then
          write(*,100)'  D0regrm_coli: analytic continuation doubtful'
          write(*,111)' D0regrm_coli: p12 = ',p12
          write(*,111)' D0regrm_coli: p23 = ',p23
          write(*,111)' D0regrm_coli: p34 = ',p34
          write(*,111)' D0regrm_coli: p14 = ',p14
          write(*,111)' D0regrm_coli: p13 = ',p13
          write(*,111)' D0regrm_coli: p24 = ',p24
          write(*,111)' D0regrm_coli: m12 = ',m12
          write(*,111)' D0regrm_coli: m22 = ',m22
          write(*,111)' D0regrm_coli: m32 = ',m32
          write(*,111)' D0regrm_coli: m42 = ',m42
          write(*,*)' D0regrm_coli: test= ',
     &        dimag(r13).eq.cd0.and.dreal(r13).gt.0d0,
     &      dimag(l23-r31*l12).gt.0d0,dimag(l34-r31*l14).gt.0d0,
     &      dimag(m12*(1d0/r13-r31)).gt.0d0,
     &      dimag(m12*(1d0/r13-r31)*r31).gt.0d0,
     &      dimag(m32-m12*r31*r31).gt.0d0,
     &      dimag(2d0*m32+l13*r31).gt.0d0
          write(*,*)' D0regrm_coli: test2= ',
     &      m12*(1d0/r13-r31)*r31,
     &      m32-m12*r31*r31,
     &      2d0*m32+l13*r31
        endif
      endif
#endif      

      ir12 = sign(1d0,1d0-abs(r12*r21))
      ir13 = sign(1d0,1d0-abs(r13*r31))
      ir14 = sign(1d0,1d0-abs(r14*r41))
      ir23 = sign(1d0,1d0-abs(r23*r32))
      ir24 = sign(1d0,1d0-abs(r24*r42))
      ir34 = sign(1d0,1d0-abs(r34*r43))

      a   =  m42*(l34/r42-l23 + (l12-l14/r42)*r31)
      b   =  m12*m42*(1d0/r13-r31)*(1d0/r42-r24)+l12*l34-l14*l23
      c   =  m12*(l34*r24-l23 + (l12-l14*r24)/r13)
      d   =  l23-r31*l12 - (l34-r31*l14)*r24
      det =  sqrt(b*b-4d0*a*c)
      x(1,4) = (-b+sqrt(b*b-4d0*a*c))/(2d0*a)
      x(2,4) = (-b-sqrt(b*b-4d0*a*c))/(2d0*a)
      if(abs(x(1,4)).gt.abs(x(2,4))) then
        x(2,4) = c/(a*x(1,4))
      else
        x(1,4) = c/(a*x(2,4))
      end if

      ix(1,4) = -sign(1d0,dreal(d))
      ix(2,4) = +sign(1d0,dreal(d))

      write(*,*) 'D0regrm_coli l ',l12,l13,l14
      write(*,*) 'D0regrm_coli l ',l23,l24,l34
      write(*,*) 'D0regrm_coli r ',r12,r13,r14
      write(*,*) 'D0regrm_coli r ',r21,r31,r41
      write(*,*) 'D0regrm_coli r ',r23,r24,r34
      write(*,*) 'D0regrm_coli r ',r32,r42,r43
      write(*,*) 'D0regrm_coli x ',x(1,4),x(2,4)
      write(*,*) 'D0regrm_coli ir ',ir12,ir13,ir14
      write(*,*) 'D0regrm_coli ir ',ir23,ir24,ir34
      write(*,*) 'D0regrm_coli ix ',ix(1,4),ix(2,4)

      x(1,1) = x(1,4)/r24
      x(2,1) = x(2,4)/r24
      x(1,2) = x(1,4)/r24*r13
      x(2,2) = x(2,4)/r24*r13
      x(1,3) = x(1,4)*r13
      x(2,3) = x(2,4)*r13
 
      ix(1,1) =  sign(1d0,ix(1,4)*dreal(r24)) ! needed for correct
      ix(2,1) =  sign(1d0,ix(2,4)*dreal(r24)) ! continuation
      ix(1,3) =  sign(1d0,ix(1,4)*dreal(r13)) ! needed for correct
      ix(2,3) =  sign(1d0,ix(2,4)*dreal(r13)) ! continuation
      ix(1,2) =  sign(1d0,ix(1,4)*dreal(r24)*dreal(r13)) ! needed for correct
      ix(2,2) =  sign(1d0,ix(2,4)*dreal(r24)*dreal(r13)) ! continuation

      D0regrm_coli = dcmplx(0d0)


c *** improves numerical stability in extreme cases

c check from here

      do i=1,2
c m12*P(1,x/r24,0,0) / (x/r24)
        pol12 = m12/x(i,1)+l12+m22*x(i,1)
c m12*P(1,0,0,x) / x
        pol14 = m12/x(i,4)+l14+m42*x(i,4)
c m32*P(0,r13*x/r24,1,0) / (r13*x/r24)
        pol23 = m32/x(i,2)+l23+m22*x(i,2)
c m32*P(0,0,1,r13*x) / (r13*x)
        pol34 = m32/x(i,3)+l34+m42*x(i,3)

c m4*m3*Q(1,0,0,x) / x
        q4  = (l34-r31*l14)+m12*(1/r13-r31)/x(i,4)
c m2*m3*Q(1,x/r24,0,0) / (x/r24)
        q1  = (l23-r31*l12)+m12*(1/r13-r31)/x(i,1)
c m2*m3 * Q(0,r13*x/r24,1,0)  / (r13*x/r24)
        q2  = (l23-r31*l12)+m32*(1d0-r31*r13)/x(i,2)
c m4*m3 * Q(0,0,1,r13*x)  / (r13)
        q3  = (l34-r31*l14)+m32*(1d0-r31*r13)/x(i,3)


c m3*m2 qbar(0,0,1,r13*x)/(m2*m3*d)
        qb2  = ((l23-r24*l34)-m42*(r24-1d0/r42)*x(i,3))
c m1*m2 qbar(1,0,0,x)/(m2*m3*d)
        qb1  = ((l12-r24*l14)-m42*(r24-1d0/r42)*x(i,4))

        qb1od=qb1/d
        qb2od=qb2/d



        write(*,*) 'cD0regrm qb2od ',qb2od,pol34/q3,pol23/q2
        write(*,*) 'cD0regrm qb1od ',qb1od,pol12/q1,pol14/q4
        write(*,*) 'cD0regrm q1 ',q1,pol12/qb2od
        write(*,*) 'cD0regrm q4 ',q4,pol14/qb2od

        t12 = abs(p12)/max(abs(m12/x(i,1)),abs(l12),abs(m22*x(i,1)))
        t14 = abs(p14)/
     &        max(abs(m12/x(i,4)),abs(l14),abs(m42*x(i,4)))
        t23 = abs(p23)/
     &        max(abs(m32/(x(i,2))),abs(l23),abs(m22*x(i,2)))
        t34 = abs(p34)/
     &        max(abs(m32/(x(i,3))),abs(l34),abs(m42*x(i,3)))
        t1  = abs(q1)/
     &      max(abs(l23),abs(m12*r31/x(i,1)),abs(r31*l12),
     &      abs(m12/r13/x(i,1)))
        t4  = abs(q4)/
     &      max(abs(l34),abs(m12*r31/x(i,4)),abs(r31*l14),
     &      abs(m12/r13/x(i,4)))
        t2  = abs(q2)/max(abs(l23),abs(r31*l12),
     &                     abs(m32/x(i,2)),abs(r31*r13/x(i,2)))
        t3  = abs(q3)/max(abs(l34),abs(r31*l14),
     &                     abs(m32/x(i,3)),abs(m32*r31*r13/x(i,3)))

        td= abs(d)/
     &      max(abs(l23),abs(r31*l12),abs(l34*r24),abs(r31*l14*r24))
        tqb1=abs(qb1od*d)/max(abs(l12),abs(r24*l14),
     &      abs(m42*r24*x(i,4)),abs(m42/r42)*x(i,4))
        tqb2=abs(qb2od*d)/max(abs(l23),abs(r24*l34),
     &      abs(m42*r24*x(i,3)),abs(m42/r42)*x(i,3))

        

        if (min(t12,t1).gt.min(t14,t4).and.min(t12,t1).gt.tb2od) 
     &      qb2od =pol12/q1
        if (min(t14,t4).gt.min(t12,t1).and.min(t14,t4).gt.tb2od) 
     &      qb2od =pol14/q4
        if (min(t23,t2).gt.min(t34,t3).and.min(t23,t2).gt.tb1od) 
     &      qb1od =pol23/q2
        if (min(t34,t3).gt.min(t23,t2).and.min(t34,t3).gt.tb1od) 
     &      qb1od =pol34/q3
        if (min(t12,tb2od).gt.t1) q1=pol12/qb2od
        if (min(t14,tb2od).gt.t4) q4=pol14/qb2od

        write(*,*) 'cD0regrm qb2od ',qb2od,pol34/q3,pol23/q2
        write(*,*) 'cD0regrm qb1od ',qb1od,pol12/q1,pol14/q4
        write(*,*) 'cD0regrm q1 ',q1,pol12/qb2od
        write(*,*) 'cD0regrm q4 ',q4,pol14/qb2od

        eta1=eta2s_coli(-x(i,1),r31,-ix(i,1),ir13,-ix(i,2))

        if(eta1.ne.cd0)then
           l1=cln_coli(q1,-1d0)

           D0regrm_coli = D0regrm_coli - (2*i-3) * eta1*l1 

           write(*,*) 'D0regrm_coli eta1 i ',eta1*l1
        endif

        eta2=eta2s_coli(-x(i,4),r31,-ix(i,4),ir13,-ix(i,3))

        if(eta2.ne.cd0)then
           l2=cln_coli(q4,-1d0) + cln_coli(r24,ir24)

           D0regrm_coli = D0regrm_coli + (2*i-3) * eta2*l2 

           write(*,*) 'D0regrm_coli eta1 i ',eta2*l2
        endif

        eta3=eta2s_coli(-x(i,4),1d0/r24,-ix(i,4),-ir24,-ix(i,1))

        if(eta3.ne.cd0)then
           iqbod1 = -1d0 ! pure guess
           l3=cln_coli(qb1od,iqbod1) + cln_coli(x(i,1),ix(i,1)) 
     &         - log(m12)

#ifdef CHECK
          if(conscheck)then
            if(dimag(r24).ne.0d0) then
              v=dimag(x(i,4))/dimag(r24)
              u=dimag(x(i,4)/r24)/dimag(1d0/r24)
              l3t=cln_coli(m12+l12*v+l14*u+m22*v*v
     &            +m42*u*u+l24*u*v,-1d0)
     &            -cln_coli(m12*(1d0/r13-r31)+(l23-r31*l12)*v
     &            +(l34-r31*l14)*u,-1d0)
              
              write(*,*) 'D0regrm_coli comp1 l3= ',l3,l3t,l3t/l3

            elseif(dimag(x(i,4)).ne.0d0) then
              l2=cln_coli(q1b,-1d0)-cln_coli(d,-1d0)

              write(*,*) 'D0regrm_coli comp2 l2= ',l2,l1,l2/l1

            elseif(dreal(x(i,4)).gt.0d0) then
              l2=cln_coli(pol14,-1d0)-cln_coli(q4,-1d0)

              write(*,*) 'D0regrm_coli comp2 l2= ',l2,l1,l2/l1

            elseif(dreal(x(i,1)).gt.0d0)then
                l2=cln_coli(pol12,-1d0)-cln_coli(q1,-1d0)

              write(*,*) 'D0regrm_coli comp3 l2= ',l2,l1,l2/l1

c@ continue here

            else
              if (flag3) then
                write(*,*) 
                write(*,*) 'D0regrm_coli: case for l2 not implemented'
                flag3 = .false.
              end if
              write(errout,*) 
              write(errout,*) 
     &            'D0regrm_coli: case for l2 not implemented '
              write(errout,*) 
     &            'D0regrm_coli: args ',p12,p23,p34,p14,p13,p24,
     &            m12,m22,m32,m42
              write(errout,*) 'r24 = ',r24
              write(errout,*) 'x(i,4) = ',x(i,4)
              
              write(errout,*) 'x(i,1) = ',x(i,1)
            endif 
            
c              write(*,*) 'D0regrm_coli l2= ',l2,l1,l2/l1

            if(abs(l2/l1-1d0).gt.testacc) then
              if (flag2) then
                write(*,*) 
                write(*,*) 'D0regrm_coli: eta terms possibly wrong'
                flag2 = .false.
              end if
              write(errout,*) 
              write(errout,*) 'D0regrm_coli: eta terms possibly wrong'
              write(errout,*) 
     &            'D0regrm_coli args ',p12,p23,p34,p14,p13,p24,
     &            m12,m22,m32,m42
              write(errout,*) 
     &            'D0regrm_coli chi ',ch0,ch1,ch2,ch3
              write(errout,*) 'l1 = ',l1
              write(errout,*) 'l2 = ',l2
              write(errout,*) 'diff   = ',abs((l2-l1)/l1)
            endif 
          endif
#endif

           D0regrm_coli = D0regrm_coli + (2*i-3) * eta3*l3 

           write(*,*) 'D0regrm_coli eta1 i ',eta3*l3
        endif




        eta4=eta2s_coli(-x(i,4),1d0/r24,-ix(i,4),-ir24,-ix(i,1))

        if(eta4.ne.0d0)then
          ch0   = m42*(r31*x(i,4))**2+l34*r31*x(i,4)+m32
          test0 = abs(ch0)/max(abs(m42*(r31*x(i,4))**2),
     &        abs(l34*r31*x(i,4)),abs(m32))
          ch1   = m42*(x(i,4))**2+l14*x(i,4)+m12
          test1 = abs(ch1)/
     &        max(abs(m42*(x(i,4))**2),abs(l14*x(i,4)),abs(m12))
          test01=min(test0,test1)
          ch2   = m22*(r31*x(i,1))**2+l23*r31*x(i,1)+m32
          test2 = abs(ch2)/max(abs(m22*(r31*x(i,1))**2),
     &        abs(l23*r31*x(i,1)),abs(m32))
          ch3   = m22*(x(i,1))**2+l12*x(i,1)+m12
          test3 = abs(ch3)/
     &        max(abs(m22*(x(i,1))**2),abs(l12*x(i,1)),abs(m12))
          test23 = min(test2,test3)
          if(dreal(x(i,4)).gt.0d0)then
            ipop(i)=dreal(ch0-ch1)
          elseif(dreal(x(i,1)).gt.0d0)then
            ipop(i)=dreal(ch2-ch3)
          endif
          if (test23.gt.test01) then
            l1 = cln_coli(ch2/ch3,ipop(i))
          else
            l1 = cln_coli(ch0/ch1,ipop(i))      
          end if


#ifdef CHECK
          if(conscheck)then
            if(dimag(r24).ne.0d0) then
              v=dimag(x(i,4))/dimag(r24)
              u=dimag(x(i,4)/r24)/dimag(1d0/r24)
              l2=cln_coli(m32+l23*r31*v+l34*r31*u+m22*(r31*v)**2
     &            +m42*(r31*u)**2+l24*r31*r31*u*v,-1d0)
     &            -cln_coli(m12+l12*v+l14*u+m22*v*v
     &            +m42*u*u+l24*u*v,-1d0)
              
              write(*,*) 'D0regrm_coli comp1 l2= ',l2,l1,l2/l1
              
            elseif(dimag(x(i,4)).ne.0d0.or.dreal(x(i,4)).gt.0d0) then
              if(test01.gt.test23)then
                l2=cln_coli(ch0,-1d0)-cln_coli(ch1,-1d0)
              elseif(test0.gt.test1)then
                l2=cln_coli(ch0,-1d0)-cln_coli(ch0*ch3/ch2,-1d0)
              else
                l2=cln_coli(ch1*ch2/ch3,-1d0)-cln_coli(ch1,-1d0)
              endif

              write(*,*) 'D0regrm_coli comp2 l2= ',l2,l1,l2/l1

            elseif(dreal(x(i,1)).gt.0d0)then
              if(test23.gt.test01)then
                l2=cln_coli(ch2,-1d0)-cln_coli(ch3,-1d0)
              elseif(test2.gt.test3)then
                l2=cln_coli(ch2,-1d0)-cln_coli(ch2*ch1/ch0,-1d0)
              else
                l2=cln_coli(ch3*ch0/ch1,-1d0)-cln_coli(ch3,-1d0)
              endif

              write(*,*) 'D0regrm_coli comp3 l2= ',l2,l1,l2/l1

            else
              if (flag3) then
                write(*,*) 
                write(*,*) 'D0regrm_coli: case for l2 not implemented'
                flag3 = .false.
              end if
              write(errout,*) 
              write(errout,*) 
     &            'D0regrm_coli: case for l2 not implemented '
              write(errout,*) 
     &            'D0regrm_coli: args ',p12,p23,p34,p14,p13,p24,
     &            m12,m22,m32,m42
              write(errout,*) 'r24 = ',r24
              write(errout,*) 'x(i,4) = ',x(i,4)
              
              write(errout,*) 'x(i,1) = ',x(i,1)
            endif 
            
c              write(*,*) 'D0regrm_coli l2= ',l2,l1,l2/l1

            if(abs(l2/l1-1d0).gt.testacc) then
              if (flag2) then
                write(*,*) 
                write(*,*) 'D0regrm_coli: eta terms possibly wrong'
                flag2 = .false.
              end if
              write(errout,*) 
              write(errout,*) 'D0regrm_coli: eta terms possibly wrong'
              write(errout,*) 
     &            'D0regrm_coli args ',p12,p23,p34,p14,p13,p24,
     &            m12,m22,m32,m42
              write(errout,*) 
     &            'D0regrm_coli chi ',ch0,ch1,ch2,ch3
              write(errout,*) 'l1 = ',l1
              write(errout,*) 'l2 = ',l2
              write(errout,*) 'diff   = ',abs((l2-l1)/l1)
            endif 
          endif
#endif


        if(eta.ne.0d0)then
          D0regrm_coli = D0regrm_coli + (2*i-3) * (
     &        + eta*(l1-log(m32/m12)) ) 
        endif
            
          write(*,*) 'D0regrm_coli eta i ',
     &        eta*(l1-log(m32/m12))
          
        endif
        
        D0regrm_coli = D0regrm_coli + (2*i-3) * (
     &      - cspcos_coli(-x(i,3),r43,-ix(i,3),ir34)
     &      - cspcos_coli(-x(i,3),1d0/r34,-ix(i,3),-ir34)
     &      + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14)
     &      + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14)
     &      - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12)
     &      - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12)
     &      + cspcos_coli(-x(i,2),r23,-ix(i,2),ir23)
     &      + cspcos_coli(-x(i,2),1d0/r32,-ix(i,2),-ir23))
        
        write(*,*) 'D0regrm_coli arg1 i ',
     &      -x(i,3),r43,-ix(i,4),ir34
        write(*,*) 'D0regrm_coli arg2 i ',
     &      -x(i,3),1d0/r34,-ix(i,4),-ir34
        write(*,*) 'D0regrm_coli arg3 i ',
     &      -x(i,4),r41,-ix(i,4),ir14
        write(*,*) 'D0regrm_coli arg4 i ',
     &      -x(i,4),1d0/r14,-ix(i,4),-ir14
        write(*,*) 'D0regrm_coli arg5 i ',
     &      -x(i,1),r21,-ix(i,1),ir12
        write(*,*) 'D0regrm_coli arg6 i ',
     &      -x(i,1),1d0/r12,-ix(i,1),-ir12
        write(*,*) 'D0regrm_coli arg7 i ',
     &      -x(i,2),r23,-ix(i,2),ir23
        write(*,*) 'D0regrm_coli arg8 i ',
     &      -x(i,2),1d0/r32,-ix(i,2),-ir23
        
        write(*,*) 'D0regrm_coli contr i ',
     &      - cspcos_coli(-x(i,3),r43,-ix(i,3),ir34),
     &      - cspcos_coli(-x(i,3),1d0/r34,-ix(i,3),-ir34),
     &      + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14),
     &      + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14),
     &      - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12),
     &      - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12),
     &      + cspcos_coli(-x(i,2),r23,-ix(i,2),ir23),
     &      + cspcos_coli(-x(i,2),1d0/r32,-ix(i,2),-ir23)
        write(*,*) 'D0regrm_coli sum i ',
     &      - cspcos_coli(-x(i,3),r43,-ix(i,3),ir34)
     &      - cspcos_coli(-x(i,3),1d0/r34,-ix(i,3),-ir34)
     &      + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14)
     &      + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14)
     &      - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12)
     &      - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12)
     &      + cspcos_coli(-x(i,2),r23,-ix(i,2),ir23)
     &      + cspcos_coli(-x(i,2),1d0/r32,-ix(i,2),-ir23)
        
      enddo
      
      write(*,*) 'D0regrm_coli sum ',D0regrm_coli
      
      D0regrm_coli = D0regrm_coli/det
      
      write(*,*) 'D0regrm_coli=   ',D0regrm_coli

#ifdef CHECK

      if (conscheck.and.(
     &    dimag(x(1,4)).eq.0d0.or.dimag(x(2,4)).eq.0d0)) then

        do j=1,3
          if(j.eq.1) then
            ix(1,4) = +sign(1d0,dreal(d))
            ix(2,4) = -sign(1d0,dreal(d))
          elseif(j.eq.2) then
            ix(1,4) = +sign(1d0,dreal(d))
            ix(2,4) = +sign(1d0,dreal(d))
          elseif(j.eq.3) then
            ix(1,4) = -sign(1d0,dreal(d))
            ix(2,4) = -sign(1d0,dreal(d))
          endif
          ix(1,1) =  sign(1d0,ix(1,4)*dreal(r24)) ! needed for correct
          ix(2,1) =  sign(1d0,ix(2,4)*dreal(r24)) ! continuation


          D0regrm_check = dcmplx(0d0)
          do i=1,2
            eta=eta2s_coli(-x(i,4),1d0/r24,-ix(i,4),-ir24,-ix(i,1))
            
            if(eta.ne.0d0)then
              ch0   = m42*(r31*x(i,4))**2+l34*r31*x(i,4)+m32
              test0 = abs(ch0)/max(abs(m42*(r31*x(i,4))**2),
     &            abs(l34*r31*x(i,4)),abs(m32))
              ch1   = m42*(x(i,4))**2+l14*x(i,4)+m12
              test1 = abs(ch1)/
     &            max(abs(m42*(x(i,4))**2),abs(l14*x(i,4)),abs(m12))
              test01=min(test0,test1)
              ch2   = m22*(r31*x(i,1))**2+l23*r31*x(i,1)+m32
              test2 = abs(ch2)/max(abs(m22*(r31*x(i,1))**2),
     &            abs(l23*r31*x(i,1)),abs(m32))
              ch3   = m22*(x(i,1))**2+l12*x(i,1)+m12
              test3 = abs(ch3)/
     &            max(abs(m22*(x(i,1))**2),abs(l12*x(i,1)),abs(m12))
              test23 = min(test2,test3)
              if(dreal(x(i,4)).gt.0d0)then
                ipop(i)=dreal(ch0-ch1)
              elseif(dreal(x(i,1)).gt.0d0)then
                ipop(i)=dreal(ch2-ch3)
              endif
              if (test23.gt.test01) then
                l1 = cln_coli(ch2/ch3,ipop(i))
              else
                l1 = cln_coli(ch0/ch1,ipop(i))      
              end if
              
              if(conscheck)then
                if(dimag(r24).ne.0d0) then
                  v=dimag(x(i,4))/dimag(r24)
                  u=dimag(x(i,4)/r24)/dimag(1d0/r24)
                  l2=cln_coli(m32+l23*r31*v+l34*r31*u+m22*(r31*v)**2
     &                +m42*(r31*u)**2+l24*r31*r31*u*v,-1d0)
     &                -cln_coli(m12+l12*v+l14*u+m22*v*v
     &                +m42*u*u+l24*u*v,-1d0)
                  
                  write(*,*) 'D0regrm_check comp1 l2= ',l2,l1,l2/l1
                  
                elseif(dimag(x(i,4)).ne.0d0.or.dreal(x(i,4)).gt.0d0)then
                  if(test01.gt.test23)then
                    l2=cln_coli(ch0,-1d0)-cln_coli(ch1,-1d0)
                  elseif(test0.gt.test1)then
                    l2=cln_coli(ch0,-1d0)-cln_coli(ch0*ch3/ch2,-1d0)
                  else
                    l2=cln_coli(ch1*ch2/ch3,-1d0)-cln_coli(ch1,-1d0)
                  endif
                  
                  write(*,*) 'D0regrm_check comp2 l2= ',l2,l1,l2/l1
                  
                elseif(dreal(x(i,1)).gt.0d0)then
                  if(test23.gt.test01)then
                    l2=cln_coli(ch2,-1d0)-cln_coli(ch3,-1d0)
                  elseif(test2.gt.test3)then
                    l2=cln_coli(ch2,-1d0)-cln_coli(ch2*ch1/ch0,-1d0)
                  else
                    l2=cln_coli(ch3*ch0/ch1,-1d0)-cln_coli(ch3,-1d0)
                  endif
                  
                  write(*,*) 'D0regrm_check comp3 l2= ',l2,l1,l2/l1
                  
                else
                  if (flag3) then
                    write(*,*) 
                    write(*,*) '
     &                  D0regrm_check: case for l2 not implemented'
                    flag3 = .false.
                  end if
                  write(errout,*) 
                  write(errout,*) 
     &                'D0regrm_check: case for l2 not implemented '
                  write(errout,*) 
     &                'D0regrm_check: args ',p12,p23,p34,p14,p13,p24,
     &                m12,m22,m32,m42
                  write(errout,*) 'r24 = ',r24
                  write(errout,*) 'x(i,4) = ',x(i,4)
                  
                  write(errout,*) 'x(i,1) = ',x(i,1)
                endif 
                
c     write(*,*) 'D0regrm_check l2= ',l2,l1,l2/l1
                
                if(abs(l2/l1-1d0).gt.testacc) then
                  if (flag2) then
                    write(*,*) 
                    write(*,*) 'D0regrm_check: eta terms possibly wrong'
                    flag2 = .false.
                  end if
                  write(errout,*) 
                  write(errout,*) 
     &                'D0regrm_check: eta terms possibly wrong'
                  write(errout,*) 
     &                'D0regrm_check args ',p12,p23,p34,p14,p13,p24,
     &                m12,m22,m32,m42
                  write(errout,*) 
     &                'D0regrm_check chi ',ch0,ch1,ch2,ch3
                  write(errout,*) 'l1 = ',l1
                  write(errout,*) 'l2 = ',l2
                  write(errout,*) 'diff   = ',abs((l2-l1)/l1)
                endif 
              endif
              
              
              if(eta.ne.0d0)then
                D0regrm_check = D0regrm_check + (2*i-3) * (
     &              + eta*(l1-log(m32/m12)) ) 
              endif 
              
              write(*,*) 'D0regrm_check eta i ',
     &            eta*(l1-log(m32/m12))
              
            endif
            
            D0regrm_check = D0regrm_check + (2*i-3) * (
     &          - cspcos_coli(-x(i,4),r43*r31,-ix(i,4),ir34)
     &          - cspcos_coli(-x(i,4),r31/r34,-ix(i,4),-ir34)
     &          + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14)
     &          + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14)
     &          - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12)
     &          - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12)
     &          + cspcos_coli(-x(i,1),r31*r23,-ix(i,1),ir23)
     &          + cspcos_coli(-x(i,1),r31/r32,-ix(i,1),-ir23))
            
            write(*,*) 'D0regrm_check arg1 i ',
     &          -x(i,4),r43*r31,-ix(i,4),ir34
            write(*,*) 'D0regrm_check arg2 i ',
     &          -x(i,4),r31/r34,-ix(i,4),-ir34
            write(*,*) 'D0regrm_check arg3 i ',
     &          -x(i,4),r41,-ix(i,4),ir14
            write(*,*) 'D0regrm_check arg4 i ',
     &          -x(i,4),1d0/r14,-ix(i,4),-ir14
            write(*,*) 'D0regrm_check arg5 i ',
     &          -x(i,1),r21,-ix(i,1),ir12
            write(*,*) 'D0regrm_check arg6 i ',
     &          -x(i,1),1d0/r12,-ix(i,1),-ir12
            write(*,*) 'D0regrm_check arg7 i ',
     &          -x(i,1),r31*r23,-ix(i,1),ir23
            write(*,*) 'D0regrm_check arg8 i ',
     &          -x(i,1),r31/r32,-ix(i,1),-ir23
            
            write(*,*) 'D0regrm_check contr i ',
     &          - cspcos_coli(-x(i,4),r43*r31,-ix(i,4),ir34),
     &          - cspcos_coli(-x(i,4),r31/r34,-ix(i,4),-ir34),
     &          + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14),
     &          + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14),
     &          - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12),
     &          - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12),
     &          + cspcos_coli(-x(i,1),r31*r23,-ix(i,1),ir23),
     &          + cspcos_coli(-x(i,1),r31/r32,-ix(i,1),-ir23)
            write(*,*) 'D0regrm_check sum i ',
     &          - cspcos_coli(-x(i,4),r43*r31,-ix(i,4),ir34)
     &          - cspcos_coli(-x(i,4),r31/r34,-ix(i,4),-ir34)
     &          + cspcos_coli(-x(i,4),r41,-ix(i,4),ir14)
     &          + cspcos_coli(-x(i,4),1d0/r14,-ix(i,4),-ir14)
     &          - cspcos_coli(-x(i,1),r21,-ix(i,1),ir12)
     &          - cspcos_coli(-x(i,1),1d0/r12,-ix(i,1),-ir12)
     &          + cspcos_coli(-x(i,1),r31*r23,-ix(i,1),ir23)
     &          + cspcos_coli(-x(i,1),r31/r32,-ix(i,1),-ir23)
            
          enddo
          
          write(*,*) 'D0regrm_check sum ',D0regrm_check
          
          D0regrm_check = D0regrm_check/det
          
          write(*,*) 'D0regrm_check =  ',D0regrm_check

          if(cdabs(D0regrm_coli/D0regrm_check-1d0).gt.testacc) then
            write(testout,*) ' D0regrm_coli inconsistent ',j
            write(testout,111) ' D0regrm_coli: p12 = ',p12
            write(testout,111) ' D0regrm_coli: p23 = ',p23
            write(testout,111) ' D0regrm_coli: p34 = ',p34
            write(testout,111) ' D0regrm_coli: p14 = ',p14
            write(testout,111) ' D0regrm_coli: p24 = ',p24
            write(testout,111) ' D0regrm_coli: p13 = ',p13
            write(testout,111) ' D0regrm_coli: m12 = ',m12
            write(testout,111) ' D0regrm_coli: m22 = ',m22
            write(testout,111) ' D0regrm_coli: m32 = ',m32
            write(testout,111) ' D0regrm_coli: m42 = ',m42
            write(testout,111) ' D0regrm_coli  = ',D0regrm_coli
            write(testout,111) ' D0regrm_check = ',D0regrm_check
            write(testout,101) ' diff        = ',
     &          cdabs(D0regrm_coli/D0regrm_check-1d0)
            
            write(*,*)   ' D0regrm_coli inconsistent ',i
            write(*,111) ' D0regrm_coli: p12 = ',p12
            write(*,111) ' D0regrm_coli: p23 = ',p23
            write(*,111) ' D0regrm_coli: p34 = ',p34
            write(*,111) ' D0regrm_coli: p14 = ',p14
            write(*,111) ' D0regrm_coli: p24 = ',p24
            write(*,111) ' D0regrm_coli: p13 = ',p13
            write(*,111) ' D0regrm_coli: m12 = ',m12
            write(*,111) ' D0regrm_coli: m22 = ',m22
            write(*,111) ' D0regrm_coli: m32 = ',m32
            write(*,111) ' D0regrm_coli: m42 = ',m42
            write(*,111) ' D0regrm_coli  = ',D0regrm_coli
            write(*,111) ' D0regrm_check = ',D0regrm_check
            write(*,101) ' diff        = ',
     &          cdabs(D0regrm_coli/D0regrm_check-1d0)
            stop
          endif
          
        enddo
      endif
#endif

      end



************************************************************************
      function cD0gen(q12,q23,q34,q14,q13,q24,m12,m22,m32,m42)
************************************************************************
*  general scalar 4-point function                                     *
*  regular case  based on general result of                            *
*        A.Denner, U.Nierste and R.Scharf, Nucl. Phys. B367 (1991) 637 *
*                                                                      *
*                     m22                                              *
*       q12  ---------------------  q23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*              m12| 1     3 | m32                                      *
*                 |         |                                          *
*                 |    4    |                                          *
*       q14  ---------------------  q34                                *
*                     m42                                              *
*                                                                      *
*  imaginary part sometimes wrong  ????????                            *
*----------------------------------------------------------------------*
*  28.05.02 Ansgar Denner       last changed 16.08.08 Ansgar Denner    *
************************************************************************
      implicit   none
      complex*16 q12,q23,q34,q14,q13,q24
      complex*16 m12,m22,m32,m42
      complex*16 k12,k13,k14,k23,k24,k34
      complex*16 m1,m2,m3,m4
      real*8     ir12,ir14,ir23,ir24,ir34,ir13,ir13or24
      real*8     ix(2,4),is(4),gamma(2)
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 a,b,c,d,det
      real*8     rd
      complex*16 x(2,4),s(4)
      complex*16 q1,q2,q3,q4,p12,p34,p14,p23,l1,l2,l3,l4
      real*8     t1,t2,t3,t4,t12,t34,t14,t23
      complex*16 cD0gen,cspcoe,olietae,ettile,cdln,cD0m0,cD0p0
      real*8     swap
      integer    k,j
      complex*16 cd0
      data       cd0 /(0d0,0d0)/
 
c      write(*,*) 'cD0gen in ',q12,q23,q34,q14,q13,q24,m12,m22,m32,m42

      if (m12.eq.0d0) then
        cD0gen = cD0m0(q23,q12,q14,q34,q13,q24,m32,m22,cd0,m42)
      else if (m22.eq.0d0) then
        cD0gen = cD0m0(q14,q12,q23,q34,q24,q13,m42,m12,cd0,m32)
      else if (m32.eq.0d0) then
        cD0gen = cD0m0(q12,q23,q34,q14,q13,q24,m12,m22,cd0,m42)
      else if (m42.eq.0d0) then
        cD0gen = cD0m0(q23,q34,q14,q12,q24,q13,m22,m32,cd0,m12)
      else if (q12.eq.0d0.and.dimag(m12/m22).eq.0d0) then
        cD0gen = cD0p0(q34,q23,q12,q14,q24,q13,m42,m32,m22,m12)
      else if (q23.eq.0d0.and.dimag(m32/m22).eq.0d0) then
        cD0gen = cD0p0(q14,q12,q23,q34,q24,q13,m42,m12,m22,m32)
      else if (q34.eq.0d0.and.dimag(m32/m42).eq.0d0) then
        cD0gen = cD0p0(q12,q23,q34,q14,q13,q24,m12,m22,m32,m42)
      else if (q14.eq.0d0.and.dimag(m12/m42).eq.0d0) then
        cD0gen = cD0p0(q23,q34,q14,q12,q24,q13,m22,m32,m42,m12)
      else if (q24.eq.0d0.and.dimag(m42/m22).eq.0d0) then
        cD0gen = cD0p0(q13,q34,q24,q12,q14,q23,m12,m32,m42,m22)
      else if (q13.eq.0d0.and.dimag(m12/m32).eq.0d0) then
        cD0gen = cD0p0(q24,q12,q13,q34,q14,q23,m42,m22,m12,m32)
      else

        if(dimag(m12).ne.0d0.or.dimag(m22).ne.0d0.or.
     &       dimag(m32).ne.0d0.or.dimag(m42).ne.0d0) then
          write(*,*) 'D0 for 4 nonzero masses and finite width'
          write(*,*) 'not implemented'
          write(*,*) 'cD0reg: q12 = ',q12
          write(*,*) 'cD0reg: q23 = ',q23
          write(*,*) 'cD0reg: q34 = ',q34
          write(*,*) 'cD0reg: q14 = ',q14
          write(*,*) 'cD0reg: q24 = ',q24
          write(*,*) 'cD0reg: q13 = ',q13
          write(*,*) 'cD0reg: m12 = ',m12
          write(*,*) 'cD0reg: m22 = ',m22
          write(*,*) 'cD0reg: m32 = ',m32
          write(*,*) 'cD0reg: m42 = ',m42
          stop
        end if

      m1 = sqrt(m12) 
      m2 = sqrt(m22) 
      m3 = sqrt(m32)
      m4 = sqrt(m42) 
      k12 = (m12+m22-q12)/(m1*m2)
      k13 = (m12+m32-q13)/(m1*m3)
      k14 = (m12+m42-q14)/(m1*m4)
      k23 = (m22+m32-q23)/(m2*m3)
      k24 = (m22+m42-q24)/(m2*m4)
      k34 = (m32+m42-q34)/(m3*m4)
      
      if (abs(k13).lt.2d0) then
         if (abs(k24).ge.2d0) then
            swap = k13
            k13  = k24
            k24  = swap
            swap = k14 
            k14  = k23
            k23  = swap
         else if (abs(k12).ge.2d0) then 
            swap = k13
            k13  = k12
            k12  = swap
            swap = k24 
            k24  = k34
            k34  = swap
         else if (abs(k14).ge.2d0) then 
            swap = k13
            k13  = k14
            k14  = swap
            swap = k24 
            k24  = k23
            k23  = swap
         else if (abs(k23).ge.2d0) then 
            swap = k13
            k13  = k23
            k23  = swap
            swap = k24 
            k24  = k14
            k14  = swap
         else if (abs(k34).ge.2d0) then 
            swap = k13
            k13  = k34
            k34  = swap
            swap = k24 
            k24  = k12
            k12  = swap
         end if
      end if

      if (k12.ne.0d0) then
        r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2)))
      else
        r12 = dcmplx(0d0,1d0)
      end if
      if (k14.ne.0d0) then
        r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2)))
      else
        r14 = dcmplx(0d0,1d0)
      end if
      if (k23.ne.0d0) then
        r23 = k23/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k23**2)))
      else
        r23 = dcmplx(0d0,1d0)
      end if
      if (k34.ne.0d0) then
        r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2)))
      else
        r34 = dcmplx(0d0,1d0)
      end if
      if (k13.ne.0d0) then
        r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2)))
      else
        r13 = dcmplx(0d0,1d0)
      end if
c      r13 = 1d0/r13
      if (k24.ne.0d0) then
        r24 = k24/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k24**2)))
      else
        r24 = dcmplx(0d0,1d0)
      end if
c      r24 = 1d0/r24
      a   =  k34/r24-k23 + (k12-k14/r24)*r13
      b   =  (1d0/r13-r13)*(1d0/r24-r24)+k12*k34-k14*k23
      c   =  k34*r24-k23 + (k12-k14*r24)/r13
      d   =  k23-r13*k12 - (k34-r13*k14)*r24
      det = sqrt(b*b-4d0*a*c)
      x(1,4) = (-b+det)/(2d0*a)
      x(2,4) = (-b-det)/(2d0*a)

c      write(*,*) 'cD0gen a ',a,k34,r24,k23,k12,k14,r13

c>      x(1,4) = (-b+sqrt(b*b-4d0*a*c))/(2d0*a)
c>      x(2,4) = (-b-sqrt(b*b-4d0*a*c))/(2d0*a)
      if(abs(x(1,4)).gt.abs(x(2,4))) then
        x(2,4) = c/(a*x(1,4))
      else
        x(1,4) = c/(a*x(2,4))
      end if

      if(dreal(k12).lt.-2d0) then
        ir12 = sign(1d1,1d0-abs(r12))
      else
        ir12 = 0d0
      end if
      if(dreal(k14).lt.-2d0) then
        ir14 = sign(1d1,1d0-abs(r14))
      else
        ir14 = 0d0
      end if
      if(dreal(k23).lt.-2d0) then
        ir23 = sign(1d1,1d0-abs(r23))
      else
        ir23 = 0d0
      end if
      if(dreal(k34).lt.-2d0) then
        ir34 = sign(1d1,1d0-abs(r34))
      else
        ir34 = 0d0
      end if
c      if(k24.le.-2d0) then
        ir24 = sign(1d1,1d0-abs(r24))
c      else
c        ir24 = 0d0
c      end if
c      if(k13.le.-2d0) then
        ir13 = sign(1d1,1d0-abs(r13))
c      else
c        ir13 = 0d0
c      end if
      ir13or24 = sign(1d1,dreal(r24)*ir13-dreal(r13)*ir24)
      

       rd = dreal(d)
c imaginary parts of x(i,4) also needed for x(i,4) < 0 in ettile!!
c      if(dreal(x(1,4)).gt.0d0) then
         ix(1,4) = sign(1d0,-rd)
c      else
c         ix(1,4) = 0d0
c      end if
c      if(dreal(x(2,4)).gt.0d0) then
         ix(2,4) = sign(1d0,rd)
c      else
c         ix(2,4) = 0d0
c      end if
c --<

c>      write(*,*) 'cD0gen q ',q12,q23,q34
c>      write(*,*) 'cD0gen q ',q14,q13,q24
c>      write(*,*) 'cD0gen m ',m12,m22,m32,m42
c>      write(*,*) 'cD0gen k ',k12,k13,k14
c>      write(*,*) 'cD0gen k ',k23,k24,k34
c>      write(*,*) 'cD0gen r ',r12,r13,r14
c>      write(*,*) 'cD0gen r ',r23,r24,r34
c>      write(*,*) 'cD0gen ir ',ir12,ir13,ir14
c>      write(*,*) 'cD0gen ir ',ir23,ir24,ir34
c>      write(*,*) 'cD0gen a ',a,b,c,d
c>      write(*,*) 'cD0gen x ',x(1,4),x(2,4)
c>      write(*,*) 'cD0gen ix ',ix(1,4),ix(2,4)

      x(1,1) = x(1,4)/r24
      x(2,1) = x(2,4)/r24
      x(1,2) = x(1,4)/r24*r13
      x(2,2) = x(2,4)/r24*r13
      x(1,3) = x(1,4)*r13
      x(2,3) = x(2,4)*r13
      s(1)  = r12
      s(2)  = r23
      s(3)  = r34
      s(4)  = r14

      is(1)  = ir12
      is(2)  = ir23
      is(3)  = ir34
      is(4)  = ir14

      if(dreal(x(1,1)).gt.0d0) then
         ix(1,1) = sign(1d0,-rd/dreal(r24))
      else
         ix(1,1) = 0d0
      end if
      if(dreal(x(2,1)).gt.0d0) then
         ix(2,1) = sign(1d0,rd/dreal(r24))
      else
         ix(2,1) = 0d0
      end if
      if(dreal(x(1,3)).gt.0d0) then
         ix(1,3) = sign(1d0,-rd*dreal(r13))
      else
         ix(1,3) = 0d0
      end if
      if(dreal(x(2,3)).gt.0d0) then
         ix(2,3) = sign(1d0,rd*dreal(r13))
      else
         ix(2,3) = 0d0
      end if
      if(dreal(x(1,2)).gt.0d0) then
         ix(1,2) = sign(1d0,-rd*dreal(r13/r24))
      else
         ix(1,2) = 0d0
      end if
      if(dreal(x(2,2)).gt.0d0) then
         ix(2,2) = sign(1d0,rd*dreal(r13/r24))
      else
         ix(2,2) = 0d0
      end if
c      gamma(1) = dsign(1d0,dreal(a*(x(1,4)-x(2,4))))
      gamma(1) = dsign(1d0,dreal(det))
      gamma(2) = -gamma(1)

      cD0gen = dcmplx(0d0)
c      do 20 k=2,1,-1
      do 20 k=1,2
        do 10 j=1,4
          cD0gen = cD0gen + (-1)**(j+k) * (
     &           cspcoe(-x(k,j),s(j),-ix(k,j),is(j))
     &         + cspcoe(-x(k,j),1d0/s(j),-ix(k,j),-is(j)) )

c>          write(*,*) 'cD0gen ',k,j,x(k,j),s(j),ix(k,j),is(j)
c>          write(*,*) 'cD0gen ',
c>     &         cspcoe(-x(k,j),s(j),-ix(k,j),is(j)) 
c>     &         , cspcoe(-x(k,j),1d0/s(j),-ix(k,j),-is(j)) 
c>          write(*,*) 'cD0gen ', (-1)**(j+k) * (
c>     &         cspcoe(-x(k,j),s(j),-ix(k,j),is(j)) 
c>     &         + cspcoe(-x(k,j),1d0/s(j),-ix(k,j),-is(j)) )
c>          write(*,*) 'cD0gen= ',k,j, cD0gen

10      continue

c        write(*,*) 'cD0gen sp = ',cD0gen

c *** improves numerical stability in extreme cases

        p12 = r24/x(k,4)+k12+x(k,4)/r24
        p14 = 1/x(k,4)+k14+x(k,4)
        p23 = r24/(x(k,4)*r13)+k23+x(k,4)*r13/r24
        p34 = 1/(r13*x(k,4))+k34+x(k,4)*r13
        q1  = (k34-r13*k14)+(1/r13-r13)/x(k,4)
        q2  = k23-r13*k12+(1/r13-r13)*r24/x(k,4)
        q3  = ((k23-r24*k34)-(r24-1d0/r24)*r13*x(k,4))/d
        q4  = ((k12-r24*k14)-(r24-1d0/r24)*x(k,4))/d

c        write(*,*) 'cD0gen ',q1,p34/q3,p14/q4
c        write(*,*) 'cD0gen ',q2,p12/q4,p23/q3
c        write(*,*) 'cD0gen ',q3,p23/q2,p34/q1
c        write(*,*) 'cD0gen ',q4,p12/q2,p14/q1
c        write(*,*) 'cD0gen ',p12,q4*q2
c        write(*,*) 'cD0gen ',p14,q4*q1
c        write(*,*) 'cD0gen ',p23,q3*q2
c        write(*,*) 'cD0gen ',p34,q3*q1

        t12 = abs(p12)/max(abs(r24/x(k,4)),abs(k12),abs(x(k,4)/r24))
        t14 = abs(p14)/
     &        max(abs(1/x(k,4)),abs(k14),abs(x(k,4)))
        t23 = abs(p23)/
     &        max(abs(r24/(r13*x(k,4))),abs(k23),abs(x(k,4)*r13/r24))
        t34 = abs(p34)/
     &        max(abs(1/(r13*x(k,4))),abs(k34),abs(x(k,4)*r13))
        t1  = abs(q1)/
     &      max(abs(k34*x(k,4)),abs(r13),abs(r13*k14*x(k,4)),abs(1/r13))
        t2  = abs(q2)/max(abs(k23),abs(r13*k12),
     &                     abs(r24*r13/x(k,4)),abs(r24/r13/x(k,4)))
        t3  = abs(q3)/max(abs(k23),abs(r24*k34),
     &                     abs(r24*r13*x(k,4)),abs(r13/r24*x(k,4)))
        t4  = abs(q4)/
     &      max(abs(k12/x(k,4)),abs(r24),abs(r24*k14/x(k,4)),abs(1/r24))

        if (min(t12,t2).gt.min(t14,t1).and.min(t12,t2).gt.t4) q4 =p12/q2
        if (min(t14,t1).gt.min(t12,t2).and.min(t14,t1).gt.t4) q4 =p14/q1
        if (min(t23,t2).gt.min(t34,t1).and.min(t23,t2).gt.t3) q3 =p23/q2
        if (min(t34,t1).gt.min(t23,t2).and.min(t34,t1).gt.t3) q3 =p34/q1
        if (min(t34,t3).gt.min(t14,t4).and.min(t34,t3).gt.t1) q1 =p34/q3
        if (min(t14,t4).gt.min(t34,t3).and.min(t14,t4).gt.t1) q1 =p14/q4
        if (min(t12,t4).gt.min(t23,t3).and.min(t12,t4).gt.t2) q2 =p12/q4
        if (min(t23,t3).gt.min(t12,t4).and.min(t23,t3).gt.t1) q2 =p23/q3
        if (min(t4,t2).gt.t12) p12=q4*q2
        if (min(t4,t1).gt.t14) p14=q4*q1
        if (min(t3,t2).gt.t23) p23=q3*q2
        if (min(t3,t1).gt.t34) p34=q3*q1


c ***
      
      if (dimag(r13).eq.0d0) then

        l1 = cdln(q1,-1d0)
        l3 = cdln(q3,sign(1d0,gamma(k)*dreal(r13)*ir24))
        l2 = cdln(q2,-1d0)
        l4 = cdln(q4,sign(1d0,gamma(k)*ir24))

c        l1 = cdln((k13-2*r13)/x(k,4)+k34-r13*k14,-1d0)
c        l3 = cdln(((k24-2*r24)*r13*x(k,4)+k23-r24*k34)/d,
c     &            dsign(1d0,gamma(k)*dreal(r13)*ir24))
c        l2 = cdln((k13-2*r13)*r24/x(k,4)+k23-r13*k12,-1d0)
c        l4 = cdln(((k24-2*r24)*x(k,4)+k12-r24*k14)/d,
c     &            dsign(1d0,gamma(k)*ir24))
        
        cD0gen = cD0gen - (-1)**k *( 
     &            ettile(-x(k,4),r13,-ix(k,4),ir13) *
     &            (  cdln(r13*x(k,4),sign(1d0,dreal(r13))*ix(k,4)) 
     &               + l1 + l3)
     &          + ettile(-x(k,4),1/r24,-ix(k,4),-ir24) *
     &            (  cdln(x(k,4)/r24,sign(1d0,dreal(r24))*ix(k,4)) 
     &               + l2 + l4)
     &          - (ettile(-x(k,4),r13/r24,-ix(k,4),ir13or24) 
     &             + olietae(r13,1/r24,ir13,-ir24,ir13or24)) *
     &            (cdln(x(k,4)*r13/r24,sign(1d0,dreal(r13/r24))*ix(k,4)) 
     &               + l3 + l2)
     &          + ettile(-x(k,4),-r13/r24,-ix(k,4),-ir13or24) 
     &            * olietae(r13,1/r24,ir13,-ir24,ir13or24) 
     &             )

c>        write(*,*) 'cD0gen a ',
c>     &            ettile(-x(k,4),r13,-ix(k,4),ir13) *
c>     &            (  cdln(r13*x(k,4),sign(1d0,dreal(r13))*ix(k,4)) 
c>     &               + l1 + l3)
c>        write(*,*) 'cD0gen b ',
c>     &          + ettile(-x(k,4),1/r24,-ix(k,4),-ir24) *
c>     &            (  cdln(x(k,4)/r24,sign(1d0,dreal(r24))*ix(k,4)) 
c>     &               + l2 + l4)
c>        write(*,*) 'cD0gen b ',
c>     &          + ettile(-x(k,4),1/r24,-ix(k,4),-ir24) ,
c>     &              cdln(x(k,4)/r24,sign(1d0,dreal(r24))*ix(k,4)) 
c>     &              , + l2, + l4
c>        write(*,*) 'cD0gen c ',
c>     &          - (ettile(-x(k,4),r13/r24,-ix(k,4),ir13or24) 
c>     &             + olietae(r13,1/r24,ir13,-ir24,ir13or24)) *
c>     &            (cdln(x(k,4)*r13/r24,sign(1d0,dreal(r13/r24))*ix(k,4)) 
c>     &               + l3 + l2)
c>        write(*,*) 'cD0gen d ',
c>     &          + ettile(-x(k,4),-r13/r24,-ix(k,4),-ir13or24) 
c>     &            * olietae(r13,1/r24,ir13,-ir24,ir13or24) 
c>
c>        write(*,*) 'cD0gen real r= ',cD0gen
       
      else

c branch not yet fully tested!!
        l1 = cdln(p34,-sign(1d0,gamma(k)*dimag(r13*x(k,4)*b)))
        l2 = cdln(p12,-sign(1d0,gamma(k)*dimag(x(k,4)/r24*b)))
        l3 = cdln(p23,-sign(1d0,gamma(k)*dimag(r13*x(k,4)/r24*b)))
c        l1 = cdln(1d0/(r13*x(k,4))+k34+r13*x(k,4),
c     &             -dsign(1d0,gamma(k)*dimag(r13*x(k,4)*b)))
c        l2 = cdln(r24/x(k,4)+k12+x(k,4)/r24,
c     &             -dsign(1d0,gamma(k)*dimag(x(k,4)/r24*b)))
c        l3 = cdln(r24/(r13*x(k,4))+k23+r13*x(k,4)/r24,
c     &             -dsign(1d0,gamma(k)*dimag(r13*x(k,4)/r24*b)))
c        write(*,*) 'cD0gen l ',l1,l2,l3
        cD0gen = cD0gen - (-1)**k *( 
     &            olietae(-x(k,4),r13,-ix(k,4),ir13,0d0) *
     &            (  cdln(r13*x(k,4),sign(1d0,dreal(r13))*ix(k,4)) 
     &               + l1)
     &          + olietae(-x(k,4),1/r24,-ix(k,4),-ir24,0d0) *
     &            (  cdln(x(k,4)/r24,sign(1d0,ix(k,4)/dreal(r24))) 
     &               + l2)
     &          - (olietae(-x(k,4),r13/r24,-ix(k,4),ir13-ir24,0d0) 
     &             + olietae(r13,1/r24,ir13,-ir24,ir13-ir24)) *
     &           (cdln(x(k,4)*r13/r24,sign(1d0,dreal(r13/r24))*ix(k,4)) 
     &               + l3 + l2)
c--> changed 21.06.05
c     &          + (olietae(-x(k,4),-r13/r24,-ix(k,4),-ir13+ir24,0d0) 
c     &             * olietae(r13,1/r24,ir13,-ir24,ir13-ir24)) *
c     &               (1-gamma(k)*sign(1d0,dreal(b)))
     &             )
        if ((1-gamma(k)*sign(1d0,dreal(b))).ne.0.and.
     &       olietae(r13,1/r24,ir13,-ir24,ir13-ir24).ne.0d0) then  
        cD0gen = cD0gen - (-1)**k *( 
     &          + (olietae(-x(k,4),-r13/r24,-ix(k,4),-ir13+ir24,0d0) 
     &             * olietae(r13,1/r24,ir13,-ir24,ir13-ir24)) *
     &               (1-gamma(k)*sign(1d0,dreal(b)))
     &             )
        end if
c <-- changed 21.06.05

c>         write(*,*) 'cD0gen eta a',
c>     &            olietae(-x(k,4),r13,-ix(k,4),ir13,0d0) *
c>     &            (  cdln(r13*x(k,4),sign(1d0,dreal(r13))*ix(k,4)) 
c>     &               + l1)
c>         write(*,*) 'cD0gen eta b',
c>     &          + olietae(-x(k,4),1/r24,-ix(k,4),-ir24,0d0) *
c>     &            (  cdln(x(k,4)/r24,sign(1d0,ix(k,4)/dreal(r24))) 
c>     &               + l2)
c>         write(*,*) 'cD0gen eta c',
c>     &          - (olietae(-x(k,4),r13/r24,-ix(k,4),ir13-ir24,0d0) 
c>     &             + olietae(r13,1/r24,ir13,-ir24,ir13-ir24)) *
c>     &           (cdln(x(k,4)*r13/r24,sign(1d0,dreal(r13/r24))*ix(k,4)) 
c>     &               + l3 + l2)
c>        if ((1-gamma(k)*sign(1d0,dreal(b))).ne.0.and.
c>     &      olietae(r13,1/r24,ir13,-ir24,ir13-ir24).ne.0d0) then  
c>         write(*,*) 'cD0gen eta d',
c>     &          + (olietae(-x(k,4),-r13/r24,-ix(k,4),-ir13+ir24,0d0) 
c>     &             * olietae(r13,1/r24,ir13,-ir24,ir13-ir24)) *
c>     &               (1-gamma(k)*sign(1d0,dreal(b)))
c>        end if
c>        write(*,*) 'cD0gen imr ',cD0gen
       
        end if

20    continue

c      write(*,*) 'cD0gen sum ',cD0gen

      cD0gen = cD0gen/(m1*m2*m3*m4*sqrt(b*b-4d0*a*c))

c      write(*,*) 'cD0gen =   ',cD0gen    

      end if
      end

************************************************************************
      function cD0p0(q12i,q23i,q34i,q14i,q13i,q24i,m12i,m22i,m32i,m42i)
************************************************************************
*  general scalar 4-point function                                     *
*  case with one zero external invariant                               *
*  based on equation (23) of                                           *
*        A.Denner, U.Nierste and R.Scharf, Nucl. Phys. B367 (1991) 637 *
*                                                                      *
*                     m22                                              *
*       q12  ---------------------  q23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*              m12| 1     3 | m32                                      *
*                 |         |                                          *
*                 |    4    |                                          *
*       q14  ---------------------  q34 = 0                            *
*                     m42                                              *
*                                                                      *
*----------------------------------------------------------------------*
*  24.01.07 Ansgar Denner       last changed 16.06.08 Ansgar Denner    *
************************************************************************
      implicit   none
      complex*16 q12,q23,q34,q14,q13,q24
      complex*16 m12,m22,m32,m42
      complex*16 q12i,q23i,q34i,q14i,q13i,q24i
      complex*16 m12i,m22i,m32i,m42i
      complex*16 k12,k13,k14,k23,k24,k34
      complex*16 l12,l13,l14,l23,l24,l34
      complex*16 m1,m2,m3,m4
      real*8     ir12,ir14,ir23,ir34,ir13
      real*8     iy(2),iz(2)
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 r21,r31,r41,r32,r42,r43
      complex*16 a,b,c,d,det
      real*8     rd
      complex*16 x(2),y(2),z(2),eta,v
      complex*16 p12,p34,p14,p23
      real*8     t12,t34,t14,t23
      complex*16 cD0p0,cspcoe,olietae,cdln,cD0m0,cD0p0o
      integer    k
      complex*16 cd0
      data       cd0 /(0d0,0d0)/

      integer    ltest
      common /ltest/ ltest
      integer    flag
      data       flag /0/
      save       flag
      integer    errflag
      save       errflag
      data       errflag /0/
      integer    errout
      data       errout /90/
      real*8     testaccd0
      common /testaccd0/ testaccd0
 
c      write(*,*) 'cD0gen in ',q12,q23,q34,q14,q13,q24,m12,m22,m32,m42

      if (flag.eq.0) then
        write(*,*) 
        write(*,*) 'cD0p0 not yet tested in practice'
        flag = 1
      end if 

      if (m12i.eq.0d0) then
        cD0p0=cD0m0(q23i,q12i,q14i,q34i,q13i,q24i,m32i,m22i,cd0,m42i)
      else if (m22i.eq.0d0) then
        cD0p0=cD0m0(q14i,q12i,q23i,q34i,q24i,q13i,m42i,m12i,cd0,m32i)
      else if (m32i.eq.0d0) then
        cD0p0=cD0m0(q12i,q23i,q34i,q14i,q13i,q24i,m12i,m22i,cd0,m42i)
      else if (m42i.eq.0d0) then
        cD0p0=cD0m0(q23i,q34i,q14i,q12i,q24i,q13i,m22i,m32i,cd0,m12i)
      else

        if(dimag(m42i/m32i).ne.0d0.or.q34i.ne.0d0) then
          write(*,*) 'cD0p0 called improperly:'
          write(*,*) 'cD0p0: q12 = ',q12i
          write(*,*) 'cD0p0: q23 = ',q23i
          write(*,*) 'cD0p0: q34 = ',q34i
          write(*,*) 'cD0p0: q14 = ',q14i
          write(*,*) 'cD0p0: q24 = ',q24i
          write(*,*) 'cD0p0: q13 = ',q13i
          write(*,*) 'cD0p0: m12 = ',m12i
          write(*,*) 'cD0p0: m22 = ',m22i
          write(*,*) 'cD0p0: m32 = ',m32i
          write(*,*) 'cD0p0: m42 = ',m42i
          stop
        end if

        q13=q12i
        q23=q23i
        q24=q34i
        q14=q14i
        q12=q13i
        q34=q24i
        m12=m12i
        m32=m22i
        m22=m32i
        m42=m42i

        l12 = (m12+m22-q12)
        l13 = (m12+m32-q13)
        l14 = (m12+m42-q14)
        l23 = (m22+m32-q23)
        l24 = (m22+m42-q24)
        l34 = (m32+m42-q34)      
        
        if (l12.ne.0d0) then
          r12 = l12/(2d0*m22)*(1d0+sqrt(1d0-4d0*m12*m22/l12**2))
        else
          r12 = sqrt(m12/m22)*dcmplx(0d0,1d0)
        end if
        r21= r12*m22/m12
        if (l14.ne.0d0) then
          r14 = l14/(2d0*m42)*(1d0+sqrt(1d0-4d0*m12*m42/l14**2))
        else
          r14 = sqrt(m12/m42)*dcmplx(0d0,1d0)
        end if
        r41 = r14*m42/m12
        if (l34.ne.0d0) then
          r34 = l34/(2d0*m42)*(1d0+sqrt(1d0-4d0*m32*m42/l34**2))
        else
          r34 = sqrt(m12/m42)*dcmplx(0d0,1d0)
        end if
        r43 = r34*m42/m32
        if (l23.ne.0d0) then
          r23 = l23/(2d0*m32)*(1d0+sqrt(1d0-4d0*m22*m32/l23**2))
        else
          r23 = sqrt(m22/m32)*dcmplx(0d0,1d0)
        end if
        r32 = r23*m32/m22
        if (l13.ne.0d0) then
          r13 = l13/(2d0*m32)*(1d0+sqrt(1d0-4d0*m12*m32/l13**2))
        else
          r13 = sqrt(m12/m32)*dcmplx(0d0,1d0)
        end if
        r31 = r13*m32/m12
c       r13 = 1d0/r13
        r24=1d0
        r42=m42/m22
        a   =  m42*((l34/r42-l23) + (l12-l14/r42)*r31)
        b   =  m32*(1d0/r31-r13)*m42*(1d0/r42-r24)+l12*l34-l14*l23
        c   =  m12*((l34*r24-l23) + (l12-l14*r24)/r13)
        d   =  m12*(l23-r31*l12 - (l34-r31*l14)*r24)
        det = sqrt(b*b-4d0*a*c)
        y(1) = (-b+det)/(2d0*a)
        y(2) = (-b-det)/(2d0*a)
        if(abs(y(1)).gt.abs(y(2))) then
          y(2) = c/(a*y(1))
        else
          y(1) = c/(a*y(2))
        end if

        if(dreal(l12).lt.0d0) then
          ir12 = sign(1d1,1d0-abs(r12*r12*m22/m12))
        else
          ir12 = 0d0
        end if
        if(dreal(l14).lt.0d0) then
          ir14 = sign(1d1,1d0-abs(r14*r14*m42/m12))
        else
          ir14 = 0d0
        end if
        if(dreal(l23).lt.0d0) then
          ir23 = sign(1d1,1d0-abs(r23*r23*m32/m22))
        else
          ir23 = 0d0
        end if
        if(dreal(l34).lt.0d0) then
          ir34 = sign(1d1,1d0-abs(r34*r34*m42/m32))
        else
          ir34 = 0d0
        end if
c       if(k13.le.-2d0) then
        ir13 = sign(1d1,1d0-abs(r13*r13*m32/m12))
c       else
c         ir13 = 0d0
c       end if

c>        write(*,*) 'cD0p0 q ',q12,q23,q34
c>        write(*,*) 'cD0p0 q ',q14,q13,q24
c>        write(*,*) 'cD0p0 m ',m12,m22,m32,m42
c>        write(*,*) 'cD0p0 l ',l12,l13,l14
c>        write(*,*) 'cD0p0 l ',l23,l24,l34
c>        write(*,*) 'cD0p0 r ',r12,r13,r14
c>        write(*,*) 'cD0p0 r ',r21,r31,r41
c>        write(*,*) 'cD0p0 r ',r23,r24,r34
c>        write(*,*) 'cD0p0 r ',r32,r42,r43
c>        write(*,*) 'cD0p0 ir',ir12,ir13,ir14
c>        write(*,*) 'cD0p0 ir',ir23,ir24,ir34
c>        write(*,*) 'cD0p0 a ',a*m42*m2*m3,b*m1*m2*m3*m4,c*m12*m2*m3,
c>     &      d*m12*m2*m3
c>        write(*,*) 'cD0p0 a ',a,b,c,d
c>        write(*,*) 'cD0p0 x ',x(1),x(2)
c>        write(*,*) 'cD0p0 ix',ix(1),ix(2)
        
        z(1) = y(1)*r31
        z(2) = y(2)*r31
        
c     result should be independent of imaginary parts of y and z!!!!
        rd = dreal(d)
        if(dreal(y(1)).gt.0d0) then
          iy(1) = sign(1d0,-rd)
        else
          iy(1) = 0d0
        end if
        if(dreal(y(2)).gt.0d0) then
          iy(2) = sign(1d0,rd)
        else
          iy(2) = 0d0
        end if
        if(dreal(z(1)).gt.0d0) then
          iz(1) = sign(1d0,-rd*dreal(r13))
        else
          iz(1) = 0d0
        end if
        if(dreal(z(2)).gt.0d0) then
          iz(2) = -sign(1d0,rd*dreal(r13))
        else
          iz(2) = 0d0
        end if
        
c>        write(*,*) 'cD0p0 y ',y(1),y(2)
c>        write(*,*) 'cD0p0 iy',iy(1),iy(2)
c>        write(*,*) 'cD0p0 z',z(1),z(2)
c>        write(*,*) 'cD0p0 iz',iz(1),iz(2)
        
        cD0p0 = dcmplx(0d0)
c       do k=2,1,-1
        do k=1,2
          cD0p0 = cD0p0 - (-1)**k * (
     &        cspcoe(-y(k),1d0/r12,-iy(k),-ir12)
     &        +cspcoe(-y(k),r21,-iy(k),ir12)
     &        -cspcoe(-y(k),1d0/r14,-iy(k),-ir14)
     &        -cspcoe(-y(k),r41,-iy(k),ir14)
     &        +cspcoe(-z(k),1d0/r34,-iz(k),-ir34)
     &        +cspcoe(-z(k),r43,-iz(k),ir34)
     &        -cspcoe(-z(k),1d0/r32,-iz(k),-ir23)
     &        -cspcoe(-z(k),r23,-iz(k),ir23)
     &        )          
          
c>          write(*,*) 'cD0p0 args yz',
c>     &        -y(k),-iy(k),-z(k),-iz(k)
c>          write(*,*) 'cD0p0 args r',
c>     &        1d0/r12,r21,1d0/(r14),r41,
c>     &        1d0/r34,r43,1d0/(r32),r23
c>          write(*,*) 'cD0p0 args ',
c>     &        1d0+y(k)/(r12),
c>     &        1d0+y(k)*r21,
c>     &        1d0+y(k)/r14,
c>     &        1d0+y(k)*r41,
c>     &        1d0+z(k)/r34,
c>     &        1d0+z(k)*r43,
c>     &        1d0+z(k)/r32,
c>     &        1d0+z(k)*r23
c>          write(*,*) 'cD0p0 csp ',
c>     &        cspcoe(-y(k),1d0/r12,-iy(k),-ir12),
c>     &        +cspcoe(-y(k),r21,-iy(k),ir12),
c>     &        -cspcoe(-y(k),1d0/r14,-iy(k),-ir14),
c>     &        -cspcoe(-y(k),r41,-iy(k),ir14),
c>     &        +cspcoe(-z(k),1d0/r34,-iz(k),-ir34),
c>     &        +cspcoe(-z(k),r43,-iz(k),ir34),
c>     &        -cspcoe(-z(k),1d0/r32,-iz(k),-ir23),
c>     &        -cspcoe(-z(k),r23,-iz(k),ir23)
c>          
c>          write(*,*) 'cD0p0 cspsum ',cD0p0
                    
c     *** improves numerical stability in extreme cases
          
          eta=olietae(-y(k),-1/z(k),-iy(k),iz(k),-ir13)
          
          if(eta.ne.0d0)then
            p14 = m12+m42*y(k)*y(k)+l14*y(k)
            p12 = m12+m22*y(k)*y(k)+l12*y(k)
            p34 = m32+m42*z(k)*z(k)+l34*z(k)
            p23 = m32+m22*z(k)*z(k)+l23*z(k)
            
            t14=abs(p14)/max(abs(m12),abs(m42*y(k)*y(k)),abs(l14*y(k)))
            t12=abs(p12)/max(abs(m12),abs(m22*y(k)*y(k)),abs(l12*y(k)))
            t34=abs(p34)/max(abs(m32),abs(m42*z(k)*z(k)),abs(l34*z(k)))
            t23=abs(p23)/max(abs(m32),abs(m22*z(k)*z(k)),abs(l23*z(k)))
            
            if (min(t12,t14).gt.min(t34,t23)) then
              if(dimag(p14/p12).ne.0d0)then
                v= cdln(p14/p12,0d0)
              else
                v=cdln(p14,-1d0)-cdln(p12,-1d0)
              endif
            else
              if(dimag(p34/p23).ne.0d0)then
                v= cdln(p34/p23,0d0)
              else
                v=cdln(p34,-1d0)-cdln(p23,-1d0)
              endif
            endif

c>            write(*,*) 'cD0p0 v',p14/p12,p34/p23              
c>            write(*,*) 'cD0p0 v',cdln(p14,-1d0)-cdln(p12,-1d0),
c>     &          cdln(p34,-1d0)-cdln(p23,-1d0)
            
            cD0p0 = cD0p0 
     &          - (-1)**k *eta*v

c>            write(*,*) 'cD0p0 etacon ', - (-1)**k *v,eta
c>            write(*,*) 'cD0p0 etasum ',cD0p0
          endif
        enddo
        
c        write(*,*) 'cD0p0 sum ',cD0p0

        cD0p0 = cD0p0/det

c        write(*,*) 'cD0p0 =   ',cD0p0    
        
c end new version

        if (ltest.lt.10) return

c begin old version

        m1 = sqrt(m12) 
        m2 = sqrt(m22) 
        m3 = sqrt(m32)
        m4 = sqrt(m42) 
        l12 = (m12+m22-q12)
        l13 = (m12+m32-q13)
        l14 = (m12+m42-q14)
        l23 = (m22+m32-q23)
        l24 = (m22+m42-q24)
        l34 = (m32+m42-q34)
        k12 = (m12+m22-q12)/(m1*m2)
        k13 = (m12+m32-q13)/(m1*m3)
        k14 = (m12+m42-q14)/(m1*m4)
        k23 = (m22+m32-q23)/(m2*m3)
        k24 = (m22+m42-q24)/(m2*m4)
        k34 = (m32+m42-q34)/(m3*m4)
        
        
        if (k12.ne.0d0) then
          r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2)))
        else
          r12 = dcmplx(0d0,1d0)
        end if
        if (k14.ne.0d0) then
          r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2)))
        else
          r14 = dcmplx(0d0,1d0)
        end if
        if (k23.ne.0d0) then
          r23 = k23/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k23**2)))
        else
          r23 = dcmplx(0d0,1d0)
        end if
        if (k34.ne.0d0) then
          r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2)))
        else
          r34 = dcmplx(0d0,1d0)
        end if
        if (k13.ne.0d0) then
          r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2)))
        else
          r13 = dcmplx(0d0,1d0)
        end if
c     r13 = 1d0/r13
        r24=m4/m2
        a   =  k34/r24-k23 + (k12-k14/r24)*r13
        b   =  (1d0/r13-r13)*(1d0/r24-r24)+k12*k34-k14*k23
        c   =  k34*r24-k23 + (k12-k14*r24)/r13
        d   =  k23-r13*k12 - (k34-r13*k14)*r24
        det = sqrt(b*b-4d0*a*c)
        x(1) = (-b+det)/(2d0*a)
        x(2) = (-b-det)/(2d0*a)
c>      x(1) = (-b+sqrt(b*b-4d0*a*c))/(2d0*a)
c>      x(2) = (-b-sqrt(b*b-4d0*a*c))/(2d0*a)
        if(abs(x(1)).gt.abs(x(2))) then
          x(2) = c/(a*x(1))
        else
          x(1) = c/(a*x(2))
        end if
        
        if(dreal(k12).lt.-2d0) then
          ir12 = sign(1d1,1d0-abs(r12))
        else
          ir12 = 0d0
        end if
        if(dreal(k14).lt.-2d0) then
          ir14 = sign(1d1,1d0-abs(r14))
        else
          ir14 = 0d0
        end if
        if(dreal(k23).lt.-2d0) then
          ir23 = sign(1d1,1d0-abs(r23))
        else
          ir23 = 0d0
        end if
        if(dreal(k34).lt.-2d0) then
          ir34 = sign(1d1,1d0-abs(r34))
        else
          ir34 = 0d0
        end if
c       if(k13.le.-2d0) then
        ir13 = sign(1d1,1d0-abs(r13))
c       else
c         ir13 = 0d0
c       end if

c>        write(*,*) 'cD0p0o q ',q12,q23,q34
c>        write(*,*) 'cD0p0o q ',q14,q13,q24
c>        write(*,*) 'cD0p0o m ',m12,m22,m32,m42
c>        write(*,*) 'cD0p0o k ',k12,k13,k14
c>        write(*,*) 'cD0p0o k ',k23,k24,k34
c>        write(*,*) 'cD0p0o k ',l12,l13,l14
c>        write(*,*) 'cD0p0o k ',l23,l24,l34
c>        write(*,*) 'cD0p0o r ',r12,r13,r14
c>        write(*,*) 'cD0p0o r ',r23,r24,r34
c>        write(*,*) 'cD0p0o r ',r12*m1/m2,r13*m1/m3,r14*m1/m4
c>        write(*,*) 'cD0p0o r ',r12*m2/m1,r13*m3/m1,r14*m4/m1
c>        write(*,*) 'cD0p0o r ',r23*m2/m3,r24*m2/m4,r34*m3/m4
c>        write(*,*) 'cD0p0o r ',r23*m3/m2,r24*m4/m2,r34*m4/m3
c>        write(*,*) 'cD0p0o ir',ir12,ir13,ir14
c>        write(*,*) 'cD0p0o ir',ir23,ir24,ir34
c>        write(*,*) 'cD0p0o a ',a,b,c,d
c>        write(*,*) 'cD0p0o x ',x(1),x(2)
c>        write(*,*) 'cD0p0o ix',ix(1),ix(2)
        
        y(1) = x(1)*m1/m4
        y(2) = x(2)*m1/m4
        z(1) = x(1)*m3/m4*r13
        z(2) = x(2)*m3/m4*r13
        
c     result should be independent of imaginary parts of y and z!!!!
        rd = dreal(d)
        if(dreal(y(1)).gt.0d0) then
          iy(1) = sign(1d0,-rd)
        else
          iy(1) = 0d0
        end if
        if(dreal(y(2)).gt.0d0) then
          iy(2) = sign(1d0,rd)
        else
          iy(2) = 0d0
        end if
        if(dreal(z(1)).gt.0d0) then
          iz(1) = sign(1d0,-rd*dreal(r13))
        else
          iz(1) = 0d0
        end if
        if(dreal(z(2)).gt.0d0) then
          iz(2) = -sign(1d0,rd*dreal(r13))
        else
          iz(2) = 0d0
        end if
        
c>        write(*,*) 'cD0p0o y ',y(1),y(2)
c>        write(*,*) 'cD0p0o iy',iy(1),iy(2)
c>        write(*,*) 'cD0p0o z',z(1),z(2)
c>        write(*,*) 'cD0p0o iz',iz(1),iz(2)
        
        cD0p0o = dcmplx(0d0)
c     do k=2,1,-1
        do k=1,2
          cD0p0o = cD0p0o - (-1)**k * (
     &        cspcoe(-y(k),m2/(m1*r12),-iy(k),-ir12)
     &        +cspcoe(-y(k),m2/m1*r12,-iy(k),ir12)
     &        -cspcoe(-y(k),m4/(m1*r14),-iy(k),-ir14)
     &        -cspcoe(-y(k),m4/m1*r14,-iy(k),ir14)
     &        +cspcoe(-z(k),m4/(m3*r34),-iz(k),-ir34)
     &        +cspcoe(-z(k),m4/m3*r34,-iz(k),ir34)
     &        -cspcoe(-z(k),m2/(m3*r23),-iz(k),-ir23)
     &        -cspcoe(-z(k),m2/m3*r23,-iz(k),ir23)
     &      )          
          
c>          write(*,*) 'cD0p0o args yz',
c>     &        -y(k),-iy(k),-z(k),-iz(k)
c>          write(*,*) 'cD0p0o args r',
c>     &        m2/(m1*r12),m2/m1*r12,m4/(m1*r14),m4/m1*r14,
c>     &        m4/(m3*r34),m4/m3*r34,m2/(m3*r23),m2/m3*r23
c>          write(*,*) 'cD0p0o args ',
c>     &        1d0+y(k)*m2/(m1*r12),
c>     &        1d0+y(k)*m2/m1*r12,
c>     &        1d0+y(k)*m4/(m1*r14),
c>     &        1d0+y(k)*m4/m1*r14,
c>     &        1d0+z(k)*m4/(m3*r34),
c>     &        1d0+z(k)*m4/m3*r34,
c>     &        1d0+z(k)*m2/(m3*r23),
c>     &        1d0+z(k)*m2/m3*r23
c>          write(*,*) 'cD0p0o csp ',
c>     &        cspcoe(-y(k),m2/(m1*r12),-iy(k),-ir12),
c>     &        +cspcoe(-y(k),m2/m1*r12,-iy(k),ir12),
c>     &        -cspcoe(-y(k),m4/(m1*r14),-iy(k),-ir14),
c>     &        -cspcoe(-y(k),m4/m1*r14,-iy(k),ir14),
c>     &        +cspcoe(-z(k),m4/(m3*r34),-iz(k),-ir34),
c>     &        +cspcoe(-z(k),m4/m3*r34,-iz(k),ir34),
c>     &        -cspcoe(-z(k),m2/(m3*r23),-iz(k),-ir23),
c>     &        -cspcoe(-z(k),m2/m3*r23,-iz(k),ir23)
c>
c>          write(*,*) 'cD0p0o cspsum ',cD0p0o
                    
c     *** improves numerical stability in extreme cases
          
          eta=olietae(-y(k),-1/z(k),-iy(k),iz(k),-ir13)
          
          if(eta.ne.0d0)then
            p14 = m12+m42*y(k)*y(k)+l14*y(k)
            p12 = m12+m22*y(k)*y(k)+l12*y(k)
            p34 = m32+m42*z(k)*z(k)+l34*z(k)
            p23 = m32+m22*z(k)*z(k)+l23*z(k)
            
            t14=abs(p14)/max(abs(m12),abs(m42*y(k)*y(k)),abs(l14*y(k)))
            t12=abs(p12)/max(abs(m12),abs(m22*y(k)*y(k)),abs(l12*y(k)))
            t34=abs(p34)/max(abs(m32),abs(m42*z(k)*z(k)),abs(l34*z(k)))
            t23=abs(p23)/max(abs(m32),abs(m22*z(k)*z(k)),abs(l23*z(k)))
            
            if (min(t12,t14).gt.min(t34,t23)) then
              if(dimag(p14/p12).ne.0d0)then
                v= cdln(p14/p12,0d0)
              else
                v=cdln(p14,-1d0)-cdln(p12,-1d0)
              endif
            else
              if(dimag(p34/p23).ne.0d0)then
                v= cdln(p34/p23,0d0)
              else
                v=cdln(p34,-1d0)-cdln(p23,-1d0)
              endif
            endif
            
c>            write(*,*) 'cD0p0o v',p14/p12,p34/p23              
c>            write(*,*) 'cD0p0o v',cdln(p14,-1d0)-cdln(p12,-1d0),
c>     &          cdln(p34,-1d0)-cdln(p23,-1d0)
            
            cD0p0o = cD0p0o 
     &          - (-1)**k *eta*v

c>            write(*,*) 'cD0p0o etacon ', - (-1)**k *v,eta
c>            write(*,*) 'cD0p0o etasum ',cD0p0o
          endif
        enddo
        
c        write(*,*) 'cD0p0o sum ',cD0p0o

        cD0p0o = cD0p0o/(m1*m2*m3*m4*det)

c        write(*,*) 'cD0p0o =   ',cD0p0o    
        
c end old version
      
        if (abs((cd0p0-cd0p0o)/cd0p0).gt.testaccd0) then
          if (errflag.eq.0) then
            write(*,*) 
            write(*,*) 'cD0p0 possibly wrong '
            errflag = 1
          end if
          write(errout,*) 
          write(errout,*) 'cD0p0 possibly wrong '
          write(errout,*) 'cD0p0 args ',
     &        q12,q23,q34,q14,q13,q24,m12,m22,m32,m42
          write(errout,*) 'cD0p0o = ',cD0p0o
          write(errout,*) 'cD0p0  = ',cD0p0
          write(errout,*) 'diff   = ',abs((cD0p0-cD0p0o)/cD0p0)
        end if
     
      endif
        
      end

************************************************************************
      function cD0msp0(q12i,q23i,q34i,q14i,q13i,q24i,
     &    m12i,m22i,m32i,m42i)
************************************************************************
*  general scalar 4-point function                                     *
*  case with one zero external invariant                               *
*  and two associated small real masses                                *
*  based on equation (23) of                                           *
*        A.Denner, U.Nierste and R.Scharf, Nucl. Phys. B367 (1991) 637 *
*                                                                      *
*                    m22                                               *
*       q12  ---------------------  q23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*              m12| 1     3 | m32 small                                *
*                 |         |                                          *
*                 |    4    |                                          *
*       q14  ---------------------  q34 = 0                            *
*                  m42 small                                           *
*                                                                      *
*----------------------------------------------------------------------*
*  24.01.07 Ansgar Denner       last changed 26.01.07 Ansgar Denner    *
************************************************************************
      implicit   none
      complex*16 q12,q23,q34,q14,q13,q24
      complex*16 m12,m22,m32,m42
      complex*16 q12i,q23i,q34i,q14i,q13i,q24i
      complex*16 m12i,m22i,m32i,m42i
      complex*16 k12,k13,k14,k23,k24,k34
      complex*16 l12,l13,l14,l23,l24,l34
      complex*16 m1,m2,m3,m4
      real*8     ir12,ir14,ir23,ir34,ir13
      real*8     iy(2),iz(2)
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 r21,r31,r41,r32,r42,r43
      complex*16 a,b,c,d
      real*8     rd
      complex*16 x(2),y(2),z(2),eta,v
      complex*16 p12,p34,p14,p23
      real*8     t12,t34,t14,t23
      complex*16 cD0msp0,cspcoe,olietae,cdln,cD0m0,cD0msp0o
      complex*16 cd0
      data       cd0 /(0d0,0d0)/
 
      real*8     pi,pi2_6
      common /pietc/  pi,pi2_6
      real*8     lambda,lambda2
      common /lambda/ lambda,lambda2
      real*8     minfscale,minfscale2,shiftms2
      common /cmsing/ minfscale,minfscale2,shiftms2

      integer    ltest
      common /ltest/ ltest
      integer    flag
      data       flag /0/
      save       flag
      integer    errflag
      save       errflag
      data       errflag /0/
      integer    errout
      data       errout /90/
      real*8     testaccd0
      common /testaccd0/ testaccd0

 2    if (flag.eq.0) then
        write(*,*) 
        write(*,*) 'cD0msp0 not yet tested in practice'
        flag = 1
      end if 

      if (m12i.eq.0d0) then
        cD0msp0 = cD0m0(q23i,q12i,q14i,q34i,q13i,q24i,
     &      m32i,m22i,cd0,m42i)
      else if (m22i.eq.0d0) then
        cD0msp0 = cD0m0(q14i,q12i,q23i,q34i,q24i,q13i,
     &      m42i,m12i,cd0,m32i)
      else if (m32i.eq.0d0) then
        cD0msp0 = cD0m0(q12i,q23i,q34i,q14i,q13i,q24i,
     &      m12i,m22i,cd0,m42i)
      else if (m42i.eq.0d0) then
        cD0msp0 = cD0m0(q23i,q34i,q14i,q12i,q24i,q13i,
     &      m22i,m32i,cd0,m12i)
      else

        if(dimag(m42i).ne.0d0.or.dimag(m32i).ne.0d0.or.
     &       (q34i).ne.0d0) then
          write(*,*) 'cD0msp0 called improperly:'
          write(*,*) 'cD0msp0: q12 = ',q12i
          write(*,*) 'cD0msp0: q23 = ',q23i
          write(*,*) 'cD0msp0: q34 = ',q34i
          write(*,*) 'cD0msp0: q14 = ',q14i
          write(*,*) 'cD0msp0: q24 = ',q24i
          write(*,*) 'cD0msp0: q13 = ',q13i
          write(*,*) 'cD0msp0: m12 = ',m12i
          write(*,*) 'cD0msp0: m22 = ',m22i
          write(*,*) 'cD0msp0: m32 = ',m32i
          write(*,*) 'cD0msp0: m42 = ',m42i
          stop
        end if

      q12=q13i
      q23=q23i
      q34=q24i
      q14=q14i
      q13=q12i
      q24=q34i
      m12=m12i
      m22=m32i*minfscale2
      m32=m22i
      m42=m42i*minfscale2


c new version
 100  continue

      l12 = (m12    -q12)
      l13 = (m12+m32-q13)
      l14 = (m12    -q14)
      l23 = (    m32-q23)
      l24 = (       -q24)
      l34 = (m32    -q34)
      
      r12 = l12/m22
      r21 = l12/m12  
      r14 = l14/m42
      r41 = l14/m12
      r23 = l23/m32
      r32 = l23/m22
      r34 = l34/m42
      r43 = l34/m32
      if (l13.ne.0d0) then
        r13 = l13/(2d0*m32)*(1d0+sqrt(1d0-4d0*m12*m32/l13**2))
      else
        r13 = sqrt(m12/m32)*dcmplx(0d0,1d0)
      end if
      r31 = r13*m32/m12
c      r13 = 1d0/r13
      r24=1d0
      r42=m42/m22
      a   =  m42*(l34/r42-l23 + (l12-l14/r42)*r31)
      b   =  l12*l34-l14*l23
      c   =  m12*(l34*r24-l23 + (l12-l14*r24)/r13)
      d   =  m12*(l23-r31*l12 - (l34-r31*l14)*r24)
      y(1) = -c/b
      y(2) = -b/a

      if(dreal(r12).lt.0d0) then
        ir12 = -1d1
      else
        ir12 = 0d0
      end if
      if(dreal(r14).lt.0d0) then
        ir14 = -1d1
      else
        ir14 = 0d0
      end if
      if(dreal(r23).lt.0d0) then
        ir23 = -1d1
      else
        ir23 = 0d0
      end if
      if(dreal(r34).lt.0d0) then
        ir34 = -1d1
      else
        ir34 = 0d0
      end if
c      if(k13.le.-2d0) then
        ir13 = sign(1d1,1d0-abs(r13*r13*m32/m12))
c      else
c        ir13 = 0d0
c      end if

c>      write(*,*) 'cD0msp0 q ',q12,q23,q34
c>      write(*,*) 'cD0msp0 q ',q14,q13,q24
c>      write(*,*) 'cD0msp0 m ',m12,m22,m32,m42
c>      write(*,*) 'cD0msp0 k ',l12,l13,l14
c>      write(*,*) 'cD0msp0 k ',l23,l24,l34
c>      write(*,*) 'cD0msp0 r ',r12,r13,r14
c>      write(*,*) 'cD0msp0 r ',r21,r31,r41
c>      write(*,*) 'cD0msp0 r ',r23,r24,r34
c>      write(*,*) 'cD0msp0 r ',r32,r42,r43
c>      write(*,*) 'cD0msp0 ir',ir12,ir13,ir14
c>      write(*,*) 'cD0msp0 ir',ir23,ir24,ir34
c>      write(*,*) 'cD0msp0 a ',a,b,c,d

      z(1) = y(1)*r31
      z(2) = y(2)*r31

      rd = dreal(d)
c result should be independent of imaginary parts of y and z!!!!
      if(dreal(y(1)).gt.0d0) then
         iy(1) = sign(1d0,-rd)
      else
         iy(1) = 0d0
      end if
      if(dreal(y(2)).gt.0d0) then
         iy(2) = sign(1d0,rd)
      else
         iy(2) = 0d0
      end if
      if(dreal(z(1)).gt.0d0) then
         iz(1) = sign(1d0,-rd*dreal(r13))
      else
         iz(1) = 0d0
      end if
      if(dreal(z(2)).gt.0d0) then
         iz(2) = -sign(1d0,rd*dreal(r13))
      else
         iz(2) = 0d0
      end if

c>      write(*,*) 'cD0msp0 y ',y(1),y(2)
c>      write(*,*) 'cD0msp0 iy',iy(1),iy(2)
c>      write(*,*) 'cD0msp0 z',z(1),z(2)
c>      write(*,*) 'cD0msp0 iz',iz(1),iz(2)

      cD0msp0 = dcmplx(0d0)
      cD0msp0 = cD0msp0 + (
     &      cspcoe(-y(1),cd0,-iy(1),-ir12)
     &      +cspcoe(-y(1),r21,-iy(1),ir12)
     &      -cspcoe(-y(1),cd0,-iy(1),-ir14)
     &      -cspcoe(-y(1),r41,-iy(1),ir14)
     &      +cspcoe(-z(1),cd0,-iz(1),-ir34)
     &      +cspcoe(-z(1),r43,-iz(1),ir34)
     &      -cspcoe(-z(1),cd0,-iz(1),-ir23)
     &      -cspcoe(-z(1),r23,-iz(1),ir23)
     &      )          

c>      write(*,*) 'cD0msp0 args yz',
c>     &      -y(1),-iy(1),-z(1),-iz(1)
c>      write(*,*) 'cD0msp0 args r',
c>     &      1d0/(r12),r21,1d0/(r14),r41,
c>     &      1d0/(r34),r43,1d0/(r32),r23
c>      write(*,*) 'cD0msp0 args ',
c>     &      1d0+y(1)/(r12),
c>     &      1d0+y(1)*r21,
c>     &      1d0+y(1)/(r14),
c>     &      1d0+y(1)*r41,
c>     &      1d0+z(1)/(r34),
c>     &      1d0+z(1)*r43,
c>     &      1d0+z(1)/(r32),
c>     &      1d0+z(1)*r23
c>      write(*,*) 'cD0msp0 csp ',
c>     &      cspcoe(-y(1),cd0,-iy(1),-ir12),
c>     &      +cspcoe(-y(1),r21,-iy(1),ir12),
c>     &      -cspcoe(-y(1),cd0,-iy(1),-ir14),
c>     &      -cspcoe(-y(1),r41,-iy(1),ir14),
c>     &      +cspcoe(-z(1),cd0,-iz(1),-ir34),
c>     &      +cspcoe(-z(1),r43,-iz(1),ir34),
c>     &      -cspcoe(-z(1),cd0,-iz(1),-ir23),
c>     &      -cspcoe(-z(1),r23,-iz(1),ir23)
c>      write(*,*) 'cD0msp0 cspsum ',cD0msp0


c *** improves numerical stability in extreme cases

        eta=olietae(-y(1),-1/z(1),-iy(1),iz(1),-ir13)

        if(eta.ne.0d0)then
          p14 = m12+l14*y(1)
          p12 = m12+l12*y(1)
          p34 = m32+l34*z(1)
          p23 = m32+l23*z(1)
          
          t14 = abs(p14)/max(abs(m12),abs(l14*y(1)))
          t12 = abs(p12)/max(abs(m12),abs(l12*y(1)))
          t34 = abs(p34)/max(abs(m32),abs(l34*z(1)))
          t23 = abs(p23)/max(abs(m32),abs(l23*z(1)))
          
          if (min(t12,t14).gt.min(t34,t23)) then
            if(dimag(p14/p12).ne.0d0)then
              v= cdln(p14/p12,0d0)
            else
              v=cdln(p14,-1d0)-cdln(p12,-1d0)
            endif
          else
            if(dimag(p34/p23).ne.0d0)then
              v= cdln(p34/p23,0d0)
            else
              v=cdln(p34,-1d0)-cdln(p23,-1d0)
            endif
          endif

c>          write(*,*) 'cD0msp0 v',p14/p12,p34/p23              
c>          write(*,*) 'cD0msp0 v',cdln(p14,-1d0)-cdln(p12,-1d0),
c>     &        cdln(p34,-1d0)-cdln(p23,-1d0)

          cD0msp0 = cD0msp0 
     &        +eta*v

c>          write(*,*) 'cD0msp0 etacon ', - (-1)**k *v,eta
c>          write(*,*) 'cD0msp0 etasum ',cD0msp0
        endif


        cD0msp0 = cD0msp0 - (
     &      cspcoe(-y(2),1d0/r12,-iy(2),-ir12)     
     &      -.5d0*(cdln(-y(2),-iy(2))+cdln(r21,ir12))**2
     &      -cspcoe(cd0,1d0/r21,iy(2),-ir12)
     &      -cspcoe(-y(2),1d0/(r14),-iy(2),-ir14)
     &      +.5d0*(cdln(-y(2),-iy(2))+cdln(r41,ir14))**2
     &      +cspcoe(cd0,1d0/r41,iy(2),-ir14)
     &      +cspcoe(-z(2),1d0/r34,-iz(2),-ir34)
     &      -.5d0*(cdln(-z(2),-iz(2))+cdln(r43,ir34))**2
     &      -cspcoe(cd0,1d0/r43,iz(2),-ir34)
     &      -cspcoe(-z(2),1d0/(r32),-iz(2),-ir23)
     &      +.5d0*(cdln(-z(2),-iz(2))+cdln(r23,ir23))**2
     &      +cspcoe(cd0,1d0/r23,iz(2),-ir23)
     &      ) 

c>      write(*,*) 'cD0msp0 args yz',
c>     &      -y(2),-iy(2),-z(2),-iz(2)
c>      write(*,*) 'cD0msp0 args r',
c>     &      1d0/(r12),r21,1d0/(r14),r41,
c>     &      1d0/(r34),r43,1d0/(r32),r23
c>      write(*,*) 'cD0msp0 args ',
c>     &      1d0+y(1)/(r12),
c>     &      1d0+y(1)*r21,
c>     &      1d0+y(1)/(r14),
c>     &      1d0+y(1)*r41,
c>     &      1d0+z(1)/(r34),
c>     &      1d0+z(1)*r43,
c>     &      1d0+z(1)/(r32),
c>     &      1d0+z(1)*r23
c>      write(*,*) 'cD0msp0 csp ',
c>     &      cspcoe(-y(2),1d0/(r12),-iy(2),-ir12),     
c>     &      -.5d0*(cdln(-y(2),-iy(2))+cdln(r21,ir12))**2
c>     &      -cspcoe(cd0,1d0/(r21),iy(2),-ir12),
c>     &      -cspcoe(-y(2),1d0/(r14),-iy(2),-ir14),
c>     &      +.5d0*(cdln(-y(2),-iy(2))+cdln(r41,ir14))**2
c>     &      +cspcoe(cd0,1d0/(r41),iy(2),-ir14),
c>     &      +cspcoe(-z(2),1d0/(r34),-iz(2),-ir34),
c>     &      -.5d0*(cdln(-z(2),-iz(2))+cdln(r43,ir34))**2
c>     &      -cspcoe(cd0,1d0/(r43),iz(2),-ir34),
c>     &      -cspcoe(-z(2),1d0/(r32),-iz(2),-ir23),
c>     &      +.5d0*(cdln(-z(2),-iz(2))+cdln(r23,ir23))**2
c>     &      +cspcoe(cd0,1d0/(r23),iz(2),-ir23)
c>
c>      write(*,*) 'cD0msp0 cspsum ',cD0msp0


c *** improves numerical stability in extreme cases

        eta=olietae(-y(2),-1/z(2),-iy(2),iz(2),-ir13)

        if(eta.ne.0d0)then
          p14 = m42*y(2)*y(2)+l14*y(2)
          p12 = m22*y(2)*y(2)+l12*y(2)
          p34 = m42*z(2)*z(2)+l34*z(2)
          p23 = m22*z(2)*z(2)+l23*z(2)
          
          t14 = abs(p14)/max(abs(m42*y(2)*y(2)),abs(l14*y(2)))
          t12 = abs(p12)/max(abs(m22*y(2)*y(2)),abs(l12*y(2)))
          t34 = abs(p34)/max(abs(m42*z(2)*z(2)),abs(l34*z(2)))
          t23 = abs(p23)/max(abs(m22*z(2)*z(2)),abs(l23*z(2)))
          
          if (min(t12,t14).gt.min(t34,t23)) then
            if(dimag(p14/p12).ne.0d0)then
              v= cdln(p14/p12,0d0)
            else
              v=cdln(p14,-1d0)-cdln(p12,-1d0)
            endif
          else
            if(dimag(p34/p23).ne.0d0)then
              v= cdln(p34/p23,0d0)
            else
              v=cdln(p34,-1d0)-cdln(p23,-1d0)
            endif
          endif

c>          write(*,*) 'cD0msp0 v',p14/p12,p34/p23              
c>          write(*,*) 'cD0msp0 v',cdln(p14,-1d0)-cdln(p12,-1d0),
c>     &        cdln(p34,-1d0)-cdln(p23,-1d0)

          cD0msp0 = cD0msp0 
     &        - eta*v

c>          write(*,*) 'cD0msp0 etacon ', - (-1)**k *v,eta
c>          write(*,*) 'cD0msp0 etasum ',cD0msp0
        endif

c      write(*,*) 'cD0msp0 sum ',cD0msp0

      cD0msp0 = cD0msp0/b

c      write(*,*) 'cD0msp0 =   ',cD0msp0    

c end new version

      if (ltest.lt.10) return

c begin old version

      m1 = sqrt(m12) 
      m2 = sqrt(m22) 
      m3 = sqrt(m32)
      m4 = sqrt(m42) 
      l12 = (m12    -q12)
      l13 = (m12+m32-q13)
      l14 = (m12    -q14)
      l23 = (    m32-q23)
      l24 = (       -q24)
      l34 = (m32    -q34)
      k12 = (m12    -q12)/(m1*m2)
      k13 = (m12+m32-q13)/(m1*m3)
      k14 = (m12    -q14)/(m1*m4)
      k23 = (   +m32-q23)/(m2*m3)
      k24 = (       -q24)/(m2*m4)
      k34 = (m32    -q34)/(m3*m4)
      

      r12 = k12
      r14 = k14
      r23 = k23
      r34 = k34
      if (k13.ne.0d0) then
        r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2)))
      else
        r13 = dcmplx(0d0,1d0)
      end if
c      r13 = 1d0/r13
      r24=m4/m2
      a   =  k34/r24-k23 + (k12-k14/r24)*r13
      b   =  k12*k34-k14*k23
      c   =  k34*r24-k23 + (k12-k14*r24)/r13
      d   =  k23-r13*k12 - (k34-r13*k14)*r24
      x(1) = -c/b
      x(2) = -b/a

      if(dreal(r12).lt.0d0) then
        ir12 = -1d1
      else
        ir12 = 0d0
      end if
      if(dreal(r14).lt.0d0) then
        ir14 = -1d1
      else
        ir14 = 0d0
      end if
      if(dreal(r23).lt.0d0) then
        ir23 = -1d1
      else
        ir23 = 0d0
      end if
      if(dreal(r34).lt.0d0) then
        ir34 = -1d1
      else
        ir34 = 0d0
      end if
c      if(k13.le.-2d0) then
        ir13 = sign(1d1,1d0-abs(r13))
c      else
c        ir13 = 0d0
c      end if

c>      write(*,*) 'cD0msp0o q ',q12,q23,q34
c>      write(*,*) 'cD0msp0o q ',q14,q13,q24
c>      write(*,*) 'cD0msp0o m ',m12,m22,m32,m42
c>      write(*,*) 'cD0msp0o k ',k12,k13,k14
c>      write(*,*) 'cD0msp0o k ',k23,k24,k34
c>      write(*,*) 'cD0msp0o r ',r12,r13,r14
c>      write(*,*) 'cD0msp0o r ',r23,r24,r34
c>      write(*,*) 'cD0msp0o r ',r12*m1/m2,r13*m1/m3,r14*m1/m4
c>      write(*,*) 'cD0msp0o r ',r12*m2/m1,r13*m3/m1,r14*m4/m1
c>      write(*,*) 'cD0msp0o r ',r23*m2/m3,r24*m2/m4,r34*m3/m4
c>      write(*,*) 'cD0msp0o r ',r23*m3/m2,r24*m4/m2,r34*m4/m3
c>      write(*,*) 'cD0msp0o ir',ir12,ir13,ir14
c>      write(*,*) 'cD0msp0o ir',ir23,ir24,ir34
c>      write(*,*) 'cD0msp0o a ',a,b,c,d
c>      write(*,*) 'cD0msp0o a ',a*m42*m2*m3,b*m1*m2*m3*m4,c*m12*m2*m3,
c>     &    d*m12*m2*m3
c>      write(*,*) 'cD0msp0o x ',x(1),x(2)
c>      write(*,*) 'cD0msp0o ix',ix(1),ix(2)

      y(1) = x(1)*m1/m4
      y(2) = x(2)*m1/m4
      z(1) = x(1)*m3/m4*r13
      z(2) = x(2)*m3/m4*r13

c result should be independent of imaginary parts of y and z!!!!
      if(dreal(y(1)).gt.0d0) then
         iy(1) = sign(1d0,-rd)
      else
         iy(1) = 0d0
      end if
      if(dreal(y(2)).gt.0d0) then
         iy(2) = sign(1d0,rd)
      else
         iy(2) = 0d0
      end if
      if(dreal(z(1)).gt.0d0) then
         iz(1) = sign(1d0,-rd*dreal(r13))
      else
         iz(1) = 0d0
      end if
      if(dreal(z(2)).gt.0d0) then
         iz(2) = -sign(1d0,rd*dreal(r13))
      else
         iz(2) = 0d0
      end if

c>      write(*,*) 'cD0msp0o y ',y(1),y(2)
c>      write(*,*) 'cD0msp0o iy',iy(1),iy(2)
c>      write(*,*) 'cD0msp0o z',z(1),z(2)
c>      write(*,*) 'cD0msp0o iz',iz(1),iz(2)

      cD0msp0o = dcmplx(0d0)
      cD0msp0o = cD0msp0o + (
     &      cspcoe(-y(1),cd0,-iy(1),-ir12)
     &      +cspcoe(-y(1),m2/m1*r12,-iy(1),ir12)
     &      -cspcoe(-y(1),cd0,-iy(1),-ir14)
     &      -cspcoe(-y(1),m4/m1*r14,-iy(1),ir14)
     &      +cspcoe(-z(1),cd0,-iz(1),-ir34)
     &      +cspcoe(-z(1),m4/m3*r34,-iz(1),ir34)
     &      -cspcoe(-z(1),cd0,-iz(1),-ir23)
     &      -cspcoe(-z(1),m2/m3*r23,-iz(1),ir23)
     &      )          

c>      write(*,*) 'cD0msp0o args yz',
c>     &      -y(1),-iy(1),-z(1),-iz(1)
c>      write(*,*) 'cD0msp0o args r',
c>     &      m2/(m1*r12),m2/m1*r12,m4/(m1*r14),m4/m1*r14,
c>     &      m4/(m3*r34),m4/m3*r34,m2/(m3*r23),m2/m3*r23
c>      write(*,*) 'cD0msp0o args ',
c>     &      1d0+y(1)*m2/(m1*r12),
c>     &      1d0+y(1)*m2/m1*r12,
c>     &      1d0+y(1)*m4/(m1*r14),
c>     &      1d0+y(1)*m4/m1*r14,
c>     &      1d0+z(1)*m4/(m3*r34),
c>     &      1d0+z(1)*m4/m3*r34,
c>     &      1d0+z(1)*m2/(m3*r23),
c>     &      1d0+z(1)*m2/m3*r23
c>      write(*,*) 'cD0msp0o csp ',
c>     &      cspcoe(-y(1),cd0,-iy(1),-ir12),
c>     &      +cspcoe(-y(1),m2/m1*r12,-iy(1),ir12),
c>     &      -cspcoe(-y(1),cd0,-iy(1),-ir14),
c>     &      -cspcoe(-y(1),m4/m1*r14,-iy(1),ir14),
c>     &      +cspcoe(-z(1),cd0,-iz(1),-ir34),
c>     &      +cspcoe(-z(1),m4/m3*r34,-iz(1),ir34),
c>     &      -cspcoe(-z(1),cd0,-iz(1),-ir23),
c>     &      -cspcoe(-z(1),m2/m3*r23,-iz(1),ir23)
c>
c>      write(*,*) 'cD0msp0o cspsum ',cD0msp0o


c *** improves numerical stability in extreme cases

        eta=olietae(-y(1),-1/z(1),-iy(1),iz(1),-ir13)

        if(eta.ne.0d0)then
          p14 = m12+l14*y(1)
          p12 = m12+l12*y(1)
          p34 = m32+l34*z(1)
          p23 = m32+l23*z(1)
          
          t14 = abs(p14)/max(abs(m12),abs(l14*y(1)))
          t12 = abs(p12)/max(abs(m12),abs(l12*y(1)))
          t34 = abs(p34)/max(abs(m32),abs(l34*z(1)))
          t23 = abs(p23)/max(abs(m32),abs(l23*z(1)))
          
          if (min(t12,t14).gt.min(t34,t23)) then
            if(dimag(p14/p12).ne.0d0)then
              v= cdln(p14/p12,0d0)
            else
              v=cdln(p14,-1d0)-cdln(p12,-1d0)
            endif
          else
            if(dimag(p34/p23).ne.0d0)then
              v= cdln(p34/p23,0d0)
            else
              v=cdln(p34,-1d0)-cdln(p23,-1d0)
            endif
          endif

c>          write(*,*) 'cD0msp0o v',p14/p12,p34/p23              
c>          write(*,*) 'cD0msp0o v',cdln(p14,-1d0)-cdln(p12,-1d0),
c>     &        cdln(p34,-1d0)-cdln(p23,-1d0)
c>
          cD0msp0o = cD0msp0o 
     &        +eta*v

c>          write(*,*) 'cD0msp0o etacon ', - (-1)**k *v,eta
c>          write(*,*) 'cD0msp0o etasum ',cD0msp0o
        endif


        cD0msp0o = cD0msp0o - (
     &      cspcoe(-y(2),m2/(m1*r12),-iy(2),-ir12)     
     &      -.5d0*(cdln(-y(2),-iy(2))+cdln(m2/m1*r12,ir12))**2
     &      -cspcoe(cd0,m1/(m2*r12),iy(2),-ir12)
     &      -cspcoe(-y(2),m4/(m1*r14),-iy(2),-ir14)
     &      +.5d0*(cdln(-y(2),-iy(2))+cdln(m4/m1*r14,ir14))**2
     &      +cspcoe(cd0,m1/(m4*r14),iy(2),-ir14)
     &      +cspcoe(-z(2),m4/(m3*r34),-iz(2),-ir34)
     &      -.5d0*(cdln(-z(2),-iz(2))+cdln(m4/m3*r34,ir34))**2
     &      -cspcoe(cd0,m3/(m4*r34),iz(2),-ir34)
     &      -cspcoe(-z(2),m2/(m3*r23),-iz(2),-ir23)
     &      +.5d0*(cdln(-z(2),-iz(2))+cdln(m2/m3*r23,ir23))**2
     &      +cspcoe(cd0,m3/(m2*r23),iz(2),-ir23)
     &      ) 

c>      write(*,*) 'cD0msp0o args yz',
c>     &      -y(2),-iy(2),-z(2),-iz(2)
c>      write(*,*) 'cD0msp0o args r',
c>     &      m2/(m1*r12),m2/m1*r12,m4/(m1*r14),m4/m1*r14,
c>     &      m4/(m3*r34),m4/m3*r34,m2/(m3*r23),m2/m3*r23
c>      write(*,*) 'cD0msp0o args ',
c>     &      1d0+y(2)*m2/(m1*r12),
c>     &      1d0+y(2)*m2/m1*r12,
c>     &      1d0+y(2)*m4/(m1*r14),
c>     &      1d0+y(2)*m4/m1*r14,
c>     &      1d0+z(2)*m4/(m3*r34),
c>     &      1d0+z(2)*m4/m3*r34,
c>     &      1d0+z(2)*m2/(m3*r23),
c>     &      1d0+z(2)*m2/m3*r23
c>      write(*,*) 'cD0msp0o csp ',
c>     &      cspcoe(-y(2),m2/(m1*r12),-iy(2),-ir12),     
c>     &      -.5d0*(cdln(-y(2),-iy(2))+cdln(m2/m1*r12,ir12))**2
c>     &      -cspcoe(cd0,m1/(m2*r12),iy(2),-ir12),
c>     &      -cspcoe(-y(2),m4/(m1*r14),-iy(2),-ir14),
c>     &      +.5d0*(cdln(-y(2),-iy(2))+cdln(m4/m1*r14,ir14))**2
c>     &      +cspcoe(cd0,m1/(m4*r14),iy(2),-ir14),
c>     &      +cspcoe(-z(2),m4/(m3*r34),-iz(2),-ir34),
c>     &      -.5d0*(cdln(-z(2),-iz(2))+cdln(m4/m3*r34,ir34))**2
c>     &      -cspcoe(cd0,m3/(m4*r34),iz(2),-ir34),
c>     &      -cspcoe(-z(2),m2/(m3*r23),-iz(2),-ir23),
c>     &      +.5d0*(cdln(-z(2),-iz(2))+cdln(m2/m3*r23,ir23))**2
c>     &      +cspcoe(cd0,m3/(m2*r23),iz(2),-ir23)
c>
c>      write(*,*) 'cD0msp0o cspsum ',cD0msp0o


c *** improves numerical stability in extreme cases

        eta=olietae(-y(2),-1/z(2),-iy(2),iz(2),-ir13)

        if(eta.ne.0d0)then
          p14 = m42*y(2)*y(2)+l14*y(2)
          p12 = m22*y(2)*y(2)+l12*y(2)
          p34 = m42*z(2)*z(2)+l34*z(2)
          p23 = m22*z(2)*z(2)+l23*z(2)
          
          t14 = abs(p14)/max(abs(m42*y(2)*y(2)),abs(l14*y(2)))
          t12 = abs(p12)/max(abs(m22*y(2)*y(2)),abs(l12*y(2)))
          t34 = abs(p34)/max(abs(m42*z(2)*z(2)),abs(l34*z(2)))
          t23 = abs(p23)/max(abs(m22*z(2)*z(2)),abs(l23*z(2)))
          
          if (min(t12,t14).gt.min(t34,t23)) then
            if(dimag(p14/p12).ne.0d0)then
              v= cdln(p14/p12,0d0)
            else
              v=cdln(p14,-1d0)-cdln(p12,-1d0)
            endif
          else
            if(dimag(p34/p23).ne.0d0)then
              v= cdln(p34/p23,0d0)
            else
              v=cdln(p34,-1d0)-cdln(p23,-1d0)
            endif
          endif

c>          write(*,*) 'cD0msp0o v',p14/p12,p34/p23              
c>          write(*,*) 'cD0msp0o v',cdln(p14,-1d0)-cdln(p12,-1d0),
c>     &        cdln(p34,-1d0)-cdln(p23,-1d0)

          cD0msp0o = cD0msp0o 
     &        - eta*v

c>          write(*,*) 'cD0msp0o etacon ', - (-1)**k *v,eta
c>          write(*,*) 'cD0msp0o etasum ',cD0msp0o
        endif

c      write(*,*) 'cD0msp0o sum ',cD0msp0o

      cD0msp0o = cD0msp0o/(m1*m2*m3*m4*b)

c      write(*,*) 'cD0msp0o =   ',cD0msp0o    

c end old version
      
      if (abs((cD0msp0-cD0msp0o)/cD0msp0).gt.testaccd0) then
        if (errflag.eq.0) then
          write(*,*) 
          write(*,*) 'cD0msp0 possibly wrong '
          errflag = 1
        end if
        write(errout,*) 
        write(errout,*) 'cD0msp0 possibly wrong '
        write(errout,*) 'cD0msp0 args ',
     &      q12,q23,q34,q14,q13,q24,m12,m22,m32,m42
        write(errout,*) 'cD0msp0o = ',cD0msp0o
        write(errout,*) 'cD0msp0  = ',cD0msp0
        write(errout,*) 'diff   = ',abs((cD0msp0-cD0msp0o)/cD0msp0)
      end if

      end if
      end

#endif
