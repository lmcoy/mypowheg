#include "global_coli.h"

#define QCONT

************************************************************************
*
*  BEWARE OF REAL ARGUMENTS OF CDLN AND CSPCOE, ...
*
*  dependence on msshift not yet tested at all
*
************************************************************************
*                                                                      *
*     Scalar 4- point functions                                        *
*                                                                      *
************************************************************************
*                                                                      *
*     last changed  16.12.08  Ansgar Denner                            *
*                                                                      *
************************************************************************
* Subroutines:                                                         *
* Functions:                                                           *
* D0_coli                                                              *
* D0reg_coli, D0m0_coli, D02m0_coli, D03m0_coli,D04m0_coli             *
* D0ms0ir1_coli, D0ms0ir1_coli, D0ms0ir1m0_coli                        *
* D0ms1ir0_coli, D0ms1ir1_coli, D0ms1ir2_coli                          *
* D0ms2ir0_coli, D0ms2ir1_coli, D0ms2ir2_coli, D0ms2ir3_coli           *
* D0ms3ir2_coli, D0ms4ir4_coli                                         *
* D0_cold                                                              *
************************************************************************

************************************************************************
      function D0_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     SCALAR 4-POINT FUNCTION                                          *
*                                                                      *
*                     m22                                              *
*       m12  ---------------------  p23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*              m12| 1     3 | m32                                      *
*                 |         |                                          *
*                 |    4    |                                          *
*       p14  ---------------------  p34                                *
*                     m42                                              *
*                                                                      *
************************************************************************
*     04.08.08 Ansgar Denner        last changed  10.12.08             *
************************************************************************
      implicit   none
      complex*16 D0_coli
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      complex*16 p2(4,4),m2(4)

      

      complex*16 elimminf2_coli
      complex*16 D0reg_coli 
      complex*16 D0m0_coli,D02m0_coli,D03m0_coli,D04m0_coli
      complex*16 D0ms0ir1_coli,D0ms0ir2_coli,D0ms0ir1m0_coli
      complex*16 D0ms1ir0_coli,D0ms1ir1_coli,D0ms1ir2_coli
      complex*16 D0ms2ir0_coli,D0ms2ir1_coli,D0ms2ir2_coli,
     &    D0ms2ir3_coli
      complex*16 D0ms3ir2_coli,D0ms4ir4_coli

      integer i,j,k
      integer nsm,nsoft,ncoll
      logical smallm2(4),smallp2(4,4),soft(4,4,4),coll(4,4),onsh(4,4)

      data coll /16*.false./, onsh /16*.false./,
     &    soft /64*.false./ ,smallm2 /4*.false./, smallp2 /16*.false./

#include "params_coli.h"
#include "checkparams_coli.h"

#ifdef CHECK
 101  format(a22,g25.17)
#endif
 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))

c      write(*,*)' D0_coli in',p12,p23,p34,p14,p13,p24,m12,m22,m32,m42


      p2(1,2) = p12
      p2(1,3) = p13
      p2(1,4) = p14    
      p2(2,3) = p23
      p2(2,4) = p24
      p2(3,4) = p34
      m2(1)   = m12
      m2(2)   = m22
      m2(3)   = m32
      m2(4)   = m42

c determine infinitesimal parameters
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

c      write(*,*)' D0_coli small',
c     &    ps12,ps23,ps34,ps14,ps13,ps24,ms12,ms22,ms32,ms42

      smallm2(1)=ms12.eq.cd0
      smallm2(2)=ms22.eq.cd0
      smallm2(3)=ms32.eq.cd0
      smallm2(4)=ms42.eq.cd0
      smallp2(1,2)=ps12.eq.cd0
      smallp2(1,3)=ps13.eq.cd0
      smallp2(1,4)=ps14.eq.cd0
      smallp2(2,3)=ps23.eq.cd0
      smallp2(2,4)=ps24.eq.cd0
      smallp2(3,4)=ps34.eq.cd0

      nsm=0
      do i=1,4
        if(smallm2(i)) nsm=nsm+1
      enddo

c determine on-shell momenta
      onsh(1,2)=ps12.eq.ms22
      onsh(1,3)=ps13.eq.ms32
      onsh(1,4)=ps14.eq.ms42
      onsh(2,1)=ps12.eq.ms12
      onsh(2,3)=ps23.eq.ms32
      onsh(2,4)=ps24.eq.ms42
      onsh(3,1)=ps13.eq.ms12
      onsh(3,2)=ps23.eq.ms22
      onsh(3,4)=ps34.eq.ms42
      onsh(4,1)=ps14.eq.ms12
      onsh(4,2)=ps24.eq.ms22
      onsh(4,3)=ps34.eq.ms32

c determine collinear singularities
      ncoll=0
      do i=1,3
      do j=i,4
        coll(i,j)=smallm2(i).and.smallm2(j).and.smallp2(i,j)
        if(coll(i,j).and.(.not.
     &      ( (p2(i,j).eq.cd0.and.m2(i).eq.m2(j)).or.
     &        (m2(i).eq.0d0.and.onsh(i,j)).or.  
     &        (m2(j).eq.0d0.and.onsh(j,i)) ))) then
          write(*,100)' D0_coli: structure of collinear singularity ',
     &        '    not supported'
          write(*,111)' D0_coli: p12 = ',p12
          write(*,111)' D0_coli: p23 = ',p23
          write(*,111)' D0_coli: p34 = ',p34
          write(*,111)' D0_coli: p14 = ',p14
          write(*,111)' D0_coli: p24 = ',p24
          write(*,111)' D0_coli: p13 = ',p13
          write(*,111)' D0_coli: m12 = ',m12
          write(*,111)' D0_coli: m22 = ',m22
          write(*,111)' D0_coli: m32 = ',m32
          write(*,111)' D0_coli: m42 = ',m42
          write(*,*)  ' D0_coli: i,j = ',i,j
          write(*,*)  ' D0_coli: t1  = ',p2(i,j).eq.cd0,m2(i).eq.m2(j)
          write(*,*)  ' D0_coli: t2  = ',m2(i).eq.cd0,onsh(i,j)
          write(*,*)  ' D0_coli: t3  = ',m2(j).eq.cd0,onsh(j,i)
          stop
        endif
        coll(j,i)=coll(i,j)
        if(coll(i,j)) ncoll=ncoll+1
       
c        if(coll(i,j)) write(*,*) 'D0_coli coll ',i,j

      enddo
      enddo

c determine soft singularities
      nsoft=0 
      do i=1,4
        do j=1,3
          if(i.ne.j)then
            do k=j,4
              if(k.ne.i.and.k.ne.j)then
                soft(i,j,k)=smallm2(i).and.onsh(i,j).and.onsh(i,k)
                soft(i,k,j)=soft(i,j,k)
                if(soft(i,j,k)) nsoft=nsoft+1

c                if(soft(i,j,k)) write(*,*) 'D0_coli soft ',i,j,k

              endif
            enddo
          endif
        enddo
      enddo

c      write(*,*) ' D0_coli: ncoll,nsoft,nsm= ',ncoll,nsoft,nsm

c regular cases
      if(ncoll.eq.0.and.nsoft.eq.0)then
        if(nsm.eq.0)then
          D0_coli = D0reg_coli(ps12,ps23,ps34,ps14,ps13,ps24,
     &        ms12,ms22,ms32,ms42)
        elseif (nsm.eq.1)then
          if (ms12.eq.cd0) then
            D0_coli = D0m0_coli(ps34,ps14,ps12,ps23,ps13,ps24,
     &          ms32,ms42,cd0,ms22)
          elseif (ms22.eq.cd0) then
            D0_coli = D0m0_coli(ps14,ps12,ps23,ps34,ps24,ps13,
     &          ms42,ms12,cd0,ms32)
          elseif (ms32.eq.cd0) then
            D0_coli = D0m0_coli(ps12,ps23,ps34,ps14,ps13,ps24,
     &          ms12,ms22,cd0,ms42)
          elseif (ms42.eq.cd0) then
            D0_coli = D0m0_coli(ps23,ps34,ps14,ps12,ps24,ps13,
     &          ms22,ms32,cd0,ms12)
#ifdef CHECK
          else
            write(*,100)' D0_coli: inconsistency 0ir 0ms 1m0'
#endif
          endif
        elseif(nsm.eq.2)then
          if (ms12.eq.cd0) then
            if(ms22.eq.cd0)then
              D0_coli = D02m0_coli(ps14,ps12,ps23,ps34,ps24,ps13,
     &            ms42,cd0,cd0,ms32)
            elseif(ms32.eq.cd0)then
              D0_coli = D02m0_coli(ps12,ps13,ps34,ps24,ps23,ps14,
     &            ms22,cd0,cd0,ms42)
            elseif(ms42.eq.cd0)then
              D0_coli = D02m0_coli(ps34,ps14,ps12,ps23,ps13,ps24,
     &            ms32,cd0,cd0,ms22)
#ifdef CHECK
            else
              write(*,100)' D0_coli: inconsistency 0ir 0ms 2m0'
#endif
            endif
          elseif(ms22.eq.cd0) then
            if(ms32.eq.cd0)then
              D0_coli = D02m0_coli(ps12,ps23,ps34,ps14,ps13,ps24,
     &            ms12,cd0,cd0,ms42)
            elseif(ms42.eq.cd0)then
              D0_coli = D02m0_coli(ps12,ps24,ps34,ps13,ps14,ps23,
     &            ms12,cd0,cd0,ms32)
#ifdef CHECK
            else
              write(*,100)' D0_coli: inconsistency 0ir 0ms 2m0'
#endif
            endif
          elseif(ms32.eq.cd0) then
            if(ms42.eq.cd0)then
              D0_coli = D02m0_coli(ps23,ps34,ps14,ps12,ps24,ps13,
     &            ms22,cd0,cd0,ms12)
#ifdef CHECK
            else
              write(*,100)' D0_coli: inconsistency 0ir 0ms 2m0'
#endif
            endif
#ifdef CHECK
          else
            write(*,100)' D0_coli: inconsistency 0ir 0ms 2m0'
#endif
          endif
        elseif (nsm.eq.3)then
          if(ms42.ne.cd0
     &        .and.ms12.eq.cd0.and.ms22.eq.cd0.and.ms32.eq.cd0
     &        )then
            D0_coli = D03m0_coli(ps14,ps12,ps23,ps34,ps24,ps13,
     &          ms42,cd0,cd0,cd0)
          elseif(ms32.ne.cd0
     &        .and.ms12.eq.cd0.and.ms22.eq.cd0.and.ms42.eq.cd0
     &          )then
            D0_coli = D03m0_coli(ps34,ps14,ps12,ps23,ps13,ps24,
     &          ms32,cd0,cd0,cd0)
          elseif(ms22.ne.cd0
     &        .and.ms12.eq.cd0.and.ms42.eq.cd0.and.ms32.eq.cd0
     &          ) then
            D0_coli = D03m0_coli(ps23,ps34,ps14,ps12,ps24,ps13,
     &          ms22,cd0,cd0,cd0)
          elseif(ms12.ne.cd0
     &        .and.ms42.eq.cd0.and.ms22.eq.cd0.and.ms32.eq.cd0
     &          ) then
            D0_coli = D03m0_coli(ps12,ps23,ps34,ps14,ps13,ps24,
     &          ms12,cd0,cd0,cd0)
#ifdef CHECK
          else
            write(*,100)' D0_coli: inconsistency 0ir 0ms 3m0'
#endif
          endif
        elseif (nsm.eq.4)then
          if(ms12.eq.cd0.and.ms22.eq.cd0.and.ms32.eq.cd0.and.
     &        ms42.eq.cd0)then
            D0_coli = D04m0_coli(ps12,ps23,ps34,ps14,ps13,ps24,
     &          cd0,cd0,cd0,cd0)
#ifdef CHECK
          else
            write(*,100)' D0_coli: inconsistency 0ir 0ms 4m0'
#endif
          endif
#ifdef CHECK
        else
          write(*,100)' D0_coli: inconsistency 0ir 0ms'
#endif
        endif

c soft singular cases without collinear singularities
      elseif(ncoll.eq.0.and.nsoft.ge.1)then
        if(soft(1,2,4))then
          if(ms32.ne.cd0)then
            D0_coli = D0ms0ir1_coli(p12,ps23,ps34,p14,ps13,ps24,
     &          m12,m22,ms32,m42)
          elseif(soft(3,4,2))then
            D0_coli = D0ms0ir2_coli(p12,p23,p34,p14,ps13,ps24,
     &          m12,m22,m32,m42)
          else
            D0_coli = D0ms0ir1m0_coli(p12,ps23,ps34,p14,ps13,ps24,
     &          m12,m22,ms32,m42)
          endif            
        elseif(soft(2,3,1))then
          if(ms42.ne.cd0)then
            D0_coli = D0ms0ir1_coli(p23,ps34,ps14,p12,ps24,ps13,
     &          m22,m32,ms42,m12)
          elseif(soft(4,1,3))then
            D0_coli = D0ms0ir2_coli(p23,p34,p14,p12,ps24,ps13,
     &          m22,m32,m42,m12)
          else
            D0_coli = D0ms0ir1m0_coli(p23,ps34,ps14,p12,ps24,ps13,
     &          m22,m32,ms42,m12)
          endif            
        elseif(soft(3,4,2))then
          if(ms12.ne.cd0)then
            D0_coli = D0ms0ir1_coli(p34,ps14,ps12,p23,ps13,ps24,
     &          m32,m42,ms12,m22)
          elseif(soft(1,2,4))then
            D0_coli = D0ms0ir2_coli(p34,p14,p12,p23,ps13,ps24,
     &          m32,m42,m12,m22)
          else
            D0_coli = D0ms0ir1m0_coli(p34,ps14,ps12,p23,ps13,ps24,
     &          m32,m42,ms12,m22)
          endif            
        elseif(soft(4,1,3))then
          if(ms22.ne.cd0)then
            D0_coli = D0ms0ir1_coli(ps14,p12,p23,ps34,ps24,ps13,
     &          m42,m12,ms22,m32)
          elseif(soft(2,3,1))then
            D0_coli = D0ms0ir2_coli(p14,p12,p23,p34,ps24,ps13,
     &          m42,m12,m22,m32)
          else
            D0_coli = D0ms0ir1m0_coli(p14,ps12,ps23,p34,ps24,ps13,
     &          m42,m12,ms22,m32)
          endif            
        else if(soft(1,2,3))then
          if(ms42.ne.cd0)then
            D0_coli = D0ms0ir1_coli(p12,ps24,ps34,p13,ps14,ps23,
     &          m12,m22,ms42,m32)
          elseif(soft(4,3,2))then
            D0_coli = D0ms0ir2_coli(p12,p24,p34,p13,ps14,ps23,
     &          m12,m22,m42,m32)
          else
            D0_coli = D0ms0ir1m0_coli(p12,ps24,ps34,p13,ps14,ps23,
     &          m12,m22,ms42,m32)
          endif            
        elseif(soft(2,3,4))then
          if(ms12.ne.cd0)then
            D0_coli = D0ms0ir1_coli(p23,ps13,ps14,p24,ps12,ps34,
     &          m22,m32,ms12,m42)
          elseif(soft(1,4,3))then
            D0_coli = D0ms0ir2_coli(p23,p13,p14,p24,ps12,ps34,
     &          m22,m32,m12,m42)
          else
            D0_coli = D0ms0ir1m0_coli(p23,ps13,ps14,p24,ps12,ps34,
     &          m22,m32,ms12,m42)
          endif            
        elseif(soft(3,4,1))then
          if(ms22.ne.cd0)then
            D0_coli = D0ms0ir1_coli(p34,ps24,ps12,p13,ps23,ps14,
     &          m32,m42,ms22,m12)
          elseif(soft(2,1,4))then
            D0_coli = D0ms0ir2_coli(p34,p24,p12,p13,ps23,ps14,
     &          m32,m42,m22,m12)
          else
            D0_coli = D0ms0ir1m0_coli(p34,ps24,ps12,p13,ps23,ps14,
     &          m32,m42,ms22,m12)
          endif            
        elseif(soft(4,1,2))then
          if(ms32.ne.cd0)then
            D0_coli = D0ms0ir1_coli(ps14,p13,p23,ps24,ps34,ps12,
     &          m42,m12,ms32,m22)
          elseif(soft(3,2,1))then
            D0_coli = D0ms0ir2_coli(p14,p13,p23,p24,ps34,ps12,
     &          m42,m12,m32,m22)
          else
            D0_coli = D0ms0ir1m0_coli(p14,ps13,ps23,p24,ps34,ps12,
     &          m42,m12,ms32,m22)
          endif            
        elseif(soft(1,3,4))then
          if(ms22.ne.cd0)then
            D0_coli = D0ms0ir1_coli(p13,ps23,ps24,p14,ps12,ps34,
     &          m12,m32,ms22,m42)
          elseif(soft(2,4,3))then
            D0_coli = D0ms0ir2_coli(p13,p23,p24,p14,ps12,ps34,
     &          m12,m32,m22,m42)
          else
            D0_coli = D0ms0ir1m0_coli(p13,ps23,ps24,p14,ps12,ps34,
     &          m12,m32,ms22,m42)
          endif            
        elseif(soft(2,4,1))then
          if(ms32.ne.cd0)then
            D0_coli = D0ms0ir1_coli(p24,ps34,ps13,p12,ps23,ps14,
     &          m22,m42,ms32,m12)
          elseif(soft(3,1,2))then
            D0_coli = D0ms0ir2_coli(p24,p34,p13,p12,ps23,ps14,
     &          m22,m42,m32,m12)
          else
            D0_coli = D0ms0ir1m0_coli(p24,ps34,ps13,p12,ps23,ps14,
     &          m22,m42,ms32,m12)
          endif            
        elseif(soft(3,1,2))then
          if(ms42.ne.cd0)then
            D0_coli = D0ms0ir1_coli(p13,ps14,ps24,p23,ps34,ps12,
     &          m32,m12,ms42,m22)
          elseif(soft(4,2,1))then
            D0_coli = D0ms0ir2_coli(p13,p14,p24,p23,ps34,ps12,
     &          m32,m12,m42,m22)
          else
            D0_coli = D0ms0ir1m0_coli(p13,ps14,ps24,p23,ps34,ps12,
     &          m32,m12,ms42,m22)
          endif            
        elseif(soft(4,2,3))then
          if(ms12.ne.cd0)then
            D0_coli = D0ms0ir1_coli(ps24,p12,p13,ps34,ps14,ps23,
     &          m42,m22,ms12,m32)
          elseif(soft(1,3,2))then
            D0_coli = D0ms0ir2_coli(p24,p12,p13,p34,ps14,ps23,
     &          m42,m22,m12,m32)
          else
            D0_coli = D0ms0ir1m0_coli(p24,ps12,ps13,p34,ps14,ps23,
     &          m42,m22,ms12,m32)
          endif            
#ifdef CHECK
        else
            write(*,100)' D0_coli: inconsistency 0ms'
#endif
        endif

c single collinear case
      elseif(ncoll.eq.1)then
        if(coll(1,2))then
          if(soft(1,2,4))then
            if(soft(2,1,3))then
              D0_coli = D0ms1ir2_coli(p12,ps23,ps34,ps14,ps13,ps24,
     &            m12,m22,ms32,ms42)
            else
              D0_coli = D0ms1ir1_coli(p12,ps23,ps34,ps14,ps13,ps24,
     &            m12,m22,ms32,ms42)
            endif
          elseif(soft(2,1,3))then
              D0_coli = D0ms1ir1_coli(p12,ps14,ps34,ps23,ps24,ps13,
     &            m22,m12,ms42,ms32)
          elseif(soft(1,2,3))then
            if(soft(2,1,4))then
              D0_coli = D0ms1ir2_coli(p12,ps24,ps34,ps13,ps14,ps23,
     &            m12,m22,ms42,ms32)
            else
              D0_coli = D0ms1ir1_coli(p12,ps24,ps34,ps13,ps14,ps23,
     &            m12,m22,ms42,ms32)
            endif
          elseif(soft(2,1,4))then
              D0_coli = D0ms1ir1_coli(p12,ps13,ps34,ps24,ps23,ps14,
     &            m22,m12,ms32,ms42)
          else
            D0_coli = D0ms1ir0_coli(ps34,ps14,p12,ps23,ps13,ps24,
     &          ms32,ms42,m12,m22)
          endif
        elseif(coll(2,3))then
          if(soft(2,3,1))then
            if(soft(3,2,4))then
              D0_coli = D0ms1ir2_coli(p23,ps34,ps14,ps12,ps24,ps13,
     &            m22,m32,ms42,ms12)
            else
              D0_coli = D0ms1ir1_coli(p23,ps34,ps14,ps12,ps24,ps13,
     &            m22,m32,ms42,ms12)
            endif
          elseif(soft(3,2,4))then
              D0_coli = D0ms1ir1_coli(p23,ps12,ps14,ps34,ps13,ps24,
     &            m32,m22,ms12,ms42)
          elseif(soft(2,3,4))then
            if(soft(3,2,1))then
              D0_coli = D0ms1ir2_coli(p23,ps24,ps14,ps13,ps34,ps12,
     &            m32,m22,ms42,ms12)
            else
              D0_coli = D0ms1ir1_coli(p23,ps13,ps14,ps24,ps12,ps34,
     &            m22,m32,ms12,ms42)
            endif
          elseif(soft(3,2,1))then
              D0_coli = D0ms1ir1_coli(p23,ps24,ps14,ps13,ps34,ps12,
     &            m32,m22,ms42,ms12)
          else
            D0_coli = D0ms1ir0_coli(ps14,ps12,p23,ps34,ps24,ps13,
     &          ms42,ms12,m22,m32)
          endif
        elseif(coll(3,4))then
          if(soft(3,4,1))then
            if(soft(4,3,2))then
              D0_coli = D0ms1ir2_coli(p34,ps24,ps12,ps13,ps23,ps14,
     &            m32,m42,ms22,ms12)
            else
              D0_coli = D0ms1ir1_coli(p34,ps24,ps12,ps13,ps23,ps14,
     &            m32,m42,ms22,ms12)
            endif
          elseif(soft(4,3,2))then
              D0_coli = D0ms1ir1_coli(p34,ps13,ps12,ps24,ps14,ps23,
     &            m42,m32,ms12,ms22)
          elseif(soft(3,4,2))then
            if(soft(4,3,1))then
              D0_coli = D0ms1ir2_coli(p34,ps14,ps12,ps23,ps13,ps24,
     &            m32,m42,ms12,ms22)
            else
              D0_coli = D0ms1ir1_coli(p34,ps14,ps12,ps23,ps13,ps24,
     &            m32,m42,ms12,ms22)
            endif
          elseif(soft(4,3,1))then
              D0_coli = D0ms1ir1_coli(p34,ps23,ps12,ps14,ps24,ps13,
     &            m42,m32,ms22,ms12)
          else
            D0_coli = D0ms1ir0_coli(ps12,ps23,p34,ps14,ps13,ps24,
     &          ms12,ms22,m32,m42)
          endif
        elseif(coll(1,4))then
          if(soft(1,4,3))then
            if(soft(4,1,2))then
              D0_coli = D0ms1ir2_coli(p14,ps24,ps23,ps13,ps12,ps34,
     &            m12,m42,ms22,ms32)
            else
              D0_coli = D0ms1ir1_coli(p14,ps24,ps23,ps13,ps12,ps34,
     &            m12,m42,ms22,ms32)
            endif
          elseif(soft(4,1,2))then
              D0_coli = D0ms1ir1_coli(p14,ps13,ps23,ps24,ps34,ps12,
     &            m42,m12,ms32,ms22)
          elseif(soft(1,4,2))then
            if(soft(4,1,3))then
              D0_coli = D0ms1ir2_coli(p14,ps34,ps23,ps12,ps13,ps24,
     &            m12,m42,ms32,ms22)
            else
              D0_coli = D0ms1ir1_coli(p14,ps34,ps23,ps12,ps13,ps24,
     &            m12,m42,ms32,ms22)
            endif
          elseif(soft(4,1,3))then
              D0_coli = D0ms1ir1_coli(p14,ps12,ps23,ps34,ps24,ps13,
     &            m42,m12,ms22,ms32)
          else
            D0_coli = D0ms1ir0_coli(ps23,ps34,p14,ps12,ps24,ps13,
     &          ms22,ms32,m42,m12)
          endif
        elseif(coll(1,3))then
          if(soft(1,3,4))then
            if(soft(3,1,2))then
              D0_coli = D0ms1ir2_coli(p13,ps23,ps24,ps14,ps12,ps34,
     &            m12,m32,ms22,ms42)
            else
              D0_coli = D0ms1ir1_coli(p13,ps23,ps24,ps14,ps12,ps34,
     &            m12,m32,ms22,ms42)
            endif
          elseif(soft(3,1,2))then
              D0_coli = D0ms1ir1_coli(p13,ps14,ps24,ps23,ps34,ps12,
     &            m32,m12,ms42,ms22)
          elseif(soft(1,3,2))then
            if(soft(3,1,4))then
              D0_coli = D0ms1ir2_coli(p13,ps34,ps24,ps12,ps14,ps23,
     &            m12,m32,ms42,ms22)
            else
              D0_coli = D0ms1ir1_coli(p13,ps34,ps24,ps12,ps14,ps23,
     &            m12,m32,ms42,ms22)
            endif
          elseif(soft(3,1,4))then
              D0_coli = D0ms1ir1_coli(p13,ps12,ps24,ps34,ps23,ps14,
     &            m32,m12,ms22,ms42)
          else
            D0_coli = D0ms1ir0_coli(ps24,ps34,p13,ps12,ps23,ps14,
     &          ms22,ms42,m32,m12)
          endif
        elseif(coll(2,4))then
          if(soft(2,4,1))then
            if(soft(4,2,3))then
              D0_coli = D0ms1ir2_coli(p24,ps34,ps13,ps12,ps23,ps14,
     &            m22,m42,ms32,ms12)
            else
              D0_coli = D0ms1ir1_coli(p24,ps34,ps13,ps12,ps23,ps14,
     &            m22,m42,ms32,ms12)
            endif
          elseif(soft(4,2,3))then
              D0_coli = D0ms1ir1_coli(p24,ps12,ps13,ps34,ps14,ps23,
     &            m42,m22,ms12,ms32)
          elseif(soft(2,4,3))then
            if(soft(4,2,1))then
              D0_coli = D0ms1ir2_coli(p24,ps23,ps13,ps14,ps34,ps12,
     &            m42,m22,ms32,ms12)
            else
              D0_coli = D0ms1ir1_coli(p24,ps14,ps13,ps23,ps12,ps34,
     &            m22,m42,ms12,ms32)
            endif
          elseif(soft(4,2,1))then
              D0_coli = D0ms1ir1_coli(p24,ps23,ps13,ps14,ps34,ps12,
     &            m42,m22,ms32,ms12)
          else
            D0_coli = D0ms1ir0_coli(ps13,ps12,p24,ps34,ps23,ps14,
     &          ms32,ms12,m22,m42)
          endif
#ifdef CHECK
        else
          write(*,100)' D0_coli: inconsistency 1ms'
#endif
        endif

c double collinear case
      elseif(ncoll.eq.2)then
        if(coll(1,2))then
          if(coll(3,4))then
            D0_coli = D0ms2ir0_coli(p12,ps23,p34,ps14,ps13,ps24,
     &          m12,m22,m32,m42)
          elseif(coll(2,3))then
            if(soft(1,2,4))then
              if(soft(3,2,4))then
                D0_coli = D0ms2ir3_coli(p12,p23,ps34,ps14,ps13,ps24,
     &              m12,m22,m32,ms42)
              else
                D0_coli = D0ms2ir2_coli(p12,p23,ps34,ps14,ps13,ps24,
     &              m12,m22,m32,ms42)
              endif
            elseif(soft(3,2,4))then
              D0_coli = D0ms2ir2_coli(p23,p12,ps14,ps34,ps13,ps24,
     &            m32,m22,m12,ms42)
            else
              D0_coli = D0ms2ir1_coli(p23,ps34,ps14,p12,ps24,ps13,
     &              m22,m32,ms42,m12)
            endif
          elseif(coll(1,4))then
            if(soft(4,1,3))then
              if(soft(2,1,3))then
                D0_coli = D0ms2ir3_coli(p14,p12,ps23,ps34,ps24,ps13,
     &              m42,m12,m22,ms32)
              else
                D0_coli = D0ms2ir2_coli(p14,p12,ps23,ps34,ps24,ps13,
     &              m42,m12,m22,ms32)
              endif
            elseif(soft(2,1,3))then
              D0_coli = D0ms2ir2_coli(p12,p14,ps34,ps23,ps24,ps13,
     &            m22,m12,m42,m32)
            else
              D0_coli = D0ms2ir1_coli(p12,ps23,ps34,p14,ps13,ps24,
     &              m12,m22,ms32,m42)
            endif
          elseif(coll(2,4))then
            if(soft(1,2,3))then
              if(soft(4,2,3))then
                D0_coli = D0ms2ir3_coli(p12,p24,ps34,ps13,ps14,ps23,
     &              m12,m22,m42,ms32)
              else
                D0_coli = D0ms2ir2_coli(p12,p24,ps34,ps13,ps14,ps23,
     &              m12,m22,m42,ms32)
              endif
            elseif(soft(4,2,3))then
              D0_coli = D0ms2ir2_coli(p24,p12,ps13,ps34,ps14,ps23,
     &            m42,m22,m12,ms32)
            else
              D0_coli = D0ms2ir1_coli(p24,ps34,ps13,p12,ps23,ps14,
     &              m22,m42,ms32,m12)
            endif
          elseif(coll(1,3))then
            if(soft(3,1,4))then
              if(soft(2,1,4))then
                D0_coli = D0ms2ir3_coli(p13,p12,ps24,ps34,ps23,ps14,
     &              m32,m12,m22,ms42)
              else
                D0_coli = D0ms2ir2_coli(p13,p12,ps24,ps34,ps23,ps14,
     &              m32,m12,m22,ms42)
              endif
c 8.12.08
c>            elseif(soft(2,1,3))then
c>              D0_coli = D0ms2ir2_coli(p24,p34,ps13,ps12,ps23,ps14,
c>     &            m22,m42,m32,m12)
            elseif(soft(2,1,4))then
              D0_coli = D0ms2ir2_coli(p12,p13,ps34,ps24,ps23,ps14,
     &            m22,m12,m32,m42)
            else
              D0_coli = D0ms2ir1_coli(p12,ps24,ps34,p13,ps14,ps23,
     &              m12,m22,ms42,m32)
            endif
          else
            write(*,100)' D0_coli: inconsistency 2ms 12'
          endif
        elseif(coll(3,4))then

c          write(*,*)' D0_coli: ',coll(4,1),coll(3,2),coll(4,2),
c     &        coll(3,1) 

          if(coll(4,1))then
            if(soft(3,4,2))then
              if(soft(1,4,2))then
                D0_coli = D0ms2ir3_coli(p34,p14,ps12,ps23,ps13,ps24,
     &              m32,m42,m12,ms22)
              else
                D0_coli = D0ms2ir2_coli(p34,p14,ps12,ps23,ps13,ps24,
     &              m32,m42,m12,ms22)
              endif
            elseif(soft(1,4,2))then
              D0_coli = D0ms2ir2_coli(p14,p34,ps23,ps12,ps13,ps24,
     &            m12,m42,m32,ms22)
            else
              D0_coli = D0ms2ir1_coli(p14,ps12,ps23,p34,ps24,ps13,
     &              m42,m12,ms22,m32)
            endif
          elseif(coll(3,2))then
            if(soft(2,3,1))then
              if(soft(4,3,1))then
                D0_coli = D0ms2ir3_coli(p23,p34,ps14,ps12,ps24,ps13,
     &              m22,m32,m42,ms12)
              else
                D0_coli = D0ms2ir2_coli(p23,p34,ps14,ps12,ps24,ps13,
     &              m22,m32,m42,ms12)
              endif
            elseif(soft(4,3,1))then
              D0_coli = D0ms2ir2_coli(p34,p23,ps12,ps14,ps24,ps13,
     &            m42,m32,m22,m12)
            else
              D0_coli = D0ms2ir1_coli(p34,ps14,ps12,p23,ps13,ps24,
     &              m32,m42,ms12,m22)
            endif
          elseif(coll(4,2))then
            if(soft(3,4,1))then
              if(soft(2,4,1))then
                D0_coli = D0ms2ir3_coli(p34,p24,ps12,ps13,ps23,ps14,
     &              m32,m42,m22,ms12)
              else
                D0_coli = D0ms2ir2_coli(p34,p24,ps12,ps13,ps23,ps14,
     &              m32,m42,m22,ms12)
              endif
            elseif(soft(2,4,1))then
              D0_coli = D0ms2ir2_coli(p24,p34,ps13,ps12,ps23,ps14,
     &            m22,m42,m32,ms12)
            else
              D0_coli = D0ms2ir1_coli(p24,ps12,ps13,p34,ps14,ps23,
     &              m42,m22,ms12,m32)
            endif
          elseif(coll(3,1))then
            if(soft(1,3,2))then
              if(soft(4,3,2))then
                D0_coli = D0ms2ir3_coli(p13,p34,ps24,ps12,ps14,ps23,
     &              m12,m32,m42,ms22)
              else
                D0_coli = D0ms2ir2_coli(p13,p34,ps24,ps12,ps14,ps23,
     &              m12,m32,m42,ms22)
              endif
c 8.12.08
c            elseif(soft(4,3,1))then
c              D0_coli = D0ms2ir2_coli(p24,p12,ps13,ps34,ps14,ps23,
c     &            m42,m22,m12,m32)
            elseif(soft(4,3,2))then
              D0_coli = D0ms2ir2_coli(p34,p13,ps12,ps24,ps14,ps23,
     &            m42,m32,m12,ms22)
            else
              D0_coli = D0ms2ir1_coli(p34,ps24,ps12,p13,ps23,ps14,
     &              m32,m42,ms22,m12)
            endif
          else
            write(*,100)' D0_coli: inconsistency 2ms 34'
          endif

        elseif(coll(2,3))then
          if(coll(4,1))then
            D0_coli = D0ms2ir0_coli(p23,ps34,p14,ps12,ps24,ps13,
     &          m22,m32,m42,m12)
C>          elseif(coll(3,4))then
C>            if(soft(2,3,1))then
C>              if(soft(4,3,1))then
C>                D0_coli = D0ms2ir3_coli(p23,p34,ps14,ps12,ps24,ps13,
C>     &              m22,m32,m42,ms13)
C>              else
C>                D0_coli = D0ms2ir2_coli(p23,p34,ps14,ps12,ps24,ps13,
C>     &              m22,m32,m42,ms13)
C>              endif
C>            elseif(soft(4,3,1))then
C>              D0_coli = D0ms2ir2_coli(p34,p23,ps12,ps14,ps24,ps13,
C>     &            m42,m32,m22,ms13)
C>            else
C>              D0_coli = D0ms2ir1_coli(p34,ps14,ps12,p23,ps13,ps24,
C>     &              m32,m42,ms13,m22)
C>            endif
C>          elseif(coll(2,1))then
C>            if(soft(1,2,4))then
C>              if(soft(3,2,4))then
C>                D0_coli = D0ms2ir3_coli(p12,p23,ps34,ps14,ps13,ps24,
C>     &              m12,m22,m32,ms43)
C>              else
C>                D0_coli = D0ms2ir2_coli(p12,p23,ps34,ps14,ps13,ps24,
C>     &              m12,m22,m32,ms43)
C>              endif
C>            elseif(soft(3,2,4))then
C>              D0_coli = D0ms2ir2_coli(p23,p12,ps14,ps34,ps13,ps24,
C>     &            m32,m22,m12,ms42)
C>            else
C>              D0_coli = D0ms2ir1_coli(p23,ps34,ps14,p12,ps24,ps13,
C>     &              m22,m32,ms42,m12)
C>            endif
          elseif(coll(3,1))then
            if(soft(2,3,4))then
              if(soft(1,3,4))then
                D0_coli = D0ms2ir3_coli(p23,p13,ps14,ps24,ps12,ps34,
     &              m22,m32,m12,ms42)
              else
                D0_coli = D0ms2ir2_coli(p23,p13,ps14,ps24,ps12,ps34,
     &              m22,m32,m12,ms42)
              endif
            elseif(soft(1,3,4))then
              D0_coli = D0ms2ir2_coli(p13,p23,ps24,ps14,ps12,ps34,
     &            m12,m32,m22,ms42)
            else
              D0_coli = D0ms2ir1_coli(p13,ps14,ps24,p23,ps34,ps12,
     &              m32,m12,ms42,m22)
            endif
          elseif(coll(2,4))then
            if(soft(4,2,1))then
              if(soft(3,2,1))then
                D0_coli = D0ms2ir3_coli(p24,p23,ps13,ps14,ps34,ps12,
     &              m42,m22,m32,ms12)
              else
                D0_coli = D0ms2ir2_coli(p24,p23,ps13,ps14,ps34,ps12,
     &              m42,m22,m32,ms12)
              endif
c 8.12.08
c            elseif(soft(3,2,4))then
c              D0_coli = D0ms2ir2_coli(p13,p14,ps24,ps23,ps34,ps12,
c     &            m32,m12,m42,m22)
            elseif(soft(3,2,1))then
              D0_coli = D0ms2ir2_coli(p23,p24,ps14,ps13,ps34,ps12,
     &            m32,m22,m42,ms12)
            else
              D0_coli = D0ms2ir1_coli(p23,ps13,ps14,p24,ps12,ps34,
     &              m22,m32,ms12,m42)
            endif
          else
            write(*,100)' D0_coli: inconsistency 3ms 23'
          endif
        elseif(coll(4,1))then
C>          if(coll(1,2))then
C>            if(soft(4,1,3))then
C>              if(soft(2,1,3))then
C>                D0_coli = D0ms2ir3_coli(p14,p12,ps23,ps43,ps24,ps13,
C>     &              m42,m12,m22,ms22)
C>              else
C>                D0_coli = D0ms2ir2_coli(p14,p12,ps23,ps43,ps24,ps13,
C>     &              m42,m12,m22,ms22)
C>              endif
C>            elseif(soft(2,1,3))then
C>              D0_coli = D0ms2ir2_coli(p12,p14,ps43,ps23,ps24,ps13,
C>     &            m22,m12,m42,ms22)
C>            else
C>              D0_coli = D0ms2ir1_coli(p12,ps23,ps43,p14,ps13,ps24,
C>     &              m12,m22,ms32,m42)
C>            endif
C>          elseif(coll(4,3))then
C>            if(soft(3,4,2))then
C>              if(soft(1,4,2))then
C>                D0_coli = D0ms2ir3_coli(p43,p14,ps12,ps23,ps13,ps24,
C>     &              m32,m42,m12,ms22)
C>              else
C>                D0_coli = D0ms2ir2_coli(p43,p14,ps12,ps23,ps13,ps24,
C>     &              m32,m42,m12,ms22)
C>              endif
C>            elseif(soft(1,4,2))then
C>              D0_coli = D0ms2ir2_coli(p14,p43,ps23,ps12,ps13,ps24,
C>     &            m12,m42,m32,ms22)
C>            else
C>              D0_coli = D0ms2ir1_coli(p14,ps12,ps23,p43,ps24,ps13,
C>     &              m42,m12,ms22,m32)
C>            endif
          if(coll(1,3))then
            if(soft(4,1,2))then
              if(soft(3,1,2))then
                D0_coli = D0ms2ir3_coli(p14,p13,ps23,ps24,ps34,ps12,
     &              m42,m12,m32,ms22)
              else
                D0_coli = D0ms2ir2_coli(p14,p13,ps23,ps24,ps34,ps12,
     &              m42,m12,m32,ms22)
              endif
            elseif(soft(3,1,2))then
              D0_coli = D0ms2ir2_coli(p13,p14,ps24,ps23,ps34,ps12,
     &            m32,m12,m42,ms22)
            else
              D0_coli = D0ms2ir1_coli(p13,ps23,ps24,p14,ps12,ps34,
     &              m12,m32,ms22,m42)
            endif
          elseif(coll(4,2))then
            if(soft(2,4,3))then
              if(soft(1,4,3))then
                D0_coli = D0ms2ir3_coli(p24,p14,ps13,ps23,ps12,ps34,
     &              m22,m42,m12,ms32)
              else
                D0_coli = D0ms2ir2_coli(p24,p14,ps13,ps23,ps12,ps34,
     &              m22,m42,m12,ms32)
              endif
            elseif(soft(1,4,3))then
              D0_coli = D0ms2ir2_coli(p14,p24,ps23,ps13,ps12,ps34,
     &            m12,m42,m22,ms32)
            else
              D0_coli = D0ms2ir1_coli(p14,ps13,ps23,p24,ps34,ps12,
     &              m42,m12,ms32,m22)
            endif
          else
            write(*,100)' D0_coli: inconsistency 2ms 41'
          endif

        elseif(coll(1,3))then
          if(coll(2,4))then
            D0_coli = D0ms2ir0_coli(p13,ps23,p24,ps14,ps12,ps34,
     &          m12,m32,m22,m42)
C>          elseif(coll(3,2))then
C>            if(soft(1,3,4))then
C>              if(soft(2,3,4))then
C>                D0_coli = D0ms2ir3_coli(p13,p23,ps24,ps14,ps12,ps34,
C>     &              m12,m32,m22,ms42)
C>              else
C>                D0_coli = D0ms2ir2_coli(p13,p23,ps24,ps14,ps12,ps34,
C>     &              m12,m32,m22,ms42)
C>              endif
C>            elseif(soft(2,3,4))then
C>              D0_coli = D0ms2ir2_coli(p23,p13,ps14,ps24,ps12,ps34,
C>     &            m22,m32,m12,ms42)
C>            else
C>              D0_coli = D0ms2ir1_coli(p23,ps24,ps14,p13,ps34,ps12,
C>     &              m32,m22,ms42,m12)
C>            endif
C>          elseif(coll(1,4))then
C>            if(soft(4,1,2))then
C>              if(soft(3,1,2))then
C>                D0_coli = D0ms2ir3_coli(p14,p13,ps23,ps24,ps34,ps12,
C>     &              m42,m12,m32,ms23)
C>              else
C>                D0_coli = D0ms2ir2_coli(p14,p13,ps23,ps24,ps34,ps12,
C>     &              m42,m12,m32,ms23)
C>              endif
C>            elseif(soft(3,1,2))then
C>              D0_coli = D0ms2ir2_coli(p13,p14,ps24,ps23,ps34,ps12,
C>     &            m32,m12,m42,ms22)
C>            else
C>              D0_coli = D0ms2ir1_coli(p13,ps23,ps24,p14,ps12,ps34,
C>     &              m12,m32,ms23,m42)
C>            endif
C>          elseif(coll(3,4))then
C>            if(soft(1,3,2))then
C>              if(soft(4,3,2))then
C>                D0_coli = D0ms2ir3_coli(p13,p34,ps24,ps12,ps14,ps23,
C>     &              m12,m32,m42,ms23)
C>              else
C>                D0_coli = D0ms2ir2_coli(p13,p34,ps24,ps12,ps14,ps23,
C>     &              m12,m32,m42,ms23)
C>              endif
C>            elseif(soft(4,3,2))then
C>              D0_coli = D0ms2ir2_coli(p34,p13,ps12,ps24,ps14,ps23,
C>     &            m42,m32,m12,ms23)
C>            else
C>              D0_coli = D0ms2ir1_coli(p34,ps24,ps12,p13,ps23,ps14,
C>     &              m32,m42,ms23,m12)
C>            endif
C>          elseif(coll(1,2))then
C>            if(soft(2,1,4))then
C>              if(soft(3,1,4))then
C>                D0_coli = D0ms2ir3_coli(p12,p13,ps34,ps24,ps23,ps14,
C>     &              m22,m12,m32,ms42)
C>              else
C>                D0_coli = D0ms2ir2_coli(p12,p13,ps34,ps24,ps23,ps14,
C>     &              m22,m12,m32,ms42)
C>              endif
C>            elseif(soft(3,1,2))then
C>              D0_coli = D0ms2ir2_coli(p34,p24,ps12,ps13,ps23,ps14,
C>     &            m32,m42,m22,m12)
C>            else
C>              D0_coli = D0ms2ir1_coli(p13,ps34,ps24,p12,ps14,ps23,
C>     &              m12,m32,ms42,m22)
C>            endif
          else
            write(*,100)' D0_coli: inconsistency 3ms 13'
          endif
C>        elseif(coll(2,4))then
C>          if(coll(4,1))then
C>            if(soft(2,4,3))then
C>              if(soft(1,4,3))then
C>                D0_coli = D0ms2ir3_coli(p24,p14,ps13,ps23,ps12,ps34,
C>     &              m22,m42,m12,ms22)
C>              else
C>                D0_coli = D0ms2ir2_coli(p24,p14,ps13,ps23,ps12,ps34,
C>     &              m22,m42,m12,ms22)
C>              endif
C>            elseif(soft(1,4,3))then
C>              D0_coli = D0ms2ir2_coli(p14,p24,ps23,ps13,ps12,ps34,
C>     &            m12,m42,m22,ms22)
C>            else
C>              D0_coli = D0ms2ir1_coli(p14,ps13,ps23,p24,ps34,ps12,
C>     &              m42,m12,ms32,m22)
C>            endif
C>          elseif(coll(2,3))then
C>            if(soft(3,2,1))then
C>              if(soft(4,2,1))then
C>                D0_coli = D0ms2ir3_coli(p23,p24,ps14,ps13,ps34,ps12,
C>     &              m32,m22,m42,ms12)
C>              else
C>                D0_coli = D0ms2ir2_coli(p23,p24,ps14,ps13,ps34,ps12,
C>     &              m32,m22,m42,ms12)
C>              endif
C>            elseif(soft(4,2,1))then
C>              D0_coli = D0ms2ir2_coli(p24,p23,ps13,ps14,ps34,ps12,
C>     &            m42,m22,m32,ms12)
C>            else
C>              D0_coli = D0ms2ir1_coli(p24,ps14,ps13,p23,ps12,ps34,
C>     &              m22,m42,ms12,m32)
C>            endif
C>          elseif(coll(4,3))then
C>            if(soft(2,4,1))then
C>              if(soft(3,4,1))then
C>                D0_coli = D0ms2ir3_coli(p24,p43,ps13,ps12,ps23,ps14,
C>     &              m22,m42,m32,ms12)
C>              else
C>                D0_coli = D0ms2ir2_coli(p24,p43,ps13,ps12,ps23,ps14,
C>     &              m22,m42,m32,ms12)
C>              endif
C>            elseif(soft(3,4,1))then
C>              D0_coli = D0ms2ir2_coli(p43,p24,ps12,ps13,ps23,ps14,
C>     &            m32,m42,m22,ms12)
C>            else
C>              D0_coli = D0ms2ir1_coli(p43,ps13,ps12,p24,ps14,ps23,
C>     &              m42,m32,ms12,m22)
C>            endif
C>          elseif(coll(2,1))then
C>            if(soft(1,2,3))then
C>              if(soft(4,2,3))then
C>                D0_coli = D0ms2ir3_coli(p12,p24,ps34,ps13,ps14,ps23,
C>     &              m12,m22,m42,ms32)
C>              else
C>                D0_coli = D0ms2ir2_coli(p12,p24,ps34,ps13,ps14,ps23,
C>     &              m12,m22,m42,ms32)
C>              endif
C>            elseif(soft(4,2,1))then
C>              D0_coli = D0ms2ir2_coli(p43,p13,ps12,ps24,ps14,ps23,
C>     &            m42,m32,m12,m22)
C>            else
C>              D0_coli = D0ms2ir1_coli(p24,ps34,ps13,p12,ps23,ps14,
C>     &              m22,m42,ms32,m12)
C>            endif
C>          else
C>            write(*,100)' D0_coli: inconsistency 2ms 24'
C>        endif

        else
          write(*,100)' D0_coli: inconsistency 2ms'
        endif

c triple collinear case
      elseif(ncoll.eq.3)then

c        write(*,*) 'D0_coli   ncoll=3'

        if(coll(1,2))then

c          write(*,*) 'D0_coli   coll12'
        
          if(coll(1,4))then
            if(coll(2,3))then
              D0_coli = D0ms3ir2_coli(p12,p23,ps34,p14,ps13,ps24,
     &            m12,m22,m32,m42)
            elseif(coll(3,4))then
              D0_coli = D0ms3ir2_coli(p14,p12,ps23,p34,ps24,ps13,
     &            m42,m12,m22,m32)
#ifdef CHECK
            else 
              write(*,100)' D0_coli: inconsistency 3ms 12 14'
#endif
            endif
          elseif(coll(2,3))then
            if(coll(3,4))then
              D0_coli = D0ms3ir2_coli(p23,p34,ps14,p12,ps24,ps13,
     &            m22,m32,m42,m12)
#ifdef CHECK
            else 
              write(*,100)' D0_coli: inconsistency 3ms 12 33'
#endif
            endif
          elseif(coll(1,3))then
            if(coll(2,4))then
              D0_coli = D0ms3ir2_coli(p12,p13,ps34,p24,ps23,ps14,
     &            m22,m12,m32,m42)
            elseif(coll(3,4))then
              D0_coli = D0ms3ir2_coli(p13,p34,ps24,p12,ps14,ps23,
     &            m12,m32,m42,m22)
#ifdef CHECK
            else 
              write(*,100)' D0_coli: inconsistency 3ms 12 13'
#endif
            endif
          elseif(coll(2,4))then
            if(coll(3,4))then
              D0_coli = D0ms3ir2_coli(p24,p34,ps13,p12,ps23,ps14,
     &            m22,m42,m32,m12)
#ifdef CHECK
            else 
              write(*,100)' D0_coli: inconsistency 3ms 12 24'
#endif
            endif
#ifdef CHECK
          else 
            write(*,100)' D0_coli: inconsistency 3ms 12 24'
#endif
          endif

        elseif(coll(3,4))then

c          write(*,*) 'D0_coli   coll34'

          if(coll(2,3))then

c            write(*,*) 'D0_coli   coll23'

            if(coll(1,4))then

c            write(*,*) 'D0_coli   coll14'

              D0_coli = D0ms3ir2_coli(p34,p14,ps12,p23,ps13,ps24,
     &            m32,m42,m12,m22)
#ifdef CHECK
            else 
              write(*,100)' D0_coli: inconsistency 3ms 34 23'
#endif
            endif
          elseif(coll(1,3))then
            if(coll(2,4))then
              D0_coli = D0ms3ir2_coli(p34,p13,ps12,p24,ps14,ps23,
     &            m42,m32,m12,m22)
#ifdef CHECK
            else 
              write(*,100)' D0_coli: inconsistency 3ms 34 13'
#endif
            endif
#ifdef CHECK
          else 
            write(*,100)' D0_coli: inconsistency 3ms 34'
#endif
          endif

        elseif(coll(1,4))then
          if(coll(1,3))then
            if(coll(2,4))then
              D0_coli = D0ms3ir2_coli(p14,p13,ps23,p24,ps34,ps12,
     &            m42,m12,m32,m22)
            elseif(coll(2,3))then
              D0_coli = D0ms3ir2_coli(p13,p23,ps24,p14,ps12,ps34,
     &            m12,m32,m22,m42)
#ifdef CHECK
            else 
              write(*,100)' D0_coli: inconsistency 3ms 14 13'
#endif
            endif
          elseif(coll(2,4))then
            if(coll(2,3))then
              D0_coli = D0ms3ir2_coli(p24,p23,ps13,p14,ps34,ps12,
     &            m42,m22,m32,m12)
#ifdef CHECK
            else 
              write(*,100)' D0_coli: inconsistency 3ms 14 24'
#endif
            endif
#ifdef CHECK
          else 
            write(*,100)' D0_coli: inconsistency 3ms 14'
#endif
          endif

        elseif(coll(2,3))then
          if(coll(2,4))then
            if(coll(1,3))then
              D0_coli = D0ms3ir2_coli(p23,p24,ps14,p13,ps34,ps12,
     &            m32,m22,m42,m12)
#ifdef CHECK
            else 
              write(*,100)' D0_coli: inconsistency 3ms 23 24'
#endif
            endif
#ifdef CHECK
          else 
            write(*,100)' D0_coli: inconsistency 3ms 23'
#endif
          endif

#ifdef CHECK
        else 
          write(*,100)' D0_coli: inconsistency 3ms'
#endif
        endif

c quartic collinear case
      elseif(ncoll.eq.4)then
        if(coll(1,2).and.coll(2,3)
#ifdef CHECK
     &      .and.coll(3,4).and.coll(1,4)
#endif
     &      )then 
          D0_coli = D0ms4ir4_coli(p12,p23,p34,p14,ps13,ps24,
     &              m12,m22,m32,m42)
        elseif(coll(1,2).and.coll(1,3)
#ifdef CHECK
     &        .and.coll(3,4).and.coll(2,4)
#endif
     &        )then 
          D0_coli = D0ms4ir4_coli(p12,p24,p34,p13,ps14,ps23,
     &              m12,m22,m42,m32)
        elseif(coll(2,3).and.coll(1,3)
#ifdef CHECK
     &        .and.coll(1,4).and.coll(2,4)
#endif
     &        )then 
          D0_coli = D0ms4ir4_coli(p13,p23,p24,p14,ps12,ps34,
     &              m12,m32,m22,m42)
#ifdef CHECK
        else 
          write(*,100)' D0_coli: inconsistency 4ms'
#endif
        endif   

#ifdef CHECK
      else
        write(*,100)' D0_coli: case with more than 4 collinear ',
     &      'singularities not implemented'
#endif
      endif

      end

************************************************************************
      function D0reg_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*  general scalar 4-point function                                     *
*  regular case  based on general result of                            *
*        A.Denner, U.Nierste and R.Scharf, Nucl. Phys. B367 (1991) 637 *
*                                                                      *
*                     m22                                              *
*       p12  ---------------------  p23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*              m12| 1     3 | m32                                      *
*                 |         |                                          *
*                 |    4    |                                          *
*       p14  ---------------------  p34                                *
*                     m42                                              *
*                                                                      *
*----------------------------------------------------------------------*
*  30.04.08 Ansgar Denner       last changed 30.04.08 Ansgar Denner    *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 l12,l13,l14,l23,l24,l34
      complex*16 r12,r13,r14,r23,r24,r34
      complex*16 r21,r31,r41,r32,r42,r43
      complex*16 D0reg_coli,D0regrp_coli
c      integer    k
c      complex*16 cd0
c      data       cd0 /(0d0,0d0)/
 
#ifdef CHECK
c      complex*16 D0reg_check
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      complex*16 elimminf2_coli
c      integer    j
c      logical flag(0:1)
c      data    flag /2*.true./
c      save    flag  
#endif

c      complex*16 cln_coli,cspenc_coli

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) ' D0reg_coli in ',p12,p23,p34,p14,p13,p24,m12,m22,m32,m42

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
          write(*,100)' D0reg_coli called improperly:'
          write(*,111)' D0reg_coli: p12 = ',p12
          write(*,111)' D0reg_coli: p23 = ',p23
          write(*,111)' D0reg_coli: p34 = ',p34
          write(*,111)' D0reg_coli: p14 = ',p14
          write(*,111)' D0reg_coli: p13 = ',p13
          write(*,111)' D0reg_coli: p24 = ',p24
          write(*,111)' D0reg_coli: m12 = ',m12
          write(*,111)' D0reg_coli: m22 = ',m22
          write(*,111)' D0reg_coli: m32 = ',m32
          write(*,111)' D0reg_coli: m42 = ',m42
          stop
        end if
      end if
#endif

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


c      r24=1d0/r42
c      r42 = r24*m42/m22


      if(dimag(r31).ne.0d0.and.
     &    dimag(m22*r31*r31).le.0d0.and.dimag(m42*r31*r31).le.0d0.and.
     &    dimag(r31*l23).le.0d0.and.dimag(r31*l34).le.0d0.and.
     &    dimag(r31*r31*l24).le.0d0 .or.
     &    dimag(r31).eq.0d0.and.dreal(r31).gt.0d0 ) then
        D0reg_coli=
     &      D0regrp_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
      elseif(dimag(r42).ne.0d0.and.
     &      dimag(m32*r42*r42).le.0d0.and.dimag(m12*r42*r42).le.0d0
     &      .and.dimag(r42*l34).le.0d0.and.dimag(r24*l14).le.0d0.and.
     &      dimag(r42*r42*l13).le.0d0.or.
     &      dimag(r42).eq.0d0.and.dreal(r42).gt.0d0 ) then
c        write(*,*) 'D0reg_coli: used rp perm 2341'
        D0reg_coli=
     &      D0regrp_coli(p23,p34,p14,p12,p24,p13,m22,m32,m42,m12)
      elseif(dimag(r21).ne.0d0.and.
     &    dimag(m32*r21*r21).le.0d0.and.dimag(m42*r21*r21).le.0d0.and.
     &    dimag(r21*l23).le.0d0.and.dimag(r21*l24).le.0d0.and.
     &    dimag(r21*r21*l34).le.0d0 .or.
     &    dimag(r21).eq.0d0.and.dreal(r21).gt.0d0 ) then
c        write(*,*) 'D0reg_coli: used rp perm 1324'
        D0reg_coli=
     &      D0regrp_coli(p13,p23,p24,p14,p12,p34,m12,m32,m22,m42)
      elseif(dimag(r43).ne.0d0.and.
     &      dimag(m22*r43*r43).le.0d0.and.dimag(m12*r43*r43).le.0d0
     &      .and.dimag(r42*l24).le.0d0.and.dimag(r34*l14).le.0d0.and.
     &      dimag(r43*r43*l12).le.0d0.or.
     &      dimag(r43).eq.0d0.and.dreal(r43).gt.0d0 ) then
        write(*,*) 'D0reg_coli: used rp perm 3241'
        D0reg_coli=
     &      D0regrp_coli(p23,p34,p14,p12,p24,p13,m22,m32,m42,m12)
      elseif(dimag(r41).ne.0d0.and.
     &    dimag(m22*r41*r41).le.0d0.and.dimag(m32*r41*r41).le.0d0.and.
     &    dimag(r41*l24).le.0d0.and.dimag(r41*l34).le.0d0.and.
     &    dimag(r41*r41*l23).le.0d0 .or.
     &    dimag(r41).eq.0d0.and.dreal(r41).gt.0d0 ) then
        write(*,*) 'D0reg_coli: used rp perm 1243'
        D0reg_coli=
     &      D0regrp_coli(p12,p24,p34,p13,p14,p23,m12,m22,m42,m32)
      elseif(dimag(r32).ne.0d0.and.
     &      dimag(m42*r32*r32).le.0d0.and.dimag(m12*r32*r32).le.0d0
     &      .and.dimag(r32*l34).le.0d0.and.dimag(r23*l13).le.0d0.and.
     &      dimag(r32*r32*l14).le.0d0.or.
     &      dimag(r32).eq.0d0.and.dreal(r32).gt.0d0 ) then
        write(*,*) 'D0reg_coli: used rp perm 2431'
        D0reg_coli=
     &      D0regrp_coli(p24,p34,p13,p12,p23,p14,m22,m42,m32,m12)
       elseif(dimag(r42).ne.0d0.and.
     &      dimag(m32*r24*r24).le.0d0.and.dimag(m12*r24*r24).le.0d0
     &      .and.dimag(r24*l23).le.0d0.and.dimag(r42*l12).le.0d0.and.
     &      dimag(r24*r24*l13).le.0d0.or.
     &      dimag(r24).eq.0d0.and.dreal(r24).gt.0d0 ) then

        write(*,*) 'D0reg_coli: used rp perm 4321'
        D0reg_coli=
     &      D0regrp_coli(p34,p23,p12,p14,p24,p13,m42,m32,m22,m12)
        

      else
        write(*,100)' D0reg_coli: case not implemented'
        write(*,111)' D0reg_coli: p12 = ',p12
        write(*,111)' D0reg_coli: p23 = ',p23
        write(*,111)' D0reg_coli: p34 = ',p34
        write(*,111)' D0reg_coli: p14 = ',p14
        write(*,111)' D0reg_coli: p13 = ',p13
        write(*,111)' D0reg_coli: p24 = ',p24
        write(*,111)' D0reg_coli: m12 = ',m12
        write(*,111)' D0reg_coli: m22 = ',m22
        write(*,111)' D0reg_coli: m32 = ',m32
        write(*,111)' D0reg_coli: m42 = ',m42
        write(*,111)' D0reg_coli: r13 = ',r13,r31
        write(*,111)' D0reg_coli: r24 = ',r24,r42
        write(*,111)' D0reg_coli: r12 = ',r12,r21
        write(*,111)' D0reg_coli: r14 = ',r14,r41
        write(*,111)' D0reg_coli: r23 = ',r23,r32
        write(*,111)' D0reg_coli: r34 = ',r34,r43
        stop
      endif



      end




************************************************************************
      function D0m0_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*  scalar 4-point function  for m32 = 0                                *
*  regular case                                                        *
*                                                                      *
*                     m22                                              *
*       p12  ---------------------  p23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*              m12| 1     3 | 0d0                                      *
*                 |         |                                          *
*                 |    4    |                                          *
*       p14  ---------------------  p34                                *
*                     m42                                              *
*                                                                      *
*----------------------------------------------------------------------*
*  29.03.92 Ansgar Denner       last changed  16.02.07 Ansgar Denner   *
************************************************************************
      implicit   none
      complex*16 D0m0_coli
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42

      integer    ltest
      common /ltest/ ltest
 
      complex*16 k12,k13,k14,k23,k24,k34
      complex*16 l12,l13,l14,l23,l24,l34
      complex*16 m1,m2,mm3,m4
      real*8     ir12,ir14,ir24
      real*8     ix1(2),ix4(2),iqbard(2)
      real*8     test0,test1,test2,test3,test4,test5
      real*8     test01,test23,test45
      real*8     u,v
      complex*16 r12,r14,r24,r21,r41,r42
      complex*16 a,b,c,d,det
      complex*16 x1(2),x4(2)
      complex*16 ch1,ch2,ch3,l1,ch4,ch5,l2
      complex*16 v1,v4,eta
      complex*16 cspcos_coli,cln_coli,eta2s_coli
      integer    i,j

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
      complex*16 D0m0_check
      logical flag(0:2)
      data    flag /3*.true./
      save    flag  
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) ' D0m0_coli in ',
c     &    p12,p23,p34,p14,p13,p24,m12,m22,m32,m42

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

        if(m32.ne.0d0.or.m12.eq.cd0.or.m22.eq.cd0.or.m42.eq.cd0
     &      .or.p23.eq.m22.and.p13.eq.m12
     &      .or.p23.eq.m22.and.p34.eq.m42
     &      .or.p13.eq.m12.and.p34.eq.m42
     &      .or.ps23.eq.cd0.and.ms22.eq.cd0
     &      .or.ps34.eq.cd0.and.ms42.eq.cd0
     &      .or.ps13.eq.cd0.and.ms12.eq.cd0
     &      .or.ps12.eq.cd0.and.ms12.eq.cd0.and.ms22.eq.cd0
c     &      .or.ps14.eq.cd0.and.ms12.eq.cd0.and.ms42.eq.cd0
     &      .or.ps24.eq.cd0.and.ms22.eq.cd0.and.ms22.eq.cd0
     &      .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &      .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &      .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0m0_coli called improperly:'
          write(*,111)' D0m0_coli: p12 = ',p12
          write(*,111)' D0m0_coli: p23 = ',p23
          write(*,111)' D0m0_coli: p34 = ',p34
          write(*,111)' D0m0_coli: p14 = ',p14
          write(*,111)' D0m0_coli: p13 = ',p13
          write(*,111)' D0m0_coli: p24 = ',p24
          write(*,111)' D0m0_coli: m12 = ',m12
          write(*,111)' D0m0_coli: m22 = ',m22
          write(*,111)' D0m0_coli: m32 = ',m32
          write(*,111)' D0m0_coli: m42 = ',m42
          write(*,*)' D0m0_coli: test= ',
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

c>      if (m12.eq.0d0) then
c>        D0m0_coli = D02m0_coli(p12,p13,p34,p24,p23,p14,m22,cd0,cd0,m42)
c>      else if (m22.eq.0d0) then
c>        D0m0_coli = D02m0_coli(p12,p23,p34,p14,p13,p24,m12,cd0,cd0,m42)
c>      else if (m42.eq.0d0) then
c>        D0m0_coli = D02m0_coli(p14,p34,p23,p12,p13,p24,m12,cd0,cd0,m22)
c>      else

c     (permutation 1 -> 3 -> 2 -> 4 -> 1 for  m -> k )
c        m1 -> m2 in normalization

        l12 = (m22+m12-p12)
        l13 = (m12    -p13)
        l14 = (m42+m12-p14)
        l23 = (m22    -p23)
        l24 = (m22+m42-p24)
        l34 = (    m42-p34)
        r12 = l12/(2d0*m22)*(1d0+sqrt(1d0-4d0*m12*m22/l12**2))
        r21 = r12*m22/m12

c>        write(*,*) 'D0m0_coli kij ',l12/sqrt(m12*m22),l13/sqrt(m12*m42)
c>     &      ,l14/sqrt(m12*m42),l23/sqrt(m42*m22)
c>     &      ,l24/sqrt(m42*m22),l34/sqrt(m42*m42)
        
        if (abs(l23*l14).lt.abs(l12*l34)) then
          r14 = l14/(2d0*m42)*(1d0+sqrt(1d0-4d0*m12*m42/l14**2))
          r24 = l24/(2d0*m42)*(1d0+sqrt(1d0-4d0*m22*m42/l24**2))
        else
          r14 = 2d0*m12/(l14*(1d0+sqrt(1d0-4d0*m12*m42/l14**2)))
          r24 = 2d0*m22/(l24*(1d0+sqrt(1d0-4d0*m22*m42/l24**2)))
        end if
        r41 = r14*m42/m12
        r42 = r24*m42/m22

c  for testing
        r24 = 1/r42
        r42 = r24*m42/m22
                
        if(dreal(l12).lt.-0d0) then
          ir12 = sign(1d1,1d0-abs(r12*r12*m22/m12))
        else
          ir12 = 0d0
        end if
        if(dreal(l14).lt.-0d0) then
          ir14 = sign(1d1,1d0-abs(r14*r14*m42/m12))
        else
          ir14 = 0d0
        end if
        if(dreal(l24).lt.-0d0) then
          ir24 = sign(1d1,1d0-abs(r24*r24*m42/m22))
        else
          ir24 = 0d0
        end if
        a   =  m42*(l34/r42 - l23)
        b   =  l13*m22*(1d0/r24-r42) + l12*l34 -l14*l23
        c   =  l13*(l12-r24*l14) - m12*l23 + m12*r24*l34
        d   =  l23 - r24*l34
        det =  sqrt(
     &      l12*l12*l34*l34 + l14*l14*l23*l23 + l24*l24*l13*l13
     &      - 2d0*(l12*l23*l34*l14 + l12*l24*l34*l13 + l14*l24*l23*l13)
     &      + 4d0*(m42*l12*l23*l13 + m12*l23*l34*l24 + m22*l14*l34*l13)
     &      - 4d0*(m12*m42*l23*l23 + m12*m22*l34*l34 + m22*m42*l13*l13)) 
        x4(1) = (-b+det)/(2d0*a)
        x4(2) = (-b-det)/(2d0*a)
        if(abs(x4(1)).gt.abs(x4(2))) then
          x4(2) = c/(a*x4(1))
        else
          x4(1) = c/(a*x4(2))
        end if
        x1(1) = x4(1)/r24
        x1(2) = x4(2)/r24
        ix4(1) = -sign(1d0,dreal(d))
        ix4(2) = +sign(1d0,dreal(d))
        ix1(1) =  sign(1d0,ix4(1)*dreal(r24))
        ix1(2) =  sign(1d0,ix4(2)*dreal(r24))

c>        write(*,*) 'D0m0_coli l ',l12,l13,l14
c>        write(*,*) 'D0m0_coli l ',l23,l24,l34
c>        write(*,*) 'D0m0_coli r ',r12,r14,r24
c>        write(*,*) 'D0m0_coli r ',r21,r41,r42
c>        write(*,*) 'D0m0_coli x4 ',x4(1),x4(2)
c>        write(*,*) 'D0m0_coli x1 ',x1(1),x1(2)
c>        write(*,*) 'D0m0_coli ir ',ir12,ir14,ir24
c>        write(*,*) 'D0m0_coli ix4 ',ix4(1),ix4(2)
c>        write(*,*) 'D0m0_coli ix1 ',ix1(1),ix1(2)

        D0m0_coli = dcmplx(0d0)
        dO i=2,1,-1
          eta=eta2s_coli(-x4(i),1d0/r24,-ix4(i),-ir24,-ix1(i))
c changed 16.10.08
c          iqbard(i) = dreal(-m22*(r42-1d0/r24)*ix4(i)/d)
          iqbard(i) = l12-r24*l14+m22*x4(i)*(1d0/r24-r42)-l23+r24*l34
          if(dreal(x4(i)).gt.0d0)then
            iqbard(i) = dreal(m12/x4(i)+l14+m42*x4(i) - (l34+l13/x4(i)))
          else
            iqbard(i) = dreal(m12/x1(i)+l12+m22*x1(i) - (l23+l13/x1(i)))
          endif

          if(eta.ne.0d0)then
            test0 = abs(d)/max(abs(l23),abs(l34*r24))
            ch1   = l12-r24*l14-m22*(r42-1d0/r24)*x4(i)
            test1 = abs(ch1)/max(abs(l12),abs(m22*r42*x4(i)),
     &          abs(m22*r24*l14),abs(1/r24*x4(i)))
            ch2    = m12/x1(i)+l12+m22*x1(i)
            test01=min(test0,test1)
            test2 = abs(ch2)/max(abs(m12/x1(i)),abs(l12),
     &          abs(m22*x1(i)))
            ch3   = l23+l13/x1(i)
            test3 = abs(ch3)/max(abs(l23),abs(l13/x1(i)))
            test23 = min(test2,test3)
            ch4    = m12/x4(i)+l14+m42*x4(i)
            test4 = abs(ch4)/max(abs(m12/x4(i)),abs(l14),
     &          abs(m42*x4(i)))
            ch5   = l34+l13/x4(i)
            test5 = abs(ch5)/max(abs(l34),abs(l13/x4(i)))
            test45 = min(test4,test5)
            if (test23.gt.test01.and.test23.gt.test45) then
              l1 = cln_coli(ch2/ch3,iqbard)
            elseif(test01.gt.test45) then
              l1 = cln_coli(ch1/d,iqbard)      
            else
              l1 = cln_coli(ch4/ch5,iqbard)      
            end if

c>            write(*,*) 'D0m0_coli l1= ',l1,
c>     &          cln_coli(ch2/ch3,iqbard),
c>     &          cln_coli(ch1/d,iqbard),   
c>     &          cln_coli(ch4/ch5,iqbard)      
c>            write(*,*) 'D0m0_coli l2= ',
c>     &          cln_coli(ch1,-1d0)-cln_coli(d,-1d0),
c>     &          cln_coli(ch2*d/ch3,-1d0)-cln_coli(d,-1d0),
c>     &          cln_coli(ch4*d/ch5,-1d0)-cln_coli(d,-1d0),
c>     &          cln_coli(ch1,-1d0)-cln_coli(ch1*ch3/ch2,-1d0),
c>     &          cln_coli(ch1,-1d0)-cln_coli(ch1*ch5/ch4,-1d0)
c>            write(*,*) 'D0m0_coli l2= ',
c>     &          cln_coli(ch4,-1d0)-cln_coli(ch5,-1d0),
c>     &          cln_coli(ch5*ch1/d,-1d0)-cln_coli(ch5,-1d0),
c>     &          cln_coli(ch5*ch2/ch3,-1d0)-cln_coli(ch5,-1d0),
c>     &          cln_coli(ch4,-1d0)-cln_coli(ch4*d/ch1,-1d0),
c>     &          cln_coli(ch4,-1d0)-cln_coli(ch4*ch3/ch2,-1d0)
c>            write(*,*) 'D0m0_coli l2= ',
c>     &          cln_coli(ch2,-1d0)-cln_coli(ch3,-1d0),
c>     &          cln_coli(ch3*ch1/d,-1d0)-cln_coli(ch3,-1d0),
c>     &          cln_coli(ch3*ch4/ch5,-1d0)-cln_coli(ch3,-1d0),
c>     &          cln_coli(ch2,-1d0)-cln_coli(ch2*d/ch1,-1d0),
c>     &          cln_coli(ch2,-1d0)-cln_coli(ch2*ch5/ch4,-1d0)

#ifdef CHECK
            if(conscheck)then
              if(dimag(r24).ne.0d0) then
                v=dimag(x4(i))/dimag(r24)
                u=dimag(x4(i)/r24)/dimag(1d0/r24)
                l2=cln_coli(m12+l12*v+l14*u+m22*v*v+m42*u*u+l24*u*v,
     &              -1d0)
     &              -cln_coli(l13+l23*v+l34*u,-1d0)

c                write(*,*) 'D0m0_coli comp l2= ',l2,l1,l2/l1

              elseif(dimag(x4(i)).ne.0d0) then
                if(test1.lt.test0)then
                  if(test1.gt.test23.and.test1.gt.test45)then
                    l2=cln_coli(ch1,-1d0)-cln_coli(d,-1d0)
                  elseif(test23.gt.test45)then
                    l2=cln_coli(ch2*d/ch3,-1d0)-cln_coli(d,-1d0)
                  else
                    l2=cln_coli(ch4*d/ch5,-1d0)-cln_coli(d,-1d0)
                  endif
                else
                  if(test1.gt.test23.and.test1.gt.test45)then
                    l2=cln_coli(l12-r24*l14+m42*x4(i)*(1d0/r42-r24),
     &                  -1d0)
     &                  -cln_coli(l23-r24*l34,-1d0)
                  elseif(test23.gt.test45)then
                    l2=cln_coli(ch1,-1d0)-cln_coli(ch1*ch3/ch2,-1d0)
                  else
                    l2=cln_coli(ch1,-1d0)-cln_coli(ch1*ch5/ch4,-1d0)
                  endif
                endif
              elseif(dreal(x4(i)).gt.0d0)then
                if(test4.lt.test5)then
                  if(test4.gt.test23.and.test4.gt.test01)then
                    l2=cln_coli(ch4,-1d0)-cln_coli(ch5,-1d0)
                  elseif(test23.lt.test01)then
                    l2=cln_coli(ch5*ch1/d,-1d0)-cln_coli(ch5,-1d0)
                  else
                    l2=cln_coli(ch5*ch2/ch3,-1d0)-cln_coli(ch5,-1d0)
                  endif
                else
                  if(test4.gt.test23.and.test4.gt.test01)then
                    l2=cln_coli(ch4,-1d0)-cln_coli(ch5,-1d0)
                  elseif(test23.lt.test01)then
                    l2=cln_coli(ch4,-1d0)-cln_coli(ch4*d/ch1,-1d0)
                  else
                    l2=cln_coli(ch4,-1d0)-cln_coli(ch4*ch3/ch2,-1d0)
                  endif
                endif
              elseif(dreal(x1(i)).gt.0d0)then
                if(test2.lt.test3)then
                  if(test2.gt.test45.and.test2.gt.test01)then
                    l2=cln_coli(ch2,-1d0)-cln_coli(ch3,-1d0)
                  elseif(test45.lt.test01)then
                    l2=cln_coli(ch3*ch1/d,-1d0)-cln_coli(ch3,-1d0)
                  else
                    l2=cln_coli(ch3*ch4/ch5,-1d0)-cln_coli(ch3,-1d0)
                  endif
                else
                  if(test2.gt.test45.and.test2.gt.test01)then
                    l2=cln_coli(ch2,-1d0)-cln_coli(ch3,-1d0)
                  elseif(test45.lt.test01)then
                    l2=cln_coli(ch2,-1d0)-cln_coli(ch2*d/ch1,-1d0)
                  else
                    l2=cln_coli(ch2,-1d0)-cln_coli(ch2*ch5/ch4,-1d0)
                  endif
                endif
              else
                if (flag3) then
                  write(*,*) 
                  write(*,*) 'D0m0_coli: case for l2 not implemented '
                  flag3 = .false.
                end if
                write(errout,*) 
                write(errout,*) 
     &              'D0m0_coli: case for l2 not implemented '
                write(errout,*) 
     &              'D0m0_coli args ',p12,p23,p34,p14,p13,p24,
     &              m12,m22,m32,m42
                write(errout,*) 'r24 = ',r24
                write(errout,*) 'x4(i) = ',x4(i)
                write(errout,*) 'x1(i) = ',x1(i)
              endif 
             
c              write(*,*) 'D0m0_coli l2= ',l2,l1,l2/l1

              if(abs(l2/l1-1d0).gt.testacc) then
                if (flag2) then
                  write(*,*) 
                  write(*,*) 'D0m0_coli: eta terms possibly wrong'
                  flag2 = .false.
                end if
                write(errout,*) 
                write(errout,*) 'D0m0_coli: eta terms  possibly wrong'
                write(errout,*) 
     &              'D0m0_coli args ',p12,p23,p34,p14,p13,p24,
     &              m12,m22,m32,m42
                write(errout,*) 'l1 = ',l1
                write(errout,*) 'l2 = ',l2
                write(errout,*) 'diff   = ',abs((l2-l1)/l1)
              endif 
            endif
#endif
          endif

          if (l13.ne.0d0) then

#ifdef ALLCHECK
c ppttbb dim
            if(flag(1))then
              write(*,100)' D0m0_coli: ',
     &            '    general case',
     &            '    not yet tested in physical process '
              write(*,111)' D0m0_coli: p12 = ',p12
              write(*,111)' D0m0_coli: p23 = ',p23
              write(*,111)' D0m0_coli: p34 = ',p34
              write(*,111)' D0m0_coli: p14 = ',p14
              write(*,111)' D0m0_coli: p24 = ',p24
              write(*,111)' D0m0_coli: p13 = ',p13
              write(*,111)' D0m0_coli: m12 = ',m12
              write(*,111)' D0m0_coli: m22 = ',m22
              write(*,111)' D0m0_coli: m32 = ',m32
              write(*,111)' D0m0_coli: m42 = ',m42
              flag(1)=.false.
            endif
#endif

            D0m0_coli = D0m0_coli + (2*i-3) * (
     &           cspcos_coli(-x4(i),r41,-ix4(i),ir14)
     &           + cspcos_coli(-x4(i),1d0/r14,-ix4(i),-ir14)
     &           - cspcos_coli(-x4(i),l34/l13,-ix4(i),dreal(l34-l13))
     &           - cspcos_coli(-x1(i),r21,-ix1(i),ir12)
     &           - cspcos_coli(-x1(i),1d0/r12,-ix1(i),-ir12)
     &           + cspcos_coli(-x1(i),l23/l13,-ix1(i),dreal(l23-l13))
     &           - eta*(l1+cln_coli(l13/m12,-1d0)) 
     &           )

c>            write(*,*) 'D0m0_coli contr ',
c>     &       (2*i-3) *cspcos_coli(-x4(i),r41,-ix4(i),ir14)
c>            write(*,*) 'D0m0_coli contr ',
c>     &       + (2*i-3) *cspcos_coli(-x4(i),1d0/r14,-ix4(i),-ir14)
c>            write(*,*) 'D0m0_coli contr ',
c>     &       - (2*i-3) *cspcos_coli(-x1(i),r21,-ix1(i),ir12)
c>            write(*,*) 'D0m0_coli contr ',
c>     &       - (2*i-3) *cspcos_coli(-x1(i),1d0/r12,-ix1(i),-ir12)
c>            write(*,*) 'D0m0_coli contr ',- (2*i-3) *
c>     &       cspcos_coli(-x4(i),l34/l13,-ix4(i),dreal(l34-l13))
c>            write(*,*) 'D0m0_coli contr ',+ (2*i-3) *
c>     &       cspcos_coli(-x1(i),l23/l13,-ix1(i),dreal(l23-l13))
c>            write(*,*) 'D0m0_coli contr ',- (2*i-3) *
c>     &       eta2s_coli(-x4(i),1d0/r24,-ix4(i),-ir24,-ix1(i))*
c>     &       (l1+cln_coli(l13/m12,-1d0)) 
c>            write(*,*) 'D0m0_coli sum ',D0m0_coli

          else

#ifdef ALLCHECK
c ppttbb dim
            if(flag(2))then
              write(*,100)' D0m0_coli: ',
     &            '    case l13=0',
     &            '    not yet tested in physical process'
              write(*,111)' D0m0_coli: p12 = ',p12
              write(*,111)' D0m0_coli: p23 = ',p23
              write(*,111)' D0m0_coli: p34 = ',p34
              write(*,111)' D0m0_coli: p14 = ',p14
              write(*,111)' D0m0_coli: p24 = ',p24
              write(*,111)' D0m0_coli: p13 = ',p13
              write(*,111)' D0m0_coli: m12 = ',m12
              write(*,111)' D0m0_coli: m22 = ',m22
              write(*,111)' D0m0_coli: m32 = ',m32
              write(*,111)' D0m0_coli: m42 = ',m42
              flag(2)=.false.
            endif
#endif

            D0m0_coli = D0m0_coli + (2*i-3) * (
     &           cspcos_coli(-x4(i),r41,-ix4(i),ir14)
     &           + cspcos_coli(-x4(i),1d0/r14,-ix4(i),-ir14)
     &           - cspcos_coli(-x1(i),r21,-ix1(i),ir12)
     &           - cspcos_coli(-x1(i),1d0/r12,-ix1(i),-ir12)
     &           - eta2s_coli(-x4(i),1d0/r24,-ix4(i),-ir24,-ix1(i))* l1
     &           + (cln_coli(-x4(i),-ix4(i)))**2/2d0
     &           - (cln_coli(-x1(i),-ix1(i)))**2/2d0
     &           + cln_coli(-x4(i),-ix4(i))*cln_coli(l34/m12,-1d0)
     &           - cln_coli(-x1(i),-ix1(i))*cln_coli(l23/m12,-1d0)
     &           )

c>            write(*,*) 'D0m0_coli contr ',
c>     &       (2*i-3) *  cspcos_coli(-x4(i),r41,-ix4(i),ir14)
c>            write(*,*) 'D0m0_coli contr ',
c>     &       (2*i-3) * cspcos_coli(-x4(i),1d0/r14,-ix4(i),-ir14)
c>            write(*,*) 'D0m0_coli contr ',
c>     &       -(2*i-3) * cspcos_coli(-x1(i),r21,-ix1(i),ir12)
c>            write(*,*) 'D0m0_coli contr ',
c>     &       -(2*i-3) * cspcos_coli(-x1(i),1d0/r12,-ix1(i),-ir12)
c>            write(*,*) 'D0m0_coli contr ',
c>     &       -(2*i-3)*eta2s_coli(-x4(i),1d0/r24,-ix4(i),-ir24,-ix1(i))
c>     &          * l1
c>            write(*,*) 'D0m0_coli contr ',
c>     &       (2*i-3) * (cln_coli(-x4(i),-ix4(i)))**2/2d0
c>            write(*,*) 'D0m0_coli contr ',
c>     &       -(2*i-3) * (cln_coli(-x1(i),-ix1(i)))**2/2d0
c>            write(*,*) 'D0m0_coli contr ',
c>     &       (2*i-3) *  cln_coli(-x4(i),-ix4(i))*cln_coli(l34/m12,-1d0)
c>            write(*,*) 'D0m0_coli contr ',
c>     &       -(2*i-3) * cln_coli(-x1(i),-ix1(i))*cln_coli(l23/m12,-1d0)

          end if
        enddo

        D0m0_coli = D0m0_coli/det

c        write(*,*) 'D0m0_coli =  ',D0m0_coli


#ifdef CHECK

       if (conscheck) then

         do j=1,3
           if(j.eq.1) then
             ix4(1) = +sign(1d0,dreal(d))
             ix4(2) = -sign(1d0,dreal(d))
           elseif(j.eq.2) then
             ix4(1) = +sign(1d0,dreal(d))
             ix4(2) = +sign(1d0,dreal(d))
           elseif(j.eq.3) then
             ix4(1) = -sign(1d0,dreal(d))
             ix4(2) = -sign(1d0,dreal(d))
           endif
           ix1(1) =  sign(1d0,ix4(1)*dreal(r24))  ! needed for correct
           ix1(2) =  sign(1d0,ix4(2)*dreal(r24))  ! continuation

c>           write(*,*) 'D0m0_check l ',l12,l13,l14
c>           write(*,*) 'D0m0_check l ',l23,l24,l34
c>           write(*,*) 'D0m0_check r ',r12,r14,r24
c>           write(*,*) 'D0m0_check r ',r21,r41,r42
c>           write(*,*) 'D0m0_check x4 ',x4(1),x4(2)
c>           write(*,*) 'D0m0_check x1 ',x1(1),x1(2)
c>           write(*,*) 'D0m0_check ir ',ir12,ir14,ir24
c>           write(*,*) 'D0m0_check ix4 ',ix4(1),ix4(2)
c>           write(*,*) 'D0m0_check ix1 ',ix1(1),ix1(2)

            D0m0_check = dcmplx(0d0)
           do i=2,1,-1
             eta=eta2s_coli(-x4(i),1d0/r24,-ix4(i),-ir24,-ix1(i))
c changed 16.10.08
c             iqbard(i) = dreal(-m22*(r42-1d0/r24)*ix4(i)/d)
             iqbard(i) = l12-r24*l14+m22*x4(i)*(1d0/r24-r42)-l23+r24*l34
             if(dreal(x4(i)).gt.0d0)then
               iqbard(i)=dreal(m12/x4(i)+l14+m42*x4(i)-(l34+l13/x4(i)))
             else
               iqbard(i)=dreal(m12/x1(i)+l12+m22*x1(i)-(l23+l13/x1(i)))
             endif
             if(eta.ne.0d0)then
               test0 = abs(d)/max(abs(l23),abs(l34*r24))
               ch1   = l12-r24*l14-m22*(r42-1d0/r24)*x4(i)
               test1 = abs(ch1)/max(abs(l12),abs(r42*x4(i)),
     &             abs(r24*l14),abs(1/r24*x4(i)))
               ch2    = m12/x1(i)+l12+m22*x1(i)
               test01=min(test0,test1)
               test2 = abs(ch2)/max(abs(m12/x1(i)),abs(l12),
     &             abs(m22*x1(i)))
               ch3   = l23+l13/x1(i)
               test3 = abs(ch3)/max(abs(l23),abs(l13/x1(i)))
               test23 = min(test2,test3)
               ch4    = m12/x4(i)+l14+m42*x4(i)
               test4 = abs(ch4)/max(abs(m12/x4(i)),abs(l14),
     &             abs(m42*x4(i)))
               ch5   = l34+l13/x4(i)
               test5 = abs(ch5)/max(abs(l34),abs(l13/x4(i)))
               test45 = min(test4,test5)
               if (test23.gt.test01.and.test23.gt.test45) then
                 l1=cln_coli(ch2/ch3,iqbard)
               elseif(test01.gt.test45) then
                 l1 = cln_coli(ch1/d,iqbard)      
               else
                 l1=cln_coli(ch4/ch5,iqbard)      
               end if

c>               write(*,*) 'D0m0_check l1= ',l1,
c>     &             cln_coli(ch2/ch3,iqbard),
c>     &             cln_coli(ch1/d,iqbard),   
c>     &             cln_coli(ch4/ch5,iqbard)      
c>               write(*,*) 'D0m0_check l2= ',
c>     &             cln_coli(ch1,-1d0)-cln_coli(d,-1d0),
c>     &             cln_coli(ch2*d/ch3,-1d0)-cln_coli(d,-1d0),
c>     &             cln_coli(ch4*d/ch5,-1d0)-cln_coli(d,-1d0),
c>     &             cln_coli(ch1,-1d0)-cln_coli(ch1*ch3/ch2,-1d0),
c>     &             cln_coli(ch1,-1d0)-cln_coli(ch1*ch5/ch4,-1d0)
c>               write(*,*) 'D0m0_check l2= ',
c>     &             cln_coli(ch4,-1d0)-cln_coli(ch5,-1d0),
c>     &             cln_coli(ch5*ch1/d,-1d0)-cln_coli(ch5,-1d0),
c>     &             cln_coli(ch5*ch2/ch3,-1d0)-cln_coli(ch5,-1d0),
c>     &             cln_coli(ch4,-1d0)-cln_coli(ch4*d/ch1,-1d0),
c>     &             cln_coli(ch4,-1d0)-cln_coli(ch4*ch3/ch2,-1d0)
c>               write(*,*) 'D0m0_check l2= ',
c>     &             cln_coli(ch2,-1d0)-cln_coli(ch3,-1d0),
c>     &             cln_coli(ch3*ch1/d,-1d0)-cln_coli(ch3,-1d0),
c>     &             cln_coli(ch3*ch4/ch5,-1d0)-cln_coli(ch3,-1d0),
c>     &             cln_coli(ch2,-1d0)-cln_coli(ch2*d/ch1,-1d0),
c>     &             cln_coli(ch2,-1d0)-cln_coli(ch2*ch5/ch4,-1d0)

               if(dimag(r24).ne.0d0) then
                 v=dimag(x4(i))/dimag(r24)
                 u=dimag(x4(i)/r24)/dimag(1d0/r24)
                 l2=cln_coli(m12+l12*v+l14*u+m22*v*v+m42*u*u+l24*u*v,
     &               -1d0)
     &               -cln_coli(l13+l23*v+l34*u,-1d0)
                 
c                 write(*,*) 'D0m0_check comp l2= ',l2,l1,l2/l1

               elseif(dimag(x4(i)).ne.0d0) then
                 if(test1.lt.test0)then
                   if(test1.gt.test23.and.test1.gt.test45)then
                     l2=cln_coli(ch1,-1d0)-cln_coli(d,-1d0)
                   elseif(test23.gt.test45)then
                     l2=cln_coli(ch2*d/ch3,-1d0)-cln_coli(d,-1d0)
                   else
                     l2=cln_coli(ch4*d/ch5,-1d0)-cln_coli(d,-1d0)
                   endif
                 else
                   if(test1.gt.test23.and.test1.gt.test45)then
                     l2=cln_coli(l12-r24*l14+m42*x4(i)*(1d0/r42-r24),
     &                   -1d0)
     &                   -cln_coli(l23-r24*l34,-1d0)
                   elseif(test23.gt.test45)then
                     l2=cln_coli(ch1,-1d0)-cln_coli(ch1*ch3/ch2,-1d0)
                   else
                     l2=cln_coli(ch1,-1d0)-cln_coli(ch1*ch5/ch4,-1d0)
                   endif
                 endif
               elseif(dreal(x4(i)).gt.0d0)then
                 if(test4.lt.test5)then
                   if(test4.gt.test23.and.test4.gt.test01)then
                     l2=cln_coli(ch4,-1d0)-cln_coli(ch5,-1d0)
                   elseif(test23.lt.test01)then
                     l2=cln_coli(ch5*ch1/d,-1d0)-cln_coli(ch5,-1d0)
                   else
                     l2=cln_coli(ch5*ch2/ch3,-1d0)-cln_coli(ch5,-1d0)
                   endif
                 else
                   if(test4.gt.test23.and.test4.gt.test01)then
                     l2=cln_coli(ch4,-1d0)-cln_coli(ch5,-1d0)
                   elseif(test23.lt.test01)then
                     l2=cln_coli(ch4,-1d0)-cln_coli(ch4*d/ch1,-1d0)
                   else
                     l2=cln_coli(ch4,-1d0)-cln_coli(ch4*ch3/ch2,-1d0)
                   endif
                 endif
               elseif(dreal(x1(i)).gt.0d0)then
                 if(test2.lt.test3)then
                   if(test2.gt.test45.and.test2.gt.test01)then
                     l2=cln_coli(ch2,-1d0)-cln_coli(ch3,-1d0)
                   elseif(test45.lt.test01)then
                     l2=cln_coli(ch3*ch1/d,-1d0)-cln_coli(ch3,-1d0)
                   else
                     l2=cln_coli(ch3*ch4/ch5,-1d0)-cln_coli(ch3,-1d0)
                   endif
                 else
                   if(test2.gt.test45.and.test2.gt.test01)then
                     l2=cln_coli(ch2,-1d0)-cln_coli(ch3,-1d0)
                   elseif(test45.lt.test01)then
                     l2=cln_coli(ch2,-1d0)-cln_coli(ch2*d/ch1,-1d0)
                   else
                     l2=cln_coli(ch2,-1d0)-cln_coli(ch2*ch5/ch4,-1d0)
                   endif
                 endif
               else
                 if (flag3) then
                   write(*,*) 
                   write(*,*) 'D0m0_check: case for l2 not implemented '
                   flag3 = .false.
                 end if
                 write(errout,*) 
                 write(errout,*) 
     &               'D0m0_check: case for l2 not implemented '
                 write(errout,*) 
     &               'D0m0_check args ',p12,p23,p34,p14,p13,p24,
     &               m12,m22,m32,m42
                 write(errout,*) 'r24 = ',r24
                 write(errout,*) 'x4(i) = ',x4(i)
                 write(errout,*) 'x1(i) = ',x1(i)
               endif 
               
c     write(*,*) 'D0m0_check l2= ',l2,l1,l2/l1

               if(abs(l2/l1-1d0).gt.testacc) then
                 if (flag2) then
                   write(*,*) 
                   write(*,*) 'D0m0_check: eta terms possibly wrong ',i
                   flag2 = .false.
                 end if
                 write(errout,*) 
                 write(errout,*) 
     &               'D0m0_check: eta terms  possibly wrong',i
                 write(errout,*) 
     &               'D0m0_check args ',p12,p23,p34,p14,p13,p24,
     &               m12,m22,m32,m42
                 write(errout,*) 'l1 = ',l1
                 write(errout,*) 'l2 = ',l2
                 write(errout,*) 'diff   = ',abs((l2-l1)/l1)
               endif 
             endif

             if (l13.ne.0d0) then

               D0m0_check = D0m0_check + (2*i-3) * (
     &             cspcos_coli(-x4(i),r41,-ix4(i),ir14)
     &             + cspcos_coli(-x4(i),1d0/r14,-ix4(i),-ir14)
     &             - cspcos_coli(-x4(i),l34/l13,-ix4(i),dreal(l34-l13))
     &             - cspcos_coli(-x1(i),r21,-ix1(i),ir12)
     &             - cspcos_coli(-x1(i),1d0/r12,-ix1(i),-ir12)
     &             + cspcos_coli(-x1(i),l23/l13,-ix1(i),dreal(l23-l13))
     &             - eta*(l1+cln_coli(l13/m12,-1d0)) 
     &             )

c>            write(*,*) 'D0m0_check contr ',
c>     &       (2*i-3) *cspcos_coli(-x4(i),r41,-ix4(i),ir14)
c>            write(*,*) 'D0m0_check contr ',
c>     &       + (2*i-3) *cspcos_coli(-x4(i),1d0/r14,-ix4(i),-ir14)
c>            write(*,*) 'D0m0_check contr ',
c>     &       - (2*i-3) *cspcos_coli(-x1(i),r21,-ix1(i),ir12)
c>            write(*,*) 'D0m0_check contr ',
c>     &       - (2*i-3) *cspcos_coli(-x1(i),1d0/r12,-ix1(i),-ir12)
c>            write(*,*) 'D0m0_check contr ',- (2*i-3) *
c>     &       cspcos_coli(-x4(i),l34/l13,-ix4(i),dreal(l34-l13))
c>            write(*,*) 'D0m0_check contr ', (2*i-3) *
c>     &       cspcos_coli(-x1(i),l23/l13,-ix1(i),dreal(l23-l13))
c>            write(*,*) 'D0m0_check contr ',- (2*i-3) *
c>     &       eta2s_coli(-x4(i),1d0/r24,-ix4(i),-ir24,-ix1(i))*
c>     &       (l1+cln_coli(l13/m12,-1d0)) 
c>            write(*,*) 'D0m0_check sum ',D0m0_check

             else
               
               D0m0_check = D0m0_check + (2*i-3) * (
     &             cspcos_coli(-x4(i),r41,-ix4(i),ir14)
     &             + cspcos_coli(-x4(i),1d0/r14,-ix4(i),-ir14)
     &             - cspcos_coli(-x1(i),r21,-ix1(i),ir12)
     &             - cspcos_coli(-x1(i),1d0/r12,-ix1(i),-ir12)
     &             - eta2s_coli(-x4(i),1d0/r24,-ix4(i),-ir24,-ix1(i))*l1
     &             + (cln_coli(-x4(i),-ix4(i)))**2/2d0
     &             - (cln_coli(-x1(i),-ix1(i)))**2/2d0
     &             + cln_coli(-x4(i),-ix4(i))*cln_coli(l34/m12,-1d0)
     &             - cln_coli(-x1(i),-ix1(i))*cln_coli(l23/m12,-1d0)
     &             )

c>            write(*,*) 'D0m0_check contr ',
c>     &       (2*i-3) *  cspcos_coli(-x4(i),r41,-ix4(i),ir14)
c>            write(*,*) 'D0m0_check contr ',
c>     &       (2*i-3) * cspcos_coli(-x4(i),1d0/r14,-ix4(i),-ir14)
c>            write(*,*) 'D0m0_check contr ',
c>     &       -(2*i-3) * cspcos_coli(-x1(i),r21,-ix1(i),ir12)
c>            write(*,*) 'D0m0_check contr ',
c>     &       -(2*i-3) * cspcos_coli(-x1(i),1d0/r12,-ix1(i),-ir12)
c>            write(*,*) 'D0m0_check contr ',-(2*i-3)*
c>     &       eta2s_coli(-x4(i),1d0/r24,-ix4(i),-ir24,-ix1(i))* l1
c>            write(*,*) 'D0m0_check contr ',
c>     &       (2*i-3) * (cln_coli(-x4(i),-ix4(i)))**2/2d0
c>            write(*,*) 'D0m0_check contr ',
c>     &       -(2*i-3) * (cln_coli(-x1(i),-ix1(i)))**2/2d0
c>            write(*,*) 'D0m0_check contr ',
c>     &       (2*i-3) *  cln_coli(-x4(i),-ix4(i))*cln_coli(l34/m12,-1d0)
c>            write(*,*) 'D0m0_check contr ',
c>     &       -(2*i-3) * cln_coli(-x1(i),-ix1(i))*cln_coli(l23/m12,-1d0)

             end if
           enddo

c           write(*,*) 'D0m0_check sum ',D0m0_check

           D0m0_check = D0m0_check/det

c           write(*,*) 'D0m0_check =  ',D0m0_check

           if(cdabs(D0m0_coli/D0m0_check-1d0).gt.testacc) then
             write(testout,*) ' D0m0_coli inconsistent ',i
             write(testout,111) ' D0m0_coli: p12 = ',p12
             write(testout,111) ' D0m0_coli: p23 = ',p23
             write(testout,111) ' D0m0_coli: p34 = ',p34
             write(testout,111) ' D0m0_coli: p14 = ',p14
             write(testout,111) ' D0m0_coli: p24 = ',p24
             write(testout,111) ' D0m0_coli: p13 = ',p13
             write(testout,111) ' D0m0_coli: m12 = ',m12
             write(testout,111) ' D0m0_coli: m22 = ',m22
             write(testout,111) ' D0m0_coli: m32 = ',m32
             write(testout,111) ' D0m0_coli: m42 = ',m42
             write(testout,111) ' D0m0_coli  = ',D0m0_coli
             write(testout,111) ' D0m0_check = ',D0m0_check
             write(testout,101) ' diff        = ',
     &           cdabs(D0m0_coli/D0m0_check-1d0)
             
             write(*,*)   ' D0m0_coli inconsistent ',i
             write(*,111) ' D0m0_coli: p12 = ',p12
             write(*,111) ' D0m0_coli: p23 = ',p23
             write(*,111) ' D0m0_coli: p34 = ',p34
             write(*,111) ' D0m0_coli: p14 = ',p14
             write(*,111) ' D0m0_coli: p24 = ',p24
             write(*,111) ' D0m0_coli: p13 = ',p13
             write(*,111) ' D0m0_coli: m12 = ',m12
             write(*,111) ' D0m0_coli: m22 = ',m22
             write(*,111) ' D0m0_coli: m32 = ',m32
             write(*,111) ' D0m0_coli: m42 = ',m42
             write(*,111) ' D0m0_coli  = ',D0m0_coli
             write(*,111) ' D0m0_check = ',D0m0_check
             write(*,101) ' diff        = ',
     &           cdabs(D0m0_coli/D0m0_check-1d0)
             stop
           endif
           
         enddo

         m2  = sqrt(m22)
         m1  = sqrt(m12)
         m4  = sqrt(m42)
         mm3 = m4

         v4 = m1/m4
         v1 = m1/m2
         
c     (permutation 1 -> 3 -> 2 -> 4 -> 1 for  m -> k )
c     m1 -> m2 in normalization
         k12 = (m22+m12-p12)/(m1*m2)
         k13 = (m12    -p13)/(m1*m4)
         k14 = (m42+m12-p14)/(m1*m4)
         k23 = (m22    -p23)/(m2*m4)
         k24 = (m22+m42-p24)/(m2*m4)
         k34 = (    m42-p34)/(m42)
         r12 = k12/2d0*(1d0+sqrt(1d0-4d0/k12**2))
         
         if (k23.eq.0.and.k34.eq.0d0) then
           write(*,*) 
     &         'D0m0_coli called with improper arguments: k23=k34=0'
           write(*,*) 'ir-singular case not implemented in D0m0_coli'
         else if (abs(k23*k14).lt.abs(k12*k34)) then
           r14 = k14/2d0*(1d0+sqrt(1d0-4d0/k14**2))
           r24 = k24/2d0*(1d0+sqrt(1d0-4d0/k24**2))
         else
           r14 = 2d0/(k14*(1d0+sqrt(1d0-4d0/k14**2)))
           r24 = 2d0/(k24*(1d0+sqrt(1d0-4d0/k24**2)))
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
         if(dreal(k24).lt.-2d0) then
           ir24 = sign(1d1,1d0-abs(r24))
         else
           ir24 = 0d0
         end if
         a   =  k34/r24 - k23
         b   =  k13*(1d0/r24-r24) + k12*k34 -k14*k23
         c   =  k13*(k12-r24*k14) - k23 + r24*k34
         d   =  k23 - r24*k34
         det =  sqrt(
     &       k12*k12*k34*k34 + k14*k14*k23*k23 + k24*k24*k13*k13
     &       - 2d0*(k12*k23*k34*k14 + k12*k24*k34*k13 + k14*k24*k23*k13)
     &       + 4d0*(k12*k23*k13 + k23*k34*k24 + k14*k34*k13)
     &       - 4d0*(k23*k23 + k34*k34 + k13*k13)) 
         x4(1) = (-b+det)/(2d0*a)
         x4(2) = (-b-det)/(2d0*a)
         if(abs(x4(1)).gt.abs(x4(2))) then
           x4(2) = c/(a*x4(1))
         else
           x4(1) = c/(a*x4(2))
         end if
         x1(1) = x4(1)/r24
         x1(2) = x4(2)/r24
         ix4(1) = -sign(1d0,dreal(d))
         ix4(2) = +sign(1d0,dreal(d))
         ix1(1) =  sign(1d0,ix4(1)*dreal(r24))
         ix1(2) =  sign(1d0,ix4(2)*dreal(r24))

c>        write(*,*) 'D0m0_check k ',k12,k13,k14
c>        write(*,*) 'D0m0_check k ',k23,k24,k34
c>        write(*,*) 'D0m0_check r ',r12,r14,r24
c>        write(*,*) 'D0m0_check x4 ',x4(1),x4(2)
c>        write(*,*) 'D0m0_check x1 ',x1(1),x1(2)
       
         D0m0_check = dcmplx(0d0)
         dO 21 i=2,1,-1
           eta=eta2s_coli(-x4(i)*v4,1d0/r24*v1/v4,-ix4(i),-ir24,-ix1(i))
           if(eta.ne.0d0) then
             ch1   = k12-r24*k14-(r24-1d0/r24)*x4(i)
             test1 = abs(ch1)/max(abs(k12),abs(r24*x4(i)),
     &           abs(r24*k14),abs(1/r24*x4(i)))
             ch2    = r24/x4(i)+k12+x4(i)/r24
             test2 = abs(ch2)/max(abs(r24/x4(i)),abs(k12),
     &           abs(x4(i)/r24))
             ch3   = k23+k13*r24/x4(i)
             test3 = abs(ch3)/max(abs(k23),abs(k13*r24/x4(i)))
             test2 = min(test2,test3)
            
             ch4    = 1d0/x4(i)+k14+x4(i)
             test4 = abs(ch4)/max(abs(1d0/x4(i)),abs(k14),
     &           abs(x4(i)))
             ch5   = k34+k13/x4(i)
             test5 = abs(ch5)/max(abs(k34),abs(k13/x4(i)))
             test4 = min(test4,test5)
             if (test2.gt.test1.and.test2.gt.test4) then
               l1 = cln_coli(ch2*m1/(ch3*mm3),
     &             iqbard)
             elseif(test1.gt.test4) then
               l1 = cln_coli(ch1*m1/(d*mm3),
     &             iqbard)      
             else
               l1 = cln_coli(ch4*m1/(ch5*mm3),
     &             iqbard)      
             end if
           endif
c     <---
           if (k13.ne.0d0) then
             D0m0_check = D0m0_check + (2*i-3) * (
     &           cspcos_coli(-x4(i)*v4,r14/v4,-ix4(i),ir14)
     &           + cspcos_coli(-x4(i)*v4,1d0/r14/v4,-ix4(i),-ir14)
     &           - cspcos_coli(-x1(i)*v1,r12/v1,-ix1(i),ir12)
     &           - cspcos_coli(-x1(i)*v1,1d0/r12/v1,-ix1(i),-ir12)
     &           - cspcos_coli(-x4(i)*v4,k34/k13/v4,-ix4(i),
     &           dreal((k34*m4-k13*m1)*mm3))
     &           + cspcos_coli(-x1(i)*v1,k23/k13/v1,-ix1(i),
     &           dreal((k23*m2-k13*m1)*mm3))
     &           - eta2s_coli(-x4(i)*v4,1d0/r24*v1/v4,-ix4(i),
     &           -ir24,-ix1(i))*
     &           (l1+cln_coli(k13*mm3/m1,-1d0)) 
     &           )

c>            write(*,*) 'D0m0_check contr ',
c>     &       (2*i-3) *cspcos_coli(-x4(i)*v4,r14/v4,-ix4(i),ir14)
c>            write(*,*) 'D0m0_check contr ',
c>     &       + (2*i-3) *cspcos_coli(-x4(i)*v4,1d0/r14/v4,-ix4(i),-ir14)
c>            write(*,*) 'D0m0_check contr ',
c>     &       - (2*i-3) *cspcos_coli(-x1(i)*v1,r12/v1,-ix1(i),ir12)
c>            write(*,*) 'D0m0_check contr ',
c>     &       - (2*i-3) *cspcos_coli(-x1(i)*v1,1d0/r12/v1,-ix1(i),-ir12)
c>            write(*,*) 'D0m0_check contr ',
c>     &   - (2*i-3) *cspcos_coli(-x4(i)*v4,k34/k13/v4,-ix4(i),
c>     &               dreal((k34*m4-k13*m1)*mm3))
c>            write(*,*) 'D0m0_check contr ',
c>     &   + (2*i-3) *cspcos_coli(-x1(i)*v1,k23/k13/v1,-ix1(i),
c>     &               dreal((k23*m2-k13*m1)*mm3))
c>            write(*,*) 'D0m0_check contr ',
c>     &          - (2*i-3) *
c>     &          eta2s_coli(-x4(i)*v4,1d0/r24*v1/v4,-ix4(i),-ir24,-ix1(i))*
c>     &       (l1+cln_coli(k13*mm3/m1,-1d0)) 

           else
             D0m0_check = D0m0_check + (2*i-3) * (
     &           cspcos_coli(-x4(i)*v4,r14/v4,-ix4(i),ir14)
     &           + cspcos_coli(-x4(i)*v4,1d0/r14/v4,-ix4(i),-ir14)
     &           - cspcos_coli(-x1(i)*v1,r12/v1,-ix1(i),ir12)
     &           - cspcos_coli(-x1(i)*v1,1d0/r12/v1,-ix1(i),-ir12)
     &           -eta2s_coli(-x4(i)*v4,1d0/r24/v4,-ix4(i),-ir24,-ix1(i))
     &           * l1
     &           + (cln_coli(-x4(i)*v4,-ix4(i)))**2/2d0
     &           - (cln_coli(-x1(i)*v1,-ix1(i)))**2/2d0
     &           + cln_coli(-x4(i)*v4,-ix4(i))
     &           *cln_coli(k34*mm3*m4/m12,-1d0)
     &           - cln_coli(-x1(i)*v1,-ix1(i))
     &           *cln_coli(k23*m2*mm3/m12,-1d0)
     &           )

c>            write(*,*) 'D0m0_check contr ',
c>     &     (2*i-3) * cspcos_coli(-x4(i)*v4,r14/v4,-ix4(i),ir14)
c>            write(*,*) 'D0m0_check contr ',
c>     &     (2*i-3) * cspcos_coli(-x4(i)*v4,1d0/r14/v4,-ix4(i),-ir14)
c>            write(*,*) 'D0m0_check contr ',
c>     &   - (2*i-3) * cspcos_coli(-x1(i)*v1,r12/v1,-ix1(i),ir12)
c>            write(*,*) 'D0m0_check contr ',
c>     &   - (2*i-3) * cspcos_coli(-x1(i)*v1,1d0/r12/v1,-ix1(i),-ir12)
c>            write(*,*) 'D0m0_check contr ',
c>     &   -(2*i-3) * eta2s_coli(-x4(i)*v4,1d0/r24/v4,-ix4(i),-ir24,-ix1(i))*l1
c>            write(*,*) 'D0m0_check contr ',
c>     &     (2*i-3) *(cln_coli(-x4(i)*v4,-ix4(i)))**2/2d0
c>            write(*,*) 'D0m0_check contr ',
c>     &   - (2*i-3) * (cln_coli(-x1(i)*v1,-ix1(i)))**2/2d0
c>            write(*,*) 'D0m0_check contr ',
c>     &    (2*i-3) * cln_coli(-x4(i)*v4,-ix4(i))*cln_coli(k34*mm3*m4/m12,-1d0)
c>            write(*,*) 'D0m0_check contr ',
c>     &   - (2*i-3) * cln_coli(-x1(i)*v1,-ix1(i))*cln_coli(k23*m2*mm3/m12,-1d0)

           end if

 21      continue

c        write(*,*) 'D0m0_check sum ',D0m0_check

         D0m0_check = D0m0_check/(det*m1*m2*mm3*m4)
c        write(*,*) 'D0m0_check =  ',D0m0_check

         if(cdabs(D0m0_coli/D0m0_check-1d0).gt.testacc) then
           write(testout,100) ' D0m0_coli inconsistent 4'
           write(testout,111) ' D0m0_coli: p12 = ',p12
           write(testout,111) ' D0m0_coli: p23 = ',p23
           write(testout,111) ' D0m0_coli: p34 = ',p34
           write(testout,111) ' D0m0_coli: p14 = ',p14
           write(testout,111) ' D0m0_coli: p24 = ',p24
           write(testout,111) ' D0m0_coli: p13 = ',p13
           write(testout,111) ' D0m0_coli: m12 = ',m12
           write(testout,111) ' D0m0_coli: m22 = ',m22
           write(testout,111) ' D0m0_coli: m32 = ',m32
           write(testout,111) ' D0m0_coli: m42 = ',m42
           write(testout,111) ' D0m0_coli  = ',D0m0_coli
           write(testout,111) ' D0m0_check = ',D0m0_check
           write(testout,101) ' diff        = ',
     &         cdabs(D0m0_coli/D0m0_check-1d0)

           write(*,100)   ' D0m0_coli inconsistent 4'
           write(*,111) ' D0m0_coli: p12 = ',p12
           write(*,111) ' D0m0_coli: p23 = ',p23
           write(*,111) ' D0m0_coli: p34 = ',p34
           write(*,111) ' D0m0_coli: p14 = ',p14
           write(*,111) ' D0m0_coli: p24 = ',p24
           write(*,111) ' D0m0_coli: p13 = ',p13
           write(*,111) ' D0m0_coli: m12 = ',m12
           write(*,111) ' D0m0_coli: m22 = ',m22
           write(*,111) ' D0m0_coli: m32 = ',m32
           write(*,111) ' D0m0_coli: m42 = ',m42
           write(*,111) ' D0m0_coli  = ',D0m0_coli
           write(*,111) ' D0m0_check = ',D0m0_check
           write(*,101) ' diff        = ',
     &         cdabs(D0m0_coli/D0m0_check-1d0)
           stop
         endif
        
       endif
#endif

      end

************************************************************************
      function D02m0_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*  scalar 4-point function  for m22 = m32 = 0                          *
*  regular case based on                                               *
*        A.Denner, U.Nierste and R.Scharf, Nucl. Phys. B367 (1991) 637 *
*                                                                      *
*                     0d0                                              *
*       p12  ---------------------  p23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*              m12| 1     3 | 0d0                                      *
*                 |         |                                          *
*                 |    4    |                                          *
*       p14  ---------------------  p34                                *
*                     m42                                              *
*                                                                      *
*----------------------------------------------------------------------*
*  14.10.08 Ansgar Denner       last changed 14.10.08 Ansgar Denner    *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 l12,l13,l14,l23,l24,l34
      complex*16 m1,m4
      real*8     ir14
      real*8     ix(2)
      complex*16 r14,r41
      complex*16 a,b,c,d
      complex*16 x(2)
      complex*16 omxy3,omxy4
      complex*16 D02m0_coli
      complex*16 cspcon_coli,cspcos_coli,cln_coli 
c      complex*16 test
      integer    i
c      complex*16 cd0
c      data       cd0 /(0d0,0d0)/
 
#ifdef CHECK
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      complex*16 elimminf2_coli
      complex*16 D02m0_check
      complex*16 k12,k13,k14,k23,k24,k34
      complex*16 mm3,mm2,v
c      complex*16 r24,r12
c      real*8     ir24,ir12
      integer    j
      logical flag(0:4)
      data    flag /5*.true./
      save    flag  
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) 'D02m0_coli in ',p12,p23,p34,p14,p13,p24,m12,m22,m32,m42

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

        if(m22.ne.cd0.or.m32.ne.cd0.or.m12.eq.cd0.or.m42.eq.cd0
     &      .or.p12.eq.m12.and.p24.eq.m42
     &      .or.p13.eq.m12.and.p34.eq.m42
     &      .or.ps23.eq.cd0
     &      .or.ps12.eq.cd0.and.ms12.eq.cd0
     &      .or.ps13.eq.cd0.and.ms12.eq.cd0
     &      .or.ps24.eq.cd0.and.ms42.eq.cd0
     &      .or.ps34.eq.cd0.and.ms42.eq.cd0
     &      .or.ps14.eq.cd0.and.ms12.eq.cd0.and.ms42.eq.cd0
     &      .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &      .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &      .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0
     &      ) then
          write(*,100)'D02m0_coli called improperly:'
          write(*,111)'D02m0_coli: p12 = ',p12,ps12
          write(*,111)'D02m0_coli: p23 = ',p23,ps23
          write(*,111)'D02m0_coli: p34 = ',p34,ps34
          write(*,111)'D02m0_coli: p14 = ',p14,ps14
          write(*,111)'D02m0_coli: p13 = ',p13,ps13
          write(*,111)'D02m0_coli: p24 = ',p24,ps24
          write(*,111)'D02m0_coli: m12 = ',m12,ms12
          write(*,111)'D02m0_coli: m22 = ',m22,ms22
          write(*,111)'D02m0_coli: m32 = ',m32,ms32
          write(*,111)'D02m0_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif

c     if (m12.eq.cd0) then
c     D02m0_coli = cD0m000(p24,p23,p13,p14,p34,p12,m42,cd0,cd0,cd0)
c     else if (m42.eq.cd0) then
c     D02m0_coli = cD0m000(p12,p23,p34,p14,p13,p24,m12,cd0,cd0,cd0)
c     else


      m1 = 0d0
      m4 = 0d0
c     write(*,*) 'D02m0_coli in',q12,q23,q34,q14,q13,q24,m12,m4
c     (permutation 1 <-> 2  for  m -> k )
      l12 = (m12    -p12)
      l13 = (m12    -p13)
      l14 = (m12+m42-p14)
      l23 = (       -p23)
      l24 = (    m42-p24)
      l34 = (    m42-p34)
      
      r14 = l14/(2d0*m42)*(1d0+sqrt(1d0-4d0*m12*m42/l14**2))
      r41 = r14*m42/m12
      a   =  l34*l24-l23*m42
      b   =  l13*l24+l12*l34-l14*l23
      c   =  l13*l12-l23*m12
      d   =  l23
      x(1) = (-b+sqrt(b*b-4d0*a*c))/(2d0*a)
      x(2) = (-b-sqrt(b*b-4d0*a*c))/(2d0*a)
      if(abs(x(1)).gt.abs(x(2))) then
        x(2) = c/(a*x(1))
      else
        x(1) = c/(a*x(2))
      end if
      if(dreal(l14).lt.0d0) then
        ir14 = sign(1d1,1d0-abs(r14*r14*m42/m12))
      else
        ir14 = 0d0
      end if
      ix(1) = -sign(1d0,dreal(d))
      ix(2) = +sign(1d0,dreal(d))
      ix(1) = 1d0
      ix(2) = -1d0
      
c>        write(*,*) 'D02m0_coli p  ',p12,p23,p34
c>        write(*,*) 'D02m0_coli p  ',p14,p13,p24
c>        write(*,*) 'D02m0_coli m  ',m12,m42
c>        write(*,*) 'D02m0_coli l  ',l12,l13,l14
c>        write(*,*) 'D02m0_coli l  ',l23,l24,l34
c>        write(*,*) 'D02m0_coli r  ',r14,r41
c>        write(*,*) 'D02m0_coli ir ',ir14
c>        write(*,*) 'D02m0_coli a  ',a,b,c,d
c>        write(*,*) 'D02m0_coli x ',x(1),x(2)
c>        write(*,*) 'D02m0_coli ix',ix(1),ix(2)

      D02m0_coli = dcmplx(0d0)
      
      if(l13.ne.0d0.and.l12.ne.0d0) then

#ifdef ALLCHECK
c ppenj
        if (flag(1)) then 
          write(*,100)' D02m0_coli: ',
     &        '    case with l13=/=0 and l12=/=0',
     &        '    not yet tested in physical process '
          write(*,111)' D02m0_coli: p12 = ',p12
          write(*,111)' D02m0_coli: p23 = ',p23
          write(*,111)' D02m0_coli: p34 = ',p34
          write(*,111)' D02m0_coli: p14 = ',p14
          write(*,111)' D02m0_coli: p24 = ',p24
          write(*,111)' D02m0_coli: p13 = ',p13
          write(*,111)' D02m0_coli: m12 = ',m12
          write(*,111)' D02m0_coli: m22 = ',m22
          write(*,111)' D02m0_coli: m32 = ',m32
          write(*,111)' D02m0_coli: m42 = ',m42
          flag(1)=.false.
        end if
#endif

        do i=1,2
          omxy3 = 1d0+x(i)*l34/l13
          omxy4 = 1d0+x(i)*l24/l12
          if (abs(omxy3).lt.abs(omxy4).and.abs(omxy3).lt.1d0) then
            omxy3 = (m12+l14*x(i)+m42*x(i)*x(i))*l23
     &          /(l12*omxy4*l13)
          else if (abs(omxy4).lt..1d0) then
            omxy4 = (m12+l14*x(i)+m42*x(i)*x(i))*l23
     &          /(l12*omxy3*l13)
          end if
          D02m0_coli = D02m0_coli + (2*i-3) * (
     &        cspcos_coli(-x(i),r41,-ix(i),ir14) 
     &        + cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
c     &           - cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13))
     &        - cspcon_coli(-x(i),l34/l13,-x(i)*l34/l13,
     &        omxy3,-ix(i),dreal(l34-l13))
c     &           - cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12))
     &        - cspcon_coli(-x(i),l24/l12,-x(i)*l24/l12,
     &        omxy4,-ix(i),dreal(l24-l12))
     &        + cln_coli(-x(i),-ix(i))*
     &        (cln_coli(l12/m12,-1d0)
     &        +cln_coli(l13/l23,dreal(l13-l23))) )

c>            write(*,*) 'D02m0_coli contr ',
c>     &           (2*i-3) *cspcos_coli(-x(i),r41,-ix(i),ir14) 
c>            write(*,*) 'D02m0_coli contr ', 
c>     &           (2*i-3) *cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
c>            write(*,*) 'D02m0_coli contr ', 
c>     &          -(2*i-3) *cspcos_coli(-x(i),l34/l13,-ix(i),
c>     &                   dreal(l34-l13)), 
c>     &           - (2*i-3) * cspcon_coli(-x(i),l34/l13,-x(i)*l34/l13,
c>     &             omxy3,-ix(i),dreal(l34-l13),dreal(l13))
c>            write(*,*) 'D02m0_coli contr ', 
c>     &         -(2*i-3) * cspcos_coli(-x(i),l24/l12,-ix(i),
c>     &                   dreal(l24-l12)),
c>     &           - (2*i-3) * cspcon_coli(-x(i),l24/l12,-x(i)*l24/l12,
c>     &             omxy4,-ix(i),dreal(l24-l12),dreal(l12))
c>            write(*,*) 'D02m0_coli contr ',(2*i-3)*cln_coli(-x(i),-ix(i))*
c>     &           (cln_coli(l12/m12,-1d0)
c>     &           +cln_coli(l13/l23,dreal(l13-l23)))

        end do
      else if(l13.ne.0d0.and.l12.eq.0d0) then

#ifdef ALLCHECK
c ppttbb dim
        if (flag(2)) then 
          write(*,100)' D02m0_coli: ',
     &        '    case with l13=/=0 and l12=0',
     &        '    not yet tested at all '
          write(*,111)' D02m0_coli: p12 = ',p12
          write(*,111)' D02m0_coli: p23 = ',p23
          write(*,111)' D02m0_coli: p34 = ',p34
          write(*,111)' D02m0_coli: p14 = ',p14
          write(*,111)' D02m0_coli: p24 = ',p24
          write(*,111)' D02m0_coli: p13 = ',p13
          write(*,111)' D02m0_coli: m12 = ',m12
          write(*,111)' D02m0_coli: m22 = ',m22
          write(*,111)' D02m0_coli: m32 = ',m32
          write(*,111)' D02m0_coli: m42 = ',m42
          flag(2)=.false.
        end if
#endif

        do i=1,2
          D02m0_coli = D02m0_coli + (2*i-3) * (
     &        cspcos_coli(-x(i),r41,-ix(i),ir14) 
     &        + cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
     &        - cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13)) 
     &        + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
     &        + cln_coli(l24/m12,-1d0)
     &        +cln_coli(l13/l23,dreal(l13-l23))))
          
c>        write(*,*) 'D02m0_coli contr ',
c>     &         (2*i-3) *cspcos_coli(-x(i),r41,-ix(i),ir14) 
c>        write(*,*) 'D02m0_coli contr ',
c>     &       (2*i-3) *cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
c>        write(*,*) 'D02m0_coli contr ', 
c>     &       -  (2*i-3) *cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13)) 
c>        write(*,*) 'D02m0_coli contr ', 
c>     &       (2*i-3) * cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
c>     &       +cln_coli(k24,-1d0))
c>        write(*,*) 'D02m0_coli contr ', + cln_coli(-x(i),-ix(i))*
c>     &             cln_coli(l13/l23,dreal(l13-l23))

        end do
      else if(l13.eq.0d0.and.l12.ne.0d0) then

#ifdef ALLCHECK
c ppttbb dim
        if (flag(3)) then 
          write(*,100)' D02m0_coli: ',
     &        '    case with l13=0 and l12=/=0',
     &        '    not yet tested at all '
          write(*,111)' D02m0_coli: p12 = ',p12
          write(*,111)' D02m0_coli: p23 = ',p23
          write(*,111)' D02m0_coli: p34 = ',p34
          write(*,111)' D02m0_coli: p14 = ',p14
          write(*,111)' D02m0_coli: p24 = ',p24
          write(*,111)' D02m0_coli: p13 = ',p13
          write(*,111)' D02m0_coli: m12 = ',m12
          write(*,111)' D02m0_coli: m22 = ',m22
          write(*,111)' D02m0_coli: m32 = ',m32
          write(*,111)' D02m0_coli: m42 = ',m42
          flag(3)=.false.
        end if
#endif

        do i=1,2
          D02m0_coli = D02m0_coli + (2*i-3) * (
     &        cspcos_coli(-x(i),r41,-ix(i),ir14) 
     &        + cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
     &        - cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12))
     &        + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
     &        + cln_coli(l34/m12,-1d0)
     &        +cln_coli(l12/l23,dreal(l12-l23))))
          
c>        write(*,*) 'D02m0_coli contr ',
c>     &         (2*i-3) *cspcos_coli(-x(i),r41,-ix(i),ir14) 
c>        write(*,*) 'D02m0_coli contr ',
c>     &       (2*i-3) *cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
c>        write(*,*) 'D02m0_coli contr ', 
c>     &       -  (2*i-3) *cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12))
c>        write(*,*) 'D02m0_coli contr ', 
c>     &       (2*i-3) * cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
c>     &       +cln_coli(l34/m12,-1d0))
c>        write(*,*) 'D02m0_coli contr ', + cln_coli(-x(i),-ix(i))*
c>     &             cln_coli(l12/l23,dreal(l12-l23))

        end do
      else 

#ifdef CHECK
        if (flag(4)) then 
          write(*,100)' D02m0_coli: ',
     &        '    case with l13=0 and l12=0',
     &        '    not yet tested at all '
          write(*,111)' D02m0_coli: p12 = ',p12
          write(*,111)' D02m0_coli: p23 = ',p23
          write(*,111)' D02m0_coli: p34 = ',p34
          write(*,111)' D02m0_coli: p14 = ',p14
          write(*,111)' D02m0_coli: p24 = ',p24
          write(*,111)' D02m0_coli: p13 = ',p13
          write(*,111)' D02m0_coli: m12 = ',m12
          write(*,111)' D02m0_coli: m22 = ',m22
          write(*,111)' D02m0_coli: m32 = ',m32
          write(*,111)' D02m0_coli: m42 = ',m42
          flag(4)=.false.
        end if
#endif
        do i=1,2
          D02m0_coli = D02m0_coli + (2*i-3) * (
     &        cspcos_coli(-x(i),r41,-ix(i),ir14) 
     &        + cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
     &        + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))
     &        + cln_coli(l24/m12,-1d0)
     &        +cln_coli(l34/l23,dreal(l34-l23))))
        end do
      end if 

c     write(*,*) 'D02m0_coli sum',D02m0_coli     

      D02m0_coli = D02m0_coli/(a*(x(1)-x(2)))

c     write(*,*) 'D02m0_coli = ',D02m0_coli      

#ifdef CHECK
      if(conscheck) then  

        do j=1,3
          if(j.eq.1) then
            ix(1)=-1
            ix(2)=-1 
          elseif(j.eq.2) then
            ix(1)=-1
            ix(2)=1 
          elseif(j.eq.3) then
            ix(1)=1
            ix(2)=1 
          endif

          D02m0_check = dcmplx(0d0)
          
          if(l13.ne.0d0.and.l12.ne.0d0) then

            do i=1,2
              omxy3 = 1d0+x(i)*l34/l13
              omxy4 = 1d0+x(i)*l24/l12
              if (abs(omxy3).lt.abs(omxy4).and.abs(omxy3).lt.1d0) then
                omxy3 = (m12+l14*x(i)+m42*x(i)*x(i))*l23
     &              /(l12*omxy4*l13)
              else if (abs(omxy4).lt..1d0) then
                omxy4 = (m12+l14*x(i)+m42*x(i)*x(i))*l23
     &              /(l12*omxy3*l13)
              end if
              D02m0_check = D02m0_check + (2*i-3) * (
     &            cspcos_coli(-x(i),r41,-ix(i),ir14) 
     &            + cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
c     &           - cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13))
     &            - cspcon_coli(-x(i),l34/l13,-x(i)*l34/l13,
     &            omxy3,-ix(i),dreal(l34-l13))
c     &           - cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12))
     &            - cspcon_coli(-x(i),l24/l12,-x(i)*l24/l12,
     &            omxy4,-ix(i),dreal(l24-l12))
     &            + cln_coli(-x(i),-ix(i))*
     &            (cln_coli(l12/m12,-1d0)
     &            +cln_coli(l13/l23,dreal(l13-l23))) )

c>            write(*,*) 'D02m0_check contr ',
c>     &           (2*i-3) *cspcos_coli(-x(i),r41,-ix(i),ir14) 
c>            write(*,*) 'D02m0_check contr ', 
c>     &           (2*i-3) *cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
c>            write(*,*) 'D02m0_check contr ', 
c>     &          -(2*i-3) *cspcos_coli(-x(i),l34/l13,-ix(i),
c>     &                   dreal(l34-l13)), 
c>     &           - (2*i-3) * cspcon_coli(-x(i),l34/l13,-x(i)*l34/l13,
c>     &             omxy3,-ix(i),dreal(l34-l13),dreal(l13))
c>            write(*,*) 'D02m0_check contr ', 
c>     &         -(2*i-3) * cspcos_coli(-x(i),l24/l12,-ix(i),
c>     &                   dreal(l24-l12)),
c>     &           - (2*i-3) * cspcon_coli(-x(i),l24/l12,-x(i)*l24/l12,
c>     &             omxy4,-ix(i),dreal(l24-l12),dreal(l12))
c>            write(*,*) 'D02m0_check contr ',(2*i-3)*cln_coli(-x(i),-ix(i))*
c>     &           (cln_coli(l12/m12,-1d0)
c>     &           +cln_coli(l13/l23,dreal(l13-l23)))

            end do
          else if(l13.ne.0d0.and.l12.eq.0d0) then


            do i=1,2
              D02m0_check = D02m0_check + (2*i-3) * (
     &            cspcos_coli(-x(i),r41,-ix(i),ir14) 
     &            + cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
     &            - cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13)) 
     &            + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
     &            + cln_coli(l24/m12,-1d0)
     &            +cln_coli(l13/l23,dreal(l13-l23))))
              
c>        write(*,*) 'D02m0_check contr ',
c>     &         (2*i-3) *cspcos_coli(-x(i),r41,-ix(i),ir14) 
c>        write(*,*) 'D02m0_check contr ',
c>     &       (2*i-3) *cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
c>        write(*,*) 'D02m0_check contr ', 
c>     &       -  (2*i-3) *cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13)) 
c>        write(*,*) 'D02m0_check contr ', 
c>     &       (2*i-3) * cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
c>     &       +cln_coli(k24,-1d0))
c>        write(*,*) 'D02m0_check contr ', + cln_coli(-x(i),-ix(i))*
c>     &             cln_coli(l13/l23,dreal(l13-l23))

            end do
          else if(l13.eq.0d0.and.l12.ne.0d0) then

            do i=1,2
              D02m0_check = D02m0_check + (2*i-3) * (
     &            cspcos_coli(-x(i),r41,-ix(i),ir14) 
     &            + cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
     &            - cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12))
     &            + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
     &            + cln_coli(l34/m12,-1d0)
     &            +cln_coli(l12/l23,dreal(l12-l23))))
              
c>        write(*,*) 'D02m0_check contr ',
c>     &         (2*i-3) *cspcos_coli(-x(i),r41,-ix(i),ir14) 
c>        write(*,*) 'D02m0_check contr ',
c>     &       (2*i-3) *cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
c>        write(*,*) 'D02m0_check contr ', 
c>     &       -  (2*i-3) *cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12))
c>        write(*,*) 'D02m0_check contr ', 
c>     &       (2*i-3) * cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
c>     &       +cln_coli(l34/m12,-1d0))
c>        write(*,*) 'D02m0_check contr ', + cln_coli(-x(i),-ix(i))*
c>     &             cln_coli(l12/l23,dreal(l12-l23))

            end do
          else 
            do i=1,2
              D02m0_check = D02m0_check + (2*i-3) * (
     &            cspcos_coli(-x(i),r41,-ix(i),ir14) 
     &            + cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
     &            + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))
     &            + cln_coli(l24/m12,-1d0)
     &            +cln_coli(l34/l23,dreal(l34-l23))))
            end do
          end if 

c     write(*,*) 'D02m0_check sum',D02m0_check     

          D02m0_check = D02m0_check/(a*(x(1)-x(2)))


          if(cdabs(D02m0_coli/D02m0_check-1d0).gt.testacc) then
            write(testout,*) ' D02m0_coli inconsistent',i
            write(testout,111) ' D02m0_coli: p12 = ',p12
            write(testout,111) ' D02m0_coli: p23 = ',p23
            write(testout,111) ' D02m0_coli: p34 = ',p34
            write(testout,111) ' D02m0_coli: p14 = ',p14
            write(testout,111) ' D02m0_coli: p24 = ',p24
            write(testout,111) ' D02m0_coli: p13 = ',p13
            write(testout,111) ' D02m0_coli: m12 = ',m12
            write(testout,111) ' D02m0_coli: m22 = ',m22
            write(testout,111) ' D02m0_coli: m32 = ',m32
            write(testout,111) ' D02m0_coli: m42 = ',m42
            write(testout,111) ' D02m0_coli  = ',D02m0_coli
            write(testout,111) ' D02m0_check = ',D02m0_check
            write(testout,101) ' diff        = ',
     &          cdabs(D02m0_coli/D02m0_check-1d0)

            write(*,*) ' D02m0_coli inconsistent',i
            write(*,111) ' D02m0_coli: p12 = ',p12
            write(*,111) ' D02m0_coli: p23 = ',p23
            write(*,111) ' D02m0_coli: p34 = ',p34
            write(*,111) ' D02m0_coli: p14 = ',p14
            write(*,111) ' D02m0_coli: p24 = ',p24
            write(*,111) ' D02m0_coli: p13 = ',p13
            write(*,111) ' D02m0_coli: m12 = ',m12
            write(*,111) ' D02m0_coli: m22 = ',m22
            write(*,111) ' D02m0_coli: m32 = ',m32
            write(*,111) ' D02m0_coli: m42 = ',m42
            write(*,111) ' D02m0_coli  = ',D02m0_coli
            write(*,111) ' D02m0_check = ',D02m0_check
            write(*,101) ' diff        = ',
     &          cdabs(D02m0_coli/D02m0_check-1d0)
            stop
          endif
        enddo

        m1  = sqrt(m12)
        m4  = sqrt(m42)

c     write(*,*) 'D02m0_check in',p12,p23,p34,p14,p13,p24,m12,m4
c     (permutation 1 <-> 2  for  m -> k )
        mm3  = dreal(m1)
        mm2  = dreal(m1)
        k12 = (m12    -p12)/(m1*mm2)
        k13 = (m12    -p13)/(m1*mm3)
        k14 = (m12+m42-p14)/(m1*m4)
        k23 = (       -p23)/(mm2*mm3)
        k24 = (    m42-p24)/(mm2*m4)
        k34 = (    m42-p34)/(mm3*m4)
        
        r14 = k14/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k14**2)))
        a   =  k34*k24-k23
        b   =  k13*k24+k12*k34-k14*k23
        c   =  k13*k12-k23
        d   =  k23
        x(1) = (-b+sqrt(b*b-4d0*a*c))/(2d0*a)
        x(2) = (-b-sqrt(b*b-4d0*a*c))/(2d0*a)
        if(abs(x(1)).gt.abs(x(2))) then
          x(2) = c/(a*x(1))
        else
          x(1) = c/(a*x(2))
        end if
        if(dreal(k14).lt.-2d0) then
          ir14 = sign(1d1,1d0-abs(r14))
        else
          ir14 = 0d0
        end if
        ix(1) = -sign(1d0,dreal(d))
        ix(2) = +sign(1d0,dreal(d))
        
c>        write(*,*) 'D02m0_checkm p  ',p12,p23,p34
c>        write(*,*) 'D02m0_checkm p  ',p14,p13,p24
c>        write(*,*) 'D02m0_checkm m  ',m12,m42
c>        write(*,*) 'D02m0_checkm k  ',k12,k13,k14
c>        write(*,*) 'D02m0_checkm k  ',k23,k24,k34
c>        write(*,*) 'D02m0_checkm r  ',r14
c>        write(*,*) 'D02m0_checkm ir ',ir14
c>        write(*,*) 'D02m0_checkm a  ',a,b,c,d
c>        write(*,*) 'D02m0_checkm x ',x(1),x(2)
c>        write(*,*) 'D02m0_checkm ix',ix(1),ix(2)

        D02m0_check = dcmplx(0d0)

        v = m1/m4
        
        if(k13.ne.0d0.and.k12.ne.0d0) then
          do i=1,2
            omxy3 = 1d0+x(i)*k34/k13
            omxy4 = 1d0+x(i)*k24/k12
            if (abs(omxy3).lt.abs(omxy4).and.abs(omxy3).lt.1d0) then
              omxy3 = (1d0+k14*x(i)+x(i)*x(i))*k23
     &            /(k12*omxy4*k13)
            else if (abs(omxy4).lt..1d0) then
              omxy4 = (1d0+k14*x(i)+x(i)*x(i))*k23
     &            /(k12*omxy3*k13)
            end if
            D02m0_check = D02m0_check + (2*i-3) * (
     &          cspcos_coli(-x(i)*v,r14/v,-ix(i),ir14) 
     &          + cspcos_coli(-x(i)*v,1d0/r14/v,-ix(i),-ir14)
c     &           - cspcos_coli(-x(i)*v,k34/k13/v,-ix(i),
c     &                   dreal((k34-k13*v)*m4*mm3))
     &          - cspcon_coli(-x(i)*v,k34/k13/v,-x(i)*k34/k13,omxy3,
     &          -ix(i),dreal((k34-k13*v)*m4*mm3))
c     &           - cspcos_coli(-x(i)*v,k24/k12/v,-ix(i),
c     &                   dreal((k24-k12*v)*m4*mm2))
     &          - cspcon_coli(-x(i)*v,k24/k12/v,-x(i)*k24/k12,omxy4,
     &          -ix(i),dreal((k24-k12*v)*m4*mm2))
     &          + cln_coli(-x(i)*v,-ix(i))*
     &          (cln_coli(k12/m1*mm2,-1d0)
     &          +cln_coli(k13*m1/k23/mm2,dreal((k13*m1-k23*mm2)*mm3))) )

c>            write(*,*) 'D02m0_check contr ',
c>     &           (2*i-3) *cspcos_coli(-x(i)*v,r14/v,-ix(i),ir14) 
c>            write(*,*) 'D02m0_check contr ', 
c>     &           (2*i-3) *cspcos_coli(-x(i)*v,1d0/r14/v,-ix(i),-ir14)
c>            write(*,*) 'D02m0_check contr ', 
c>     &          -(2*i-3) *cspcos_coli(-x(i)*v,k34/k13/v,-ix(i),
c>     &                   dreal((k34-k13*v)*m4*mm3)) 
c>            write(*,*) 'D02m0_check contr ', 
c>     &         -(2*i-3) * cspcos_coli(-x(i)*v,k24/k12/v,-ix(i),
c>     &                   dreal((k24-k12*v)*m4*mm2))
c>            write(*,*) 'D02m0_check contr ',(2*i-3)*cln_coli(-x(i)*v,-ix(i))*
c>     &           (cln_coli(k12/m1*mm2,-1d0)
c>     &           +cln_coli(k13*m1/k23/mm2,dreal(k13*m1-k23*mm2)*mm3))

          end do
        else if(k13.ne.0d0.and.k12.eq.0d0) then
          do i=1,2
            D02m0_check = D02m0_check + (2*i-3) * (
     &          cspcos_coli(-x(i)*v,r14/v,-ix(i),ir14) 
     &          + cspcos_coli(-x(i)*v,1d0/r14/v,-ix(i),-ir14)
     &          - cspcos_coli(-x(i)*v,k34/k13/v,-ix(i),
     &          dreal(k34-k13*v)*m4*mm3) 
     &          + cln_coli(-x(i)*v,-ix(i))*(cln_coli(-x(i)*v,-ix(i))/2d0
     &          + cln_coli(k24*mm2*m4/m12,-1d0)+cln_coli(k13*m1/k23/mm2,
     &          dreal(k13*m1-k23*mm2)*mm3)))
            
c>        write(*,*) 'D02m0_check contr ',
c>     &         (2*i-3) *cspcos_coli(-x(i),r14,-ix(i),ir14) 
c>        write(*,*) 'D02m0_check contr ',
c>     &       (2*i-3) *cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
c>        write(*,*) 'D02m0_check contr ', 
c>     &       -  (2*i-3) *cspcos_coli(-x(i),k34/k13,-ix(i),dreal(k34-k13)) 
c>        write(*,*) 'D02m0_check contr ', 
c>     &       (2*i-3) * cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
c>     &       +cln_coli(k24,-1d0))
c>        write(*,*) 'D02m0_check contr ', + cln_coli(-x(i),-ix(i))*
c>     &             cln_coli(k13/k23,dreal(k13-k23))

          end do
        else if(k13.eq.0d0.and.k12.ne.0d0) then
          do i=1,2
            D02m0_check = D02m0_check + (2*i-3) * (
     &          cspcos_coli(-x(i)*v,r14/v,-ix(i),ir14) 
     &          + cspcos_coli(-x(i)*v,1d0/r14/v,-ix(i),-ir14)
     &          - cspcos_coli(-x(i)*v,k24/k12/v,-ix(i),dreal(k24-k12))
     &          + cln_coli(-x(i)*v,-ix(i))*(cln_coli(-x(i)*v,-ix(i))/2d0
     &          + cln_coli(k34*mm3*m4/m12,-1d0)
     &          +cln_coli(k12*m1/k23/mm3,dreal(k12*m1-k23*mm3)*mm2)))
            
c>        write(*,*) 'D02m0_check contr ',
c>     &         (2*i-3) *cspcos_coli(-x(i),r14,-ix(i),ir14) 
c>        write(*,*) 'D02m0_check contr ',
c>     &       (2*i-3) *cspcos_coli(-x(i),1d0/r14,-ix(i),-ir14)
c>        write(*,*) 'D02m0_check contr ', 
c>     &       -  (2*i-3) *cspcos_coli(-x(i),k24/k12,-ix(i),dreal(k24-k12))
c>        write(*,*) 'D02m0_check contr ', 
c>     &       (2*i-3) * cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
c>     &       +cln_coli(k34,-1d0))
c>        write(*,*) 'D02m0_check contr ', + cln_coli(-x(i),-ix(i))*
c>     &             cln_coli(k12/k23,dreal(k12-k23))

          end do
        else 
          do i=1,2
            D02m0_check = D02m0_check + (2*i-3) * (
     &          cspcos_coli(-x(i)*v,r14/v,-ix(i),ir14) 
     &          + cspcos_coli(-x(i)*v,1d0/r14/v,-ix(i),-ir14)
     &          + cln_coli(-x(i)*v,-ix(i))*(cln_coli(-x(i)*v,-ix(i))
     &          + cln_coli(k24*m4*mm2/m12,-1d0)
     &          +cln_coli(k34*m4/k23/mm2,dreal(k34*m4-k23*mm2)*mm3)))
          end do
        end if 
c     write(*,*) 'D02m0_check sum',D02m0_check     

        D02m0_check = D02m0_check/(m1*mm2*mm3*m4*a*(x(1)-x(2)))

c     write(*,*) 'D02m0_check = ',D02m0_check        

        if(cdabs(D02m0_coli/D02m0_check-1d0).gt.testacc) then
          write(testout,100) ' D02m0_coli inconsistent 4'
          write(testout,111) ' D02m0_coli: p12 = ',p12
          write(testout,111) ' D02m0_coli: p23 = ',p23
          write(testout,111) ' D02m0_coli: p34 = ',p34
          write(testout,111) ' D02m0_coli: p14 = ',p14
          write(testout,111) ' D02m0_coli: p24 = ',p24
          write(testout,111) ' D02m0_coli: p13 = ',p13
          write(testout,111) ' D02m0_coli: m12 = ',m12
          write(testout,111) ' D02m0_coli: m22 = ',m22
          write(testout,111) ' D02m0_coli: m32 = ',m32
          write(testout,111) ' D02m0_coli: m42 = ',m42
          write(testout,111) ' D02m0_coli  = ',D02m0_coli
          write(testout,111) ' D02m0_check = ',D02m0_check
          write(testout,101) ' diff        = ',
     &        cdabs(D02m0_coli/D02m0_check-1d0)

          write(*,100)   ' D02m0_coli inconsistent 4'
          write(*,111) ' D02m0_coli: p12 = ',p12
          write(*,111) ' D02m0_coli: p23 = ',p23
          write(*,111) ' D02m0_coli: p34 = ',p34
          write(*,111) ' D02m0_coli: p14 = ',p14
          write(*,111) ' D02m0_coli: p24 = ',p24
          write(*,111) ' D02m0_coli: p13 = ',p13
          write(*,111) ' D02m0_coli: m12 = ',m12
          write(*,111) ' D02m0_coli: m22 = ',m22
          write(*,111) ' D02m0_coli: m32 = ',m32
          write(*,111) ' D02m0_coli: m42 = ',m42
          write(*,111) ' D02m0_coli  = ',D02m0_coli
          write(*,111) ' D02m0_check = ',D02m0_check
          write(*,101) ' diff        = ',
     &        cdabs(D02m0_coli/D02m0_check-1d0)
          stop
        endif
        
      endif
#endif

      end

************************************************************************
      function D03m0_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)            
************************************************************************
*  scalar 4-point function  for m22 = m32 = m42 = 0                    *
*  regular case                                                        *
*                                                                      *
*                     0d0                                              *
*       p12  ---------------------  p23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*             m12 | 1     3 | 0d0                                      *
*                 |         |                                          *
*                 |    4    |                                          *
*       p14  ---------------------  p34                                *
*                     0d0                                              *
*                                                                      *
*----------------------------------------------------------------------*
*  21.08.08 Ansgar Denner       last changed 21.08.08 Ansgar Denner    *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 l12,l13,l14,l23,l24,l34
      real*8     acc
      complex*16 a,b,c,d
      complex*16 x(2)
      real*8     ix(2)
      complex*16 D03m0_coli,cspcon_coli,cln_coli
      integer    i

#ifdef CHECK
      complex*16 D03m0_check
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      logical    flag(0:4)
      data       flag /5*.true./
      integer    j
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) ' D03m0_coli in ',p12,p23,p34,p14,p13,p24,m12
      
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

        if(m22.ne.cd0.or.m32.ne.cd0.or.m42.ne.cd0.or.ms12.eq.cd0
     &      .or.ps23.eq.cd0.or.ps34.eq.cd0.or.ps24.eq.cd0
     &      .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &      .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &      .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0
     &      ) then
          write(*,100)' D03m0_coli called improperly:'
          write(*,111)' D03m0_coli: p12 = ',p12,ps12
          write(*,111)' D03m0_coli: p23 = ',p23,ps23
          write(*,111)' D03m0_coli: p34 = ',p34,ps34
          write(*,111)' D03m0_coli: p14 = ',p14,ps14
          write(*,111)' D03m0_coli: p13 = ',p13,ps13
          write(*,111)' D03m0_coli: p24 = ',p24,ps24
          write(*,111)' D03m0_coli: m12 = ',m12,ms12
          write(*,111)' D03m0_coli: m22 = ',m22,ms22
          write(*,111)' D03m0_coli: m32 = ',m32,ms32
          write(*,111)' D03m0_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif

c      if (m12.eq.0d0) then
c        D03m0_coli = cD00000(p12,p23,p34,p14,p13,p24,cd0,cd0,cd0,cd0)
c      else

      acc = 1d-7
c     eps = 1d-15
      l12 = (m12    -p12)
      l13 = (m12    -p13)
      l14 = (m12    -p14)
      l23 = (       -p23)
      l24 = (       -p24)
      l34 = (       -p34)

      a   =  l34*l24
      b   =  l13*l24+l12*l34-l14*l23
      c   =  l13*l12-l23*m12
      d   =  l23*m12
      x(1) = (-b+sqrt(b*b-4d0*a*c))/(2d0*a)
      x(2) = (-b-sqrt(b*b-4d0*a*c))/(2d0*a)
      if(abs(x(1)).gt.abs(x(2))) then
        x(2) = c/(a*x(1))
      else
        x(1) = c/(a*x(2))
      end if

c        write(*,*) ' D03m0_coli a ',a,b,c,d
c        write(*,*) ' D03m0_coli x ',x(1),x(2)
c>      write(*,*) ' D03m0_coli p ',p12,p23,p34
c>      write(*,*) ' D03m0_coli p ',p14,p13,p24
c>      write(*,*) ' D03m0_coli m ',m12,m22,m32,m42
c>      write(*,*) ' D03m0_coli k ',k12,k13,k14
c>      write(*,*) ' D03m0_coli k ',k23,k24,k34
c>      write(*,*) ' D03m0_coli a ',a,b,c,d
c>      write(*,*) ' D03m0_coli x ',x(1),x(2)
c>      write(*,*) ' D03m0_coli arg1 ',1,x(1),k14,x(1)*k14
c>      write(*,*) ' D03m0_coli arg1 ',2,x(2),k14,x(2)*k14
c>
c----> added 23.05.03 to avoid log(0) in nearly singular cases
      if (abs(m12+x(1)*l14).lt.acc*abs(m12)
     &    .or.abs(m12+x(2)*l14).lt.acc*abs(m12)
     &    .or.abs(l13+x(1)*l34).lt.acc*abs(l13).and.l13.ne.0d0
     &    .or.abs(l13+x(2)*l34).lt.acc*abs(l13).and.l13.ne.0d0
     &    .or.abs(l12+x(1)*l24).lt.acc*abs(l12).and.l12.ne.0d0
     &    .or.abs(l12+x(2)*l24).lt.acc*abs(l12).and.l12.ne.0d0) then
        l14 = (m12    -p12)
        l13 = (m12    -p13)
        l12 = (m12    -p14)
        l34 = (       -p23)
        l24 = (       -p24)
        l23 = (       -p34)
        a   =  l34*l24
        b   =  l13*l24+l12*l34-l14*l23
        c   =  l13*l12-l23*m12
        d   =  l23*m12
        x(1) = (-b+sqrt(b*b-4d0*a*c))/(2d0*a)
        x(2) = (-b-sqrt(b*b-4d0*a*c))/(2d0*a)
        if(abs(x(1)).gt.abs(x(2))) then
          x(2) = c/a/x(1)
        else
          x(1) = c/a/x(2)
        end if
c       write(*,*)  'permutation!'
c      write(*,*) ' D03m0_coli k ',k12,k13,k14
c      write(*,*) ' D03m0_coli k ',k23,k24,k34
c      write(*,*) ' D03m0_coli a ',a,b,c,d
c      write(*,*) ' D03m0_coli x ',x(1),x(2)
c      write(*,*) ' D03m0_coli arg1 ',1,x(1),k14,x(1)*k14
c      write(*,*) ' D03m0_coli arg1 ',2,x(2),k14,x(2)*k14
        if (abs(m12+x(1)*l14).lt.acc*abs(m12)
     &      .or.abs(m12+x(2)*l14).lt.acc*abs(m12)
     &      .or.abs(l13+x(1)*l34).lt.acc*abs(l13).and.l13.ne.0d0
     &      .or.abs(l13+x(2)*l34).lt.acc*abs(l13).and.l13.ne.0d0
     &      .or.abs(l12+x(1)*l24).lt.acc*abs(l12).and.l12.ne.0d0
     &      .or.abs(l12+x(2)*l24).lt.acc*abs(l12).and.l12.ne.0d0) then
          l12 = (m12    -p12)
          l14 = (m12    -p13)
          l13 = (m12    -p14)
          l24 = (       -p23)
          l23 = (       -p24)
          l34 = (       -p34)
          a   =  l34*l24
          b   =  l13*l24+l12*l34-l14*l23
          c   =  l13*l12-l23
          d   =  l23
          x(1) = (-b+sqrt(b*b-4d0*a*c))/(2d0*a)
          x(2) = (-b-sqrt(b*b-4d0*a*c))/(2d0*a)
          if(abs(x(1)).gt.abs(x(2))) then
            x(2) = c/a/x(1)
          else
            x(1) = c/a/x(2)
          end if
        end if
      end if 
      ix(1) = -sign(1d0,dreal(d))
      ix(2) = +sign(1d0,dreal(d))
      ix(1) = -1d0
      ix(2) = +1d0

      D03m0_coli = dcmplx(0d0)
      if(l13.ne.0d0.and.l12.ne.0d0) then

#ifdef CHECK
        if (flag(1)) then 
          write(*,100)' D03m0_coli: ',
     &        '    case with l13=/=0 and l12=/=0',
     &        '    not yet tested in physical process '
          write(*,111)' D03m0_coli: p12 = ',p12
          write(*,111)' D03m0_coli: p23 = ',p23
          write(*,111)' D03m0_coli: p34 = ',p34
          write(*,111)' D03m0_coli: p14 = ',p14
          write(*,111)' D03m0_coli: p24 = ',p24
          write(*,111)' D03m0_coli: p13 = ',p13
          write(*,111)' D03m0_coli: m12 = ',m12
          write(*,111)' D03m0_coli: m22 = ',m22
          write(*,111)' D03m0_coli: m32 = ',m32
          write(*,111)' D03m0_coli: m42 = ',m42
          flag(1)=.false.
        end if
#endif

        do i=1,2

        D03m0_coli = D03m0_coli + (2*i-3) * (
     &       cspcon_coli(-x(i),l14/m12,
     &        -x(i)*l14/m12,1d0+x(i)*l14/m12,-ix(i),-1d0)
     &     - cspcon_coli(-x(i),l34/l13,
     &        -x(i)*l34/l13,1d0+x(i)*l34/l13,-ix(i),dreal(l34-l13)) 
     &     - cspcon_coli(-x(i),l24/l12,
     &        -x(i)*l24/l12,1d0+x(i)*l24/l12,-ix(i),dreal(l24-l12))
     &     + cln_coli(-x(i),-ix(i))
     &        *(cln_coli(l12,-1d0)+cln_coli(l13,-1d0)
     &        -cln_coli(l23,-1d0)-log(m12) ) )
        end do
      else if(l13.eq.0d0.and.l12.ne.0d0) then

#ifdef ALLCHECK
c ppttbb dim
        if (flag(2)) then 
          write(*,100)' D03m0_coli: ',
     &        '    case with l13=0 and l12=/=0',
     &        '    not yet tested in physical process '
          write(*,111)' D03m0_coli: p12 = ',p12
          write(*,111)' D03m0_coli: p23 = ',p23
          write(*,111)' D03m0_coli: p34 = ',p34
          write(*,111)' D03m0_coli: p14 = ',p14
          write(*,111)' D03m0_coli: p24 = ',p24
          write(*,111)' D03m0_coli: p13 = ',p13
          write(*,111)' D03m0_coli: m12 = ',m12
          write(*,111)' D03m0_coli: m22 = ',m22
          write(*,111)' D03m0_coli: m32 = ',m32
          write(*,111)' D03m0_coli: m42 = ',m42
          flag(2)=.false.
        end if
#endif

        do i=1,2
        D03m0_coli = D03m0_coli + (2*i-3) * (
     &       cspcon_coli(-x(i),l14/m12,
     &        -x(i)*l14/m12,1d0+x(i)*l14/m12,-ix(i),-1d0)
     &     - cspcon_coli(-x(i),l24/l12,
     &        -x(i)*l24/l12,1d0+x(i)*l24/l12,-ix(i),dreal(l24-l12))
     &     + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
     &        +cln_coli(l12,-1d0)+cln_coli(l34,-1d0)
     &        -cln_coli(l23,-1d0)-log(m12))  )
       end do

      else if(l13.ne.0d0.and.l12.eq.0d0) then

#ifdef ALLCHECK
c ppttbb dim
        if (flag(3)) then 
          write(*,100)' D03m0_coli: ',
     &        '    case with l13=/=0 and l12=0',
     &        '    not yet tested in physical process '
          write(*,111)' D03m0_coli: p12 = ',p12
          write(*,111)' D03m0_coli: p23 = ',p23
          write(*,111)' D03m0_coli: p34 = ',p34
          write(*,111)' D03m0_coli: p14 = ',p14
          write(*,111)' D03m0_coli: p24 = ',p24
          write(*,111)' D03m0_coli: p13 = ',p13
          write(*,111)' D03m0_coli: m12 = ',m12
          write(*,111)' D03m0_coli: m22 = ',m22
          write(*,111)' D03m0_coli: m32 = ',m32
          write(*,111)' D03m0_coli: m42 = ',m42
          flag(3)=.false.
        end if
#endif

        do i=1,2
        D03m0_coli = D03m0_coli + (2*i-3) * (
     &       cspcon_coli(-x(i),l14/m12,
     &        -x(i)*l14/m12,1d0+x(i)*l14/m12,-ix(i),-1d0)
     &     - cspcon_coli(-x(i),l34/l13,
     &        -x(i)*l34/l13,1d0+x(i)*l34/l13,-ix(i),dreal(l34-l13)) 
     &     + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
     &         +cln_coli(l24,-1d0)+cln_coli(l13,-1d0)
     &        -cln_coli(l23,-1d0)-log(m12) )  )
        end do
      else if(l13.eq.0d0.and.l12.eq.0d0) then

#ifdef ALLCHECK
c ppttbb dim
        if (flag(4)) then 
          write(*,100)' D03m0_coli: ',
     &        '    case with l13=0 and l12=0',
     &        '    not yet tested in physical process '
          write(*,111)' D03m0_coli: p12 = ',p12
          write(*,111)' D03m0_coli: p23 = ',p23
          write(*,111)' D03m0_coli: p34 = ',p34
          write(*,111)' D03m0_coli: p14 = ',p14
          write(*,111)' D03m0_coli: p24 = ',p24
          write(*,111)' D03m0_coli: p13 = ',p13
          write(*,111)' D03m0_coli: m12 = ',m12
          write(*,111)' D03m0_coli: m22 = ',m22
          write(*,111)' D03m0_coli: m32 = ',m32
          write(*,111)' D03m0_coli: m42 = ',m42
          flag(4)=.false.
        end if
#endif

        do i=1,2
        D03m0_coli = D03m0_coli + (2*i-3) * (
     &         cspcon_coli(-x(i),l14/m12,
     &        -x(i)*l14/m12,1d0+x(i)*l14/m12,-ix(i),-1d0)
     &         + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))
     &         +cln_coli(l24,-1d0)+cln_coli(l34,-1d0)
     &        -cln_coli(l23,-1d0)-log(m12)  )  )
        end do
      end if

      D03m0_coli = D03m0_coli/(a*(x(1)-x(2)))



#ifdef CHECK
      if(conscheck) then
        do j=1,4
          if(j.eq.1) then
            ix(1)=-1
            ix(2)=-1 
          elseif(j.eq.2) then
            ix(1)=-1
            ix(2)=1 
          elseif(j.eq.3) then
            ix(1)=1
            ix(2)=1 
          endif
         
          D03m0_check = dcmplx(0d0)
          if(l13.ne.0d0.and.l12.ne.0d0) then
            do i=1,2

              D03m0_check = D03m0_check + (2*i-3) * (
     &            cspcon_coli(-x(i),l14/m12,
     &            -x(i)*l14/m12,1d0+x(i)*l14/m12,-ix(i),-1d0)
     &            - cspcon_coli(-x(i),l34/l13,
     &            -x(i)*l34/l13,1d0+x(i)*l34/l13,-ix(i),dreal(l34-l13)) 
     &            - cspcon_coli(-x(i),l24/l12,
     &            -x(i)*l24/l12,1d0+x(i)*l24/l12,-ix(i),dreal(l24-l12))
     &            + cln_coli(-x(i),-ix(i))
     &            *(cln_coli(l12,-1d0)+cln_coli(l13,-1d0)
     &            -cln_coli(l23,-1d0)-log(m12) ) )
            end do
          else if(l13.eq.0d0.and.l12.ne.0d0) then
            do i=1,2
              D03m0_check = D03m0_check + (2*i-3) * (
     &            cspcon_coli(-x(i),l14/m12,
     &            -x(i)*l14/m12,1d0+x(i)*l14/m12,-ix(i),-1d0)
     &            - cspcon_coli(-x(i),l24/l12,
     &            -x(i)*l24/l12,1d0+x(i)*l24/l12,-ix(i),dreal(l24-l12))
     &            + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
     &            +cln_coli(l12,-1d0)+cln_coli(l34,-1d0)
     &            -cln_coli(l23,-1d0)-log(m12) ) )
            end do
          else if(l13.ne.0d0.and.l12.eq.0d0) then
            do i=1,2
              D03m0_check = D03m0_check + (2*i-3) * (
     &            cspcon_coli(-x(i),l14/m12,
     &            -x(i)*l14/m12,1d0+x(i)*l14/m12,-ix(i),-1d0)
     &            - cspcon_coli(-x(i),l34/l13,
     &            -x(i)*l34/l13,1d0+x(i)*l34/l13,-ix(i),dreal(l34-l13)) 
     &            + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
     &            +cln_coli(l24,-1d0)+cln_coli(l13,-1d0)
     &            -cln_coli(l23,-1d0)-log(m12) ) )
            end do
          else if(l13.eq.0d0.and.l12.eq.0d0) then
            do i=1,2
              D03m0_check = D03m0_check + (2*i-3) * (
     &            cspcon_coli(-x(i),l14/m12,
     &            -x(i)*l14/m12,1d0+x(i)*l14/m12,-ix(i),-1d0)
     &            + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))
     &            +cln_coli(l24,-1d0)+cln_coli(l34,-1d0)
     &            -cln_coli(l23,-1d0)-log(m12)  )  )
            end do
          end if
          
          D03m0_check = D03m0_check/(a*(x(1)-x(2)))
 
          if(cdabs(D03m0_coli/D03m0_check-1d0).gt.testacc) then
            write(testout,100) ' D03m0_coli inconsistent'
            write(testout,111) ' D03m0_coli: p12 = ',p12
            write(testout,111) ' D03m0_coli: p23 = ',p23
            write(testout,111) ' D03m0_coli: p34 = ',p34
            write(testout,111) ' D03m0_coli: p14 = ',p14
            write(testout,111) ' D03m0_coli: p24 = ',p24
            write(testout,111) ' D03m0_coli: p13 = ',p13
            write(testout,111) ' D03m0_coli: m12 = ',m12
            write(testout,111) ' D03m0_coli: m22 = ',m22
            write(testout,111) ' D03m0_coli: m32 = ',m32
            write(testout,111) ' D03m0_coli: m42 = ',m42
            write(testout,111) ' D03m0_coli  = ',D03m0_coli
            write(testout,111) ' D03m0_check = ',D03m0_check
            write(testout,101) ' diff        = ',
     &          cdabs(D03m0_coli/D03m0_check-1d0)

            write(*,*)   ' D03m0_coli inconsistent'
            write(*,111) ' D03m0_coli: p12 = ',p12
            write(*,111) ' D03m0_coli: p23 = ',p23
            write(*,111) ' D03m0_coli: p34 = ',p34
            write(*,111) ' D03m0_coli: p14 = ',p14
            write(*,111) ' D03m0_coli: p24 = ',p24
            write(*,111) ' D03m0_coli: p13 = ',p13
            write(*,111) ' D03m0_coli: m12 = ',m12
            write(*,111) ' D03m0_coli: m22 = ',m22
            write(*,111) ' D03m0_coli: m32 = ',m32
            write(*,111) ' D03m0_coli: m42 = ',m42
            write(*,111) ' D03m0_coli  = ',D03m0_coli
            write(*,111) ' D03m0_check = ',D03m0_check
            write(*,101) ' diff        = ',
     &          cdabs(D03m0_coli/D03m0_check-1d0)
            stop
          endif
        enddo
      endif
#endif

      end


************************************************************************
      function D04m0_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*  scalar 4-point function  for m12 = m22 = m32 = m42 = 0              *
*  regular case                                                        *
*                                                                      *
*                     0d0                                              *
*       m12  ---------------------  p23                                *
*                 |    2    |                                          *
*                 |         |                                          *
*             0d0 | 1     3 | 0d0                                      *
*                 |         |                                          *
*                 |    4    |                                          *
*       p14  ---------------------  p34                                *
*                     0d0                                              *
*                                                                      *
*----------------------------------------------------------------------*
*  20.01.95 Ansgar Denner       last changed 07.05.08 Ansgar Denner    *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 l12,l13,l14,l23,l24,l34
      complex*16 m2
      real*8     eps
      complex*16 a,b,c
      complex*16 x(2)
      real*8     ix(2)
      complex*16 D04m0_coli,cspcos_coli,cln_coli
      integer    i
      common /infeps/ eps
 
#ifdef CHECK
      complex*16 D04m0_check
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      integer    j
      logical    flag(0:4)
      data       flag /5*.true./

#endif

#include "params_coli.h"
#include "checkparams_coli.h"

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

        if(m12.ne.0d0.and.m22.ne.0d0.or.m32.ne.0d0.or.m42.ne.0d0
     &      .or.ps12.eq.cd0.or.ps23.eq.cd0.or.ps34.eq.cd0
     &      .or.ps14.eq.cd0.or.ps13.eq.cd0.or.ps24.eq.cd0
     &      .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &      .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &      .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0
     &      ) then
          write(*,*) ' D04m0_coli called improperly:'
          write(*,111)' D04m0_coli: p12 = ',p12,ps12
          write(*,111)' D04m0_coli: p23 = ',p23,ps23
          write(*,111)' D04m0_coli: p34 = ',p34,ps34
          write(*,111)' D04m0_coli: p14 = ',p14,ps14
          write(*,111)' D04m0_coli: p13 = ',p13,ps13
          write(*,111)' D04m0_coli: p24 = ',p24,ps24
          write(*,111)' D04m0_coli: m12 = ',m12,ms12
          write(*,111)' D04m0_coli: m22 = ',m22,ms22
          write(*,111)' D04m0_coli: m32 = ',m32,ms32
          write(*,111)' D04m0_coli: m42 = ',m42,ms42
          stop
        end if
      end if

      if (flag(0)) then 
        write(*,100)' D04m0_coli: ',
     &      '    not yet tested at all '
        write(*,111)' D04m0_coli: p12 = ',p12
        write(*,111)' D04m0_coli: p23 = ',p23
        write(*,111)' D04m0_coli: p34 = ',p34
        write(*,111)' D04m0_coli: p14 = ',p14
        write(*,111)' D04m0_coli: p24 = ',p24
        write(*,111)' D04m0_coli: p13 = ',p13
        write(*,111)' D04m0_coli: m12 = ',m12
        write(*,111)' D04m0_coli: m22 = ',m22
        write(*,111)' D04m0_coli: m32 = ',m32
        write(*,111)' D04m0_coli: m42 = ',m42
        flag(0)=.false.
      end if
#endif

      l12 = (       -p12)
      l13 = (       -p13)
      l14 = (       -p14)
      l23 = (       -p23)
      l24 = (       -p24)
      l34 = (       -p34)
      a   =  l34*l24
      b   =  l13*l24+l12*l34-l14*l23
      c   =  l13*l12
      x(1) = (-b+sqrt(b*b-4d0*a*c))/(2d0*a)
      x(2) = (-b-sqrt(b*b-4d0*a*c))/(2d0*a)
      if(abs(x(1)).gt.abs(x(2))) then
        x(2) = c/(a*x(1))
      else
        x(1) = c/(a*x(2))
      end if
      ix(1) = -1d0
      ix(2) = +1d0
 
      D04m0_coli = dcmplx(0d0)

      if(l13.ne.0d0.and.l12.ne.0d0) then

#ifdef CHECK
        if (flag(1)) then 
          write(*,100)' D04m0_coli: ',
     &        '    case l13=/=0 and l12=/=0 ',
     &        '    not yet tested at all '
          write(*,111)' D04m0_coli: p12 = ',p12
          write(*,111)' D04m0_coli: p23 = ',p23
          write(*,111)' D04m0_coli: p34 = ',p34
          write(*,111)' D04m0_coli: p14 = ',p14
          write(*,111)' D04m0_coli: p24 = ',p24
          write(*,111)' D04m0_coli: p13 = ',p13
          write(*,111)' D04m0_coli: m12 = ',m12
          write(*,111)' D04m0_coli: m22 = ',m22
          write(*,111)' D04m0_coli: m32 = ',m32
          write(*,111)' D04m0_coli: m42 = ',m42
          flag(1)=.false.
        end if
#endif

        do i=1,2
          D04m0_coli = D04m0_coli + (2*i-3) * (
     &        - cln_coli(-x(i),-ix(i))**2/2d0
     &        - cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13))
     &        - cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12))
     &        + cln_coli(-x(i),-ix(i))
     &        *(cln_coli(l12,-1d0)+cln_coli(l13,-1d0)
     &        -cln_coli(l14,-1d0)-cln_coli(l23,-1d0))  )

c>          write(*,*) 'D04m0_coli',
c>     &        - cln_coli(-x(i),-ix(i))**2/2d0,
c>     &        - cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13)),
c>     &        - cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12)),
c>     &        + cln_coli(-x(i),-ix(i))*(cln_coli(l12,-1d0)+cln_coli(l13,-1d0)
c>     &        -cln_coli(l14,-1d0)-cln_coli(l23,-1d0))

        enddo
      else if(l13.eq.0d0.and.l12.ne.0d0) then

#ifdef CHECK
        if (flag(2)) then 
          write(*,100)' D04m0_coli: ',
     &        '    case l13=/=0 and l12=/=0 ',
     &        '    not yet tested at all '
          write(*,111)' D04m0_coli: p12 = ',p12
          write(*,111)' D04m0_coli: p23 = ',p23
          write(*,111)' D04m0_coli: p34 = ',p34
          write(*,111)' D04m0_coli: p14 = ',p14
          write(*,111)' D04m0_coli: p24 = ',p24
          write(*,111)' D04m0_coli: p13 = ',p13
          write(*,111)' D04m0_coli: m12 = ',m12
          write(*,111)' D04m0_coli: m22 = ',m22
          write(*,111)' D04m0_coli: m32 = ',m32
          write(*,111)' D04m0_coli: m42 = ',m42
          flag(2)=.false.
        end if
#endif

        do i=1,2
          D04m0_coli = D04m0_coli + (2*i-3) * (
     &        - cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12))
     &        + cln_coli(-x(i),-ix(i))
     &        *(cln_coli(l12,-1d0)+cln_coli(l34,-1d0)
     &        -cln_coli(l23,-1d0)-cln_coli(l14,-1d0))  )
        end do
      else if(l13.ne.0d0.and.l12.eq.0d0) then

#ifdef CHECK
        if (flag(3)) then 
          write(*,100)' D04m0_coli: ',
     &        '    case l13=/=0 and l12=0 ',
     &        '    not yet tested at all '
          write(*,111)' D04m0_coli: p12 = ',p12
          write(*,111)' D04m0_coli: p23 = ',p23
          write(*,111)' D04m0_coli: p34 = ',p34
          write(*,111)' D04m0_coli: p14 = ',p14
          write(*,111)' D04m0_coli: p24 = ',p24
          write(*,111)' D04m0_coli: p13 = ',p13
          write(*,111)' D04m0_coli: m12 = ',m12
          write(*,111)' D04m0_coli: m22 = ',m22
          write(*,111)' D04m0_coli: m32 = ',m32
          write(*,111)' D04m0_coli: m42 = ',m42
          flag(3)=.false.
        end if
#endif

        do i=1,2
          D04m0_coli = D04m0_coli + (2*i-3) * (
     &        - cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13))
     &        + cln_coli(-x(i),-ix(i))
     &        *(cln_coli(l24,-1d0)+cln_coli(l13,-1d0)
     &        -cln_coli(l23,-1d0)-cln_coli(l14,-1d0))  )
        end do
      else if(l13.eq.0d0.and.l12.eq.0d0) then

#ifdef CHECK
        if (flag(4)) then 
          write(*,100)' D04m0_coli: ',
     &        '    case l13=0 and l12=0 ',
     &        '    not yet tested at all '
          write(*,111)' D04m0_coli: p12 = ',p12
          write(*,111)' D04m0_coli: p23 = ',p23
          write(*,111)' D04m0_coli: p34 = ',p34
          write(*,111)' D04m0_coli: p14 = ',p14
          write(*,111)' D04m0_coli: p24 = ',p24
          write(*,111)' D04m0_coli: p13 = ',p13
          write(*,111)' D04m0_coli: m12 = ',m12
          write(*,111)' D04m0_coli: m22 = ',m22
          write(*,111)' D04m0_coli: m32 = ',m32
          write(*,111)' D04m0_coli: m42 = ',m42
          flag(4)=.false.
        end if
#endif

        do i=1,2
          D04m0_coli = D04m0_coli + (2*i-3) * (
     &        + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
     &        +cln_coli(l24,-1d0)+cln_coli(l34,-1d0)
     &        -cln_coli(l23,-1d0)-cln_coli(l14,-1d0))  )
        end do
      end if

c      write(*,*) ' D04m0_coli',D04m0_coli

      D04m0_coli = D04m0_coli/(m2*m2*a*(x(1)-x(2)))

c      write(*,*) ' D04m0_coli',D04m0_coli


#ifdef CHECK
      if(conscheck) then
        do j=1,4
          if(j.eq.1) then
            ix(1)=-1
            ix(2)=-1 
          elseif(j.eq.2) then
            ix(1)=-1
            ix(2)=1 
          elseif(j.eq.3) then
            ix(1)=1
            ix(2)=1 
          endif
         
          D04m0_check = dcmplx(0d0)

          if(l13.ne.0d0.and.l12.ne.0d0) then
            do i=1,2
              D04m0_check = D04m0_check + (2*i-3) * (
     &            - cln_coli(-x(i),-ix(i))**2/2d0
     &            - cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13))
     &            - cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12))
     &            + cln_coli(-x(i),-ix(i))
     &            *(cln_coli(l12,-1d0)+cln_coli(l13,-1d0)
     &            -cln_coli(l14,-1d0)-cln_coli(l23,-1d0))  )

c>          write(*,*) 'D04m0_check',
c>     &        - cln_coli(-x(i),-ix(i))**2/2d0,
c>     &        - cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13)),
c>     &        - cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12)),
c>     &        + cln_coli(-x(i),-ix(i))*(cln_coli(l12,-1d0)+cln_coli(l13,-1d0)
c>     &        -cln_coli(l14,-1d0)-cln_coli(l23,-1d0))

            enddo
          else if(l13.eq.0d0.and.l12.ne.0d0) then
            do i=1,2
              D04m0_check = D04m0_check + (2*i-3) * (
     &            - cspcos_coli(-x(i),l24/l12,-ix(i),dreal(l24-l12))
     &            + cln_coli(-x(i),-ix(i))
     &            *(cln_coli(l12,-1d0)+cln_coli(l34,-1d0)
     &            -cln_coli(l23,-1d0)-cln_coli(l14,-1d0))  )
            end do
          else if(l13.ne.0d0.and.l12.eq.0d0) then
            do i=1,2
              D04m0_check = D04m0_check + (2*i-3) * (
     &            - cspcos_coli(-x(i),l34/l13,-ix(i),dreal(l34-l13))
     &            + cln_coli(-x(i),-ix(i))
     &            *(cln_coli(l24,-1d0)+cln_coli(l13,-1d0)
     &            -cln_coli(l23,-1d0)-cln_coli(l14,-1d0))  )
            end do
          else if(l13.eq.0d0.and.l12.eq.0d0) then
            do i=1,2
              D04m0_check = D04m0_check + (2*i-3) * (
     &            + cln_coli(-x(i),-ix(i))*(cln_coli(-x(i),-ix(i))/2d0
     &            +cln_coli(l24,-1d0)+cln_coli(l34,-1d0)
     &            -cln_coli(l23,-1d0)-cln_coli(l14,-1d0))  )
            end do
          end if


          D04m0_check = D04m0_check/(m2*m2*a*(x(1)-x(2)))

c      write(*,*) ' D04m0_check =   ',D04m0_check,j    
 
          if(cdabs(D04m0_coli/D04m0_check-1d0).gt.testacc) then
            write(testout,100)   ' D04m0_coli inconsistent'
            write(testout,111) ' D04m0_coli: p12 = ',p12
            write(testout,111) ' D04m0_coli: p23 = ',p23
            write(testout,111) ' D04m0_coli: p34 = ',p34
            write(testout,111) ' D04m0_coli: p14 = ',p14
            write(testout,111) ' D04m0_coli: p24 = ',p24
            write(testout,111) ' D04m0_coli: p13 = ',p13
            write(testout,111) ' D04m0_coli: m12 = ',m12
            write(testout,111) ' D04m0_coli: m22 = ',m22
            write(testout,111) ' D04m0_coli: m32 = ',m32
            write(testout,111) ' D04m0_coli: m42 = ',m42
            write(testout,111) ' D04m0_coli  = ',D04m0_coli
            write(testout,111) ' D04m0_check = ',D04m0_check
            write(testout,101) ' diff        = ',
     &          cdabs(D04m0_coli/D04m0_check-1d0)

            write(*,*)   ' D04m0_coli inconsistent'
            write(*,111) ' D04m0_coli: p12 = ',p12
            write(*,111) ' D04m0_coli: p23 = ',p23
            write(*,111) ' D04m0_coli: p34 = ',p34
            write(*,111) ' D04m0_coli: p14 = ',p14
            write(*,111) ' D04m0_coli: p24 = ',p24
            write(*,111) ' D04m0_coli: p13 = ',p13
            write(*,111) ' D04m0_coli: m12 = ',m12
            write(*,111) ' D04m0_coli: m22 = ',m22
            write(*,111) ' D04m0_coli: m32 = ',m32
            write(*,111) ' D04m0_coli: m42 = ',m42
            write(*,111) ' D04m0_coli  = ',D04m0_coli
            write(*,111) ' D04m0_check = ',D04m0_check
            write(*,101) ' diff        = ',
     &          cdabs(D04m0_coli/D04m0_check-1d0)
            stop
          endif
        enddo
      endif
#endif


      end


************************************************************************
      function D0ms0ir1m0_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     4-point function                                                 *
*     with 0 mass singularities and up to 1 soft singularities         *
*                                                                      *
*     assumes                                                          *
*     m12 small, p12=m22, p14=m42, m32=0                               *
*                                                                      *
*                        m22                                           *
*       p12=m22  ---------------------  p23                            *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | 0                                    *
*                     |         |                                      *
*                     |    4    |                                      *
*       p14=m42  ---------------------  p34                            *
*                         m42                                          *
*                                                                      *
************************************************************************
*     22.08.08 Ansgar Denner        last changed  10.10.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms0ir1m0_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
c      complex*16 D0ms0ir1m0_check
      logical flag(0:6)
      data    flag /7*.true./
      save    flag  
#endif

      real*8     ir24
      complex*16 m2,m4,k24,r24
      complex*16 y,logxs,logy
      complex*16 cln_coli,cspenc_coli

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
#ifdef CHECK
 101  format(a22,g25.17)
      if(argcheck)then
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

        if(ms12.ne.cd0.or.ms32.ne.cd0.or.ms22.eq.cd0.or.ms42.eq.cd0
     &        .or.ps12.ne.ms22.or.ps14.ne.ms42 
     &        .or.ps23.eq.ms22.and.ps34.eq.ms42
     &        .or.ps24.eq.cd0.or.ps13.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &        .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0ms0ir1m0_coli called improperly:'
          write(*,111)' D0ms0ir1m0_coli: p12 = ',p12,ps12
          write(*,111)' D0ms0ir1m0_coli: p23 = ',p23,ps23
          write(*,111)' D0ms0ir1m0_coli: p34 = ',p34,ps34
          write(*,111)' D0ms0ir1m0_coli: p14 = ',p14,ps14
          write(*,111)' D0ms0ir1m0_coli: p13 = ',p13,ps13
          write(*,111)' D0ms0ir1m0_coli: p24 = ',p24,ps24
          write(*,111)' D0ms0ir1m0_coli: m12 = ',m12,ms12
          write(*,111)' D0ms0ir1m0_coli: m22 = ',m22,ms22
          write(*,111)' D0ms0ir1m0_coli: m32 = ',m32,ms32
          write(*,111)' D0ms0ir1m0_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif


c according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 (ib)
c          and (ii)    (D0ir0)
c Denner Dittmaier (3.4/3.5)

      m2 = sqrt(m22)
      m4 = sqrt(m42)
      k24 = (m22+m42-p24)/(m2*m4)
      if (k24.ne.0d0) then
c     |r24| > 0 yields wrong results !!!!!!!    (25.06.03)
c     analytic continuation requires |r24| < 0
        r24 = 2d0/(k24*(1d0+sqrt(dcmplx(1d0-4d0/k24**2))))
      else
        r24 = dcmplx(0d0,1d0)
      end if
      if(dreal(k24).lt.-2d0) then
        ir24 = sign(1d1,1d0-abs(r24))
      else
        ir24 = 0d0
      end if

      if(p24.ne.(m2-m4)**2) then
        logxs = cln_coli(r24,ir24)
        if(p23.ne.m32.and.p34.ne.m42)then 
          y    = m2*(p34-m42)/(m4*(p23-m22))
          logy = log(m2/m4)+cln_coli(p34-m42,1d0)-cln_coli(p23-m22,1d0)
          D0ms0ir1m0_coli =  logxs*(
     &         2d0*cln_coli(1d0-r24*r24,-dreal(r24)*ir24)
     &         - logxs/2d0
     &         + cln_coli(p13/(p23-m22),dreal(p23-m22))
     &         + cln_coli(p13/(p34-m42),dreal(p34-m42))) 
     &         + pi2_6 + logy*logy/2d0
     &         + cspenc_coli(r24*r24,dreal(r24)*ir24)
     &         - cspenc_coli(r24*y,1d0) - cln_coli(1d0-r24*y,-1d0)
     &         *(logxs+logy)
     &         - cspenc_coli(r24/y,1d0) - cln_coli(1d0-r24/y,-1d0)
     &         *(logxs-logy)  

          if(m12.eq.cd0)then

#ifdef CHECK
            if(flag(1))then
              write(*,100)' D0ms0ir1m0_coli: (3.4) m32=0',
     &            '    general case',
     &            '    not yet tested at all '
              write(*,111)' D0ms0ir1m0_coli: p12 = ',p12
              write(*,111)' D0ms0ir1m0_coli: p23 = ',p23
              write(*,111)' D0ms0ir1m0_coli: p34 = ',p34
              write(*,111)' D0ms0ir1m0_coli: p14 = ',p14
              write(*,111)' D0ms0ir1m0_coli: p24 = ',p24
              write(*,111)' D0ms0ir1m0_coli: p13 = ',p13
              write(*,111)' D0ms0ir1m0_coli: m12 = ',m12
              write(*,111)' D0ms0ir1m0_coli: m22 = ',m22
              write(*,111)' D0ms0ir1m0_coli: m32 = ',m32
              write(*,111)' D0ms0ir1m0_coli: m42 = ',m42
              flag(1)=.false.
            endif
#endif

            D0ms0ir1m0_coli = D0ms0ir1m0_coli 
     &          -logxs*log(muir2/(m2*m4))
#ifdef SING
     &          -logxs*delta1ir
#endif
          else

#ifdef CHECK
            if(flag(2))then
              write(*,100)' D0ms0ir1m0_coli: (3.5) m32=0',
     &            '    general case',
     &            '    not yet tested at all '
              write(*,111)' D0ms0ir1m0_coli: p12 = ',p12
              write(*,111)' D0ms0ir1m0_coli: p23 = ',p23
              write(*,111)' D0ms0ir1m0_coli: p34 = ',p34
              write(*,111)' D0ms0ir1m0_coli: p14 = ',p14
              write(*,111)' D0ms0ir1m0_coli: p24 = ',p24
              write(*,111)' D0ms0ir1m0_coli: p13 = ',p13
              write(*,111)' D0ms0ir1m0_coli: m12 = ',m12
              write(*,111)' D0ms0ir1m0_coli: m22 = ',m22
              write(*,111)' D0ms0ir1m0_coli: m32 = ',m32
              write(*,111)' D0ms0ir1m0_coli: m42 = ',m42
              flag(2)=.false.
            endif
#endif

            D0ms0ir1m0_coli = D0ms0ir1m0_coli 
     &          -logxs*log(m12/(m2*m4))
          endif

          D0ms0ir1m0_coli = D0ms0ir1m0_coli*
     &         *r24/(m2*m4*p13*(1d0-r24*r24))

        else if (p23.eq.m22) then

          D0ms0ir1m0_coli = ( logxs*(
     &         2d0*cln_coli(1d0-r24*r24,-dreal(r24)*ir24)
     &         - logxs
     &         + 2d0*cln_coli(p13/(p34-m42),dreal(p34-m42)) )
     &         - pi2_6 
     &         + cspenc_coli(r24*r24,dreal(r24)*ir24) )
          if(m12.eq.cd0)then

#ifdef CHECK
            if(flag(3))then
              write(*,100)' D0ms0ir1m0_coli: (3.4) m32=0',
     &            '    case with p23=m22',
     &            '    not yet tested at all '
              write(*,111)' D0ms0ir1m0_coli: p12 = ',p12
              write(*,111)' D0ms0ir1m0_coli: p23 = ',p23
              write(*,111)' D0ms0ir1m0_coli: p34 = ',p34
              write(*,111)' D0ms0ir1m0_coli: p14 = ',p14
              write(*,111)' D0ms0ir1m0_coli: p24 = ',p24
              write(*,111)' D0ms0ir1m0_coli: p13 = ',p13
              write(*,111)' D0ms0ir1m0_coli: m12 = ',m12
              write(*,111)' D0ms0ir1m0_coli: m22 = ',m22
              write(*,111)' D0ms0ir1m0_coli: m32 = ',m32
              write(*,111)' D0ms0ir1m0_coli: m42 = ',m42
              flag(3)=.false.
            endif
#endif

            D0ms0ir1m0_coli = D0ms0ir1m0_coli 
     &          -logxs*log(muir2/m42)
#ifdef SING
     &          -logxs*delta1ir
#endif
          else

#ifdef CHECK
            if(flag(4))then
              write(*,100)' D0ms0ir1m0_coli: (3.5) m32=0',
     &            '    case with p23=m22',
     &            '    not yet tested at all '
              write(*,111)' D0ms0ir1m0_coli: p12 = ',p12
              write(*,111)' D0ms0ir1m0_coli: p23 = ',p23
              write(*,111)' D0ms0ir1m0_coli: p34 = ',p34
              write(*,111)' D0ms0ir1m0_coli: p14 = ',p14
              write(*,111)' D0ms0ir1m0_coli: p24 = ',p24
              write(*,111)' D0ms0ir1m0_coli: p13 = ',p13
              write(*,111)' D0ms0ir1m0_coli: m12 = ',m12
              write(*,111)' D0ms0ir1m0_coli: m22 = ',m22
              write(*,111)' D0ms0ir1m0_coli: m32 = ',m32
              write(*,111)' D0ms0ir1m0_coli: m42 = ',m42
              flag(4)=.false.
            endif
#endif

            D0ms0ir1m0_coli = D0ms0ir1m0_coli 
     &          -logxs*log(m12/m42)
          endif

          D0ms0ir1m0_coli = D0ms0ir1m0_coli*
     &         *r24/(m2*m4*p13*(1d0-r24*r24))

        else if (p34.eq.m42) then

          D0ms0ir1m0_coli = ( logxs*(
     &         2d0*cln_coli(1-r24*r24,-dreal(r24)*ir24)
     &         - logxs 
     &         + 2d0*cln_coli(p13/(p23-m22),dreal(p23-m22)) )
     &         - pi2_6 
     &         + cspenc_coli(r24*r24,dreal(r24)*ir24) )

          D0ms0ir1m0_coli = D0ms0ir1m0_coli*
     &         *r24/(m2*m4*p13*(1d0-r24*r24))

          if(m12.eq.cd0)then

#ifdef CHECK
            if(flag(5))then
              write(*,100)' D0ms0ir1m0_coli: (3.4) m32=0',
     &            '    case with p34=m42',
     &            '    not yet tested at all '
              write(*,111)' D0ms0ir1m0_coli: p12 = ',p12
              write(*,111)' D0ms0ir1m0_coli: p23 = ',p23
              write(*,111)' D0ms0ir1m0_coli: p34 = ',p34
              write(*,111)' D0ms0ir1m0_coli: p14 = ',p14
              write(*,111)' D0ms0ir1m0_coli: p24 = ',p24
              write(*,111)' D0ms0ir1m0_coli: p13 = ',p13
              write(*,111)' D0ms0ir1m0_coli: m12 = ',m12
              write(*,111)' D0ms0ir1m0_coli: m22 = ',m22
              write(*,111)' D0ms0ir1m0_coli: m32 = ',m32
              write(*,111)' D0ms0ir1m0_coli: m42 = ',m42
              flag(5)=.false.
            endif
#endif

            D0ms0ir1m0_coli = D0ms0ir1m0_coli 
     &          -logxs*log(muir2/m22)
#ifdef SING
     &          -logxs*delta1ir
#endif
          else

#ifdef CHECK
            if(flag(6))then
              write(*,100)' D0ms0ir1m0_coli: (3.5) m32=0',
     &            '    case with p34=m42',
     &            '    not yet tested at all '
              write(*,111)' D0ms0ir1m0_coli: p12 = ',p12
              write(*,111)' D0ms0ir1m0_coli: p23 = ',p23
              write(*,111)' D0ms0ir1m0_coli: p34 = ',p34
              write(*,111)' D0ms0ir1m0_coli: p14 = ',p14
              write(*,111)' D0ms0ir1m0_coli: p24 = ',p24
              write(*,111)' D0ms0ir1m0_coli: p13 = ',p13
              write(*,111)' D0ms0ir1m0_coli: m12 = ',m12
              write(*,111)' D0ms0ir1m0_coli: m22 = ',m22
              write(*,111)' D0ms0ir1m0_coli: m32 = ',m32
              write(*,111)' D0ms0ir1m0_coli: m42 = ',m42
              flag(6)=.false.
            endif
#endif

            D0ms0ir1m0_coli = D0ms0ir1m0_coli 
     &          -logxs*log(m12/m22)
          endif

        else 
          write(*,100)' D0ms0ir1m0_coli called improperly'
          write(*,111)' D0ms0ir1m0_coli: p12 = ',p12
          write(*,111)' D0ms0ir1m0_coli: p23 = ',p23
          write(*,111)' D0ms0ir1m0_coli: p34 = ',p34
          write(*,111)' D0ms0ir1m0_coli: p14 = ',p14
          write(*,111)' D0ms0ir1m0_coli: p13 = ',p13
          write(*,111)' D0ms0ir1m0_coli: p24 = ',p24
          write(*,111)' D0ms0ir1m0_coli: m12 = ',m12
          write(*,111)' D0ms0ir1m0_coli: m22 = ',m22
          write(*,111)' D0ms0ir1m0_coli: m32 = ',m32
          write(*,111)' D0ms0ir1m0_coli: m42 = ',m42
          stop
        endif
      else
        write(*,100)' D0ms0ir1m0_coli: case not implemented'
        if (p24.eq.(m2-m4)**2) write(*,100)' p24 = (m2-m4)^2 '
        write(*,111)' D0ms0ir1m0_coli: p12 = ',p12
        write(*,111)' D0ms0ir1m0_coli: p23 = ',p23
        write(*,111)' D0ms0ir1m0_coli: p34 = ',p34
        write(*,111)' D0ms0ir1m0_coli: p14 = ',p14
        write(*,111)' D0ms0ir1m0_coli: p13 = ',p13
        write(*,111)' D0ms0ir1m0_coli: p24 = ',p24
        write(*,111)' D0ms0ir1m0_coli: m12 = ',m12
        write(*,111)' D0ms0ir1m0_coli: m22 = ',m22
        write(*,111)' D0ms0ir1m0_coli: m32 = ',m32
        write(*,111)' D0ms0ir1m0_coli: m42 = ',m42
        stop
      end if

      end

************************************************************************
      function D0ms0ir1_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     4-point function                                                 *
*     with 0 mass singularities and up to 1 soft singularities         *
*                                                                      *
*     assumes                                                          *
*     m12 small, p12=m22, p14=m42, m32=/=0                             *
*                                                                      *
*                        m22                                           *
*       p12=m22  ---------------------  p23                            *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | m32                                  *
*                     |         |                                      *
*                     |    4    |                                      *
*       p14=m42  ---------------------  p34                            *
*                         m42                                          *
*                                                                      *
************************************************************************
*     22.08.08 Ansgar Denner        last changed  23.09.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms0ir1_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      complex*16 D0ms0ir1_check
      logical flag(0:2),errflag
      data    flag /3*.true./ , errflag /.true./
      save    flag,errflag
#endif

      complex*16 m1,m2,m3,m4
      complex*16 k24,k23,k34,r24,r23,r34
      real*8     ir24,ir23,ir34
      complex*16 logxs,logx2,logx3
      complex*16 cln_coli,cspenc_coli

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

c      write(*,*) 'D0ms0ir1_coli in ',
c     &    p12,p23,p34,p14,p13,p24,m12,m22,m32,m42

 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
#ifdef CHECK
 101  format(a22,g25.17)
      if(argcheck)then
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

        if(ms12.ne.cd0.or.ms32.eq.cd0.or.ms22.eq.cd0.or.ms42.eq.cd0
     &        .or.ps12.ne.ms22.or.ps14.ne.ms42 
     &        .or.ps24.eq.cd0.or.ps13.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &        .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0ms0ir1_coli called improperly:'
          write(*,111)' D0ms0ir1_coli: p12 = ',p12,ps12
          write(*,111)' D0ms0ir1_coli: p23 = ',p23,ps23
          write(*,111)' D0ms0ir1_coli: p34 = ',p34,ps34
          write(*,111)' D0ms0ir1_coli: p14 = ',p14,ps14
          write(*,111)' D0ms0ir1_coli: p13 = ',p13,ps13
          write(*,111)' D0ms0ir1_coli: p24 = ',p24,ps24
          write(*,111)' D0ms0ir1_coli: m12 = ',m12,ms12
          write(*,111)' D0ms0ir1_coli: m22 = ',m22,ms22
          write(*,111)' D0ms0ir1_coli: m32 = ',m32,ms32
          write(*,111)' D0ms0ir1_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif

c according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 (ib)
c          and (i)    (D0irr)
c Denner Dittmaier (3.4/3.5)

      m1 = sqrt(m12)
      m2 = sqrt(m22)
      m3 = sqrt(m32)
      m4 = sqrt(m42)
      k24 = (m22+m42-p24)/(m2*m4)
      if (k24.ne.0d0) then
c     |r24| > 0 yields wrong results !!!!!!!  
c     analytic continuation requires |r24| < 1
        r24 = 2d0/(k24*(1d0+sqrt(dcmplx(1d0-4d0/k24**2))))
      else
        r24 = dcmplx(0d0,1d0)
      end if
      if(dreal(k24).lt.-2d0) then
        ir24 = sign(1d1,1d0-abs(r24))
      else
        ir24 = 0d0
      end if
      
      k23 = (m22+m32-p23)/(m2*m3)
      k34 = (m32+m42-p34)/(m3*m4)
      if (k23.ne.0d0) then
        r23 = k23/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k23**2)))
      else
        r23 = dcmplx(0d0,1d0)
      end if
      if (dreal(k34).ne.0d0) then
        r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2)))
      else
        r34 = dcmplx(0d0,1d0)
      end if
      if(dreal(k23).lt.-2d0) then
        ir23= sign(1d1,1d0-abs(r23))
      else
        ir23= 0d0
      end if
      if(dreal(k34).lt.-2d0) then
        ir34= sign(1d1,1d0-abs(r34))
      else
        ir34= 0d0
      end if
      if(p24.ne.(m2-m4)**2) then

c      write(*,*) 'D0ms0ir1_coli  r ',r24,r23,r34,1/r34
c      write(*,*) 'D0ms0ir1_coli l args',1d0-r24*r24,1d0-r24*r34*r23,
c     &    1d0-r24/r34*r23,1d0-r24*r34/r23,1d0-r24/r34/r23

        logxs = cln_coli(r24,ir24)
        logx2 = cln_coli(r23,ir23)
        logx3 = cln_coli(r34,ir34)
        D0ms0ir1_coli = 
     &      2d0*logxs*cln_coli(1d0-r24*r24,-dreal(r24)*ir24)
     &      + 3d0*pi2_6
     &       + cspenc_coli(r24*r24,dreal(r24)*ir24)
     &       + (logx2)**2 + (logx3)**2
     &       - cspenc_coli(r24*r34*r23,1d0)
     &       - cln_coli(1d0-r24*r34*r23,-1d0)*(logxs+logx2+logx3)
     &       - cspenc_coli(r24/r34*r23,1d0) 
     &       - cln_coli(1d0-r24/r34*r23,-1d0)*(logxs+logx2-logx3)
     &       - cspenc_coli(r24*r34/r23,1d0)
     &       - cln_coli(1d0-r24*r34/r23,-1d0)*(logxs-logx2+logx3)
     &       - cspenc_coli(r24/(r34*r23),1d0) 
     &       - cln_coli(1d0-r24/(r34*r23),-1d0) 
     &                  * (logxs-logx2-logx3) 

        if(m12.eq.cd0)then

#ifdef ALLCHECK
c ppttbb dim
          if(flag(1))then
            write(*,100)' D0ms0ir1_coli: (3.4)',
     &          '    not yet tested in physical process '
            write(*,111)' D0ms0ir1_coli: p12 = ',p12
            write(*,111)' D0ms0ir1_coli: p23 = ',p23
            write(*,111)' D0ms0ir1_coli: p34 = ',p34
            write(*,111)' D0ms0ir1_coli: p14 = ',p14
            write(*,111)' D0ms0ir1_coli: p24 = ',p24
            write(*,111)' D0ms0ir1_coli: p13 = ',p13
            write(*,111)' D0ms0ir1_coli: m12 = ',m12
            write(*,111)' D0ms0ir1_coli: m22 = ',m22
            write(*,111)' D0ms0ir1_coli: m32 = ',m32
            write(*,111)' D0ms0ir1_coli: m42 = ',m42
            flag(1)=.false.
          endif
#endif

          D0ms0ir1_coli =   D0ms0ir1_coli +
     &        2d0*logxs*cln_coli((m32-p13)/(m3*sqrt(muir2)),-1d0) 
#ifdef SING
     &        -logxs*delta1ir
#endif
        else

#ifdef CHECK
         if(flag(2))then
            write(*,100)' D0ms0ir1_coli: (3.5)',
     &          '    not yet tested at all '
            write(*,111)' D0ms0ir1_coli: p12 = ',p12
            write(*,111)' D0ms0ir1_coli: p23 = ',p23
            write(*,111)' D0ms0ir1_coli: p34 = ',p34
            write(*,111)' D0ms0ir1_coli: p14 = ',p14
            write(*,111)' D0ms0ir1_coli: p24 = ',p24
            write(*,111)' D0ms0ir1_coli: p13 = ',p13
            write(*,111)' D0ms0ir1_coli: m12 = ',m12
            write(*,111)' D0ms0ir1_coli: m22 = ',m22
            write(*,111)' D0ms0ir1_coli: m32 = ',m32
            write(*,111)' D0ms0ir1_coli: m42 = ',m42
            flag(2)=.false.
          endif
#endif

          D0ms0ir1_coli =   D0ms0ir1_coli +
     &        2d0*logxs*cln_coli((m32-p13)/(m3*m1),-1d0) 
        endif

        D0ms0ir1_coli = D0ms0ir1_coli
     &       *r24/(m2*m4*(p13-m32)*(1-r24*r24))

c        write(*,*) ' D0ms0ir1_coli 2',D0ms0ir1_coli

      else
        write(*,100)' D0ms0ir1_coli: case not implemented'
        if (p24.eq.(m2-m4)**2) write(*,100)' p24 = (m2-m4)^2 '
        if (m32.eq.0) write(*,100)' m32 = 0'
        write(*,111)' D0ms0ir1: p12 = ',p12
        write(*,111)' D0ms0ir1: p23 = ',p23
        write(*,111)' D0ms0ir1: p34 = ',p34
        write(*,111)' D0ms0ir1: p14 = ',p14
        write(*,111)' D0ms0ir1: p24 = ',p24
        write(*,111)' D0ms0ir1: p13 = ',p13
        write(*,111)' D0ms0ir1: m12 = ',m12
        write(*,111)' D0ms0ir1: m22 = ',m22
        write(*,111)' D0ms0ir1: m32 = ',m32
        write(*,111)' D0ms0ir1: m42 = ',m42
        stop
      end if

#ifdef CHECK

      if (conscheck) then

        m2 = sqrt(m22)
        m4 = sqrt(m42)
        k24 = (m22+m42-p24)/(m2*m4)
        if (k24.ne.0d0) then
c     |r24| > 0 yields wrong results !!!!!!!    (25.06.03)
c     analytic continuation requires |r24| < 0
          r24 = 2d0/(k24*(1d0+sqrt(dcmplx(1d0-4d0/k24**2))))
        else
          r24 = dcmplx(0d0,1d0)
        end if
        if(dreal(k24).lt.-2d0) then
          ir24 = sign(1d1,1d0-abs(r24))
        else
          ir24 = 0d0
        end if
        
        m3 = sqrt(m32)
        k23 = (m22+m32-p23)/(m2*m3)
        k34 = (m32+m42-p34)/(m3*m4)
        if (k23.ne.0d0) then
          r23 = k23/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k23**2)))
        else
          r23 = dcmplx(0d0,1d0)
        end if
c     r23 = 1d0/r23
        if (k34.ne.0d0) then
          r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2)))
        else
          r34 = dcmplx(0d0,1d0)
        end if
c     r34 = 1d0/r34
        if(dreal(k23).lt.-2d0) then
          ir23= sign(1d1,1d0-abs(r23))
        else
          ir23= 0d0
        end if
        if(dreal(k34).lt.-2d0) then
          ir34= sign(1d1,1d0-abs(r34))
        else
          ir34= 0d0
        end if
        if(p24.ne.(m2-m4)**2) then
          logxs = cln_coli(r24,ir24)
          logx2 = cln_coli(r23,ir23)
          logx3 = cln_coli(r34,ir34)
          D0ms0ir1_check = 
     &        (2d0*logxs*(cln_coli(1d0-r24*r24,-dreal(r24)*ir24)
     &        + cln_coli((m32-p13)/(m3*sqrt(muir2)),-1d0) ) + 3d0*pi2_6
#ifdef SING
     &        -logxs*delta1ir
#endif
     &        + cspenc_coli(r24*r24,dreal(r24)*ir24)
     &        + (logx2)**2 + (logx3)**2
     &        - cspenc_coli(r24*r34*r23,1d0) 
     &        - cln_coli(1d0-r24*r34*r23,-1d0)
     &        *(logxs+logx2+logx3)
     &        - cspenc_coli(r24/r34*r23,1d0) 
     &        - cln_coli(1d0-r24/r34*r23,-1d0)
     &        *(logxs+logx2-logx3)
     &        - cspenc_coli(r24*r34/r23,1d0) 
     &        - cln_coli(1d0-r24*r34/r23,-1d0)
     &        *(logxs-logx2+logx3)
     &        - cspenc_coli(r24/(r34*r23),1d0) 
     &        - cln_coli(1d0-r24/(r34*r23),-1d0) * (logxs-logx2-logx3))
     &        *r24/(m2*m4*(p13-m32)*(1-r24*r24))
        else
          write(*,100)' D0ms0ir1_coli: case not implemented'
          if (p24.eq.(m2-m4)**2) write(*,100)' p24 = (m2-m4)^2 '
          if (m32.eq.0) write(*,100)' m32 = 0'
          write(*,111)' D0ms0ir1_coli: p12 = ',p12
          write(*,111)' D0ms0ir1_coli: p23 = ',p23
          write(*,111)' D0ms0ir1_coli: p34 = ',p34
          write(*,111)' D0ms0ir1_coli: p14 = ',p14
          write(*,111)' D0ms0ir1_coli: p24 = ',p24
          write(*,111)' D0ms0ir1_coli: p13 = ',p13
          write(*,111)' D0ms0ir1_coli: m12 = ',m12
          write(*,111)' D0ms0ir1_coli: m22 = ',m22
          write(*,111)' D0ms0ir1_coli: m32 = ',m32
          write(*,111)' D0ms0ir1_coli: m42 = ',m42
          stop
        end if
        
        if (abs(D0ms0ir1_coli/D0ms0ir1_check-1d0).gt.testacc) then
          if (errflag) then
            write(*,100)
            write(*,100) ' D0ms0ir1_coli inconsistent'
            write(*,111) ' D0ms0ir1_coli: p12 = ',p12
            write(*,111) ' D0ms0ir1_coli: p23 = ',p23
            write(*,111) ' D0ms0ir1_coli: p34 = ',p34
            write(*,111) ' D0ms0ir1_coli: p14 = ',p14
            write(*,111) ' D0ms0ir1_coli: p13 = ',p13
            write(*,111) ' D0ms0ir1_coli: p24 = ',p24
            write(*,111) ' D0ms0ir1_coli: m12 = ',m12
            write(*,111) ' D0ms0ir1_coli: m22 = ',m22
            write(*,111) ' D0ms0ir1_coli: m32 = ',m32
            write(*,111) ' D0ms0ir1_coli: m42 = ',m42
            write(*,111) ' D0ms0ir1_coli  = ',D0ms0ir1_coli
            write(*,111) ' D0ms0ir1_check = ',D0ms0ir1_check
            write(*,101) ' diff        = ',
     &        cdabs(D0ms0ir1_coli/D0ms0ir1_check-1d0)
            errflag = .true.
          end if
          write(testout,100) ' D0ms0ir1_coli inconsistent'
          write(testout,111) ' D0ms0ir1_coli: p12 = ',p12
          write(testout,111) ' D0ms0ir1_coli: p23 = ',p23
          write(testout,111) ' D0ms0ir1_coli: p34 = ',p34
          write(testout,111) ' D0ms0ir1_coli: p14 = ',p14
          write(testout,111) ' D0ms0ir1_coli: p13 = ',p13
          write(testout,111) ' D0ms0ir1_coli: p24 = ',p24
          write(testout,111) ' D0ms0ir1_coli: m12 = ',m12
          write(testout,111) ' D0ms0ir1_coli: m22 = ',m22
          write(testout,111) ' D0ms0ir1_coli: m32 = ',m32
          write(testout,111) ' D0ms0ir1_coli: m42 = ',m42
          write(testout,111) ' D0ms0ir1_coli  = ',D0ms0ir1_coli
          write(testout,111) ' D0ms0ir1_check = ',D0ms0ir1_check
          write(testout,101) ' diff        = ',
     &        cdabs(D0ms0ir1_coli/D0ms0ir1_check-1d0)
        end if
      
      endif  
#endif

      end

************************************************************************
      function D0ms0ir2_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     4-point function                                                 *
*     with 0 mass singularities and up to 2 soft singularities         *
*                                                                      *
*     assumes                                                          *
*     m12 small, p12=m22, p14=m42, m32=0                               *
*                                                                      *
*                        m22                                           *
*       p12=m22  ---------------------  p23=m22                        *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | m32 small                            *
*                     |         |                                      *
*                     |    4    |                                      *
*       p14=m42  ---------------------  p34=m42                        *
*                         m42                                          *
*                                                                      *
************************************************************************
*     22.08.08 Ansgar Denner        last changed  10.10.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms0ir2_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      logical flag(0:8)
      data    flag /9*.true./
      save    flag  
#endif

      real*8     ir24
      complex*16 ma2(4)
      complex*16 m2,m4,k24,r24
      complex*16 cln_coli
      integer    nsoft
      integer    i,j,k
      logical    soft(4),onsh(4,4)

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
#ifdef CHECK
 101  format(a22,g25.17)
      if(argcheck)then
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

        if(ms12.ne.cd0.or.ms32.ne.cd0.or.ms22.eq.cd0.or.ms42.eq.cd0
     &        .or.ps12.ne.ms22.or.ps14.ne.ms42 
     &        .or.ps23.ne.ms22.or.ps34.ne.ms42
     &        .or.ps24.eq.cd0.or.ps13.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &        .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0ms0ir2_coli called improperly:'
          write(*,111)' D0ms0ir2_coli: p12 = ',p12,ps12
          write(*,111)' D0ms0ir2_coli: p23 = ',p23,ps23
          write(*,111)' D0ms0ir2_coli: p34 = ',p34,ps34
          write(*,111)' D0ms0ir2_coli: p14 = ',p14,ps14
          write(*,111)' D0ms0ir2_coli: p13 = ',p13,ps13
          write(*,111)' D0ms0ir2_coli: p24 = ',p24,ps24
          write(*,111)' D0ms0ir2_coli: m12 = ',m12,ms12
          write(*,111)' D0ms0ir2_coli: m22 = ',m22,ms22
          write(*,111)' D0ms0ir2_coli: m32 = ',m32,ms32
          write(*,111)' D0ms0ir2_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif

      ma2(1)=m12
      ma2(2)=m22
      ma2(3)=m32
      ma2(4)=m42

c determine on-shell momenta
      onsh(1,2)=p12.eq.m22
      onsh(1,3)=p13.eq.m32
      onsh(1,4)=p14.eq.m42
      onsh(2,1)=p12.eq.m12
      onsh(2,3)=p23.eq.m32
      onsh(2,4)=p24.eq.m42
      onsh(3,1)=p13.eq.m12
      onsh(3,2)=p23.eq.m22
      onsh(3,4)=p34.eq.m42
      onsh(4,1)=p14.eq.m22
      onsh(4,2)=p24.eq.m22
      onsh(4,3)=p34.eq.m32

c count/determine soft singularities
      nsoft=0 
      do i=1,4
        do j=1,3
          if(i.ne.j)then
            do k=j,4
              if(k.ne.i.and.k.ne.j)then
                soft(i)=ma2(i).eq.cd0.and.onsh(i,j).and.onsh(i,k)
c                soft(i,k,j)=soft(i,j,k)
                if(soft(i)) nsoft=nsoft+1
              endif
            enddo
          endif
        enddo
      enddo

#ifdef CHECK
      if(flag(0))then
        write(*,100)' D0ms0ir2_coli: (3.8),(3.9),(3.10)',
     &        '    not yet tested at all '
        write(*,111)' D0ms0ir2_coli: p12 = ',p12
        write(*,111)' D0ms0ir2_coli: p23 = ',p23
        write(*,111)' D0ms0ir2_coli: p34 = ',p34
        write(*,111)' D0ms0ir2_coli: p14 = ',p14
        write(*,111)' D0ms0ir2_coli: p24 = ',p24
        write(*,111)' D0ms0ir2_coli: p13 = ',p13
        write(*,111)' D0ms0ir2_coli: m12 = ',m12
        write(*,111)' D0ms0ir2_coli: m22 = ',m22
        write(*,111)' D0ms0ir2_coli: m32 = ',m32
        write(*,111)' D0ms0ir2_coli: m42 = ',m42
        flag(0)=.false.
      endif
#endif

c according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 (ib)
c          and (iii)    (D0ir)
c Denner Dittmaier (3.8),(3.9),(3.10)

      m2 = sqrt(m22)
      m4 = sqrt(m42)
      k24 = (m22+m42-p24)/(m2*m4)
      if (k24.ne.0d0) then
c     |r24| > 0 yields wrong results !!!!!!! 
c     analytic continuation requires |r24| < 0
        r24 = 2d0/(k24*(1d0+sqrt(dcmplx(1d0-4d0/k24**2))))
      else
        r24 = dcmplx(0d0,1d0)
      end if
      if(dreal(k24).lt.-2d0) then
        ir24 = sign(1d1,1d0-abs(r24))
      else
        ir24 = 0d0
      end if
      
      if(p24.ne.(m2-m4)**2) then
        if(m12.eq.cd0.and.m32.eq.cd0) then 

#ifdef CHECK
          if(m32.ne.cd0.or.m42.ne.0d0.or.p34.ne.0d0.or.
     &        p14.eq.m12.or.p23.eq.m22)
     &        then
            write(*,100)' D0ms0ir2_coli:  (3.8)',
     &          '    case with 2 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms0ir2_coli: p12 = ',p12
            write(*,111)' D0ms0ir2_coli: p23 = ',p23
            write(*,111)' D0ms0ir2_coli: p34 = ',p34
            write(*,111)' D0ms0ir2_coli: p14 = ',p14
            write(*,111)' D0ms0ir2_coli: p24 = ',p24
            write(*,111)' D0ms0ir2_coli: p13 = ',p13
            write(*,111)' D0ms0ir2_coli: m12 = ',m12
            write(*,111)' D0ms0ir2_coli: m22 = ',m22
            write(*,111)' D0ms0ir2_coli: m32 = ',m32
            write(*,111)' D0ms0ir2_coli: m42 = ',m42
            stop
          endif
          
          if(flag(1))then
            write(*,100)' D0ms0ir2_coli:  (3.8)',
     &          '    case with 2 zero masses:',
     &          '    not yet tested at all'
            write(*,111)' D0ms0ir2_coli: p12 = ',p12
            write(*,111)' D0ms0ir2_coli: p23 = ',p23
            write(*,111)' D0ms0ir2_coli: p34 = ',p34
            write(*,111)' D0ms0ir2_coli: p14 = ',p14
            write(*,111)' D0ms0ir2_coli: p24 = ',p24
            write(*,111)' D0ms0ir2_coli: p13 = ',p13
            write(*,111)' D0ms0ir2_coli: m12 = ',m12
            write(*,111)' D0ms0ir2_coli: m22 = ',m22
            write(*,111)' D0ms0ir2_coli: m32 = ',m32
            write(*,111)' D0ms0ir2_coli: m42 = ',m42
            flag(1)=.false.
          endif
#endif

          D0ms0ir2_coli = 2d0*cln_coli(r24,ir24)
     &      *cln_coli(-p13/muir2,-1d0)
#ifdef SING
     &        +delta1ir*2d0*cln_coli(r24,ir24)
#endif
          D0ms0ir2_coli = D0ms0ir2_coli 
     &        *r24/(m2*m4*p13*(1-r24*r24))


        elseif(m12.eq.cd0.or.m32.eq.cd0) then 
#ifdef CHECK
          if(m32.ne.cd0.or.m42.ne.0d0.or.p34.ne.0d0.or.
     &        p14.eq.m12.or.p23.eq.m22)
     &        then
            write(*,100)' D0ms0ir2_coli:  (3.9)',
     &          '    case with 1 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms0ir2_coli: p12 = ',p12
            write(*,111)' D0ms0ir2_coli: p23 = ',p23
            write(*,111)' D0ms0ir2_coli: p34 = ',p34
            write(*,111)' D0ms0ir2_coli: p14 = ',p14
            write(*,111)' D0ms0ir2_coli: p24 = ',p24
            write(*,111)' D0ms0ir2_coli: p13 = ',p13
            write(*,111)' D0ms0ir2_coli: m12 = ',m12
            write(*,111)' D0ms0ir2_coli: m22 = ',m22
            write(*,111)' D0ms0ir2_coli: m32 = ',m32
            write(*,111)' D0ms0ir2_coli: m42 = ',m42
            stop
          endif
          
          if(flag(2))then
            write(*,100)' D0ms0ir2_coli:  (3.9)',
     &          '    case with 1 zero masses:',
     &          '    not yet tested at all'
            write(*,111)' D0ms0ir2_coli: p12 = ',p12
            write(*,111)' D0ms0ir2_coli: p23 = ',p23
            write(*,111)' D0ms0ir2_coli: p34 = ',p34
            write(*,111)' D0ms0ir2_coli: p14 = ',p14
            write(*,111)' D0ms0ir2_coli: p24 = ',p24
            write(*,111)' D0ms0ir2_coli: p13 = ',p13
            write(*,111)' D0ms0ir2_coli: m12 = ',m12
            write(*,111)' D0ms0ir2_coli: m22 = ',m22
            write(*,111)' D0ms0ir2_coli: m32 = ',m32
            write(*,111)' D0ms0ir2_coli: m42 = ',m42
            flag(2)=.false.
          endif
#endif
 
          D0ms0ir2_coli = cln_coli(r24,ir24)*
     &        (cln_coli(-p13/muir2,-1d0)+cln_coli(-p13/(m12+m32),-1d0))
     &        *r24/(m2*m4*p13*(1-r24*r24))
#ifdef SING
     &        +delta1ir*
     &        cln_coli(r24,ir24)*r24/(m2*m4*p13*(1-r24*r24))
#endif
        else
#ifdef CHECK
          if(m32.ne.cd0.or.m42.ne.0d0.or.p34.ne.0d0.or.
     &        p14.eq.m12.or.p23.eq.m22)
     &        then
            write(*,100)' D0ms0ir2_coli:  (3.10)',
     &          '    case with 0 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms0ir2_coli: p12 = ',p12
            write(*,111)' D0ms0ir2_coli: p23 = ',p23
            write(*,111)' D0ms0ir2_coli: p34 = ',p34
            write(*,111)' D0ms0ir2_coli: p14 = ',p14
            write(*,111)' D0ms0ir2_coli: p24 = ',p24
            write(*,111)' D0ms0ir2_coli: p13 = ',p13
            write(*,111)' D0ms0ir2_coli: m12 = ',m12
            write(*,111)' D0ms0ir2_coli: m22 = ',m22
            write(*,111)' D0ms0ir2_coli: m32 = ',m32
            write(*,111)' D0ms0ir2_coli: m42 = ',m42
            stop
          endif
          
          if(flag(3))then
            write(*,100)' D0ms0ir2_coli:  (3.10)',
     &          '    case with 0 zero masses:',
     &          '    not yet tested at all'
            write(*,111)' D0ms0ir2_coli: p12 = ',p12
            write(*,111)' D0ms0ir2_coli: p23 = ',p23
            write(*,111)' D0ms0ir2_coli: p34 = ',p34
            write(*,111)' D0ms0ir2_coli: p14 = ',p14
            write(*,111)' D0ms0ir2_coli: p24 = ',p24
            write(*,111)' D0ms0ir2_coli: p13 = ',p13
            write(*,111)' D0ms0ir2_coli: m12 = ',m12
            write(*,111)' D0ms0ir2_coli: m22 = ',m22
            write(*,111)' D0ms0ir2_coli: m32 = ',m32
            write(*,111)' D0ms0ir2_coli: m42 = ',m42
            flag(3)=.false.
          endif
#endif
          
          D0ms0ir2_coli = cln_coli(r24,ir24)*
     &        (cln_coli(-p13/m12,-1d0)+cln_coli(-p13/m32,-1d0))
     &        *r24/(m2*m4*p13*(1-r24*r24))

        endif
      else
        write(*,100)' D0ms0ir2_coli:',
     &        '    not implemented for det=0 '
        write(*,111)' D0ms0ir2_coli: p12 = ',p12
        write(*,111)' D0ms0ir2_coli: p23 = ',p23
        write(*,111)' D0ms0ir2_coli: p34 = ',p34
        write(*,111)' D0ms0ir2_coli: p14 = ',p14
        write(*,111)' D0ms0ir2_coli: p24 = ',p24
        write(*,111)' D0ms0ir2_coli: p13 = ',p13
        write(*,111)' D0ms0ir2_coli: m12 = ',m12
        write(*,111)' D0ms0ir2_coli: m22 = ',m22
        write(*,111)' D0ms0ir2_coli: m32 = ',m32
        write(*,111)' D0ms0ir2_coli: m42 = ',m42

      end if

      end

************************************************************************
      function D0ms1ir0_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     4-point function                                                 *
*     with 1 mass singularities and up to 0 soft singularities         *
*                                                                      *
*     assumes                                                          *
*     p34, m32, m42 small                                              *
*     p14=/=m12, p23=/=m22                                             *
*                                                                      *
*                        m22                                           *
*          p12  ---------------------  p23=/=m22                       *
*                     |    2    |                                      *
*                     |         |                                      *
*                 m12 | 1     3 | m32 small                            *
*                     |         |                                      *
*                     |    4    |                                      *
*     p14/==m12  ---------------------  p34 small                      *
*                      m42 small                                       *
*                                                                      *
************************************************************************
*     22.08.08 Ansgar Denner        last changed  06.11.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms1ir0_coli

#ifdef CHECK
      complex*16 elimminf2_coli,D0ms1ir0_check
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      real*8     test1,test2,test3
      logical flag(0:5),flag3,flag2
      data    flag /6*.true./, flag3 /.true./,  flag2 /.true./
      save    flag  
#endif

      complex*16 cln_coli,cspenc_coli,cspcos_coli,eta2s_coli
      complex*16 cspco0_coli
      complex*16 det,n12
      complex*16 m1,m2,m3,m4
      complex*16 mm12,mm22,mm32,mm42,q12,q13,q14,q23,q24,q34
      complex*16 h12,h23,h34,h14
      complex*16 k12,k13,k14,k23,k24,k34
      complex*16 l12,l13,l14,l23,l24,l34
      complex*16 r12,r13,r14,r23,r24,r34,r21,r31,r41,r32,r42,r43
      complex*16 invr12
      complex*16 a,b,c,d
      complex*16 x(2),y(2),z(2),x4(2),x1(2)
      complex*16 lqbard
      real*8     rd
      real*8     ir12,ir14,ir23,ir34,ir13,ir24,ik14,ik24,imp4,imp3
      real*8     iy(2),iz(2),ix4(2),ix1(2),ix(2),iqbard(2)
      real*8     rt12,rt34,rt14,rt23
      integer    i,j

      real*8     v,u
      complex*16 mm3,mm2
      complex*16 eta
      complex*16 ch1,ch2,ch3,l1(2),lx1,lx4,sp3,sp4

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
#ifdef CHECK
 101  format(a22,g25.17)
      if(argcheck)then
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

        if(ms32.ne.cd0.or.ms42.ne.cd0.or.ps34.ne.cd0
     &      .or.ps12.eq.cd0.and.ms12.eq.cd0.and.ms22.eq.cd0
     &      .or.ps14.eq.ms12.or.ps23.eq.ms22
c     &      .or.ps24.eq.cd0.or.ps13.eq.cd0
     &      .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &      .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &      .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0ms1ir0_coli called improperly:'
          write(*,111)' D0ms1ir0_coli: p12 = ',p12,ps12
          write(*,111)' D0ms1ir0_coli: p23 = ',p23,ps23
          write(*,111)' D0ms1ir0_coli: p34 = ',p34,ps34
          write(*,111)' D0ms1ir0_coli: p14 = ',p14,ps14
          write(*,111)' D0ms1ir0_coli: p13 = ',p13,ps13
          write(*,111)' D0ms1ir0_coli: p24 = ',p24,ps24
          write(*,111)' D0ms1ir0_coli: m12 = ',m12,ms12
          write(*,111)' D0ms1ir0_coli: m22 = ',m22,ms22
          write(*,111)' D0ms1ir0_coli: m32 = ',m32,ms32
          write(*,111)' D0ms1ir0_coli: m42 = ',m42,ms42
          write(*,*)' D0ms1ir0_coli: test= ',
     &        ms32.ne.cd0,ms42.ne.cd0,ps34.ne.cd0,
     &      ps12.eq.cd0.and.ms12.eq.cd0.and.ms22.eq.cd0,
     &      ps14.eq.ms12,ps23.eq.ms22,
c     &      ps24.eq.cd0,ps13.eq.cd0,
     &      dimag(p12).ne.0d0,dimag(p23).ne.0d0,
     &      dimag(p34).ne.0d0,dimag(p14).ne.0d0,
     &      dimag(p24).ne.0d0,dimag(p13).ne.0d0
          stop
        end if
      end if
#endif

      if(m32.eq.cd0.and.m42.eq.cd0)then
* 0 soft singularities and 2 zero masses
*
*                  m22             
*        p12 ----------------  p23    
*                |  2   |         
*            m12 |1    3| 0       
*                |   4  |         
*        p14 ----------------  0
*                   0             
*

#ifdef CHECK
        if(m32.ne.cd0.or.m42.ne.0d0.or.p34.ne.0d0.or.
     &      p14.eq.m12.or.p23.eq.m22)
     &      then
          write(*,100)' D0ms1ir0_coli:  (3.12)',
     &        '    case with 2 zero masses:',
     &        '    wrong arguments'
          write(*,111)' D0ms1ir0_coli: p12 = ',p12
          write(*,111)' D0ms1ir0_coli: p23 = ',p23
          write(*,111)' D0ms1ir0_coli: p34 = ',p34
          write(*,111)' D0ms1ir0_coli: p14 = ',p14
          write(*,111)' D0ms1ir0_coli: p13 = ',p13
          write(*,111)' D0ms1ir0_coli: p24 = ',p24
          write(*,111)' D0ms1ir0_coli: m12 = ',m12
          write(*,111)' D0ms1ir0_coli: m22 = ',m22
          write(*,111)' D0ms1ir0_coli: m32 = ',m32
          write(*,111)' D0ms1ir0_coli: m42 = ',m42
          stop
        endif
#endif
       
#ifdef ALLCHECK 
c ppttbb dim
        if(flag(1))then
          write(*,100)' D0ms1ir0_coli:  (3.12)',
     &        '    case with 2 zero masses:',
     &        '    not yet tested in physical process'
          write(*,111)' D0ms1ir0_coli: p12 = ',p12
          write(*,111)' D0ms1ir0_coli: p23 = ',p23
          write(*,111)' D0ms1ir0_coli: p34 = ',p34
          write(*,111)' D0ms1ir0_coli: p14 = ',p14
          write(*,111)' D0ms1ir0_coli: p13 = ',p13
          write(*,111)' D0ms1ir0_coli: p24 = ',p24
          write(*,111)' D0ms1ir0_coli: m12 = ',m12
          write(*,111)' D0ms1ir0_coli: m22 = ',m22
          write(*,111)' D0ms1ir0_coli: m32 = ',m32
          write(*,111)' D0ms1ir0_coli: m42 = ',m42
          flag(1)=.false.
        endif
#endif

        if(m12.ne.cd0)then
          mm22=m22
          mm12=m12
          q14=p14
          q23=p23
          q13=p13
          q24=p24
        else
          mm22=m12
          mm12=m22
          q14=p23
          q23=p14
          q13=p24
          q24=p13
        endif          

c according to Ellis, Zanderighi, JHEP 0802 (2008) 002   box13
c p2^2 -> p14,  p3^2 -> p12,  p4^2 -> p23,  s12 -> p143, s23 -> p24
c m3^2 -> m12, m4_^2 -> m22 
c Denner Dittmaier (3.12)

        det=(mm12-q13)*(mm22-q24)-(mm12-q14)*(mm22-q23)

        if(det.eq.cd0)then
          write(*,100)' D0ms1ir0_coli:  (3.12)',
     &        '    case with 2 zero masses',
     &        '    zero Cayley determinant not implemented'
          write(*,111)' D0ms1ir0_coli: p12 = ',p12
          write(*,111)' D0ms1ir0_coli: p23 = ',p23
          write(*,111)' D0ms1ir0_coli: p34 = ',p34
          write(*,111)' D0ms1ir0_coli: p14 = ',p14
          write(*,111)' D0ms1ir0_coli: p13 = ',p13
          write(*,111)' D0ms1ir0_coli: p24 = ',p24
          write(*,111)' D0ms1ir0_coli: m12 = ',m12
          write(*,111)' D0ms1ir0_coli: m22 = ',m22
          write(*,111)' D0ms1ir0_coli: m32 = ',m32
          write(*,111)' D0ms1ir0_coli: m42 = ',m42
          stop
        endif

        if(p12.ne.cd0)then
          l12 = (mm12+mm22-p12)      
           if (l12.ne.0d0) then
c            r12 = l12/(2d0*mm22)*(1d0+sqrt(1d0-4d0*mm12*mm22/l12**2))
            n12 = .5d0*l12*(1d0+sqrt(1d0-4d0*mm12*mm22/l12**2))
          else
            n12 = sqrt(mm12*mm22)*dcmplx(0d0,1d0)
          end if
        else
          if(dreal(mm12).gt.dreal(mm22))then
            n12=mm12
          else
            n12=mm22
          endif
        endif
        if(n12.ne.cd0) then
          r12=mm12/n12
          r21=mm22/n12
          if(m12.ne.cd0) invr12 = 1d0/r12
          ir12 = sign(1d0,1d0-abs(r12*r21))
        else
          invr12 = cd0
          r21 = (mm12-p12)/mm12
          ir12 = -1d0
        endif   
        
        if(abs(p12-mm12).gt.abs(mm12))then
c case for mm22=0 and mm12=0
          D0ms1ir0_coli =
#ifdef SING
     &      delta1ir*(cln_coli((mm12-q14)/(mm12-q13),dreal(q13-q14))
     &      +cln_coli((mm22-q23)/(mm22-q24),dreal(q24-q23)))
#endif
     &      -2d0*cspenc_coli(1d0-(mm12-q14)/(mm12-q13),
     &        dreal(q14-q13))
     &      -2d0*cspenc_coli(1d0-(mm22-q23)/(mm22-q24),
     &        dreal(q23-q24))
     &      +2d0*cspcos_coli((mm12-q14)/(mm12-q13),
     &        (mm22-q23)/(mm22-q24),
     &              -dreal(q14-q13),-dreal(q23-q24))
     &      +cspcos_coli((mm22-q24)/(mm12-q14),r12,
     &              -dreal(mm12-q14-mm22+q24),ir12)
     &      -cspcos_coli((mm12-q14)/(mm22-q24),r21,
     &              dreal(mm12-q14-mm22+q24),ir12)
     &      +cspcos_coli((mm12-q13)/(mm22-q23),r21,
     &              -dreal(mm22-q23-mm12+q13),ir12)
     &      -cspcos_coli((mm22-q23)/(mm12-q13),r12,
     &              dreal(mm22-q23-mm12+q13),ir12)
     &      +2d0*cln_coli((mm12-q13)/muir2,-1d0)
     &      *cln_coli((mm22-q24)/muir2,-1d0)
     &      -(cln_coli((mm12-q14)/muir2,-1d0)**2)
     &      -(cln_coli((mm22-q23)/muir2,-1d0)**2)
     &      +cln_coli((mm12-q14)/(mm22-q24),dreal(mm12-q14-mm22+q24))
     &      *cln_coli(n12/muir2,-ir12)
     &      +cln_coli((mm22-q23)/(mm12-q13),dreal(mm22-q23-mm12+q13))
     &      *cln_coli(n12/muir2,-ir12)
     &      +.5d0*(cln_coli((mm12-q13)/(mm22-q23),
     &      -dreal(mm22-q23-mm12+q13)))**2
     &      +.5d0*(cln_coli((mm22-q24)/(mm12-q14),
     &      -dreal(mm12-q14-mm22+q24)))**2

        write(*,*) 'D0ms1ir0_coli a = ',D0ms1ir0_coli

        else
c case for mm22=0 and p12=mm12
          D0ms1ir0_coli =
#ifdef SING
     &      delta1ir*(cln_coli((mm12-q14)/(mm12-q13),dreal(q13-q14))
     &      +cln_coli((mm22-q23)/(mm22-q24),dreal(q24-q23)))
#endif
     &      -2d0*cspenc_coli(1d0-(mm12-q14)/(mm12-q13),
     &      dreal(q14-q13))
     &      -2d0*cspenc_coli(1d0-(mm22-q23)/(mm22-q24),
     &      dreal(q23-q24))
     &      +2d0*cspcos_coli((mm12-q14)/(mm12-q13),
     &      (mm22-q23)/(mm22-q24),
     &              -dreal(q14-q13),-dreal(q23-q24))
     &      -cspcos_coli((mm12-q14)/(mm22-q24),invr12,
     &              dreal(mm12-q14-mm22+q24),-ir12)
     &      -cspcos_coli((mm12-q14)/(mm22-q24),r21,
     &              dreal(mm12-q14-mm22+q24),ir12)
     &      +cspcos_coli((mm12-q13)/(mm22-q23),r21,
     &              -dreal(mm22-q23-mm12+q13),ir12)
     &      +cspcos_coli((mm12-q13)/(mm22-q23),invr12,
     &              -dreal(mm22-q23-mm12+q13),ir12)
     &      +2d0*cln_coli((mm12-q13)/muir2,-1d0)
     &      *cln_coli((mm22-q24)/muir2,-1d0)
     &      -(cln_coli((mm12-q14)/muir2,-1d0)**2)
     &      -(cln_coli((mm22-q23)/muir2,-1d0)**2)
     &      +cln_coli((mm12-q14)/(mm22-q24),dreal(mm12-q14-mm22+q24))
     &      *cln_coli(mm12/muir2,-ir12)
     &      +cln_coli((mm22-q23)/(mm12-q13),dreal(mm22-q23-mm12+q13))
     &      *cln_coli(mm12/muir2,-ir12)
     &      +(cln_coli((mm22-q23)/(mm12-q13),
     &      -dreal(mm12-q13-mm22+q23)))**2
  
c        write(*,*) 'D0ms1ir0_coli b = ',D0ms1ir0_coli
        endif

        D0ms1ir0_coli = D0ms1ir0_coli/det

#ifdef CHECK
        if(conscheck)then
          if(mm12.ne.cd0.and.mm22.ne.cd0)then
            D0ms1ir0_check =
#ifdef SING
     &      delta1ir*(cln_coli((mm12-q14)/(mm12-q13),dreal(q13-q14))
     &      +cln_coli((mm22-q23)/(mm22-q24),dreal(q24-q23)))
#endif
     &      -2d0*cspenc_coli(1d0-(mm12-q14)/(mm12-q13),
     &          dreal(q14-q13))
     &      -2d0*cspenc_coli(1d0-(mm22-q23)/(mm22-q24),
     &          dreal(q23-q24))
     &      +2d0*cspcos_coli((mm12-q14)/(mm12-q13),
     &          (mm22-q23)/(mm22-q24),
     &              -dreal(q14-q13),-dreal(q23-q24))
     &      -cspcos_coli((mm12-q14)/(mm22-q24),1d0/r12,
     &              dreal(mm12-q14-mm22+q24),-ir12)
     &      -cspcos_coli((mm12-q14)/(mm22-q24),r21,
     &              dreal(mm12-q14-mm22+q24),ir12)
     &      -cspcos_coli((mm22-q23)/(mm12-q13),1d0/r21,
     &              dreal(mm22-q23-mm12+q13),-ir12)
     &      -cspcos_coli((mm22-q23)/(mm12-q13),r12,
     &              dreal(mm22-q23-mm12+q13),ir12)
     &      +2d0*cln_coli((mm12-q13)/muir2,-1d0)
     &      *cln_coli((mm22-q24)/muir2,-1d0)
     &      -(cln_coli((mm12-q14)/muir2,-1d0)**2)
     &      -(cln_coli((mm22-q23)/muir2,-1d0)**2)
     &      +cln_coli((mm12-q14)/(mm22-q24),
     &          dreal(mm12-q14-mm22+q24))
     &      *log(mm12/muir2)
     &      +cln_coli((mm22-q23)/(mm12-q13),
     &          dreal(mm22-q23-mm12+q13))
     &      *log(mm22/muir2)
     &      -.5d0*(cln_coli(r12,ir12))**2
     &      -.5d0*(cln_coli(r21,ir12))**2
 
            write(*,*) 'D0ms1ir0_coli c = ',D0ms1ir0_coli

            D0ms1ir0_check = D0ms1ir0_check/det

            if(cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0).gt.testacc) then
              write(testout,100) ' D0ms1ir0_coli inconsistent a'
              write(testout,111) ' D0ms1ir0_coli: p12 = ',p12
              write(testout,111) ' D0ms1ir0_coli: p23 = ',p23
              write(testout,111) ' D0ms1ir0_coli: p34 = ',p34
              write(testout,111) ' D0ms1ir0_coli: p14 = ',p14
              write(testout,111) ' D0ms1ir0_coli: p24 = ',p24
              write(testout,111) ' D0ms1ir0_coli: p13 = ',p13
              write(testout,111) ' D0ms1ir0_coli: m12 = ',m12
              write(testout,111) ' D0ms1ir0_coli: m22 = ',m22
              write(testout,111) ' D0ms1ir0_coli: m32 = ',m32
              write(testout,111) ' D0ms1ir0_coli: m42 = ',m42
              write(testout,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
              write(testout,111) ' D0ms1ir0_check = ',D0ms1ir0_check
              write(testout,101) ' diff        = ',
     &            cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
              
              write(*,100)   ' D0ms1ir0_coli inconsistent'
              write(*,111) ' D0ms1ir0_coli: p12 = ',p12
              write(*,111) ' D0ms1ir0_coli: p23 = ',p23
              write(*,111) ' D0ms1ir0_coli: p34 = ',p34
              write(*,111) ' D0ms1ir0_coli: p14 = ',p14
              write(*,111) ' D0ms1ir0_coli: p24 = ',p24
              write(*,111) ' D0ms1ir0_coli: p13 = ',p13
              write(*,111) ' D0ms1ir0_coli: m12 = ',m12
              write(*,111) ' D0ms1ir0_coli: m22 = ',m22
              write(*,111) ' D0ms1ir0_coli: m32 = ',m32
              write(*,111) ' D0ms1ir0_coli: m42 = ',m42
              write(*,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
              write(*,111) ' D0ms1ir0_check = ',D0ms1ir0_check
              write(*,101) ' diff        = ',
     &            cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
              stop
            endif
          endif
        endif
#endif

      elseif(m32.eq.cd0.or.m42.eq.cd0)then
* 0 soft singularities and 1 zero masses
*
*                  m22             
*        p12 ----------------  p23    
*                |  2   |         
*            m12 |1    3| 0       
*                |   4  |         
*        p14 ----------------  p34=m42
*                  m42             
*

#ifdef CHECK
        if(.not.(m32.eq.cd0.and.m42.eq.p34.or.
     &      m42.eq.cd0.or.m32.eq.p34).or.
     &      p14.eq.m12.or.p23.eq.m22)
     &      then
          write(*,100)' D0ms1ir0_coli:  (3.13)',
     &        '    case with 1 zero masses:',
     &        '    wrong arguments'
          write(*,111)' D0ms1ir0_coli: p12 = ',p12
          write(*,111)' D0ms1ir0_coli: p23 = ',p23
          write(*,111)' D0ms1ir0_coli: p34 = ',p34
          write(*,111)' D0ms1ir0_coli: p14 = ',p14
          write(*,111)' D0ms1ir0_coli: p24 = ',p24
          write(*,111)' D0ms1ir0_coli: p13 = ',p13
          write(*,111)' D0ms1ir0_coli: m12 = ',m12
          write(*,111)' D0ms1ir0_coli: m22 = ',m22
          write(*,111)' D0ms1ir0_coli: m32 = ',m32
          write(*,111)' D0ms1ir0_coli: m42 = ',m42
          stop
        endif
#endif
        
        if(m32.eq.cd0)then
          if(m12.ne.cd0)then
            mm42=m42
            mm12=m12
            mm22=m22
            q13=p13
            q23=p23
            q14=p14
            q24=p24
         else
            mm42=m42
            mm12=m22
            mm22=m12
            q13=p23
            q23=p13
            q14=p24
            q24=p14        
         endif
        else
          if(m22.ne.cd0)then
            mm42=m32
            mm12=m22
            mm22=m12
            q13=p24
            q23=p14
            q14=p23
            q24=p13
          else
            mm42=m32
            mm12=m12
            mm22=m22
            q13=p14
            q23=p24
            q14=p13
            q24=p23
          endif
        endif

        if (mm12.ne.0d0.and.mm22.ne.0d0) then
c cD0ms0e
c Denner Dittmaier (3.13)

#ifdef CHECK
          if(flag(2))then
            write(*,100)' D0ms1ir0_coli:  (3.13), m12,m22 =/=0',
     &          '    case with 1 zero masses:',
     &          '    not yet tested at all'
            write(*,111)' D0ms1ir0_coli: p12 = ',p12
            write(*,111)' D0ms1ir0_coli: p23 = ',p23
            write(*,111)' D0ms1ir0_coli: p34 = ',p34
            write(*,111)' D0ms1ir0_coli: p14 = ',p14
            write(*,111)' D0ms1ir0_coli: p24 = ',p24
            write(*,111)' D0ms1ir0_coli: p13 = ',p13
            write(*,111)' D0ms1ir0_coli: m12 = ',m12
            write(*,111)' D0ms1ir0_coli: m22 = ',m22
            write(*,111)' D0ms1ir0_coli: m32 = ',m32
            write(*,111)' D0ms1ir0_coli: m42 = ',m42
            write(*,111)' D0ms1ir0_coli: q23 = ',q23
            write(*,111)' D0ms1ir0_coli: q14 = ',q14
            write(*,111)' D0ms1ir0_coli: q24 = ',q24
            write(*,111)' D0ms1ir0_coli: q13 = ',q13
            write(*,111)' D0ms1ir0_coli: mm12= ',mm12
            write(*,111)' D0ms1ir0_coli: mm22= ',mm22
            write(*,111)' D0ms1ir0_coli: mm32= ',mm32
            write(*,111)' D0ms1ir0_coli: mm42= ',mm42
            flag(2)=.false.
          endif
#endif

          m2  = sqrt(mm22)
          m1  = sqrt(mm12)
          m4  = sqrt(mm42)

c   (permutation 1 -> 3 -> 2 -> 4 -> 1 for  m -> k )
          k12 = (mm22+mm12-p12)/(m1*m2)
          k13 = (m12    -q13)/(m1*m4)
          k14 = (    mm12-q14)/(m1*m4*coliminfscale)
          k23 = (mm22    -q23)/(m2*m4)
c         k14 = (    mm12-q14)/(m1*m4)
c         k23 = (mm22    -q23)/(m2*m4*coliminfscale)
          k24 = (mm22    -q24)/(m2*m4*coliminfscale)

          r12 = k12/2d0*(1d0+sqrt(1d0-4d0/k12**2))
          if(dreal(k12).lt.-2D0) then
            ir12 = sign(1d1,1d0-abs(r12))
          else
            ir12 = 0D0
          end if
       
C>      x1 = -k14
C>      x2 = -k13*k24/k23
C>      ix1 = 1d0
C>      ix2 = 0d0

c      write(80,*) 'D0ms1ir0_coli k ',k12,k14,k24
c      write(80,*) 'D0ms1ir0_coli k ',k13,k23
c      write(80,*) 'D0ms1ir0_coli r ',r12,k14,k24
c     write(*,*) 'D0ms1ir0_coli x4 ',x1,x2
c     write(*,*) 'D0ms1ir0_coli x1 ',x1/k24,x2/k24

C>      D0ms1ir0_coli = -cspcos_coli(-x1,k14,-ix1,-1d0) 
C>     &           - cspcos_coli(-x1,1d0/k14,-ix1,1d0)
C>     &           + cspcos_coli(-x1/k24,r12,-ix1/dreal(k24),ir12)
C>     &           + cspcos_coli(-x1/k24,1d0/r12,-ix1/dreal(k24),-ir12)
C>     &           + pi2_6
C>     &           - cspcos_coli(-x1/k24,k23/k13,-ix1/dreal(k24),0d0)
C>     &           + eta2s_coli(-x1,1d0/k24,-ix1,1d0,ix1)
C>     &             *(cln_coli( (k24/x1+k12+x1/k24)/(k23+k13*k24/x1))
C>     &               + cln_coli(k23))
C>     &           + cspcos_coli(-x2,k14,-ix2,-1d0) 
C>     &           + cspcos_coli(-x2,1d0/k14,-ix2,1d0)
C>     &           - cspcos_coli(-x2/k24,r12,-ix2/dreal(k24),ir12)
C>     &           - cspcos_coli(-x2/k24,1d0/r12,-ix2/dreal(k24),-ir12)
C>     &           - pi2_6
C>     &           + cspcos_coli(-x2/k24,k23/k13,-ix2/dreal(k24),0d0)
C>     &           - eta2s_coli(-x2,1d0/k24,-ix2,1d0,ix2)
C>     &             *(cln_coli(-k24*(k14+x2)/k23)+cln_coli(k13))
C>
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &       -cspcos_coli(-x1,k14,-ix1,-1d0) 
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           - cspcos_coli(-x1,1d0/k14,-ix1,1d0)
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           + cspcos_coli(-x1/k24,r12,-ix1/dreal(k24),ir12)
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           + cspcos_coli(-x1/k24,1d0/r12,-ix1/dreal(k24),-ir12)
C>c      write(*,*) 'D0ms1ir0_coli contr ',
C>c     &           + pi2_6
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           - cspcos_coli(-x1/k24,k23/k13,-ix1/dreal(k24),0d0)
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           + eta2s_coli(-x1,1/k24,-ix1,1d0,ix1)
C>     &             *(cln_coli( (k24/x1+k12+x1/k24)/(k23+k13*k24/x1))
C>     &               + cln_coli(k23))
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           + cspcos_coli(-x2,k14,-ix2,-1d0) 
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           + cspcos_coli(-x2,1d0/k14,-ix2,1d0)
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           - cspcos_coli(-x2/k24,r12,-ix2/k24,ir12)
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           - cspcos_coli(-x2/k24,1d0/r12,-ix2/dreal(k24),-ir12)
C>c      write(*,*) 'D0ms1ir0_coli contr ',
C>c     &           - pi2_6
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           + cspcos_coli(-x2/k24,k23/k13,-ix2/dreal(k24),0d0)
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           - eta2s_coli(-x2,dcmplx(1/k24),-ix2,1d0,ix2)
C>     &             *(cln_coli(-k24*(k14+x2)/k23)+cln_coli(k13))
C>
C>      write(*,*) 'D0ms1ir0_coli sum ',D0ms1ir0_coli
C>      D0ms1ir0_coli = D0ms1ir0_coli/(m22*m3*m4*(k14*k23-k13*k24))
C>
C>      write(*,*) 'D0ms1ir0_coli =   ',D0ms1ir0_coli
C>
C>      D0ms1ir0_coli = -cspcos_coli(k14,k14,-1d0,-1d0) 
C>c     &           - cspcos_coli(k14),dcmplx(1d0/k14),-1d0,1d0)
C>     &           + cspcos_coli(k14/k24,r12,-1d0/k24,ir12)
C>     &           + cspcos_coli(k14/k24,1d0/r12,-1d0/k24,-ir12)
C>c     &           + pi2_6
C>     &           - cspcos_coli(k14/k24,k23/k13,-1d0/k24,0d0)
C>c     &           + eta2s_coli(k14,1d0/k24,-1d0,1d0,ix1)
C>c     &             *(cln_coli( (k24/x1+k12+x1/k24)/(k23+k13*k24/x1))
C>c     &               + cln_coli(k23))
C>     &           + cspcos_coli(k13*k24/k23,dcmplx(k14),0d0,-1d0) 
C>     &           + cspcos_coli(k13*k24/k23,dcmplx(1d0/k14),0d0,1d0)
C>     &           - cspcos_coli(k13*k24/k23/k24,r12,0d0/k24,ir12)
C>     &           - cspcos_coli(k13*k24/k23/k24,1d0/r12,0d0/k24,-ir12)
C>c     &           - pi2_6
C>c     &           + cspcos_coli(k13*k24/k23/k24,k23/k13,0d0/k24,0d0)
C>     &           - eta2s_coli(k13*k24/k23,dcmplx(1/k24),0d0,1d0,ix2)
C>     &             *(cln_coli(-k24*(k14+x2)/k23)+cln_coli(k13))
C>
C>      write(*,*) 'D0ms1ir0_coli1 contr ',
C>     &       -cspcos_coli(dcmplx(k14),dcmplx(k14),-1d0,-1d0),
C>     &       pi2_6 + 2d0*cln_coli(k14,-1d0)**2 
C>c      write(*,*) 'D0ms1ir0_coli1 contr ',
C>c     &           - cspcos_coli(dcmplx(k14),dcmplx(1d0/k14),-1d0,1d0)
C>      write(*,*) 'D0ms1ir0_coli1 contr ',
C>     &           + cspcos_coli(dcmplx(k14)/k24,r12,-1d0/k24,ir12)
C>      write(*,*) 'D0ms1ir0_coli1 contr ',
C>     &           + cspcos_coli(dcmplx(k14)/k24,1d0/r12,-1d0/k24,-ir12)
C>c      write(*,*) 'D0ms1ir0_coli1 contr ',
C>c     &           + pi2_6
C>      write(*,*) 'D0ms1ir0_coli1 contr ',
C>     &           - cspcos_coli(dcmplx(k14)/k24,k23/k13,-1d0/k24,0d0),
C>     &           - cspcos_coli(dcmplx(k14)/k24,k23/k13,k14-k24,0d0)
C>c      write(*,*) 'D0ms1ir0_coli1 contr ',
C>c     &           + eta2s_coli(dcmplx(k14),dcmplx(1/k24),-1d0,1d0,ix1)
C>c     &             *(cln_coli( (k24/x1+k12+x1/k24)/(k23+k13*k24/x1))
C>c     &               + cln_coli(k23))
C>      write(*,*) 'D0ms1ir0_coli1 contr ',
C>     &           + cspcos_coli(k13*k24/k23,dcmplx(k14),0d0,-1d0),
C>     &     - pi2_6 - 0.5d0*(cln_coli(k13*k24/k23, )+cln_coli(k14,-1d0))**2
C>      write(*,*) 'D0ms1ir0_coli1 contr ',
C>     &           + cspcos_coli(k13*k24/k23,dcmplx(1d0/k14),0d0,1d0),
C>     &           + cspcos_coli(k13/k23,dcmplx(k24/k14),0d0,k24-k14)
C>     &           + eta2s_coli(k13*k24/k23,dcmplx(1/k24),0d0,1d0,0d0)
C>     &            *cdlog(1-k13*k24/(k23*k14))
C>      write(*,*) 'D0ms1ir0_coli1 contr ',
C>     &           - cspcos_coli(k13*k24/k23/k24,r12,0d0/k24,ir12)
C>      write(*,*) 'D0ms1ir0_coli1 contr ',
C>     &           - cspcos_coli(k13*k24/k23/k24,1d0/r12,0d0/k24,-ir12)
C>c      write(*,*) 'D0ms1ir0_coli1 contr ',
C>c     &           - pi2_6
C>c      write(*,*) 'D0ms1ir0_coli1 contr ',
C>c     &           + cspcos_coli(k13*k24/k23/k24,k23/k13,0d0/k24,0d0)
C>      write(*,*) 'D0ms1ir0_coli1 contr ',
C>     &           - eta2s_coli(k13*k24/k23,dcmplx(1/k24),0d0,1d0,0d0)
C>     &             *(cln_coli(-k24*(k14+x2)/k23)+cln_coli(k13))
C>
C>      write(*,*) 'D0ms1ir0_coli1 sum ',D0ms1ir0_coli
C>      D0ms1ir0_coli = D0ms1ir0_coli/(m22*m3*m4*(k14*k23-k13*k24))
C>
C>      write(*,*) 'D0ms1ir0_coli1 =   ',D0ms1ir0_coli

          D0ms1ir0_coli = 
     &        2d0*cspcos_coli(k24/k14,k13/k23,
     &          dreal(k24-k14),dreal(k13-k23)) 
     &        + cspcos_coli(k14/k24,r12,dreal(k14-k24),ir12)
     &        + cspcos_coli(k14/k24,1d0/r12,dreal(k14-k24),-ir12)
     &        - cspcos_coli(k13/k23,r12,dreal(k13-k23),ir12)
     &        - cspcos_coli(k13/k23,1d0/r12,dreal(k13-k23),-ir12)
     &        - 2d0*cln_coli(k14,-1d0)*
     &        (cln_coli(k13/k14,dreal(k13-k14))
     &        - cln_coli(k23/k24,dreal(k23-k24)))
c OMISSION correct??? not tested  19.05.04
cc     -olieta(k13*k24/k23,1/k24)* ...
c      if(dreal(k24).lt.0d0.and.dimag(k13*k24/k23).gt.0d0) then
c        D0ms1ir0_coli = D0ms1ir0_coli - dcmplx(0d0,-2d0*pi)*
c     &       (cdlog(-k24*(k14*k23-k13*k24)/(k23*k23))
c     &       -cdlog((k14*k23-k13*k24)/(k23*k14))
c     &       + cdlog(k23/k14) - cln_coli(k24,-1d0) )
c     &       + 12d0*pi2_6   
c      end if
       
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &           2d0*cspcos_coli(dcmplx(k24/k14),k13/k23,dreal(k24-k14),0d0) 
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &     + cspcos_coli(dcmplx(k14/k24),r12,dreal(k14-k24),ir12)
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &     + cspcos_coli(dcmplx(k14/k24),1d0/r12,dreal(k14-k24),-ir12)
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &     - cspcos_coli(k13/k23,r12,dreal(k13-k23),ir12)
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &     - cspcos_coli(k13/k23,1d0/r12,dreal(k13-k23),-ir12)
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &     - 2d0*cln_coli(k14,-1d0)*
C>     &     (cdlog(k13/k14)-cdlog(k23/k24))
C>c     -olieta(k13*k24/k23,1/k24)* ...
C>      if(k24.lt.0d0.and.dimag(k13*k24/k23).gt.0d0) then
C>      write(*,*) 'D0ms1ir0_coli contr ',
C>     &                      - dcmplx(0d0,-2d0*pi)*
C>     &       (cdlog(-k24*(k14*k23-k13*k24)/(k23*k23))
C>     &       -cdlog((k14*k23-k13*k24)/(k23*k14))
C>     &       + cdlog(k23/k14) - cln_coli(k24,-1d0) )
C>     &       + 12d0*pi2_6   
C>      end if
C>
C>      write(*,*) 'D0ms1ir0_coli sum ',D0ms1ir0_coli

          D0ms1ir0_coli = D0ms1ir0_coli/
     &        (m2*coliminfscale*m42*m1*(k14*k23-k13*k24))

c      write(*,*) 'D0ms1ir0_coli =   ',D0ms1ir0_coli              

#ifdef CHECK
          if (conscheck) then
            
            m1  = sqrt(mm12)
            m2  = sqrt(mm22)
            m4  = sqrt(mm42)
            mm3  = m1
        
            k12 = (mm12+mm22-p12)/(m1*m2)
            k13 = (mm12     -q13)/(m1*mm3)
            k14 = (mm12     -q14)/(m1*m4*coliminfscale)
            k23 = (mm22     -q23)/(m2*mm3)
            k24 = (mm22     -q24)/(m2*m4*coliminfscale)
c     k34 = 0d0
        
            r12 = k12/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k12**2)))
c     1/r24 = k24, r24 = 0
c     a   =   -k23
c     b   =   k13*k24-k14*k23
c     c   =   k12*k13-k13/k24*k14-k23
c     d   =   k23
c     det =  k14*k23 - k24*k13
c     x4(1) = (-b+det)/2d0/a
c     x4(2) = (-b-det)/2d0/a
            x4(1) = (k13*k24-k14*k23)/k23
            x4(2) = -(k12*k13-k13/k24*k14-k23)/(k23*x4(1))
            x1(1) = x4(1)*k24
            x1(2) = x4(2)*k24
            
            if(dreal(k12).lt.-2d0) then
              ir12 = sign(1d1,1d0-abs(r12))
            else
              ir12 = 0d0
            end if
            if(dreal(k14).lt.-2d0) then
              ik14 = -1d0
            else
              ik14 = 0d0
            end if
            if(dreal(k24).lt.-2d0) then
              ik24 = -1d0
            else              
              ik24 = 0d0
            end if
            ix4(1) =  sign(1d0,dreal(x4(1)))
            ix4(2) =  -ix4(1)
            ix1(1) =  sign(1d0,ix4(1)/dreal(k24))
            ix1(2) =  sign(1d0,ix4(2)/dreal(k24))

C>      WRITE(testout,*) 'D0ms1ir0_check Q  ',Q12,Q23,Q34
C>      WRITE(testout,*) 'D0ms1ir0_check Q  ',Q14,Q24,Q13
C>      WRITE(testout,*) 'D0ms1ir0_check m  ',m12,m2,m4
C>      WRITE(testout,*) 'D0ms1ir0_check K  ',K12,K13,K14
C>      WRITE(testout,*) 'D0ms1ir0_check K  ',K23,K24
C>      write(testout,*) 'D0ms1ir0_check r  ',r12,1/k14,1/k24
C>      write(testout,*) 'D0ms1ir0_check ir ',ir12
C>      write(testout,*) 'D0ms1ir0_check x4 ',x4(1),x4(2)
C>      write(testout,*) 'D0ms1ir0_check ix4',ix4(1),ix4(2)
C>      write(testout,*) 'D0ms1ir0_check x1 ',x1(1),x1(2)
C>      write(testout,*) 'D0ms1ir0_check ix1',ix1(1),ix1(2)

            lx1 = cln_coli(-x1(1),-ix1(1))
            lx4 = cln_coli(-x4(1),-ix4(1))
            D0ms1ir0_check =
     &          - cspcos_coli(-x4(1),1d0/k14,-ix4(1),-ik14)
     &          +(lx4+cln_coli(k14,ik14))**2/2d0
     &          -(lx1+cln_coli(r12,ir12))**2/2d0
     &          -(lx1+cln_coli(1d0/r12,-ir12))**2/2d0
     &          +(lx1+cln_coli(k23/k13,-dreal(k13)))**2/2d0
c to be checked
c     &          + ettile(-x4(1),k24,-ix4(1),ik24) *
     &          +eta2s_coli(-x4(1),k24,-ix4(1),ik24,-ix4(1)*dreal(k24))
     &          *(cln_coli(k24/k23*x4(1),ix4(1)*dreal(k24/k23))
     &          +cln_coli(k13,-1d0))

C>        write(testout,*) 'D0ms1ir0_check contr nn ', 
C>     &       lx4,cln_coli(dcmplx(k14),Ik14),x4(1),ix4(1)
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>     &       - cspcos_coli(-x4(1),1d0/k14,-ix4(1),-ik14)
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>     &       + pi2_6 +(lx4+cln_coli(k14,ik14))**2/2d0
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>     &       - pi2_6 -(lx1+cln_coli(r12,ir12))**2/2d0
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>     &       - pi2_6 -(lx1+cln_coli(1d0/r12,-ir12))**2/2d0
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>     &       + pi2_6 +(lx1+cln_coli(k23/k13,-k13))**2/2d0
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>c     &       + ettile(-dcmplx(x4(1)),dcmplx(k24),-ix4(1),ik24)
C>     &          +eta2s_coli(-x4(1),k24,-ix4(1),ik24,-ix4(1)*dreal(k24))
C>     &         *(cln_coli(k24/k23*x4(1),k24*ix4(1)/k23)
C>     &                   +cln_coli(k13,-1d0))

            ch1    = k12-k14/k24+k24*x4(2)
            test1 = abs(ch1)/
     &          max(abs(k12),abs(k24*x4(2)),abs(k14/k24))
            ch2    = 1d0/k24/x4(2)+k12+x4(2)*k24
            test2 = abs(ch2)/max(abs(1d0/k24/x4(2)),
     &          abs(k12),abs(x4(2)*k24))
            ch3    = k23+k13/k24/x4(2)
            test3 = abs(ch3)/max(abs(k23),abs(k13/k24/x4(2)))
            test2 = min(test2,test3)
            if (test1.lt.test2) then
              l1(2) = cln_coli(ch2/ch3,ix4(2)*dreal(k24/k23))
            else
              l1(2) = cln_coli(ch1/k23,ix4(2)*dreal(k24/k23))
            end if
            
            D0ms1ir0_check = D0ms1ir0_check 
     &          + pi2_6
     &          + cspcos_coli(-x4(2),k14,-ix4(2),ik14)
     &          - cspcos_coli(-x1(2),r12,-ix1(2),ir12)
     &          - cspcos_coli(-x1(2),1d0/r12,-ix1(2),-ir12)
     &          + cspcos_coli(-x1(2),k23/k13,-ix1(2),-dreal(k13))
c     &          - ettile(-x4(2),k24,-ix4(2),ik24) *
     &          - eta2s_coli(-x4(2),k24,-ix4(2),ik24,-ix4(2)*dreal(k24))
     &          *(l1(2)+cln_coli(k13,-1d0))
            
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>     &       + pi2_6
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>     &       + cspcos_coli(-x4(2),k14,-ix4(2),ik14)
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>     &       - cspcos_coli(-x1(2),r12,-ix1(2),ir12)
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>     &       - cspcos_coli(-x1(2),1d0/r12,-ix1(2),-ir12)
C>        write(testout,*) 'D0ms1ir0_check contr n ', 
C>     &       + cspcos_coli(-x1(2),k23/k13,-ix1(2),-dreal(k13))
C>        write(testout,*) 'D0ms1ir0_check contr n ', l1(2),
C>C     &       - ettile(-x4(2),k24,-ix4(2),ik24)
C>     &          - eta2s_coli(-x4(2),k24,-ix4(2),ik24,-ix4(2)*dreal(k24))
C>     &          *(l1(2)+cln_coli(k13,-1d0))
C>
C>        write(testout,*) 'D0ms1ir0_check sum ', D0ms1ir0_check
            D0ms1ir0_check = -D0ms1ir0_check/
     &          (m1*m2*mm3*m4*(k13*k24-k14*k23)*coliminfscale)
c        write(testout,*) 'D0ms1ir0_check = ', D0ms1ir0_check       

            if(cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0).gt.testacc) then
              write(testout,100)   ' D0ms1ir0_coli inconsistent b'
              write(testout,111) ' D0ms1ir0_coli: p12 = ',p12
              write(testout,111) ' D0ms1ir0_coli: p23 = ',p23
              write(testout,111) ' D0ms1ir0_coli: p34 = ',p34
              write(testout,111) ' D0ms1ir0_coli: p14 = ',p14
              write(testout,111) ' D0ms1ir0_coli: p24 = ',p24
              write(testout,111) ' D0ms1ir0_coli: p13 = ',p13
              write(testout,111) ' D0ms1ir0_coli: m12 = ',m12
              write(testout,111) ' D0ms1ir0_coli: m22 = ',m22
              write(testout,111) ' D0ms1ir0_coli: m32 = ',m32
              write(testout,111) ' D0ms1ir0_coli: m42 = ',m42
              write(testout,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
              write(testout,111) ' D0ms1ir0_check = ',D0ms1ir0_check
              write(testout,101) ' diff        = ',
     &            cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
              
              write(*,100)   ' D0ms1ir0_coli inconsistent b'
              write(*,111) ' D0ms1ir0_coli: p12 = ',p12
              write(*,111) ' D0ms1ir0_coli: p23 = ',p23
              write(*,111) ' D0ms1ir0_coli: p34 = ',p34
              write(*,111) ' D0ms1ir0_coli: p14 = ',p14
              write(*,111) ' D0ms1ir0_coli: p24 = ',p24
              write(*,111) ' D0ms1ir0_coli: p13 = ',p13
              write(*,111) ' D0ms1ir0_coli: m12 = ',m12
              write(*,111) ' D0ms1ir0_coli: m22 = ',m22
              write(*,111) ' D0ms1ir0_coli: m32 = ',m32
              write(*,111) ' D0ms1ir0_coli: m42 = ',m42
              write(*,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
              write(*,111) ' D0ms1ir0_check = ',D0ms1ir0_check
              write(*,101) ' diff        = ',
     &            cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
              stop
            endif
          endif
#endif

        else if (mm12.ne.0d0) then
c cD0ms00e
c Denner Dittmaier (3.13)

#ifdef ALLCHECK
c ppenj
          if(flag(4))then
            write(*,100)' D0ms1ir0_coli:  (3.13) m12 =/=0, m22=0',
     &          '    case with 1 zero masses:',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms1ir0_coli: p12 = ',p12
            write(*,111)' D0ms1ir0_coli: p23 = ',p23
            write(*,111)' D0ms1ir0_coli: p34 = ',p34
            write(*,111)' D0ms1ir0_coli: p14 = ',p14
            write(*,111)' D0ms1ir0_coli: p24 = ',p24
            write(*,111)' D0ms1ir0_coli: p13 = ',p13
            write(*,111)' D0ms1ir0_coli: m12 = ',m12
            write(*,111)' D0ms1ir0_coli: m22 = ',m22
            write(*,111)' D0ms1ir0_coli: m32 = ',m32
            write(*,111)' D0ms1ir0_coli: m42 = ',m42
c            write(*,111)' D0ms1ir0_coli: q12 = ',q12
            write(*,111)' D0ms1ir0_coli: q23 = ',q23
c            write(*,111)' D0ms1ir0_coli: q34 = ',q34
            write(*,111)' D0ms1ir0_coli: q14 = ',q14
            write(*,111)' D0ms1ir0_coli: q24 = ',q24
            write(*,111)' D0ms1ir0_coli: q13 = ',q13
            write(*,111)' D0ms1ir0_coli: mm12= ',mm12
            write(*,111)' D0ms1ir0_coli: mm22= ',mm22
            write(*,111)' D0ms1ir0_coli: mm32= ',mm32
            write(*,111)' D0ms1ir0_coli: mm42= ',mm42
            flag(4)=.false.
          endif
#endif

          m1  = sqrt(mm12)
          m4  = sqrt(mm42)

c     (permutation 3 <-> 2 for  m -> k )
          mm3 = m1
          mm2 = m1
          k12 = (mm12    -q13)/(m1*mm3)
          k13 = (mm12    -p12)/(m1*mm2)
          k14 = (mm12    -q14)/(m1*m4)
          k23 = (       -q23)/(mm2*mm3)
          k34 = (       -q24)/(mm2*m4)
          
          x(1) = (k12*k34-k14*k23)/k23
          x(2) = -(k13*k12-k23)/(k12*k34-k14*k23)
          ix(1) = dreal(x(1))
          ix(2) = -dreal(x(1))

c>          write(*,*) 'D0ms1ir0_coli: q ',p12,q23,p34
c>          write(*,*) 'D0ms1ir0_coli: q ',q14,q13,q24
c>          write(*,*) 'D0ms1ir0_coli: m ',m12,m4
c>          write(*,*) 'D0ms1ir0_coli: k ',k13,k12,k14
c>          write(*,*) 'D0ms1ir0_coli: k ',k23,k34,k12
c>          write(*,*) 'D0ms1ir0_coli: x ',x(1),x(2)
c>          write(*,*) 'D0ms1ir0_coli: ix4',ix(1),ix(2)
        
c       only one massless line can lead to IR singularities
c       if(k34.ne.0d0) then
          D0ms1ir0_coli = pi2_6  
     &        - cspcos_coli(-1d0/x(2),1d0/k14,ix(2),1d0)
     &        - cspcos_coli(-x(1),1d0/k14,-ix(1),1d0)
     &        + cspcos_coli(-1d0/x(2),k13/k34,ix(2),dreal(k13-k34)) 
     &        + (cln_coli(-x(2)*coliminfscale2,-ix(2)) 
     &        - cln_coli(-x(1),-ix(1)))*
     &        (cln_coli(k12/k23,dreal(k12-k23))
     &        +cln_coli(k34/k14,dreal(k34-k14)) )
c       else 
c         D0ms1ir0_coli = 
c    &               cspcos_coli(-1d0/x(2),1d0/k14,ix(2),1d0)
c    &             + cspcos_coli(-x(1),1d0/k14,-ix(1),1d0)
c    &           + (cln_coli(-x(2)*coliminfscale2,-ix(2)) - cln_coli(-x(1),-ix(1)))*
c    &             (cln_coli(k12/k23,dreal(k12-k23))+cln_coli(k13/k14,k13-k14))
c    &              - 0.5d0*(cln_coli(-x(2),-ix(2)) + cln_coli(-x(1),-ix(1))) 
c       end if 

c>      write(*,*) 'D0ms1ir0_coli contr',
c>     &             pi2_6   
c>      write(*,*) 'D0ms1ir0_coli contr',
c>     &             - cspcos_coli(-x(1),dcmplx(1d0/k14),-ix(1),1d0)
c>      write(*,*) 'D0ms1ir0_coli contr',
c>     &             - cspcos_coli(-1d0/x(2),dcmplx(1d0/k14),ix(2),1d0)
c>      write(*,*) 'D0ms1ir0_coli contr',
c>     &             + cspcos_coli(-1d0/x(2),dcmplx(k13/k34),ix(2),dreal(k13-k34)) 
c>      write(*,*) 'D0ms1ir0_coli contr',
c>     &           + (cln_coli(-x(2)*coliminfscale2,-ix(2)) - cln_coli(-x(1),-ix(1)))*
c>     &             (cln_coli(k12/k23,dreal(k12-k23))
c>     &             +cln_coli(k34/k14,dreal(k34-k14)))
C>
C>      write(*,*) 'D0ms1ir0_coli sum ',D0ms1ir0_coli     

          D0ms1ir0_coli = 
     &        D0ms1ir0_coli/(mm2*m4*mm3*m1*(k14*k23-k12*k34))

c      write(*,*) 'D0ms1ir0_coli = ',D0ms1ir0_coli     

#ifdef CHECK
          if (conscheck) then

            m1 = sqrt(mm12)
            m4 = sqrt(mm42)
            k14 = (mm12 - q14)/(m1*m4*coliminfscale)
            k24 = -q24/(m1*m4*coliminfscale)
            k12 = (mm12 - p12)/mm12
            k13 = (mm12 - q13)/mm12
            k23 = -q23/mm12
            x(1) = (k13*k24-k14*k23)/(k23)
            x(2) = (k23-k12*k13)/(k13*k24-k14*k23)
            ix(1) = dreal(x(1))
            ix(2) = -ix(1)
      
C>        write(*,*) 'D0ms1ir0_check: k  ',k12,k13,k14
C>        write(*,*) 'D0ms1ir0_check: k  ',k23,k24
C>        write(*,*) 'D0ms1ir0_check: x4 ',x(1),x(2)
C>        write(*,*) 'D0ms1ir0_check: ix4',x(1),-x(1)
C>
C>        write(*,*) 'D0ms1ir0_check: contr ',
C>     &       - cspcos_coli(-dcmplx(x(1)),dcmplx(k14),-x(1),-1d0)
C>        write(*,*) 'D0ms1ir0_check: contr ',
C>     &       - cspcos_coli(-dcmplx(x(1)),dcmplx(1d0/k14),-x(1),1d0)
C>        write(*,*) 'D0ms1ir0_check: contr ',
C>     &       + cspcos_coli(-dcmplx(x(1)),dcmplx(k24/k12),-x(1),dreal(k24-k12))
C>        write(*,*) 'D0ms1ir0_check: contr ',
C>     &       - cln_coli(-x(1),-x(1)) *
C>     &       (cln_coli(k12,-1d0) + cln_coli(k13/k23,k13-k23))
C>        write(*,*) 'D0ms1ir0_check: contr ',
C>     &       + cspcos_coli(-dcmplx(x(2)),dcmplx(k14),x(1),-1d0)
C>        write(*,*) 'D0ms1ir0_check: contr ',
C>     &       + cspcos_coli(-dcmplx(x(2)),dcmplx(1d0/k14),x(1),1d0)
C>        write(*,*) 'D0ms1ir0_check: contr ',
C>     &       - cspcos_coli(-dcmplx(x(2)),dcmplx(k24/k12),x(1),dreal(k24-k12))
C>        write(*,*) 'D0ms1ir0_check: contr ',
C>     &       cln_coli(-x(2),x(1))*
C>     &       (cln_check(k12,-1d0) + cln_coli(k13/k23,k13-k23))
C>
C>        D0ms1ir0_check = 
C>     &       - cspcos_coli(-dcmplx(x(1)),dcmplx(k14),-x(1),-1d0)
C>     &       - cspcos_coli(-dcmplx(x(1)),dcmplx(1d0/k14),-x(1),1d0)
C>     &       + cspcos_coli(-dcmplx(x(2)),dcmplx(k14),x(1),-1d0)
C>     &       + cspcos_coli(-dcmplx(x(2)),dcmplx(1d0/k14),x(1),1d0)
C>     &       + cspcos_coli(-dcmplx(x(1)),dcmplx(k24/k12),-x(1),dreal(k24-k12))
C>     &       - cspcos_coli(-dcmplx(x(2)),dcmplx(k24/k12),x(1),dreal(k24-k12))
C>     &       + (cln_coli(-x(2),x(1)) - cln_coli(-x(1),-x(1))) *
C>     &       (cln_coli(k12,-1d0) + cln_coli(k13/k23,k13-k23))
C>        write(*,*) 'D0ms1ir0_check 1  ',D0ms1ir0_check
C>        D0ms1ir0_check = D0ms1ir0_check/((q14-m12)*q23-(q13-m12)*q24)

c         write(*,*) 'D0ms1ir0_check 1=  ',D0ms1ir0_check

            D0ms1ir0_check = pi2_6 
     &          - cspcos_coli(-x(1),1d0/k14,-ix(1),1d0)
     &          - cspcos_coli(-1d0/x(2),1d0/k14,-ix(1),1d0)
     &          + cspcos_coli(-1d0/x(2),k12/k24,-ix(1),dreal(k12-k24))
     &          + (cln_coli(-x(2),ix(1)) - cln_coli(-x(1),-ix(1))) *
     &          (cln_coli(k24/k14,dreal(k24-k14)) 
     &          + cln_coli(k13/k23,dreal(k13-k23)))

c>        write(*,*) 'D0ms1ir0_check cont  ',
c>     &          pi2_6 
c>        write(*,*) 'D0ms1ir0_check cont  ',
c>     &       - cspcos_coli(-x(1),dcmplx(1d0/k14),-ix(1),1d0)
c>        write(*,*) 'D0ms1ir0_check cont  ',
c>     &       - cspcos_coli(-1d0/x(2),1d0/k14,-ix(1),1d0)
c>        write(*,*) 'D0ms1ir0_check cont  ',
c>     &       + cspcos_coli(-1d0/x(2),k12/k24,-ix(1),dreal(k12-k24))
c>        write(*,*) 'D0ms1ir0_check cont  ',
c>     &       + (cln_coli(-x(2),ix(1)) - cln_coli(-x(1),-ix(1)1)) *
c>     &       (cln_coli(k24/k14,dreal(k24-k14)) 
c>     &        + cln_coli(k13/k23,dreal(k13-k23)))
        
c        write(*,*) 'D0ms1ir0_check sum  ',D0ms1ir0_check

            D0ms1ir0_check = 
     &          D0ms1ir0_check/((q14-mm12)*q23-(q13-mm12)*q24)

c        write(*,*) 'D0ms1ir0_check=  ',D0ms1ir0_check

c        write(*,*) 'D0ms1ir0_check 1',D0ms1ir0_check


            if(cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0).gt.testacc) then
              write(testout,100) ' D0ms1ir0_coli inconsistent c'
              write(testout,111) ' D0ms1ir0_coli: p12 = ',p12
              write(testout,111) ' D0ms1ir0_coli: p23 = ',p23
              write(testout,111) ' D0ms1ir0_coli: p34 = ',p34
              write(testout,111) ' D0ms1ir0_coli: p14 = ',p14
              write(testout,111) ' D0ms1ir0_coli: p24 = ',p24
              write(testout,111) ' D0ms1ir0_coli: p13 = ',p13
              write(testout,111) ' D0ms1ir0_coli: m12 = ',m12
              write(testout,111) ' D0ms1ir0_coli: m22 = ',m22
              write(testout,111) ' D0ms1ir0_coli: m32 = ',m32
              write(testout,111) ' D0ms1ir0_coli: m42 = ',m42
              write(testout,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
              write(testout,111) ' D0ms1ir0_check = ',D0ms1ir0_check
              write(testout,101) ' diff        = ',
     &            cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
              
              write(*,100) ' D0ms1ir0_coli inconsistent c'
              write(*,111) ' D0ms1ir0_coli: p12 = ',p12
              write(*,111) ' D0ms1ir0_coli: p23 = ',p23
              write(*,111) ' D0ms1ir0_coli: p34 = ',p34
              write(*,111) ' D0ms1ir0_coli: p14 = ',p14
              write(*,111) ' D0ms1ir0_coli: p24 = ',p24
              write(*,111) ' D0ms1ir0_coli: p13 = ',p13
              write(*,111) ' D0ms1ir0_coli: m12 = ',m12
              write(*,111) ' D0ms1ir0_coli: m22 = ',m22
              write(*,111) ' D0ms1ir0_coli: m32 = ',m32
              write(*,111) ' D0ms1ir0_coli: m42 = ',m42
              write(*,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
              write(*,111) ' D0ms1ir0_check = ',D0ms1ir0_check
              write(*,101) ' diff        = ',
     &            cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
              stop
            endif
          endif
#endif

        else
          write(*,*) 'D0ms1ir0_coli with m12=0 and m32=0',
     &        ' or m22=0 and m42=0',
     &        ' not yet implemented'
          write(*,*) 'D0ms1ir0_coli: p12 = ',p12
          write(*,*) 'D0ms1ir0_coli: p23 = ',p23
          write(*,*) 'D0ms1ir0_coli: p34 = ',p34
          write(*,*) 'D0ms1ir0_coli: p14 = ',p14
          write(*,*) 'D0ms1ir0_coli: p24 = ',p24
          write(*,*) 'D0ms1ir0_coli: p13 = ',p13
          write(*,*) 'D0ms1ir0_coli: m12 = ',m12
          write(*,*) 'D0ms1ir0_coli: m22 = ',m22
          write(*,*) 'D0ms1ir0_coli: m32 = ',m32
          write(*,*) 'D0ms1ir0_coli: m42 = ',m42
          stop
        end if

      else 
c cD0mse
c Denner Dittmaier (3.14)
* 0 soft singularities and 0 zero masses
*
*                  m22             
*        p12 ----------------  p23    
*                |  2   |         
*            m12 |1    3| m32=m42       
*                |   4  |         
*        p14 ----------------  0
*                  m42             
*

        if(m12.ne.cd0.and.m22.ne.cd0)then

#ifdef CHECK
          if(m32.ne.m42.or.p34.ne.0d0.or.
     &        p14.eq.m12.or.p23.eq.m22)
     &        then
            write(*,100)' D0ms1ir0_coli:  (3.14)  m12,m22 =/=0',
     &          '    case with 0 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms1ir0_coli: p12 = ',p12
            write(*,111)' D0ms1ir0_coli: p23 = ',p23
            write(*,111)' D0ms1ir0_coli: p34 = ',p34
            write(*,111)' D0ms1ir0_coli: p14 = ',p14
            write(*,111)' D0ms1ir0_coli: p24 = ',p24
            write(*,111)' D0ms1ir0_coli: p13 = ',p13
            write(*,111)' D0ms1ir0_coli: m12 = ',m12
            write(*,111)' D0ms1ir0_coli: m22 = ',m22
            write(*,111)' D0ms1ir0_coli: m32 = ',m32
            write(*,111)' D0ms1ir0_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef ALLCHECK
c ppenj
          if(flag(3))then
            write(*,100)' D0ms1ir0_coli:  (3.14) m12,m22 =/=0',
     &          '    case with 0 zero masses:',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms1ir0_coli: p12 = ',p12
            write(*,111)' D0ms1ir0_coli: p23 = ',p23
            write(*,111)' D0ms1ir0_coli: p34 = ',p34
            write(*,111)' D0ms1ir0_coli: p14 = ',p14
            write(*,111)' D0ms1ir0_coli: p24 = ',p24
            write(*,111)' D0ms1ir0_coli: p13 = ',p13
            write(*,111)' D0ms1ir0_coli: m12 = ',m12
            write(*,111)' D0ms1ir0_coli: m22 = ',m22
            write(*,111)' D0ms1ir0_coli: m32 = ',m32
            write(*,111)' D0ms1ir0_coli: m42 = ',m42
            flag(3)=.false.
          endif
#endif

c cD0mse
c     Denner Dittmaier (3.14)
          q12=p13
          q23=p23
          q34=p24
          q14=p14
          q13=p12
          q24=p34
          mm12=m12
          mm22=m32*coliminfscale2
          mm32=m22
          mm42=m42*coliminfscale2
          
          l12 = (mm12    -q12)
          l13 = (mm12+mm32-q13)
          l14 = (mm12    -q14)
          l23 = (    mm32-q23)
          l24 = (       -q24)
          l34 = (mm32    -q34)
          
          r12 = l12/mm22
          r21 = l12/mm12  
          r14 = l14/mm42
          r41 = l14/mm12
          r23 = l23/mm32
          r32 = l23/mm22
          r34 = l34/mm42
          r43 = l34/mm32
          if (l13.ne.0d0) then
            r13 = l13/(2d0*mm32)*(1d0+sqrt(1d0-4d0*mm12*mm32/l13**2))
          else
            r13 = sqrt(mm12/mm32)*dcmplx(0d0,1d0)
          end if
          r31 = r13*mm32/mm12
c     r13 = 1d0/r13
          r24=1d0
          r42=mm42/mm22
          a   =  mm42*(l34/r42-l23 + (l12-l14/r42)*r31)
          b   =  l12*l34-l14*l23
          c   =  mm12*(l34*r24-l23 + (l12-l14*r24)/r13)
          d   =  mm12*(l23-r31*l12 - (l34-r31*l14)*r24)
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
c     if(k13.le.-2d0) then
          ir13 = sign(1d1,1d0-abs(r13*r13*mm32/mm12))
c     else
c     ir13 = 0d0
c     end if

c>          write(*,*) 'D0ms1ir0_coli q ',q12,q23,q34
c>          write(*,*) 'D0ms1ir0_coli q ',q14,q13,q24
c>          write(*,*) 'D0ms1ir0_coli m ',mm12,mm22,mm32,mm42
c>          write(*,*) 'D0ms1ir0_coli k ',l12,l13,l14
c>          write(*,*) 'D0ms1ir0_coli k ',l23,l24,l34
c>          write(*,*) 'D0ms1ir0_coli r ',r12,r13,r14
c>          write(*,*) 'D0ms1ir0_coli r ',r21,r31,r41
c>          write(*,*) 'D0ms1ir0_coli r ',r23,r24,r34
c>          write(*,*) 'D0ms1ir0_coli r ',r32,r42,r43
c>          write(*,*) 'D0ms1ir0_coli ir',ir12,ir13,ir14
c>          write(*,*) 'D0ms1ir0_coli ir',ir23,ir34
c>          write(*,*) 'D0ms1ir0_coli a ',a,b,c,d

          z(1) = y(1)*r31
          z(2) = y(2)*r31

          rd = dreal(d)
c     result should be independent of imaginary parts of y and z!!!!
c          if(dreal(y(1)).gt.0d0) then
            iy(1) = sign(1d0,-rd)
c          else
c            iy(1) = 0d0
c          end if
c          if(dreal(y(2)).gt.0d0) then
            iy(2) = sign(1d0,rd)
c          else
c            iy(2) = 0d0
c          end if
          iz(1)=sign(1d0,iy(1)*dreal(r13))
          iz(2)=sign(1d0,iy(2)*dreal(r13))

c>          write(*,*) 'D0ms1ir0_coli y ',y(1),y(2)
c>          write(*,*) 'D0ms1ir0_coli iy',iy(1),iy(2)
c>          write(*,*) 'D0ms1ir0_coli z',z(1),z(2)
c>          write(*,*) 'D0ms1ir0_coli iz',iz(1),iz(2)

          D0ms1ir0_coli = dcmplx(0d0)
          D0ms1ir0_coli = D0ms1ir0_coli + (
     &        cspcos_coli(-y(1),cd0,-iy(1),-ir12)
     &        +cspcos_coli(-y(1),r21,-iy(1),ir12)
     &        -cspcos_coli(-y(1),cd0,-iy(1),-ir14)
     &        -cspcos_coli(-y(1),r41,-iy(1),ir14)
     &        +cspcos_coli(-z(1),cd0,-iz(1),-ir34)
     &        +cspcos_coli(-z(1),r43,-iz(1),ir34)
     &        -cspcos_coli(-z(1),cd0,-iz(1),-ir23)
     &        -cspcos_coli(-z(1),r23,-iz(1),ir23)
     &        )          

c>          write(*,*) 'D0ms1ir0_coli args yz',
c>     &        -y(1),-iy(1),-z(1),-iz(1)
c>          write(*,*) 'D0ms1ir0_coli args r',
c>     &        r21,r41,
c>     &        r43,r23
c>          write(*,*) 'D0ms1ir0_coli args ',
c>     &        1d0+y(1)*r21,
c>     &        1d0+y(1)*r41,
c>     &        1d0+z(1)*r43,
c>     &        1d0+z(1)*r23


c>          write(*,*) 'D0ms1ir0_coli arg12 ',
c>     &       -y(1),r21,-iy(1),ir12 
c>          write(*,*) 'D0ms1ir0_coli arg14 ',
c>     &       -y(1),r41,-iy(1),ir14
c>          write(*,*) 'D0ms1ir0_coli csp ',
c>     &        cspcos_coli(-y(1),cd0,-iy(1),-ir12),
c>     &        +cspcos_coli(-y(1),r21,-iy(1),ir12),
c>     &        -cspcos_coli(-y(1),cd0,-iy(1),-ir14),
c>     &        -cspcos_coli(-y(1),r41,-iy(1),ir14),
c>     &        +cspcos_coli(-z(1),cd0,-iz(1),-ir34),
c>     &        +cspcos_coli(-z(1),r43,-iz(1),ir34),
c>     &        -cspcos_coli(-z(1),cd0,-iz(1),-ir23),
c>     &        -cspcos_coli(-z(1),r23,-iz(1),ir23)
c>          write(*,*) 'D0ms1ir0_coli cspsum ',D0ms1ir0_coli


c     *** improves numerical stability in extreme cases

          eta=eta2s_coli(-y(1),-1/z(1),-iy(1),iz(1),-ir13)

          if(eta.ne.0d0)then
            h14 = mm12+l14*y(1)
            h12 = mm12+l12*y(1)
            h34 = mm32+l34*z(1)
            h23 = mm32+l23*z(1)
            
            rt14 = abs(h14)/max(abs(mm12),abs(l14*y(1)))
            rt12 = abs(h12)/max(abs(mm12),abs(l12*y(1)))
            rt34 = abs(h34)/max(abs(mm32),abs(l34*z(1)))
            rt23 = abs(h23)/max(abs(mm32),abs(l23*z(1)))
            
            if (min(rt12,rt14).gt.min(rt34,rt23)) then
              if(dimag(h14/h12).ne.0d0)then
                lqbard= cln_coli(h14/h12,0d0)
              else
                lqbard=cln_coli(h14,-1d0)-cln_coli(h12,-1d0)
              endif
            else
              if(dimag(h34/h23).ne.0d0)then
                lqbard= cln_coli(h34/h23,0d0)
              else
                lqbard=cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
              endif
            endif
                        
c>            write(*,*) 'D0ms1ir0_coli lqbard',h14/h12,h34/h23              
c>            write(*,*) 'D0ms1ir0_coli lqbard',lqbard,
c>     &          cln_coli(h14,-1d0)-cln_coli(h12,-1d0),
c>     &          cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
            
            D0ms1ir0_coli = D0ms1ir0_coli + eta*lqbard

c>            write(*,*) 'D0ms1ir0_coli etacon ',lqbard*eta
c>            write(*,*) 'D0ms1ir0_coli etasum ',D0ms1ir0_coli

          endif          
          
          D0ms1ir0_coli = D0ms1ir0_coli - (
     &        cspcos_coli(-y(2),1d0/r12,-iy(2),-ir12)     
     &        -.5d0*(cln_coli(-y(2),-iy(2))+cln_coli(r21,ir12))**2
     &        -cspcos_coli(cd0,1d0/r21,iy(2),-ir12)
     &        -cspcos_coli(-y(2),1d0/(r14),-iy(2),-ir14)
     &        +.5d0*(cln_coli(-y(2),-iy(2))+cln_coli(r41,ir14))**2
     &        +cspcos_coli(cd0,1d0/r41,iy(2),-ir14)
     &        +cspcos_coli(-z(2),1d0/r34,-iz(2),-ir34)
     &        -.5d0*(cln_coli(-z(2),-iz(2))+cln_coli(r43,ir34))**2
     &        -cspcos_coli(cd0,1d0/r43,iz(2),-ir34)
     &        -cspcos_coli(-z(2),1d0/(r32),-iz(2),-ir23)
     &        +.5d0*(cln_coli(-z(2),-iz(2))+cln_coli(r23,ir23))**2
     &        +cspcos_coli(cd0,1d0/r23,iz(2),-ir23)
     &        ) 
          
c>          write(*,*) 'D0ms1ir0_coli args yz',
c>     &        -y(2),-iy(2),-z(2),-iz(2)
c>          write(*,*) 'D0ms1ir0_coli args r',
c>     &        1d0/(r12),r21,1d0/(r14),r41,
c>     &        1d0/(r34),r43,1d0/(r32),r23
c>          write(*,*) 'D0ms1ir0_coli args ',
c>     &        1d0+y(1)/(r12),
c>     &        1d0+y(1)*r21,
c>     &        1d0+y(1)/(r14),
c>     &        1d0+y(1)*r41,
c>     &        1d0+z(1)/(r34),
c>     &        1d0+z(1)*r43,
c>     &        1d0+z(1)/(r32),
c>     &        1d0+z(1)*r23
c>          write(*,*) 'D0ms1ir0_coli arg25 ',
c>     &        -z(2),1d0/(r34),-iz(2),-ir34
c>          write(*,*) 'D0ms1ir0_coli arg26 ',
c>     &        cd0,1d0/(r43),iz(2),-ir34
c>          write(*,*) 'D0ms1ir0_coli arg27 ',
c>     &        -z(2),1d0/(r32),-iz(2),-ir23
c>          write(*,*) 'D0ms1ir0_coli arg28 ',
c>     &        cd0,1d0/(r23),iz(2),-ir23

c>          write(*,*) 'D0ms1ir0_coli csp ',
c>     &        cspcos_coli(-y(2),1d0/(r12),-iy(2),-ir12),     
c>     &        -.5d0*(cln_coli(-y(2),-iy(2))+cln_coli(r21,ir12))**2
c>     &        -cspcos_coli(cd0,1d0/(r21),iy(2),-ir12),
c>     &        -cspcos_coli(-y(2),1d0/(r14),-iy(2),-ir14),
c>     &        +.5d0*(cln_coli(-y(2),-iy(2))+cln_coli(r41,ir14))**2
c>     &        +cspcos_coli(cd0,1d0/(r41),iy(2),-ir14),
c>     &        +cspcos_coli(-z(2),1d0/(r34),-iz(2),-ir34),
c>     &        -.5d0*(cln_coli(-z(2),-iz(2))+cln_coli(r43,ir34))**2
c>     &        -cspcos_coli(cd0,1d0/(r43),iz(2),-ir34),
c>     &        -cspcos_coli(-z(2),1d0/(r32),-iz(2),-ir23),
c>     &        +.5d0*(cln_coli(-z(2),-iz(2))+cln_coli(r23,ir23))**2
c>     &        +cspcos_coli(cd0,1d0/(r23),iz(2),-ir23)
c>
c>          write(*,*) 'D0ms1ir0_coli cspsum ',D0ms1ir0_coli


c     *** improves numerical stability in extreme cases

          eta=eta2s_coli(-y(2),-1/z(2),-iy(2),iz(2),-ir13)
         
          if(eta.ne.0d0)then
            h14 = mm42*y(2)*y(2)+l14*y(2)
            h12 = mm22*y(2)*y(2)+l12*y(2)
            h34 = mm42*z(2)*z(2)+l34*z(2)
            h23 = mm22*z(2)*z(2)+l23*z(2)
            
            rt14 = abs(h14)/max(abs(mm42*y(2)*y(2)),abs(l14*y(2)))
            rt12 = abs(h12)/max(abs(mm22*y(2)*y(2)),abs(l12*y(2)))
            rt34 = abs(h34)/max(abs(mm42*z(2)*z(2)),abs(l34*z(2)))
            rt23 = abs(h23)/max(abs(mm22*z(2)*z(2)),abs(l23*z(2)))
            
            if (min(rt12,rt14).gt.min(rt34,rt23)) then
              if(dimag(h14/h12).ne.0d0)then
                lqbard= cln_coli(h14/h12,0d0)
              else
                lqbard=cln_coli(h14,-1d0)-cln_coli(h12,-1d0)
              endif
            else
              if(dimag(h34/h23).ne.0d0)then
                lqbard= cln_coli(h34/h23,0d0)
              else
                lqbard=cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
              endif
            endif
            
c>            write(*,*) 'D0ms1ir0_coli lqbard',h14/h12,h34/h23              
c>            write(*,*) 'D0ms1ir0_coli lqbard',
c>     &          cln_coli(h14,-1d0)-cln_coli(h12,-1d0),
c>     &          cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
            
            D0ms1ir0_coli = D0ms1ir0_coli 
     &          - eta*lqbard

c>            write(*,*) 'D0ms1ir0_coli etacon ', -lqbard*eta
c>            write(*,*) 'D0ms1ir0_coli etasum ',D0ms1ir0_coli

          endif
          
c          write(*,*) 'D0ms1ir0_coli sum ',D0ms1ir0_coli

          D0ms1ir0_coli = D0ms1ir0_coli/b

c          write(*,*) 'D0ms1ir0_coli =   ',D0ms1ir0_coli    
          
#ifdef CHECK
          if(conscheck) then

            do j=1,3
              if(j.eq.1) then
                if(dreal(y(1)).gt.0d0) then
                  iy(1) = -sign(1d0,-rd)
                else
                  iy(1) = 0d0
                end if
                if(dreal(y(2)).gt.0d0) then
                  iy(2) = -sign(1d0,rd)
                else
                  iy(2) = 0d0
                end if
              elseif(j.eq.2) then
                if(dreal(y(1)).gt.0d0) then
                  iy(1) = -sign(1d0,-rd)
                else
                  iy(1) = 0d0
                end if
                if(dreal(y(2)).gt.0d0) then
                  iy(2) = sign(1d0,rd)
                else
                  iy(2) = 0d0
                end if
              elseif(j.eq.3) then
c                if(dreal(y(1)).gt.0d0) then
                  iy(1) = sign(1d0,-rd)
c                else
c                  iy(1) = 0d0
c                end if
c                if(dreal(y(2)).gt.0d0) then
                  iy(2) = -sign(1d0,rd)
c                else
c                  iy(2) = 0d0
c                end if
              endif
              iz(1)=sign(1d0,iy(1)*dreal(r13))
              iz(2)=sign(1d0,iy(2)*dreal(r13))


c>          write(*,*) 'D0ms1ir0_check y ',y(1),y(2)
c>          write(*,*) 'D0ms1ir0_check iy',iy(1),iy(2)
c>          write(*,*) 'D0ms1ir0_check z',z(1),z(2)
c>          write(*,*) 'D0ms1ir0_check iz',iz(1),iz(2)

              D0ms1ir0_check = dcmplx(0d0)
              D0ms1ir0_check = D0ms1ir0_check + (
     &            cspcos_coli(-y(1),cd0,-iy(1),-ir12)
     &            +cspcos_coli(-y(1),r21,-iy(1),ir12)
     &            -cspcos_coli(-y(1),cd0,-iy(1),-ir14)
     &            -cspcos_coli(-y(1),r41,-iy(1),ir14)
     &            +cspcos_coli(-z(1),cd0,-iz(1),-ir34)
     &            +cspcos_coli(-z(1),r43,-iz(1),ir34)
     &            -cspcos_coli(-z(1),cd0,-iz(1),-ir23)
     &            -cspcos_coli(-z(1),r23,-iz(1),ir23)
     &            )          

c>          write(*,*) 'D0ms1ir0_coli args yz',
c>     &        -y(1),-iy(1),-z(1),-iz(1)
c>          write(*,*) 'D0ms1ir0_coli args r',
c>     &        r21,r41,r43,r23
c>          write(*,*) 'D0ms1ir0_coli args ir',
c>     &        ir12,ir14,ir34,ir23
c>          write(*,*) 'D0ms1ir0_coli args ',
c>     &        1d0+y(1)*r21,
c>     &        1d0+y(1)*r41,
c>     &        1d0+z(1)*r43,
c>     &        1d0+z(1)*r23

c>              write(*,*) 'D0ms1ir0_check arg12 ',
c>     &            -y(1),r21,-iy(1),ir12 
c>              write(*,*) 'D0ms1ir0_check arg14 ',
c>     &            -y(1),r41,-iy(1),ir14
c>              write(*,*) 'D0ms1ir0_check csp ',
c>     &            cspcos_coli(-y(1),cd0,-iy(1),-ir12),
c>     &            +cspcos_coli(-y(1),r21,-iy(1),ir12),
c>     &            -cspcos_coli(-y(1),cd0,-iy(1),-ir14),
c>     &            -cspcos_coli(-y(1),r41,-iy(1),ir14),
c>     &            +cspcos_coli(-z(1),cd0,-iz(1),-ir34),
c>     &            +cspcos_coli(-z(1),r43,-iz(1),ir34),
c>     &            -cspcos_coli(-z(1),cd0,-iz(1),-ir23),
c>     &            -cspcos_coli(-z(1),r23,-iz(1),ir23)
c>              write(*,*) 'D0ms1ir0_check cspsum ',D0ms1ir0_check

c     *** improves numerical stability in extreme cases

              eta=eta2s_coli(-y(1),-1/z(1),-iy(1),iz(1),-ir13)

              if(eta.ne.0d0)then
                h14 = mm12+l14*y(1)
                h12 = mm12+l12*y(1)
                h34 = mm32+l34*z(1)
                h23 = mm32+l23*z(1)
                
                rt14 = abs(h14)/max(abs(mm12),abs(l14*y(1)))
                rt12 = abs(h12)/max(abs(mm12),abs(l12*y(1)))
                rt34 = abs(h34)/max(abs(mm32),abs(l34*z(1)))
                rt23 = abs(h23)/max(abs(mm32),abs(l23*z(1)))
                
                if (min(rt12,rt14).gt.min(rt34,rt23)) then
                  if(dimag(h14/h12).ne.0d0)then
                    lqbard= cln_coli(h14/h12,0d0)
                  else
                    lqbard=cln_coli(h14,-1d0)-cln_coli(h12,-1d0)
                  endif
                else
                  if(dimag(h34/h23).ne.0d0)then
                    lqbard= cln_coli(h34/h23,0d0)
                  else
                    lqbard=cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
                  endif
                endif
                        
c>            write(*,*) 'D0ms1ir0_check lqbard',h14/h12,h34/h23              
c>            write(*,*) 'D0ms1ir0_check lqbard',
c>     &          cln_coli(h14,-1d0)-cln_coli(h12,-1d0),
c>     &          cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
            
                D0ms1ir0_check = D0ms1ir0_check + eta*lqbard

c>                write(*,*) 'D0ms1ir0_check etacon ',lqbard*eta
c>                write(*,*) 'D0ms1ir0_check etasum ',D0ms1ir0_check

              endif          
              
              D0ms1ir0_check = D0ms1ir0_check - (
     &            cspcos_coli(-y(2),1d0/r12,-iy(2),-ir12)     
     &            -.5d0*(cln_coli(-y(2),-iy(2))+cln_coli(r21,ir12))**2
     &            -cspcos_coli(cd0,1d0/r21,iy(2),-ir12)
     &            -cspcos_coli(-y(2),1d0/(r14),-iy(2),-ir14)
     &            +.5d0*(cln_coli(-y(2),-iy(2))+cln_coli(r41,ir14))**2
     &            +cspcos_coli(cd0,1d0/r41,iy(2),-ir14)
     &            +cspcos_coli(-z(2),1d0/r34,-iz(2),-ir34)
     &            -.5d0*(cln_coli(-z(2),-iz(2))+cln_coli(r43,ir34))**2
     &            -cspcos_coli(cd0,1d0/r43,iz(2),-ir34)
     &            -cspcos_coli(-z(2),1d0/(r32),-iz(2),-ir23)
     &            +.5d0*(cln_coli(-z(2),-iz(2))+cln_coli(r23,ir23))**2
     &            +cspcos_coli(cd0,1d0/r23,iz(2),-ir23)
     &            ) 
              
c>          write(*,*) 'D0ms1ir0_check args yz',
c>     &        -y(2),-iy(2),-z(2),-iz(2)
c>          write(*,*) 'D0ms1ir0_check args r',
c>     &        1d0/(r12),r21,1d0/(r14),r41,
c>     &        1d0/(r34),r43,1d0/(r32),r23
c>          write(*,*) 'D0ms1ir0_check args ',
c>     &        1d0+y(1)*r21,
c>     &        1d0+y(1)*r41,
c>     &        1d0+z(1)*r43,
c>     &        1d0+z(1)*r23
c>              write(*,*) 'D0ms1ir0_check arg25 ',
c>     &            -z(2),1d0/(r34),-iz(2),-ir34
c>              write(*,*) 'D0ms1ir0_check arg26 ',
c>     &            cd0,1d0/(r43),iz(2),-ir34
c>              write(*,*) 'D0ms1ir0_check arg27 ',
c>     &            -z(2),1d0/(r32),-iz(2),-ir23
c>              write(*,*) 'D0ms1ir0_check arg28 ',
c>     &            cd0,1d0/(r23),iz(2),-ir23

c>              write(*,*) 'D0ms1ir0_check csp ',
c>     &            cspcos_coli(-y(2),1d0/(r12),-iy(2),-ir12),     
c>     &            -.5d0*(cln_coli(-y(2),-iy(2))+cln_coli(r21,ir12))**2
c>     &            -cspcos_coli(cd0,1d0/(r21),iy(2),-ir12),
c>     &            -cspcos_coli(-y(2),1d0/(r14),-iy(2),-ir14),
c>     &            +.5d0*(cln_coli(-y(2),-iy(2))+cln_coli(r41,ir14))**2
c>     &            +cspcos_coli(cd0,1d0/(r41),iy(2),-ir14),
c>     &            +cspcos_coli(-z(2),1d0/(r34),-iz(2),-ir34),
c>     &            -.5d0*(cln_coli(-z(2),-iz(2))+cln_coli(r43,ir34))**2
c>     &            -cspcos_coli(cd0,1d0/(r43),iz(2),-ir34),
c>     &            -cspcos_coli(-z(2),1d0/(r32),-iz(2),-ir23),
c>     &            +.5d0*(cln_coli(-z(2),-iz(2))+cln_coli(r23,ir23))**2
c>     &            +cspcos_coli(cd0,1d0/(r23),iz(2),-ir23)
c>
c>              write(*,*) 'D0ms1ir0_check cspsum ',D0ms1ir0_check


c     *** improves numerical stability in extreme cases

              eta=eta2s_coli(-y(2),-1/z(2),-iy(2),iz(2),-ir13)
              
              if(eta.ne.0d0)then
                h14 = mm42*y(2)*y(2)+l14*y(2)
                h12 = mm22*y(2)*y(2)+l12*y(2)
                h34 = mm42*z(2)*z(2)+l34*z(2)
                h23 = mm22*z(2)*z(2)+l23*z(2)
                
                rt14 = abs(h14)/max(abs(mm42*y(2)*y(2)),abs(l14*y(2)))
                rt12 = abs(h12)/max(abs(mm22*y(2)*y(2)),abs(l12*y(2)))
                rt34 = abs(h34)/max(abs(mm42*z(2)*z(2)),abs(l34*z(2)))
                rt23 = abs(h23)/max(abs(mm22*z(2)*z(2)),abs(l23*z(2)))
                
                if (min(rt12,rt14).gt.min(rt34,rt23)) then
                  if(dimag(h14/h12).ne.0d0)then
                    lqbard= cln_coli(h14/h12,0d0)
                  else
                    lqbard=cln_coli(h14,-1d0)-cln_coli(h12,-1d0)
                  endif
                else
                  if(dimag(h34/h23).ne.0d0)then
                    lqbard= cln_coli(h34/h23,0d0)
                  else
                    lqbard=cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
                  endif
                endif
                
c>            write(*,*) 'D0ms1ir0_coli lqbard',h14/h12,h34/h23              
c>            write(*,*) 'D0ms1ir0_coli lqbard',
c>     &          cln_coli(h14,-1d0)-cln_coli(h12,-1d0),
c>     &          cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
            
                D0ms1ir0_check = D0ms1ir0_check 
     &              - eta*lqbard

c                write(*,*) 'D0ms1ir0_check etacon ', -lqbard*eta
c                write(*,*) 'D0ms1ir0_check etasum ',D0ms1ir0_check

              endif
              
c              write(*,*) 'D0ms1ir0_check sum ',D0ms1ir0_check

              D0ms1ir0_check = D0ms1ir0_check/b

c              write(*,*) 'D0ms1ir0_check =   ',D0ms1ir0_check    

              if(cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0).gt.testacc)then
                write(testout,*) ' D0ms1ir0_coli inconsistent',i
                write(testout,111) ' D0ms1ir0_coli: p12 = ',p12
                write(testout,111) ' D0ms1ir0_coli: p23 = ',p23
                write(testout,111) ' D0ms1ir0_coli: p34 = ',p34
                write(testout,111) ' D0ms1ir0_coli: p14 = ',p14
                write(testout,111) ' D0ms1ir0_coli: p24 = ',p24
                write(testout,111) ' D0ms1ir0_coli: p13 = ',p13
                write(testout,111) ' D0ms1ir0_coli: m12 = ',m12
                write(testout,111) ' D0ms1ir0_coli: m22 = ',m22
                write(testout,111) ' D0ms1ir0_coli: m32 = ',m32
                write(testout,111) ' D0ms1ir0_coli: m42 = ',m42
                write(testout,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
                write(testout,111) ' D0ms1ir0_check = ',D0ms1ir0_check
                write(testout,101) ' diff        = ',
     &              cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
                
                write(*,*) ' D0ms1ir0_coli inconsistent',i
                write(*,111) ' D0ms1ir0_coli: p12 = ',p12
                write(*,111) ' D0ms1ir0_coli: p23 = ',p23
                write(*,111) ' D0ms1ir0_coli: p34 = ',p34
                write(*,111) ' D0ms1ir0_coli: p14 = ',p14
                write(*,111) ' D0ms1ir0_coli: p24 = ',p24
                write(*,111) ' D0ms1ir0_coli: p13 = ',p13
                write(*,111) ' D0ms1ir0_coli: m12 = ',m12
                write(*,111) ' D0ms1ir0_coli: m22 = ',m22
                write(*,111) ' D0ms1ir0_coli: m32 = ',m32
                write(*,111) ' D0ms1ir0_coli: m42 = ',m42
                write(*,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
                write(*,111) ' D0ms1ir0_check = ',D0ms1ir0_check
                write(*,101) ' diff        = ',
     &              cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
c                stop
              endif
            enddo

            m1 = sqrt(mm12) 
            m2 = sqrt(mm22) 
            m3 = sqrt(mm32)
            m4 = sqrt(mm42) 
            l12 = (mm12    -q12)
            l13 = (mm12+mm32-q13)
            l14 = (mm12    -q14)
            l23 = (    mm32-q23)
            l24 = (       -q24)
            l34 = (mm32    -q34)
            k12 = (mm12    -q12)/(m1*m2)
            k13 = (mm12+mm32-q13)/(m1*m3)
            k14 = (mm12    -q14)/(m1*m4)
            k23 = (   +mm32-q23)/(m2*m3)
            k24 = (       -q24)/(m2*m4)
            k34 = (mm32    -q34)/(m3*m4)
            
            r12 = k12
            r14 = k14
            r23 = k23
            r34 = k34
            if (k13.ne.0d0) then
              r13 = k13/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k13**2)))
            else
              r13 = dcmplx(0d0,1d0)
            end if
c     r13 = 1d0/r13
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
c           if(k13.le.-2d0) then
              ir13 = sign(1d1,1d0-abs(r13))
c           else
c             ir13 = 0d0
c           end if

c>            write(*,*) 'D0ms1ir0_check q ',q12,q23,q34
c>            write(*,*) 'D0ms1ir0_check q ',q14,q13,q24
c>            write(*,*) 'D0ms1ir0_check m ',mm12,mm22,mm32,mm42
c>            write(*,*) 'D0ms1ir0_check k ',k12,k13,k14
c>            write(*,*) 'D0ms1ir0_check k ',k23,k24,k34
c>            write(*,*) 'D0ms1ir0_check r ',r12,r13,r14
c>            write(*,*) 'D0ms1ir0_check r ',r23,r24,r34
c>            write(*,*) 'D0ms1ir0_check r ',r12*m1/m2,r13*m1/m3,r14*m1/m4
c>            write(*,*) 'D0ms1ir0_check r ',r12*m2/m1,r13*m3/m1,r14*m4/m1
c>            write(*,*) 'D0ms1ir0_check r ',r23*m2/m3,r24*m2/m4,r34*m3/m4
c>            write(*,*) 'D0ms1ir0_check r ',r23*m3/m2,r24*m4/m2,r34*m4/m3
c>            write(*,*) 'D0ms1ir0_check ir',ir12,ir13,ir14
c>            write(*,*) 'D0ms1ir0_check ir',ir23,ir34
c>            write(*,*) 'D0ms1ir0_check a ',a,b,c,d
c>            write(*,*) 'D0ms1ir0_check a ',a*mm42*m2*m3,b*m1*m2*m3*m4,
c>     &          c*mm12*m2*m3,d*m12*m2*m3
c>            write(*,*) 'D0ms1ir0_check x ',x(1),x(2)
c>            write(*,*) 'D0ms1ir0_check ix',ix(1),ix(2)

            y(1) = x(1)*m1/m4
            y(2) = x(2)*m1/m4
            z(1) = x(1)*m3/m4*r13
            z(2) = x(2)*m3/m4*r13

c     result should be independent of imaginary parts of y and z!!!!
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
              iz(1)=sign(1d0,iy(1)*dreal(r13))
              iz(2)=sign(1d0,iy(2)*dreal(r13))

c>            write(*,*) 'D0ms1ir0_check y ',y(1),y(2)
c>            write(*,*) 'D0ms1ir0_check iy',iy(1),iy(2)
c>            write(*,*) 'D0ms1ir0_check z',z(1),z(2)
c>            write(*,*) 'D0ms1ir0_check iz',iz(1),iz(2)

            D0ms1ir0_check = dcmplx(0d0)
            D0ms1ir0_check = D0ms1ir0_check + (
     &          cspcos_coli(-y(1),cd0,-iy(1),-ir12)
     &          +cspcos_coli(-y(1),m2/m1*r12,-iy(1),ir12)
     &          -cspcos_coli(-y(1),cd0,-iy(1),-ir14)
     &          -cspcos_coli(-y(1),m4/m1*r14,-iy(1),ir14)
     &          +cspcos_coli(-z(1),cd0,-iz(1),-ir34)
     &          +cspcos_coli(-z(1),m4/m3*r34,-iz(1),ir34)
     &          -cspcos_coli(-z(1),cd0,-iz(1),-ir23)
     &          -cspcos_coli(-z(1),m2/m3*r23,-iz(1),ir23)
     &          )          

c>            write(*,*) 'D0ms1ir0_check args yz',
c>     &          -y(1),-iy(1),-z(1),-iz(1)
c>            write(*,*) 'D0ms1ir0_check args r',
c>     &          m2/(m1*r12),m2/m1*r12,m4/(m1*r14),m4/m1*r14,
c>     &          m4/(m3*r34),m4/m3*r34,m2/(m3*r23),m2/m3*r23
c>            write(*,*) 'D0ms1ir0_check args ',
c>     &          1d0+y(1)*m2/(m1*r12),
c>     &          1d0+y(1)*m2/m1*r12,
c>     &          1d0+y(1)*m4/(m1*r14),
c>     &          1d0+y(1)*m4/m1*r14,
c>     &          1d0+z(1)*m4/(m3*r34),
c>     &          1d0+z(1)*m4/m3*r34,
c>     &          1d0+z(1)*m2/(m3*r23),
c>     &          1d0+z(1)*m2/m3*r23
c>            write(*,*) 'D0ms1ir0_check csp ',
c>     &          cspcos_coli(-y(1),cd0,-iy(1),-ir12),
c>     &          +cspcos_coli(-y(1),m2/m1*r12,-iy(1),ir12),
c>     &          -cspcos_coli(-y(1),cd0,-iy(1),-ir14),
c>     &          -cspcos_coli(-y(1),m4/m1*r14,-iy(1),ir14),
c>     &          +cspcos_coli(-z(1),cd0,-iz(1),-ir34),
c>     &          +cspcos_coli(-z(1),m4/m3*r34,-iz(1),ir34),
c>     &          -cspcos_coli(-z(1),cd0,-iz(1),-ir23),
c>     &          -cspcos_coli(-z(1),m2/m3*r23,-iz(1),ir23)
c>
c>            write(*,*) 'D0ms1ir0_check cspsum ',D0ms1ir0_check


c     *** improves numerical stability in extreme cases

            eta=eta2s_coli(-y(1),-1/z(1),-iy(1),iz(1),-ir13)
            
            if(eta.ne.0d0)then
              h14 = mm12+l14*y(1)
              h12 = mm12+l12*y(1)
              h34 = mm32+l34*z(1)
              h23 = mm32+l23*z(1)
              
              rt14 = abs(h14)/max(abs(mm12),abs(l14*y(1)))
              rt12 = abs(h12)/max(abs(mm12),abs(l12*y(1)))
              rt34 = abs(h34)/max(abs(mm32),abs(l34*z(1)))
              rt23 = abs(h23)/max(abs(mm32),abs(l23*z(1)))
              
              if (min(rt12,rt14).gt.min(rt34,rt23)) then
                if(dimag(h14/h12).ne.0d0)then
                  lqbard= cln_coli(h14/h12,0d0)
                else
                  lqbard=cln_coli(h14,-1d0)-cln_coli(h12,-1d0)
                endif
              else
                if(dimag(h34/h23).ne.0d0)then
                  lqbard= cln_coli(h34/h23,0d0)
                else
                  lqbard=cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
                endif
              endif
              
c>              write(*,*) 'D0ms1ir0_check lqbard',h14/h12,h34/h23
c>              write(*,*) 'D0ms1ir0_check lqbard',
c>     &            cln_coli(h14,-1d0)-cln_coli(h12,-1d0),
c>     &            cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
              
              D0ms1ir0_check = D0ms1ir0_check 
     &            +eta*lqbard

c              write(*,*) 'D0ms1ir0_check etacon ', lqbard*eta
c              write(*,*) 'D0ms1ir0_check etasum ',D0ms1ir0_check

            endif            
            
            D0ms1ir0_check = D0ms1ir0_check - (
     &          cspcos_coli(-y(2),m2/(m1*r12),-iy(2),-ir12)     
     &          -.5d0*
     &          (cln_coli(-y(2),-iy(2))+cln_coli(m2/m1*r12,ir12))**2
     &          -cspcos_coli(cd0,m1/(m2*r12),iy(2),-ir12)
     &          -cspcos_coli(-y(2),m4/(m1*r14),-iy(2),-ir14)
     &          +.5d0*
     &          (cln_coli(-y(2),-iy(2))+cln_coli(m4/m1*r14,ir14))**2
     &          +cspcos_coli(cd0,m1/(m4*r14),iy(2),-ir14)
     &          +cspcos_coli(-z(2),m4/(m3*r34),-iz(2),-ir34)
     &          -.5d0*
     &          (cln_coli(-z(2),-iz(2))+cln_coli(m4/m3*r34,ir34))**2
     &          -cspcos_coli(cd0,m3/(m4*r34),iz(2),-ir34)
     &          -cspcos_coli(-z(2),m2/(m3*r23),-iz(2),-ir23)
     &          +.5d0*
     &          (cln_coli(-z(2),-iz(2))+cln_coli(m2/m3*r23,ir23))**2
     &          +cspcos_coli(cd0,m3/(m2*r23),iz(2),-ir23)
     &          ) 
            
c>            write(*,*) 'D0ms1ir0_check args yz',
c>     &          -y(2),-iy(2),-z(2),-iz(2)
c>            write(*,*) 'D0ms1ir0_check args r',
c>     &          m2/(m1*r12),m2/m1*r12,m4/(m1*r14),m4/m1*r14,
c>     &          m4/(m3*r34),m4/m3*r34,m2/(m3*r23),m2/m3*r23
c>            write(*,*) 'D0ms1ir0_check args ',
c>     &          1d0+y(2)*m2/(m1*r12),
c>     &          1d0+y(2)*m2/m1*r12,
c>     &          1d0+y(2)*m4/(m1*r14),
c>     &          1d0+y(2)*m4/m1*r14,
c>     &          1d0+z(2)*m4/(m3*r34),
c>     &          1d0+z(2)*m4/m3*r34,
c>     &          1d0+z(2)*m2/(m3*r23),
c>     &          1d0+z(2)*m2/m3*r23
c>            write(*,*) 'D0ms1ir0_check csp ',
c>     &          cspcos_coli(-y(2),m2/(m1*r12),-iy(2),-ir12),     
c>     &          -.5d0*
c>     &          (cln_coli(-y(2),-iy(2))+cln_coli(m2/m1*r12,ir12))**2
c>     &          -cspcos_coli(cd0,m1/(m2*r12),iy(2),-ir12),
c>     &          -cspcos_coli(-y(2),m4/(m1*r14),-iy(2),-ir14),
c>     &          +.5d0*
c>     &          (cln_coli(-y(2),-iy(2))+cln_coli(m4/m1*r14,ir14))**2
c>     &          +cspcos_coli(cd0,m1/(m4*r14),iy(2),-ir14),
c>     &          +cspcos_coli(-z(2),m4/(m3*r34),-iz(2),-ir34),
c>     &          -.5d0*
c>     &          (cln_coli(-z(2),-iz(2))+cln_coli(m4/m3*r34,ir34))**2
c>     &          -cspcos_coli(cd0,m3/(m4*r34),iz(2),-ir34),
c>     &          -cspcos_coli(-z(2),m2/(m3*r23),-iz(2),-ir23),
c>     &          +.5d0*
c>     &          (cln_coli(-z(2),-iz(2))+cln_coli(m2/m3*r23,ir23))**2
c>     &          +cspcos_coli(cd0,m3/(m2*r23),iz(2),-ir23)
c>            
c>            write(*,*) 'D0ms1ir0_check cspsum ',D0ms1ir0_check
            

c     *** improves numerical stability in extreme cases

            eta=eta2s_coli(-y(2),-1/z(2),-iy(2),iz(2),-ir13)
            
            if(eta.ne.0d0)then
              h14 = mm42*y(2)*y(2)+l14*y(2)
              h12 = mm22*y(2)*y(2)+l12*y(2)
              h34 = mm42*z(2)*z(2)+l34*z(2)
              h23 = mm22*z(2)*z(2)+l23*z(2)
              
              rt14 = abs(h14)/max(abs(mm42*y(2)*y(2)),abs(l14*y(2)))
              rt12 = abs(h12)/max(abs(mm22*y(2)*y(2)),abs(l12*y(2)))
              rt34 = abs(h34)/max(abs(mm42*z(2)*z(2)),abs(l34*z(2)))
              rt23 = abs(h23)/max(abs(mm22*z(2)*z(2)),abs(l23*z(2)))
              
              if (min(rt12,rt14).gt.min(rt34,rt23)) then
                if(dimag(h14/h12).ne.0d0)then
                  lqbard= cln_coli(h14/h12,0d0)
                else
                  lqbard=cln_coli(h14,-1d0)-cln_coli(h12,-1d0)
                endif
              else
                if(dimag(h34/h23).ne.0d0)then
                  lqbard= cln_coli(h34/h23,0d0)
                else
                  lqbard=cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
                endif
              endif
              
c>              write(*,*) 'D0ms1ir0_check lqbard',h14/h12,h34/h23
c>              write(*,*) 'D0ms1ir0_check lqbard',
c>     &            cln_coli(h14,-1d0)-cln_coli(h12,-1d0),
c>     &            cln_coli(h34,-1d0)-cln_coli(h23,-1d0)
        
              D0ms1ir0_check = D0ms1ir0_check 
     &            - eta*lqbard

c>              write(*,*) 'D0ms1ir0_check etacon ', -lqbard*eta
c>              write(*,*) 'D0ms1ir0_check etasum ',D0ms1ir0_check

            endif
      
c            write(*,*) 'D0ms1ir0_check sum ',D0ms1ir0_check

            D0ms1ir0_check = D0ms1ir0_check/(m1*m2*m3*m4*b)

c            write(*,*) 'D0ms1ir0_check =   ',D0ms1ir0_check         
      
          if(cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0).gt.testacc) then
            write(testout,100) ' D0ms1ir0_coli inconsistent d'
            write(testout,111) ' D0ms1ir0_coli: p12 = ',p12
            write(testout,111) ' D0ms1ir0_coli: p23 = ',p23
            write(testout,111) ' D0ms1ir0_coli: p34 = ',p34
            write(testout,111) ' D0ms1ir0_coli: p14 = ',p14
            write(testout,111) ' D0ms1ir0_coli: p24 = ',p24
            write(testout,111) ' D0ms1ir0_coli: p13 = ',p13
            write(testout,111) ' D0ms1ir0_coli: m12 = ',m12
            write(testout,111) ' D0ms1ir0_coli: m22 = ',m22
            write(testout,111) ' D0ms1ir0_coli: m32 = ',m32
            write(testout,111) ' D0ms1ir0_coli: m42 = ',m42
            write(testout,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
            write(testout,111) ' D0ms1ir0_check = ',D0ms1ir0_check
            write(testout,101) ' diff        = ',
     &          cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
            
            write(*,100) ' D0ms1ir0_coli inconsistent d'
            write(*,111) ' D0ms1ir0_coli: p12 = ',p12
            write(*,111) ' D0ms1ir0_coli: p23 = ',p23
            write(*,111) ' D0ms1ir0_coli: p34 = ',p34
            write(*,111) ' D0ms1ir0_coli: p14 = ',p14
            write(*,111) ' D0ms1ir0_coli: p24 = ',p24
            write(*,111) ' D0ms1ir0_coli: p13 = ',p13
            write(*,111) ' D0ms1ir0_coli: m12 = ',m12
            write(*,111) ' D0ms1ir0_coli: m22 = ',m22
            write(*,111) ' D0ms1ir0_coli: m32 = ',m32
            write(*,111) ' D0ms1ir0_coli: m42 = ',m42
            write(*,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
            write(*,111) ' D0ms1ir0_check = ',D0ms1ir0_check
            write(*,101) ' diff        = ',
     &          cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
c            stop
          endif
        endif
#endif

      else if(m12.ne.cd0.or.m22.ne.cd0)then

#ifdef CHECK
          if(m32.ne.m42.or.p34.ne.0d0.or.
     &        p14.eq.m12.or.p23.eq.m22)
     &        then
            write(*,100)' D0ms1ir0_coli:  (3.14)  m12=/=0, m22=0',
     &          '    case with 0 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms1ir0_coli: p12 = ',p12
            write(*,111)' D0ms1ir0_coli: p23 = ',p23
            write(*,111)' D0ms1ir0_coli: p34 = ',p34
            write(*,111)' D0ms1ir0_coli: p14 = ',p14
            write(*,111)' D0ms1ir0_coli: p24 = ',p24
            write(*,111)' D0ms1ir0_coli: p13 = ',p13
            write(*,111)' D0ms1ir0_coli: m12 = ',m12
            write(*,111)' D0ms1ir0_coli: m22 = ',m22
            write(*,111)' D0ms1ir0_coli: m32 = ',m32
            write(*,111)' D0ms1ir0_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef ALLCHECK
c ppenj
          if(flag(5))then
            write(*,100)' D0ms1ir0_coli:  (3.14)  m12=/=0, m22=0',
     &          '    case with 0 zero masses:',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms1ir0_coli: p12 = ',p12
            write(*,111)' D0ms1ir0_coli: p23 = ',p23
            write(*,111)' D0ms1ir0_coli: p34 = ',p34
            write(*,111)' D0ms1ir0_coli: p14 = ',p14
            write(*,111)' D0ms1ir0_coli: p24 = ',p24
            write(*,111)' D0ms1ir0_coli: p13 = ',p13
            write(*,111)' D0ms1ir0_coli: m12 = ',m12
            write(*,111)' D0ms1ir0_coli: m22 = ',m22
            write(*,111)' D0ms1ir0_coli: m32 = ',m32
            write(*,111)' D0ms1ir0_coli: m42 = ',m42
            flag(5)=.false.
          endif
#endif

c shift zero mass to 3 and small masses to 1 and 4
        if(m22.eq.cd0)then
          mm12=m42
          mm22=m12
          mm42=m32
          q13=p24
          q23=p12
          q24=p13
          q34=p23
          q12=p14
        else
          mm12=m32
          mm22=m22
          mm42=m42
          q13=p13
          q23=p12
          q24=p24
          q34=p14
          q12=p23
        endif

c     cD0ms0m
c     Denner Dittmaier (3.14) with  zero mass 
        mm12=mm12*coliminfscale2
        mm42=mm42*coliminfscale2


        l12 = (mm22   -q12)
        l13 = (       -q13)
        l14 = (mm42+mm12    )
        l23 = (mm22   -q23)
        l24 = (mm22   -q24)
        l34 = (       -q34)
        
        r12 = l12/mm22
        r21 = l12/mm12      
        r24 = mm22/l24
        r42 = mm42/l24      
        r14 = 1d0
        r41 = mm42/mm12
        
        if(dreal(l12).lt.-0d0) then
          ir12 = -1d0
        else
          ir12 = 0d0
        end if
        if(dreal(l24).lt.-0d0) then
          ir24 = 1d0
        else
          ir24 = 0d0
        end if
        
        y(1) = -l12/l24
        y(2) = -l13/l34
        z(1) = y(1)/r24
        z(2) = y(2)/r24

c>        d=l23-l34/l24*m22
c>        iy(1) = -sign(1d0,dreal(d))
c>        iy(2) = +sign(1d0,dreal(d))
c>        iz(1) =  sign(1d0,iy(1)*dreal(r24))
c>        iz(2) =  sign(1d0,iy(2)*dreal(r24))

        iz(1) = -ir12           ! ensures correct continution
        iy(2) = dreal(l34-l13)  ! ensures correct continution
        iz(2) =  sign(1d0,iy(2)*dreal(r24))
        iy(1) =  sign(1d0,iz(1)*dreal(r24))

c>        write(*,*) 'D0ms1ir0_coli q ',q12,q23,q34,q14,q13,q24
c>        write(*,*) 'D0ms1ir0_coli m ',m12,m22,m32,m42
c>        write(*,*) 'D0ms1ir0_coli l ',l12,l13,l14
c>        write(*,*) 'D0ms1ir0_coli l ',l23,l24,l34
c>        write(*,*) 'D0ms1ir0_coli r ',r12,r14,r24
c>        write(*,*) 'D0ms1ir0_coli r ',r21,r41,r42
c>        write(*,*) 'D0ms1ir0_coli y ',y(1),y(2)
c>        write(*,*) 'D0ms1ir0_coli z ',z(1),z(2)
c>        write(*,*) 'D0ms1ir0_coli ir ',ir12,ir14,ir24
c>        write(*,*) 'D0ms1ir0_coli iy ',iy(1),iy(2)
c>        write(*,*) 'D0ms1ir0_coli iz ',iz(1),iz(2)
        
        D0ms1ir0_coli = dcmplx(0d0)
        dO 11 i=2,1,-1
          eta=eta2s_coli(-y(i),1d0/r24,-iy(i),-ir24,-iz(i))

          if(i.eq.1)then
            if(eta.ne.0d0)then
c changed 17.10.08
              if(dreal(y(i)).gt.0d0)then
                iqbard(i)=dreal(m12/y(i)+l14+m42*y(i) - (l34+l13/y(i)))
              else
                iqbard(i) = dreal(- (l23+l13/z(i)))
              endif
              l1(1) = cln_coli((mm12-l14*l12/l24+mm42*(l12/l24)**2)
     &            /(l13-l34*l12/l24),iqbard(i))

c>            write(*,*) 'D0ms1ir0_coli eta1 ',
c>     &          (mm12-l14*l12/l24+m42*(l12/l24)**2)
c>     &          /(l13-l34*l12/l24)            
c>            write(*,*) 'D0ms1ir0_coli l1(1)= ',l1(1)
c>            write(*,*) 'D0ms1ir0_coli l2args= ',
c>     &          (mm12+l14*y(i)+mm42*y(i)*y(i))/(l13+l34*y(i))
c>     &                *(l23-r24*l34),
c>     &          l23-r24*l34,
c>     &          mm12+l14*y(i)+mm42*y(i)*y(i),
c>     &          l13+l34*y(i),
c>     &          (mm12+l14*y(i)+mm42*y(i)*y(i))/(l13+l34*y(i))
c>     &          *(l13+l23*z(i)),
c>     &          l13+l23*z(i)
c>            write(*,*) 'D0ms1ir0_coli l2= ',
c>     &          cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))/(l13+l34*y(i))
c>     &          *(l23-r24*l34),-1d0)
c>     &          -cln_coli(l23-r24*l34,-1d0),
c>     &          cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
c>     &          -cln_coli(l13+l34*y(i),-1d0),
c>     &          cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))/(l13+l34*y(i))
c>     &          *(l13+l23*z(i)),-1d0)
c>     &          -cln_coli(l13+l23*z(i),-1d0)

#ifdef CHECK
              if(conscheck)then
                if(dimag(r24).ne.0d0) then
                  v=dimag(y(i))/dimag(r24)
                  u=dimag(y(i)/r24)/dimag(1d0/r24)
                  l1(2)=cln_coli(l12*v+mm22*v*v+l24*u*v,-1d0)
     &                -cln_coli(l13+l23*v+l34*u,-1d0)

c                  write(*,*) 'D0ms1ir0_coli comp1 l1(2)= ',
c     &                l1(2),l1(1),l1(2)/l1(1)

                elseif(dimag(y(i)).ne.0d0) then
                  l1(2)=cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))
     &                /(l13+l34*y(i))*(l23-r24*l34),-1d0)
     &                -cln_coli(l23-r24*l34,-1d0)
                elseif(dreal(y(i)).gt.0d0)then
                  l1(2)=cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
     &                -cln_coli(l13+l34*y(i),-1d0)
                elseif(dreal(z(i)).gt.0d0)then
                  l1(2)=cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))
     &                /(l13+l34*y(i))*(l13+l23*z(i)),-1d0)
     &                -cln_coli(l13+l23*z(i),-1d0)
                else
                  if (flag3) then
                    write(*,*) 
                    write(*,*) 
     &                  'D0ms1ir0_coli: case for l1(2) not implemented '
                    flag3 = .false.
                  end if
                  write(errout,*) 
                  write(errout,*) 
     &                'D0ms1ir0_coli: case for l1(2) not implemented '
                  write(errout,*) 
     &                'D0ms1ir0_coli args ',q12,q23,q34,q14,q13,q24,
     &                mm12,mm22,m32,mm42
                  write(errout,*) 'r24 = ',r24
                  write(errout,*) 'y(i) = ',y(i)
                  write(errout,*) 'z(i) = ',z(i)
                endif 
                
                if(abs(l1(2)/l1(1)-1d0).gt.testacc) then
                  if (flag2) then
                    write(*,*) 
                    write(*,*)
     &                  'D0ms1ir0_coli: eta terms possibly wrong a'
                    flag2 = .false.
                  end if
                  write(errout,*) 
                  write(errout,*) 
     &                'D0ms1ir0_coli: eta terms  possibly wrong a'
                  write(errout,*) 
     &                'D0ms1ir0_coli args ',q12,q23,q34,q14,q13,q24,
     &                mm12,mm22,m32,mm42
                  write(errout,*) 'l1(1) = ',l1(1)
                  write(errout,*) 'l1(2) = ',l1(2)
                  write(errout,*) 'diff   = ',abs((l1(2)-l1(1))/l1(1))
                endif 
              endif
#endif
            endif
            sp3 =  cspcos_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13))

c           sp4 = cspcos_coli(-z(i),1d0/r12,-iz(i),-ir12)
c>            imp4=-ir12*sign(1d0,dreal(-z(i)))
c>            if(imp4.eq.0d0) imp4 =-iz(i)*sign(1d0,dreal(r12))
c>            sp4 = eta2s_coli(-z(i),1d0/r12,-iz(i),-ir12,imp4)
c>     &          *cln_coli((1d0+y(i))*(1d0+y(i)*mm42/mm12)
c>     &          *(1d0+l23*z(i)/l13)
c>     &          /((1d0+y(i)*l34/l13)*z(i)*r21),-imp4)
c>
c>            write(*,*) 'D0ms1ir0_coli sp4= ',sp4
c>
c>            sp4=cspco0_coli(-z(i),1d0/r12,
c>     &          (1d0+y(i))*(1d0+y(i)*mm42/mm12)*(1d0+l23*z(i)/l13)
c>     &          /((1d0+y(i)*l34/l13)*z(i)*r21),
c>     &          -iz(i),-ir12)

            sp4=0        ! enforce by choice of iz(i)

c            write(*,*) 'D0ms1ir0_coli sp4= ',sp4

          else
c changed 17.10.08
            if(dreal(y(i)).gt.0d0)then
              iqbard(i) = dreal(mm12+l14*y(i)+mm42*y(i)*y(i))
            else
              iqbard(i) = dreal(l12+m22*z(i) - (l23+l13/z(i)))
            endif
            if(eta.ne.0d0)then
              l1(1) = cln_coli((l12-l24*l13/l34)/(l23-mm22*l34/l24),
     &            iqbard(i))

c>            write(*,*) 'D0ms1ir0_coli eta2 ',
c>     &          (l12-l24*l13/l34)/(l23-mm22*l34/l24)
c>            write(*,*) 'D0ms1ir0_coli l1(1)= ',l1(1)
c>            write(*,*) 'D0ms1ir0_coli l2= ',
c>     &          cln_coli(l12+mm42*y(i)/r42,-1d0)
c>     &          -cln_coli(l23-r24*l34,-1d0),
c>     &          cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
c>     &          -cln_coli((1d0+y(2))*(1d0+y(2)*mm42/mm12)
c>     &          *(1d0+l23*z(2)/l13)/((1+z(2)/r12)*z(2)*r21),-1d0),
c>     &          cln_coli(l12*z(i)+mm22*z(i)*z(i),-1d0)
c>     &          -cln_coli(l13+l23*z(i),-1d0)

#ifdef CHECK
              if(conscheck)then
                if(dimag(r24).ne.0d0) then
                  v=dimag(y(i))/dimag(r24)
                  u=dimag(y(i)/r24)/dimag(1d0/r24)
                  if(v.ne.0)then
                    l1(2)=cln_coli(l12*v+mm22*v*v+l24*u*v,-1d0)
     &                  -cln_coli(l13+l23*v+l34*u,-1d0)

c                    write(*,*) 'D0ms1ir0_coli comp21 l1(2)= ',
c     &                  l1(2),l1(1),l1(2)/l1(1)

                  else
                    l1(2)=cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
     &                  -cln_coli((1d0+y(2))*(1d0+y(2)*mm42/mm12)
     &                  *(l13+l23*z(2))/
     &                  ((1+z(2)/r12)*z(2)*r21),-1d0)

c                    write(*,*) 'D0ms1ir0_coli comp22 l1(2)= ',
c     &                  l1(2),l1(1),l1(2)/l1(1)

                  endif
                elseif(dimag(y(i)).ne.0d0) then
                  l1(2)=cln_coli(l12+l24*y(i),-1d0)
     &                -cln_coli(l23-r24*l34,-1d0)
                elseif(dreal(y(i)).gt.0d0)then
                  l1(2)=cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
     &                -cln_coli((1d0+y(2))*(1d0+y(2)*mm42/mm12)
     &                *(l13+l23*z(2))/((1+z(2)/r12)*z(2)*r21),-1d0)
                elseif(dreal(z(i)).gt.0d0)then
                  l1(2)=cln_coli(l12*z(i)+mm22*z(i)*z(i),-1d0)
     &                -cln_coli(l13+l23*z(i),-1d0)
                else
                  if (flag3) then
                    write(*,*) 
                    write(*,*) 
     &                  'D0ms1ir0_coli: case for l1(2) not implemented '
                    flag3 = .false.
                  end if
                  write(errout,*) 
                  write(errout,*) 
     &                'D0ms1ir0_coli: case for l1(2) not implemented'
                  write(errout,*) 'D0ms1ir0_coli args ',q12,q23,q34,q14,q13,
     &                q24,mm12,mm22,m32,mm42
                  write(errout,*) 'r24 = ',r24
                  write(errout,*) 'y(i) = ',y(i)
                  write(errout,*) 'z(i) = ',z(i)
                endif 
                
                if(abs(l1(2)/l1(1)-1d0).gt.testacc) then
                  if (flag2) then
                    write(*,*) 
                    write(*,*)
     &                  'D0ms1ir0_coli: eta terms possibly wrong b'
                    flag2 = .false.
                  end if
                  write(errout,*) 
                  write(errout,*) 
     &                'D0ms1ir0_coli: eta terms possibly wrong b'
                  write(errout,*) 
     &                'D0ms1ir0_coli args ',q12,q23,q34,q14,q13,
     &                q24,mm12,mm22,m32,mm42
                  write(errout,*) 'l1(1) = ',l1(1)
                  write(errout,*) 'l1(2) = ',l1(2)
                  write(errout,*) 'diff   = ',abs((l1(2)-l1(1))/l1(1))
                endif 
              endif
#endif
            endif
c>            imp3=dreal(l34-l13)*sign(1d0,dreal(-y(i)))
c>            sp3 = eta2s_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13),imp3)
c>     &          *cln_coli((1d0+y(i))*(1d0+y(i)*mm42/mm12)
c>     &          *(1d0+l23*z(i)/l13)
c>     &          /((1d0+z(i)/r12)*z(i)*r21),-imp3)

c>            write(*,*) 'D0ms1ir0_coli sp3 ', 
c>     &           eta2s_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13),imp3),
c>     &          (1d0+y(i))*(1d0+y(i)*mm42/mm12)
c>     &          *(1d0+l23*z(i)/l13)
c>     &          /((1d0+z(i)/r12)*z(i)*r21)
c>            write(*,*) 'D0ms1ir0_coli sp3 ',sp3,-imp3

c>            sp3=cspco0_coli(-y(i),l34/l13,
c>     &          (1d0+y(i))*(1d0+y(i)*mm42/mm12)*(1d0+l23*z(i)/l13)
c>     &          /((1d0+z(i)/r12)*z(i)*r21),
c>     &          -iy(i),dreal(l34-l13))

            sp3 = 0 ! enforce by choicxe of iy(i)

c            Write(*,*) 'D0ms1ir0_coli sp3 ',sp3

            sp4 = cspcos_coli(-z(i),1d0/r12,-iz(i),-ir12)

c>            write(*,*) 'D0ms1ir0_coli test1',1d0+y(i),1d0+y(i)*mm42/mm12,
c>     &          (1d0+l23*z(i)/l13),(1d0+z(i)/r12),z(i),r21
c>            write(*,*) 'D0ms1ir0_coli test1',1d0+y(i),1d0+y(i)*mm42/mm12,
c>     &          (1d0+y(i))*(1d0+y(i)*mm42/mm12)
c>     &          *(1d0+l23*z(i)/l13)
c>     &          /((1d0+z(i)/r12)*z(i)*r21),
c>     &          (1d0+y(i))*(1d0+y(i)*mm42/mm12)/
c>     &          ((1d0+y(i))*(1d0+y(i)*mm42/mm12)
c>     &          *(1d0+l23*z(i)/l13)
c>     &          /((1d0+z(i)/r12)*z(i)*r21))
c>            write(*,*) 'D0ms1ir0_coli test2',
c>     &          (l12-l24*l13/l34)/(l23-mm22*l34/l24),l13/mm12,
c>     &          (l12-l24*l13/l34)/(l23-mm22*l34/l24)*l13/mm12

          endif

c>          write(*,*) 'D0ms1ir0_coli contr ',
c>     &        (2*i-3) *cspcos_coli(-y(i),mm42/mm12,-iy(i),0d0)
c>          write(*,*) 'D0ms1ir0_coli contr ',
c>     &        + (2*i-3) *cspcos_coli(-y(i),cd1,-iy(i),0d0)
c>          write(*,*) 'D0ms1ir0_coli contr ',
c>c     &           - cspcos_coli(-z(i),r21,-iz(i),ir12)
c>     &        + (2*i-3) *.5d0*(cln_coli(-z(i),-iz(i))
c>     &        +cln_coli(r21,ir12))**2
c>     &        + (2*i-3) *cspcos_coli(dcmplx(1d0),cd0,0d0,0d0)
c>          write(*,*) 'D0ms1ir0_coli contr ',
c>c     &           - cspcos_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13))
c>     &        -(2*i-3) *sp3
c>          write(*,*) 'D0ms1ir0_coli contr ',
c>     &        - (2*i-3) *sp4
c>          write(*,*) 'D0ms1ir0_coli contr ',
c>     &        + (2*i-3) 
c>     &        *cspcos_coli(-z(i),l23/l13,-iz(i),dreal(l23-l13))
c>          write(*,*) 'D0ms1ir0_coli contr ',
c>     &        - (2*i-3) *eta*(l1(1)+cln_coli(l13/mm12,-1d0))           
c>          write(*,*) 'D0ms1ir0_coli eta ',
c>     &        eta,(l1(1)+cln_coli(l13/mm12,-1d0)) 
c>          write(*,*) 'D0ms1ir0_coli eta ',
c>     &        l1(1),cln_coli(l13/mm12,-1d0)
c>          write(*,*) 'D0ms1ir0_coli eta2 ',
c>     &        (l12-l24*l13/l34),(l23-mm22*l34/l24),
c>     &        (l12-l24*l13/l34)/(l23-mm22*l34/l24)
          
          D0ms1ir0_coli = D0ms1ir0_coli + (2*i-3) * (
     &        cspcos_coli(-y(i),mm42/mm12,-iy(i),0d0)

     &        + cspcos_coli(-y(i),cd1,-iy(i),0d0)
c     &           - cspcos_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13))
     &        -sp3
     &        + cspcos_coli(cd1,cd0,0d0,0d0)
     &        + .5d0*(cln_coli(-z(i),-iz(i))+cln_coli(r21,ir12))**2
c     &           - cspcos_coli(-z(i),r21,-iz(i),ir12)
     &        - sp4
     &        + cspcos_coli(-z(i),l23/l13,-iz(i),dreal(l23-l13))
     &        - eta*(l1(1)+cln_coli(l13/mm12,-1d0)) 
     &        )

c          write(*,*) 'D0ms1ir0_coli sum ',D0ms1ir0_coli

 11     continue
        
c        write(*,*) 'D0ms1ir0_coli sum ',D0ms1ir0_coli
         
        D0ms1ir0_coli = -D0ms1ir0_coli/(l12*l34-l13*l24)

#ifdef CHECK
        if (conscheck) then

          do j=0,3
            if(j.eq.0) then
              iy(1) = -sign(1d0,dreal(d))
              iy(2) = +sign(1d0,dreal(d))
            elseif(j.eq.1) then
              iy(1) = +sign(1d0,dreal(d))
              iy(2) = -sign(1d0,dreal(d))
            elseif(j.eq.2) then
              iy(1) = +sign(1d0,dreal(d))
              iy(2) = +sign(1d0,dreal(d))
            elseif(j.eq.3) then
              iy(1) = -sign(1d0,dreal(d))
              iy(2) = -sign(1d0,dreal(d))
            endif
            iz(1) =  sign(1d0,iy(1)*dreal(r24))
            iz(2) =  sign(1d0,iy(2)*dreal(r24))
            
c>            write(*,*) 'D0ms1ir0_check m ',m12,m22,m32,m42
c>            write(*,*) 'D0ms1ir0_check l ',l12,l13,l14
c>            write(*,*) 'D0ms1ir0_check l ',l23,l24,l34
c>            write(*,*) 'D0ms1ir0_check r ',r12,r14,r24
c>            write(*,*) 'D0ms1ir0_check r ',r21,r41,r42
c>            write(*,*) 'D0ms1ir0_check y ',y(1),y(2)
c>            write(*,*) 'D0ms1ir0_check z ',z(1),z(2)
c>            write(*,*) 'D0ms1ir0_check ir ',ir12,ir14,ir24
c>            write(*,*) 'D0ms1ir0_check iy ',iy(1),iy(2)
c>            write(*,*) 'D0ms1ir0_check iz ',iz(1),iz(2)

c change iz and iy if required by continution 
            imp4=-ir12*sign(1d0,dreal(-z(1)))
            if(imp4.eq.0d0) imp4 =-iz(1)*sign(1d0,dreal(r12))
            if(eta2s_coli(-z(1),1d0/r12,-iz(1),-ir12,imp4).ne.cd0
     &          .and.dimag((1d0+y(1))*(1d0+y(1)*mm42/mm12)
     &          *(1d0+l23*z(1)/l13)/((1d0+y(1)*l34/l13)*z(1)*r21))
     &          *imp4.gt.0d0) then
              iz(1)=-iz(1)
              iy(1)=-iy(1)
            endif
            imp3=dreal(l34-l13)*sign(1d0,dreal(-y(2)))
            if(eta2s_coli(-y(2),l34/l13,-iy(2),dreal(l34-l13),imp3)
     &          .ne.cd0
     &          .and.dimag((1d0+y(2))*(1d0+y(2)*mm42/mm12)
     &              *(1d0+l23*z(2)/l13)/((1d0+z(2)/r12)*z(2)*r21))
     &          *imp3.gt.0d0)then
              iz(2)=-iz(2)
              iy(2)=-iy(2)
            endif
            
            D0ms1ir0_check = dcmplx(0d0)
            dO 12 i=2,1,-1
              eta=eta2s_coli(-y(i),1d0/r24,-iy(i),-ir24,-iz(i))

              if(i.eq.1)then
                if(eta.ne.0d0)then
c     changed 17.10.08
                  if(dreal(y(i)).gt.0d0)then
                    iqbard(i)=dreal(m12/y(i)+l14+m42*y(i) 
     &                  - (l34+l13/y(i)))
                  else
                    iqbard(i) = dreal(- (l23+l13/z(i)))
                  endif
                  l1(1) = cln_coli((mm12-l14*l12/l24+mm42*(l12/l24)**2)
     &                /(l13-l34*l12/l24),iqbard(i))

c>            write(*,*) 'D0ms1ir0_check eta1 ',
c>     &          (mm12-l14*l12/l24+m42*(l12/l24)**2)
c>     &          /(l13-l34*l12/l24)            
c>            write(*,*) 'D0ms1ir0_check l1(1)= ',l1(1)
c>            write(*,*) 'D0ms1ir0_check l2args= ',
c>     &          (mm12+l14*y(i)+mm42*y(i)*y(i))/(l13+l34*y(i))
c>     &                *(l23-r24*l34),
c>     &          l23-r24*l34,
c>     &          mm12+l14*y(i)+mm42*y(i)*y(i),
c>     &          l13+l34*y(i),
c>     &          (mm12+l14*y(i)+mm42*y(i)*y(i))/(l13+l34*y(i))
c>     &          *(l13+l23*z(i)),
c>     &          l13+l23*z(i)
c>            write(*,*) 'D0ms1ir0_check l2= ',
c>     &          cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))/(l13+l34*y(i))
c>     &          *(l23-r24*l34),-1d0)
c>     &          -cln_coli(l23-r24*l34,-1d0),
c>     &          cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
c>     &          -cln_coli(l13+l34*y(i),-1d0),
c>     &          cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))/(l13+l34*y(i))
c>     &          *(l13+l23*z(i)),-1d0)
c>     &          -cln_coli(l13+l23*z(i),-1d0)

                  if(dimag(r24).ne.0d0) then
                    v=dimag(y(i))/dimag(r24)
                    u=dimag(y(i)/r24)/dimag(1d0/r24)
                    l1(2)=cln_coli(l12*v+mm22*v*v+l24*u*v,-1d0)
     &                  -cln_coli(l13+l23*v+l34*u,-1d0)

c     write(*,*) 'D0ms1ir0_check comp1 l1(2)= ',
c     &                l1(2),l1(1),l1(2)/l1(1)

                  elseif(dimag(y(i)).ne.0d0) then
                    l1(2)=cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))
     &                  /(l13+l34*y(i))*(l23-r24*l34),-1d0)
     &                  -cln_coli(l23-r24*l34,-1d0)
                  elseif(dreal(y(i)).gt.0d0)then
                    l1(2)=cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
     &                  -cln_coli(l13+l34*y(i),-1d0)
                  elseif(dreal(z(i)).gt.0d0)then
                    l1(2)=cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))
     &                  /(l13+l34*y(i))*(l13+l23*z(i)),-1d0)
     &                  -cln_coli(l13+l23*z(i),-1d0)
                  else
                    if (flag3) then
                      write(*,*) 
                      write(*,*) 
     &                  'D0ms1ir0_check: case for l1(2) not implemented'
                      flag3 = .false.
                    end if
                    write(errout,*) 
                    write(errout,*) 
     &                  'D0ms1ir0_check: case for l1(2) not implemented'
                    write(errout,*) 
     &                  'D0ms1ir0_check args ',q12,q23,q34,q14,q13,q24,
     &                  mm12,mm22,m32,mm42
                    write(errout,*) 'r24 = ',r24
                    write(errout,*) 'y(i) = ',y(i)
                    write(errout,*) 'z(i) = ',z(i)
                  endif 
                  
                  if(abs(l1(2)/l1(1)-1d0).gt.testacc) then
                    if (flag2) then
                      write(*,*) 
                      write(*,*)
     &                    'D0ms1ir0_check: eta terms possibly wrong a'
                      flag2 = .false.
                    end if
                    write(errout,*) 
                    write(errout,*) 
     &                  'D0ms1ir0_check: eta terms  possibly wrong a'
                    write(errout,*) 
     &                  'D0ms1ir0_check args ',q12,q23,q34,q14,q13,q24,
     &                  mm12,mm22,m32,mm42
                    write(errout,*) 'l1(1) = ',l1(1)
                    write(errout,*) 'l1(2) = ',l1(2)
                    write(errout,*) 'diff   = ',abs((l1(2)-l1(1))/l1(1))
                  endif 
                endif
                sp3 =  cspcos_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13))
c     sp4 = cspcos_coli(-z(i),1d0/r12,-iz(i),-ir12)
                imp4=-ir12*sign(1d0,dreal(-z(i)))
                if(imp4.eq.0d0) imp4 =-iz(i)*sign(1d0,dreal(r12))
                sp4 = eta2s_coli(-z(i),1d0/r12,-iz(i),-ir12,imp4)
     &              *cln_coli((1d0+y(i))*(1d0+y(i)*mm42/mm12)
     &              *(l13+l23*z(i))
     &              /((l13+y(i)*l34)*z(i)*r21),-imp4)

c                write(*,*) 'D0ms1ir0_check sp4= ',sp4

                sp4=cspco0_coli(-z(i),1d0/r12,
     &              (1d0+y(i))*(1d0+y(i)*mm42/mm12)*(l13+l23*z(i))
     &              /((l13+y(i)*l34)*z(i)*r21),
     &              -iz(i),-ir12)

c                write(*,*) 'D0ms1ir0_check sp4= ',sp4

              else
c     changed 17.10.08
                if(dreal(y(i)).gt.0d0)then
                  iqbard(i) = dreal(mm12+l14*y(i)+mm42*y(i)*y(i))
                else
                  iqbard(i) = dreal(l12+m22*z(i) - (l23+l13/z(i)))
                endif
                if(eta.ne.0d0)then
                  l1(1) = cln_coli((l12-l24*l13/l34)/(l23-mm22*l34/l24),
     &                iqbard(i))

c>            write(*,*) 'D0ms1ir0_check eta2 ',
c>     &          (l12-l24*l13/l34)/(l23-mm22*l34/l24)
c>            write(*,*) 'D0ms1ir0_check l1(1)= ',l1(1)
c>            write(*,*) 'D0ms1ir0_check l2= ',
c>     &          cln_coli(l12+mm42*y(i)/r42,-1d0)
c>     &          -cln_coli(l23-r24*l34,-1d0),
c>     &          cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
c>     &          -cln_coli((1d0+y(2))*(1d0+y(2)*mm42/mm12)
c>     &          *(1d0+l23*z(2)/l13)/((1+z(2)/r12)*z(2)*r21),-1d0),
c>     &          cln_coli(l12*z(i)+mm22*z(i)*z(i),-1d0)
c>     &          -cln_coli(l13+l23*z(i),-1d0)

                  if(dimag(r24).ne.0d0) then
                    v=dimag(y(i))/dimag(r24)
                    u=dimag(y(i)/r24)/dimag(1d0/r24)
                    if(v.ne.0)then
                      l1(2)=cln_coli(l12*v+mm22*v*v+l24*u*v,-1d0)
     &                    -cln_coli(l13+l23*v+l34*u,-1d0)

c     write(*,*) 'D0ms1ir0_check comp21 l1(2)= ',
c     &                  l1(2),l1(1),l1(2)/l1(1)

                    else
                      l1(2)=cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
     &                    -cln_coli((1d0+y(2))*(1d0+y(2)*mm42/mm12)
     &                    *(l13+l23*z(2))/
     &                    ((1+z(2)/r12)*z(2)*r21),-1d0)

c     write(*,*) 'D0ms1ir0_check comp22 l1(2)= ',
c     &                  l1(2),l1(1),l1(2)/l1(1)

                    endif
                  elseif(dimag(y(i)).ne.0d0) then
                    l1(2)=cln_coli(l12+l24*y(i),-1d0)
     &                  -cln_coli(l23-r24*l34,-1d0)
                  elseif(dreal(y(i)).gt.0d0)then
                    l1(2)=cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
     &                  -cln_coli((1d0+y(2))*(1d0+y(2)*mm42/mm12)
     &                  *(l13+l23*z(2))/((1+z(2)/r12)*z(2)*r21),-1d0)
                  elseif(dreal(z(i)).gt.0d0)then
                    l1(2)=cln_coli(l12*z(i)+mm22*z(i)*z(i),-1d0)
     &                  -cln_coli(l13+l23*z(i),-1d0)
                  else
                    if (flag3) then
                      write(*,*) 
                      write(*,*) 
     &                 'D0ms1ir0_check: case for l1(2) not implemented '
                      flag3 = .false.
                    end if
                    write(errout,*) 
                    write(errout,*) 
     &                  'D0ms1ir0_check: case for l1(2) not implemented'
                    write(errout,*) 'D0ms1ir0_check args ',
     &                  q12,q23,q34,q14,q13,
     &                  q24,mm12,mm22,m32,mm42
                    write(errout,*) 'r24 = ',r24
                    write(errout,*) 'y(i) = ',y(i)
                    write(errout,*) 'z(i) = ',z(i)
                  endif 
                  
                  if(abs(l1(2)/l1(1)-1d0).gt.testacc) then
                    if (flag2) then
                      write(*,*) 
                      write(*,*)
     &                    'D0ms1ir0_check: eta terms possibly wrong b'
                      flag2 = .false.
                    end if
                    write(errout,*) 
                    write(errout,*) 
     &                  'D0ms1ir0_check: eta terms possibly wrong b'
                    write(errout,*) 
     &                  'D0ms1ir0_check args ',q12,q23,q34,q14,q13,
     &                  q24,mm12,mm22,m32,mm42
                    write(errout,*) 'l1(1) = ',l1(1)
                    write(errout,*) 'l1(2) = ',l1(2)
                    write(errout,*) 'diff   = ',abs((l1(2)-l1(1))/l1(1))
                  endif 
                endif
                imp3=dreal(l34-l13)*sign(1d0,dreal(-y(i)))
                
                sp3=eta2s_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13),imp3)
     &              *cln_coli((1d0+y(i))*(1d0+y(i)*mm42/mm12)
     &              *(1d0+l23*z(i)/l13)
     &              /((1d0+z(i)/r12)*z(i)*r21),-imp3)

c>                write(*,*) 'D0ms1ir0_check sp3 ', 
c>     &              eta2s_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13),imp3)
c>     &              ,(1d0+y(i))*(1d0+y(i)*mm42/mm12)
c>     &              *(1d0+l23*z(i)/l13)
c>     &              /((1d0+z(i)/r12)*z(i)*r21)
c>                write(*,*) 'D0ms1ir0_check sp3 ',sp3,-imp3

                sp3=cspco0_coli(-y(i),l34/l13,
     &              (1d0+y(i))*(1d0+y(i)*mm42/mm12)*(1d0+l23*z(i)/l13)
     &              /((1d0+z(i)/r12)*z(i)*r21),
     &              -iy(i),dreal(l34-l13))

c                write(*,*) 'D0ms1ir0_check sp3 ',sp3

                sp4 = cspcos_coli(-z(i),1d0/r12,-iz(i),-ir12)
              endif

c>              write(*,*) 'D0ms1ir0_check contr ',
c>     &            (2*i-3) *cspcos_coli(-y(i),mm42/mm12,-iy(i),0d0)
c>              write(*,*) 'D0ms1ir0_check contr ',
c>     &            + (2*i-3) *cspcos_coli(-y(i),cd1,-iy(i),0d0)
c>              write(*,*) 'D0ms1ir0_check contr ',
c>c     &           - cspcos_coli(-z(i),r21,-iz(i),ir12)
c>     &            + (2*i-3) *.5d0*(cln_coli(-z(i),-iz(i))
c>     &            +cln_coli(r21,ir12))**2
c>     &            + (2*i-3) *cspcos_coli(dcmplx(1d0),cd0,0d0,0d0)
c>              write(*,*) 'D0ms1ir0_check contr ',
c>c     &           - cspcos_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13))
c>     &            -(2*i-3) *sp3
c>              write(*,*) 'D0ms1ir0_check contr ',
c>     &            - (2*i-3) *sp4
c>              write(*,*) 'D0ms1ir0_check contr ',
c>     &            + (2*i-3) 
c>     &            *cspcos_coli(-z(i),l23/l13,-iz(i),dreal(l23-l13))
c>              write(*,*) 'D0ms1ir0_check contr ',
c>     &            - (2*i-3) *eta*(l1(1)+cln_coli(l13/mm12,-1d0))           
c>              write(*,*) 'D0ms1ir0_check eta ',
c>     &            eta,(l1(1)+cln_coli(l13/mm12,-1d0)) 
c>              write(*,*) 'D0ms1ir0_check eta ',
c>     &            l1(1),cln_coli(l13/mm12,-1d0)
c>              write(*,*) 'D0ms1ir0_check eta2 ',
c>     &            (l12-l24*l13/l34),(l23-mm22*l34/l24),
c>     &            (l12-l24*l13/l34)/(l23-mm22*l34/l24)
              
              D0ms1ir0_check = D0ms1ir0_check + (2*i-3) * (
     &            cspcos_coli(-y(i),mm42/mm12,-iy(i),0d0)

     &            + cspcos_coli(-y(i),cd1,-iy(i),0d0)
c     &           - cspcos_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13))
     &            -sp3
     &            + cspcos_coli(cd1,cd0,0d0,0d0)
     &            + .5d0*(cln_coli(-z(i),-iz(i))+cln_coli(r21,ir12))**2
c     &           - cspcos_coli(-z(i),r21,-iz(i),ir12)
     &            - sp4
     &            + cspcos_coli(-z(i),l23/l13,-iz(i),dreal(l23-l13))
     &            - eta*(l1(1)+cln_coli(l13/mm12,-1d0)) 
     &            )

c              write(*,*) 'D0ms1ir0_check sum ',D0ms1ir0_check

 12         continue
            
c            write(*,*) 'D0ms1ir0_check sum ',D0ms1ir0_check
            
            D0ms1ir0_check = -D0ms1ir0_check/(l12*l34-l13*l24)

            if(cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0).gt.testacc) then
              write(testout,*) ' D0ms1ir0_coli inconsistent',i
              write(testout,111) ' D0ms1ir0_coli: p12 = ',p12
              write(testout,111) ' D0ms1ir0_coli: p23 = ',p23
              write(testout,111) ' D0ms1ir0_coli: p34 = ',p34
              write(testout,111) ' D0ms1ir0_coli: p14 = ',p14
              write(testout,111) ' D0ms1ir0_coli: p24 = ',p24
              write(testout,111) ' D0ms1ir0_coli: p13 = ',p13
              write(testout,111) ' D0ms1ir0_coli: m12 = ',m12
              write(testout,111) ' D0ms1ir0_coli: m22 = ',m22
              write(testout,111) ' D0ms1ir0_coli: m32 = ',m32
              write(testout,111) ' D0ms1ir0_coli: m42 = ',m42
              write(testout,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
              write(testout,111) ' D0ms1ir0_check = ',D0ms1ir0_check
              write(testout,101) ' diff        = ',
     &            cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
              
              write(*,*) ' D0ms1ir0_coli inconsistent',i
              write(*,111) ' D0ms1ir0_coli: p12 = ',p12
              write(*,111) ' D0ms1ir0_coli: p23 = ',p23
              write(*,111) ' D0ms1ir0_coli: p34 = ',p34
              write(*,111) ' D0ms1ir0_coli: p14 = ',p14
              write(*,111) ' D0ms1ir0_coli: p24 = ',p24
              write(*,111) ' D0ms1ir0_coli: p13 = ',p13
              write(*,111) ' D0ms1ir0_coli: m12 = ',m12
              write(*,111) ' D0ms1ir0_coli: m22 = ',m22
              write(*,111) ' D0ms1ir0_coli: m32 = ',m32
              write(*,111) ' D0ms1ir0_coli: m42 = ',m42
              write(*,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
              write(*,111) ' D0ms1ir0_check = ',D0ms1ir0_check
              write(*,101) ' diff        = ',
     &            cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
c     stop
            endif
 
          enddo

          l12 = (mm22    -q12)
          l13 = (        -q13)
          l14 = (mm42 +mm12     )
          l23 = (mm22    -q23)
          l24 = (mm22    -q24)
          l34 = (        -q34)

          r12 = l12/mm22
          r21 = l12/mm12      
          r24 = l24/mm42
          r42 = l24/mm22      
          r14 = 1d0
          r41 = mm42/mm12
          
          if(dreal(l12).lt.-0d0) then
            ir12 = -1d0
          else
            ir12 = 0d0
          end if
          if(dreal(l24).lt.-0d0) then
            ir24 = -1d0
          else
            ir24 = 0d0
          end if
          ir14=0d0

          y(1)=(l24*l13*l14-mm42*l12*l13-mm12*l24*l34)/
     &        ((l12*l34-l13*l24)*mm42)
          y(2)=(l12*l34-l13*l24)*l24/((l23*l24-mm22*l34)*mm42)
          z(1) = y(1)/r24
          z(2) = y(2)/r24

          d=-l34*l24
          iy(1) = -sign(1d0,dreal(d))
          iy(2) = +sign(1d0,dreal(d))
          iz(1) =  sign(1d0,iy(1)*dreal(r24))
          iz(2) =  sign(1d0,iy(2)*dreal(r24))

c>          write(*,*) 'D0ms1ir0_check l ',l12,l13,l14
c>          write(*,*) 'D0ms1ir0_check l ',l23,l24,l34
c>          write(*,*) 'D0ms1ir0_check r ',r12,r14,r24
c>          write(*,*) 'D0ms1ir0_check r ',r21,r41,r42
c>          write(*,*) 'D0ms1ir0_check y ',y(1),y(2)
c>          write(*,*) 'D0ms1ir0_check z ',z(1),z(2)

          D0ms1ir0_check = dcmplx(0d0)

          eta=eta2s_coli(-y(1),1d0/r24,-iy(1),-ir24,-iz(1))
c>            ch1   = l12-r24*l14-mm22*(r42-1d0/r24)*y(i)
c>            test1 = abs(ch1)/max(abs(l12),abs(r42*y(i)),
c>     &          abs(r24*l14),abs(1/r24*y(i)))
c>            ch2    = mm12/z(i)+l12+mm22*z(i)
c>            test2 = abs(ch2)/max(abs(mm12/z(i)),abs(l12),
c>     &          abs(mm22*z(i)))
c>            ch3   = l23+l13/z(i)
c>            test3 = abs(ch3)/max(abs(l23),abs(l13/z(i)))
c>            test2 = min(test2,test3)
c>            if (test1.lt.test2) then
c>              l1 = cln_coli(ch2/ch3,dreal(-mm22*(r42-1d0/r24)*iy(i)/d))
c>            else
c>              l1 = cln_coli(ch1/d,dreal(-mm22*(r42-1d0/r24)*iy(i)/d))      
c>            end if
          if(eta.ne.0d0)then
            l1(1) = cln_coli((mm12+l12*z(1))/l13,dreal(iy(1)/l34))

c>              write(*,*) 'D0ms1ir0_check l1(1)= ',l1(1)
c>              write(*,*) 'D0ms1ir0_check l2= ',
c>     &            cln_coli(l12-r24*l14-mm42*y(1)*r24,-1d0)
c>     &            -cln_coli(-r24*l34,-1d0),
c>     &            cln_coli(mm12+l14*y(1)+mm42*y(1)*y(1),-1d0)
c>     &            -cln_coli(l13+l34*y(1),-1d0),
c>     &            cln_coli(mm12+l12*z(1),-1d0)-cln_coli(l13,-1d0)

            if(conscheck)then
              if(dimag(r24).ne.0d0) then
                v=dimag(y(1))/dimag(r24)
                u=dimag(y(1)/r24)/dimag(1d0/r24)
                l1(2)=cln_coli(mm12+l12*v+l14*u+mm42*u*u+l24*u*v,-1d0)
     &              -cln_coli(l13+l34*u,-1d0)

c                write(*,*) 'D0ms1ir0_check comp1 l1(2)= ',l1(2),
c     &              l1(1),l1(2)/l1(1)              

              elseif(dimag(y(1)).ne.0d0) then
                l1(2)=cln_coli(l12-r24*l14-mm42*y(1)*r24,-1d0)
     &              -cln_coli(-r24*l34,-1d0)
              elseif(dreal(y(1)).gt.0d0)then
                l1(2)=cln_coli(mm12+l14*y(1)+mm42*y(1)*y(1),-1d0)
     &              -cln_coli(l13+l34*y(1),-1d0)
              elseif(dreal(z(1)).gt.0d0)then
                l1(2)=cln_coli(mm12+l12*z(1),-1d0)-cln_coli(l13,-1d0)
              else
                if (flag3) then
                  write(*,*) 
                  write(*,*) 
     &                'D0ms1ir0_check: case for l1(2) not implemented '
                  flag3 = .false.
                end if
                write(errout,*) 
                write(errout,*) 
     &              'D0ms1ir0_check: case for l1(2) not implemented '
                write(errout,*) 
     &              'D0ms1ir0_check args ',q12,q23,q34,q14,q13,q24,
     &              mm12,mm22,m32,mm42
                write(errout,*) 'r24 = ',r24
                write(errout,*) 'y(1) = ',y(1)
                write(errout,*) 'z(1) = ',z(1)
              endif 
              
              if(abs(l1(2)/l1(1)-1d0).gt.testacc) then
                if (flag2) then
                  write(*,*) 
                  write(*,*) 'D0ms1ir0_check: eta terms possibly wrong '
                  flag2 = .false.
                end if
                write(errout,*) 
                write(errout,*) 
     &              'D0ms1ir0_check: eta terms  possibly wrong 4'
                write(errout,*) 
     &              'D0ms1ir0_check args ',q12,q23,q34,q14,q13,q24,
     &              mm12,mm22,m32,mm42
                write(errout,*) 'l1(1) = ',l1(1)
                write(errout,*) 'l1(2) = ',l1(2)
                write(errout,*) 'diff   = ',abs((l1(2)-l1(1))/l1(1))
              endif 
            endif
          else
            l1(1)=0d0
          endif

          D0ms1ir0_check = D0ms1ir0_check 
     &        - cspcos_coli(-y(1),r41,-iy(1),ir14)
     &        - cspcos_coli(-y(1),cd1,-iy(1),-ir14)
     &        + cspcos_coli(-z(1),r21,-iz(1),ir12)
     &        + pi2_6
     &        + cspcos_coli(-y(1),l34/l13,-iy(1),dreal(l34-l13))
     &        - pi2_6
     &        + eta*(l1(1)+cln_coli(l13/mm12,-1d0)) 

c>          write(*,*) 'D0ms1ir0_check sum ',D0ms1ir0_check
c>
c>          write(*,*) 'D0ms1ir0_check contr 1',
c>     &           - cspcos_coli(-y(1),r41,-iy(1),ir14),
c>     &           - cspcos_coli(-y(1),cd1,-iy(1),-ir14),
c>     &           + cspcos_coli(-z(1),r21,-iz(1),ir12),
c>     &           + pi2_6,
c>     &           + cspcos_coli(-y(1),l34/l13,-iy(1),dreal(l34-l13)),
c>     &           - pi2_6,
c>     &           + eta*(l1(1)+cln_coli(l13/mm12,-1d0))

          eta=eta2s_coli(-y(2),1d0/r24,-iy(2),-ir24,-iz(2))
          if(eta.ne.0d0)then
            l1(1) = cln_coli(mm42*y(2)/l34,dreal(iy(2)/l34))

c>              write(*,*) 'D0ms1ir0_check l1(1)= ',l1(1)
c>              write(*,*) 'D0ms1ir0_check l1(2)= ',
c>     &            cln_coli(-mm42*y(2)*r24,-1d0)
c>     &            -cln_coli(-r24*l34,-1d0),
c>     &            cln_coli(mm42*y(2)*y(2),-1d0)
c>     &            -cln_coli(l34*y(2),-1d0),
c>     &            cln_coli(l12*z(2)+mm22*z(2)*z(2),-1d0)
c>     &            -cln_coli(l13+l23*z(2),-1d0)

            if(conscheck)then
              if(dimag(r24).ne.0d0) then
                v=dimag(y(2))/dimag(r24)
                u=dimag(y(2)/r24)/dimag(1d0/r24)
                l1(2)=cln_coli(mm42*u*u+l24*u*v,-1d0)
     &              -cln_coli(l34*u,-1d0)

c     write(*,*) 'D0ms1ir0_check comp2 l1(2)= ',l1(2),l1(1),l1(2)/l1(1)

              elseif(dimag(y(2)).ne.0d0) then
                l1(2)=cln_coli(-mm42*y(2)*r24,-1d0)
     &              -cln_coli(-r24*l34,-1d0)
              elseif(dreal(y(2)).gt.0d0)then
                l1(2)=cln_coli(mm42*y(2)*y(2),-1d0)
     &              -cln_coli(l34*y(2),-1d0)
              elseif(dreal(z(2)).gt.0d0)then
                l1(2)=cln_coli(l12*z(2)+mm22*z(2)*z(2),-1d0)
     &              -cln_coli(l13+l23*z(2),-1d0)
              else
                if (flag3) then
                  write(*,*) 
                  write(*,*) 
     &                'D0ms1ir0_check: case for l1(2) not implemented '
                  flag3 = .false.
                end if
                write(errout,*) 
                write(errout,*) 
     &              'D0ms1ir0_check: case for l1(2) not implemented '
                write(errout,*) 
     &              'D0ms1ir0_check args ',q12,q23,q34,q14,q13,q24,
     &              mm12,mm22,m32,mm42
                write(errout,*) 'r24 = ',r24
                write(errout,*) 'y(2) = ',y(2)
                write(errout,*) 'z(2) = ',z(2)
              endif 
              
              if(abs(l1(2)/l1(1)-1d0).gt.testacc) then
                if (flag2) then
                  write(*,*) 
                  write(*,*) 
     &                'D0ms1ir0_check: eta terms possibly wrong '
                  flag2 = .false.
                end if
                write(errout,*) 
                write(errout,*) 
     &              'D0ms1ir0_check: eta terms  possibly wrong '
                write(errout,*) 
     &              'D0ms1ir0_check args ',q12,q23,q34,q14,q13,q24,
     &              mm12,mm22,m32,mm42
                write(errout,*) 'l1(1) = ',l1(1)
                write(errout,*) 'l1(2) = ',l1(2)
                write(errout,*) 'diff   = ',abs((l1(2)-l1(1))/l1(1))
              endif 
            endif
          else
            l1(1)=0d0
          endif

c>          write(*,*) 'D0ms1ir0_check eta 2',
c>     &          eta,l1(1),mm42*y(2)/l34
c>          write(*,*) 'D0ms1ir0_check eta 2',
c>     &          -y(2),1d0/r24,-iy(2),-ir24,-iz(2)

          D0ms1ir0_check = D0ms1ir0_check 
c     &           + cspcos_coli(-y(2),r41,-iy(2),ir14)
c     &           + cspcos_coli(-y(2),1d0/r14,-iy(i),-ir14)
     &        -pi2_6-.5d0*
     &        (cln_coli(-y(2),-iy(2))+cln_coli(r41,ir14))**2
     &        -pi2_6-.5d0*(cln_coli(-y(2),-iy(2)))**2
c     &        - cspcos_coli(-z(2),r21,-iz(2),ir12)
     &        +pi2_6+.5d0*
     &        (cln_coli(-z(2),-iz(2))+cln_coli(r21,ir12))**2
     &        - cspcos_coli(-z(2),1d0/r12,-iz(2),-ir12)
c     &           - cspcos_coli(-y(2),l34/l13,-iy(2),dreal(l34-l13))
     &        +pi2_6+.5d0*
     &        (cln_coli(-y(2),-iy(2))
     &        +cln_coli(l34/l13,dreal(l34-l13)))**2
     &        + cspcos_coli(-z(2),l23/l13,-iz(2),dreal(l23-l13))
     &        - eta*(l1(1)+cln_coli(l13/mm12,-1d0)) 

c>          write(*,*) 'D0ms1ir0_check contr 2 ',
c>c     &           + cspcos_coli(-y(2),r41,-iy(2),ir14)
c>c     &           + cspcos_coli(-y(2),1d0/r14,-iy(i),-ir14)
c>     &      -pi2_6-.5d0*(cln_coli(-y(2),-iy(2)))**2,
c>     &      +pi2_6+.5d0*
c>     &      (cln_coli(-y(2),-iy(2))+cln_coli(l34/l13,dreal(l34-l13)))**2,
c>c     &        - cspcos_coli(-z(2),r21,-iz(2),ir12)
c>     &      +pi2_6+.5d0*
c>     &      (cln_coli(-z(2),-iz(2))+cln_coli(r21,ir12))**2,
c>     &      - cspcos_coli(-z(2),1d0/r12,-iz(2),-ir12),
c>c     &           - cspcos_coli(-y(2),l34/l13,-iy(2),dreal(l34-l13))
c>     &        -pi2_6-.5d0*
c>     &      (cln_coli(-y(2),-iy(2))+cln_coli(r41,ir14))**2,
c>     &      + cspcos_coli(-z(2),l23/l13,-iz(2),dreal(l23-l13)),
c>     &      - eta*(l1(1)+cln_coli(l13/mm12,-1d0)) 
c>          write(*,*) 'D0ms1ir0_check sum ',D0ms1ir0_check

          D0ms1ir0_check = D0ms1ir0_check/(l12*l34-l13*l24)

c          write(*,*) 'D0ms1ir0_check =  ',D0ms1ir0_check             

c     end test version
          
          if(cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0).gt.testacc) then
            write(testout,100) ' D0ms1ir0_coli inconsistent e'
            write(testout,111) ' D0ms1ir0_coli: p12 = ',p12
            write(testout,111) ' D0ms1ir0_coli: p23 = ',p23
            write(testout,111) ' D0ms1ir0_coli: p34 = ',p34
            write(testout,111) ' D0ms1ir0_coli: p14 = ',p14
            write(testout,111) ' D0ms1ir0_coli: p24 = ',p24
            write(testout,111) ' D0ms1ir0_coli: p13 = ',p13
            write(testout,111) ' D0ms1ir0_coli: m12 = ',m12
            write(testout,111) ' D0ms1ir0_coli: m22 = ',m22
            write(testout,111) ' D0ms1ir0_coli: m32 = ',m32
            write(testout,111) ' D0ms1ir0_coli: m42 = ',m42
            write(testout,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
            write(testout,111) ' D0ms1ir0_check = ',D0ms1ir0_check
            write(testout,101) ' diff        = ',
     &          cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
            
            write(*,100) ' D0ms1ir0_coli inconsistent e '
            write(*,111) ' D0ms1ir0_coli: p12 = ',p12
            write(*,111) ' D0ms1ir0_coli: p23 = ',p23
            write(*,111) ' D0ms1ir0_coli: p34 = ',p34
            write(*,111) ' D0ms1ir0_coli: p14 = ',p14
            write(*,111) ' D0ms1ir0_coli: p24 = ',p24
            write(*,111) ' D0ms1ir0_coli: p13 = ',p13
            write(*,111) ' D0ms1ir0_coli: m12 = ',m12
            write(*,111) ' D0ms1ir0_coli: m22 = ',m22
            write(*,111) ' D0ms1ir0_coli: m32 = ',m32
            write(*,111) ' D0ms1ir0_coli: m42 = ',m42
            write(*,111) ' D0ms1ir0_coli  = ',D0ms1ir0_coli
            write(*,111) ' D0ms1ir0_check = ',D0ms1ir0_check
            write(*,101) ' diff        = ',
     &          cdabs(D0ms1ir0_coli/D0ms1ir0_check-1d0)
c            stop
          endif
        endif
#endif


      endif
      endif

      end

************************************************************************
      function D0ms1ir1_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     4-point function                                                 *
*     with 1 mass singularities and up to 1 soft singularities         *
*                                                                      *
*     assumes                                                          *
*     p12, m12, m22 small                                              *
*     p13, p24, p14=m42, p23=/=m32, p34 finite                         *
*                                                                      *
*                     m22 small                                        *
*     p12 small  ---------------------  p23=/=m32                      *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | m32                                  *
*                     |         |                                      *
*                     |    4    |                                      *
*       p14=m42  ---------------------  p34                            *
*                         m42                                          *
*                                                                      *
************************************************************************
*     22.08.08 Ansgar Denner        last changed  31.08.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms1ir1_coli

#ifdef CHECK
      complex*16 D0ms1ir1_check
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      logical flag(0:8)
      data    flag /9*.true./
      save    flag  
#endif

      complex*16 m2,m3,m4
      complex*16 k24,k23,k34
      complex*16 l34
      complex*16 r34,r43inv
      real*8     ir34
      complex*16 y,z,logx3,logy,logxs

      complex*16 cln_coli,cspcon_coli,cspcos_coli,cspenc_coli

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
#ifdef CHECK
 101  format(a22,g25.17)
      if(argcheck)then
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

        if(ms12.ne.cd0.or.ms22.ne.cd0
     &        .or.ms42.eq.cd0.or.ps12.ne.cd0.or.ps23.eq.ms32
     &        .or.ps23.eq.cd0.and.ms32.eq.cd0.or.ps14.ne.ms42
     &        .or.ps24.eq.cd0.or.ps13.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &        .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0ms1ir1_coli called improperly:'
          write(*,111)' D0ms1ir1_coli: p12 = ',p12,ps12
          write(*,111)' D0ms1ir1_coli: p23 = ',p23,ps23
          write(*,111)' D0ms1ir1_coli: p34 = ',p34,ps34
          write(*,111)' D0ms1ir1_coli: p14 = ',p14,ps14
          write(*,111)' D0ms1ir1_coli: p13 = ',p13,ps13
          write(*,111)' D0ms1ir1_coli: p24 = ',p24,ps24
          write(*,111)' D0ms1ir1_coli: m12 = ',m12,ms12
          write(*,111)' D0ms1ir1_coli: m22 = ',m22,ms22
          write(*,111)' D0ms1ir1_coli: m32 = ',m32,ms32
          write(*,111)' D0ms1ir1_coli: m42 = ',m42,ms42
          write(*,*)' D0ms1ir1_coli: tests',
     &        ms12.ne.cd0,ms22.ne.cd0
     &        ,ms42.eq.cd0,ps12.ne.cd0,ps23.eq.ms32
     &        ,ps23.eq.cd0.and.ms32.eq.cd0,ps14.ne.ms42
     &        ,ps24.eq.cd0,ps13.eq.cd0
     &        ,dimag(p12).ne.0d0,dimag(p23).ne.0d0
     &        ,dimag(p34).ne.0d0,dimag(p14).ne.0d0
     &        ,dimag(p24).ne.0d0,dimag(p13).ne.0d0
          stop
        end if
      end if
#endif

      if(m12.eq.cd0)then
        if(m22.eq.cd0)then
* 1 soft singularities and 2 zero masses
*
*                   0             
*    p12=0   ----------------  p23    
*                |  2   |         
*              0 |1    3| m32       
*                |   4  |         
*    p14=m42 ----------------  p34
*                   m42             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.0d0.or.p12.ne.0d0.or.
     &        p14.ne.m42.or.p23.eq.m32)
     &        then
            write(*,100)' D0ms1ir1_coli:  (3.17)',
     &          '    case with 1 soft singularity and 2 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms1ir1_coli: p12 = ',p12
            write(*,111)' D0ms1ir1_coli: p23 = ',p23
            write(*,111)' D0ms1ir1_coli: p34 = ',p34
            write(*,111)' D0ms1ir1_coli: p14 = ',p14
            write(*,111)' D0ms1ir1_coli: p13 = ',p13
            write(*,111)' D0ms1ir1_coli: p24 = ',p24
            write(*,111)' D0ms1ir1_coli: m12 = ',m12
            write(*,111)' D0ms1ir1_coli: m22 = ',m22
            write(*,111)' D0ms1ir1_coli: m32 = ',m32
            write(*,111)' D0ms1ir1_coli: m42 = ',m42
            stop
          endif
#endif

#ifdef ALLCHECK          
c ppttbb dim
          if(flag(1))then
            write(*,100)' D0ms1ir1_coli:  (3.17)',
     &          '    case with 1 soft singularity and  2 zero masses:',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms1ir1_coli: p12 = ',p12
            write(*,111)' D0ms1ir1_coli: p23 = ',p23
            write(*,111)' D0ms1ir1_coli: p34 = ',p34
            write(*,111)' D0ms1ir1_coli: p14 = ',p14
            write(*,111)' D0ms1ir1_coli: p13 = ',p13
            write(*,111)' D0ms1ir1_coli: p24 = ',p24
            write(*,111)' D0ms1ir1_coli: m12 = ',m12
            write(*,111)' D0ms1ir1_coli: m22 = ',m22
            write(*,111)' D0ms1ir1_coli: m32 = ',m32
            write(*,111)' D0ms1ir1_coli: m42 = ',m42
            flag(1)=.false.
          endif
#endif

c     according to Ellis, Zanderighi, JHEP 0802 (2008) 002   box12
c     p3^2 -> p34, p4^2 -> p23, m3^2 -> m42, m4^2-> m32
c     s23 -> p13, s12 -> p24
c     Denner Dittmaier (3.17)


          l34 = m32+m42-p34
          if(l34.ne.cd0)then
            r34 = l34/(2d0*m42)*(1d0+sqrt(1d0-4d0*m32*m42/l34**2))
c          r43 = r34*m42/m32
            r43inv= m32/(r34*m42)
            if(m32.ne.cd0)then
              ir34 = sign(1d0,1d0-abs(r34/r43inv))
            else
              ir34 = -1d0
            endif
          else
            r34 = dcmplx(0d0,1d0)*sqrt(m32/m42)
            r43inv= -r34
          endif

c>          if(r43inv.ne.cd0)then
c>          write(*,*) 'D0ms1ir1_coli test1 ',
c>     &        -cspcos_coli((m32-p23)/(m42-p24),1d0/r43inv,
c>     &        dreal(p24-m42-p23+m32),ir34)
c>     &        -.5d0*(cln_coli(r43inv,-ir34))**2
c>     &        +cln_coli((m32-p23)/(m42-p24),dreal(p24-p23))*log(m32/m42)
c>          write(*,*) 'D0ms1ir1_coli test2 ',
c>     &        cspcos_coli((m42-p24)/(m32-p23),r43inv,
c>     &        -dreal(p24-m42-p23+m32),-ir34)
c>     &        +.5d0*(cln_coli((m32-p23)/(m42-p24),
c>     &                dreal(p24-m42-p23+m32)))**2
c>     &        +cln_coli((m32-p23)/(m42-p24),dreal(p24-m42-p23+m32))*
c>     &        cln_coli(r34,ir34)
c>           write(*,*) 'D0ms1ir1_coli test3 ',
c>     &        -cspcos_coli((m32-p23)/(m42-p24),1d0/r43inv,
c>     &        dreal(p24-m42-p23+m32),ir34)
c>     &        -.5d0*(cln_coli(r43inv,-ir34))**2
c>     &        +cln_coli((m32-p23)/(m42-p24),dreal(p24-p23))*log(m32/m42)
c>     &        -cspcos_coli((m32-p23)/(m42-p24),1d0/r34,
c>     &        dreal(p24-m42-p23+m32),-ir34)
c>     &        -.5d0*(cln_coli(r34,ir34))**2
c>          write(*,*) 'D0ms1ir1_coli test4 ',
c>     &        cspcos_coli((m42-p24)/(m32-p23),r43inv,
c>     &        -dreal(p24-m42-p23+m32),-ir34)
c>     &        +(cln_coli((m32-p23)/(m42-p24),
c>     &                dreal(p24-m42-p23+m32)))**2
c>     &        +cspcos_coli((m42-p24)/(m32-p23),r34,
c>     &        -dreal(p24-m42-p23+m32),ir34)
c>         endif

          D0ms1ir1_coli=
#ifdef SING
     &        .5d0*delta2ir
     &        -delta1ir*(cln_coli((m42-p24)/sqrt(m42*muir2),-1d0)
     &        +cln_coli((m32-p13)/(m32-p23),dreal(p23-p13))) 
#endif
     &        -2d0*
     &        cspcon_coli((m32-p23),1d0/(m32-p13),(m32-p23)/(m32-p13),
     &        (p23-p13)/(m32-p13),-1d0,1d0)
     &        +2d0*cln_coli((m42-p24)/sqrt(m42*muir2),-1d0)
     &        *cln_coli((m32-p13)/sqrt(m42*muir2),-1d0)
     &        -(cln_coli((m32-p23)/sqrt(m42*muir2),-1d0))**2
     &        -pi2_6
     &        +cspcos_coli((m42-p24)/(m32-p23),r34,
     &        -dreal(p24-m42-p23+m32),ir34)
     &        +cspcos_coli((m42-p24)/(m32-p23),r43inv,
     &        -dreal(p24-m42-p23+m32),-ir34)
     &        +(cln_coli((m32-p23)/(m42-p24),
     &                dreal(p24-m42-p23+m32)))**2
c     &        -cspcos_coli((m32-p23)/(m42-p24),1d0/r34,
c     &        dreal(p24-m42-p23+m32),-ir34)
c     &        -.5d0*(cln_coli(r34,ir34))**2
c     &        -cspcos_coli((m32-p23)/(m42-p24),r43,
c     &        dreal(p24-m42-p23+m32),ir34)
c     &        -.5d0*(cln_coli(r43,ir34))**2
c     &        +cln_coli((m32-p23)/(m42-p24),dreal(p24-p23))*log(m32/m42)


          D0ms1ir1_coli=D0ms1ir1_coli/((p13-m32)*(p24-m42))
         

        else                    ! m12=0, m22=/=0
          if(m32.ne.cd0)then

* 1 soft singularities and 1 zero masses
*
*                  m22             
*    p12=m22 ----------------  p23    
*                |  2   |         
*              0 |1    3| m32       
*                |   4  |         
*    p14=m42 ----------------  p34
*                   m42             
*

#ifdef CHECK
            if(m12.ne.cd0.or.m22.eq.0d0.or.p12.ne.m22.or.
     &          p14.ne.m42.or.m32.eq.cd0)
     &          then
              write(*,100)' D0ms1ir1_coli:  (3.18) m32=/=0',
     &            '    case with 1 soft singularity and 1 zero mass:',
     &            '    wrong arguments'
              write(*,111)' D0ms1ir1_coli: p12 = ',p12
              write(*,111)' D0ms1ir1_coli: p23 = ',p23
              write(*,111)' D0ms1ir1_coli: p34 = ',p34
              write(*,111)' D0ms1ir1_coli: p14 = ',p14
              write(*,111)' D0ms1ir1_coli: p13 = ',p13
              write(*,111)' D0ms1ir1_coli: p24 = ',p24
              write(*,111)' D0ms1ir1_coli: m12 = ',m12
              write(*,111)' D0ms1ir1_coli: m22 = ',m22
              write(*,111)' D0ms1ir1_coli: m32 = ',m32
              write(*,111)' D0ms1ir1_coli: m42 = ',m42
              stop
            endif
            
            if(flag(2))then
              write(*,100)' D0ms1ir1_coli:  (3.18) m32=/=0 ',
     &            '    case with 1 soft singularity and 1 zero mass:',
     &        '    not yet tested at all'
              write(*,111)' D0ms1ir1_coli: p12 = ',p12
              write(*,111)' D0ms1ir1_coli: p23 = ',p23
              write(*,111)' D0ms1ir1_coli: p34 = ',p34
              write(*,111)' D0ms1ir1_coli: p14 = ',p14
              write(*,111)' D0ms1ir1_coli: p13 = ',p13
              write(*,111)' D0ms1ir1_coli: p24 = ',p24
              write(*,111)' D0ms1ir1_coli: m12 = ',m12
              write(*,111)' D0ms1ir1_coli: m22 = ',m22
              write(*,111)' D0ms1ir1_coli: m32 = ',m32
              write(*,111)' D0ms1ir1_coli: m42 = ',m42
              flag(2)=.false.
            endif
#endif

c     cD0ire
c     according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349  (ia)
c     Ellis, Zanderighi, JHEP 0802 (2008) 002    box16 m2 small
c     Denner, Dittmaier (3.18)

            m2 = sqrt(m22*coliminfscale)
            m3 = sqrt(m32)
            m4 = sqrt(m42)

c     xs = 1/k24   x2 = 1/k23   x3 = r34
            k24 = (m42-p24)/(m2*m4)
            k23 = (m32-p23)/(m2*m3)
            k34 = (m32+m42-p34)/(m3*m4)
            if (k34.ne.0d0) then
              r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2)))
            else
              r34 = dcmplx(0d0,1d0)
            end if
            if (dreal(k34).lt.-2d0) then
              ir34= sign(1d1,1d0-abs(r34))
            else
              ir34= 0d0
            end if
            
            y   = m4*(p23-m32)/(m3*(p24-m42))
            logy = log(m4/m3)+cln_coli(p23-m32,1d0)
     &          -cln_coli(p24-m42,1d0)
            logx3= cln_coli(r34,ir34) 
            D0ms1ir1_coli = (- 2d0*cln_coli(k24,-1d0)
     &          *cln_coli((m32-p13)/(m3*sqrt(muir2)),-1d0)
     &          + (cln_coli(k23,-1d0))**2 + logx3*logx3 + 3d0*pi2_6
#ifdef SING
     &          + cln_coli(k24,-1d0)*delta1ir
     &          + 0.25d0*colishiftms2
#endif
     &          - cspenc_coli(y*r34,1d0) 
     &          - cln_coli(1d0-y*r34,-1d0)*(logy+logx3)
     &          - cspenc_coli(y/r34,1d0) 
     &          - cln_coli(1d0-y/r34,-1d0)*(logy-logx3)
     &          )

              D0ms1ir1_coli =  D0ms1ir1_coli/((p24-m42)*(m32-p13))
            
c     write(*,*) 'D0ms1ir1_coli 2',D0ms1ir1_coli

#ifdef CHECK
            if (conscheck) then

              m2 = sqrt(m22*coliminfscale)
              m3 = sqrt(m32)
              m4 = sqrt(m42)
c     xs = 1/k24   x2 = 1/k23   x3 = r34
              k24 = (m42-p24)/(m2*m4)
              k23 = (m32-p23)/(m2*m3)
              k34 = (m32+m42-p34)/(m3*m4)
              if (k34.ne.0d0) then
                r34 = k34/2d0*(1d0+sqrt(dcmplx(1d0-4d0/k34**2)))
              else
                r34 = dcmplx(0d0,1d0)
              end if
              if(dreal(k34).lt.-2d0) then
                ir34= sign(1d1,1d0-abs(r34))
              else
                ir34= 0d0
              end if
              
              y    = m4*(p23-m32)/(m3*(p24-m42))
              logy = log(m4/m3)+cln_coli(p23-m32,1d0)
     &            -cln_coli(p24-m42,1d0)
              logx3 = cln_coli(r34,ir34) 
              D0ms1ir1_check = (- 2d0*cln_coli(k24,-1d0)
     &            *cln_coli((m32-p13)/(m3*sqrt(muir2)),-1d0)
     &            + (cln_coli(k23,-1d0))**2 + logx3**2 + 3d0*pi2_6
#ifdef SING
     &          + cln_coli(k24,-1d0)*delta1ir
     &          + 0.25d0*colishiftms2
#endif
     &            - cspenc_coli(y*r34,1d0) 
     &            - cln_coli(1d0-y*r34,-1d0)*(logy+logx3)
     &            - cspenc_coli(y/r34,1d0) 
     &            - cln_coli(1d0-y/r34,-1d0)*(logy-logx3)
     &            )

              D0ms1ir1_check =  D0ms1ir1_check/((p24-m42)*(m32-p13))

c     write(*,*) 'D0ms1ir1_coli 1',D0ms1ir1_check

              if(cdabs(D0ms1ir1_coli/D0ms1ir1_check-1d0).gt.testacc)
     &            then
                write(testout,*) ' D0ms1ir1_coli inconsistent '
                write(testout,100)' D0ms1ir1_coli:  (3.18) m22=/=0',
     &            '    case with 1 soft singularity and 1 zero mass:'
                write(testout,111) ' D0ms1ir1_coli: p12 = ',p12
                write(testout,111) ' D0ms1ir1_coli: p23 = ',p23
                write(testout,111) ' D0ms1ir1_coli: p34 = ',p34
                write(testout,111) ' D0ms1ir1_coli: p14 = ',p14
                write(testout,111) ' D0ms1ir1_coli: p24 = ',p24
                write(testout,111) ' D0ms1ir1_coli: p13 = ',p13
                write(testout,111) ' D0ms1ir1_coli: m12 = ',m12
                write(testout,111) ' D0ms1ir1_coli: m22 = ',m22
                write(testout,111) ' D0ms1ir1_coli: m32 = ',m32
                write(testout,111) ' D0ms1ir1_coli: m42 = ',m42
                write(testout,111) ' D0ms1ir1_coli  = ',D0ms1ir1_coli
                write(testout,111) ' D0ms1ir1_check = ',D0ms1ir1_check
                write(testout,101) ' diff        = ',
     &              cdabs(D0ms1ir1_coli/D0ms1ir1_check-1d0)
                
                write(*,*)   ' D0ms1ir1_coli inconsistent '
                write(*,100)' D0ms1ir1_coli:  (3.18) m22=/=0',
     &            '    case with 1 soft singularity and 1 zero mass:'
                write(*,111) ' D0ms1ir1_coli: p12 = ',p12
                write(*,111) ' D0ms1ir1_coli: p23 = ',p23
                write(*,111) ' D0ms1ir1_coli: p34 = ',p34
                write(*,111) ' D0ms1ir1_coli: p14 = ',p14
                write(*,111) ' D0ms1ir1_coli: p24 = ',p24
                write(*,111) ' D0ms1ir1_coli: p13 = ',p13
                write(*,111) ' D0ms1ir1_coli: m12 = ',m12
                write(*,111) ' D0ms1ir1_coli: m22 = ',m22
                write(*,111) ' D0ms1ir1_coli: m32 = ',m32
                write(*,111) ' D0ms1ir1_coli: m42 = ',m42
                write(*,111) ' D0ms1ir1_coli  = ',D0ms1ir1_coli
                write(*,111) ' D0ms1ir1_check = ',D0ms1ir1_check
                write(*,101) ' diff        = ',
     &              cdabs(D0ms1ir1_coli/D0ms1ir1_check-1d0)
                stop
              endif
            endif
#endif
          
          else                  ! m32=0

* 1 soft singularities and 2 zero masses
*
*                  m22             
*    p12=m22 ----------------  p23    
*                |  2   |         
*              0 |1    3| 0       
*                |   4  |         
*    p14=m42 ----------------  p34
*                   m42             
*

#ifdef CHECK
            if(m12.ne.cd0.or.m22.eq.0d0.or.p12.ne.m22.or.
     &          p14.ne.m42.or.m32.ne.cd0.or.
     &          m22.eq.p23.and.m42.eq.p34)
     &          then
              write(*,100)' D0ms1ir1_coli:  (3.18) m32=0',
     &            '    case with 1 soft singularity and 2 zero masses:',
     &            '    wrong arguments'
              write(*,111)' D0ms1ir1_coli: p12 = ',p12
              write(*,111)' D0ms1ir1_coli: p23 = ',p23
              write(*,111)' D0ms1ir1_coli: p34 = ',p34
              write(*,111)' D0ms1ir1_coli: p14 = ',p14
              write(*,111)' D0ms1ir1_coli: p13 = ',p13
              write(*,111)' D0ms1ir1_coli: p24 = ',p24
              write(*,111)' D0ms1ir1_coli: m12 = ',m12
              write(*,111)' D0ms1ir1_coli: m22 = ',m22
              write(*,111)' D0ms1ir1_coli: m32 = ',m32
              write(*,111)' D0ms1ir1_coli: m42 = ',m42
            endif
            
#endif

c     cD0ir0e
c     according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349  (iia)
c     Ellis, Zanderighi, JHEP 0802 (2008) 002    box16 m2 small, m3=0
c     Denner, Dittmaier (3.18) with m2=0

            m2 = sqrt(m22*coliminfscale2)
            m4 = sqrt(m42)
            k24 = (m42-p24)/(m2*m4)
            logxs = -cln_coli(k24,-1d0)
c     xs = 1/k24   

            if (p34.ne.m42) then

#ifdef CHECK
              if(flag(3))then
                write(*,100)' D0ms1ir1_coli: (3.18) m32=0, p34=/=m42',
     &              '    case with 1 soft singularity ',
     &              '    and 2 zero masses:',
     &              '    not yet tested at all'
                write(*,111)' D0ms1ir1_coli: p12 = ',p12
                write(*,111)' D0ms1ir1_coli: p23 = ',p23
                write(*,111)' D0ms1ir1_coli: p34 = ',p34
                write(*,111)' D0ms1ir1_coli: p14 = ',p14
                write(*,111)' D0ms1ir1_coli: p13 = ',p13
                write(*,111)' D0ms1ir1_coli: p24 = ',p24
                write(*,111)' D0ms1ir1_coli: m12 = ',m12
                write(*,111)' D0ms1ir1_coli: m22 = ',m22
                write(*,111)' D0ms1ir1_coli: m32 = ',m32
                write(*,111)' D0ms1ir1_coli: m42 = ',m42
                flag(3)=.false.
c                stop
              endif
#endif

              y    = m2*(p34-m42)/(m4*p23)
              z    = 1d0/(k24*y)
              logy = log(m2/m4)+cln_coli(p34-m42,1d0)
     &            -cln_coli(p23,1d0)
              D0ms1ir1_coli = ( logxs*( -logxs/2d0 
     &            - log(muir2/(m2*m4))
     &            +cln_coli(p13/p23,dreal(p23))
     &            +cln_coli(p13/(p34-m42),dreal(p34-m42)) )
     &            + logy**2/2d0 + pi2_6 
#ifdef SING
     &            - logxs*delta1ir
     &            + 0.25d0*colishiftms2
#endif
     &            - cspenc_coli(z,1d0) 
     &            - cln_coli(1d0-z,-1d0)*(logxs-logy) )

              D0ms1ir1_coli =  D0ms1ir1_coli/((p24-m42)*(-p13))

            else

#ifdef CHECK
              if(flag(4))then
                write(*,100)' D0ms1ir1_coli:  (3.18) m32=0, p34=m42 ',
     &              '    case with 1 soft singularity and 1 zero mass:',
     &              '    not yet tested at all'
                write(*,111)' D0ms1ir1_coli: p12 = ',p12
                write(*,111)' D0ms1ir1_coli: p23 = ',p23
                write(*,111)' D0ms1ir1_coli: p34 = ',p34
                write(*,111)' D0ms1ir1_coli: p14 = ',p14
                write(*,111)' D0ms1ir1_coli: p13 = ',p13
                write(*,111)' D0ms1ir1_coli: p24 = ',p24
                write(*,111)' D0ms1ir1_coli: m12 = ',m12
                write(*,111)' D0ms1ir1_coli: m22 = ',m22
                write(*,111)' D0ms1ir1_coli: m32 = ',m32
                write(*,111)' D0ms1ir1_coli: m42 = ',m42
                flag(4)=.false.
              endif
#endif

c     according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 (iic??)
c     Ellis, Zanderighi, JHEP 0802 (2008) 002    box16 m2 small, m3=0
c     Denner, Dittmaier (3.18) 

c        write(*,*) 'D0ms1ir1_coli input',p13,p23,p24,m22,m42
c        write(*,*) 'D0ms1ir1_coli pars ',muir2,coliminfscale2
c        write(*,*) 'D0ms1ir1_coli logxs',logxs


              D0ms1ir1_coli = ( logxs*(
     &            - logxs - log(muir2/(m22*coliminfscale2))
     &            + 2d0*cln_coli(p13/p23,dreal(p23)) )
     &            - pi2_6 )

              D0ms1ir1_coli =  D0ms1ir1_coli/((p24-m42)*(-p13))

            end if        

          endif                 !(m32.ne.cd0)
        endif                   !(m22.eq.cd0)
      elseif(m12.ne.cd0)then
        if(m22.eq.cd0)then      ! m12=/=0, m22=0
* 0 soft singularities and 1 zero masses
*
*                   0             
*    p12=m12 ----------------  p23    
*                |  2   |         
*            m12 |1    3| m32       
*                |   4  |         
*    p14=m42 ----------------  p34
*                   m42             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.ne.0d0.or.p12.ne.m12.or.
     &        p23.eq.m32)
     &        then
            write(*,100)' D0ms1ir1_coli:  (3.19)',
     &          '    case with 0 soft singularity and 1 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms1ir1_coli: p12 = ',p12
            write(*,111)' D0ms1ir1_coli: p23 = ',p23
            write(*,111)' D0ms1ir1_coli: p34 = ',p34
            write(*,111)' D0ms1ir1_coli: p14 = ',p14
            write(*,111)' D0ms1ir1_coli: p13 = ',p13
            write(*,111)' D0ms1ir1_coli: p24 = ',p24
            write(*,111)' D0ms1ir1_coli: m12 = ',m12
            write(*,111)' D0ms1ir1_coli: m22 = ',m22
            write(*,111)' D0ms1ir1_coli: m32 = ',m32
            write(*,111)' D0ms1ir1_coli: m42 = ',m42
            stop
          endif
          
          if(flag(5))then
            write(*,100)' D0ms1ir1_coli:  (3.19)',
     &          '    case with 0 soft singularity and 1 zero masses:',
c     &        '    not yet tested at all'
     &          '    not yet implemented'
            write(*,111)' D0ms1ir1_coli: p12 = ',p12
            write(*,111)' D0ms1ir1_coli: p23 = ',p23
            write(*,111)' D0ms1ir1_coli: p34 = ',p34
            write(*,111)' D0ms1ir1_coli: p14 = ',p14
            write(*,111)' D0ms1ir1_coli: p13 = ',p13
            write(*,111)' D0ms1ir1_coli: p24 = ',p24
            write(*,111)' D0ms1ir1_coli: m12 = ',m12
            write(*,111)' D0ms1ir1_coli: m22 = ',m22
            write(*,111)' D0ms1ir1_coli: m32 = ',m32
            write(*,111)' D0ms1ir0_coli: m42 = ',m42
            flag(5)=.false.
            stop
          endif
#endif

c     Denner, Dittmaier (3.19)

          D0ms1ir1_coli=cd0
      
        else                    ! m12=/=0, m22=/=0
*     1 soft singularities and 1 zero masses
*     
*                 m22=m12             
*     p12=0   ----------------  p23    
*                |  2   |         
*            m12 |1    3| m32       
*                |   4  |         
*     p14=m42 ----------------  p34
*                  m42             
*     

#ifdef CHECK
          if(m12.ne.m22.or.m22.eq.0d0.or.p12.ne.cd0)
     &        then
            write(*,100)' D0ms1ir1_coli:  (3.20)',
     &          '    case with 0 soft singularity and 0 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms1ir1_coli: p12 = ',p12
            write(*,111)' D0ms1ir1_coli: p23 = ',p23
            write(*,111)' D0ms1ir1_coli: p34 = ',p34
            write(*,111)' D0ms1ir1_coli: p14 = ',p14
            write(*,111)' D0ms1ir1_coli: p13 = ',p13
            write(*,111)' D0ms1ir1_coli: p24 = ',p24
            write(*,111)' D0ms1ir1_coli: m12 = ',m12
            write(*,111)' D0ms1ir1_coli: m22 = ',m22
            write(*,111)' D0ms1ir1_coli: m32 = ',m32
            write(*,111)' D0ms1ir1_coli: m42 = ',m42
            stop
          endif
          
          if(flag(6))then
            write(*,100)' D0ms1ir1_coli:  (3.20)',
     &          '    case with 0 soft singularity and 0 zero masses:',
c     &        '    not yet tested at all'
     &          '    not yet implemented'
            write(*,111)' D0ms1ir1_coli: p12 = ',p12
            write(*,111)' D0ms1ir1_coli: p23 = ',p23
            write(*,111)' D0ms1ir1_coli: p34 = ',p34
            write(*,111)' D0ms1ir1_coli: p14 = ',p14
            write(*,111)' D0ms1ir1_coli: p13 = ',p13
            write(*,111)' D0ms1ir1_coli: p24 = ',p24
            write(*,111)' D0ms1ir1_coli: m12 = ',m12
            write(*,111)' D0ms1ir1_coli: m22 = ',m22
            write(*,111)' D0ms1ir1_coli: m32 = ',m32
            write(*,111)' D0ms1ir0_coli: m42 = ',m42
            flag(6)=.false.
            stop
          endif
#endif

c     Denner, Dittmaier (3.20)

          D0ms1ir1_coli=cd0
          

        endif
      endif                     ! m22.eq.cd0


      end

************************************************************************
      function D0ms1ir2_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     4-point function                                                 *
*     with 1 mass singularities and up to 2 soft singularities         *
*                                                                      *
*     assumes                                                          *
*     p12, m12, m22 small                                              *
*     p13, p24, p14=m42, p23=m32, p34 finite                           *
*                                                                      *
*                     m22 small                                        *
*     p12 small  ---------------------  p23=m32                        *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | m32                                  *
*                     |         |                                      *
*                     |    4    |                                      *
*       p14=m42  ---------------------  p34                            *
*                         m42                                          *
*                                                                      *
************************************************************************
*     22.08.08 Ansgar Denner        last changed  09.12.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms1ir2_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      logical flag(0:3)
      data    flag /4*.true./
      save    flag  
#endif

      complex*16 cln_coli
      complex*16 m2(4)
      complex*16 l34,r34
      real*8     ir34
c      complex*16 cln_coli,cspcon_coli,cspenc_coli
      integer    nsoft
      integer    i,j,k
      logical    soft(4),onsh(4,4)

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
#ifdef CHECK
 101  format(a22,g25.17)
      if(argcheck)then
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

        if(ms12.ne.cd0.or.ms22.ne.cd0.or.ms32.eq.cd0
     &        .or.ms42.eq.cd0.or.ps12.ne.cd0.or.ps23.ne.ms32
     &        .or.ps14.ne.ms42
     &        .or.ps24.eq.cd0.or.ps13.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &        .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0ms1ir2_coli called improperly:'
          write(*,111)' D0ms1ir2_coli: p12 = ',p12,ps12
          write(*,111)' D0ms1ir2_coli: p23 = ',p23,ps23
          write(*,111)' D0ms1ir2_coli: p34 = ',p34,ps34
          write(*,111)' D0ms1ir2_coli: p14 = ',p14,ps14
          write(*,111)' D0ms1ir2_coli: p13 = ',p13,ps13
          write(*,111)' D0ms1ir2_coli: p24 = ',p24,ps24
          write(*,111)' D0ms1ir2_coli: m12 = ',m12,ms12
          write(*,111)' D0ms1ir2_coli: m22 = ',m22,ms22
          write(*,111)' D0ms1ir2_coli: m32 = ',m32,ms32
          write(*,111)' D0ms1ir2_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif

      m2(1)=m12
      m2(2)=m22
      m2(3)=m32
      m2(4)=m42

c determine on-shell momenta
      onsh(1,2)=p12.eq.m22
      onsh(1,3)=p13.eq.m32
      onsh(1,4)=p14.eq.m42
      onsh(2,1)=p12.eq.m12
      onsh(2,3)=p23.eq.m32
      onsh(2,4)=p24.eq.m42
      onsh(3,1)=p13.eq.m12
      onsh(3,2)=p23.eq.m22
      onsh(3,4)=p34.eq.m42
      onsh(4,1)=p14.eq.m22
      onsh(4,2)=p24.eq.m22
      onsh(4,3)=p34.eq.m32

c count/determine soft singularities
      nsoft=0 
      do i=1,4
        do j=1,3
          if(i.ne.j)then
            do k=j,4
              if(k.ne.i.and.k.ne.j)then
                soft(i)=m2(i).eq.cd0.and.onsh(i,j).and.onsh(i,k)
c                soft(i,k,j)=soft(i,j,k)
                if(soft(i)) nsoft=nsoft+1
              endif
            enddo
          endif
        enddo
      enddo

c>c      if(flag(0))then
c>        write(*,100)' D0ms1ir2_coli: not yet implemented '
c>        write(*,111)' D0ms1ir2_coli: p12 = ',p12
c>        write(*,111)' D0ms1ir2_coli: p23 = ',p23
c>        write(*,111)' D0ms1ir2_coli: p34 = ',p34
c>        write(*,111)' D0ms1ir2_coli: p14 = ',p14
c>        write(*,111)' D0ms1ir2_coli: p24 = ',p24
c>        write(*,111)' D0ms1ir2_coli: p13 = ',p13
c>        write(*,111)' D0ms1ir2_coli: m12 = ',m12
c>        write(*,111)' D0ms1ir2_coli: m22 = ',m22
c>        write(*,111)' D0ms1ir2_coli: m32 = ',m32
c>        write(*,111)' D0ms1ir2_coli: m42 = ',m42
c>c        flag(0)=.false.
c>        stop
c>c      endif

      if(m12.eq.cd0.and.m22.eq.cd0)then
* 2 soft singularities and 2 zero masses
*
*                   0             
*    p12=0   ----------------  p23=m32    
*                |  2   |         
*              0 |1    3| m32       
*                |   4  |         
*    p14=m42 ----------------  p34
*                   m42             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.0d0.or.p12.ne.0d0.or.
     &        p14.ne.m42.or.p23.ne.m32)
     &        then
            write(*,100)' D0ms1ir2_coli:  (3.22)',
     &          '    case with 2 soft singularities and 2 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms1ir1_coli: p12 = ',p12
            write(*,111)' D0ms1ir2_coli: p23 = ',p23
            write(*,111)' D0ms1ir2_coli: p34 = ',p34
            write(*,111)' D0ms1ir2_coli: p14 = ',p14
            write(*,111)' D0ms1ir2_coli: p13 = ',p13
            write(*,111)' D0ms1ir2_coli: p24 = ',p24
            write(*,111)' D0ms1ir2_coli: m12 = ',m12
            write(*,111)' D0ms1ir2_coli: m22 = ',m22
            write(*,111)' D0ms1ir2_coli: m32 = ',m32
            write(*,111)' D0ms1ir2_coli: m42 = ',m42
            stop
          endif
#endif
      
#ifdef ALLCHECK    
c ppttbb dim
          if(flag(1))then
            write(*,100)' D0ms1ir2_coli:  (3.22)',
     &          '    case with 2 soft singularities and 2 zero masses:',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms1ir2_coli: p12 = ',p12
            write(*,111)' D0ms1ir2_coli: p23 = ',p23
            write(*,111)' D0ms1ir2_coli: p34 = ',p34
            write(*,111)' D0ms1ir2_coli: p14 = ',p14
            write(*,111)' D0ms1ir2_coli: p13 = ',p13
            write(*,111)' D0ms1ir2_coli: p24 = ',p24
            write(*,111)' D0ms1ir2_coli: m12 = ',m12
            write(*,111)' D0ms1ir2_coli: m22 = ',m22
            write(*,111)' D0ms1ir2_coli: m32 = ',m32
            write(*,111)' D0ms1ir2_coli: m42 = ',m42
            flag(1)=.false.
          endif
#endif

c     according to Ellis, Zanderighi, JHEP 0802 (2008) 002   box11
c     p3^2 -> p34, m3^2 -> m32, m4^2-> m42
c     s12 -> p13, s23 -> p24
c     Denner Dittmaier (3.22)

          l34 = m32+m42-p34
          r34 = l34/(2d0*m42)*(1d0+sqrt(1d0-4d0*m32*m42/l34**2))
          r34= sqrt(m42/m32)*r34
          ir34 = sign(1d0,1d0-abs(r34*r34))

          D0ms1ir2_coli =
#ifdef SING
     &        delta2ir
     &        -delta1ir*
     &        (cln_coli((m32-p13)/sqrt(muir2*m32),-1d0)
     &        +cln_coli((m42-p24)/sqrt(muir2*m42),-1d0))
#endif
     &        +2d0*cln_coli((m32-p13)/sqrt(muir2*m32),-1d0)
     &        *cln_coli((m42-p24)/sqrt(muir2*m42),-1d0)
     &        -(cln_coli(r34,-1d0))**2
     &        -4d0*pi2_6
 
          D0ms1ir2_coli = D0ms1ir2_coli/((p13-m32)*(p24-m42))

         

        elseif(m12.ne.cd0.or.m22.ne.cd0)then

* 2 soft singularities and 1 zero masses
*
*                  m22             
*    p12=m22 ----------------  p23=m32    
*                |  2   |         
*              0 |1    3| m32       
*                |   4  |         
*    p14=m42 ----------------  p34
*                   m42             
*

#ifdef CHECK
            if(m12+m22.eq.cd0.or.m22*m12.ne.0d0.or.p12.ne.(m22+m12).or.
     &          p14.ne.m42.or.p23.ne.m32)
     &          then
              write(*,100)' D0ms1ir2_coli:  (3.23)',
     &            '    case with 2 soft singularities and 1 zero mass:',
     &            '    wrong arguments'
              write(*,111)' D0ms1ir2_coli: p12 = ',p12
              write(*,111)' D0ms1ir2_coli: p23 = ',p23
              write(*,111)' D0ms1ir2_coli: p34 = ',p34
              write(*,111)' D0ms1ir2_coli: p14 = ',p14
              write(*,111)' D0ms1ir2_coli: p13 = ',p13
              write(*,111)' D0ms1ir2_coli: p24 = ',p24
              write(*,111)' D0ms1ir2_coli: m12 = ',m12
              write(*,111)' D0ms1ir2_coli: m22 = ',m22
              write(*,111)' D0ms1ir2_coli: m32 = ',m32
              write(*,111)' D0ms1ir2_coli: m42 = ',m42
              stop
            endif
            
            if(flag(2))then
              write(*,100)' D0ms1ir2_coli:  (3.23)',
     &            '    case with 2 soft singularities and 1 zero mass:',
c     &            '    not yet tested at all'
     &            '    not yet implemented'
              write(*,111)' D0ms1ir2_coli: p12 = ',p12
              write(*,111)' D0ms1ir2_coli: p23 = ',p23
              write(*,111)' D0ms1ir2_coli: p34 = ',p34
              write(*,111)' D0ms1ir2_coli: p14 = ',p14
              write(*,111)' D0ms1ir2_coli: p13 = ',p13
              write(*,111)' D0ms1ir2_coli: p24 = ',p24
              write(*,111)' D0ms1ir2_coli: m12 = ',m12
              write(*,111)' D0ms1ir2_coli: m22 = ',m22
              write(*,111)' D0ms1ir2_coli: m32 = ',m32
              write(*,111)' D0ms1ir2_coli: m42 = ',m42
              flag(2)=.false.
            endif
#endif

c     Ellis, Zanderighi, JHEP 0802 (2008) 002    box??
c new case
c     Denner, Dittmaier (3.23)


            D0ms1ir2_coli =  cd0

            D0ms1ir2_coli =  D0ms1ir2_coli/((p24-m42)*(m32-p13))
            
c     write(*,*) 'D0ms1ir2_coli 2',D0ms1ir2_coli


        else

* 2 soft singularities and 0 zero masses
*
*                 m22=m12             
*      p12=0 ----------------  p23=m32    
*                |  2   |         
*      m12 small |1    3| m32       
*                |   4  |         
*    p14=m42 ----------------  p34
*                   m42             
*

#ifdef CHECK
            if(m12.ne.m22.or.m22.eq.0d0.or.p12.ne.cd0.or.
     &          p14.ne.m42.or.p23.ne.m32)
     &          then
              write(*,100)' D0ms1ir2_coli:  (3.24)',
     &            '    case with 2 soft singularities and 0 zero mass:',
     &            '    wrong arguments'
              write(*,111)' D0ms1ir2_coli: p12 = ',p12
              write(*,111)' D0ms1ir2_coli: p23 = ',p23
              write(*,111)' D0ms1ir2_coli: p34 = ',p34
              write(*,111)' D0ms1ir2_coli: p14 = ',p14
              write(*,111)' D0ms1ir2_coli: p13 = ',p13
              write(*,111)' D0ms1ir2_coli: p24 = ',p24
              write(*,111)' D0ms1ir2_coli: m12 = ',m12
              write(*,111)' D0ms1ir2_coli: m22 = ',m22
              write(*,111)' D0ms1ir2_coli: m32 = ',m32
              write(*,111)' D0ms1ir2_coli: m42 = ',m42
              stop
            endif
            
            if(flag(3))then
              write(*,100)' D0ms1ir2_coli:  (3.23)',
     &            '    case with 2 soft singularities and 0 zero mass:',
c     &            '    not yet tested at all'
     &            '    not yet implemented'
              write(*,111)' D0ms1ir2_coli: p12 = ',p12
              write(*,111)' D0ms1ir2_coli: p23 = ',p23
              write(*,111)' D0ms1ir2_coli: p34 = ',p34
              write(*,111)' D0ms1ir2_coli: p14 = ',p14
              write(*,111)' D0ms1ir2_coli: p13 = ',p13
              write(*,111)' D0ms1ir2_coli: p24 = ',p24
              write(*,111)' D0ms1ir2_coli: m12 = ',m12
              write(*,111)' D0ms1ir2_coli: m22 = ',m22
              write(*,111)' D0ms1ir2_coli: m32 = ',m32
              write(*,111)' D0ms1ir2_coli: m42 = ',m42
              flag(3)=.false.
            endif
#endif

c     Ellis, Zanderighi, JHEP 0802 (2008) 002    box??
c new case
c     Denner, Dittmaier (3.24)


            D0ms1ir2_coli =  cd0

            D0ms1ir2_coli =  D0ms1ir2_coli/((p24-m42)*(m32-p13))
            
c     write(*,*) 'D0ms1ir2_coli 2',D0ms1ir2_coli

      endif

      end

************************************************************************
      function D0ms2ir0_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     4-point function                                                 *
*     with 2 mass singularities and up to 0 soft singularities         *
*                                                                      *
*     assumes                                                          *
*     p12, p34, m12, m22, m32, m42 small                               *
*     p13, p24, p14, p23 finite                                        *
*                                                                      *
*                      m22 small                                       *
*     p12 small  ---------------------  p23                            *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | m32 small                            *
*                     |         |                                      *
*                     |    4    |                                      *
*           p14  ---------------------  p34 small                      *
*                      m42 small                                       *
*                                                                      *
************************************************************************
*     22.08.08 Ansgar Denner        last changed  16.12.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms2ir0_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      logical flag(0:7)
      data    flag /8*.true./
      save    flag  
#endif

      complex*16 m2(4),mm12,q12,q34
      complex*16 cln_coli,cspenc_coli,cspcon_coli,cspcos_coli
      integer    nm0
      integer    i
#ifdef CHECK
      integer    nsoft
      integer    j,k
      logical    soft(4),onsh(4,4)
#endif

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

 100  format(((a)))
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
#ifdef CHECK
 101  format(a22,g25.17)
      if(argcheck)then
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

        if(ms12.ne.cd0.or.ms22.ne.cd0.or.ms32.ne.cd0
     &      .or.ms42.ne.cd0.or.ps12.ne.cd0.or.ps34.ne.cd0 
     &      .or.ps23.eq.cd0.or.ps14.eq.cd0 
     &      .or.ps24.eq.cd0.or.ps13.eq.cd0
     &      .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &      .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &      .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0ms2ir0_coli called improperly:'
          write(*,111)' D0ms2ir0_coli: p12 = ',p12,ps12
          write(*,111)' D0ms2ir0_coli: p23 = ',p23,ps23
          write(*,111)' D0ms2ir0_coli: p34 = ',p34,ps34
          write(*,111)' D0ms2ir0_coli: p14 = ',p14,ps14
          write(*,111)' D0ms2ir0_coli: p13 = ',p13,ps13
          write(*,111)' D0ms2ir0_coli: p24 = ',p24,ps24
          write(*,111)' D0ms2ir0_coli: m12 = ',m12,ms12
          write(*,111)' D0ms2ir0_coli: m22 = ',m22,ms22
          write(*,111)' D0ms2ir0_coli: m32 = ',m32,ms32
          write(*,111)' D0ms2ir0_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif

      m2(1)=m12
      m2(2)=m22
      m2(3)=m32
      m2(4)=m42

#ifdef CHECK
c determine on-shell momenta
      onsh(1,2)=p12.eq.m22
      onsh(1,3)=p13.eq.m32
      onsh(1,4)=p14.eq.m42
      onsh(2,1)=p12.eq.m12
      onsh(2,3)=p23.eq.m32
      onsh(2,4)=p24.eq.m42
      onsh(3,1)=p13.eq.m12
      onsh(3,2)=p23.eq.m22
      onsh(3,4)=p34.eq.m42
      onsh(4,1)=p14.eq.m22
      onsh(4,2)=p24.eq.m22
      onsh(4,3)=p34.eq.m32

c count/determine soft singularities
      nsoft=0 
      do i=1,4
        do j=1,3
          if(i.ne.j)then
            do k=j,4
              if(k.ne.i.and.k.ne.j)then
                soft(i)=m2(i).eq.cd0.and.onsh(i,j).and.onsh(i,k)
c                soft(i,k,j)=soft(i,j,k)
                if(soft(i)) nsoft=nsoft+1
              endif
            enddo
          endif
        enddo
      enddo
#endif

c count/determine zero masses
      nm0=0
      do i=1,4
        if(m2(i).eq.cd0) nm0=nm0+1
      enddo

c>c      if(flag(0))then
c>        write(*,100)' D0ms2ir0_coli: not yet implemented '
c>        write(*,111)' D0ms2ir0_coli: p12 = ',p12
c>        write(*,111)' D0ms2ir0_coli: p23 = ',p23
c>        write(*,111)' D0ms2ir0_coli: p34 = ',p34
c>        write(*,111)' D0ms2ir0_coli: p14 = ',p14
c>        write(*,111)' D0ms2ir0_coli: p24 = ',p24
c>        write(*,111)' D0ms2ir0_coli: p13 = ',p13
c>        write(*,111)' D0ms2ir0_coli: m12 = ',m12
c>        write(*,111)' D0ms2ir0_coli: m22 = ',m22
c>        write(*,111)' D0ms2ir0_coli: m32 = ',m32
c>        write(*,111)' D0ms2ir0_coli: m42 = ',m42
c>c        flag(0)=.false.
c>        stop
c>c      endif

#ifdef CHECK
      if(nsoft.eq.0)then
#endif

        if(nm0.eq.4)then
* 0 soft singularities and 4 zero masses
*
*                   0             
*          0 ----------------  p23    
*                |  2   |         
*              0 |1    3| 0       
*                |   4  |         
*        p14 ----------------  0
*                   0             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.0d0.or.m32.ne.cd0.or.m42.ne.0d0
     &        .or.p12.ne.cd0.or.p34.ne.cd0.or.
     &        p14.eq.cd0.or.p23.eq.cd0)
     &        then
            write(*,100)' D0ms2ir0_coli:  (3.49)',
     &          '    case with 0 soft singularities',
     &          '    and 4 zero masses: wrong arguments'
            write(*,111)' D0ms2ir0_coli: p12 = ',p12
            write(*,111)' D0ms2ir0_coli: p23 = ',p23
            write(*,111)' D0ms2ir0_coli: p34 = ',p34
            write(*,111)' D0ms2ir0_coli: p14 = ',p14
            write(*,111)' D0ms2ir0_coli: p24 = ',p24
            write(*,111)' D0ms2ir0_coli: p13 = ',p13
            write(*,111)' D0ms2ir0_coli: m12 = ',m12
            write(*,111)' D0ms2ir0_coli: m22 = ',m22
            write(*,111)' D0ms2ir0_coli: m32 = ',m32
            write(*,111)' D0ms2ir0_coli: m42 = ',m42
            stop
          endif
#endif
      
#ifdef ALLCHECK
c ppttbb dim 
          if(flag(1))then
            write(*,100)' D0ms2ir0_coli:  (3.49)',
     &          '    case with 0 soft singularities',
     &          '    and 4 zero masses: ',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms2ir0_coli: p12 = ',p12
            write(*,111)' D0ms2ir0_coli: p23 = ',p23
            write(*,111)' D0ms2ir0_coli: p34 = ',p34
            write(*,111)' D0ms2ir0_coli: p14 = ',p14
            write(*,111)' D0ms2ir0_coli: p24 = ',p24
            write(*,111)' D0ms2ir0_coli: p13 = ',p13
            write(*,111)' D0ms2ir0_coli: m12 = ',m12
            write(*,111)' D0ms2ir0_coli: m22 = ',m22
            write(*,111)' D0ms2ir0_coli: m32 = ',m32
            write(*,111)' D0ms2ir0_coli: m42 = ',m42
            flag(1)=.false.
c            stop
          endif
#endif

c according to Ellis, Zanderighi, JHEP 0802 (2008) 002   box3
c p_4^2 -> p14, p_2^2 -> p23, s12 -> p13, s23 -> p24
c and Berger, Klasen, Tait, PRD62 (2000) 095014
c Denner Dittmaier (3.49)

          D0ms2ir0_coli =
#ifdef SING
     &        2d0*delta1ir*(cln_coli(-p14/muir2,-1d0)
     &        +cln_coli(-p23/muir2,-1d0)
     &        -cln_coli(-p24/muir2,-1d0)
     &        -cln_coli(-p13/muir2,-1d0))
#endif
     &        -(cln_coli(-p14/muir2,-1d0))**2
     &        -(cln_coli(-p23/muir2,-1d0))**2
     &        +(cln_coli(-p24/muir2,-1d0))**2
     &        +(cln_coli(-p13/muir2,-1d0))**2
     &        -2d0*cspenc_coli(1d0-p23/p13,dreal(p23-p13))
     &        -2d0*cspenc_coli(1d0-p23/p24,dreal(p23-p24))
     &        -2d0*cspenc_coli(1d0-p14/p13,dreal(p14-p13))
     &        -2d0*cspenc_coli(1d0-p14/p24,dreal(p14-p24))
     &        +2d0*cspcos_coli(p23/p13,p14/p24,
     &        dreal(p13-p23),dreal(p24-p14))
     &        -(cln_coli(p13/p24,p24-p13))**2
 
          D0ms2ir0_coli = D0ms2ir0_coli/(p13*p24-p14*p23)

        elseif(nm0.eq.3)then
* 0 soft singularities and 3 zero masses
*
*    small          0            
*    p12=m12 ----------------  p23    
*                |  2   |         
*      m12 small |1    3| 0       
*                |   4  |         
*        p14 ----------------  0
*                   0             
*

#ifdef CHECK
          if(.not.(m12.eq.cd0.and.m32.eq.cd0.and.m42.eq.0d0.and.
     &        m22.ne.cd0.and.p12.eq.m22).and.
     &        .not.(m22.eq.cd0.and.m32.eq.cd0.and.m42.eq.0d0.and.
     &        m12.ne.cd0.and.p12.eq.m12).and.
     &        .not.(m12.eq.cd0.and.m22.eq.cd0.and.m42.eq.0d0.and.
     &        m32.ne.cd0.and.p34.eq.m32).and.
     &        .not.(m12.eq.cd0.and.m22.eq.cd0.and.m32.eq.0d0.and.
     &        m42.ne.cd0.and.p34.eq.m42).or.
     &        p14.eq.cd0.or.p23.eq.cd0)
     &        then
            write(*,100)' D0ms2ir0_coli:  (3.50)',
     &          '    case with 0 soft singularities',
     &          '    and 3 zero masses: wrong arguments'
            write(*,111)' D0ms2ir0_coli: p12 = ',p12
            write(*,111)' D0ms2ir0_coli: p23 = ',p23
            write(*,111)' D0ms2ir0_coli: p34 = ',p34
            write(*,111)' D0ms2ir0_coli: p14 = ',p14
            write(*,111)' D0ms2ir0_coli: p24 = ',p24
            write(*,111)' D0ms2ir0_coli: p13 = ',p13
            write(*,111)' D0ms2ir0_coli: m12 = ',m12
            write(*,111)' D0ms2ir0_coli: m22 = ',m22
            write(*,111)' D0ms2ir0_coli: m32 = ',m32
            write(*,111)' D0ms2ir0_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(2))then
            write(*,100)' D0ms2ir0_coli:  (3.50)',
     &          '    case with 0 soft singularities and 3 zero masses:',
c     &          '    not yet tested at all'
     &          '    not yet implemented'
            write(*,111)' D0ms2ir0_coli: p12 = ',p12
            write(*,111)' D0ms2ir0_coli: p23 = ',p23
            write(*,111)' D0ms2ir0_coli: p34 = ',p34
            write(*,111)' D0ms2ir0_coli: p14 = ',p14
            write(*,111)' D0ms2ir0_coli: p24 = ',p24
            write(*,111)' D0ms2ir0_coli: p13 = ',p13
            write(*,111)' D0ms2ir0_coli: m12 = ',m12
            write(*,111)' D0ms2ir0_coli: m22 = ',m22
            write(*,111)' D0ms2ir0_coli: m32 = ',m32
            write(*,111)' D0ms2ir0_coli: m42 = ',m42
            flag(2)=.false.
            stop
          endif
#endif

          if(m22.eq.cd0)then
            mm12=m12
            q12=p12
            q34=p34
          elseif(m12.eq.cd0)then
            mm12=m22
            q12=p12
            q34=p34
          elseif(m42.eq.cd0)then
            mm12=m32
            q12=p34
            q34=p12
          else
            mm12=m42
            q12=p34
            q34=p12
          endif

c according to Ellis, Zanderighi, JHEP 0802 (2008) 002   box10
c and Berger, Klasen, Tait, PRD62 (2000) 095014
c Denner Dittmaier (3.50)
c not yet implemented, following formula wrong!!!

          D0ms2ir0_coli =
#ifdef SING
     &        delta2ir
     &        -delta1ir*(cln_coli((-p24)/muir2,-1d0)
     &        +cln_coli((m32-p13)/(m32-p34),-1d0)
     &        -cln_coli((m32-p24)/(m32-p23),-1d0)
     &        +log(m42/muir2))
     &        -.5d0*colishiftms2 
#endif
     &        +2d0*cln_coli(-p24/muir2,-1d0)
     &        *cln_coli((m32-p13)/m32,-1d0)
     &        -cln_coli((m32-p23)/muir2,-1d0)
     &        *cln_coli((m32-p23)/m32,-1d0)
     &        -cln_coli((m32-p34)/muir2,-1d0)
     &        *cln_coli((m32-p34)/m32,-1d0)
     &        +.5d0*(cln_coli(-p24/muir2,-1d0))**2
     &        -.5d0*(cln_coli(-p24/m32,-1d0))**2
     &        -2d0*(
     &        cspcon_coli(m32-p34,1d0/(m32-p13),(m32-p34)/(m32-p13),
     &        (p34-p13)/(m32-p13),-1d0,1d0)
     &        +cspcon_coli(m32-p23,1d0/(m32-p13),(m32-p23)/(m32-p13),
     &        (p23-p13)/(m32-p13),-1d0,1d0) )
     &        -cspcon_coli((m32-p34)/m32,-(m32-p23)/p24,
     &        (m32-p34)*(m32-p23)/(-m32*p24),
     &        ((m32-p34)*(m32-p23)+-m32*p24)/(-m32*p24),
     &        -1d0,dreal(p24-p23))
     &        -0.25d0*pi2 
 
          D0ms2ir0_coli = D0ms2ir0_coli/((p13-m32)*p24)

        elseif(nm0.eq.2)then
          if(m12.eq.m22.and.m42.eq.m32)then

* 0 soft singularities and 2 zero masses on collinear leg
*
*                 m22=m12            
*          0  ----------------  p23    
*                |  2   |         
*      m12 small |1    3| 0       
*                |   4  |         
*        p14 ----------------  0
*                   0             
*

#ifdef CHECK
            if(p12.ne.cd0.or.p34.ne.cd0.or.
     &          m12.ne.cd0.and.m32.ne.cd0.or.
     &          m12.eq.cd0.and.m32.eq.cd0.or.
     &          p14.eq.cd0.or.p23.eq.cd0)
     &          then
              write(*,100)' D0ms2ir0_coli:  (3.51)',
     &            '    case with 0 soft singularities',
     &            '    and 2 zero masses on collinear leg: ',
     &            '    wrong arguments'
              write(*,111)' D0ms2ir0_coli: p12 = ',p12
              write(*,111)' D0ms2ir0_coli: p23 = ',p23
              write(*,111)' D0ms2ir0_coli: p34 = ',p34
              write(*,111)' D0ms2ir0_coli: p14 = ',p14
              write(*,111)' D0ms2ir0_coli: p24 = ',p24
              write(*,111)' D0ms2ir0_coli: p13 = ',p13
              write(*,111)' D0ms2ir0_coli: m12 = ',m12
              write(*,111)' D0ms2ir0_coli: m22 = ',m22
              write(*,111)' D0ms2ir0_coli: m32 = ',m32
              write(*,111)' D0ms2ir0_coli: m42 = ',m42
              stop
            endif
#endif
            
#ifdef CHECK
            if(flag(3))then
              write(*,100)' D0ms2ir0_coli: (3.51)',
     &            '    case with 0 soft singularities',
     &            '    and 2 zero masses on collinear leg',
     &            '    not yet implemented'
              write(*,111)' D0ms2ir0_coli: p12 = ',p12
              write(*,111)' D0ms2ir0_coli: p23 = ',p23
              write(*,111)' D0ms2ir0_coli: p34 = ',p34
              write(*,111)' D0ms2ir0_coli: p14 = ',p14
              write(*,111)' D0ms2ir0_coli: p24 = ',p24
              write(*,111)' D0ms2ir0_coli: p13 = ',p13
              write(*,111)' D0ms2ir0_coli: m12 = ',m12
              write(*,111)' D0ms2ir0_coli: m22 = ',m22
              write(*,111)' D0ms2ir0_coli: m32 = ',m32
              write(*,111)' D0ms2ir0_coli: m42 = ',m42
              flag(3)=.false.
              stop
            endif
#endif            

            if(m32.eq.cd0)then
              mm12=m12
              q12=p12
              q34=p34
            else
              mm12=m32
              q12=p34
              q34=p12
            endif

c according to Ellis, Zanderighi, JHEP 0802 (2008) 002   
c           box13 with p32=0 and m3=m4 small
c Denner Dittmaier (3.51)
c not yet implemented, following formula wrong!!!

            D0ms2ir0_coli =
#ifdef SING
     &          delta2ir
     &          -delta1ir*(cln_coli((-p24)/muir2,-1d0)
     &          +cln_coli((m32-p13)/(m32-p34),-1d0)
     &          -cln_coli((m32-p24)/(m32-p23),-1d0)
     &          +log(m42/muir2))
     &          -.5d0*colishiftms2 
#endif
     &          +2d0*cln_coli(-p24/muir2,-1d0)
     &          *cln_coli((m32-p13)/m32,-1d0)
     &          -cln_coli((m32-p23)/muir2,-1d0)
     &          *cln_coli((m32-p23)/m32,-1d0)
     &          -cln_coli((m32-p34)/muir2,-1d0)
     &          *cln_coli((m32-p34)/m32,-1d0)
     &          +.5d0*(cln_coli(-p24/muir2,-1d0))**2
     &          -.5d0*(cln_coli(-p24/m32,-1d0))**2
     &          -2d0*(
     &          cspcon_coli(m32-p34,1d0/(m32-p13),(m32-p34)/(m32-p13),
     &          (p34-p13)/(m32-p13),-1d0,1d0)
     &          +cspcon_coli(m32-p23,1d0/(m32-p13),(m32-p23)/(m32-p13),
     &          (p23-p13)/(m32-p13),-1d0,1d0) )
     &          -cspcon_coli((m32-p34)/m32,-(m32-p23)/p24,
     &          (m32-p34)*(m32-p23)/(-m32*p24),
     &          ((m32-p34)*(m32-p23)+-m32*p24)/(-m32*p24),
     &          -1d0,dreal(p24-p23))
     &          -0.25d0*pi2 
            
            D0ms2ir0_coli = D0ms2ir0_coli/((p13-m32)*p24)
            
          elseif(m12.eq.m32.and.m22.eq.m42)then

* 0 soft singularities and 2 opposite zero masses  
*
*    small          0            
*    p12=m12 ----------------  p23    
*                |  2   |         
*      m12 small |1    3| m32 small       
*                |   4  |         
*        p14 ----------------  p34=m32
*                   0          small
*

#ifdef CHECK
            if(p12.eq.cd0.or.p34.eq.cd0.or.
     &          m12.ne.cd0.and.m22.ne.cd0.or.
     &          m12.eq.cd0.and.m22.eq.cd0.or.
     &          p14.eq.cd0.or.p23.eq.cd0)
     &          then
              write(*,100)' D0ms2ir0_coli:  (3.52)',
     &            '    case with 0 soft singularities',
     &            '    and 2 zero masses on opposite lines: ',
     &            '    wrong arguments'
              write(*,111)' D0ms2ir0_coli: p12 = ',p12
              write(*,111)' D0ms2ir0_coli: p23 = ',p23
              write(*,111)' D0ms2ir0_coli: p34 = ',p34
              write(*,111)' D0ms2ir0_coli: p14 = ',p14
              write(*,111)' D0ms2ir0_coli: p24 = ',p24
              write(*,111)' D0ms2ir0_coli: p13 = ',p13
              write(*,111)' D0ms2ir0_coli: m12 = ',m12
              write(*,111)' D0ms2ir0_coli: m22 = ',m22
              write(*,111)' D0ms2ir0_coli: m32 = ',m32
              write(*,111)' D0ms2ir0_coli: m42 = ',m42
              stop
            endif
#endif
            
#ifdef CHECK
            if(flag(4))then
              write(*,100)' D0ms2ir0_coli: (3.52)',
     &            '    case with 0 soft singularities',
     &            '    and 2 zero masses on opposite lines:',
     &            '    not yet implemented'
              write(*,111)' D0ms2ir0_coli: p12 = ',p12
              write(*,111)' D0ms2ir0_coli: p23 = ',p23
              write(*,111)' D0ms2ir0_coli: p34 = ',p34
              write(*,111)' D0ms2ir0_coli: p14 = ',p14
              write(*,111)' D0ms2ir0_coli: p24 = ',p24
              write(*,111)' D0ms2ir0_coli: p13 = ',p13
              write(*,111)' D0ms2ir0_coli: m12 = ',m12
              write(*,111)' D0ms2ir0_coli: m22 = ',m22
              write(*,111)' D0ms2ir0_coli: m32 = ',m32
              write(*,111)' D0ms2ir0_coli: m42 = ',m42
              flag(4)=.false.
              stop
            endif
#endif
            
            if(m22.eq.cd0)then
              mm12=m12
              q12=p12
              q34=p34
            else
              mm12=m32
              q12=p34
              q34=p12
            endif
            
c     Denner Dittmaier (3.51)
c     new case
            D0ms2ir0_coli = cd0
#ifdef SING
#endif
 
            D0ms2ir0_coli = D0ms2ir0_coli/((p13-m32)*p24)
            
          elseif(m12.eq.m42.and.m22.eq.m32)then

* 0 soft singularities and 2 zero masses on non-collinear leg  
*
*    small          0            
*    p12=m12 ----------------  p23    
*                |  2   |         
*      m12 small |1    3| 0
*                |   4  |         
*        p14 ----------------  p34=m42
*                m42 small       small
*

#ifdef CHECK
            if(p12.eq.cd0.or.p34.eq.cd0.or.
     &          m12.ne.cd0.and.m22.ne.cd0.or.
     &          m12.eq.cd0.and.m22.eq.cd0.or.
     &          p14.eq.cd0.or.p23.eq.cd0)
     &          then
              write(*,100)' D0ms2ir0_coli:  (3.53)',
     &            '    case with 0 soft singularities',
     &            '    and 2 zero masses on noncollinear leg: ',
     &            ' wrong arguments'
              write(*,111)' D0ms2ir0_coli: p12 = ',p12
              write(*,111)' D0ms2ir0_coli: p23 = ',p23
              write(*,111)' D0ms2ir0_coli: p34 = ',p34
              write(*,111)' D0ms2ir0_coli: p14 = ',p14
              write(*,111)' D0ms2ir0_coli: p24 = ',p24
              write(*,111)' D0ms2ir0_coli: p13 = ',p13
              write(*,111)' D0ms2ir0_coli: m12 = ',m12
              write(*,111)' D0ms2ir0_coli: m22 = ',m22
              write(*,111)' D0ms2ir0_coli: m32 = ',m32
              write(*,111)' D0ms2ir0_coli: m42 = ',m42
              stop
            endif
#endif
            
#ifdef CHECK
            if(flag(5))then
              write(*,100)' D0ms2ir0_coli: (3.53)',
     &            '    case with 0 soft singularities',
     &            '    and 2 zero masses on non-collinear leg',
     &            '    not yet implemented'
              write(*,111)' D0ms2ir0_coli: p12 = ',p12
              write(*,111)' D0ms2ir0_coli: p23 = ',p23
              write(*,111)' D0ms2ir0_coli: p34 = ',p34
              write(*,111)' D0ms2ir0_coli: p14 = ',p14
              write(*,111)' D0ms2ir0_coli: p24 = ',p24
              write(*,111)' D0ms2ir0_coli: p13 = ',p13
              write(*,111)' D0ms2ir0_coli: m12 = ',m12
              write(*,111)' D0ms2ir0_coli: m22 = ',m22
              write(*,111)' D0ms2ir0_coli: m32 = ',m32
              write(*,111)' D0ms2ir0_coli: m42 = ',m42
              flag(5)=.false.
              stop
            endif
#endif
            
            if(m22.eq.cd0)then
              mm12=m12
              q12=p12
              q34=p34
            else
              mm12=m42
              q12=p34
              q34=p12
            endif
            
c     Denner Dittmaier (3.53)
c     new case
          D0ms2ir0_coli = cd0
#ifdef SING
#endif
          
          D0ms2ir0_coli = D0ms2ir0_coli/((p13-m32)*p24)
          
          
#ifdef CHECK
        else
          write(*,100)' D0ms2ir0_coli: case with more than 0 soft 
     &        singularities and 2 zero masses: wrong arguments '
          write(*,111)' D0ms2ir0_coli: nm0 = ',nm0
          write(*,111)' D0ms2ir0_coli: p23 = ',p23
          write(*,111)' D0ms2ir0_coli: p34 = ',p34
          write(*,111)' D0ms2ir0_coli: p14 = ',p14
          write(*,111)' D0ms2ir0_coli: p24 = ',p24
          write(*,111)' D0ms2ir0_coli: p13 = ',p13
          write(*,111)' D0ms2ir0_coli: m12 = ',m12
          write(*,111)' D0ms2ir0_coli: m22 = ',m22
          write(*,111)' D0ms2ir0_coli: m32 = ',m32
          write(*,111)' D0ms2ir0_coli: m42 = ',m42    
          stop
#endif
        endif
        
      elseif(nm0.eq.1)then
* 0 soft singularities and 1 zero masses
*
*    small       m22 small            
*    p12=m22 ----------------  p23    
*                |  2   |         
*              0 |1    3| m32 small       
*                |   4  |         
*        p14 ----------------  0
*                 m42=m32             
*

#ifdef CHECK
          if(.not.(m12.ne.cd0.and.m32.ne.cd0.and.m42.ne.0d0.and.
     &        m22.eq.cd0.and.p12.eq.m12.and.p34.eq.cd0).and.
     &        .not.(m22.ne.cd0.and.m32.ne.cd0.and.m42.ne.0d0.and.
     &        m12.eq.cd0.and.p12.eq.m22.and.p34.eq.cd0).and.
     &        .not.(m12.ne.cd0.and.m22.ne.cd0.and.m42.ne.0d0.and.
     &        m32.eq.cd0.and.p34.eq.m42.and.p12.eq.cd0).and.
     &        .not.(m12.ne.cd0.and.m22.ne.cd0.and.m32.ne.0d0.and.
     &        m42.eq.cd0.and.p34.eq.m32.and.p12.eq.cd0).or.
     &        p14.eq.cd0.or.p23.eq.cd0)
     &        then
            write(*,100)' D0ms2ir0_coli:  (3.54)',
     &          '    case with 0 soft singularities',
     &          '    and 1 zero masses: wrong arguments'
            write(*,111)' D0ms2ir0_coli: p12 = ',p12
            write(*,111)' D0ms2ir0_coli: p23 = ',p23
            write(*,111)' D0ms2ir0_coli: p34 = ',p34
            write(*,111)' D0ms2ir0_coli: p14 = ',p14
            write(*,111)' D0ms2ir0_coli: p24 = ',p24
            write(*,111)' D0ms2ir0_coli: p13 = ',p13
            write(*,111)' D0ms2ir0_coli: m12 = ',m12
            write(*,111)' D0ms2ir0_coli: m22 = ',m22
            write(*,111)' D0ms2ir0_coli: m32 = ',m32
            write(*,111)' D0ms2ir0_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(6))then
            write(*,100)' D0ms2ir0_coli:  (3.54)',
     &          '    case with 0 soft singularities and 1 zero masses:',
c     &          '     not yet tested at all'
     &          '    not yet implemented'
            write(*,111)' D0ms2ir0_coli: p12 = ',p12
            write(*,111)' D0ms2ir0_coli: p23 = ',p23
            write(*,111)' D0ms2ir0_coli: p34 = ',p34
            write(*,111)' D0ms2ir0_coli: p14 = ',p14
            write(*,111)' D0ms2ir0_coli: p24 = ',p24
            write(*,111)' D0ms2ir0_coli: p13 = ',p13
            write(*,111)' D0ms2ir0_coli: m12 = ',m12
            write(*,111)' D0ms2ir0_coli: m22 = ',m22
            write(*,111)' D0ms2ir0_coli: m32 = ',m32
            write(*,111)' D0ms2ir0_coli: m42 = ',m42
            flag(6)=.false.
            stop
          endif
#endif

          if(m12.ne.cd0)then
            mm12=m12
            q12=p12
            q34=p34
          elseif(m22.ne.cd0)then
            mm12=m22
            q12=p12
            q34=p34
          elseif(m32.ne.cd0)then
            mm12=m32
            q12=p34
            q34=p12
          else
            mm12=m42
            q12=p34
            q34=p12
          endif

c Denner Dittmaier (3.54)
c new case
c not yet implemented, following formula wrong!!!


          D0ms2ir0_coli =
#ifdef SING
     &        delta2ir
     &        -delta1ir*(cln_coli((-p24)/muir2,-1d0)
     &        +cln_coli((m32-p13)/(m32-p34),-1d0)
     &        -cln_coli((m32-p24)/(m32-p23),-1d0)
     &        +log(m42/muir2))
     &        -.5d0*colishiftms2 
#endif
     &        +2d0*cln_coli(-p24/muir2,-1d0)
     &        *cln_coli((m32-p13)/m32,-1d0)
     &        -cln_coli((m32-p23)/muir2,-1d0)
     &        *cln_coli((m32-p23)/m32,-1d0)
     &        -cln_coli((m32-p34)/muir2,-1d0)
     &        *cln_coli((m32-p34)/m32,-1d0)
     &        +.5d0*(cln_coli(-p24/muir2,-1d0))**2
     &        -.5d0*(cln_coli(-p24/m32,-1d0))**2
     &        -2d0*(
     &        cspcon_coli(m32-p34,1d0/(m32-p13),(m32-p34)/(m32-p13),
     &        (p34-p13)/(m32-p13),-1d0,1d0)
     &        +cspcon_coli(m32-p23,1d0/(m32-p13),(m32-p23)/(m32-p13),
     &        (p23-p13)/(m32-p13),-1d0,1d0) )
     &        -cspcon_coli((m32-p34)/m32,-(m32-p23)/p24,
     &        (m32-p34)*(m32-p23)/(-m32*p24),
     &        ((m32-p34)*(m32-p23)+-m32*p24)/(-m32*p24),
     &        -1d0,dreal(p24-p23))
     &        -0.25d0*pi2 
 
          D0ms2ir0_coli = D0ms2ir0_coli/((p13-m32)*p24)






        elseif(nm0.eq.0)then
* 0 soft singularities and 0 zero masses
*
*                 m22=m12             
*          0 ----------------  p23    
*                |  2   |         
*       m12 small|1    3| m32 small       
*                |   4  |         
*        p14 ----------------  0
*                 m42=m32             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.eq.0d0.or.m32.eq.cd0.or.m42.eq.0d0
     &        .or.p12.ne.cd0.or.p34.ne.cd0.or.
     &        p14.eq.cd0.or.p23.eq.cd0)
     &        then
            write(*,100)' D0ms2ir0_coli: (3.55)',
     &          '    case with 0 soft singularities',
     &          '    and 0 zero masses: wrong arguments'
            write(*,111)' D0ms2ir0_coli: p12 = ',p12
            write(*,111)' D0ms2ir0_coli: p23 = ',p23
            write(*,111)' D0ms2ir0_coli: p34 = ',p34
            write(*,111)' D0ms2ir0_coli: p14 = ',p14
            write(*,111)' D0ms2ir0_coli: p24 = ',p24
            write(*,111)' D0ms2ir0_coli: p13 = ',p13
            write(*,111)' D0ms2ir0_coli: m12 = ',m12
            write(*,111)' D0ms2ir0_coli: m22 = ',m22
            write(*,111)' D0ms2ir0_coli: m32 = ',m32
            write(*,111)' D0ms2ir0_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(7))then
            write(*,100)' D0ms2ir0_coli:  (3.55)',
     &          '    case with 0 soft singularities and 0 zero masses:',
c     &          '    not yet tested at all'
     &          '    not yet implemented'
            write(*,111)' D0ms2ir0_coli: p12 = ',p12
            write(*,111)' D0ms2ir0_coli: p23 = ',p23
            write(*,111)' D0ms2ir0_coli: p34 = ',p34
            write(*,111)' D0ms2ir0_coli: p14 = ',p14
            write(*,111)' D0ms2ir0_coli: p24 = ',p24
            write(*,111)' D0ms2ir0_coli: p13 = ',p13
            write(*,111)' D0ms2ir0_coli: m12 = ',m12
            write(*,111)' D0ms2ir0_coli: m22 = ',m22
            write(*,111)' D0ms2ir0_coli: m32 = ',m32
            write(*,111)' D0ms2ir0_coli: m42 = ',m42
            flag(7)=.false.
            stop
          endif
#endif

c Denner Dittmaier (3.55)
c new case
          D0ms2ir0_coli =cd0
#ifdef SING
#endif
          D0ms2ir0_coli = D0ms2ir0_coli/((p13-m32)*p24)

        endif

#ifdef CHECK
      else
        write(*,100)' D0ms2ir0_coli: case with more than 0 soft 
     &      singularities not implemented '
        write(*,111)' D0ms2ir0_coli: nsoft = ',nsoft
        write(*,111)' D0ms2ir0_coli: p23 = ',p23
        write(*,111)' D0ms2ir0_coli: p34 = ',p34
        write(*,111)' D0ms2ir0_coli: p14 = ',p14
        write(*,111)' D0ms2ir0_coli: p24 = ',p24
        write(*,111)' D0ms2ir0_coli: p13 = ',p13
        write(*,111)' D0ms2ir0_coli: m12 = ',m12
        write(*,111)' D0ms2ir0_coli: m22 = ',m22
        write(*,111)' D0ms2ir0_coli: m32 = ',m32
        write(*,111)' D0ms2ir0_coli: m42 = ',m42    
        stop

      endif
#endif

      end

************************************************************************
      function D0ms2ir1_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     4-point function                                                 *
*     with 2 mass singularities and up to 1 soft singularities         *
*                                                                      *
*     assumes                                                          *
*     p12, p14, m12, m22, m42 small                                    *
*     p13, p24, p23, p34, m32 finite                                   *
*                                                                      *
*                     m22 small                                        *
*     p12 small  ---------------------  p23                            *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | m32                                  *
*                     |         |                                      *
*                     |    4    |                                      *
*     p14 small  ---------------------  p34                            *
*                       m42 small                                      *
*                                                                      *
************************************************************************
*     27.08.08 Ansgar Denner        last changed  16.12.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms2ir1_coli

#ifdef CHECK
      complex*16 elimminf2_coli,D0ms2ir1_check
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      logical flag(0:7),flag3,flag2
      data    flag /8*.true./,flag3 /.true./, flag2 /.true./
      save    flag,flag2,flag3
#endif

      complex*16 q12,q23,q34,q14,q13,q24,mm12,mm22,mm32,mm42
      complex*16 l12,l23,l34,l14,l13,l24
      complex*16 r12,r14,r24,r21,r41,r42
      complex*16 d
      complex*16 m2(4)
      complex*16 y(2),z(2),eta,l1,l2,x1,x2
      complex*16 logxs
      complex*16 cln_coli,cspcon_coli,cspenc_coli,cspcos_coli,eta2s_coli
      real*8     ir12,ir24,ir14,iy(2),iz(2),u,v,ix1,ix2
      integer    nsoft
      integer    i,j,k
      logical    soft(4),onsh(4,4)

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"


c      write(*,*) ' D0ms2ir1_coli in ',p12,p23,p34,p14,p13,p24,
c     &    m12,m22,m32,m42

#ifdef CHECK
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if(argcheck)then
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

        if(ms12.ne.cd0.or.ms22.ne.cd0.or.ms42.ne.cd0
     &        .or.ps12.ne.cd0.or.ps14.ne.cd0 
     &        .or.ps34.eq.ms32.or.ps23.eq.ms32
     &        .or.ps24.eq.cd0.or.ps13.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &        .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0ms2ir1_coli called improperly:'
          write(*,111)' D0ms2ir1_coli: p12 = ',p12,ps12
          write(*,111)' D0ms2ir1_coli: p23 = ',p23,ps23
          write(*,111)' D0ms2ir1_coli: p34 = ',p34,ps34
          write(*,111)' D0ms2ir1_coli: p14 = ',p14,ps14
          write(*,111)' D0ms2ir1_coli: p13 = ',p13,ps13
          write(*,111)' D0ms2ir1_coli: p24 = ',p24,ps24
          write(*,111)' D0ms2ir1_coli: m12 = ',m12,ms12
          write(*,111)' D0ms2ir1_coli: m22 = ',m22,ms22
          write(*,111)' D0ms2ir1_coli: m32 = ',m32,ms32
          write(*,111)' D0ms2ir1_coli: m42 = ',m42,ms42
          write(*,*)' D0ms2ir1_coli: test  ',
     &        ms12.ne.cd0,ms22.ne.cd0,ms42.ne.cd0
     &        ,ps12.ne.cd0,ps14.ne.cd0 
     &        ,ps34.eq.ms32,ps23.eq.ms32
     &        ,ps24.eq.cd0,ps13.eq.cd0
     &        ,dimag(p12).ne.0d0,dimag(p23).ne.0d0
     &        ,dimag(p34).ne.0d0,dimag(p14).ne.0d0
     &        ,dimag(p24).ne.0d0,dimag(p13).ne.0d0
          stop
        end if
      end if
#endif

      m2(1)=m12
      m2(2)=m22
      m2(3)=m32
      m2(4)=m42

c determine on-shell momenta
      onsh(1,2)=p12.eq.m22
      onsh(1,3)=p13.eq.m32
      onsh(1,4)=p14.eq.m42
      onsh(2,1)=p12.eq.m12
      onsh(2,3)=p23.eq.m32
      onsh(2,4)=p24.eq.m42
      onsh(3,1)=p13.eq.m12
      onsh(3,2)=p23.eq.m22
      onsh(3,4)=p34.eq.m42
      onsh(4,1)=p14.eq.m22
      onsh(4,2)=p24.eq.m22
      onsh(4,3)=p34.eq.m32

c count/determine soft singularities
      nsoft=0 
      do i=1,4
        do j=1,3
          if(i.ne.j)then
            do k=j,4
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

        if(m22.eq.cd0.and.m42.eq.cd0)then
          if(m12.eq.cd0)then
* 1 soft singularities and 3 zero masses
*     
*                   0             
*          0 ----------------  p23    
*                |  2   |         
*              0 |1    3| m32       
*                |   4  |         
*          0 ----------------  p34
*                   0             
*

#ifdef CHECK
            if(m12.ne.cd0.or.m22.ne.0d0.or.m42.ne.0d0
     &          .or.p12.ne.cd0.or.p14.ne.cd0.or.
     &          p34.eq.m32.or.p23.eq.m32)
     &          then
              write(*,100)' D0ms2ir1_coli:  (3.26)',
     &            '    case with 1 soft singularities',
     &            '    and 3 zero masses: wrong arguments'
              write(*,111)' D0ms2ir1_coli: p12 = ',p12
              write(*,111)' D0ms2ir1_coli: p23 = ',p23
              write(*,111)' D0ms2ir1_coli: p34 = ',p34
              write(*,111)' D0ms2ir1_coli: p14 = ',p14
              write(*,111)' D0ms2ir1_coli: p24 = ',p24
              write(*,111)' D0ms2ir1_coli: p13 = ',p13
              write(*,111)' D0ms2ir1_coli: m12 = ',m12
              write(*,111)' D0ms2ir1_coli: m22 = ',m22
              write(*,111)' D0ms2ir1_coli: m32 = ',m32
              write(*,111)' D0ms2ir1_coli: m42 = ',m42
              stop
            endif
#endif
            
#ifdef ALLCHECK
c ppttbb dim
            if(flag(1))then
              write(*,100)' D0ms2ir1_coli:  (3.26)',
     &            '    case with 1 soft singularities',
     &            '    and 3 zero masses ',
     &            '    not yet tested in physical process'
              write(*,111)' D0ms2ir1_coli: p12 = ',p12
              write(*,111)' D0ms2ir1_coli: p23 = ',p23
              write(*,111)' D0ms2ir1_coli: p34 = ',p34
              write(*,111)' D0ms2ir1_coli: p14 = ',p14
              write(*,111)' D0ms2ir1_coli: p24 = ',p24
              write(*,111)' D0ms2ir1_coli: p13 = ',p13
              write(*,111)' D0ms2ir1_coli: m12 = ',m12
              write(*,111)' D0ms2ir1_coli: m22 = ',m22
              write(*,111)' D0ms2ir1_coli: m32 = ',m32
              write(*,111)' D0ms2ir1_coli: m42 = ',m42
              flag(1)=.false.
            endif
#endif

c     according to Ellis, Zanderighi, JHEP 0802 (2008) 002   box8
c p_3^2 -> p23, p_4^2 -> p34, m_4^2 -> m32, s12 -> p24, s23 -> p13
c     and Berger, Klasen, Tait, PRD62 (2000) 095014
c     Denner Dittmaier (3.26)
            D0ms2ir1_coli =
#ifdef SING
     &          delta2ir
     &          -delta1ir*(cln_coli((-p24)/muir2,-1d0)
     &          +cln_coli((m32-p13)/(m32-p34),dreal(p34-p13))
     &          +cln_coli((m32-p13)/(m32-p23),dreal(p23-p13)))
c?     &          -.5d0*colishiftms2 
#endif
c changed to allow for zero m32
c>     &          +2d0*cln_coli(-p24/muir2,-1d0)
c>     &          *cln_coli((m32-p13)/m32,-1d0)
c>     &          -cln_coli((m32-p23)/muir2,-1d0)
c>     &          *cln_coli((m32-p23)/m32,-1d0)
c>     &          -cln_coli((m32-p34)/muir2,-1d0)
c>     &          *cln_coli((m32-p34)/m32,-1d0)
c>     &          -.5d0*(cln_coli(-p24/m32,-1d0))**2
c>     &          -cspcon_coli((m32-p34)/m32,-(m32-p23)/p24,
c>     &          (m32-p34)*(m32-p23)/(-m32*p24),
c>     &          ((m32-p34)*(m32-p23)+m32*p24)/(m32*p24),
c>     &          -1d0,dreal(p24-p23))
     &          +cln_coli(-p24/muir2,-1d0)
     &          *(2d0*cln_coli((m32-p13)/(-p24),dreal(p24-p13+m32))
     &          -cln_coli((m32-p23)/(-p24),dreal(p24-p23+m32))
     &          -cln_coli((m32-p34)/(-p24),dreal(p24-p34+m32)))
     &          -.5d0*(cln_coli((m32-p23)/(m32-p34),dreal(p34-p23)))**2
     &          +cspcon_coli(m32/(m32-p34),-p24/(m32-p23),
     &          -m32*p24/((m32-p34)*(m32-p23)),
     &          ((m32-p34)*(m32-p23)+m32*p24)/((m32-p34)*(m32-p23)),
     &          1d0,dreal(p23-p24))
c
     &          +.5d0*(cln_coli(-p24/muir2,-1d0))**2
     &          -2d0*(
     &          cspcon_coli(m32-p34,1d0/(m32-p13),(m32-p34)/(m32-p13),
     &          (p34-p13)/(m32-p13),-1d0,1d0)
     &          +cspcon_coli(m32-p23,1d0/(m32-p13),(m32-p23)/(m32-p13),
     &          (p23-p13)/(m32-p13),-1d0,1d0) )
     &          -2d0*pi2_6 
            
            D0ms2ir1_coli = D0ms2ir1_coli/((p13-m32)*p24)

          elseif (m32.eq.p34) then
    
* 1 soft singularities and 3 zero masses
*     
*                   0             
*    p12=m12 ----------------  p23    
*                |  2   |         
*            m12 |1    3| m32       
*                |   4  |         
*    p14=m12 ----------------  p34=m32
*                   0             
*

#ifdef CHECK
            if(m12.eq.cd0.or.m22.ne.0d0.or.m42.ne.0d0
     &          .or.p12.ne.m12.or.p14.ne.m12.or.
     &          p34.ne.m32.or.p23.eq.m32)
     &          then
              write(*,100)' D0ms1ir2_coli:  (3.29)IR ',
     &            '    case with 1 soft singularity',
     &            '    and 2 opposite zero masses: ',
     &            '    wrong arguments'
              write(*,111)' D0ms2ir1_coli: p12 = ',p12
              write(*,111)' D0ms2ir1_coli: p23 = ',p23
              write(*,111)' D0ms2ir1_coli: p34 = ',p34
              write(*,111)' D0ms2ir1_coli: p14 = ',p14
              write(*,111)' D0ms2ir1_coli: p24 = ',p24
              write(*,111)' D0ms2ir1_coli: p13 = ',p13
              write(*,111)' D0ms2ir1_coli: m12 = ',m12
              write(*,111)' D0ms2ir1_coli: m22 = ',m22
              write(*,111)' D0ms2ir1_coli: m32 = ',m32
              write(*,111)' D0ms2ir1_coli: m42 = ',m42
              stop
            endif
#endif

#ifdef CHECK
            if(flag(2))then
              write(*,100)' D0ms2ir1_coli:  (3.29) IR ',
     &            '    case with 1 soft singularity ',
     &            '    and 2 opposite zero masses:',
     &            '    not yet tested at all'
              write(*,111)' D0ms2ir1_coli: p12 = ',p12
              write(*,111)' D0ms2ir1_coli: p23 = ',p23
              write(*,111)' D0ms2ir1_coli: p34 = ',p34
              write(*,111)' D0ms2ir1_coli: p14 = ',p14
              write(*,111)' D0ms2ir1_coli: p13 = ',p13
              write(*,111)' D0ms2ir1_coli: p24 = ',p24
              write(*,111)' D0ms2ir1_coli: m12 = ',m12
              write(*,111)' D0ms2ir1_coli: m22 = ',m22
              write(*,111)' D0ms2ir1_coli: m32 = ',m32
              write(*,111)' D0ms2ir1_coli: m42 = ',m42
              flag(2)=.false.
              stop
            endif
#endif

c     cD0ir0e
c     according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349  (iic)
c     Denner, Dittmaier (3.29) with IR singularity 

c            k13 = (m32-p13)/sqrt(m12*coliminfscale2*m32)
c            logxs = -cln_coli(k13,-1d0)
            logxs= -(cln_coli((m32-p13)/m32,-1d0)
     &          +cln_coli((m32-p13)/(m12*coliminfscale2),-1d0))/2d0
c     xs = 1/k24   

              D0ms2ir1_coli = ( logxs*(
     &            - logxs - log(muir2/(m12*coliminfscale2))
     &            + 2d0*cln_coli(p24/p12,dreal(p12)) )
     &            - pi2_6 )

              D0ms2ir1_coli =  D0ms2ir1_coli/((p13-m32)*(-p24))

#ifdef CHECK
            else
              write(*,100)' D0ms2ir1_coli: case not implemented'
              write(*,111)' D0ms2ir1_coli: p12 = ',p12
              write(*,111)' D0ms2ir1_coli: p23 = ',p23
              write(*,111)' D0ms2ir1_coli: p34 = ',p34
              write(*,111)' D0ms2ir1_coli: p14 = ',p14
              write(*,111)' D0ms2ir1_coli: p13 = ',p13
              write(*,111)' D0ms2ir1_coli: p24 = ',p24
              write(*,111)' D0ms2ir1_coli: m12 = ',m12
              write(*,111)' D0ms2ir1_coli: m22 = ',m22
              write(*,111)' D0ms2ir1_coli: m32 = ',m32
              write(*,111)' D0ms2ir1_coli: m42 = ',m42
              stop
#endif

            endif

          elseif(m22.eq.cd0.or.m42.eq.cd0)then
* 1 soft singularities and 2 zero masses
*
*    small        m22 small             
*    p12=m22 ----------------  p23    
*                |  2   |         
*              0 |1    3| m32       
*                |   4  |         
*          0 ----------------  p34
*                    0             
*

#ifdef CHECK
          if(m12.ne.cd0.or.p14.ne.m42.or.p12.ne.m22.or.
     &        .not.(m22.eq.cd0.and.m42.ne.cd0.and.p23.ne.m32.or.
     &              m42.eq.cd0.and.m22.ne.cd0.and.p34.ne.m32))
     &        then
            write(*,100)' D0ms1ir2_coli:  (3.27)',
     &          '    case with 1 soft singularities',
     &          '    and 2 adjacen zero masses: ',
     &          '    wrong arguments'
            write(*,111)' D0ms2ir1_coli: p12 = ',p12
            write(*,111)' D0ms2ir1_coli: p23 = ',p23
            write(*,111)' D0ms2ir1_coli: p34 = ',p34
            write(*,111)' D0ms2ir1_coli: p14 = ',p14
            write(*,111)' D0ms2ir1_coli: p24 = ',p24
            write(*,111)' D0ms2ir1_coli: p13 = ',p13
            write(*,111)' D0ms2ir1_coli: m12 = ',m12
            write(*,111)' D0ms2ir1_coli: m22 = ',m22
            write(*,111)' D0ms2ir1_coli: m32 = ',m32
            write(*,111)' D0ms2ir1_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(3))then
            write(*,100)' D0ms2ir1_coli:  (3.27)',
     &          '    case with 1 soft singularities',
     &          '    and 2 adjacent zero masses ',
     &          '    not yet tested at all'
            write(*,111)' D0ms2ir1_coli: p12 = ',p12
            write(*,111)' D0ms2ir1_coli: p23 = ',p23
            write(*,111)' D0ms2ir1_coli: p34 = ',p34
            write(*,111)' D0ms2ir1_coli: p14 = ',p14
            write(*,111)' D0ms2ir1_coli: p24 = ',p24
            write(*,111)' D0ms2ir1_coli: p13 = ',p13
            write(*,111)' D0ms2ir1_coli: m12 = ',m12
            write(*,111)' D0ms2ir1_coli: m22 = ',m22
            write(*,111)' D0ms2ir1_coli: m32 = ',m32
            write(*,111)' D0ms2ir1_coli: m42 = ',m42
            flag(3)=.false.
          endif
#endif

          if(m22.ne.cd0) then
            mm22=m22*coliminfscale2
            q23=p23
            q34=p34
          else
            mm22=m42*coliminfscale2
            q23=p34
            q34=p23
          endif

c according to Ellis, Zanderighi, JHEP 0802 (2008) 002  
c              box12 m3 small               
c Denner Dittmaier (3.27)

c>          write(*,*) 'D0ms2ir1_coli 1 ',-p24/(m32-q34),m32/(m32-q23),
c>     &        (-p24*m32)/((m32-q34)*(m32-q23)),
c>     &        ((m32-q34)*(m32-q23)+p24*m32)/((m32-q34)*(m32-q23)),
c>     &        -dreal(p24-q34),1d0
c>          write(*,*) 'D0ms2ir1_coli 2 ',(m32-q34)/(-p24),(m32-q23)/m32,
c>     &        (m32-q34)*(m32-q23)/(-p24*m32),
c>     &        ((m32-q34)*(m32-q23)+p24*m32)/(p24*m32),
c>     &        dreal(p24-q34),-1d0
c>          write(*,*) 'D0ms2ir1_coli 3 ',-p24/(m32-q34),m32/(m32-q23),
c>     &        (-p24*m32)/((m32-q34)*(m32-q23)),
c>     &        ((m32-q34)*(m32-q23)+p24*m32)/((m32-q34)*(m32-q23)),
c>     &        -dreal(p24-q34),1d0
c>          write(*,*) 'D0ms2ir1_coli 4 ',
c>     &        (m32-q34),1d0/(m32-p13),(m32-q34)/(m32-p13),
c>     &        (q34-p13)/(m32-p13),-1d0,1d0


          D0ms2ir1_coli=
#ifdef SING
     &        .5d0*delta2ir
     &        -delta1ir*(cln_coli((-p24)/sqrt(mm22*muir2),-1d0)
     &        +cln_coli((m32-p13)/(m32-q34),dreal(q34-p13))) 
     &        -.25d0*colishiftms2
#endif
     &        -2d0*
     &        cspcon_coli((m32-q34),1d0/(m32-p13),(m32-q34)/(m32-p13),
     &        (q34-p13)/(m32-p13),-1d0,1d0)
     &        +2d0*cln_coli(-p24/sqrt(mm22*muir2),-1d0)
     &        *cln_coli((m32-p13)/sqrt(mm22*muir2),-1d0)
     &        -(cln_coli((m32-q34)/sqrt(mm22*muir2),-1d0))**2
     &        -.5d0*(cln_coli((m32-q23)/mm22,-1d0))**2
     &        -2d0*pi2_6
     &        +cspcon_coli(-p24/(m32-q34),m32/(m32-q23),
     &        (-p24*m32)/((m32-q34)*(m32-q23)),
     &        ((m32-q34)*(m32-q23)+p24*m32)/((m32-q34)*(m32-q23)),
     &        -dreal(p24-q34),1d0)
     &        +.5d0*(cln_coli((m32-q34)/(-p24),dreal(p24-q34)))**2
     &        +cln_coli((m32-q34)/(-p24),dreal(p24-q34))
     &              *cln_coli((m32-q23)/mm22,-1d0)


c>              if(m32.ne.cd0)then
c>                write(*,*) 'test1=',
c>     &        -cspcon_coli((m32-q34)/(-p24),(m32-q23)/m32,
c>     &        (m32-q34)*(m32-q23)/(-p24*m32),
c>     &        ((m32-q34)*(m32-q23)+p24*m32)/(p24*m32),
c>     &        dreal(p24-q34),-1d0)
c>     &        -.5d0*(cln_coli((m32-q23)/m32,-1d0))**2
c>     &        +cln_coli((m32-q34)/(-p24),dreal(p24-q34))*log(m32/mm22)
c>                write(*,*) 'test2=',
c>     &        +cspcon_coli(-p24/(m32-q34),m32/(m32-q23),
c>     &        (-p24*m32)/((m32-q34)*(m32-q23)),
c>     &        ((m32-q34)*(m32-q23)+p24*m32)/((m32-q34)*(m32-q23)),
c>     &        -dreal(p24-q34),1d0)
c>     &        +.5d0*(cln_coli((m32-q34)/(-p24),dreal(p24-q34)))**2
c>     &        +cln_coli((m32-q34)/(-p24),dreal(p24-q34))
c>     &              *cln_coli((m32-q23)/mm22,-1d0)
c>              endif

          D0ms2ir1_coli=D0ms2ir1_coli/((p13-m32)*p24)
          

        else
* 1 soft singularities and 1 zero masses
*
*    small       m22 small                     
*    p12=m22 ----------------  p23    
*                |  2   |         
*              0 |1    3| m32       
*                |   4  |         
*    p14=m42 ----------------  p34
*    small       m42 small             
*

#ifdef CHECK
          if(m12.ne.cd0.or.p14.ne.m42.or.p12.ne.m22.or.
     &       m22.eq.cd0.or.m42.eq.cd0)
     &        then
            write(*,100)' D0ms2ir2_coli:  (3.31)',
     &          '    case with 1 soft singularities and 1 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms2ir1_coli: p12 = ',p12
            write(*,111)' D0ms2ir1_coli: p23 = ',p23
            write(*,111)' D0ms2ir1_coli: p34 = ',p34
            write(*,111)' D0ms2ir1_coli: p14 = ',p14
            write(*,111)' D0ms2ir1_coli: p24 = ',p24
            write(*,111)' D0ms2ir1_coli: p13 = ',p13
            write(*,111)' D0ms2ir1_coli: m12 = ',m12
            write(*,111)' D0ms2ir1_coli: m22 = ',m22
            write(*,111)' D0ms2ir1_coli: m32 = ',m32
            write(*,111)' D0ms2ir1_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef ALLCHECK
c ppenj
          if(flag(4))then
            write(*,100)' D0ms2ir1_coli:  (3.31)',
     &          '    case with 1 soft singularities and 1 zero masses',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms2ir1_coli: p12 = ',p12
            write(*,111)' D0ms2ir1_coli: p23 = ',p23
            write(*,111)' D0ms2ir1_coli: p34 = ',p34
            write(*,111)' D0ms2ir1_coli: p14 = ',p14
            write(*,111)' D0ms2ir1_coli: p24 = ',p24
            write(*,111)' D0ms2ir1_coli: p13 = ',p13
            write(*,111)' D0ms2ir1_coli: m12 = ',m12
            write(*,111)' D0ms2ir1_coli: m22 = ',m22
            write(*,111)' D0ms2ir1_coli: m32 = ',m32
            write(*,111)' D0ms2ir1_coli: m42 = ',m42
            flag(4)=.false.
          endif
#endif

c according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 (ib)
c          and (iib)    (D0irem and D0ir0em)
c Denner Dittmaier (3.31)
      
          if(m32.ne.cd0)then
            D0ms2ir1_coli= 
     &          -cln_coli(-p24/(sqrt(m22*m42)*coliminfscale2),-1d0)
     &          *(cln_coli((m32-p13)/muir2,-1d0)
     &          +cln_coli((m32-p13)/m32,-1d0))
     &          + (cln_coli((m32-p23)/sqrt(m22*m32*coliminfscale2),
     &          -1d0))**2 
     &          + (cln_coli((m32-p34)/sqrt(m42*m32*coliminfscale2),
     &          -1d0))**2 
     &          + 3d0*pi2_6
     &          - cspenc_coli(-(p23-m32)*(p34-m32)/(m32*p24),1d0) 
     &          + cln_coli(1d0+(p23-m32)*(p34-m32)/(m32*p24),-1d0)
     &          *(cln_coli(-p24/(sqrt(m22*m42)*coliminfscale2),-1d0) 
     &          -cln_coli((m32-p23)/sqrt(m22*m32*coliminfscale2),-1d0)
     &          -cln_coli((m32-p34)/sqrt(m42*m32*coliminfscale2),-1d0))
#ifdef SING
     &          +cln_coli(-p24/(sqrt(m22*m42)*coliminfscale2),-1d0)
     &          *delta1ir
     &          + 0.5d0*colishiftms2
#endif
 
            D0ms2ir1_coli = -D0ms2ir1_coli/((p13-m32)*p24)
           
          else
            D0ms2ir1_coli = 
     &          cln_coli(-p24/(sqrt(m22*m42)*coliminfscale2),-1d0)
     &          *(.5d0*cln_coli(-p24/(sqrt(m22*m42)*coliminfscale2),
     &          -1d0)
     &          - log(muir2/(sqrt(m22*m42)*coliminfscale2))
     &          +cln_coli(p13/p23,dreal(p23)) 
     &          + cln_coli(p13/p34,dreal(p34)) )
     &          - pi2_6
     &          - (cln_coli(sqrt(m42/m22)*p23/p34,dreal(p34)))**2/2d0 
#ifdef SING
     &          +cln_coli(-p24/(sqrt(m22*m42)*coliminfscale2),-1d0)
     &          *delta1ir
     &           - 0.5d0*colishiftms2
#endif

            D0ms2ir1_coli = D0ms2ir1_coli/(p13*p24)
           
          endif
        endif

      elseif(nsoft.eq.0)then

        if(m22.eq.cd0.and.m42.eq.cd0)then
* 1 soft singularities and 3 zero masses
*
*                   0             
*    p12=m12 ----------------  p23    
*                |  2   |         
*      m12 small |1    3| m32       
*                |   4  |         
*    p14=m12 ----------------  p34
*                   0             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.ne.0d0.or.m42.ne.0d0
     &        .or.p12.ne.m12.or.p14.ne.m12.or.
     &        p34.eq.m32.or.p23.eq.m32)
     &        then
            write(*,100)' D0ms1ir2_coli:  (3.29)',
     &          '    case with 0 soft singularities',
     &          '    and 2 zero masses: wrong arguments'
            write(*,111)' D0ms2ir1_coli: p12 = ',p12
            write(*,111)' D0ms2ir1_coli: p23 = ',p23
            write(*,111)' D0ms2ir1_coli: p34 = ',p34
            write(*,111)' D0ms2ir1_coli: p14 = ',p14
            write(*,111)' D0ms2ir1_coli: p24 = ',p24
            write(*,111)' D0ms2ir1_coli: p13 = ',p13
            write(*,111)' D0ms2ir1_coli: m12 = ',m12
            write(*,111)' D0ms2ir1_coli: m22 = ',m22
            write(*,111)' D0ms2ir1_coli: m32 = ',m32
            write(*,111)' D0ms2ir1_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(5))then
            write(*,100)' D0ms2ir1_coli:  (3.29)',
     &          '    case with 0 soft singularities and 2 zero masses ',
     &          '    not yet tested at all'
            write(*,111)' D0ms2ir1_coli: p12 = ',p12
            write(*,111)' D0ms2ir1_coli: p23 = ',p23
            write(*,111)' D0ms2ir1_coli: p34 = ',p34
            write(*,111)' D0ms2ir1_coli: p14 = ',p14
            write(*,111)' D0ms2ir1_coli: p24 = ',p24
            write(*,111)' D0ms2ir1_coli: p13 = ',p13
            write(*,111)' D0ms2ir1_coli: m12 = ',m12
            write(*,111)' D0ms2ir1_coli: m22 = ',m22
            write(*,111)' D0ms2ir1_coli: m32 = ',m32
            write(*,111)' D0ms2ir1_coli: m42 = ',m42
            flag(5)=.false.
          endif
#endif

c D0mseee00
c Denner Dittmaier (3.29)

          x1 = p13*p24/(p23*p34)
          x2 = p12/p13*coliminfscale2
          ix1 = dreal(x2)
          ix2 = -ix1

          D0ms2ir1_coli = pi2_6 
     &        + (cln_coli(-x2,-ix2) - cln_coli(-x1, -ix1))
     &        *( cln_coli(p23/(p13),dreal(p13-p23)) 
     &        + cln_coli(p34/p24,dreal(p24-p34))
     &        + 0.5d0*cln_coli(x1*x2,-dreal(x1*x2)))
      
c>      write(*,*) 'D0ms2ir1_coli x',x1,x2,ix1,ix2
c>      write(*,*) 'D0ms2ir1_coli s',
c>     &     pi2_6, 
c>     &      + (cln_coli(-x2*minfscale2,-ix2) - cln_coli(-x1, -ix1)),
c>     &      cln_coli(p23/(p13),dreal(p13-p23)), 
c>     &     + cln_coli(p34/p24,dreal(p24-p34)),
c>     &     + 0.5d0*cln_coli(x1*x2*minfscale2,-dreal(x1*x2))
c>
c>       write(*,*) 'D0ms2ir1_coli s',D0ms2ir1_coli

          D0ms2ir1_coli = D0ms2ir1_coli/(p13*p24) 

c       write(*,*) 'D0ms2ir1_coli= ',D0ms2ir1_coli



#ifdef SING
#endif

        elseif(m22.eq.cd0.or.m42.eq.cd0)then
* 1 soft singularities and 1 zero masses
*
*                m22=m12             
*          0 ----------------  p23    
*                |  2   |         
*      m12 small |1    3| m32       
*                |   4  |         
*    p14=m12 ----------------  p34
*    small           0                 
*

#ifdef CHECK
          if(m12.eq.cd0.or.p14.ne.m22.or.p12.ne.m42.or.
     &        .not.(m22.eq.m12.and.m42.eq.cd0.and.p34.ne.m32.or.
     &              m42.eq.m12.and.m22.eq.cd0.and.p23.ne.m32))
     &        then
            write(*,100)' D0ms1ir2_coli: (3.28)',
     &          '    case with 0 soft singularities and 1 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms2ir1_coli: p12 = ',p12
            write(*,111)' D0ms2ir1_coli: p23 = ',p23
            write(*,111)' D0ms2ir1_coli: p34 = ',p34
            write(*,111)' D0ms2ir1_coli: p14 = ',p14
            write(*,111)' D0ms2ir1_coli: p24 = ',p24
            write(*,111)' D0ms2ir1_coli: p13 = ',p13
            write(*,111)' D0ms2ir1_coli: m12 = ',m12
            write(*,111)' D0ms2ir1_coli: m22 = ',m22
            write(*,111)' D0ms2ir1_coli: m32 = ',m32
            write(*,111)' D0ms2ir1_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef ALLCHECK
c ppenj
          if(flag(6))then
            write(*,100)' D0ms2ir1_coli:  (3.28)',
     &          '    case with 0 soft singularities and 1 zero masses ',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms2ir1_coli: p12 = ',p12
            write(*,111)' D0ms2ir1_coli: p23 = ',p23
            write(*,111)' D0ms2ir1_coli: p34 = ',p34
            write(*,111)' D0ms2ir1_coli: p14 = ',p14
            write(*,111)' D0ms2ir1_coli: p24 = ',p24
            write(*,111)' D0ms2ir1_coli: p13 = ',p13
            write(*,111)' D0ms2ir1_coli: m12 = ',m12
            write(*,111)' D0ms2ir1_coli: m22 = ',m22
            write(*,111)' D0ms2ir1_coli: m32 = ',m32
            write(*,111)' D0ms2ir1_coli: m42 = ',m42
            flag(6)=.false.
          endif
#endif

c permutation 3 -> 4 -> 1 -> 2 -> 3   m22=0
c or          3 -> 2 -> 1 -> 4 -> 3   m42=0
          if(m42.eq.cd0) then
            mm12=m22
            mm22=m32
            mm42=m12
            q12=p23
            q23=p34
            q13=p24
            q24=p13
          else
            mm12=m42
            mm42=m12
            mm22=m32
            q12=p34
            q23=p23
            q13=p24
            q24=p13
          endif

c>          write(*,111)' D0ms2ir1_coli: q12 = ',q12
c>          write(*,111)' D0ms2ir1_coli: q23 = ',q23
c>          write(*,111)' D0ms2ir1_coli: q34 = ',q34
c>          write(*,111)' D0ms2ir1_coli: q14 = ',q14
c>          write(*,111)' D0ms2ir1_coli: q24 = ',q24
c>          write(*,111)' D0ms2ir1_coli: q13 = ',q13
c>          write(*,111)' D0ms2ir1_coli: mm12= ',mm12
c>          write(*,111)' D0ms2ir1_coli: mm22= ',mm22
c>          write(*,111)' D0ms2ir1_coli: mm32= ',mm32
c>          write(*,111)' D0ms2ir1_coli: mm42= ',mm42

c cD0ms0ee
c Denner Dittmaier (3.28)

          mm12=mm12*coliminfscale2
          mm42=mm42*coliminfscale2

          l12 = (mm22    -q12)
          l13 = (        -q13)
          l14 = (mm42+mm12   )
          l23 = (mm22    -q23)
          l24 = (mm22    -q24)
          l34 = 0d0

          r12 = l12/mm22
          r21 = l12/mm12      
          r24 = mm22/l24
          r42 = mm42/l24      
          r14 = 1d0
          r41 = mm42/mm12
          
          if(dreal(l12).lt.-0d0) then
            ir12 = -1d0
          else
            ir12 = 0d0
          end if
          if(dreal(l24).lt.-0d0) then
            ir24 = 1d0
          else
            ir24 = 0d0
          end if

          y(1) = -l12/l24
          y(2) = l13*l24/(l23*mm12)
          z(1) = y(1)/r24
          z(2) = y(2)/r24
          d=l23
          iy(1) = -sign(1d0,dreal(d))
          iy(2) = +sign(1d0,dreal(d))
c     iy and iz correlated for l1 version!!!!
c     iy(1) = -iy(1)
c     iy(2) = -iy(2)
          iz(1) =  sign(1d0,iy(1)*dreal(r24))
          iz(2) =  sign(1d0,iy(2)*dreal(r24))

c>          write(*,*) 'D0ms2ir1_coli l ',l12,l13,l14
c>          write(*,*) 'D0ms2ir1_coli l ',l23,l24,l34
c>          write(*,*) 'D0ms2ir1_coli r ',r12,r14,r24
c>          write(*,*) 'D0ms2ir1_coli r ',r21,r41,r42
c>          write(*,*) 'D0ms2ir1_coli y ',y(1),y(2)
c>          write(*,*) 'D0ms2ir1_coli z ',z(1),z(2)
c>          write(*,*) 'D0ms2ir1_coli ir ',ir12,ir14,ir24
c>          write(*,*) 'D0ms2ir1_coli iy ',iy(1),iy(2)
c>          write(*,*) 'D0ms2ir1_coli iz ',iz(1),iz(2)

          D0ms2ir1_coli = dcmplx(0d0)
          do i=2,1,-1
            eta=eta2s_coli(-y(i),1d0/r24,-iy(i),-ir24,-iz(i))

c           write(*,*) 'D0ms2ir1_coli etai ',i,eta
c           write(*,*) 'D0ms2ir1_coli etai ',
c     &          -y(i),1d0/r24,-iy(i),-ir24,-iz(i)

            if(i.eq.1)then
              D0ms2ir1_coli = D0ms2ir1_coli 
     &            - cspcos_coli(-y(i),cd1,-iy(i),0d0)
     &            - cspcos_coli(-y(i),cd1,-iy(i),0d0)
c     &           + cspcos_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13))
     &            + pi2_6
     &            - cspcos_coli(dcmplx(1d0),cd0,0d0,0d0)
     &            - .5d0*(cln_coli(-z(i),-iz(i))+cln_coli(r21,ir12))**2
c     &           + cspcos_coli(-z(i),r21,-iz(i),ir12)
     &            + cspcos_coli(-z(i),1d0/r12,-iz(i),-ir12)
     &            - cspcos_coli(-z(i),l23/l13,-iz(i),dreal(l23-l13))

c>              write(*,*) 'D0ms2ir1_coli contr ',
c>     &            -cspcos_coli(-y(i),cd1,-iy(i),0d0)
c>              write(*,*) 'D0ms2ir1_coli contr ',
c>     &            -cspcos_coli(-y(i),cd1,-iy(i),0d0)
c>              write(*,*) 'D0ms2ir1_coli contr ',
c>c     &          cspcos_coli(-z(i),r21,-iz(i),ir12)
c>     &            -.5d0*(cln_coli(-z(i),-iz(i))+cln_coli(r21,ir12))**2
c>     &            -cspcos_coli(dcmplx(1d0),cd0,0d0,0d0)
c>              write(*,*) 'D0ms2ir1_coli contr ',
c>     &            cspcos_coli(-z(i),1d0/r12,-iz(i),-ir12)
c>              write(*,*) 'D0ms2ir1_coli contr ',
c>c     &           cspcos_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13))
c>     &            pi2_6
c>              write(*,*) 'D0ms2ir1_coli contr ',
c>     &            -cspcos_coli(-z(i),l23/l13,-iz(i),dreal(l23-l13))

              if(eta.ne.0d0)then
                l1 = cln_coli((mm12-l14*l12/l24+mm42*(l12/l24)**2)
     &              /l13,dreal(iy(i)/(l24*l23-l34*mm22)))

c>              write(*,*) 'D0ms2ir1_coli eta1 ',
c>     &            (mm12-l14*l12/l24+mm42*(l12/l24)**2)/l13
c>              write(*,*) 'D0ms2ir1_coli l1= ',l1
c>              write(*,*) 'D0ms2ir1_coli l2= ',
c>     &             cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))*l23/l13,-1d0)
c>     &                -cln_coli(l23,-1d0),
c>     &             cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
c>     &                -cln_coli(l13,-1d0),
c>     &             cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))
c>     &                *(l13+l23*z(i))/l13,-1d0)
c>     &                -cln_coli(l13+l23*z(i),-1d0)

                D0ms2ir1_coli = D0ms2ir1_coli 
     &              + eta*(l1+cln_coli(l13/mm12,-1d0)) 
                
c>                write(*,*) 'D0ms2ir1_coli contr ',
c>     &              eta*(l1+cln_coli(l13/mm12,-1d0)) 
c>                write(*,*) 'D0ms2ir1_coli eta ',
c>     &              eta,(l1+cln_coli(l13/mm12,-1d0)) 
c>                write(*,*) 'D0ms2ir1_coli l1 ',
c>     &              l1,cln_coli(l13/mm12,-1d0)

#ifdef CHECK
                if(conscheck)then
                  if(dimag(r24).ne.0d0) then
                    v=dimag(y(i))/dimag(r24)
                    u=dimag(y(i)/r24)/dimag(1d0/r24)
                    l2=cln_coli(l12*v+mm22*v*v+l24*u*v,-1d0)
     &                  -cln_coli(l13+l23*v,-1d0)

c     write(*,*) 'D0ms2ir1_coli l2= ',l2,l1

                  elseif(dimag(y(i)).ne.0d0) then
                    l2=cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))*l23/l13,
     &                  -1d0)
     &                  -cln_coli(l23,-1d0)
                  elseif(dreal(y(i)).gt.0d0)then
                    l2=cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
     &                  -cln_coli(l13,-1d0)
                  elseif(dreal(z(i)).gt.0d0)then
                    l2=cln_coli((mm12+l14*y(i)+mm42*y(i)*y(i))
     &                  *(l13+l23*z(i))/l13,-1d0)
     &                  -cln_coli(l13+l23*z(i),-1d0)
                  else
                    if (flag3) then
                      write(*,*) 
                      write(*,*) 
     &                    'D0ms2ir1_coli: case for l2 not implemented '
                      flag3 = .false.
                    end if
                    write(errout,*) 
                    write(errout,*)
     &                  'D0ms2ir1_coli: case for l2 not implemented'
                    write(errout,*) 
     &                  'D0ms2ir1_coli args ',q12,q23,q34,q14,
     &                  q13,q24,mm12,mm22,mm32,mm42
                    write(errout,*) 'r24 = ',r24
                    write(errout,*) 'y(i) = ',y(i)
                    write(errout,*) 'z(i) = ',z(i)
                  endif 
                  
                  if(abs(l2/l1-1d0).gt.testacc) then
                    if (flag2) then
                      write(*,*) 
                      write(*,*) 
     &                    'D0ms2ir1_coli: eta terms possibly wrong '
                      flag2 = .false.
                    end if
                    write(errout,*) 
                    write(errout,*) 
     &                  'D0ms2ir1_coli: etaterms  possibly wrong '
                    write(errout,*) 'D0ms2ir1_coli args p ',
     &                  p12,p23,p34,p14,
     &                  p13,p24,m12,m22,m32,m42
                    write(errout,*) 'D0ms2ir1_coli args q ',
     &                  q12,q23,q34,q14,
     &                  q13,q24,mm12,mm22,mm32,mm42
                    write(errout,*) 'l1 = ',l1
                    write(errout,*) 'l2 = ',l2
                    write(errout,*) 'diff   = ',abs((l2-l1)/l1)
                  endif 
                endif
#endif
              endif

c              write(*,*) 'D0ms2ir1_coli sum ',D0ms2ir1_coli

            else

              D0ms2ir1_coli = D0ms2ir1_coli +
c     &            cspcos_coli(-y(2),cd1,-iy(2),0d0)
     &            -pi2_6-.5d0*cln_coli(-y(2),-iy(2))**2
c     &            cspcos_coli(-y(2),cd1,-iy(2),0d0)
     &            -pi2_6-.5d0*cln_coli(-y(2),-iy(2))**2
c     &           - cspcos_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13))
     &            -pi2_6
c     &           - cspcos_coli(-z(i),r21,-iz(i),ir12)
     &            +pi2_6+
     &            .5d0*(cln_coli(-z(i),-iz(i))+cln_coli(r21,ir12))**2
c     &           - cspcos_coli(-z(i),1d0/r12,-iz(i),-ir12)
     &            +pi2_6
     &            +.5d0*
     &            (cln_coli(-z(2),-iz(2))+cln_coli(1d0/r12,-ir12))**2
c     &           + cspcos_coli(-z(i),l23/l13,-iz(i),dreal(l23-l13))
     &            -pi2_6-.5d0*(cln_coli(-z(2),-iz(2))
     &            +cln_coli(l23/l13,dreal(l23-l13)))**2

c>              write(*,*) 'D0ms2ir1_coli contr ',
c>     &            -pi2_6-.5d0*cln_coli(-y(2),-iy(2))**2
c>              write(*,*) 'D0ms2ir1_coli contr ',
c>     &            -pi2_6-.5d0*cln_coli(-y(2),-iy(2))**2
c>              write(*,*) 'D0ms2ir1_coli contr ',
c>c     &           - cspcos_coli(-z(i),r21,-iz(i),ir12)
c>     &            pi2_6
c>     &            +.5d0*(cln_coli(-z(i),-iz(i))+cln_coli(r21,ir12))**2
c>              write(*,*) 'D0ms2ir1_coli contr ',
c>c     &           cspcos_coli(-z(i),1d0/r12,-iz(i),-ir12)
c>     &            +pi2_6
c>     &            +.5d0*
c>     &            (cln_coli(-z(2),-iz(2))+cln_coli(1d0/r12,-ir12))**2
c>              write(*,*) 'D0ms2ir1_coli contr ',
c>c     &           - cspcos_coli(-y(i),l34/l13,-iy(i),dreal(l34-l13))
c>     &            -pi2_6
c>              write(*,*) 'D0ms2ir1_coli contr ',
c>c     &           cspcos_coli(-z(i),l23/l13,-iz(i),dreal(l23-l13))
c>     &            -pi2_6-.5d0*(cln_coli(-z(2),-iz(2))
c>     &            +cln_coli(l23/l13,dreal(l23-l13)))**2
              
              if(eta.ne.0d0)then
                l1 = cln_coli(l24*y(2)/l23,dreal(iy(i)/(l24*l23)))

c>              write(*,*) 'D0ms2ir1_coli eta2 ',
c>     &                  l24*y(2)/l23
c>              write(*,*) 'D0ms2ir1_coli l1= ',l1
c>              write(*,*) 'D0ms2ir1_coli l2= ',
c>     &            cln_coli(l24*y(2),-1d0)-cln_coli(l23,-1d0),
c>     &            cln_coli(mm42*y(2)*y(2),-1d0)-cln_coli(l13,-1d0),
c>     &            cln_coli(mm22*z(2)*z(2),-1d0)-cln_coli(l23*z(2),-1d0)


                D0ms2ir1_coli = D0ms2ir1_coli 
     &              - eta*(l1+cln_coli(l13/mm12,-1d0)) 

c>                write(*,*) 'D0ms2ir1_coli contr ',
c>     &              - eta*(l1+cln_coli(l13/mm12,-1d0)) 
c>                write(*,*) 'D0ms2ir1_coli eta ',
c>     &              eta,(l1+cln_coli(l13/mm12,-1d0)) 
c>                write(*,*) 'D0ms2ir1_coli eta ',
c>     &              l1,cln_coli(l13/mm12,-1d0)

#ifdef CHECK
                if(conscheck)then
                  if(dimag(r24).ne.0d0) then
                    v=dimag(y(i))/dimag(r24)
                    u=dimag(y(i)/r24)/dimag(1d0/r24)
                    if(v.ne.0)then
c                      l2=cln_coli(l12*v+mm22*v*v+l24*u*v,-1d0)
c     &                    -cln_coli(l13+l23*v,-1d0)
                      l2=cln_coli(mm22*v+l24*u,-1d0)
     &                    -cln_coli(l23,-1d0)

c                      write(*,*) 'D0ms2ir1_coli compb l2= ',l2,l1,l2/l1

                    else
                      l2=cln_coli(mm12+l14*y(i)+mm42*y(i)*y(i),-1d0)
     &                    -cln_coli(l13+l34*y(i),-1d0)
c     &                    -cln_coli((1d0+y(2))*(1d0+y(2))
c     &                    *(l13+l23*z(2))/
c     &                    ((1+z(2)/r12)*z(2)*r21),-1d0)

c                      write(*,*) 'D0ms2ir1_coli compb l2= ',l2,l1,l2/l1

                    endif
                  elseif(dimag(y(2)).ne.0d0) then
                    l2=cln_coli(l24*y(2),-1d0)-cln_coli(l23,-1d0)
                  elseif(dreal(y(i)).gt.0d0)then
                    l2=cln_coli(mm42*y(2)*y(2),-1d0)-cln_coli(l13,-1d0)
                  elseif(dreal(z(i)).gt.0d0)then
                    l2=cln_coli(mm22*z(2)*z(2),-1d0)
     &                  -cln_coli(l23*z(2),-1d0)
                  else
                    if (flag3) then
                      write(*,*) 
                      write(*,*) 
     &                    'D0ms2ir1_coli: case for l2 not implemented '
                      flag3 = .false.
                    end if
                    write(errout,*) 
                    write(errout,*)
     &                  'D0ms2ir1_coli: case for l2 not implemented'
                    write(errout,*) 'D0ms2ir1_coli args p ',
     &                  p12,p23,p34,p14,p13,
     &                  p24,m12,m22,m32,m42
                    write(errout,*) 'D0ms2ir1_coli args q ',
     &                  q12,q23,q34,q14,q13,
     &                  q24,mm12,mm22,mm32,mm42
                    write(errout,*) 'r24 = ',r24
                    write(errout,*) 'y(i) = ',y(i)
                    write(errout,*) 'z(i) = ',z(i)
                  endif 
                  
                  if(abs(l2/l1-1d0).gt.testacc) then
                    if (flag2) then
                      write(*,*) 
                      write(*,*) 
     &                    'D0ms2ir1_coli: eta terms possibly wrong '
                      flag2 = .false.
                    end if
                    write(errout,*) 
                    write(errout,*) 
     &                  'D0ms2ir1_coli: etaterms  possibly wrong '
                    write(errout,*) 'D0ms2ir1_coli args ',
     &                  q12,q23,q34,q14,q13,
     &                  q24,mm12,mm22,mm32,mm42
                    write(errout,*) 'l1 = ',l1
                    write(errout,*) 'l2 = ',l2
                    write(errout,*) 'diff   = ',abs((l2-l1)/l1)
                  endif 
                endif
#endif
              endif

c              write(*,*) 'D0ms2ir1_coli sum ',D0ms2ir1_coli

            endif

          enddo

c     write(*,*) 'D0ms2ir1_coli sum ',D0ms2ir1_coli

#ifdef SING
          D0ms2ir1_coli = D0ms2ir1_coli  + 0.5d0*colishiftms2
#endif

          D0ms2ir1_coli = D0ms2ir1_coli/(l13*l24)

c     write(*,*) 'D0ms2ir1_coli =  ',D0ms2ir1_coli

#ifdef CHECK
          if (conscheck) then

          l12 = (mm22   -q12)
          l13 = (       -q13)
          l14 = (mm42+mm12  )
          l23 = (mm22   -q23)
          l24 = (mm22   -q24)
          l34 = 0d0

          r12 = l12/mm22
          r21 = l12/mm12      
          r24 = l24/mm42
          r42 = l24/mm22      
          r14 = 1d0
          r41 = mm42/mm12
          
          if(dreal(l12).lt.-0d0) then
            ir12 = -1d0
          else
            ir12 = 0d0
          end if
          if(dreal(l24).lt.-0d0) then
            ir24 = -1d0
          else
            ir24 = 0d0
          end if
          ir14=0d0

          y(1)=(l24*l14-mm42*l12)/(-l24*mm42)
          y(2)=-l13*l24/(l23*mm42)
          z(1) = y(1)/r24
          z(2) = y(2)/r24

          d=l23
          iy(1) = -sign(1d0,dreal(d))
          iy(2) = +sign(1d0,dreal(d))
          iz(1) =  sign(1d0,iy(1)*dreal(r24))
          iz(2) =  sign(1d0,iy(2)*dreal(r24))

c>          write(*,*) 'D0ms2ir1_check l ',l12,l13,l14
c>          write(*,*) 'D0ms2ir1_check l ',l23,l24,l34
c>          write(*,*) 'D0ms2ir1_check r ',r12,r14,r24
c>          write(*,*) 'D0ms2ir1_check r ',r21,r41,r42
c>          write(*,*) 'D0ms2ir1_check y ',y(1),y(2)
c>          write(*,*) 'D0ms2ir1_check z ',z(1),z(2)

          D0ms2ir1_check = dcmplx(0d0)
          
          eta=eta2s_coli(-y(1),1d0/r24,-iy(1),-ir24,-iz(1))
c>            ch1   = l12-r24*l14-mm22*(r42-1d0/r24)*y(i)
c>            test1 = abs(ch1)/max(abs(l12),abs(r42*y(i)),
c>     &          abs(r24*l14),abs(1/r24*y(i)))
c>            ch2    = mm12/z(i)+l12+mm22*z(i)
c>            test2 = abs(ch2)/max(abs(mm12/z(i)),abs(l12),
c>     &          abs(m22*z(i)))
c>            ch3   = l23+l13/z(i)
c>            test3 = abs(ch3)/max(abs(l23),abs(l13/z(i)))
c>            test2 = min(test2,test3)
c>            if (test1.lt.test2) then
c>              l1 = cln_coli(ch2/ch3,dreal(-m22*(r42-1d0/r24)*iy(i)/d))
c>            else
c>              l1 = cln_coli(ch1/d,dreal(-mm22*(r42-1d0/r24)*iy(i)/d))      
c>            end if
          if(eta.ne.0d0)then
            l1 = cln_coli((1d0+y(1))**2*mm12/l13,dreal(-iy(1)*l24/l23))

c>              write(*,*) 'D0ms2ir1_check l1= ',l1
c>              write(*,*) 'D0ms2ir1_check l2= ',
c>     &         cln_coli((mm12+l12*z(1))*l23/l13,-1d0)
c>     &            -cln_coli(l23,-1d0),
c>     &         cln_coli(mm12+l14*y(1)+mm42*y(1)*y(1),-1d0)
c>     &            -cln_coli(l13,-1d0),
c>     &         cln_coli(mm12+l12*z(1),-1d0)-cln_coli(l13,-1d0)

            if(conscheck)then
              if(dimag(r24).ne.0d0) then
                v=dimag(y(1))/dimag(r24)
                u=dimag(y(1)/r24)/dimag(1d0/r24)
                l2=cln_coli(mm12+l12*v+l14*u+mm42*u*u+l24*u*v,-1d0)
     &              -cln_coli(l13,-1d0)

                write(*,*) 'D0ms2ir1_check compa l2= ',l2,l1,l2/l1         

              elseif(dimag(y(1)).ne.0d0) then
                l2=cln_coli((mm12+l12*z(1))*l23/l13,-1d0)
     &              -cln_coli(l23,-1d0)
              elseif(dreal(y(1)).gt.0d0)then
                l2=cln_coli(mm12+l14*y(1)+mm42*y(1)*y(1),-1d0)
     &              -cln_coli(l13,-1d0)
              elseif(dreal(z(1)).gt.0d0)then
                l2=cln_coli(mm12+l12*z(1),-1d0)-cln_coli(l13,-1d0)
              else
                if (flag3) then
                  write(*,*) 
                  write(*,*) 
     &                'D0ms2ir1_check: case for l2 not implemented '
                  flag3 = .false.
                end if
                write(errout,*) 
                write(errout,*) 
     &              'D0ms2ir1_check: case for l2 not implemented '
                write(errout,*) 'D0ms2ir1_check args ',
     &              q12,q23,q34,q14,q13,q24,
     &              mm12,mm22,mm32,mm42
                write(errout,*) 'r24 = ',r24
                write(errout,*) 'y(1) = ',y(1)
                write(errout,*) 'z(1) = ',z(1)
              endif 
              
              if(abs(l2/l1-1d0).gt.testacc) then
                if (flag2) then
                  write(*,*) 
                  write(*,*) 'D0ms2ir1_check: eta terms possibly wrong '
                  flag2 = .false.
                end if
                write(errout,*) 
                write(errout,*) 
     &              'D0ms2ir1_check: eta terms  possibly wrong '
                write(errout,*) 'D0ms2ir1_check args ',
     &              q12,q23,q34,q14,q13,q24,
     &              mm12,mm22,mm32,mm42
                write(errout,*) 'l1 = ',l1
                write(errout,*) 'l2 = ',l2
                write(errout,*) 'diff   = ',abs((l2-l1)/l1)
              endif 
            endif
          else
            l1=0d0
          endif

c     write(*,*) 'D0ms2ir1_check eta 1',eta,l1,(mm12+l12*z(1))/l13
c     write(*,*) 'D0ms2ir1_check eta 2',-y(1),1d0/r24,-iy(1),-ir24,-iz(1)
          
          D0ms2ir1_check = D0ms2ir1_check 
     &        - cspcos_coli(-y(1),cd1,-iy(1),ir14)
     &        - cspcos_coli(-y(1),cd1,-iy(1),-ir14)
     &        + cspcos_coli(-z(1),r21,-iz(1),ir12)
     &        + pi2_6
     &        + pi2_6
     &        - pi2_6
     &        + eta*(l1+cln_coli(l13/mm12,-1d0)) 

c>          write(*,*) 'D0ms2ir1_check contr 1',
c>     &           - cspcos_coli(-y(1),cd1,-iy(1),ir14),
c>     &           - cspcos_coli(-y(1),cd1,-iy(1),-ir14),
c>     &           + cspcos_coli(-z(1),r21,-iz(1),ir12),
c>     &           + pi2_6,
c>     &           + pi2_6,
c>     &           - pi2_6,
c>     &           + eta*(l1+cln_coli(l13/mm12,-1d0))

          eta=eta2s_coli(-y(2),1d0/r24,-iy(2),-ir24,-iz(2))
          if(eta.ne.0d0)then
            l1 = cln_coli(-l24*y(2)/l23,dreal(-iy(2)*l24/l23))

c>              write(*,*) 'D0ms2ir1_check l1= ',l1
c>              write(*,*) 'D0ms2ir1_check l2= ',
c>     &         cln_coli(-y(2)*l24,-1d0)-cln_coli(l23,-1d0),
c>     &         cln_coli(mm42*y(2)*y(2),-1d0)-cln_coli(l13,-1d0),
c>     &         cln_coli(l12*z(2)+mm22*z(2)*z(2),-1d0)
c>     &            -cln_coli(-(l12*z(2)+mm22*z(2)*z(2))*l23/(y(2)*l24),-1d0)

            if(conscheck)then
              if(dimag(r24).ne.0d0) then
                v=dimag(y(2))/dimag(r24)
                u=dimag(y(2)/r24)/dimag(1d0/r24)
                l2=cln_coli(mm42*u*u+l24*u*v,-1d0)
     &              -cln_coli(l13+l23*v,-1d0)

                write(*,*) 'D0ms2ir1_check comp l2= ',l2,l1,l2/l1

              elseif(dimag(y(2)).ne.0d0) then
                l2=cln_coli(-y(2)*l24,-1d0)-cln_coli(l23,-1d0)
              elseif(dreal(y(2)).gt.0d0)then
                l2=cln_coli(mm42*y(2)*y(2),-1d0)-cln_coli(l13,-1d0)
              elseif(dreal(z(2)).gt.0d0)then
                l2=cln_coli(l12*z(2)+mm22*z(2)*z(2),-1d0)
     &              -cln_coli(-(l12*z(2)+mm22*z(2)*z(2))*l23/(y(2)*l24),
     &              -1d0)
              else
                if (flag3) then
                  write(*,*) 
                  write(*,*) 
     &                'D0ms2ir1_check: case for l2 not implemented '
                  flag3 = .false.
                end if
                write(errout,*) 
                write(errout,*) 
     &              'D0ms2ir1_check: case for l2 not implemented '
                write(errout,*) 'D0ms2ir1_check args ',
     &              q12,q23,q34,q14,q13,q24,
     &              mm12,mm22,mm32,mm42
                write(errout,*) 'r24 = ',r24
                write(errout,*) 'y(2) = ',y(2)
                write(errout,*) 'z(2) = ',z(2)
              endif 
              
              if(abs(l2/l1-1d0).gt.testacc) then
                if (flag2) then
                  write(*,*) 
                  write(*,*) 'D0ms2ir1_check: eta terms possibly wrong '
                  flag2 = .false.
                end if
                write(errout,*) 
                write(errout,*) 
     &              'D0ms2ir1_check: etaterms  possibly wrong '
                write(errout,*) 'D0ms2ir1_check args ',
     &              q12,q23,q34,q14,q13,q24,
     &              mm12,mm22,mm32,mm42
                write(errout,*) 'l1 = ',l1
                write(errout,*) 'l2 = ',l2
                write(errout,*) 'diff   = ',abs((l2-l1)/l1)
              endif 
            endif
          else
            l1=0d0
          endif

c     write(*,*) 'D0ms2ir1_check eta 2',eta,l1,-l24*y(2)/l23
c     write(*,*) 'D0ms2ir1_check eta 2',-y(2),1d0/r24,-iy(2),-ir24,-iz(2)

          D0ms2ir1_check = D0ms2ir1_check 
c     &           + cspcos_coli(-y(2),r41,-iy(2),ir14)
c     &           + cspcos_coli(-y(2),1d0/r14,-iy(i),-ir14)
     &      -pi2_6-.5d0*(cln_coli(-y(2),-iy(2)))**2
     &      -pi2_6-.5d0*(cln_coli(-y(2),-iy(2)))**2
c     &        - cspcos_coli(-z(2),r21,-iz(2),ir12)
     &      +pi2_6+.5d0*
     &      (cln_coli(-z(2),-iz(2))+cln_coli(r21,ir12))**2
     &      - cspcos_coli(-z(2),1d0/r12,-iz(2),-ir12)
c     &           - cspcos_coli(-y(2),l34/l13,-iy(2),dreal(l34-l13))
     &      -pi2_6
c     &        + cspcos_coli(-z(2),l23/l13,-iz(2),dreal(l23-l13))
     &      +eta2s_coli(-z(2),l23/l13,-iz(2),dreal(l23-l13),0d0)
     &      *cln_coli((l12+mm22*z(2))*mm12/(l24*l24*z(2)),-iz(2))
     &      - eta*(l1+cln_coli(l13/mm12,-1d0)) 

c>          write(*,*) 'D0ms2ir1_check contr 2 ',
c>c     &           + cspcos_coli(-y(2),r41,-iy(2),ir14)
c>     &        -pi2_6-.5d0*(cln_coli(-y(2),-iy(2)))**2
c>c     &           + cspcos_coli(-y(2),1d0/r14,-iy(i),-ir14)
c>     &        -pi2_6-.5d0*(cln_coli(-y(2),-iy(2)))**2,
c>c     &           - cspcos_coli(-y(2),l34/l13,-iy(2),dreal(l34-l13))
c>     &        -pi2_6
c>c     &        - cspcos_coli(-z(2),r21,-iz(2),ir12)
c>     &        +pi2_6+.5d0*
c>     &        (cln_coli(-z(2),-iz(2))+cln_coli(r21,ir12))**2,
c>     &        - cspcos_coli(-z(2),1d0/r12,-iz(2),-ir12),
c>c     &           - cspcos_coli(-y(2),l34/l13,-iy(2),dreal(l34-l13))
c>     &        -pi2_6-.5d0*
c>     &        (cln_coli(-y(2),-iy(2))+cln_coli(r41,ir14))**2,
c>c     &        + cspcos_coli(-z(2),l23/l13,-iz(2),dreal(l23-l13)),
c>     &         +eta2s_coli(-z(2),l23/l13,-iz(2),dreal(l23-l13),0d0)
c>     &      *cln_coli((l12+mm22*z(2))*mm12/(l24*l24*z(2)),-iz(2)),
c>     &        - eta*(l1+cln_coli(l13/mm12,-1d0)) 
c>          write(*,*) 'D0ms2ir1_check sum ',D0ms2ir1_check

#ifdef SING
          D0ms2ir1_check = D0ms2ir1_check - 0.5d0*colishiftms2
#endif

          D0ms2ir1_check = D0ms2ir1_check/(-l13*l24)

c         write(*,*) 'D0ms2ir1_check =  ',D0ms2ir1_check          

          if(cdabs(D0ms2ir1_coli/D0ms2ir1_check-1d0).gt.testacc) then
            write(testout,*) ' D0ms2ir1_coli inconsistent ',i
            write(testout,111) ' D0ms2ir1_coli: p12 = ',p12
            write(testout,111) ' D0ms2ir1_coli: p23 = ',p23
            write(testout,111) ' D0ms2ir1_coli: p34 = ',p34
            write(testout,111) ' D0ms2ir1_coli: p14 = ',p14
            write(testout,111) ' D0ms2ir1_coli: p24 = ',p24
            write(testout,111) ' D0ms2ir1_coli: p13 = ',p13
            write(testout,111) ' D0ms2ir1_coli: m12 = ',m12
            write(testout,111) ' D0ms2ir1_coli: m22 = ',m22
            write(testout,111) ' D0ms2ir1_coli: m32 = ',m32
            write(testout,111) ' D0ms2ir1_coli: m42 = ',m42
            write(testout,111) ' D0ms2ir1_coli: q12 = ',q12
            write(testout,111) ' D0ms2ir1_coli: q23 = ',q23
            write(testout,111) ' D0ms2ir1_coli: q34 = ',q34
            write(testout,111) ' D0ms2ir1_coli: q14 = ',q14
            write(testout,111) ' D0ms2ir1_coli: q24 = ',q24
            write(testout,111) ' D0ms2ir1_coli: q13 = ',q13
            write(testout,111) ' D0ms2ir1_coli: mm12= ',mm12
            write(testout,111) ' D0ms2ir1_coli: mm22= ',mm22
            write(testout,111) ' D0ms2ir1_coli: mm32= ',mm32
            write(testout,111) ' D0ms2ir1_coli: mm42= ',mm42
            write(testout,111) ' D0ms2ir1_coli  = ',D0ms2ir1_coli
            write(testout,111) ' D0ms2ir1_check = ',D0ms2ir1_check
            write(testout,101) ' diff        = ',
     &          cdabs(D0ms2ir1_coli/D0ms2ir1_check-1d0)
            
            write(*,*)   ' D0ms2ir1_coli inconsistent ',i
            write(*,111) ' D0ms2ir1_coli: p12 = ',p12
            write(*,111) ' D0ms2ir1_coli: p23 = ',p23
            write(*,111) ' D0ms2ir1_coli: p34 = ',p34
            write(*,111) ' D0ms2ir1_coli: p14 = ',p14
            write(*,111) ' D0ms2ir1_coli: p24 = ',p24
            write(*,111) ' D0ms2ir1_coli: p13 = ',p13
            write(*,111) ' D0ms2ir1_coli: m12 = ',m12
            write(*,111) ' D0ms2ir1_coli: m22 = ',m22
            write(*,111) ' D0ms2ir1_coli: m32 = ',m32
            write(*,111) ' D0ms2ir1_coli: m42 = ',m42
            write(*,111) ' D0ms2ir1_coli: q12 = ',q12
            write(*,111) ' D0ms2ir1_coli: q23 = ',q23
            write(*,111) ' D0ms2ir1_coli: q34 = ',q34
            write(*,111) ' D0ms2ir1_coli: q14 = ',q14
            write(*,111) ' D0ms2ir1_coli: q24 = ',q24
            write(*,111) ' D0ms2ir1_coli: q13 = ',q13
            write(*,111) ' D0ms2ir1_coli: mm12= ',mm12
            write(*,111) ' D0ms2ir1_coli: mm22= ',mm22
            write(*,111) ' D0ms2ir1_coli: mm32= ',mm32
            write(*,111) ' D0ms2ir1_coli: mm42= ',mm42
            write(*,111) ' D0ms2ir1_coli  = ',D0ms2ir1_coli
            write(*,111) ' D0ms2ir1_check = ',D0ms2ir1_check
            write(*,101) ' diff        = ',
     &          cdabs(D0ms2ir1_coli/D0ms2ir1_check-1d0)
            stop
          endif
        endif
#endif


        else
* 1 soft singularities and 0 zero masses
*
*                m22=m12                     
*          0 ----------------  p23    
*                |  2   |         
*      m12 small |1    3| m32       
*                |   4  |         
*          0 ----------------  p34
*                m42=m12             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.ne.m12.or.m42.ne.m12.or.
     &       p12.ne.cd0.or.p14.ne.cd0)
     &        then
            write(*,100)' D0ms2ir1_coli:  (3.30)',
     &          '    case with 0 soft singularities',
     &          '    and 0 zero masses: wrong arguments'
            write(*,111)' D0ms2ir1_coli: p12 = ',p12
            write(*,111)' D0ms2ir1_coli: p23 = ',p23
            write(*,111)' D0ms2ir1_coli: p34 = ',p34
            write(*,111)' D0ms2ir1_coli: p14 = ',p14
            write(*,111)' D0ms2ir1_coli: p24 = ',p24
            write(*,111)' D0ms2ir1_coli: p13 = ',p13
            write(*,111)' D0ms2ir1_coli: m12 = ',m12
            write(*,111)' D0ms2ir1_coli: m22 = ',m22
            write(*,111)' D0ms2ir1_coli: m32 = ',m32
            write(*,111)' D0ms2ir1_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(6))then
            write(*,100)' D0ms2ir1_coli:  (3.30)',
     &          '    case with 0 soft singularities',
     &          '    and 0 zero masses not yet implemented'
            write(*,111)' D0ms2ir1_coli: p12 = ',p12
            write(*,111)' D0ms2ir1_coli: p23 = ',p23
            write(*,111)' D0ms2ir1_coli: p34 = ',p34
            write(*,111)' D0ms2ir1_coli: p14 = ',p14
            write(*,111)' D0ms2ir1_coli: p24 = ',p24
            write(*,111)' D0ms2ir1_coli: p13 = ',p13
            write(*,111)' D0ms2ir1_coli: m12 = ',m12
            write(*,111)' D0ms2ir1_coli: m22 = ',m22
            write(*,111)' D0ms2ir1_coli: m32 = ',m32
            write(*,111)' D0ms2ir1_coli: m42 = ',m42
            flag(6)=.false.
            stop
          endif
#endif

c Denner Dittmaier (3.30)
c new case
          D0ms2ir1_coli =cd0
#ifdef SING
#endif
 
          D0ms2ir1_coli = D0ms2ir1_coli/((p13-m32)*p24)

        endif

#ifdef CHECK
      else
        write(*,100)' D0ms2ir1_coli: case with more than 1 soft 
     &      singularities not implemented '
        write(*,111)' D0ms2ir1_coli: nsoft = ',nsoft
        write(*,111)' D0ms2ir1_coli: p23 = ',p23
        write(*,111)' D0ms2ir1_coli: p34 = ',p34
        write(*,111)' D0ms2ir1_coli: p14 = ',p14
        write(*,111)' D0ms2ir1_coli: p24 = ',p24
        write(*,111)' D0ms2ir1_coli: p13 = ',p13
        write(*,111)' D0ms2ir1_coli: m12 = ',m12
        write(*,111)' D0ms2ir1_coli: m22 = ',m22
        write(*,111)' D0ms2ir1_coli: m32 = ',m32
        write(*,111)' D0ms2ir1_coli: m42 = ',m42    
        stop
#endif
      endif

      end

************************************************************************
      function D0ms2ir2_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     4-point function                                                 *
*     with 2 mass singularities and up to 2 soft singularities         *
*                                                                      *
*     assumes                                                          *
*     p12, p23, m12, m22, m32 small                                    *
*     p13, p24, p14=m42, p34=/=m42 finite                              *
*                                                                      *
*                     m22 small                                        *
*     p12 small  ---------------------  p23 small                      *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | m32 small                            *
*                     |         |                                      *
*                     |    4    |                                      *
*       p14=m42  ---------------------  p34                            *
*                         m42                                          *
*                                                                      *
************************************************************************
*     22.08.08 Ansgar Denner        last changed  16.12.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms2ir2_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      logical flag(0:8)
      data    flag /9*.true./
      save    flag  
#endif

      complex*16 m2(4)
      complex*16 cln_coli,cspenc_coli
      integer    nsoft
      integer    i,j,k
      logical    soft(4),onsh(4,4)

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
        ms42 = elimminf2_coli(m42)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps34 = elimminf2_coli(p34)
        ps14 = elimminf2_coli(p14)
        ps24 = elimminf2_coli(p24)
        ps13 = elimminf2_coli(p13)

        if(ms12.ne.cd0.or.ms22.ne.cd0.or.ms32.ne.cd0
     &        .or.ms42.eq.cd0.or.ps12.ne.cd0.or.ps23.ne.cd0 
     &        .or.ps14.ne.ms42
     &        .or.ps24.eq.cd0.or.ps13.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &        .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,100)' D0ms2ir2_coli called improperly:'
          write(*,111)' D0ms2ir2_coli: p12 = ',p12,ps12
          write(*,111)' D0ms2ir2_coli: p23 = ',p23,ps23
          write(*,111)' D0ms2ir2_coli: p34 = ',p34,ps34
          write(*,111)' D0ms2ir2_coli: p14 = ',p14,ps14
          write(*,111)' D0ms2ir2_coli: p13 = ',p13,ps13
          write(*,111)' D0ms2ir2_coli: p24 = ',p24,ps24
          write(*,111)' D0ms2ir2_coli: m12 = ',m12,ms12
          write(*,111)' D0ms2ir2_coli: m22 = ',m22,ms22
          write(*,111)' D0ms2ir2_coli: m32 = ',m32,ms32
          write(*,111)' D0ms2ir2_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif

      m2(1)=m12
      m2(2)=m22
      m2(3)=m32
      m2(4)=m42

c determine on-shell momenta
      onsh(1,2)=p12.eq.m22
      onsh(1,3)=p13.eq.m32
      onsh(1,4)=p14.eq.m42
      onsh(2,1)=p12.eq.m12
      onsh(2,3)=p23.eq.m32
      onsh(2,4)=p24.eq.m42
      onsh(3,1)=p13.eq.m12
      onsh(3,2)=p23.eq.m22
      onsh(3,4)=p34.eq.m42
      onsh(4,1)=p14.eq.m22
      onsh(4,2)=p24.eq.m22
      onsh(4,3)=p34.eq.m32

c count/determine soft singularities
      nsoft=0 
      do i=1,4
        do j=1,3
          if(i.ne.j)then
            do k=j,4
              if(k.ne.i.and.k.ne.j)then
                soft(i)=m2(i).eq.cd0.and.onsh(i,j).and.onsh(i,k)
c                soft(i,k,j)=soft(i,j,k)
                if(soft(i)) nsoft=nsoft+1
              endif
            enddo
          endif
        enddo
      enddo

      if(nsoft.eq.2)then

        if(m32.eq.cd0)then
* 2 soft singularities and 3 zero masses
*
*                   0             
*          0 ----------------  0    
*                |  2   |         
*              0 |1    3| 0       
*                |   4  |         
*        m42 ----------------  p34
*                  m42             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.0d0.or.m32.ne.0d0.or.m42.eq.cd0
     &        .or.p12.ne.cd0.or.p23.ne.cd0.or.p34.eq.cd0.or.p14.ne.m42)
     &        then
            write(*,100)' D0ms2ir2_coli:  (3.33)',
     &          '    case with 2 soft singularities',
     &          '    and 3 zero masses: wrong arguments'
            write(*,111)' D0ms2ir2_coli: p12 = ',p12
            write(*,111)' D0ms2ir2_coli: p23 = ',p23
            write(*,111)' D0ms2ir2_coli: p34 = ',p34
            write(*,111)' D0ms2ir2_coli: p14 = ',p14
            write(*,111)' D0ms2ir2_coli: p24 = ',p24
            write(*,111)' D0ms2ir2_coli: p13 = ',p13
            write(*,111)' D0ms2ir2_coli: m12 = ',m12
            write(*,111)' D0ms2ir2_coli: m22 = ',m22
            write(*,111)' D0ms2ir2_coli: m32 = ',m32
            write(*,111)' D0ms2ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef ALLCHECK
c ppttbb dim
          if(flag(1))then
            write(*,100)' D0ms2ir2_coli:  (3.33)',
     &          '    case with 2 soft singularities and 3 zero masses',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms2ir2_coli: p12 = ',p12
            write(*,111)' D0ms2ir2_coli: p23 = ',p23
            write(*,111)' D0ms2ir2_coli: p34 = ',p34
            write(*,111)' D0ms2ir2_coli: p14 = ',p14
            write(*,111)' D0ms2ir2_coli: p24 = ',p24
            write(*,111)' D0ms2ir2_coli: p13 = ',p13
            write(*,111)' D0ms2ir2_coli: m12 = ',m12
            write(*,111)' D0ms2ir2_coli: m22 = ',m22
            write(*,111)' D0ms2ir2_coli: m32 = ',m32
            write(*,111)' D0ms2ir2_coli: m42 = ',m42
            flag(1)=.false.
          endif
#endif

c Beenakker et al., NPB653 (2003) 151 
c Ellis, Ross, Terrano, NPB 178 (1981) 421              box7 
c m^2 -> m42, p4^2 -> p34, s12 -> p13, s23 -> p24
c Denner Dittmaier (3.33)

          D0ms2ir2_coli =
#ifdef SING
     &        1.5d0*delta2ir
     &        -delta1ir*(2d0*cln_coli((m42-p24)/m42,-1d0)
     &        +cln_coli(-p13/m42,-1d0)
     &        -cln_coli((m42-p34)/m42,-1d0)
     &        +1.5d0*log(m42/muir2))
#endif
     &        +log(m42/muir2)*(2d0*cln_coli((m42-p24)/m42,-1d0)
     &        +cln_coli(-p13/m42,-1d0)
     &        -cln_coli((m42-p34)/m42,-1d0)+.75d0*log(m42/muir2))
     &        +2d0*cln_coli(-p13/m42,-1d0)
     &        *cln_coli((m42-p24)/m42,-1d0)
     &        -(cln_coli((m42-p34)/m42,-1d0))**2
     &        -2d0*
     &        cspenc_coli((p34-p24)/(m42-p24),dreal(p34-p24))
     &        -4d0*pi2_6
 
          D0ms2ir2_coli = D0ms2ir2_coli/(p13*(p24-m42))

        else
* 2 soft singularities and 2 zero masses
*
*                   0             
*          0 ----------------  m32  
*                |  2   |         
*              0 |1    3| m32     
*                |   4  |         
*        m42 ----------------  p34
*                  m42             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.0d0.or.m42.eq.cd0.or.m32.eq.cd0.or.
     &        p12.ne.0d0.or.p34.eq.cd0.or.p14.ne.m42.or.p23.ne.m32)
     &        then
            write(*,100)' D0ms2ir2_coli:  (3.35)',
     &          '    case with 2 soft singularities',
     &          '    and 2 zero masses: wrong arguments'
            write(*,111)' D0ms2ir2_coli: p12 = ',p12
            write(*,111)' D0ms2ir2_coli: p23 = ',p23
            write(*,111)' D0ms2ir2_coli: p34 = ',p34
            write(*,111)' D0ms2ir2_coli: p14 = ',p14
            write(*,111)' D0ms2ir2_coli: p24 = ',p24
            write(*,111)' D0ms2ir2_coli: p13 = ',p13
            write(*,111)' D0ms2ir2_coli: m12 = ',m12
            write(*,111)' D0ms2ir2_coli: m22 = ',m22
            write(*,111)' D0ms2ir2_coli: m32 = ',m32
            write(*,111)' D0ms2ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(2))then
            write(*,100)' D0ms2ir2_coli:  (3.35)',
     &          '    case with 2 soft singularities',
     &          '    and 2 zero masses not yet tested at all'
            write(*,111)' D0ms2ir2_coli: p12 = ',p12
            write(*,111)' D0ms2ir2_coli: p23 = ',p23
            write(*,111)' D0ms2ir2_coli: p34 = ',p34
            write(*,111)' D0ms2ir2_coli: p14 = ',p14
            write(*,111)' D0ms2ir2_coli: p24 = ',p24
            write(*,111)' D0ms2ir2_coli: p13 = ',p13
            write(*,111)' D0ms2ir2_coli: m12 = ',m12
            write(*,111)' D0ms2ir2_coli: m22 = ',m22
            write(*,111)' D0ms2ir2_coli: m32 = ',m32
            write(*,111)' D0ms2ir2_coli: m42 = ',m42
            flag(2)=.false.
          endif
#endif
         
c Beenakker et al., NPB653 (2003) 151 
c Ellis, Ross, Terrano, NPB 178 (1981) 421              box11  m3 small
c Denner Dittmaier (3.35)

          D0ms2ir2_coli =
#ifdef SING
     &        delta2ir
     &        -delta1ir*
     &        (cln_coli(-p13/sqrt(muir2*m32*coliminfscale2),-1d0)
     &        +cln_coli((m42-p24)/sqrt(muir2*m42*coliminfscale2),-1d0))
     &        -0.25d0*colishiftms2
#endif
     &        +2d0*cln_coli(-p13/sqrt(muir2*m32*coliminfscale2),-1d0)
     &        *cln_coli((m42-p24)/sqrt(muir2*m42*coliminfscale2),-1d0)
     &        -(cln_coli((m42-p34)/sqrt(m42*m32*coliminfscale2),-1d0))
     &        **2
     &        -4d0*pi2_6
 
          D0ms2ir2_coli = D0ms2ir2_coli/(p13*(p24-m42))

        endif

      elseif(nsoft.eq.1)then

        if(m12.ne.cd0)then
          if(m32.ne.cd0)then

* 1 soft singularity and one zero mass opposite to finite mass
*
*                   0             
*        m12 ----------------  m32  
*                |  2   |         
*            m12 |1    3| m32       
*                |   4  |         
*        m42 ----------------  p34  
*                   m42             
*

#ifdef CHECK
            if(m12.eq.cd0.or.m22.ne.cd0.or.m32.eq.cd0.or.m42.eq.cd0.or.
     &          p12.ne.m12.or.p23.ne.m32.or.p14.ne.m42.or.p34.eq.cd0)
     &          then
              write(*,100)' D0ms2ir2_coli:  (3.40)',
     &            '    case with 1 soft singularity ',
     &            '    and 1 zero mass opposite to finite mass: ',
     &            '    wrong arguments'
              write(*,111)' D0ms2ir2_coli: p12 = ',p12
              write(*,111)' D0ms2ir2_coli: p23 = ',p23
              write(*,111)' D0ms2ir2_coli: p34 = ',p34
              write(*,111)' D0ms2ir2_coli: p14 = ',p14
              write(*,111)' D0ms2ir2_coli: p24 = ',p24
              write(*,111)' D0ms2ir2_coli: p13 = ',p13
              write(*,111)' D0ms2ir2_coli: m12 = ',m12
              write(*,111)' D0ms2ir2_coli: m22 = ',m22
              write(*,111)' D0ms2ir2_coli: m32 = ',m32
              write(*,111)' D0ms2ir2_coli: m42 = ',m42
              stop
            endif
#endif
          
#ifdef CHECK
            if(flag(3))then
              write(*,100)' D0ms2ir2_coli:  (3.40)',
     &          '    case with 1 soft singularity ',
     &            '    and 1 zero mass opposite to finite mass ',
     &            '    not yet implemented'
              write(*,111)' D0ms2ir2_coli: p12 = ',p12
              write(*,111)' D0ms2ir2_coli: p23 = ',p23
              write(*,111)' D0ms2ir2_coli: p34 = ',p34
              write(*,111)' D0ms2ir2_coli: p14 = ',p14
              write(*,111)' D0ms2ir2_coli: p24 = ',p24
              write(*,111)' D0ms2ir2_coli: p13 = ',p13
              write(*,111)' D0ms2ir2_coli: m12 = ',m12
              write(*,111)' D0ms2ir2_coli: m22 = ',m22
              write(*,111)' D0ms2ir2_coli: m32 = ',m32
              write(*,111)' D0ms2ir2_coli: m42 = ',m42
              flag(3)=.false.
            endif
#endif 
                   
c according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 (i)
c              m0=m3, m2,m4 small
c Ellis, Zanderighi, JHEP 0802 (2008) 002  
c              box16 m22,m42 small,  p22=m32                   
c Denner Dittmaier (3.40)

            D0ms2ir2_coli=cd0
#ifdef SING
#endif

            D0ms2ir2_coli=D0ms2ir2_coli/(p13*(p24-m42))
          
          else   ! m32=0
          
* 1 soft singularity and one zero mass opposite to finite mass
*
*                   0             
*        m12 ----------------  0  
*                |  2   |         
*            m12 |1    3| 0       
*                |   4  |         
*        m42 ----------------  p34  
*                   m42             
*

#ifdef CHECK
            if(m12.eq.cd0.or.m22.ne.cd0.or.m32.ne.cd0.or.m42.eq.cd0.or.
     &          p12.ne.m12.or.p23.ne.cd0.or.p14.ne.m42.or.p34.eq.cd0)
     &          then
              write(*,100)' D0ms2ir2_coli:  (3.34)',
     &          '    case with 1 soft singularity ',
     &            '    and two adjacent zero masses: wrong arguments'
              write(*,111)' D0ms2ir2_coli: p12 = ',p12
              write(*,111)' D0ms2ir2_coli: p23 = ',p23
              write(*,111)' D0ms2ir2_coli: p34 = ',p34
              write(*,111)' D0ms2ir2_coli: p14 = ',p14
              write(*,111)' D0ms2ir2_coli: p24 = ',p24
              write(*,111)' D0ms2ir2_coli: p13 = ',p13
              write(*,111)' D0ms2ir2_coli: m12 = ',m12
              write(*,111)' D0ms2ir2_coli: m22 = ',m22
              write(*,111)' D0ms2ir2_coli: m32 = ',m32
              write(*,111)' D0ms2ir2_coli: m42 = ',m42
              stop
            endif
#endif
          
#ifdef CHECK
            if(flag(4))then
              write(*,100)' D0ms2ir2_coli:  (3.34)',
     &          '    case with 1 soft singularity ',
     &            '    and 2 adjacent zero  masses not yet implemented'
              write(*,111)' D0ms2ir2_coli: p12 = ',p12
              write(*,111)' D0ms2ir2_coli: p23 = ',p23
              write(*,111)' D0ms2ir2_coli: p34 = ',p34
              write(*,111)' D0ms2ir2_coli: p14 = ',p14
              write(*,111)' D0ms2ir2_coli: p24 = ',p24
              write(*,111)' D0ms2ir2_coli: p13 = ',p13
              write(*,111)' D0ms2ir2_coli: m12 = ',m12
              write(*,111)' D0ms2ir2_coli: m22 = ',m22
              write(*,111)' D0ms2ir2_coli: m32 = ',m32
              write(*,111)' D0ms2ir2_coli: m42 = ',m42
              flag(4)=.false.
            endif
#endif
                    
c new case
c Denner Dittmaier (3.34)

            D0ms2ir2_coli=cd0
#ifdef SING
#endif
            
            D0ms2ir2_coli=D0ms2ir2_coli/(p13*(p24-m42))
          
          endif
        else !   if(m12.eq.cd0)then
          if(m32.ne.cd0)then

* 1 soft singularity and one zero mass opposite to finite mass
*
*                  m12             
*          0 ----------------  m12  
*                |  2   |         
*            m12 |1    3| 0       
*                |   4  |         
*        m42 ----------------  p34  
*                   m42             
*

#ifdef CHECK
            if(m12.eq.cd0.or.m22.ne.m12.or.m32.ne.cd0.or.m42.eq.cd0.or.
     &          p12.ne.cd0.or.p23.ne.m12.or.p14.ne.m42.or.p34.eq.cd0)
     &          then
              write(*,100)' D0ms2ir2_coli:  (3.40)',
     &          '    case with 1 soft singularity ',
     &            '    and 2 adjacent equal nonzero masses: ',
     &            '    wrong arguments'
              write(*,111)' D0ms2ir2_coli: p12 = ',p12
              write(*,111)' D0ms2ir2_coli: p23 = ',p23
              write(*,111)' D0ms2ir2_coli: p34 = ',p34
              write(*,111)' D0ms2ir2_coli: p14 = ',p14
              write(*,111)' D0ms2ir2_coli: p24 = ',p24
              write(*,111)' D0ms2ir2_coli: p13 = ',p13
              write(*,111)' D0ms2ir2_coli: m12 = ',m12
              write(*,111)' D0ms2ir2_coli: m22 = ',m22
              write(*,111)' D0ms2ir2_coli: m32 = ',m32
              write(*,111)' D0ms2ir2_coli: m42 = ',m42
              stop
            endif
#endif
          
#ifdef CHECK
            if(flag(5))then
              write(*,100)' D0ms2ir2_coli:  (3.40)',
     &          '    case with 1 soft singularity ',
     &            '    and 2 adjacent equal nonzero masses ',
     &            '    not yet implemented'
              write(*,111)' D0ms2ir2_coli: p12 = ',p12
              write(*,111)' D0ms2ir2_coli: p23 = ',p23
              write(*,111)' D0ms2ir2_coli: p34 = ',p34
              write(*,111)' D0ms2ir2_coli: p14 = ',p14
              write(*,111)' D0ms2ir2_coli: p24 = ',p24
              write(*,111)' D0ms2ir2_coli: p13 = ',p13
              write(*,111)' D0ms2ir2_coli: m12 = ',m12
              write(*,111)' D0ms2ir2_coli: m22 = ',m22
              write(*,111)' D0ms2ir2_coli: m32 = ',m32
              write(*,111)' D0ms2ir2_coli: m42 = ',m42
              flag(5)=.false.
            endif
#endif
                    
c according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 (ic)
c              m0=m3, m2,m4 small
c Ellis, Zanderighi, JHEP 0802 (2008) 002  
c              box16 m22,m42 small,  p22=m32                   
c Denner Dittmaier (3.40)

            D0ms2ir2_coli = (
     &          -2d0*cln_coli((m42-p24)/(sqrt(m22*coliminfscale2*m42)),
     &          -1d0)
     &          *cln_coli((-p13)/sqrt(m22*coliminfscale2*muir2),-1d0)
     &          +(cln_coli((m42-p34)/sqrt(m22*coliminfscale2*m42),-1d0))
     &          **2 
     &          + pi2_6 
     &          + 2d0*cspenc_coli((p24-p34)/(p24-m42),dreal(p34-m42)) ) 
     &          /((p24-m42)*(m32-p13))
            
#ifdef SING
     &          +cln_coli((m42-p24)/(sqrt(m22*coliminfscale2*m42)),-1d0)
     &          *delta1ir
     &          - 0.25d0*colishiftms2
#endif

            D0ms2ir2_coli=D0ms2ir2_coli/(p13*(p24-m42))
          
          else   ! m32=0, m12=0
          
* 1 soft singularity and one zero mass opposite to finite mass
*
*                  m22             
*        m22 ----------------  m22  
*                |  2   |         
*              0 |1    3| 0       
*                |   4  |         
*        m42 ----------------  p34  
*                   m42             
*

#ifdef CHECK
            if(m12.ne.cd0.or.m22.eq.cd0.or.m32.ne.cd0.or.m42.eq.cd0.or.
     &          p12.ne.m22.or.p23.ne.m22.or.p14.ne.m42.or.p34.eq.cd0)
     &          then
              write(*,100)' D0ms2ir2_coli:  (3.38)',
     &          '    case with 1 soft singularity ',
     &            '    and two opposite zero masses: wrong arguments'
              write(*,111)' D0ms2ir2_coli: p12 = ',p12
              write(*,111)' D0ms2ir2_coli: p23 = ',p23
              write(*,111)' D0ms2ir2_coli: p34 = ',p34
              write(*,111)' D0ms2ir2_coli: p14 = ',p14
              write(*,111)' D0ms2ir2_coli: p24 = ',p24
              write(*,111)' D0ms2ir2_coli: p13 = ',p13
              write(*,111)' D0ms2ir2_coli: m12 = ',m12
              write(*,111)' D0ms2ir2_coli: m22 = ',m22
              write(*,111)' D0ms2ir2_coli: m32 = ',m32
              write(*,111)' D0ms2ir2_coli: m42 = ',m42
              stop
            endif
#endif
          
#ifdef CHECK
            if(flag(6))then
              write(*,100)' D0ms2ir2_coli:  (3.38)',
     &          '    case with 1 soft singularity ',
     &            '    and 2 opposite zero masses not yet tested at all'
              write(*,111)' D0ms2ir2_coli: p12 = ',p12
              write(*,111)' D0ms2ir2_coli: p23 = ',p23
              write(*,111)' D0ms2ir2_coli: p34 = ',p34
              write(*,111)' D0ms2ir2_coli: p14 = ',p14
              write(*,111)' D0ms2ir2_coli: p24 = ',p24
              write(*,111)' D0ms2ir2_coli: p13 = ',p13
              write(*,111)' D0ms2ir2_coli: m12 = ',m12
              write(*,111)' D0ms2ir2_coli: m22 = ',m22
              write(*,111)' D0ms2ir2_coli: m32 = ',m32
              write(*,111)' D0ms2ir2_coli: m42 = ',m42
              flag(6)=.false.
            endif
#endif
                    
c according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 (iic)
c Ellis, Zanderighi, JHEP 0802 (2008) 002  
c              box15 m22 small,  p22=m22                   
c Denner Dittmaier (3.38)

            D0ms2ir2_coli = 
     &          cln_coli((m42-p24)/sqrt(m22*coliminfscale2*m42),-1d0)*
     &          (cln_coli((m42-p24)/sqrt(m22*coliminfscale2*m42),-1d0)
     &          -log(muir2/m42)
     &          +2d0*cln_coli(p13/(p34-m42),dreal(p34-m42)) )
     &          + pi2_6 
#ifdef SING
     &          -cln_coli((m42-p24)/sqrt(m22*coliminfscale2*m42),-1d0)*
     &          delta1ir
     &          + 0.25d0*colishiftms2
#endif
            D0ms2ir2_coli = D0ms2ir2_coli/((p24-m42)*p13)

          endif
        endif

      elseif(nsoft.eq.0)then

        if(m32.eq.cd0)then
* no soft singularity and 1 zero mass
*
*                m22 small            
*          0 ----------------  m22  
*                |  2   |         
*            m22 |1    3| 0       
*                |   4  |         
*        m42 ----------------  p34  
*                  m42             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.ne.m12.or.m42.eq.cd0.or.
     &        p34.eq.cd0.or.p23.ne.m12.or.p14.ne.m42.or.p12.ne.cd0)
     &        then
            write(*,100)' D0ms2ir2_coli:  (3.36)',
     &          '    case with 0 soft singularity ',
     &          '    and 1 zero mass: wrong arguments'
            write(*,111)' D0ms2ir2_coli: p12 = ',p12
            write(*,111)' D0ms2ir2_coli: p23 = ',p23
            write(*,111)' D0ms2ir2_coli: p34 = ',p34
            write(*,111)' D0ms2ir2_coli: p14 = ',p14
            write(*,111)' D0ms2ir2_coli: p24 = ',p24
            write(*,111)' D0ms2ir2_coli: p13 = ',p13
            write(*,111)' D0ms2ir2_coli: m12 = ',m12
            write(*,111)' D0ms2ir2_coli: m22 = ',m22
            write(*,111)' D0ms2ir2_coli: m32 = ',m32
            write(*,111)' D0ms2ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(7))then
            write(*,100)' D0ms2ir2_coli:  (3.36)',
     &          '    case with 0 soft singularity ',
     &          '    and 1 zero mass not yet implemented'
            write(*,111)' D0ms2ir2_coli: p12 = ',p12
            write(*,111)' D0ms2ir2_coli: p23 = ',p23
            write(*,111)' D0ms2ir2_coli: p34 = ',p34
            write(*,111)' D0ms2ir2_coli: p14 = ',p14
            write(*,111)' D0ms2ir2_coli: p24 = ',p24
            write(*,111)' D0ms2ir2_coli: p13 = ',p13
            write(*,111)' D0ms2ir2_coli: m12 = ',m12
            write(*,111)' D0ms2ir2_coli: m22 = ',m22
            write(*,111)' D0ms2ir2_coli: m32 = ',m32
            write(*,111)' D0ms2ir2_coli: m42 = ',m42
            flag(7)=.false.
            stop
          endif
#endif
          
c new case
c Denner Dittmaier (3.36)

          D0ms2ir2_coli=0d0
#ifdef SING
#endif

          D0ms2ir2_coli=D0ms2ir2_coli/(p13*p24)
          
        else ! m32.ne.0
* 0 soft singularity and 1 zero mass
*
*                   m22 small            
*          0 ----------------  0  
*                |  2   |         
*            m22 |1    3| m22       
*                |   4  |         
*        m42 ----------------  p34  
*                   m42             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.ne.m12.or.m32.ne.m12.or.m42.eq.cd0.or.
     &        p12.ne.cd0.or.p23.ne.cd0.or.p14.ne.m42.or.p34.eq.cd0)
     &        then
            write(*,100)' D0ms2ir2_coli:  (3.39)',
     &          '    case with 0 soft singularity ',
     &          '    and 0 zero mass: wrong arguments'
            write(*,111)' D0ms2ir2_coli: p12 = ',p12
            write(*,111)' D0ms2ir2_coli: p23 = ',p23
            write(*,111)' D0ms2ir2_coli: p34 = ',p34
            write(*,111)' D0ms2ir2_coli: p14 = ',p14
            write(*,111)' D0ms2ir2_coli: p24 = ',p24
            write(*,111)' D0ms2ir2_coli: p13 = ',p13
            write(*,111)' D0ms2ir2_coli: m12 = ',m12
            write(*,111)' D0ms2ir2_coli: m22 = ',m22
            write(*,111)' D0ms2ir2_coli: m32 = ',m32
            write(*,111)' D0ms2ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(8))then
            write(*,100)' D0ms2ir2_coli:  (3.39)',
     &          '    case with 0 soft singularity ',
     &          '    and 0 zero mass not yet implemented'
            write(*,111)' D0ms2ir2_coli: p12 = ',p12
            write(*,111)' D0ms2ir2_coli: p23 = ',p23
            write(*,111)' D0ms2ir2_coli: p34 = ',p34
            write(*,111)' D0ms2ir2_coli: p14 = ',p14
            write(*,111)' D0ms2ir2_coli: p24 = ',p24
            write(*,111)' D0ms2ir2_coli: p13 = ',p13
            write(*,111)' D0ms2ir2_coli: m12 = ',m12
            write(*,111)' D0ms2ir2_coli: m22 = ',m22
            write(*,111)' D0ms2ir2_coli: m32 = ',m32
            write(*,111)' D0ms2ir2_coli: m42 = ',m42
            flag(8)=.false.
          endif
#endif
          
c new case
c Denner Dittmaier (3.39)

          D0ms2ir2_coli=0d0

          D0ms2ir2_coli=D0ms2ir2_coli/(p13*p24)
          
        endif

#ifdef CHECK
      else
        write(*,100)' D0ms2ir2_coli: case with more than 2 soft ', 
     &      'singularities not implemented '
        write(*,111)' D0ms2ir2_coli: nsoft = ',nsoft
        write(*,111)' D0ms2ir2_coli: p23 = ',p23
        write(*,111)' D0ms2ir2_coli: p34 = ',p34
        write(*,111)' D0ms2ir2_coli: p14 = ',p14
        write(*,111)' D0ms2ir2_coli: p24 = ',p24
        write(*,111)' D0ms2ir2_coli: p13 = ',p13
        write(*,111)' D0ms2ir2_coli: m12 = ',m12
        write(*,111)' D0ms2ir2_coli: m22 = ',m22
        write(*,111)' D0ms2ir2_coli: m32 = ',m32
        write(*,111)' D0ms2ir2_coli: m42 = ',m42    
        stop
#endif
      endif

      end



************************************************************************
      function D0ms2ir3_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     mass-singular 4-point function                                   *
*     with two collinear and up to 3 soft singularities                *
*     p12, p23, m12, m22, m32   small                                  *
*     p14=m42, p34, p13, p24    finite                                 *
*                                                                      *
*                     m22 small                                        *
*     p12 small  ---------------------  p23 small                      *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | m32 small                            *
*                     |         |                                      *
*                     |    4    |                                      *
*       p14=m42  ---------------------  p34                            *
*                         m42                                          *
*                                                                      *
************************************************************************
*     22.08.08 Ansgar Denner        last changed  22.08.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms2ir3_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      logical flag(0:6)
      data    flag /7*.true./
      save    flag  
#endif

      complex*16 q13,q24,mm32,m2(4)
      complex*16 cln_coli
      integer    nsoft
      integer    i,j,k
      logical    soft(4),onsh(4,4)

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
        ms42 = elimminf2_coli(m42)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps34 = elimminf2_coli(p34)
        ps14 = elimminf2_coli(p14)
        ps24 = elimminf2_coli(p24)
        ps13 = elimminf2_coli(p13)

        if(ms12.ne.cd0.or.ms22.ne.cd0.or.ms32.ne.cd0
     &        .or.ms42.eq.cd0.or.ps12.ne.cd0.or.ps23.ne.cd0 
     &        .or.ps34.ne.ms42.or.ps14.ne.ms42
     &        .or.ps24.eq.cd0.or.ps13.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &        .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0) then
          write(*,111)' D0ms2ir3_coli called improperly:'
          write(*,111)' D0ms2ir3_coli: p12 = ',p12,ps12
          write(*,111)' D0ms2ir3_coli: p23 = ',p23,ps23
          write(*,111)' D0ms2ir3_coli: p34 = ',p34,ps34
          write(*,111)' D0ms2ir3_coli: p14 = ',p14,ps14
          write(*,111)' D0ms2ir3_coli: p13 = ',p13,ps13
          write(*,111)' D0ms2ir3_coli: p24 = ',p24,ps24
          write(*,111)' D0ms2ir3_coli: m12 = ',m12,ms12
          write(*,111)' D0ms2ir3_coli: m22 = ',m22,ms22
          write(*,111)' D0ms2ir3_coli: m32 = ',m32,ms32
          write(*,111)' D0ms2ir3_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif

      m2(1)=m12
      m2(2)=m22
      m2(3)=m32
      m2(4)=m42

c determine on-shell momenta
      onsh(1,2)=p12.eq.m22
      onsh(1,3)=p13.eq.m32
      onsh(1,4)=p14.eq.m42
      onsh(2,1)=p12.eq.m12
      onsh(2,3)=p23.eq.m32
      onsh(2,4)=p24.eq.m42
      onsh(3,1)=p13.eq.m12
      onsh(3,2)=p23.eq.m22
      onsh(3,4)=p34.eq.m42
      onsh(4,1)=p14.eq.m22
      onsh(4,2)=p24.eq.m22
      onsh(4,3)=p34.eq.m32

c count/determine soft singularities
      nsoft=0 
      do i=1,4
        do j=1,3
          if(i.ne.j)then
            do k=j,4
              if(k.ne.i.and.k.ne.j)then
                soft(i)=m2(i).eq.cd0.and.onsh(i,j).and.onsh(i,k)
c                soft(i,k,j)=soft(i,j,k)
                if(soft(i)) nsoft=nsoft+1
              endif
            enddo
          endif
        enddo
      enddo

      if(nsoft.eq.3)then

* four soft singularities
*
*                   0             
*          0 ----------------  0  
*                |  2   |         
*              0 |1    3| 0       
*                |   4  |         
*        m42 ----------------  m42  
*                   m42             
*

#ifdef CHECK
        if(m12.ne.cd0.or.m22.ne.0d0.or.m32.ne.0d0.or.m42.eq.cd0
     &      .or.p12.ne.cd0.or.p23.ne.cd0.or.p34.ne.m42.or.p14.ne.m42)
     &      then
          write(*,100)' D0ms2ir3_coli:  (3.42)',
     &          '    case with 3 soft singularities: ',
     &        '    wrong arguments'
          write(*,111)' D0ms2ir3_coli: p12 = ',p12
          write(*,111)' D0ms2ir3_coli: p23 = ',p23
          write(*,111)' D0ms2ir3_coli: p34 = ',p34
          write(*,111)' D0ms2ir3_coli: p14 = ',p14
          write(*,111)' D0ms2ir3_coli: p24 = ',p24
          write(*,111)' D0ms2ir3_coli: p13 = ',p13
          write(*,111)' D0ms2ir3_coli: m12 = ',m12
          write(*,111)' D0ms2ir3_coli: m22 = ',m22
          write(*,111)' D0ms2ir3_coli: m32 = ',m32
          write(*,111)' D0ms2ir3_coli: m42 = ',m42
          stop
        endif
#endif

#ifdef CHECK
        if(flag(1))then
          write(*,100)' D0ms2ir3_coli:  (3.42)',
     &          '    case with 3 soft singularities',
     &        '    not yet tested at all'
          write(*,111)' D0ms2ir3_coli: p12 = ',p12
          write(*,111)' D0ms2ir3_coli: p23 = ',p23
          write(*,111)' D0ms2ir3_coli: p34 = ',p34
          write(*,111)' D0ms2ir3_coli: p14 = ',p14
          write(*,111)' D0ms2ir3_coli: p24 = ',p24
          write(*,111)' D0ms2ir3_coli: p13 = ',p13
          write(*,111)' D0ms2ir3_coli: m12 = ',m12
          write(*,111)' D0ms2ir3_coli: m22 = ',m22
          write(*,111)' D0ms2ir3_coli: m32 = ',m32
          write(*,111)' D0ms2ir3_coli: m42 = ',m42
          flag(1)=.false.
        endif
#endif

c Beenakker, Kuiff, van Neerven, Smith, PRD 40 (1989) 54 
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002    box6
c Denner, Dittmaier (3.42)
        D0ms2ir3_coli = 
#ifdef SING
     &      2d0*delta2ir
     &      -delta1ir*(2d0*cln_coli((m42-p13)/sqrt(m42*muir2),-1d0)
     &                 +cln_coli(-p24/muir2,-1d0)) 
#endif
     &      +2d0*cln_coli((m42--p13)/sqrt(m42*muir2),-1d0)
     &      *cln_coli(-p24/muir2,-1d0)
     &      -5d0*pi2_6

        D0ms2ir3_coli = D0ms2ir3_coli/(p13*p24)

      elseif(nsoft.eq.2)then

        if(m12.eq.cd0.and.m32.eq.cd0)then
* 2 soft singularities and two opposite zero masses
*
*                m22 small            
*        m22  ---------------- m22  
*                |  2   |         
*             0  |1    3| 0       
*                |   4  |         
*        m42 ----------------  m42  
*                  m42             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.eq.0d0.or.m32.ne.0d0.or.m42.eq.cd0
     &        .or.p12.ne.m22.or.p23.ne.m22.or.p34.ne.m42.or.p14.ne.m42)
     &        then
            write(*,100)' D0ms2ir3_coli:  (3.45)',
     &          '    case with 2 soft singularities',
     &          '    and 2 opposite zero masses: wrong arguments'
            write(*,111)' D0ms2ir3_coli: p12 = ',p12
            write(*,111)' D0ms2ir3_coli: p23 = ',p23
            write(*,111)' D0ms2ir3_coli: p34 = ',p34
            write(*,111)' D0ms2ir3_coli: p14 = ',p14
            write(*,111)' D0ms2ir3_coli: p24 = ',p24
            write(*,111)' D0ms2ir3_coli: p13 = ',p13
            write(*,111)' D0ms2ir3_coli: m12 = ',m12
            write(*,111)' D0ms2ir3_coli: m22 = ',m22
            write(*,111)' D0ms2ir3_coli: m32 = ',m32
            write(*,111)' D0ms2ir3_coli: m42 = ',m42
            stop
          endif
#endif

#ifdef CHECK
          if(flag(2))then
            write(*,100)' D0ms2ir3_coli:  (3.45)',
     &          '    case with 2 soft singularities',
     &          '    and 2 opposite zero masses not yet tested at all'
            write(*,111)' D0ms2ir3_coli: p12 = ',p12
            write(*,111)' D0ms2ir3_coli: p23 = ',p23
            write(*,111)' D0ms2ir3_coli: p34 = ',p34
            write(*,111)' D0ms2ir3_coli: p14 = ',p14
            write(*,111)' D0ms2ir3_coli: p24 = ',p24
            write(*,111)' D0ms2ir3_coli: p13 = ',p13
            write(*,111)' D0ms2ir3_coli: m12 = ',m12
            write(*,111)' D0ms2ir3_coli: m22 = ',m22
            write(*,111)' D0ms2ir3_coli: m32 = ',m32
            write(*,111)' D0ms2ir3_coli: m42 = ',m42
            flag(2)=.false.
          endif
#endif

* according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349  (iiia)
c Ellis, Zanderighi, JHEP 0802 (2008) 002    box14 m2 small
c Denner, Dittmaier (3.45)
          D0ms2ir3_coli =
#ifdef SING
     &        2d0*cln_coli((m42-q24)/(sqrt(m22*m42)*coliminfscale),-1d0)
     &        *delta1ir
#endif
     &        -2d0*cln_coli((m42-q24)/(sqrt(m22*m42)*coliminfscale),
     &        -1d0)
     &        *cln_coli(-q13/muir2,-1d0)
          D0ms2ir3_coli = D0ms2ir3_coli/(q13*(m42-q24))

        else
* 2 soft singularities and two adjacent zero masses
*
*                   0            
*          0  ---------------- m32  
*                |  2   |         
*             0  |1    3| m32 small       
*                |   4  |         
*        m42 ----------------  m42  
*                  m42             
*

#ifdef CHECK
          if(m22.ne.cd0.or.m42.eq.cd0.or.
     &        .not.(m32.eq.cd0.and.m12.ne.cd0.or.
     &        m32.ne.cd0.and.m12.eq.cd0).or.
     &        p12.ne.m12.or.p23.ne.m32.or.p34.ne.m42.or.p14.ne.m42)
     &        then
            write(*,100)' D0ms2ir3_coli:  (3.43)',
     &          '    case with 2 soft singularities',
     &          '    and two adjacent zero masses: wrong arguments'
            write(*,111)' D0ms2ir3_coli: p12 = ',p12
            write(*,111)' D0ms2ir3_coli: p23 = ',p23
            write(*,111)' D0ms2ir3_coli: p34 = ',p34
            write(*,111)' D0ms2ir3_coli: p14 = ',p14
            write(*,111)' D0ms2ir3_coli: p24 = ',p24
            write(*,111)' D0ms2ir3_coli: p13 = ',p13
            write(*,111)' D0ms2ir3_coli: m12 = ',m12
            write(*,111)' D0ms2ir3_coli: m22 = ',m22
            write(*,111)' D0ms2ir3_coli: m32 = ',m32
            write(*,111)' D0ms2ir3_coli: m42 = ',m42
            stop
          endif
#endif

#ifdef CHECK
          if(flag(3))then
            write(*,100)' D0ms2ir3_coli:  (3.43)',
     &          '    case with 2 soft singularities',
     &          '    and 2 adjacent zero masses not yet tested at all'
            write(*,111)' D0ms2ir3_coli: p12 = ',p12
            write(*,111)' D0ms2ir3_coli: p23 = ',p23
            write(*,111)' D0ms2ir3_coli: p34 = ',p34
            write(*,111)' D0ms2ir3_coli: p14 = ',p14
            write(*,111)' D0ms2ir3_coli: p24 = ',p24
            write(*,111)' D0ms2ir3_coli: p13 = ',p13
            write(*,111)' D0ms2ir3_coli: m12 = ',m12
            write(*,111)' D0ms2ir3_coli: m22 = ',m22
            write(*,111)' D0ms2ir3_coli: m32 = ',m32
            write(*,111)' D0ms2ir3_coli: m42 = ',m42
            flag(3)=.false.
          endif
#endif

          if(m12.ne.cd0)then 
            mm32=m12*coliminfscale2
          else
            mm32=m32*coliminfscale2
          endif
         
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002    
c        box11 p32=m42, m3 small  
c Denner, Dittmaier (3.43)
          D0ms2ir3_coli =
#ifdef SING
     &        delta2ir
     &        -delta1ir*(cln_coli(-p13/sqrt(muir2*mm32),-1d0)
     &        +cln_coli((m42-p24)/sqrt(muir2*m42),-1d0))
#endif
     &        +2d0*cln_coli(-p13/sqrt(muir2*mm32),-1d0)
     &        *cln_coli((m42-p24)/sqrt(muir2*m42),-1d0)
     &        -2.5d0*pi2_6
c constant std  -3*pi2/4
 
          D0ms2ir3_coli = D0ms2ir3_coli/(p13*(p24-m42))

        endif

      elseif(nsoft.eq.1)then

        if(m12.ne.cd0.and.m32.ne.cd0)then
* 1 soft singularities opposite to finite mass m4
*
*                   0            
*        m12  ---------------- m32  
*                |  2   |         
*           m12  |1    3| m32       
*                |   4  |         
*        m42 ----------------  m42  
*                  m42             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.ne.0d0.or.m32.eq.0d0.or.m42.eq.cd0
     &        .or.p12.ne.m12.or.p23.ne.m32.or.p34.ne.m42.or.p14.ne.m42)
     &        then
            write(*,100)' D0ms2ir3_coli:  (3.47)',
     &          '    case with 1 soft singularity ',
     &          ' opposite to finite mass: wrong arguments'
            write(*,111)' D0ms2ir3_coli: p12 = ',p12
            write(*,111)' D0ms2ir3_coli: p23 = ',p23
            write(*,111)' D0ms2ir3_coli: p34 = ',p34
            write(*,111)' D0ms2ir3_coli: p14 = ',p14
            write(*,111)' D0ms2ir3_coli: p24 = ',p24
            write(*,111)' D0ms2ir3_coli: p13 = ',p13
            write(*,111)' D0ms2ir3_coli: m12 = ',m12
            write(*,111)' D0ms2ir3_coli: m22 = ',m22
            write(*,111)' D0ms2ir3_coli: m32 = ',m32
            write(*,111)' D0ms2ir3_coli: m42 = ',m42
            stop
          endif
#endif

#ifdef CHECK
          if(flag(4))then
            write(*,100)' D0ms2ir3_coli:  (3.47)',
     &          '    case with 1 soft singularity ',
     &          'opposite to finite mass not yet implemented'
            write(*,111)' D0ms2ir3_coli: p12 = ',p12
            write(*,111)' D0ms2ir3_coli: p23 = ',p23
            write(*,111)' D0ms2ir3_coli: p34 = ',p34
            write(*,111)' D0ms2ir3_coli: p14 = ',p14
            write(*,111)' D0ms2ir3_coli: p24 = ',p24
            write(*,111)' D0ms2ir3_coli: p13 = ',p13
            write(*,111)' D0ms2ir3_coli: m12 = ',m12
            write(*,111)' D0ms2ir3_coli: m22 = ',m22
            write(*,111)' D0ms2ir3_coli: m32 = ',m32
            write(*,111)' D0ms2ir3_coli: m42 = ',m42
            flag(4)=.false.
            stop
          endif
#endif

* according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 
c        (i) with m22=m32=m02, m42,m22 small
c Ellis, Zanderighi, JHEP 0802 (2008) 002    
c        box16 m2,m4 small, p22=m32=p32
c Denner, Dittmaier (3.47)
          D0ms2ir3_coli = cd0
#ifdef SING
#endif
          D0ms2ir3_coli = D0ms2ir3_coli/(q13*(m42-q24))

        else
* 1 soft singularity adjacent to finite mass m42
*
*                m22 small            
*         0  ---------------- m22  
*                |  2   |         
*           m22  |1    3| 0       
*                |   4  |         
*        m42 ----------------  m42  
*                  m42             
*

#ifdef CHECK
          if(m22.eq.cd0.or.m42.eq.cd0.or.
     &        .not.(m32.eq.m22.and.m12.eq.cd0.and.
     &              p12.eq.m22.and.p23.eq.cd0.or.
     &              m12.eq.m22.and.m32.eq.cd0.and.
     &              p23.eq.m22.and.p12.eq.cd0).or.
     &        p34.ne.m42.or.p14.ne.m42)
     &        then
            write(*,100)' D0ms2ir3_coli:  (3.44)',
     &          '    case with 1 soft singularity ',
     &        'adjacent to finite masses: wrong arguments'
            write(*,111)' D0ms2ir3_coli: p12 = ',p12
            write(*,111)' D0ms2ir3_coli: p23 = ',p23
            write(*,111)' D0ms2ir3_coli: p34 = ',p34
            write(*,111)' D0ms2ir3_coli: p14 = ',p14
            write(*,111)' D0ms2ir3_coli: p24 = ',p24
            write(*,111)' D0ms2ir3_coli: p13 = ',p13
            write(*,111)' D0ms2ir3_coli: m12 = ',m12
            write(*,111)' D0ms2ir3_coli: m22 = ',m22
            write(*,111)' D0ms2ir3_coli: m32 = ',m32
            write(*,111)' D0ms2ir3_coli: m42 = ',m42
            stop
          endif
#endif

#ifdef CHECK
          if(flag(5))then
            write(*,100)' D0ms2ir3_coli:  (3.44)',
     &          '    case with 1 soft singularity ',
     &          'adjacent to finite mass not yet implemented'
            write(*,111)' D0ms2ir3_coli: p12 = ',p12
            write(*,111)' D0ms2ir3_coli: p23 = ',p23
            write(*,111)' D0ms2ir3_coli: p34 = ',p34
            write(*,111)' D0ms2ir3_coli: p14 = ',p14
            write(*,111)' D0ms2ir3_coli: p24 = ',p24
            write(*,111)' D0ms2ir3_coli: p13 = ',p13
            write(*,111)' D0ms2ir3_coli: m12 = ',m12
            write(*,111)' D0ms2ir3_coli: m22 = ',m22
            write(*,111)' D0ms2ir3_coli: m32 = ',m32
            write(*,111)' D0ms2ir3_coli: m42 = ',m42
            flag(5)=.false.
            stop
          endif
#endif

* according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 
c        (ic) with m22=m32=m02, m42,m22 small
c Ellis, Zanderighi, JHEP 0802 (2008) 002    
c        box16 p32=m42, m3=m2, p22=0, m2 small  
c Denner, Dittmaier (3.44)
          D0ms2ir3_coli =cd0

#ifdef SING
#endif
 
          D0ms2ir3_coli = D0ms2ir3_coli/(p13*(p24-m42))

        endif

      elseif(nsoft.eq.0)then

* no soft singularity
*
*                  m22            
*          0 ----------------  0  
*                |  2   |         
*            m22 |1    3| m22      
*                |   4  |         
*        m42 ----------------  m42  
*                  m42            
*

#ifdef CHECK
        if(m12.ne.m22.or.m12.ne.m32.or.m42.eq.cd0.or.m12.eq.cd0.or.
     &      p12.ne.cd0.or.p23.ne.cd0.or.p34.ne.m42.or.p14.ne.m42)
     &      then
          write(*,100)' D0ms2ir3_coli:  (3.46)',
     &          '    soft regular case: ',
     &        '    wrong arguments'
          write(*,111)' D0ms2ir3_coli: p12 = ',p12
          write(*,111)' D0ms2ir3_coli: p23 = ',p23
          write(*,111)' D0ms2ir3_coli: p34 = ',p34
          write(*,111)' D0ms2ir3_coli: p14 = ',p14
          write(*,111)' D0ms2ir3_coli: p24 = ',p24
          write(*,111)' D0ms2ir3_coli: p13 = ',p13
          write(*,111)' D0ms2ir3_coli: m12 = ',m12
          write(*,111)' D0ms2ir3_coli: m22 = ',m22
          write(*,111)' D0ms2ir3_coli: m32 = ',m32
          write(*,111)' D0ms2ir3_coli: m42 = ',m42
          stop
        endif
#endif

#ifdef CHECK
        if(flag(6))then
          write(*,100)' D0ms2ir3_coli:  (3.46)',
     &          '    soft regular case ',
     &        '    not yet implemented'
          write(*,111)' D0ms2ir3_coli: p12 = ',p12
          write(*,111)' D0ms2ir3_coli: p23 = ',p23
          write(*,111)' D0ms2ir3_coli: p34 = ',p34
          write(*,111)' D0ms2ir3_coli: p14 = ',p14
          write(*,111)' D0ms2ir3_coli: p24 = ',p24
          write(*,111)' D0ms2ir3_coli: p13 = ',p13
          write(*,111)' D0ms2ir3_coli: m12 = ',m12
          write(*,111)' D0ms2ir3_coli: m22 = ',m22
          write(*,111)' D0ms2ir3_coli: m32 = ',m32
          write(*,111)' D0ms2ir3_coli: m42 = ',m42
          flag(6)=.false.
          stop         
        endif
#endif

c new case 
c Denner, Dittmaier (3.46)
        D0ms2ir3_coli = 0d0

#ifdef CHECK
      else
        write(*,100)' D0ms2ir3_coli: case with more than 3 soft 
     &      singularities not implemented '
        write(*,111)' D0ms2ir3_coli: nsoft = ',nsoft
        write(*,111)' D0ms2ir3_coli: p23 = ',p23
        write(*,111)' D0ms2ir3_coli: p34 = ',p34
        write(*,111)' D0ms2ir3_coli: p14 = ',p14
        write(*,111)' D0ms2ir3_coli: p24 = ',p24
        write(*,111)' D0ms2ir3_coli: p13 = ',p13
        write(*,111)' D0ms2ir3_coli: m12 = ',m12
        write(*,111)' D0ms2ir3_coli: m22 = ',m22
        write(*,111)' D0ms2ir3_coli: m32 = ',m32
        write(*,111)' D0ms2ir3_coli: m42 = ',m42    
        stop
#endif
      endif

      end


************************************************************************
      function D0ms3ir2_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     4-point function                                                 *
*     with 3 mass singularities and 2 soft singularities               *
*                                                                      *
*     assumes                                                          *
*     p12, p23, p14, m12, m22, m32, m42 small                          *
*     p13, p24, p34 finite                                             *
*                                                                      *
*                     m22 small                                        *
*     p12 small  ---------------------  p23 small                      *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | m32 small                            *
*                     |         |                                      *
*                     |    4    |                                      *
*     p14 small  ---------------------  p34                            *
*                     m42 small                                        *
*                                                                      *
************************************************************************
*     05.08.08 Ansgar Denner        last changed  10.12.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms3ir2_coli

#ifdef CHECK
      complex*16 D0ms3ir2_check
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      complex*16 l13,l23,l24,x
      real*8     ix
      logical flag(0:10)
      data    flag /11*.true./
      save    flag  
#endif

      complex*16 q13,q24,mm12,mm22,mm42,mm32,m2(4)
      complex*16 cln_coli,cspcon_coli,cspenc_coli
      integer    nsoft
      integer    i,j,k
      logical    soft(4),onsh(4,4)

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
        ms42 = elimminf2_coli(m42)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps34 = elimminf2_coli(p34)
        ps14 = elimminf2_coli(p14)
        ps24 = elimminf2_coli(p24)
        ps13 = elimminf2_coli(p13)

        if(ms12.ne.cd0.or.ms22.ne.cd0.or.ms32.ne.cd0
     &        .or.ms42.ne.cd0.or.ps12.ne.cd0.or.ps23.ne.cd0 
     &        .or.ps34.eq.cd0.or.ps14.ne.cd0
     &        .or.ps24.eq.cd0.or.ps13.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &        .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0
     &        .or.ps13.eq.cd0.or.ps24.eq.cd0) then
          write(*,100)' D0ms3ir2_coli called improperly:'
          write(*,111)' D0ms3ir2_coli: p12 = ',p12,ps12
          write(*,111)' D0ms3ir2_coli: p23 = ',p23,ps23
          write(*,111)' D0ms3ir2_coli: p34 = ',p34,ps34
          write(*,111)' D0ms3ir2_coli: p14 = ',p14,ps14
          write(*,111)' D0ms3ir2_coli: p13 = ',p13,ps13
          write(*,111)' D0ms3ir2_coli: p24 = ',p24,ps24
          write(*,111)' D0ms3ir2_coli: m12 = ',m12,ms12
          write(*,111)' D0ms3ir2_coli: m22 = ',m22,ms22
          write(*,111)' D0ms3ir2_coli: m32 = ',m32,ms32
          write(*,111)' D0ms3ir2_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif

      m2(1)=m12
      m2(2)=m22
      m2(3)=m32
      m2(4)=m42

c determine on-shell momenta
      onsh(1,2)=p12.eq.m22
      onsh(1,3)=p13.eq.m32
      onsh(1,4)=p14.eq.m42
      onsh(2,1)=p12.eq.m12
      onsh(2,3)=p23.eq.m32
      onsh(2,4)=p24.eq.m42
      onsh(3,1)=p13.eq.m12
      onsh(3,2)=p23.eq.m22
      onsh(3,4)=p34.eq.m42
      onsh(4,1)=p14.eq.m22
      onsh(4,2)=p24.eq.m22
      onsh(4,3)=p34.eq.m32

c count/determine soft singularities
      nsoft=0 
      do i=1,4
        do j=1,3
          if(i.ne.j)then
            do k=j,4
              if(k.ne.i.and.k.ne.j)then
                soft(i)=m2(i).eq.cd0.and.onsh(i,j).and.onsh(i,k)
c                soft(i,k,j)=soft(i,j,k)
                if(soft(i)) nsoft=nsoft+1
              endif
            enddo
          endif
        enddo
      enddo

c      write(*,*) 'D0ms3ir2_coli: nsoft= ',nsoft

      if(nsoft.eq.2)then

        if(m32.eq.cd0.and.m42.eq.cd0)then
* 2 soft singularities and 4 zero masses
*
*                   0             
*          0 ----------------  0    
*                |  2   |         
*              0 |1    3| 0       
*                |   4  |         
*          0 ----------------  p34
*                   0             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.0d0.or.m32.ne.0d0.or.m42.ne.cd0
     &        .or.p12.ne.cd0.or.p23.ne.cd0.or.p34.eq.cd0.or.p14.ne.cd0)
     &        then
            write(*,100)' D0ms3ir2_coli:  (3.57)',
     &          '    case with 2 soft singularities and 4 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef ALLCHECK
c ppttbb dim
          if(flag(1))then
            write(*,100)' D0ms3ir2_coli:  (3.57)',
     &          '    case with 2 soft singularities and 4 zero masses:',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            flag(1)=.false.
          endif
#endif 

c according to Ellis, Zanderighi, JHEP 0802 (2008) 002
c and Bern, Dixon, Kosower, Nucl. Phys. B412 (1994) 751
c Ellis, Ross, Terrano, NPB 178 (1981) 421      box2
c Denner Dittmaier (3.57)
          D0ms3ir2_coli =
#ifdef SING
     &        2d0*delta2ir
     &        -2d0*delta1ir*(cln_coli(-p13/muir2,-1d0)
     &        +cln_coli(-p24/muir2,-1d0)
     &        -cln_coli(-p34/muir2,-1d0))
#endif
     &        -(cln_coli(-p34/muir2,-1d0))**2
     &        +2d0*(
     &         cln_coli(-p13/muir2,-1d0)*cln_coli(-p24/muir2,-1d0)
C     &        -cspcon_coli(-p34,-1d0/p13,p34/p13,1d0-p34/p13,-1d0,1d0)
C     &        -cspcon_coli(-p34,-1d0/p24,p34/p24,1d0-p34/p24,-1d0,1d0)
     &        -cspenc_coli(1d0-p34/p13,dreal(p34-p13))
     &        -cspenc_coli(1d0-p34/p24,dreal(p34-p24))
     &        )
     &        -4d0*pi2_6

          D0ms3ir2_coli = D0ms3ir2_coli/(p13*p24)

        elseif(m32.eq.cd0.or.m42.eq.cd0)then
* 2 soft singularities and 3 zero masses
*
*                   0             
*          0 ----------------  m32  
*                |  2   |         
*              0 |1    3| m32     
*                |   4  |         
*          0 ----------------  p34
*                   0             
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.0d0.or.p12.ne.0d0.or.p34.eq.cd0
     &        .or.p14.ne.m42.or.p23.ne.m32.or.m32*m42.ne.cd0)
     &        then
            write(*,100)' D0ms3ir2_coli:',
     &          '    case with 2 soft singularities and 3 zero masses:',
     &          '    wrong arguments'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(2))then
            write(*,100)' D0ms3ir2_coli:',
     &          '    case with 2 soft singularities and 3 zero masses:',
     &          '    not yet tested at all'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            flag(2)=.false.
          endif
#endif
          
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002
c and Bern, Dixon, Kosower, Nucl. Phys. B412 (1994) 751
c Ellis, Ross, Terrano, NPB 178 (1981) 421              box7
c Denner Dittmaier (3.58)

          mm32=(m32+m42)*coliminfscale2

          D0ms3ir2_coli =
#ifdef SING
     &        1.5d0*delta2ir
     &        -delta1ir*(2d0*cln_coli(-p13/sqrt(muir2*mm32),-1d0)
     &        +cln_coli(-p24/muir2,-1d0)
     &        -cln_coli(-p34/muir2,-1d0)
     &        +.5d0*cln_coli(mm32/muir2,-1d0))
     &        +1.75*colishiftms2
#endif
     &        +2d0*cln_coli(-p13/sqrt(muir2*mm32),-1d0)
     &        *cln_coli(-p24/muir2,-1d0)
     &        -(cln_coli(-p34/muir2,-1d0))**2
     &        -cln_coli(mm32/muir2,-1d0)*(
     &         2d0*cln_coli(-p13/sqrt(muir2*mm32),-1d0)
     &        +cln_coli(-p24/muir2,-1d0)
     &        -2d0*cln_coli(-p34/muir2,-1d0)
     &        +.75d0*cln_coli(mm32/muir2,-1d0))
     &        -2d0*
     &        cspcon_coli(-p34,-1d0/p13,p34/p13,1d0-p34/p13,-1d0,1d0)
     &        -4d0*pi2_6
 
          D0ms3ir2_coli = D0ms3ir2_coli/(p13*p24)

        else  ! m32.ne.cd0.and.m42.ne.cd0
* 2 soft singularities and 2 zero masses
*
*                   0             
*          0 ----------------  m32  
*                |  2   |         
*              0 |1    3| m32     
*                |   4  |         
*        m42 ----------------  p34
*                  m42            
*

#ifdef CHECK
          if(m12.ne.cd0.or.m22.ne.0d0.or.p12.ne.0d0.or.p34.eq.cd0
     &        .or.p14.ne.m42.or.p23.ne.m32.or.m32*m42.eq.0d0)
     &        then
            write(*,100)' D0ms3ir2_coli:  (3.64)',
     &          '    case with 2 soft singularities',
     &          '    and 2 zero masses: wrong arguments'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(3))then
            write(*,100)' D0ms3ir2_coli:  (3.64)',
     &          '    case with 2 soft singularities',
     &          '    and 2 zero masses not yet testted'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            flag(3)=.false.
          endif
#endif
          
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002   box11 m32,m42 -> 0
c Denner Dittmaier (3.64)

          D0ms3ir2_coli =
#ifdef SING
     &        delta2ir
     &        -delta1ir*(cln_coli(-p13/sqrt(muir2*m32*coliminfscale2),
     &        -1d0)
     &        +cln_coli(-p24/sqrt(muir2*m42*coliminfscale2),-1d0))
     &        -colishiftms2
#endif
     &        +2d0*cln_coli(-p13/sqrt(muir2*m32*coliminfscale2),-1d0)
     &        *cln_coli(-p24/sqrt(muir2*m42*coliminfscale2),-1d0)
     &        -(cln_coli(-p34/(sqrt(m32*m42)*coliminfscale2),-1d0))**2
     &        -4d0*pi2_6
 
          D0ms3ir2_coli = D0ms3ir2_coli/(p13*p24)

        endif
        
      elseif(nsoft.eq.1)then

        if(m32.eq.cd0.and.m42.eq.cd0.and.(m12.eq.cd0.or.m22.eq.cd0))then
* 1 soft singularity and three zero masses
*
*                   0             
*         m2 ----------------  0  
*                |  2   |         
*             m2 |1    3| 0       
*                |   4  |         
*         m2 ----------------  p34  
*                   0             
*

#ifdef CHECK
          if(.not.(
     &        (m12.eq.p12.and.m12.eq.p14.and.p23.eq.cd0.and.p34.ne.cd0
     &        .and.m22.eq.cd0.and.m32.eq.cd0.and.m42.eq.cd0).or.
     &        (m22.eq.p12.and.m22.eq.p23.and.p14.eq.cd0.and.p34.ne.cd0
     &        .and.m12.eq.cd0.and.m42.eq.cd0.and.m32.eq.cd0)))
     &        then
            write(*,100)' D0ms3ir2_coli:  (3.63)',
     &          '    case with 1 soft singularity ',
     &          '    and three zero masses: wrong arguments'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(4))then
            write(*,100)' D0ms3ir2_coli:  (3.63)',
     &          '    case with 1 soft singularity ',
     &          '    and 3 zero masses not yet tested at all'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            flag(4)=.false.
          endif
#endif
          
          if(m12.ne.cd0)then 
            mm12=m12*coliminfscale2 
            q13=p13
            q24=p24
          else
            mm12=m22*coliminfscale2
            q13=p24
            q24=p13
          endif
          
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002   box9  p32=m2 -> 0
c Denner Dittmaier (3.63)

          D0ms3ir2_coli=
#ifdef SING
     &        .5d0*delta2ir
     &        -delta1ir*(cln_coli(-q13/sqrt(mm12*muir2),-1d0)
     &        +cln_coli(-q24/p34,-1d0)) 
     &        +.25d0*colishiftms2
#endif
     &        +2d0*
     &        cspcon_coli(-q24,-1d0/p34,q24/p34,1d0-q24/p34,-1d0,1d0)
     &        +(cln_coli(-q13/sqrt(mm12*muir2),-1d0)
     &        *cln_coli(-q24/p34,-1d0))**2
     &        +pi2_6

          D0ms3ir2_coli=D0ms3ir2_coli/(q13*q24)
          
        elseif((m12.eq.cd0.and.m42.eq.cd0).or.
     &         (m22.eq.cd0.and.m32.eq.cd0))then
* 1 soft singularity and two adjacent zero masses
*
*                   0             
*         m2 ----------------  0  
*                |  2   |         
*             m2 |1    3| 0       
*                |   4  |         
*          0 ----------------  p34  
*                   m2             
*

#ifdef CHECK
          if(.not.(
     &        (m12.eq.p12.and.m12.eq.m42.and.p23.eq.cd0.and.p34.ne.cd0
     &        .and.m22.eq.cd0.and.m32.eq.cd0.and.p14.eq.cd0).or.
     &        (m22.eq.p12.and.m22.eq.m32.and.p14.eq.cd0.and.p34.ne.cd0
     &        .and.m12.eq.cd0.and.m42.eq.cd0.and.p23.eq.cd0)))
     &        then
            write(*,100)' D0ms3ir2_coli:  (3.59)',
     &          '    case with 1 soft singularity ',
     &          '    and two adjacent zero masses: wrong arguments'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(5))then
            write(*,100)' D0ms3ir2_coli:  (3.59)',
     &          '    case with 1 soft singularity ',
     &          '    and 2 adjacent zero masses not yet tested at all'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            flag(5)=.false.
          endif
#endif
          
          if(m12.ne.cd0)then 
            mm12=m12*coliminfscale2 
            q13=p13
            q24=p24
          else
            mm12=m22*coliminfscale2
            q13=p24
            q24=p13
          endif
          
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002  
c              box12 m3=m4 small p32=0                   
c Denner Dittmaier (3.59)

          D0ms3ir2_coli=
#ifdef SING
     &        .5d0*delta2ir
     &        -delta1ir*(cln_coli(-q24/sqrt(mm12*muir2),-1d0)
     &        +cln_coli(q13/p34,dreal(p34-q13))) 
     &        +.25d0*colishiftms2
#endif
     &        -2d0*(
     &        cspcon_coli(-p34,-1d0/q24,p34/q24,1d0-p34/q24,-1d0,1d0)
     &        +cspcon_coli(-p34,-1d0/q13,p34/q13,1d0-p34/q13,-1d0,1d0)
     &        )
     &        +2d0*cln_coli(-q24/sqrt(mm12*muir2),-1d0)
     &        *cln_coli(-q13/sqrt(mm12*muir2),-1d0)
     &        -(cln_coli(-p34/sqrt(mm12*muir2),-1d0))**2
     &        -pi2_6
c constant differs from std: -1/2 pi2

          D0ms3ir2_coli=D0ms3ir2_coli/(q13*q24)
          
        elseif((m12.eq.cd0.and.m32.eq.cd0).or.
     &         (m22.eq.cd0.and.m42.eq.cd0))then
* 1 soft singularity and two opposite zero masses
*
*                   m22            
*        m22 ----------------  m22  
*                |  2   |         
*              0 |1    3| 0       
*                |   4  |         
*        m42 ----------------  p34  
*                   m42             
*

#ifdef CHECK
          if(.not.(
     &        (m22.eq.p12.and.m22.eq.p23.and.p14.eq.m42.and.p34.ne.cd0
     &        .and.m12.eq.cd0.and.m32.eq.cd0).or.
     &        (m12.eq.p12.and.m12.eq.p14.and.p23.eq.m32.and.p34.ne.cd0
     &        .and.m22.eq.cd0.and.m42.eq.cd0)))
     &        then
            write(*,100)' D0ms3ir2_coli:  (3.66)',
     &          '    case with 1 soft singularity ',
     &          '    and two opposite zero masses: wrong arguments'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(6))then
            write(*,100)' D0ms3ir2_coli:  (3.66)',
     &          '    case with 1 soft singularity ',
     &          '    and 2 opposite zero masses not yet tested at all'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            flag(6)=.false.
          endif
#endif

          if(m22.ne.cd0)then 
            mm22=m22*coliminfscale2 
            mm42=m42*coliminfscale2 
            q13=p13
            q24=p24
          else
            mm22=m12*coliminfscale2
            mm42=m32*coliminfscale2 
            q13=p24
            q24=p13
          endif
          
c according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 iid
c Ellis, Zanderighi, JHEP 0802 (2008) 002  
c              box15 m22,m42 small,  p22=m22                   
c Denner Dittmaier (3.66)
          D0ms3ir2_coli= 
     &        cln_coli(-q24/sqrt(mm22*mm42),-1d0)*
     &        (cln_coli(-q24/sqrt(mm22*mm42),-1d0)
     &        +cln_coli(mm42/muir2,0d0)
#ifdef SING
     &        -delta1ir
     &        + colishiftms2
#endif
     &        +2d0*cln_coli(q13/p34,dreal(p34-q13)))
     &        + pi2_6 

          D0ms3ir2_coli=D0ms3ir2_coli/(q13*q24)
          
        elseif(m12.eq.cd0.or.m22.eq.cd0)then
* 1 soft singularity and one zero mass
*
*                   m22            
*        m22 ----------------  0  
*                |  2   |         
*              0 |1    3| m22       
*                |   4  |         
*        m42 ----------------  p34  
*                   m42             
*

#ifdef CHECK
          if(.not.(
     &        (m22.eq.p12.and.m22.eq.m32.and.p14.eq.m42.and.p34.ne.cd0
     &        .and.m12.eq.cd0.and.p23.eq.cd0).or.
     &        (m12.eq.p12.and.m12.eq.m42.and.p23.eq.m32.and.p34.ne.cd0
     &        .and.m22.eq.cd0.and.p14.eq.cd0)))
     &        then
            write(*,100)' D0ms3ir2_coli:  (3.65)',
     &          '    case with 1 soft singularity and one zero mass: ',
     &          '    wrong arguments'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef ALLCHECK
c ppenj
          if(flag(7))then
            write(*,100)' D0ms3ir2_coli:  (3.65)',
     &          '    case with 1 soft singularity and 1 zero mass',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            flag(7)=.false.
          endif
#endif          

          if(m22.ne.cd0)then 
            mm22=m22*coliminfscale2 
            mm42=m42*coliminfscale2 
            q13=p13
            q24=p24
          else
            mm22=m12*coliminfscale2
            mm42=m32*coliminfscale2 
            q13=p24
            q24=p13
          endif

c  according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349 (id)
c    cD0ireem
c Ellis, Zanderighi, JHEP 0802 (2008) 002  
c              box16 m22,m42 small,  p22=0, m32=m22                   
c Denner Dittmaier (3.65)
          D0ms3ir2_coli= 2d0*cln_coli(-q24/sqrt(mm22*mm42),-1d0)
c     &     *(cln_coli(-q13/sqrt(p34*muir2),-1d0)
     &     *(cln_coli(-q13/sqrt(mm22*muir2),-1d0)
#ifdef SING
     &        -.5d0*delta1ir
#endif
     &        )
     &        - (cln_coli(-p34/sqrt(mm22*mm42),-1d0))**2 - pi2_6
     &        - 2d0*
     &        cspenc_coli(1d0-p34/q24,dreal(p34-q24))

          D0ms3ir2_coli=D0ms3ir2_coli/(q13*q24)
          
#ifdef CHECK
        else
          write(*,100)' D0ms3ir2: inconsistency 1 soft'
#endif
        endif

      elseif(nsoft.eq.0)then

        if(m32.eq.cd0.and.m42.eq.cd0)then
* no soft singularity and two zero masses
*
*                   m2             
*          0 ----------------  m2  
*                |  2   |         
*             m2 |1    3| 0       
*                |   4  |         
*         m2 ----------------  p34  
*                   0             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.ne.m12.or.
     &        p34.eq.cd0.or.p23.ne.m12.or.p14.ne.m12.or.p12.ne.cd0)
     &        then
            write(*,111)' D0ms3ir2_coli:  (3.62)',
     &          '    case with 0 soft singularity ',
     &          '    and 2 zero masses: wrong arguments'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
c eeqqg
          if(flag(8))then
            write(*,100)' D0ms3ir2_coli:  (3.62)',
     &          '    case with 0 soft singularity ',
     &          '    and 2 zero masses not yet tested at all'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            flag(8)=.false.
          endif
#endif
          
c cD0ms00eee
c Denner Dittmaier (3.62)

#ifdef SING
#endif


          mm12=m12*coliminfscale2
c     m42=m42i*minfscale2

          D0ms3ir2_coli= -2d0*pi2_6 
     &        -(cln_coli(-p24,-1d0)+cln_coli(-p13,-1d0)
     &        -cln_coli(-p34,-1d0)-cln_coli(mm12,-1d0))**2
#ifdef SING
     &        -colishiftms2
#endif

          D0ms3ir2_coli=-D0ms3ir2_coli/(p13*p24)

c     write(*,*) 'D0ms3ir2_coli = ',D0ms3ir2_coli

#ifdef CHECK
          if(conscheck) then

            l23=-p34
            l13=-p13
            l24=-p24
            x=mm12*l23/(l13*l24)
            ix=-1d0

            D0ms3ir2_check= 4d0*pi2_6 + cln_coli(-x,-ix)**2
     &          +2d0*cln_coli(-x,-ix)
     &          *(cln_coli(l24,-1d0)+cln_coli(l13,-1d0)
     &          -cln_coli(l23,-1d0)-cln_coli(mm12,-1d0))
#ifdef SING
     &          -colishiftms2
#endif
            
            D0ms3ir2_check=-D0ms3ir2_check/(l13*l24)

c     write(*,*) 'D0ms3ir2_check= ',D0ms3ir2_check   

            if(cdabs(D0ms3ir2_coli/D0ms3ir2_check-1d0).gt.testacc) then
              write(testout,100) ' D0ms3ir2_coli inconsistent'
              write(testout,111) ' D0ms3ir2_coli: p12 = ',p12
              write(testout,111) ' D0ms3ir2_coli: p23 = ',p23
              write(testout,111) ' D0ms3ir2_coli: p34 = ',p34
              write(testout,111) ' D0ms3ir2_coli: p14 = ',p14
              write(testout,111) ' D0ms3ir2_coli: p24 = ',p24
              write(testout,111) ' D0ms3ir2_coli: p13 = ',p13
              write(testout,111) ' D0ms3ir2_coli: m12 = ',m12
              write(testout,111) ' D0ms3ir2_coli: m22 = ',m22
              write(testout,111) ' D0ms3ir2_coli: m32 = ',m32
              write(testout,111) ' D0ms3ir2_coli: m42 = ',m42
              write(testout,111) ' D0ms3ir2_coli  = ',D0ms3ir2_coli
              write(testout,111) ' D0ms3ir2_check = ',D0ms3ir2_check
              write(testout,101) ' diff        = ',
     &            cdabs(D0ms3ir2_coli/D0ms3ir2_check-1d0)
              
              write(*,100) ' D0ms3ir2_coli inconsistent'
              write(*,111) ' D0ms3ir2_coli: p12 = ',p12
              write(*,111) ' D0ms3ir2_coli: p23 = ',p23
              write(*,111) ' D0ms3ir2_coli: p34 = ',p34
              write(*,111) ' D0ms3ir2_coli: p14 = ',p14
              write(*,111) ' D0ms3ir2_coli: p24 = ',p24
              write(*,111) ' D0ms3ir2_coli: p13 = ',p13
              write(*,111) ' D0ms3ir2_coli: m12 = ',m12
              write(*,111) ' D0ms3ir2_coli: m22 = ',m22
              write(*,111) ' D0ms3ir2_coli: m32 = ',m32
              write(*,111) ' D0ms3ir2_coli: m42 = ',m42
              write(*,111) ' D0ms3ir2_coli  = ',D0ms3ir2_coli
              write(*,111) ' D0ms3ir2_check = ',D0ms3ir2_check
              write(*,101) ' diff        = ',
     &            cdabs(D0ms3ir2_coli/D0ms3ir2_check-1d0)
              stop
            endif
          endif
#endif

          
        elseif(m32.eq.cd0.or.m42.eq.cd0)then
* 0 soft singularity and 1 zero mass
*
*                   m2             
*          0 ----------------  m2 
*                |  2   |         
*             m2 |1    3| 0       
*                |   4  |         
*          0 ----------------  p34  
*                   m2             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.ne.m12.or.
     &        p12.ne.cd0.or.p34.eq.cd0.or.p34.eq.m12
     &        .or..not.(p23.eq.m22.and.m32.eq.cd0
     &                  .and.m12.eq.m42.and.p14.eq.cd0
     &              .or.p14.eq.m12.and.m42.eq.cd0
     &                  .and.m22.eq.m32.and.p23.eq.cd0))
     &        then
            write(*,100)' D0ms3ir2_coli:   (3.60)',
     &          '    case with 0 soft singularity ',
     &          '    and 1 zero mass: wrong arguments'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(9))then
            write(*,100)' D0ms3ir2_coli:   (3.60)',
     &          '    case with 0 soft singularity ',
     &          '    and 1 zero mass not yet implemented'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            flag(9)=.false.
            stop
          endif
#endif
          
          if(m42.ne.cd0)then 
            mm42=m42*coliminfscale2 
          else
            mm42=m32*coliminfscale2
          endif
          
c new case
c Denner Dittmaier (3.60)

          D0ms3ir2_coli=0d0

          D0ms3ir2_coli=D0ms3ir2_coli/(p13*p24)
          
        elseif((m12.eq.cd0.and.m32.eq.cd0).or.
     &         (m22.eq.cd0.and.m42.eq.cd0))then
* 0 soft singularity and 0 zero masses
*
*                   m2            
*         0 ----------------  0  
*                |  2   |         
*             m2 |1    3| m2       
*                |   4  |         
*         0 ----------------  p34  
*                   m2             
*

#ifdef CHECK
          if(m12.eq.cd0.or.m22.ne.m12.or.m32.ne.m12.or.m42.ne.m12.or.
     &        p12.ne.cd0.or.p23.ne.cd0.or.p14.ne.cd0.or.p34.eq.cd0)
     &        then
            write(*,100)' D0ms3ir2_coli:   (3.61)',
     &          '    case with 0 soft singularity ',
     &          '    and 0 zero masses: wrong arguments'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            stop
          endif
#endif
          
#ifdef CHECK
          if(flag(10))then
            write(*,100)' D0ms3ir2_coli:   (3.61)',
     &          '    case with 0 soft singularity ',
     &          '    and 0 zero masses not implemented'
            write(*,111)' D0ms3ir2_coli: p12 = ',p12
            write(*,111)' D0ms3ir2_coli: p23 = ',p23
            write(*,111)' D0ms3ir2_coli: p34 = ',p34
            write(*,111)' D0ms3ir2_coli: p14 = ',p14
            write(*,111)' D0ms3ir2_coli: p24 = ',p24
            write(*,111)' D0ms3ir2_coli: p13 = ',p13
            write(*,111)' D0ms3ir2_coli: m12 = ',m12
            write(*,111)' D0ms3ir2_coli: m22 = ',m22
            write(*,111)' D0ms3ir2_coli: m32 = ',m32
            write(*,111)' D0ms3ir2_coli: m42 = ',m42
            flag(10)=.false.
          endif
#endif
          
c new case
c Denner Dittmaier (3.61)
          D0ms3ir2_coli= cd0

          D0ms3ir2_coli=D0ms3ir2_coli/(p13*p24)
          

#ifdef CHECK
        else
          write(*,100)' D0ms3ir2: inconsistency 1ir'
#endif
        endif

#ifdef CHECK
      else
        write(*,111)' D0ms3ir2_coli: case with more than 2 soft 
     &      singularities not implemented '
        write(*,111)' D0ms3ir2_coli: nsoft = ',nsoft
        write(*,111)' D0ms3ir2_coli: p23 = ',p23
        write(*,111)' D0ms3ir2_coli: p34 = ',p34
        write(*,111)' D0ms3ir2_coli: p14 = ',p14
        write(*,111)' D0ms3ir2_coli: p24 = ',p24
        write(*,111)' D0ms3ir2_coli: p13 = ',p13
        write(*,111)' D0ms3ir2_coli: m12 = ',m12
        write(*,111)' D0ms3ir2_coli: m22 = ',m22
        write(*,111)' D0ms3ir2_coli: m32 = ',m32
        write(*,111)' D0ms3ir2_coli: m42 = ',m42    
        stop
#endif
      endif

      end


************************************************************************
      function D0ms4ir4_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
************************************************************************
*     quartic mass-singular 4-point function                           *
*     p12, p23, p34, p14, m12, m22, m32, m42 small                     *
*     p13, p24  finite                                                 *
*                                                                      *
*                     m22 small                                        *
*     p12 small  ---------------------  p23 small                      *
*                     |    2    |                                      *
*                     |         |                                      *
*           m12 small | 1     3 | m32 small                            *
*                     |         |                                      *
*                     |    4    |                                      *
*     p14 small  ---------------------  p34 small                      *
*                     m42 small                                        *
*                                                                      *
************************************************************************
*     05.08.08 Ansgar Denner        last changed  27.07.09             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 D0ms4ir4_coli

#ifdef CHECK
      complex*16 elimminf2_coli
      complex*16 ps12,ps23,ps34,ps14,ps13,ps24
      complex*16 ms12,ms22,ms32,ms42
      logical flag(0:6)
      data    flag /7*.true./
      save    flag  
#endif

      complex*16 q13,q24,mm12,mm22,mm42,mm32,m2(4)
      complex*16 cln_coli
      integer    nsoft
      integer    i,j,k
      logical    soft(4),onsh(4,4),softt(4,4,4)

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
        ms42 = elimminf2_coli(m42)
        ps12 = elimminf2_coli(p12)
        ps23 = elimminf2_coli(p23)
        ps34 = elimminf2_coli(p34)
        ps14 = elimminf2_coli(p14)
        ps24 = elimminf2_coli(p24)
        ps13 = elimminf2_coli(p13)

        if(ms12.ne.cd0.or.ms22.ne.cd0.or.ms32.ne.cd0
     &        .or.ms42.ne.cd0.or.ps12.ne.cd0.or.ps23.ne.cd0 
     &        .or.ps34.ne.cd0.or.ps14.ne.cd0
     &        .or.ps24.eq.cd0.or.ps13.eq.cd0
     &        .or.dimag(p12).ne.0d0.or.dimag(p23).ne.0d0
     &        .or.dimag(p34).ne.0d0.or.dimag(p14).ne.0d0
     &        .or.dimag(p24).ne.0d0.or.dimag(p13).ne.0d0
     &        .or.ps13.eq.cd0.or.ps24.eq.cd0) then
          write(*,100)' D0ms4ir4_coli called improperly:'
          write(*,111)' D0ms4ir4_coli: p12 = ',p12,ps12
          write(*,111)' D0ms4ir4_coli: p23 = ',p23,ps23
          write(*,111)' D0ms4ir4_coli: p34 = ',p34,ps34
          write(*,111)' D0ms4ir4_coli: p14 = ',p14,ps14
          write(*,111)' D0ms4ir4_coli: p13 = ',p13,ps13
          write(*,111)' D0ms4ir4_coli: p24 = ',p24,ps24
          write(*,111)' D0ms4ir4_coli: m12 = ',m12,ms12
          write(*,111)' D0ms4ir4_coli: m22 = ',m22,ms22
          write(*,111)' D0ms4ir4_coli: m32 = ',m32,ms32
          write(*,111)' D0ms4ir4_coli: m42 = ',m42,ms42
          stop
        end if
      end if
#endif

      m2(1)=m12
      m2(2)=m22
      m2(3)=m32
      m2(4)=m42

c determine on-shell momenta
      onsh(1,2)=p12.eq.m22
      onsh(1,3)=p13.eq.m32
      onsh(1,4)=p14.eq.m42
      onsh(2,1)=p12.eq.m12
      onsh(2,3)=p23.eq.m32
      onsh(2,4)=p24.eq.m42
      onsh(3,1)=p13.eq.m12
      onsh(3,2)=p23.eq.m22
      onsh(3,4)=p34.eq.m42
      onsh(4,1)=p14.eq.m12
      onsh(4,2)=p24.eq.m22
      onsh(4,3)=p34.eq.m32

c count/determine soft singularities
      nsoft=0 
      do i=1,4
        soft(i)=.false.
        do j=1,3
          if(i.ne.j)then
            do k=j,4
              if(k.ne.i.and.k.ne.j)then
                softt(i,j,k)=m2(i).eq.cd0.and.onsh(i,j).and.onsh(i,k)
                soft(i)=softt(i,j,k).or.soft(i)
                if(softt(i,j,k)) nsoft=nsoft+1
              endif
            enddo
          endif
        enddo
      enddo

      if(nsoft.eq.4)then

* four soft singularities
*
*                   0             
*          0 ----------------  0  
*                |  2   |         
*              0 |1    3| 0       
*                |   4  |         
*          0 ----------------  0  
*                   0             
*

#ifdef CHECK
        if(m12.ne.cd0.or.m22.ne.0d0.or.m32.ne.0d0.or.m42.ne.cd0
     &      .or.p12.ne.cd0.or.p23.ne.cd0.or.p34.ne.cd0.or.p14.ne.cd0)
     &      then
          write(*,100)' D0ms4ir4_coli:   (3.68)',
     &          '    case with 4 soft singularities: ',
     &        '    wrong arguments'
          write(*,111)' D0ms4ir4_coli: p12 = ',p12
          write(*,111)' D0ms4ir4_coli: p23 = ',p23
          write(*,111)' D0ms4ir4_coli: p34 = ',p34
          write(*,111)' D0ms4ir4_coli: p14 = ',p14
          write(*,111)' D0ms4ir4_coli: p24 = ',p24
          write(*,111)' D0ms4ir4_coli: p13 = ',p13
          write(*,111)' D0ms4ir4_coli: m12 = ',m12
          write(*,111)' D0ms4ir4_coli: m22 = ',m22
          write(*,111)' D0ms4ir4_coli: m32 = ',m32
          write(*,111)' D0ms4ir4_coli: m42 = ',m42
          stop
        endif
#endif

#ifdef CHECK
        if(flag(1))then
          write(*,100)' D0ms4ir4_coli:   (3.68)',
     &        '    case with 4 soft singularities',
     &        '    not yet tested in physical process'
          write(*,111)' D0ms4ir4_coli: p12 = ',p12
          write(*,111)' D0ms4ir4_coli: p23 = ',p23
          write(*,111)' D0ms4ir4_coli: p34 = ',p34
          write(*,111)' D0ms4ir4_coli: p14 = ',p14
          write(*,111)' D0ms4ir4_coli: p24 = ',p24
          write(*,111)' D0ms4ir4_coli: p13 = ',p13
          write(*,111)' D0ms4ir4_coli: m12 = ',m12
          write(*,111)' D0ms4ir4_coli: m22 = ',m22
          write(*,111)' D0ms4ir4_coli: m32 = ',m32
          write(*,111)' D0ms4ir4_coli: m42 = ',m42
          flag(1)=.false.
        endif
#endif

c Bern, Dixon, Kosower, Nucl. Phys. B412 (1994) 751
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002    box1
c Denner, Dittmaier (3.68)
        D0ms4ir4_coli = 
#ifdef SING
     &      4d0*delta2ir
     &      -2d0*delta1ir*(cln_coli(-p13/muir2,-1d0)
     &                     +cln_coli(-p24/muir2,-1d0))
#endif
     &      +2d0*cln_coli(-p13/muir2,-1d0)
     &                     *cln_coli(-p24/muir2,-1d0)
     &      -10d0*pi2_6

        D0ms4ir4_coli = D0ms4ir4_coli/(p13*p24)

      elseif(nsoft.eq.3)then

* three soft singularities
*
*                   0             
*         m2 ----------------  0  
*                |  2   |         
*             m2 |1    3| 0       
*                |   4  |         
*         m2 ----------------  0  
*                   0             
*

#ifdef CHECK
        if(.not.(
     &      (m12.eq.p12.and.m12.eq.p14.and.p23.eq.cd0.and.p34.eq.cd0
     &      .and.m22.eq.cd0.and.m32.eq.cd0.and.m42.eq.cd0).or.
     &      (m22.eq.p12.and.m22.eq.p23.and.p14.eq.cd0.and.p34.eq.cd0
     &      .and.m12.eq.cd0.and.m32.eq.cd0.and.m42.eq.cd0).or.
     &      (m32.eq.p23.and.m32.eq.p34.and.p14.eq.cd0.and.p12.eq.cd0
     &      .and.m12.eq.cd0.and.m22.eq.cd0.and.m42.eq.cd0).or.
     &      (m42.eq.p34.and.m42.eq.p14.and.p23.eq.cd0.and.p12.eq.cd0
     &      .and.m12.eq.cd0.and.m22.eq.cd0.and.m32.eq.cd0)))
     &      then
          write(*,100)' D0ms4ir4_coli:   (3.69)',
     &        '    case with 3 soft singularities: ',
     &        '    wrong arguments'
          write(*,111)' D0ms4ir4_coli: p12 = ',p12
          write(*,111)' D0ms4ir4_coli: p23 = ',p23
          write(*,111)' D0ms4ir4_coli: p34 = ',p34
          write(*,111)' D0ms4ir4_coli: p14 = ',p14
          write(*,111)' D0ms4ir4_coli: p24 = ',p24
          write(*,111)' D0ms4ir4_coli: p13 = ',p13
          write(*,111)' D0ms4ir4_coli: m12 = ',m12
          write(*,111)' D0ms4ir4_coli: m22 = ',m22
          write(*,111)' D0ms4ir4_coli: m32 = ',m32
          write(*,111)' D0ms4ir4_coli: m42 = ',m42
          stop
        endif
#endif

#ifdef CHECK
        if(flag(2))then
          write(*,100)' D0ms4ir4_coli:   (3.69)',
     &        '    case with 3 soft singularities',
     &        '    not yet tested in physical process'
          write(*,111)' D0ms4ir4_coli: p12 = ',p12
          write(*,111)' D0ms4ir4_coli: p23 = ',p23
          write(*,111)' D0ms4ir4_coli: p34 = ',p34
          write(*,111)' D0ms4ir4_coli: p14 = ',p14
          write(*,111)' D0ms4ir4_coli: p24 = ',p24
          write(*,111)' D0ms4ir4_coli: p13 = ',p13
          write(*,111)' D0ms4ir4_coli: m12 = ',m12
          write(*,111)' D0ms4ir4_coli: m22 = ',m22
          write(*,111)' D0ms4ir4_coli: m32 = ',m32
          write(*,111)' D0ms4ir4_coli: m42 = ',m42
          flag(2)=.false.
        endif
#endif

        if(m12.ne.cd0)then 
          mm12=m12*coliminfscale2
          q13=p13
          q24=p24
        elseif(m22.ne.cd0)then
          mm12=m22*coliminfscale2
          q13=p24
          q24=p13
        elseif(m32.ne.cd0)then
          mm12=m32*coliminfscale2
          q13=p13
          q24=p24
        elseif(m42.ne.cd0)then
          mm12=m42*coliminfscale2
          q13=p24
          q24=p13
        endif
         
c Beenakker, Kuiff, van Neerven, Smith, PRD 40 (1989) 54 
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002    box6  m4 small 
c Denner, Dittmaier (3.69)
        D0ms4ir4_coli=
#ifdef SING
     &      2d0*delta2ir
     &      -delta1ir*(2d0*cln_coli(-q13/sqrt(mm12*muir2),-1d0)
     &                 +cln_coli(-q24/muir2,-1d0)) 
#endif
     &      +2d0*cln_coli(-q13/sqrt(mm12*muir2),-1d0)
     &      *cln_coli(-q24/muir2,-1d0)
     &      -5d0*pi2_6

        D0ms4ir4_coli=D0ms4ir4_coli/(q13*q24)


      elseif(nsoft.eq.2)then

        if(soft(1).and.soft(3).or.soft(2).and.soft(4))then

* two opposite soft singularities
*
*                  m22            
*        m22 ---------------- m22
*                |  2   |         
*              0 |1    3| 0       
*                |   4  |         
*        m42 ---------------- m42
*                  m42            
*

#ifdef CHECK
          if(.not.(
     &        (p12.eq.m22.and.m22.eq.p23.and.m12.eq.cd0.and.m32.eq.cd0
     &        .and.p14.eq.m42.and.m42.eq.p34).or.
     &        (p23.eq.m32.and.m32.eq.p34.and.m22.eq.cd0.and.m42.eq.cd0
     &        .and.p12.eq.m12.and.m12.eq.p14).or.
     &        (p12.eq.m22.and.m22.eq.p24.and.m12.eq.cd0.and.m42.eq.cd0
     &        .and.p13.eq.m32.and.m32.eq.p34).or.
     &        (p13.eq.m32.and.m32.eq.p23.and.m12.eq.cd0.and.m22.eq.cd0
     &        .and.p14.eq.m42.and.m42.eq.p24).or.
     &        (p24.eq.m42.and.m42.eq.p34.and.m22.eq.cd0.and.m32.eq.cd0
     &        .and.p12.eq.m12.and.m12.eq.p13).or.
     &        (p23.eq.m22.and.m22.eq.p24.and.m32.eq.cd0.and.m42.eq.cd0
     &        .and.p13.eq.m12.and.m12.eq.p14)))
     &        then
            write(*,100)' D0ms4ir4_coli:   (3.73)',
     &          '    case with 2 opposite soft  singularities:',
     &          '    wrong arguments'
            write(*,111)' D0ms4ir4_coli: p12 = ',p12
            write(*,111)' D0ms4ir4_coli: p23 = ',p23
            write(*,111)' D0ms4ir4_coli: p34 = ',p34
            write(*,111)' D0ms4ir4_coli: p14 = ',p14
            write(*,111)' D0ms4ir4_coli: p24 = ',p24
            write(*,111)' D0ms4ir4_coli: p13 = ',p13
            write(*,111)' D0ms4ir4_coli: m12 = ',m12
            write(*,111)' D0ms4ir4_coli: m22 = ',m22
            write(*,111)' D0ms4ir4_coli: m32 = ',m32
            write(*,111)' D0ms4ir4_coli: m42 = ',m42
            stop
          endif
#endif

          if(soft(1).and.soft(3))then
            mm22=m22
            mm42=m42
            q13=p13
            q24=p24
          else
            mm22=m32
            mm42=m12
            q13=p24
            q24=p13
          endif
 
#ifdef CHECK
          if(flag(3))then
            write(*,100)' D0ms4ir4_coli:   (3.73)',
     &          '    case with 2 opposite soft singularities ',
     &          '    not yet tested in physcial process'
            write(*,111)' D0ms4ir4_coli: p12 = ',p12
            write(*,111)' D0ms4ir4_coli: p23 = ',p23
            write(*,111)' D0ms4ir4_coli: p34 = ',p34
            write(*,111)' D0ms4ir4_coli: p14 = ',p14
            write(*,111)' D0ms4ir4_coli: p24 = ',p24
            write(*,111)' D0ms4ir4_coli: p13 = ',p13
            write(*,111)' D0ms4ir4_coli: m12 = ',m12
            write(*,111)' D0ms4ir4_coli: m22 = ',m22
            write(*,111)' D0ms4ir4_coli: m32 = ',m32
            write(*,111)' D0ms4ir4_coli: m42 = ',m42
            flag(3)=.false.
          endif
#endif
 

* according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349  (iiib)
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002    box14 m2,m4 small
c Denner, Dittmaier (3.73)
          D0ms4ir4_coli =
#ifdef SING
     &        -2d0*cln_coli(-q24/(sqrt(mm22*mm42)*coliminfscale2),-1d0)
     &        *delta1ir
#endif
     &        +2d0*cln_coli(-q24/(sqrt(mm22*mm42)*coliminfscale2),-1d0)
     &        *cln_coli(-q13/muir2,-1d0)
          D0ms4ir4_coli = D0ms4ir4_coli/(q13*q24)

        else

* two adjacent soft singularities
*
*                   0             
*         0  ---------------- mm2
*                |  2   |         
*              0 |1    3| mm2     
*                |   4  |         
*        mm2 ----------------  0 
*                  mm2            
*

#ifdef CHECK
          if(.not.(
     &        (m12.eq.cd0.and.m22.eq.cd0.and.p12.eq.cd0.and.p34.eq.cd0
     &        .and.p14.eq.m42.and.m42.eq.m32.and.m32.eq.p23).or.
     &        (m22.eq.cd0.and.m32.eq.cd0.and.p23.eq.cd0.and.p14.eq.cd0
     &        .and.p12.eq.m12.and.m12.eq.m42.and.m42.eq.p34).or.
     &        (m32.eq.cd0.and.m42.eq.cd0.and.p34.eq.cd0.and.p12.eq.cd0
     &        .and.p23.eq.m22.and.m22.eq.m12.and.m12.eq.p14).or.
     &        (m42.eq.cd0.and.m12.eq.cd0.and.p14.eq.cd0.and.p23.eq.cd0
     &        .and.p34.eq.m32.and.m32.eq.m22.and.m22.eq.p12).or.
     &        (m32.eq.cd0.and.m12.eq.cd0.and.p13.eq.cd0.and.p24.eq.cd0
     &        .and.p34.eq.m42.and.m42.eq.m22.and.m22.eq.p12).or.
     &        (m42.eq.cd0.and.m22.eq.cd0.and.p24.eq.cd0.and.p13.eq.cd0
     &        .and.p34.eq.m32.and.m32.eq.m12.and.m12.eq.p12)))
     &        then
            write(*,100)' D0ms4ir4_coli:   (3.70)',
     &          '    case with 2 adjacent soft singularities:',
     &          '    wrong arguments'
            write(*,111)' D0ms4ir4_coli: p12 = ',p12
            write(*,111)' D0ms4ir4_coli: p23 = ',p23
            write(*,111)' D0ms4ir4_coli: p34 = ',p34
            write(*,111)' D0ms4ir4_coli: p14 = ',p14
            write(*,111)' D0ms4ir4_coli: p24 = ',p24
            write(*,111)' D0ms4ir4_coli: p13 = ',p13
            write(*,111)' D0ms4ir4_coli: m12 = ',m12
            write(*,111)' D0ms4ir4_coli: m22 = ',m22
            write(*,111)' D0ms4ir4_coli: m32 = ',m32
            write(*,111)' D0ms4ir4_coli: m42 = ',m42
            stop
          endif
#endif

#ifdef CHECK
          if(flag(4))then
            write(*,100)' D0ms4ir4_coli:   (3.70)',
     &          '    case with 2 adjacent soft singularities',
     &          '    not yet tested in physical process'
            write(*,111)' D0ms4ir4_coli: p12 = ',p12
            write(*,111)' D0ms4ir4_coli: p23 = ',p23
            write(*,111)' D0ms4ir4_coli: p34 = ',p34
            write(*,111)' D0ms4ir4_coli: p14 = ',p14
            write(*,111)' D0ms4ir4_coli: p24 = ',p24
            write(*,111)' D0ms4ir4_coli: p13 = ',p13
            write(*,111)' D0ms4ir4_coli: m12 = ',m12
            write(*,111)' D0ms4ir4_coli: m22 = ',m22
            write(*,111)' D0ms4ir4_coli: m32 = ',m32
            write(*,111)' D0ms4ir4_coli: m42 = ',m42
            flag(4)=.false.
          endif
#endif 
          if(soft(1).and.soft(2))then
            q13=p13
            q24=p24
            mm32=m32*coliminfscale2
          elseif(soft(2).and.soft(3))then
            q13=p24
            q24=p13
            mm32=m42*coliminfscale2
          elseif(soft(3).and.soft(4))then
            q13=p13
            q24=p24
            mm32=m12*coliminfscale2
          elseif(soft(4).and.soft(1))then
            q13=p24
            q24=p13
            mm32=m22*coliminfscale2
          endif

c Hoepker, DESY-T-96-02
c according to Ellis, Zanderighi, JHEP 0802 (2008) 002
c     box11 p32=0 m3=m4
c Denner, Dittmaier (3.70)

          D0ms4ir4_coli=
#ifdef SING
     &        delta2ir
     &        -delta1ir*(cln_coli(-q13/sqrt(mm32*muir2),-1d0)
     &                 +cln_coli(-q24/sqrt(mm32*muir2),-1d0))
     &        +0.5d0*colishiftms2
#endif
     &        +2d0*cln_coli(-q13/sqrt(mm32*muir2),-1d0)
     &        *cln_coli(-q24/sqrt(mm32*muir2),-1d0)
     &        -4d0*pi2_6
          
          D0ms4ir4_coli=D0ms4ir4_coli/(q13*q24)
          
        endif

      elseif(nsoft.eq.1)then

* 1 soft singularity
*
*                   m2            
*         m2 ----------------  0  
*                |  2   |         
*              0 |1    3| m2      
*                |   4  |         
*         m2 ----------------  0  
*                   m2            
*

#ifdef CHECK
        if(.not.(
     &      (m12.eq.cd0.and.p12.eq.m22.and.m22.eq.m32.and.m32.eq.m42
     &      .and.m42.eq.p14.and.p23.eq.cd0.and.p34.eq.cd0).or.
     &      (m22.eq.cd0.and.p23.eq.m32.and.m32.eq.m42.and.m42.eq.m12
     &      .and.m12.eq.p12.and.p34.eq.cd0.and.p14.eq.cd0).or.
     &      (m32.eq.cd0.and.p34.eq.m42.and.m42.eq.m12.and.m12.eq.m22
     &      .and.m22.eq.p23.and.p14.eq.cd0.and.p12.eq.cd0).or.
     &      (m42.eq.cd0.and.p14.eq.m12.and.m12.eq.m22.and.m22.eq.m32
     &      .and.m32.eq.p34.and.p12.eq.cd0.and.p23.eq.cd0)))
     &      then
          write(*,100)' D0ms4ir4_coli:   (3.71)',
     &        '    case with 1 soft singularity: ',
     &        '    wrong arguments'
          write(*,111)' D0ms4ir4_coli: p12 = ',p12
          write(*,111)' D0ms4ir4_coli: p23 = ',p23
          write(*,111)' D0ms4ir4_coli: p34 = ',p34
          write(*,111)' D0ms4ir4_coli: p14 = ',p14
          write(*,111)' D0ms4ir4_coli: p24 = ',p24
          write(*,111)' D0ms4ir4_coli: p13 = ',p13
          write(*,111)' D0ms4ir4_coli: m12 = ',m12
          write(*,111)' D0ms4ir4_coli: m22 = ',m22
          write(*,111)' D0ms4ir4_coli: m32 = ',m32
          write(*,111)' D0ms4ir4_coli: m42 = ',m42
          stop
        endif
#endif

        if(m12.eq.cd0)then 
          mm32=m32*coliminfscale2
          q13=p13
          q24=p24
        elseif(m22.eq.cd0)then
          mm32=m42*coliminfscale2
          q13=p24
          q24=p13
        elseif(m32.eq.cd0)then
          mm32=m12*coliminfscale2
          q13=p13
          q24=p24
        elseif(m42.eq.cd0)then
          mm32=m22*coliminfscale2
          q13=p24
          q24=p13
        endif

#ifdef CHECK
        if(flag(5)) then         
          write(*,100)' D0ms4ir4_coli:   (3.71)',
     &        '    case with 1 soft singularity ',
     &        '    not yet tested in physical process'
          write(*,111)' D0ms4ir4_coli: p12 = ',p12
          write(*,111)' D0ms4ir4_coli: p23 = ',p23
          write(*,111)' D0ms4ir4_coli: p34 = ',p34
          write(*,111)' D0ms4ir4_coli: p14 = ',p14
          write(*,111)' D0ms4ir4_coli: p24 = ',p24
          write(*,111)' D0ms4ir4_coli: p13 = ',p13
          write(*,111)' D0ms4ir4_coli: m12 = ',m12
          write(*,111)' D0ms4ir4_coli: m22 = ',m22
          write(*,111)' D0ms4ir4_coli: m32 = ',m32
          write(*,111)' D0ms4ir4_coli: m42 = ',m42
          flag(5)=.false.
        endif
#endif
 
* according to W.Beenakker and A.Denner, Nucl. Phys. B338 (1990) 349  (ie) 
c Ellis, Zanderighi, JHEP 0802 (2008) 002
c     box16 p32=0, p22=0  m2=m3=m4
c Denner, Dittmaier (3.71)

        D0ms4ir4_coli = (2d0*cln_coli(-q24/(mm32),-1d0)
     &     *cln_coli(-q13/sqrt(mm32*muir2),-1d0)
     &     - 3d0*pi2_6 
#ifdef SING
     &      + colishiftms2
#endif
     &      )
     &     /(q24*q13)


      elseif(nsoft.eq.0)then

* no soft singularity
*
*                   m2            
*          0 ----------------  0  
*                |  2   |         
*             m2 |1    3| m2      
*                |   4  |         
*          0 ----------------  0  
*                   m2            
*

#ifdef CHECK
        if(m12.ne.m22.or.m12.ne.m32.or.m12.ne.m42.or.
     &      p12.ne.cd0.or.p23.ne.cd0.or.p34.ne.cd0.or.p14.ne.cd0)
     &      then
          write(*,*) ' D0ms4ir4_coli:   (3.72)',
     &        '    soft regular case: ',
     &        '    wrong arguments'
          write(*,111)' D0ms4ir4_coli: p12 = ',p12
          write(*,111)' D0ms4ir4_coli: p23 = ',p23
          write(*,111)' D0ms4ir4_coli: p34 = ',p34
          write(*,111)' D0ms4ir4_coli: p14 = ',p14
          write(*,111)' D0ms4ir4_coli: p24 = ',p24
          write(*,111)' D0ms4ir4_coli: p13 = ',p13
          write(*,111)' D0ms4ir4_coli: m12 = ',m12
          write(*,111)' D0ms4ir4_coli: m22 = ',m22
          write(*,111)' D0ms4ir4_coli: m32 = ',m32
          write(*,111)' D0ms4ir4_coli: m42 = ',m42
          stop
        endif
#endif

#ifdef CHECK
        if(flag(6))then
          write(*,100)' D0ms4ir4_coli:   (3.72)',
     &        '    soft regular case ',
     &        '    not yet tested in physical process'
          write(*,111)' D0ms4ir4_coli: p12 = ',p12
          write(*,111)' D0ms4ir4_coli: p23 = ',p23
          write(*,111)' D0ms4ir4_coli: p34 = ',p34
          write(*,111)' D0ms4ir4_coli: p14 = ',p14
          write(*,111)' D0ms4ir4_coli: p24 = ',p24
          write(*,111)' D0ms4ir4_coli: p13 = ',p13
          write(*,111)' D0ms4ir4_coli: m12 = ',m12
          write(*,111)' D0ms4ir4_coli: m22 = ',m22
          write(*,111)' D0ms4ir4_coli: m32 = ',m32
          write(*,111)' D0ms4ir4_coli: m42 = ',m42
          flag(6)=.false.
        endif
#endif

c new case     27.07.09
c Denner, Dittmaier (3.72)
        mm12=m12*coliminfscale2
        D0ms4ir4_coli = 
     &      cln_coli(-p13/mm12,-1d0)*cln_coli(-p24/mm12,-1d0)-3d0*pi2_6
#ifdef SING
        D0ms4ir4_coli = D0ms4ir4_coli + colishiftms2
#endif

        D0ms4ir4_coli = D0ms4ir4_coli*2d0/(p13*p24)

      endif

      end

************************************************************************
