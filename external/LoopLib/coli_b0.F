#include "global_coli.h"
 
************************************************************************
*                                                                      *
*     Scalar 1-, and 2-point functions                                 *
*                                                                      *
************************************************************************
*                                                                      *
*     last changed  21.10.08  Ansgar Denner                            *
*                                                                      *
************************************************************************
* subroutines:                                                         *
* crootse,elminf2iv_coli,cerrB0_coli                                   *
* functions:                                                           *
* A0_coli                                                              *
* B0_coli, B1_coli, B11_coli, B00_coli , B111_coli, B001_coli          *
* BB0_coli, DB0_coli, DB1_coli, DB00_coli                              *
* Bn_coli,B00n_coli, B0000n_coli                                       *
* fpve_coli, yfpve_coli, xyfpve_coli, xlogxe_coli                      *
************************************************************************
      function A0_coli(m2)
************************************************************************
*     scalar one-loop-integral A0_coli                                 *
*-----------------------------------------------------------------------
*     02.10.89 Ansgar Denner      last changed 24.10.04                *
************************************************************************
      implicit   none
      complex*16 A0_coli,m2
 
      complex*16 elimminf2_coli

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"

      if (m2.eq.0d0) then
        A0_coli = 0d0
      else if (elimminf2_coli(m2).eq.0d0) then
        A0_coli = m2 * ( 1d0 - log(m2*coliminfscale2/muuv2) )
      else
        A0_coli = m2 * ( 1d0 - log(m2/muuv2) )
      end if
#ifdef SING
        A0_coli = A0_coli + m2 * deltauv
#endif
      end

************************************************************************
      subroutine crootse_coli(p2,m12,m22,x1,x2,y1,y2,r)
************************************************************************
*     roots of quadratic equation                                      *
*     p2*x*x+(m22-m12-p2)*x+m12     = p2*(x-x1)*(x-x2)                 *
*                                     p2*(x-1+y1)*(x-1+y2)             *
*----------------------------------------------------------------------*
*     21.07.93 Ansgar Denner      last changed  21.10.08 ad            *
************************************************************************
      implicit   none
      complex*16 p2,m12,m22
      complex*16 x1,x2,y1,y2,r
      complex*16 q12,q22,m1,m2
 
c      write(*,*) 'crootse_coli:in ',p2,m12,m22

c --->   case y1=y2=0  <-> B0_coli(m2,m,0) makes it difficult!
      m1 = sqrt(m12)
      m2 = sqrt(m22)
      r=sqrt((p2-(m1+m2)**2)*(p2-(m1-m2)**2))
c     r=sqrt(p2*p2+m12*m12+m22*m22-2d0*(m12*m22+m12*p2+m22*p2))
      q12=m12-m22+p2
      q22=m22-m12+p2
      x1=(q12+r)/(2d0*p2)
      x2=(q12-r)/(2d0*p2)

c      write(*,*) 'crootse_coli: ', r,x1,x2
c      write(*,*) 'crootse_coli: ', p2,q12,q22

      if (abs(x2).gt.abs(x1)) then
        x1=m12/(p2*x2)
      elseif (abs(x1).gt.abs(x2)) then
        x2=m12/(p2*x1)
      end if
      y2=(q22+r)/(2d0*p2)
      y1=(q22-r)/(2d0*p2)
      if (abs(y2).gt.abs(y1)) then
        y1=m22/(p2*y2)
      elseif (abs(y1).gt.abs(y2)) then
        y2=m22/(p2*y1)
      end if
c      x1 = x1 + dcmplx(0d0, abs(p2*x1)/p2*epsr)
c      x2 = x2 + dcmplx(0d0,-abs(p2*x2)/p2*epsr)
c      y1 = y1 + dcmplx(0d0,-abs(p2*y1)/p2*epsr)
c      y2 = y2 + dcmplx(0d0, abs(p2*y2)/p2*epsr)

c      write(*,*) 'crootse_coli out: ',x1,x2,y1,y2,r
      end

************************************************************************
      function B0_coli(p2,m12,m22)
************************************************************************
*     scalar one loop integral B0_coli                                 *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed 24.10.08 ad             *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      complex*16 B0_coli,fpve_coli,xlogxe_coli,cln_coli,eta2_coli
      complex*16 x1,x2,y1,y2,r
      logical    finarg
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
 
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' B0_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' B0_coli: p2  = ',p2
          write(*,111) ' B0_coli: m12 = ',m12
          write(*,111) ' B0_coli: m22 = ',m22
          stop
        endif
      endif
#endif

      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)

c      write(*,*) 'B0_coli in',p2,m12,m22,q2,mm12,mm22,finarg
c      write(*,*) 'B0_coli in',minfscale

      if (finarg) then
c---->  general case
        if (abs(q2).gt.calacc*abs(mm12+mm22)) then
          call crootse_coli(q2,mm12,mm22,x1,x2,y1,y2,r)
c          write(*,*) 'B0_coli x ',y1,y2,x1,x2
          if (abs(y1).gt.0.5d0.and.abs(y2).gt.0.5d0) then
            B0_coli= (-log(mm22/muuv2) - fpve_coli(1,x1,y1,1d0) 
     &                             - fpve_coli(1,x2,y2,-1d0))
c           write(*,*) 'B0_coli reg1',B0_coli
          elseif (abs(x1).lt.10d0.and.abs(x2).lt.10d0) then
            B0_coli = 2d0 - cln_coli(q2/muuv2,-1d0)
     &           + xlogxe_coli(-x1,-1d0) + xlogxe_coli(-x2,1d0) 
     &           - xlogxe_coli(y1,-1d0)  - xlogxe_coli(y2,1d0)
c olieta2 added 15.09.02
            if (mm12.ne.0d0) B0_coli = B0_coli 
     &           - eta2_coli(-x1,-x2,mm12/q2,-1d0,1d0,-q2)
     
c             write(*,*) 'B0_coli reg2',2d0, - cln_coli(q2/muuv2,-1d0),
c     &           + xlogxe_coli(-x1,-1d0), + xlogxe_coli(-x2,1d0), 
c     &           - xlogxe_coli(y1,-1d0),  - xlogxe_coli(y2,1d0),
c     &           - eta2_coli(-x1,-x2,mm12/q2,-1d0,1d0,-q2)
c             write(*,*) 'B0_coli reg2',B0_coli
          elseif (abs(x1).gt.0.5d0.and.abs(x2).gt.0.5d0) then
            B0_coli= (-log(mm12/muuv2) - fpve_coli(1,y1,x1,-1d0) 
     &                             - fpve_coli(1,y2,x2,1d0))
c           write(*,*) 'B0_coli reg3',B0_coli
          else
            call errB0_coli(q2,mm12,mm22,B0_coli,'B0_coli   ',12)
          endif
c---->  zero momentum
        elseif (abs(mm12-mm22).gt.calacc*abs(mm12+mm22)) then
          x2 = mm12/(mm12-mm22)
          y2 = mm22/(mm22-mm12)
          if (abs(y2).gt.0.5d0) then
            B0_coli = (-log(mm22/muuv2) - fpve_coli(1,x2,y2,mm22-mm12))
c           write(*,*) 'B0_coli sing1'
          else
            B0_coli = (-log(mm12/muuv2) - fpve_coli(1,y2,x2,mm12-mm22))
c           write(*,*) 'B0_coli sing2'
          endif
        else if (mm22.ne.0d0) then
c          write(*,*) 'B0_coli sing3'
          B0_coli = -log(mm22/muuv2) + q2/(6d0*mm22)
        else
          write (*,*) 'error'
          call errB0_coli(q2,mm12,mm22,B0_coli,'B0_coli   ',20)
        endif
      else
        if (q2.eq.mm12.and.mm22.eq.0d0.and.q2.ne.0d0) then
          B0_coli = 2d0 - log(mm12*coliminfscale2/muuv2)
        else if (q2.eq.mm22.and.mm12.eq.0d0.and.q2.ne.0d0) then
          B0_coli = 2d0 - log(mm22*coliminfscale2/muuv2)
        else if (q2.eq.0d0.and.mm12.eq.mm22.and.mm12.ne.0d0) then
          B0_coli = -log(mm12*coliminfscale2/muuv2)
        else if (q2.eq.0d0.and.mm12.eq.0d0.and.mm22.ne.0d0) then
          B0_coli = 1d0 - log(mm22*coliminfscale2/muuv2)
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.ne.0d0) then
          B0_coli = 1d0 - log(mm12*coliminfscale2/muuv2)
        else if (q2.ne.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          B0_coli = 2d0 - cln_coli(-q2*coliminfscale2/muuv2,-1d0)
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          B0_coli = log(muuv2/muir2)
#ifdef SING
          B0_coli = B0_coli-delta1ir

#endif
        else
          call errB0_coli(q2,mm12,mm22,B0_coli,'B0_coli  ',28)
        end if
      end if

c      write(*,*) 'B0_coli a ',B0_coli,deltauv,muuv2

#ifdef SING
      B0_coli = B0_coli + deltauv
#endif

c      write(*,*) 'B0_coli = ',B0_coli

      end

************************************************************************
      function BB0_coli(p2,m12,m22)
************************************************************************
*     scalar one loop integral (B0_coli(q2)-B0_coli(0))/q2             *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*     appears only multiplied with (m12-m22)                           *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed 21.10.08 ad             *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      complex*16 BB0_coli,fpve_coli,yfpve_coli,cln_coli
      complex*16 x1,x2,y1,y2,r
      logical    finarg
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
 
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' BB0_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' BB0_coli: p2  = ',p2
          write(*,111) ' BB0_coli: m12 = ',m12
          write(*,111) ' BB0_coli: m22 = ',m22
          stop
        endif
      endif
#endif

c      write(*,*)  'BB0 in ',p2,m1,m2
      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)
c      write(*,*)  'BB0 in ',p2,m1,m2,q2,mm1,mm2,finarg
      if(mm12.ne.cd0.or.mm22.ne.cd0)then
        if (finarg) then
c---->  general case
          if (abs(q2).gt.calacc*abs(mm12+mm22)) then
            call crootse_coli(q2,mm12,mm22,x1,x2,y1,y2,r)

c      write(*,*)  'BB0_coli x ',x1,x2,y1,y2,r
c       write(*,*) 'BB0_coli ',q2,mm12,mm22
c       write(*,*) 'BB0_coli ',x1,x2,y1,y2,r
c       write(*,*) 'BB0_coli = ',-fpve_coli(1,x2,y2,-1d0)/q2
c       write(*,*) 'BB0_coli = ',-fpve_coli(1,x1,y1,1d0)/q2
c       write(*,*) 'BB0_coli = ',-fpve_coli(1,y1,x1,-1d0)/q2
c       write(*,*) 'BB0_coli = ',-fpve_coli(1,y2,x2,1d0)/q2
c       write(*,*) 'BB0_coli = ',(-0.5d0/x1 - yfpve_coli(2,x1,y1,1d0)/x1
c    &             -0.5d0/x2 - yfpve_coli(2,x2,y2,-1d0)/x2 )/(mm22-mm12)
c       write(*,*) 'BB0_coli = ',(-0.5d0 - yfpv(1,x1,y1)
c    &             -0.5d0/x2 - yfpve_coli(2,x2,y2,-1d0)/x2 )/(mm22-mm12)
c       write(*,*) 'BB0_coli = ',(-0.5d0/x1 - yfpv(2,x1,y1)/x1
c    &             - 0.5d0 - yfpve_coli(1,x2,y2,-1d0))/(mm22-mm12)
c       write(*,*) 'BB0_coli = ',(-0.5d0 - yfpv(1,x1,y1)
c    &             - 0.5d0 - yfpve_coli(1,x2,y2,-1d0))/(mm22-mm12)
c       write(*,*) 'BB0_coli = ',(.5d0 +
c    &      (1d0-2d0*x1)*fpve_coli(2,x1,y1,1d0))/x1/x1/q2
c       write(*,*) 'BB0_coli = ',
c    &      (2d0 + (1d0-2d0*x1)*fpve_coli(0,x1,y1,1d0))/q2

            if (abs(x1+y2).lt.calacc*abs(x1+x2).or.
     &          abs(x2+y1).lt.calacc*abs(x1+x2)) then
              if(m12*m22.ne.0d0) then
                BB0_coli = 1d0 + (m12*cln_coli(-p2/m12,-1d0)
     &              - m22*cln_coli(-p2/m22,-1d0))/(m22-m12)
              else if (m12.ne.0d0) then
                BB0_coli = 1d0 - cln_coli(-p2/m12,-1d0)
              else if (m22.ne.0d0) then
                BB0_coli = 1d0 - cln_coli(-p2/m22,-1d0)
              else
                write(*,*) 'BB0 ',p2,m12,m22,q2,mm12,mm22,finarg
                call errB0_coli(q2,m12,m22,BB0_coli,'BB0_coli  ',24)              
              end if
              write(*,*) 'BB0: case not yet tested'
c     call errB0_coli(q2,mm12,mm22,BB0_coli,'BB0_coli  ',24)
            else if (abs(x1).lt.calacc*abs(x1+x2)) then
              BB0_coli = -fpve_coli(1,x2,y2,-1d0)/q2

c     write(*,*) 'BB0_coli reg1'

            else if (abs(x2).lt.calacc*abs(x1+x2)) then 
              BB0_coli = -fpve_coli(1,x1,y1,1d0)/q2

c     write(*,*) 'BB0_coli reg2'

            else if (abs(y2).lt.calacc*abs(y1+y2)) then
              BB0_coli = -fpve_coli(1,y1,x1,-1d0)/q2

c     write(*,*) 'BB0_coli reg3'

            else if (abs(y1).lt.calacc*abs(y1+y2)) then
              BB0_coli = -fpve_coli(1,y2,x2,1d0)/q2

c     write(*,*) 'BB0_coli reg4'

c     :  mm22-mm12 = q2*(1-x1-x2)
            else if (abs(y1-x2).gt.sqrt(calacc)*abs(y1+x2).or.
     &            abs(y2-x1).gt.sqrt(calacc)*abs(y2+x1)) then
              if(abs(x1).gt.10d0) then
                BB0_coli = -0.5d0/x1 - yfpve_coli(2,x1,y1,1d0)/x1
c     write(*,*) 'BB0_coli reg15'
              else
                BB0_coli = -0.5d0 - yfpve_coli(1,x1,y1,1d0)
c     write(*,*) 'BB0_coli reg16'
              endif
              if(abs(x2).gt.10d0) then
                BB0_coli = BB0_coli - 0.5d0/x2 
     &              - yfpve_coli(2,x2,y2,-1d0)/x2

c     write(*,*) 'BB0_coli reg25'

              else
                BB0_coli = BB0_coli - 0.5d0 - yfpve_coli(1,x2,y2,-1d0)

c     write(*,*) 'BB0_coli reg26'

              endif
              BB0_coli = BB0_coli/(mm22-mm12)
            else
              if (abs(x1).gt.10d0) then
                BB0_coli = (.5d0 + (1d0-2d0*x1)*fpve_coli(2,x1,y1,1d0))
     &              /(x1*x1*q2)

c     write(*,*) 'BB0_coli reg7'

              elseif (y1.ne.0d0.and.x1.ne.0d0) then
                BB0_coli = (2d0 + (1d0-2d0*x1)*fpve_coli(0,x1,y1,1d0))
     &              /q2

c     write(*,*) 'BB0_coli reg8'

c     elseif (y2.ne.0d0.and.x2.ne.0d0) then
c     BB0_coli = (2d0 + (1d0-2d0*x2)*fpve_coli(0,x2,y2,-1d0))/q2
              else
                call errB0_coli(q2,mm12,mm22,BB0_coli,'BB0_coli  ',16)
              endif
            endif
c---- >  zero momentum
          elseif (abs(mm12-mm22).gt.calacc*abs(mm12+mm22)) then
            x2 = mm12/(mm12-mm22)
            y2 = mm22/(mm22-mm12)
            if (abs(x2).lt.10d0) then
              BB0_coli = (0.5d0 + yfpve_coli(1,x2,y2,mm22-mm12))
     &            /(mm12-mm22)
c     write(*,*) 'BB0_coli sing1'
            else
              BB0_coli = (y2*fpve_coli(2,x2,y2,mm22-mm12) + 0.5d0)/mm12
c     write(*,*) 'BB0_coli sing2'
            endif
          else if (mm12.ne.0d0) then
            BB0_coli = 1d0/6d0/mm12
c     write(*,*) 'BB0_coli sing3'
          else 
            call errB0_coli(q2,mm12,mm22,BB0_coli,'BB0_coli  ',20)
          endif
        else
          if (q2.eq.mm12.and.mm22.eq.0d0.and.q2.ne.0d0) then
            BB0_coli = 1d0/q2
          else if (q2.eq.mm22.and.mm12.eq.0d0.and.q2.ne.0d0) then
            BB0_coli = 1d0/q2
          else if (q2.eq.0d0.and.mm12.eq.mm22.and.mm12.ne.0d0) then
            BB0_coli = 1d0/(6d0*mm12)
          else if (q2.eq.0d0.and.mm12.eq.0d0.and.mm22.ne.0d0) then
            BB0_coli = 1d0/(2d0*mm22)
          else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.ne.0d0) then
            BB0_coli = 1d0/(2d0*mm12)
          else
            call errB0_coli(q2,mm12,mm22,BB0_coli,'BB0_coli  ',28)
          end if
        endif
      else
        if(m12*m22.ne.0d0) then
          BB0_coli = (1d0 + (m12*cln_coli(-p2/(m12*coliminfscale2),-1d0)
     &        - m22*cln_coli(-p2/(m22*coliminfscale2),-1d0))
     &        /(m22-m12))/p2
        else if (m12.ne.0d0) then
          BB0_coli = (1d0 - cln_coli(-p2/(m12*coliminfscale2),-1d0))/p2
        else if (m22.ne.0d0) then
          BB0_coli = (1d0 - cln_coli(-p2/(m22*coliminfscale2),-1d0))/p2
        else
c     write(*,*) 'BB0_coli ',p2,m12,m22,q2,mm1,mm2,finarg
          call errB0_coli(q2,m12,m22,BB0_coli,'BB0_coli  ',24)              
        end if
#ifdef SING
          BB0_coli = BB0_coli + delta1ir
#endif
      end if
c     write(*,*) 'BB0_coli ',BB0_coli   
      end

************************************************************************
      function DB0_coli(p2,m12,m22)
************************************************************************
*     derivative of scalar one loop integral B0_coli                   *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  28.10.08 ad            *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      complex*16 DB0_coli,fpve_coli,yfpve_coli,xyfpve_coli
      complex*16 x1,x2,y1,y2,r
      logical    finarg
 
#ifdef CHECK
      logical flag
      data    flag /.true./
      save    flag  
#endif

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
  
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' DB0_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' DB0_coli: p2  = ',p2
          write(*,111) ' DB0_coli: m12 = ',m12
          write(*,111) ' DB0_coli: m22 = ',m22
          stop
        endif
      endif
#endif

c      write(*,*) 'DB0f in ',p2,m12,m22
      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)

      if (finarg) then
c IR-singular case
        if(p2.eq.m12+m22.and.m12*m22.eq.0d0) then
          DB0_coli = -(1d0 + 0.5d0*log(muir2/(m12+m22)))/p2
#ifdef SING
          DB0_coli = DB0_coli -0.5d0*delta1ir/p2
#endif

c---->  general case
        else if (abs(q2).gt.calacc*abs(mm12+mm22)) then
          call crootse_coli(q2,mm12,mm22,x1,x2,y1,y2,r)
c            write(*,*) 'DB0_coli reg',x1,x2
c:  r = q2*(x1-x2)
c changed 7.8.02  B0_coli(mm2,0,m) was wrong on pc  (xi=0,yi=1)
          if (abs(x1-x2).gt.sqrt(calacc)*abs(x1+x2)) then

c            write(*,*) 'DB0_coli reg1',x1,x2

            if (abs(x1).gt.0.1d0.or.abs(x2).gt.0.1d0) then
              DB0_coli = (yfpve_coli(1,x2,y2,-1d0) 
     &            - yfpve_coli(1,x1,y1,1d0))/r

c              write(*,*) 'DB0_coli reg11',DB0_coli

            else  
c added 20.09.04   for small x1,x2
              DB0_coli = (xyfpve_coli(0,x2,y2,-1d0) 
     &            - xyfpve_coli(0,x1,y1,1d0))/r
     &                 -1d0/q2

c              write(*,*) 'DB0_coli reg12',DB0_coli

            end if

c            write(*,*) 'DB0_coli reg1',
c     &       yfpve_coli(1,x2,y2,-1d0), - yfpve_coli(1,x1,y1,1d0)
c            write(*,*) 'DB0_coli reg1',
c     &       x2,y2,-1d0,x1,y1,1d0
c            write(*,*) 'DB0_coli reg1',DB0_coli

          elseif (abs(x1).gt.10d0) then
            DB0_coli = -(.5d0 + (1d0-2d0*x1)*fpve_coli(2,x1,y1,1d0))
     &          /(q2*x1*x1)

c           write(*,*) 'DB0_coli reg2'

          elseif (abs(y1).gt.calacc.and.abs(x1).gt.calacc) then
            DB0_coli = -( 2d0 + (1d0-2d0*x1)*fpve_coli(0,x1,y1,1d0) )/q2

c           write(*,*) 'DB0_coli reg3'

          elseif (mm12*mm22.eq.0d0.and.q2.eq.(mm12+mm22)) then
            DB0_coli = -(1d0 + 0.5d0*log(muir2/(mm12+mm22)))/q2

c            write(*,*) 'DB0_coli reg4'

c added 30.09.04 for very small width
          else if (abs(x1-x2).gt.calacc*abs(x1+x2)) then

c            write(*,*) 'DB0_coli reg5',x1,x2

            DB0_coli = (yfpve_coli(1,x2,y2,-1d0) 
     &          - yfpve_coli(1,x1,y1,1d0))/r
c end added 30.09.04
          else
            call errB0_coli(q2,mm12,mm22,DB0_coli,'DB0_coli  ',12)
          end if
c---->  zero momentum
        elseif (abs(mm12-mm22).gt.calacc*abs(mm12+mm22)) then

c          write(*,*) 'DB0_coli sing'

          x2 = mm12/(mm12-mm22)
          y2 = mm22/(mm22-mm12)
          if (abs(x2).lt.10d0) then
            DB0_coli = (0.5d0 + yfpve_coli(1,x2,y2,mm22-mm12))
     &          /(mm12-mm22)

c            write(*,*) 'DB0_coli sing1'

          else
            DB0_coli = (y2*fpve_coli(2,x2,y2,mm22-mm12) + 0.5d0)/mm12

c            write(*,*) 'DB0_coli sing2'

          endif
        else if (mm12.ne.0d0) then
          DB0_coli = 1d0/(6d0*mm12)

c          write(*,*) 'DB0_coli sing3'

        else
          call errB0_coli(q2,mm12,mm22,DB0_coli,'DB0_coli  ',20)
        endif
        
      else
        if (q2.eq.mm12.and.mm22.eq.0d0.and.q2.ne.0d0) then
          DB0_coli = -(1d0 + .5d0*log(muir2/(mm12*coliminfscale2)))/q2
#ifdef SING
          DB0_coli = DB0_coli -0.5d0*delta1ir/q2
#endif
        else if (q2.eq.mm22.and.mm12.eq.0d0.and.q2.ne.0d0) then
          DB0_coli = -(1d0 + .5d0*log(muir2/(mm22*coliminfscale2)))/q2
#ifdef SING
          DB0_coli = DB0_coli -0.5d0*delta1ir/q2
#endif
        else if (q2.eq.0d0.and.mm12.eq.mm22.and.mm12.ne.0d0) then
          DB0_coli = 1d0/(6d0*mm12)
        else if (q2.eq.0d0.and.mm12.eq.0d0.and.mm22.ne.0d0) then
          DB0_coli = 1d0/(2d0*mm22)
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.ne.0d0) then
          DB0_coli = 1d0/(2d0*mm12)
        else if (q2.ne.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          DB0_coli = -1d0/(q2)
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          DB0_coli = cd0
#ifdef CHECK
          if(flag) then
            write(*,*) ' DB0_coli for three vanishing arguments! '
            flag=.false.
          endif
#endif
        else
          call errB0_coli(q2,mm12,mm22,DB0_coli,'DB0_coli  ',28)
        end if
      end if

c      write(*,*) 'DB0_coli =',DB0_coli
      end

************************************************************************
      function B1_coli(p2,m12,m22)
************************************************************************
*     vector one loop integral B1_coli                                 *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  24.10.08 ad            *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      logical    finarg
      complex*16 B1_coli,fpve_coli,xlogxe_coli,cln_coli,eta2_coli
      complex*16 x1,x2,y1,y2,r
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
 
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' B1_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' B1_coli: p2  = ',p2
          write(*,111) ' B1_coli: m12 = ',m12
          write(*,111) ' B1_coli: m22 = ',m22
          stop
        endif
      endif
#endif

c      integer    testout
c      testout=70
 
      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)
      if (finarg) then
c---->  general case
        if (abs(q2).gt.calacc*abs(mm12+mm22)) then
          call crootse_coli(q2,mm12,mm22,x1,x2,y1,y2,r)
c          write(testout,*) 'B1_coli : x',x1,x2,y1,y2,r
          if (abs(y1).gt.0.5d0.and.abs(y2).gt.0.5d0) then
            B1_coli = -(-log(mm22/muuv2) 
     &           - fpve_coli(2,x1,y1,1d0) - fpve_coli(2,x2,y2,-1d0))/2d0
          elseif (abs(x1).lt.10d0.and.abs(x2).lt.10d0) then
            B1_coli = 1d0 - cln_coli(q2/muuv2,-1d0)
     &           + x1*xlogxe_coli(-x1,-1d0) + x1 
     &          + x2*xlogxe_coli(-x2,1d0) + x2
     &           - (1d0+x1)*xlogxe_coli(y1,-1d0) 
     &          - (1d0+x2)*xlogxe_coli(y2,1d0)
c olieta2 added 15.09.02
            if (mm12.ne.0d0) B1_coli = B1_coli 
     &           - eta2_coli(-x1,-x2,mm12/q2,-1d0,1d0,-q2)
            B1_coli = -B1_coli/2d0
          elseif (abs(x1).gt.0.5d0.and.abs(x2).gt.0.5d0) then
            if (abs(x1).gt.10) then
              B1_coli = ( log(mm12/muuv2) 
     &            + fpve_coli(2,x1,y1,1d0) - fpve_coli(0,x1,y1,1d0)
     &            + ((1d0+x2)*fpve_coli(1,y2,x2,1d0) + .5d0)  )/2d0
            else 
              B1_coli = ( log(mm12/muuv2) 
     &            + ((1d0+x1)*fpve_coli(1,y1,x1,-1d0) + .5d0)
     &            + fpve_coli(2,x2,y2,-1d0) - fpve_coli(0,x2,y2,-1d0))/2d0
            end if
          else
            call errB0_coli(q2,mm12,mm22,B1_coli,'B1_coli   ',12)
          endif
c---->  zero momentum
        elseif (abs(mm12-mm22).gt.calacc*abs(mm12+mm22)) then
          x2 = mm12/(mm12-mm22)
          y2 = mm22/(mm22-mm12)
          if(abs(y2).gt.0.5d0) then
            B1_coli = -(-log(mm22/muuv2) 
     &          - fpve_coli(2,x2,y2,mm22-mm12))/2d0
          else
            B1_coli = -(-log(mm12/muuv2) 
     &          - (1+x2)*fpve_coli(1,y2,x2,mm12-mm22) 
     &           - .5d0)/2d0
          endif
        else if (mm22.ne.0d0) then
          B1_coli = -(-log(mm22/muuv2))/2d0
        else
          call errB0_coli(q2,mm12,mm22,B1_coli,'B1_coli   ',20)
        endif
      else
        if (q2.eq.mm12.and.mm22.eq.0d0.and.q2.ne.0d0) then
          B1_coli = -1.5d0 + 0.5d0*log(mm12*coliminfscale2/muuv2)
        else if (q2.eq.mm22.and.mm12.eq.0d0.and.q2.ne.0d0) then
          B1_coli = -0.5d0 + 0.5d0*log(mm22*coliminfscale2/muuv2)
        else if (q2.eq.0d0.and.mm12.eq.mm22.and.mm12.ne.0d0) then
          B1_coli = 0.5d0*log(mm12*coliminfscale2/muuv2)
        else if (q2.eq.0d0.and.mm12.eq.0d0.and.mm22.ne.0d0) then
          B1_coli = -0.25d0 + 0.5d0*log(mm22*coliminfscale2/muuv2)
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.ne.0d0) then
          B1_coli = -0.75d0 + 0.5d0*log(mm12*coliminfscale2/muuv2)
        else if (q2.ne.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          B1_coli = -1d0 + 0.5d0*cln_coli(-q2*coliminfscale2/muuv2,-1d0)
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          B1_coli = -.5d0*log(muuv2/muir2)
#ifdef SING
          B1_coli = B1_coli + .5d0*delta1ir
#endif
        else
          call errB0_coli(q2,mm12,mm22,B1_coli,'B1_coli  ',28)
          write(*,*) q2-mm22,mm12
        end if
      end if   

c      write(testout,*) ' B1_coli  ',p2,m12,m22
c      write(testout,*) ' B1_coli =',B1_coli

#ifdef SING
      B1_coli = B1_coli - 0.5d0 * deltauv
#endif
      end

************************************************************************
      function DB1_coli(p2,m12,m22)
************************************************************************
*     derivative of vector one loop integral B1_coli                   *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  28.10.08 ad            *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      complex*16 DB1_coli,fpve_coli,yfpve_coli,xyfpve_coli
      complex*16 x1,x2,y1,y2,r
      logical    finarg
 
#ifdef CHECK
      logical flag
      data    flag /.true./
      save    flag  
#endif

#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
 
c      integer    testout
     
c      testout = 20

#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' DB1_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' DB1_coli: p2  = ',p2
          write(*,111) ' DB1_coli: m12 = ',m12
          write(*,111) ' DB1_coli: m22 = ',m22
          stop
        endif
      endif
#endif

c      write(testout,*) 'DB1_coli in', p2,m12,m22
      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg) 
      if (finarg) then
c IR-singular case
c added 23.09.04
        if(p2.eq.m12.and.m22.eq.0d0) then
          DB1_coli = (3d0 + log(muir2/p2))/(2d0*p2)
#ifdef SING
          DB1_coli = DB1_coli + delta1ir/(2d0*p2)
#endif
c---->  general case
        else if (abs(q2).gt.calacc*abs(mm12+mm22)) then
          call crootse_coli(q2,mm12,mm22,x1,x2,y1,y2,r)

c            write(testout,*) 'DB1_coli x',x1,x2,y1,y2,r,x1+x2

          if (abs(x1-x2).gt.sqrt(calacc)*abs(x1+x2)) then
c            write(*,*) 'DB1_coli reg1',x1,x2
            if (abs(x1).gt.0.1d0.or.abs(x2).gt.0.1d0) then
              DB1_coli = -(yfpve_coli(2,x2,y2,-1d0) 
     &            - yfpve_coli(2,x1,y1,1d0))/r

c             write(testout,*) 'DB1_coli reg1a',DB1_coli

c added 20.09.04   for small x1,x2
            else 
              DB1_coli = -(x2*xyfpve_coli(0,x2,y2,-1d0)
     &            -x1*xyfpve_coli(0,x1,y1,1d0))
     &             /r    +(x1+x2-.5d0)/q2
c             writet(testout,*) 'DB1_coli reg1b',DB1_coli
            end if
c            write(*,*) 'DB1_coli reg1',DB1_coli
          elseif (abs(x1).gt.10d0) then
            DB1_coli = (2d0/3d0 + (2d0-3d0*x1)*fpve_coli(3,x1,y1,1d0))
     &              /(q2*x1*x1)
c changed 09.09.04
c            DB1_coli = dreal(DB1_coli)
          elseif (abs(y1).gt.calacc) then
            DB1_coli = ( 3d0/2d0 + (2d0-3d0*x1)*fpve_coli(1,x1,y1,1d0) )
     &          /q2
c changed 09.09.04
c            DB1_coli = dreal(DB1_coli)
          elseif (mm22.eq.0d0.and.q2.eq.mm12) then
            DB1_coli = (3d0 + log(muir2/q2))/(q2*2d0)
c---> ir singular part
c added 30.09.04 for very small width
          else if (abs(x1-x2).gt.calacc*abs(x1+x2)) then
c            write(*,*) 'DB0_coli reg5',x1,x2
            DB1_coli = -(yfpve_coli(2,x2,y2,-1d0) 
     &          - yfpve_coli(2,x1,y1,1d0))/r
c end added 30.09.04
          else
            call errB0_coli(q2,mm12,mm22,DB1_coli,'DB1_coli  ',12)
          end if
c---->  zero momentum
        elseif (abs(mm12-mm22).gt.calacc*abs(mm12+mm22)) then
          x2 = mm12/(mm12-mm22)
          y2 = mm22/(mm22-mm12)
          if (abs(x2).lt.10d0) then
            DB1_coli = -(1d0/3d0 + yfpve_coli(2,x2,y2,mm22-mm12))
     &          /(mm12-mm22)
          else
            DB1_coli = -(y2*fpve_coli(3,x2,y2,mm22-mm12) + 1d0/3d0)/mm12
          endif
        else if (mm12.ne.0d0) then
          DB1_coli  = -1d0/(12d0*mm12)
        else
          call errB0_coli(q2,mm12,mm22,DB1_coli,'DB1_coli  ',20)
        endif
      else
        if (q2.eq.mm12.and.mm22.eq.0d0.and.q2.ne.0d0) then
          DB1_coli = (3d0 + log(muir2/(mm12*coliminfscale2)))/(2d0*q2)
#ifdef SING
          DB1_coli = DB1_coli + delta1ir/(2d0*p2)
#endif
        else if (q2.eq.mm22.and.mm12.eq.0d0.and.q2.ne.0d0) then
          DB1_coli = -0.5d0/q2
        else if (q2.eq.0d0.and.mm12.eq.mm22.and.mm12.ne.0d0) then
          DB1_coli = -1d0/(12d0*mm12)
        else if (q2.eq.0d0.and.mm12.eq.0d0.and.mm22.ne.0d0) then
          DB1_coli = -1d0/(6d0*mm22)
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.ne.0d0) then
          DB1_coli = -1d0/(3d0*mm12)
        else if (q2.ne.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          DB1_coli = 1d0/(2d0*q2)
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          DB1_coli = cd0
#ifdef CHECK
          if(flag) then
            write(*,*) ' DB0_coli for three vanishing arguments! '
            flag=.false.
          endif
#endif
        else
          call errB0_coli(q2,mm12,mm22,DB1_coli,'DB1_coli  ',28)
        end if
      end if

c      write(*,*) 'DB1 ',db1f,p2,m1,m2,q2,mm1,mm2
      end

************************************************************************
      function B11_coli(p2,m12,m22)
************************************************************************
*     tensor one loop integral B11_coli                                *
*     p2 outer momentum, m1,m2 inner masses squared                    *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  24.10.08 ad            *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      complex*16 B11_coli,fpve_coli,yfpve_coli,xlogxe_coli,cln_coli,
     &    eta2_coli
      complex*16 x1,x2,y1,y2,r
      logical    finarg
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
 
 
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' B11_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' B11_coli: p2  = ',p2
          write(*,111) ' B11_coli: m12 = ',m12
          write(*,111) ' B11_coli: m22 = ',m22
          stop
        endif
      endif
#endif

      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)
      if (finarg) then
c---->  general case
        if (abs(q2).gt.calacc*abs(mm12+mm22)) then
          call crootse_coli(q2,mm12,mm22,x1,x2,y1,y2,r)
          if (abs(y1).gt.0.5d0.and.abs(y2).gt.0.5d0) then
            B11_coli = (-log(mm22/muuv2) 
     &           - fpve_coli(3,x1,y1,1d0) - fpve_coli(3,x2,y2,-1d0))/3d0
          elseif (abs(x1).lt.10d0.and.abs(x2).lt.10d0) then
            B11_coli = 2d0/3d0 - cln_coli(q2/muuv2,-1d0)
     &           + x1**2*xlogxe_coli(-x1,-1d0) + x1*x1 + x1/2d0
     &           + x2**2*xlogxe_coli(-x2, 1d0) + x2*x2 + x2/2d0
     &           - (1d0+x1+x1*x1)*xlogxe_coli(y1,-1d0) 
     &           - (1d0+x2+x2*x2)*xlogxe_coli(y2,1d0)
c olieta2 added 15.09.02
            if (mm12.ne.0d0) B11_coli = B11_coli 
     &           - eta2_coli(-x1,-x2,mm12/q2,-1d0,1d0,-q2)
            B11_coli = B11_coli/3d0
          elseif (abs(x1).gt.0.5d0.and.abs(x2).gt.0.5d0) then
            if (abs(x1).gt.10) then
              B11_coli = ( -log(mm12/muuv2)
     &              - fpve_coli(3,x1,y1,1d0) + fpve_coli(0,x1,y1,1d0)
     &            - ((x2*x2+x2+1d0)*fpve_coli(1,y2,x2,1d0) 
     &                  + x2/2d0 + 2d0/3d0)
     &            )/3d0
            else
              B11_coli = ( -log(mm12/muuv2)
     &            - ((x1*x1+x1+1d0)*fpve_coli(1,y1,x1,-1d0) 
     &                  + x1/2d0 + 2d0/3d0)
     &              + fpve_coli(3,x2,y2,-1d0) - fpve_coli(0,x2,y2,-1d0)
     &            )/3d0
            end if 
          else
            call errB0_coli(q2,mm12,mm22,B11_coli,'B11_coli  ',12)
          endif
c---->  zero momentum
        elseif (abs(mm12-mm22).gt.calacc*abs(mm12+mm22)) then
          x2 = mm12/(mm12-mm22)
          y2 = mm22/(mm22-mm12)
          if(abs(y2).gt.0.5d0) then
            B11_coli = (-log(mm22/muuv2) - fpve_coli(3,x2,y2,mm22-mm12))
     &          /3d0
          else
            B11_coli = (-log(mm12/muuv2) 
     &           + (1d0+x2+x2*x2)*yfpve_coli(0,x2,y2,mm22-mm12)
     &           + x2*x2 + x2/2d0 + 1d0/3d0)/3d0
          endif
        else if(mm22.ne.0d0) then
          B11_coli = -log(mm22/muuv2)/3d0
        else
          call errB0_coli(q2,mm12,mm22,B11_coli,'B11_coli  ',20)
        endif
      else
        if (q2.eq.mm12.and.mm22.eq.0d0.and.q2.ne.0d0) then
          B11_coli = 11d0/9d0 - log(mm12*coliminfscale2/muuv2)/3d0
        else if (q2.eq.mm22.and.mm12.eq.0d0.and.q2.ne.0d0) then
          B11_coli = 2d0/9d0 - log(mm22*coliminfscale2/muuv2)/3d0
        else if (q2.eq.0d0.and.mm12.eq.mm22.and.mm12.ne.0d0) then
          B11_coli = -log(mm12*coliminfscale2/muuv2)/3d0
        else if (q2.eq.0d0.and.mm12.eq.0d0.and.mm22.ne.0d0) then
          B11_coli = 1d0/9d0 - log(mm22*coliminfscale2/muuv2)/3d0
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.ne.0d0) then
          B11_coli = 11d0/18d0 - log(mm12*coliminfscale2/muuv2)/3d0
        else if (q2.ne.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          B11_coli = 13d0/18d0 - cln_coli(-q2*coliminfscale2/muuv2,-1d0)
     &        /3d0
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          B11_coli = log(muuv2/muir2)/3d0
#ifdef SING
          B11_coli = B11_coli - delta1ir/3d0
#endif
        else
          call errB0_coli(q2,mm12,mm22,B11_coli,'B11_coli  ',28)
          write(*,*) q2-mm22,mm12
        end if
      end if   
#ifdef SING
      B11_coli = B11_coli + deltauv/3d0
#endif
      end

************************************************************************
      function B00_coli(p2,m12,m22)
************************************************************************
*     one loop integral B00_coli                                       *
*     p2 outer momentum                                                *
*     m12,m22 inner masses squared                                     *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  21.10.08 ad            *
************************************************************************
      implicit none
      complex*16 A0_coli,B0_coli,B1_coli,B00_coli
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      logical    finarg
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
  
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' B00_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' B00_coli: p2  = ',p2
          write(*,111) ' B00_coli: m12 = ',m12
          write(*,111) ' B00_coli: m22 = ',m22
          stop
        endif
      endif
#endif

c      write(*,*) 'B00_coli ',p2,m12,m22
      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)
c      write(*,*) 'B00_coli ',p2,m12,m22,q2,mm12,mm22,finarg
      if (finarg) then
        if (abs(q2)/(abs(q2)+abs(mm12+mm22)).gt.calacc) then
          B00_coli = (q2+mm12-mm22)*B1_coli(q2,mm12,mm22) 
     &        + 2d0*mm12*B0_coli(q2,mm12,mm22)
     &        + A0_coli(mm22) + mm12+mm22-q2/3d0
          B00_coli = B00_coli/6d0
c     write(*,*) 'B00_coli ',(q2+mm12-mm22)*B1_coli(q2,mm12,mm22), 
c     &         + 2d0*mm12*B0_coli(q2,mm12,mm22),
c     &         + A0_coli(mm22), + mm12+mm22-q2/3d0
c     write(*,*) 'B00_coli ',B00_coli
        else
          B00_coli = (2d0*(mm22*B0_coli(dcmplx(0d0),mm12,mm22) 
     &        + A0_coli(mm12)) + mm12 + mm22)/8d0
        end if
      else
        B00_coli = 0d0
      end if
      end

************************************************************************
      function DB00_coli(p2,m12,m22)
************************************************************************
*     derivative of scalar one loop integral b20 = B00_coli            *
*     p2 outer momentum                                                *
*     m12,m22 inner masses squared                                     *
*     NOT REALLY TESTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*
*----------------------------------------------------------------------*
*     17.10.05 Ansgar Denner      last changed  21.10.08 ad            *
************************************************************************
      implicit none
      complex*16 DB0_coli,B1_coli,DB1_coli,DB00_coli
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      logical    finarg
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
 
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' DB00_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' DB00_coli: p2  = ',p2
          write(*,111) ' DB00_coli: m12 = ',m12
          write(*,111) ' DB00_coli: m22 = ',m22
          stop
        endif
      endif
#endif

c      write(*,*) 'B00_coli ',p2,m12,m22
      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)
c      write(*,*) 'B00_coli ',p2,m12,m22
c      if (finarg) then
c        if (abs(q2)/(abs(q2)+abs(mm12+mm22)).gt.calacc) then
          DB00_coli = B1_coli(q2,mm12,mm22) 
     &        +   (q2+mm12-mm22)*DB1_coli(q2,mm12,mm22) 
     &        + 2d0*mm12*DB0_coli(q2,mm12,mm22)
     &        -1d0/3d0
          DB00_coli = DB00_coli/6d0
c        else
c          write(*,*) 'DB00_coli ',p2,m12,m22
c          write(*,*) 'DB00_coli not implemented 1'
c          stop
c        end if
c      else
c        write(*,*) 'DB00_coli ',p2,m12,m22
c        write(*,*) 'DB00_coli not implemented 2'
c        stop
c      end if
      end

************************************************************************
      function B111_coli(p2,m12,m22)
************************************************************************
*     one loop integral B111_coli                                      *
*     p2 outer momentum, m12,m22 inner masses                          *
*----------------------------------------------------------------------*
*     05.02.04 Ansgar Denner      last changed 24.10.08 ad             *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      complex*16 B111_coli,fpve_coli,yfpve_coli,xlogxe_coli,cln_coli,
     &    eta2_coli
      complex*16 x1,x2,y1,y2,r
      logical    finarg
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
 
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' B111_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' B111_coli: p2  = ',p2
          write(*,111) ' B111_coli: m12 = ',m12
          write(*,111) ' B111_coli: m22 = ',m22
          stop
        endif
      endif
#endif

      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)
      if (finarg) then
c---->  general case
        if (abs(q2).gt.calacc*abs(mm12+mm22)) then
          call crootse_coli(q2,mm12,mm22,x1,x2,y1,y2,r)
          if (abs(y1).gt.0.5d0.and.abs(y2).gt.0.5d0) then
            B111_coli = -(-log(mm22/muuv2) 
     &           - fpve_coli(4,x1,y1,1d0) - fpve_coli(4,x2,y2,-1d0))/4d0
c           write(*,*) 'B111_coli reg1'
          elseif (abs(x1).lt.10d0.and.abs(x2).lt.10d0) then
            B111_coli = 1d0/2d0 - cln_coli(q2/muuv2,-1d0)
     &           + x1**3*xlogxe_coli(-x1,-1d0)
     &          + x1*x1*x1 + x1*x1/2d0 + x1/3d0
     &           + x2**3*xlogxe_coli(-x2,1d0) 
     &          + x2*x2*x2 + x2*x2/2d0 + x2/3d0
     &           - (1d0+x1+x1*x1+x1*x1*x1)*xlogxe_coli(y1,-1d0) 
     &           - (1d0+x2+x2*x2+x2*x2*x2)*xlogxe_coli(y2,1d0)
c olieta2 added 15.09.02
            if (mm12.ne.0d0) B111_coli = B111_coli 
     &           - eta2_coli(-x1,-x2,mm12/q2,-1d0,1d0,-q2)
            B111_coli = -B111_coli/4d0

c           write(*,*) 'B111_coli reg2'

          elseif (abs(x1).gt.0.5d0.and.abs(x2).gt.0.5d0) then
            if (abs(x1).gt.10) then
              B111_coli = ( log(mm12/muuv2)
     &              + fpve_coli(4,x1,y1,1d0) - fpve_coli(0,x1,y1,1d0)
     &            + ((x2*x2*x2+x2*x2+x2+1)*fpve_coli(1,y2,x2,1d0) 
     &            + x2*x2/2d0 + 2d0/3d0*x2 + 3d0/4d0) )/4d0
            else
              B111_coli = ( log(mm12/muuv2)
     &            + ((x1*x1*x1+x1*x1+x1+1)*fpve_coli(1,y1,x1,-1d0) 
     &            + x1*x1/2d0 + 2d0/3d0*x1 + 3d0/4d0)
     &            + fpve_coli(4,x2,y2,-1d0) 
     &            - fpve_coli(0,x2,y2,-1d0))/4d0
            end if

c           write(*,*) 'B111_coli reg3'

          else
            call errB0_coli(q2,mm12,mm22,B111_coli,'B111_coli ',12)
          endif
c---->  zero momentum
        elseif (abs(mm12-mm22).gt.calacc*abs(mm12+mm22)) then
          x2 = mm12/(mm12-mm22)
          y2 = mm22/(mm22-mm12)
          if(abs(y2).gt.0.5d0) then
            B111_coli = (log(mm22/muuv2) 
     &          + fpve_coli(4,x2,y2,mm22-mm12))/4d0

c           write(*,*) 'B111_coli sing1'

          else
            B111_coli =(log(mm12/muuv2) 
     &           - (1d0+x2+x2*x2+x2*x2*x2)*yfpve_coli(0,x2,y2,mm22-mm12)
     &           - x2*x2*x2 - x2*x2/2d0 - x2/3d0 - 1d0/4d0)/4d0

c           write(*,*) 'B111_coli sing2'

          endif
        else if (mm22.ne.0d0) then
          B111_coli = log(mm22/muuv2)/4d0

c         write(*,*) 'B111_coli sing3'

        else
          call errB0_coli(q2,mm12,mm22,B111_coli,'B111_coli ',20)
        endif
      else
        if (q2.eq.mm12.and.mm22.eq.0d0.and.q2.ne.0d0) then
          B111_coli = -25d0/24d0 + 0.25d0*log(mm12*coliminfscale2/muuv2)
        else if (q2.eq.mm22.and.mm12.eq.0d0.and.q2.ne.0d0) then
          B111_coli = -0.125d0 + 0.25d0*log(mm22*coliminfscale2/muuv2)
        else if (q2.eq.0d0.and.mm12.eq.mm22.and.mm12.ne.0d0) then
          B111_coli = 0.25d0*log(mm12*coliminfscale2/muuv2)
        else if (q2.eq.0d0.and.mm12.eq.0d0.and.mm22.ne.0d0) then
          B111_coli = -0.0625d0 + 0.25d0*log(mm22*coliminfscale2/muuv2)
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.ne.0d0) then
          B111_coli = -25d0/48d0 + 0.25d0*log(mm12*coliminfscale2/muuv2)
        else if (q2.ne.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          B111_coli = -7d0/12d0 
     &        + 0.25d0*cln_coli(-q2*coliminfscale2/muuv2,-1d0)
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          B111_coli = -.25d0 * log(muuv2/muir2)
#ifdef SING
          B111_coli = B111_coli + .25d0 * delta1ir
#else
          B111_coli = cd0
#endif
        else
          call errB0_coli(q2,mm12,mm22,B111_coli,'B111_coli  ',28)
          write(*,*) q2-mm22,mm12
        end if
      end if   
#ifdef SING
      B111_coli = B111_coli -.25d0 * deltauv
#endif
      end

************************************************************************
      function B001_coli(p2,m12,m22)
************************************************************************
*     scalar one loop integral B301 = B001_coli                        *
*     p2 outer momentum                                                *
*     m12,m22 inner masses squared                                     *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  21.10.08 ad            *
************************************************************************
      implicit none
      complex*16 A0_coli,B1_coli,B11_coli,B001_coli
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      logical    finarg
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
  
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' B001_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' B001_coli: p2  = ',p2
          write(*,111) ' B001_coli: m12 = ',m12
          write(*,111) ' B001_coli: m22 = ',m22
          stop
        endif
      endif
#endif

      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)
      if (finarg) then
        if (abs(q2)/(abs(q2)+abs(mm12+mm22)).gt.calacc) then
          B001_coli = (q2+mm12-mm22)*B11_coli(q2,mm12,mm22) 
     &        + 2d0*mm12*B1_coli(q2,mm12,mm22)
     &        - A0_coli(mm22) - (2d0*mm12+4d0*mm22-q2)/6d0
          B001_coli = B001_coli/8d0
        else
          B001_coli = (6d0*(mm12*B1_coli(dcmplx(0d0),mm12,mm22) 
     &        - A0_coli(mm22)) 
     &        - mm12 - 2d0*mm22) /36d0
        end if
      else
        B001_coli = 0d0
      end if
      end

************************************************************************
      function Bn_coli(n,p2,m12,m22)
************************************************************************
*     tensor one loop integral cB1..1 of rank n                        *
*     p2 outer momentum, m12,m22 inner masses                          *
*----------------------------------------------------------------------*
*     21.10.04 Ansgar Denner      last changed 24.10.08 ad             *
************************************************************************
      implicit   none
      integer    n,np
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      complex*16 Bn_coli,fpve_coli,yfpve_coli,xlogxe_coli,cln_coli,
     &    eta2_coli
      complex*16 x1,x2,y1,y2,r,sum1,sum2,sum1a,sum2a
      logical    finarg
      integer    i,tout
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
 
      tout = 20

c      write(*,*) 'Bn_coli in',n,p2,m12,m22
      
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' Bn_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' Bn_coli: p2  = ',p2
          write(*,111) ' Bn_coli: m12 = ',m12
          write(*,111) ' Bn_coli: m22 = ',m22
          stop
        endif
      endif
#endif

      np = n+1
      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)

c      write(*,*) 'Bn_coli small',q2,mm12,mm22,finarg

      if (finarg) then

c      write(*,*) 'Bn_coli finarg'

c---->  general case
        if (abs(q2).gt.calacc*abs(mm12+mm22)) then

          call crootse_coli(q2,mm12,mm22,x1,x2,y1,y2,r)

          if (abs(y1).gt.0.5d0.and.abs(y2).gt.0.5d0) then


            Bn_coli = (-1)**n*(-log(mm22/muuv2) 
     &           - fpve_coli(np,x1,y1,1d0) - fpve_coli(np,x2,y2,-1d0))
     &         /np

          elseif (abs(x1).lt.10d0.and.abs(x2).lt.10d0) then
            if(n.eq.0) then
              Bn_coli = 2d0 - cln_coli(q2/muuv2,-1d0)
     &            + xlogxe_coli(-x1,-1d0) + xlogxe_coli(-x2,1d0) 
     &            - xlogxe_coli(y1,-1d0) - xlogxe_coli(y2,1d0)
            else
              sum1  = 1d0
              sum2  = 1d0
              sum1a = 1d0
              sum2a = 1d0
              do i=2,n+1
                sum1 = sum1*x1+1d0 
                sum2 = sum2*x2+1d0
                sum1a = sum1a*x1 + 1d0/i 
                sum2a = sum2a*x2 + 1d0/i 
              end do
              Bn_coli = - cln_coli(q2/muuv2,-1d0)
     &            + x1**n*xlogxe_coli(-x1,-1d0)+ sum1a
     &            + x2**n*xlogxe_coli(-x2,1d0) + sum2a
     &            - sum1*xlogxe_coli(y1,-1d0) 
     &            - sum2*xlogxe_coli(y2,1d0) 
            end if
            if (mm12.ne.0d0) Bn_coli = Bn_coli 
     &          - eta2_coli(-x1,-x2,mm12/q2,-1d0,1d0,-q2)
            Bn_coli = (-1)**n*Bn_coli/np
c            write(tout,*) 'Bn_coli reg2'
          elseif (abs(x1).gt.0.5d0.and.abs(x2).gt.0.5d0) then
            if(n.eq.0) then
              Bn_coli = (-1)**np*( log(mm12/muuv2)
     &           + fpve_coli(1,y1,x1,-1d0) + fpve_coli(1,y2,x2,1d0) )/np
            else
              if (abs(x1).gt.10) then
                sum2  = 1d0
                sum2a = .5d0
                do i=2,n
                  sum2 = sum2*x2+1d0
                  sum2a = sum2a*x2 + dfloat(i)/dfloat(i+1) 
                end do
                sum2 = sum2*x2+1d0
                Bn_coli = (-1)**np*( log(mm12/muuv2)
     &              + fpve_coli(np,x1,y1,1d0) - fpve_coli(0,x1,y1,1d0)
     &              + sum2*fpve_coli(1,y2,x2,1d0) + sum2a )/np
              else
                sum1  = 1d0
                sum1a = .5d0
                do i=2,n
                  sum1 = sum1*x1+1d0 
                  sum1a = sum1a*x1 + dfloat(i)/dfloat(i+1) 
                end do
                sum1 = sum1*x1+1d0
                Bn_coli = (-1)**np*( log(mm12/muuv2)
     &              + sum1*fpve_coli(1,y1,x1,-1d0) + sum1a 
     &              + fpve_coli(np,x2,y2,-1d0) 
     &              - fpve_coli(0,x2,y2,-1d0))/np
              end if
            end if
c            write(20,*) 'Bn_coli reg3'
          else
            call errB0_coli(q2,mm12,mm22,Bn_coli,'Bn_coli ',12)
          endif
c---->  zero momentum
        elseif (abs(mm12-mm22).gt.calacc*abs(mm12+mm22)) then
          x2 = mm12/(mm12-mm22)
          y2 = mm22/(mm22-mm12)
          if(abs(y2).gt.0.5d0) then
            Bn_coli = (-1)**np*
     &          (log(mm22/muuv2) + fpve_coli(np,x2,y2,mm22-mm12))/np
c           write(*,*) 'Bn_coli sing1'
          else
            sum2  = 1d0
            sum2a = 1d0
            do i=2,n+1
              sum2 = sum2*x2+1d0
              sum2a = sum2a*x2 + 1d0/i 
            end do
            Bn_coli = (-1)**np*(log(mm12/muuv2) 
     &          - sum2*yfpve_coli(0,x2,y2,mm22-mm12) - sum2a)/np
c           write(*,*) 'Bn_coli sing2'
          end if
        else if (mm22.ne.0d0) then
          Bn_coli = (-1)**np*log(mm22/muuv2)/np
c         write(*,*) 'Bn_coli sing3'
        else
          call errB0_coli(q2,mm12,mm22,Bn_coli,'Bn_coli ',20)
        endif
      else
        if (q2.eq.mm12.and.mm22.eq.0d0.and.q2.ne.0d0) then
          sum1 = 1d0/np
          do i=1,n
             sum1 = sum1 + 1d0/i
          end do
          Bn_coli = (-1)**np*(log(mm12*coliminfscale2/muuv2) 
     &        - 2d0*sum1)/np
        else if (q2.eq.mm22.and.mm12.eq.0d0.and.q2.ne.0d0) then
          Bn_coli = (-1)**np*(log(mm22*coliminfscale2/muuv2) 
     &        - 2d0/np)/np
        else if (q2.eq.0d0.and.mm12.eq.mm22.and.mm12.ne.0d0) then
          Bn_coli = (-1)**np*(log(mm12*coliminfscale2/muuv2))/np
        else if (q2.eq.0d0.and.mm12.eq.0d0.and.mm22.ne.0d0) then
          Bn_coli = (-1)**np*(log(mm22*coliminfscale2/muuv2) 
     &        - 1d0/np)/np
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.ne.0d0) then
          sum1 = 1d0/np
          do i=1,n
             sum1 = sum1 + 1d0/i
          end do
          Bn_coli = (-1)**np*(log(mm12*coliminfscale2/muuv2) - sum1)/np
        else if (q2.ne.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          sum1 = 1d0/np
          do i=1,n
             sum1 = sum1 + 1d0/i
          end do
          Bn_coli = (-1)**np*(cln_coli(-q2*coliminfscale2/muuv2,-1d0) 
     &               - 1d0*sum1 - 1d0/np )/np
        else if (q2.eq.0d0.and.mm22.eq.0d0.and.mm12.eq.0d0) then
          Bn_coli = (-1)**np*log(muir2/muuv2)/np
#ifdef SING
          Bn_coli = Bn_coli + (-1)**np*delta1ir/np
#endif
        else
          call errB0_coli(q2,mm12,mm22,Bn_coli,'Bn_coli  ',28)
          write(*,*) q2-mm22,mm12
        end if
      end if   

c      write(*,*) 'Bn_coli a',Bn_coli,deltauv,muuv2

#ifdef SING
      Bn_coli = Bn_coli - (-1)**np*deltauv/np
#endif

c      write(*,*) 'Bn_coli =',Bn_coli

      end

************************************************************************
      function B00n_coli(n,p2,m12,m22)
************************************************************************
*     tensor one loop integral cB001..1 of rank n+2                    *
*     p2 outer momentum, m12,m22 inner masses                          *
*----------------------------------------------------------------------*
*     22.10.04 Ansgar Denner      last changed  21.10.08 ad            *
************************************************************************
      implicit none
      integer    n
      complex*16 A0_coli,Bn_coli,B00n_coli
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      logical    finarg
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
 
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' B00n_coli: imaginary momentum arguments not ',
     &             'implemented'
          write(*,111) ' B00n_coli: p2  = ',p2
          write(*,111) ' B00n_coli: m12 = ',m12
          write(*,111) ' B00n_coli: m22 = ',m22
          stop
        endif
      endif
#endif

      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)
      if (finarg) then
        if (abs(q2)/(abs(q2)+abs(mm12+mm22)).gt.calacc) then
          B00n_coli = (q2+mm12-mm22)*Bn_coli(n+1,q2,mm12,mm22) 
     &        + 2d0*mm12*Bn_coli(n,q2,mm12,mm22)
     &        +(-1)**n *( A0_coli(mm22) 
     &             +2d0*(mm22 + mm12/(n+1) - q2/(n+3))/(n+2) )
          B00n_coli = B00n_coli/(6d0+2d0*n)
        else
          B00n_coli = (mm12*Bn_coli(n,dcmplx(0d0),mm12,mm22)
     &        +(-1)**n *( A0_coli(mm22) 
     &             + (mm22 + mm12/(n+1))/(n+2)))
     &          /(4d0+2d0*n)
        end if
      else
        B00n_coli = 0d0
      end if
      end

************************************************************************
      function B0000n_coli(n,p2,m12,m22)
************************************************************************
*     tensor one loop integral cB00001..1 of rank n+2                  *
*     p2 outer momentum, m12,m22 inner masses                          *
*----------------------------------------------------------------------*
*     22.10.04 Ansgar Denner      last changed  21.10.08 ad            *
************************************************************************
      implicit none
      integer    n
      complex*16 A0_coli,B00n_coli,B0000n_coli
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      logical    finarg,finarg2
      data finarg2 /.true./
      save       finarg2
 
#include "params_coli.h"
#include "common_coli.h"
#include "checkparams_coli.h"
  
#ifdef CHECK 
 100  format(((a)))
 101  format(a22,g25.17)
 111  format(a22,2('(',g24.17,',',g24.17,') ':))
      if (argcheck) then
        if(imag(p2).ne.0d0) then
          write(*,100) ' B0000n_coli: imaginary momentum arguments',
     &             ' not implemented'
          write(*,111) ' B0000n_coli: p2  = ',p2
          write(*,111) ' B0000n_coli: m12 = ',m12
          write(*,111) ' B0000n_coli: m22 = ',m22
          stop
        endif
      endif
#endif

      if (finarg2.and.n.gt.2) then
        write(*,*) 'B0000n_coli : function has not been tested!!!!!'
        write(*,*) 'B0000n_coli : n = ',n     
        finarg2 = .true.
      end if

      call elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)
      if (finarg) then
        if (abs(q2)/(abs(q2)+abs(mm12+mm22)).gt.calacc) then
          B0000n_coli = (q2+mm12-mm22)*B00n_coli(n+1,q2,mm12,mm22) 
     &        + 2d0*mm12*B00n_coli(n,q2,mm12,mm22)
     &        +(-1)**n *( mm22*(A0_coli(mm22)+mm22/2d0)/4d0 
     &             + (mm22*mm22  
     &               + 2d0*mm12*(mm22+mm12/(n+1))/(n+2)
     &               + 2d0*q2/(n+4)*(q2/(n+5)-mm22-2d0*mm12/(n+2)))
     &              /(2*(n+3)) )
          B0000n_coli = B0000n_coli/(10d0+2d0*n)
        else
          B0000n_coli = (mm12*B00n_coli(n,dcmplx(0d0),mm12,mm22)
     &        +(-1)**n *( mm22*(A0_coli(mm22)+mm22/2d0)/4d0
     &             + (mm22*mm22  
     &               + 2d0*mm12*(mm22+mm12/(n+1))/(n+2)) /(4*(n+3)) ) )
     &          /(6d0+2d0*n)
c>          write(*,*) 'cB0000 ',(6d0+2d0*n)*B0000n_coli
c>     &          - (-1)**n *( 
c>     &             + (mm22*mm22  
c>     &               + 2d0*mm12*(mm22+mm12/(n+1))/(n+2)) /(4*(n+3)) )
c>          write(*,*) 'cB0000 ',mm12*B00n_coli(n,dcmplx(0d0),mm12,mm22)
c>     &        +(-1)**n * mm22*(A0_coli(mm22)+2d0*mm22/(n+4))/(n+4)
c>          write(*,*) 'cB0000 ',mm12*B00n_coli(n,dcmplx(0d0),mm12,mm22),
c>     &        +(-1)**n * mm22*A0_coli(mm22)/(n+4),
c>     &        +(-1)**n * mm22*(2d0*mm22/(n+4))/(n+4)
        end if
      else
        B0000n_coli = 0d0
      end if
      end

************************************************************************
      subroutine elminf2iv_coli(p2,m12,m22,q2,mm12,mm22,finarg)
************************************************************************
*     eliminate small masses squared in vertex   p2,m12,m22            *
*----------------------------------------------------------------------*
*     16.03.04 Ansgar Denner                                           *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 q2
      complex*16 mm12,mm22
      complex*16 p2s
      complex*16 mm12s,mm22s
      complex*16 elimminf2_coli
      logical    finarg
 
      p2s  = elimminf2_coli(p2)
      mm12s = elimminf2_coli(m12)
      mm22s = elimminf2_coli(m22)

      if(p2s.ne.0d0.or.mm12s.ne.0d0.or.mm22s.ne.0d0) then
        q2 = p2s
        mm12 = mm12s
        mm22 = mm22s
        finarg = .true.
      else
        q2 = p2
        mm12 = m12
        mm22 = m22
        finarg  = .false.
      end if
      end

************************************************************************
      subroutine errB0_coli(q2,m12,m22,value,name,code)
************************************************************************
*     error messages for B0 functions                                  *
*     q2 outer momentum      m12,m22 inner masses squared              *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  21.10.08 ad            *
***************************:********************************************
      implicit none
      complex*16 value
      complex*16 q2
      complex*16 m12,m22
      character  name*6
      integer    code
 
      if (code.eq.16) then
        write(*,*) name,' not well defined'
      elseif (code.eq.12) then
        write(*,*) name,' not defined'
      elseif (code.eq.4) then
        write(*,*) name,' possibly not well defined'
      elseif (code.eq.20) then
        write(*,*) name,' not defined for zero arguments'
      elseif (code.eq.24) then
        write(*,*) name,' not defined for almost zero masses'
      elseif (code.eq.28) then
        write(*,*) name,' not defined for this case of small masses'
      end if
      write(*,*) ' q2,m12,m22 = ',q2,m12,m22
      write(*,*) '  ',name,'    = ',value
      stop
      end


************************************************************************
      function fpve_coli(n,x,y,eps)
************************************************************************
*     passarino and veltmans f(n,x) function                           *
*     y = 1-x,  x  = x + i eps                                         *
*----------------------------------------------------------------------*
*     21.07.93 Ansgar Denner      last changed  16.03.04 ad            *
************************************************************************
      implicit   none
      complex*16 x,y,fpve_coli,cln_coli
      real*8     calacc,eps
      integer    n,m
      common /calacc/ calacc
 
      if(abs(x).lt.10d0) then
        if (n.eq.0) then
          fpve_coli = -cln_coli(-y/x,eps)
        elseif (x.ne.dcmplx(0d0)) then
          fpve_coli = -x**n*cln_coli(-y/x,eps)
          do m=1,n
            fpve_coli = fpve_coli - x**(n-m)/m
          enddo
        else
          fpve_coli = -1d0/n
        endif
      else
        fpve_coli = 0d0
        do m=1,10
          fpve_coli = fpve_coli + 1d0/((m+n)*x**m)
          if (abs(1d0/(x**m*fpve_coli)).lt.calacc**2) return
        enddo
      endif
      end

************************************************************************
      function yfpve_coli(n,x,y,eps)
************************************************************************
*     passarino and veltmans f(n,x) function times (1-x) = y           *
*     x  = x + i eps                                                   *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  16.03.04 ad            *
************************************************************************
      implicit   none
      complex*16 x,y,fpve_coli,yfpve_coli
      real*8     eps 
      integer    n
 
      if(abs(y).eq.0d0) then
        yfpve_coli = 0d0
      else
        yfpve_coli = y*fpve_coli(n,x,y,eps)
      endif
      end

************************************************************************
      function xyfpve_coli(n,x,y,eps)
************************************************************************
*     passarino and veltmans f(n,x) function times x*(1-x) = x*y       *
*     x  = x + i eps                                                   *
*----------------------------------------------------------------------*
*     20.09.04 Ansgar Denner      last changed  20.09.04 ad            *
************************************************************************
      implicit   none
      complex*16 x,y,fpve_coli,xyfpve_coli
      real*8     eps 
      integer    n
 
      if(abs(y).eq.0d0.or.abs(x).eq.0d0) then
        xyfpve_coli = 0d0
      else
        xyfpve_coli = x*y*fpve_coli(n,x,y,eps)
      endif
      end

************************************************************************
      function xlogxe_coli(x,eps)
************************************************************************
*     x*log(x),    x  = x + i eps                                      *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  16.03.04 ad            *
************************************************************************
      implicit   none
      complex*16 x,xlogxe_coli,cln_coli
      real*8     eps
 
      if(abs(x).eq.0d0) then
        xlogxe_coli = 0d0
      else
        xlogxe_coli = x*cln_coli(x,eps)
      endif
      end

