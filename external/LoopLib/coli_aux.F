#include "global_coli.h"

#define WIB

************************************************************************
*                                                                      *
*     auxiliary functions for COLI library                             *
*                                                                      *
************************************************************************
*                                                                      *
*     last changed  28.10.09  Ansgar Denner                            *
*                                                                      *
************************************************************************
* subroutines:                                                         *
* initcoli                                                             *
* setmuuv2_coli, setmuir2_coli, setdeltauv_coli, setdeltair_coli       *
* setshiftms2_coli                                                     *
* setminf2_coli, clearcoliminf2, setminfscale2_coli                    *
* functions:                                                           *
* cln_coli, cspenc_coli, cspenh_coli                                   *  
* csp2con_coli, csp2cos_coli, csp3con_coli, csp4con_coli               *
* eta2_coli, eta2s_coli                                                *
* minfscaledown_coli, minfscaledown2_coli                              *
* elimminf_coli, elimminf2_coli                                        *
*
*
*
* ettile,olieta,olietal,olietae,olietae2,olietae3,olietae2l,olietae3l  *
* fpve,  yfpve, xyfpve,xlogxe                                          *
* fpve,  yfpve, xyfpve,xlogxe                                          *

************************************************************************
      subroutine initcoli
************************************************************************
* fixing of default values for various variables                       *
*----------------------------------------------------------------------*
* 28.08.08  Ansgar Denner         last changed  28.08.08               *
************************************************************************
      implicit   none
      real*8     eps,calacc,impacc
      integer    sym
 
      common /impacc/ impacc
      common /calacc/ calacc
      common /sym/    sym
  
      integer  errout,testout
      common /out/ errout,testout

      logical    init
      data init /.false./
      save init

      if (init) return
      init = .true.
#ifdef PRINT
      print *,'======================================================='
      print *,' COLI 1.0beta: a package to evaluate one-loop integrals'
      print *,'written by A. Denner, PSI, Villigen, Switzerland'
      print *,'======================================================='
      print *,'THIS CODE IS PRELIMINARY!!!'
      print *,'There is no guarantee that the results are correct.'
c     print *,'This code has been used and tested for e+e- -> 4f.' 
c     print *,'Cases that were not required there' 
c     print *,'are either not implemented or not tested.'
      print *,'======================================================='
      print *,'If you use the code, you must refer to it'
      print *,'and cite the papers it is based on.'
      print *,'A. Denner and S. Dittmaier, Nucl.Phys.B734:62-115,2006' 
      print *,'A. Denner and S. Dittmaier, Nucl.Phys.B658:175-202,2003'
      print *,'  for the tensor integral reduction,'
      print *,'A. Denner, U. Nierste, and R. Scharf,'
      print *,'Nucl.Phys.B367:637-656,1991'
      print *,'A. Denner and W. Beenakker, Nucl.Phys.B338:349-370,1990'
      print *,'  for the scalar integrals.'
      print *,'========================================================'
#endif

c print information messages
C      call setinfo_coli
      call unsetinfo_coli
c mu^2 of dimensional regularization in uv-divergent loop integrals
      call setmuuv2_coli(1d0)
c photon/gluon mass in ir-divergent loop integrals
c mu^2 of dimensional regularization in ir-divergent loop integrals
      call setmuir2_coli(1d0)
c scale factor for mass-singular masses
      call setminfscale2_coli(1d0)
#ifdef SING
c uv-divergence 2/(4-D)
      call setdeltauv_coli(0d0)
c ir-divergence 2/(4-D)
      call setdeltair_coli(0d0,0d0)
c extra term for mass-singular squared logarithms
      call setshiftms2_coli(0d0)
#endif
c infinitesimal parameter
      eps    = 1d-20
c size of imaginary parts below which explicit ieps take over
      impacc = 1d-16
c should be set to the accuracy of the computer: n digits -> 1d-n
      calacc = 5d-16

c sym = 1 calculation of symmetric tensor coefficients
c sym = 0 no calculation of symmetric tensor coefficients
c     sym = 0

c errout:  output channel for error messages
      errout = 6
c testout: output channel for testing purposes
      testout = 9

      end

************************************************************************
      function cln_coli(cz,eps)
************************************************************************
*     complex logarithm of cz + i*eps                                  *
*----------------------------------------------------------------------*
*     09.01.90 Ansgar Denner         last changed 15.03.04             *
************************************************************************
      implicit   none
      real*8     eps
      complex*16 cln_coli,cz

#include "params_coli.h"
#include "checkparams_coli.h"

#ifdef CHECK
      if(argcheck) then  
         if(abs(dimag(cz)).lt.impacc*abs(dreal(cz)).and.
     &        dreal(cz).le.0d0.and.eps.eq.0d0.or.
     &        cdabs(cz).eq.0d0) then 
            write(*,*) 'cln_coli:  argument on cut '
            write(*,*) 'cln_coli:  eps = ',eps
            write(*,*) 'cln_coli:  cz  = ',cz
c            stop
         endif
      endif
#endif

      if(abs(dimag(cz)).gt.impacc*abs(dreal(cz)).or.dreal(cz).gt.0d0)
     &    then
        cln_coli=cdlog(cz)
      elseif(eps.ne.0) then
        cln_coli=cdlog(-cz)+dcmplx(0d0,pi)*dsign(1d0,eps)
      elseif(dimag(cz).ne.0d0) then
        cln_coli=cdlog(cz)
        write(*,*) 'cln_coli:  argument effectively on cut '
        write(*,*) 'cln_coli:  eps = ',eps
        write(*,*) 'cln_coli:  cz  = ',cz
      else 
        cln_coli=cdlog(-cz)+dcmplx(0d0,pi)
        write(*,*) 'cln_coli:  argument on cut '
        write(*,*) 'cln_coli:  eps = ',eps
        write(*,*) 'cln_coli:  cz  = ',cz
        stop
      endif

      end
************************************************************************
      function cspenc_coli(cz,eps)
************************************************************************
*       complex spence function  of cz + i*eps                         *
*       calculated by mapping on the area where there is a quickly     *
*       convergent series                                              *
*----------------------------------------------------------------------*
*     08.01.90 Ansgar Denner        last changed  3.02.97 ad           *
************************************************************************
      implicit   none
      real*8     eps
      complex*16 cspenc_coli,cz
      real*8     az,rz,az1
      complex*16 cz1,cspenh_coli,cln_coli

#include "params_coli.h"
#include "checkparams_coli.h"

#ifdef CHECK
      if(argcheck) then  
         if(dreal(cz).gt.1d0.and.
     &        abs(dimag(cz)).lt.impacc*abs(dreal(cz)).and.eps.eq.0)then
            write(*,*) 'cspen_coli:  argument on cut '
            write(*,*) 'cspen_coli:  eps = ',eps
            write(*,*) 'cspen_coli:  cz  = ',cz
            stop
         endif
      endif
#endif

      cz1     = 1d0-cz
      az1     = abs(cz1)
      az      = abs(cz)
      rz      = dreal(cz)

c      write(*,*) 'cspenc ',cz,az,az1,az1-1d0,az1.lt.1d0

      if (az1.lt.calacc) then
         cspenc_coli = pi2_6
      else if (rz.lt.0.5d0) then
         if (az.lt.1d0) then
            cspenc_coli =  cspenh_coli(cz,eps)
         else
            cspenc_coli = -pi2_6 - .5d0*cln_coli(-cz,-eps)**2
     &                -  cspenh_coli(1d0/cz,-eps)
         end if
      else
         if (az1.lt.1d0) then
            cspenc_coli =  pi2_6 - cln_coli(cz,eps)*cln_coli(cz1,-eps)
     &                       -  cspenh_coli(cz1,-eps)
         else
            cspenc_coli = 2d0*pi2_6 + .5d0*cln_coli(-cz1,-eps)**2
     &              - cln_coli(cz,eps)*cln_coli(cz1,-eps)
     &              +  cspenh_coli(1d0/cz1,eps)
         end if
      end if
      end

************************************************************************
      function cspenh_coli(cz,eps)
************************************************************************
*       complex spence function of cz + i*eps                          *
*       in convergence region                                          *
*       calculation of bernoulli series                                *
*----------------------------------------------------------------------*
*     09.01.90 Ansgar Denner        last changed 02.09.08              *
************************************************************************
      implicit   none
      complex*16 cspenh_coli,cln_coli,cz,x,x2
      real*8     eps
      real*8 b(11)/
     1   0.1666666666666666666666666667d0,
     2  -0.0333333333333333333333333333d0,
     3   0.0238095238095238095238095238d0,
     4  -0.0333333333333333333333333333d0,
     5   0.0757575757575757575757575758d0,
     6  -0.2531135531135531135531135531d0,
     7   1.1666666666666666666666666667d0,
     8  -7.0921568627450980392156862745d0,
     9  54.97117794486215538847117794486d0,
     +  -529.124242424242424242424242424242d0,
     1  6192.123188405797101449275362318d0  /
c     note:                 b(n)=b2n
c     b(1)=1./6.
c     b(2)=-1./30.
c     b(3)=1./42.
c     b(4)=-1./30.
c     b(5)=5./66.
c     b(6)=-691./2730.
c     b(7)=7./6.
c     b(8)=-3617./510.
c     b(9)=43867./798.
c     b(10)=-174611./330.
c     b(11)=854513./138.
c     pi=3.1415926535897932384
c     pi*pi/6.=1.6449..., pi*pi/3=3.28986...
c
      integer    j
      real*8     factor
      complex*16 power,term,csp

#include "params_coli.h"
#include "checkparams_coli.h"

#ifdef CHECK
      if(argcheck) then  
         if(abs(cz).gt.1d0+calacc.or.dreal(cz).gt.0.5d0)then
            write(*,*)'cspenh_coli:  argument not in convergence region'
            write(*,*) 'cspenh_coli:  eps = ',eps
            write(*,*) 'cspenh_coli:  cz  = ',cz
            write(*,*) 'cspenh_coli:  |cz|= ',abs(cz),abs(cz)-1d0,
     &          abs(cz).gt.1d0
            write(*,*) 'cspenh_coli: Re(z)= ',dreal(cz),dreal(cz)-.5d0,
     &          dreal(cz).gt.0.5d0
            stop
         endif
      endif
#endif

      b(11)  =    854513d0/ 138d0
      b(10)  =  - 174611d0/ 330d0
      b(9)   =     43867d0/ 798d0
      b(8)   =  -   3617d0/ 510d0
      b(7)   =         7d0/   6d0
      b(6)   =  -    691d0/2730d0
      b(5)   =         5d0/  66d0
      b(4)   =  -      1d0/  30d0
      b(3)   =         1d0/  42d0
      b(2)   =  -      1d0/  30d0
      b(1)   =         1d0/   6d0
      x      =  -cln_coli(1d0-cz,-eps)
      x2     =  x*x
      power  =  x
      factor =  1d0
      cspenh_coli =  x - x2/4d0
      do 10 j=2,22,2
         factor = factor / j / (j+1)
         power  = power * x2
         term   = b(j/2) * factor * power
         csp    = cspenh_coli + term
         if (csp.eq.cspenh_coli) return
         cspenh_coli = csp
10    continue
      if (cdabs(term/csp).gt.calacc) then
        write(80,*) 'cspenh_coli converges badly  ',cz,x
        write(80,*) 'cspenh_coli converges badly  ',csp-term,cspenh_coli
     &      ,term
      end if 
      end


************************************************************************
      function cspcon_coli(z1,z2,pz,ompz,i1,i2)
************************************************************************
*  complex spence function   sp(1-z1*z2)                               *
*                          plus continuation terms                     *
*  pz = z1*z2, ompz = 1-pz                                             * 
*  i1,i2 determine signs of infinitesimal imaginary parts of           *
*  z1,z2                                                               *
*----------------------------------------------------------------------*
*  16.10.08 Ansgar Denner      last changed 12.12.08                   *
************************************************************************
      implicit   none
      complex*16 z1,z2,pz,ompz
      real*8     i1,i2,ip
      complex*16 cspcon_coli,cspenc_coli,cln_coli,sumln
#ifdef CHECK
      complex*16 cspcon_check
#endif


#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) 'cspcon_coli in: ',z1,z2,pz,ompz,i1,i2,ip,pz.ne.1d0

#ifdef CHECK
      if(argcheck) then  
         if(abs(z1*z2-pz).gt.1d3*calacc*abs(pz).or.
     &        abs(1d0-pz-ompz).gt.
     &      1d3*calacc*max(cdabs(pz),cdabs(1d0-pz))) then
            write(*,*) 'cspcon_coli:  arguments inconsistent  '
            write(*,100) 'cspcon_coli: z1 = ',z1
            write(*,100) 'cspcon_coli: z2 = ',z2
            write(*,100) 'cspcon_coli: zp = ',pz,z1*z2
            write(*,100) 'cspcon_coli: omzp = ',ompz,1d0-pz
            write(*,100) 'cspcon_coli: test1  ',
     &          abs(z1*z2-pz),1d3*calacc*abs(pz)
            write(*,100) 'cspcon_coli: test2  ',
     &          abs(1d0-pz-ompz),
     &          1d3*calacc*max(cdabs(pz),cdabs(1d0-pz))
 100      format(a21,8g25.17)
            stop
         endif
      endif
#endif

c to ensure correct continuation for infinitesimal imaginary parts
      ip=-sign(1d0,i1)
      if(pz.eq.0d0)then
        cspcon_coli = pi2_6
      else
        sumln=cln_coli(z1,i1)+cln_coli(z2,i2)
        if(abs(pz).lt.1d0) then
          if(dreal(pz).gt.0.5d0) then
            cspcon_coli = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &          *(cln_coli(pz,ip)-sumln)
          else
            cspcon_coli = pi2_6-cspenc_coli(pz,ip)
     &          -cln_coli(1d0-pz,-ip)*sumln 
          endif
        elseif(pz.ne.1d0) then
          if(dreal(pz).lt.2d0) then
            cspcon_coli = -cspenc_coli(1d0-1d0/pz,ip)
     &          +cln_coli(1d0-1d0/pz,ip)
     &          *(cln_coli(pz,ip)-sumln)
     &          -.5d0*(sumln)**2
c>            cspcon_coli = -cspenc_coli(1d0-1d0/pz,ip)
c>     &          +cln_coli(1d0-pz,-ip)
c>     &          *(cln_coli(pz,ip)-sumln)
c>     &        -.5d0*(cln_coli(pz,ip))**2
          else
            cspcon_coli = -pi2_6+cspenc_coli(1d0/pz,-ip)
     &          -cln_coli(1d0-1d0/pz,ip)*sumln
     &          -.5d0*(sumln)**2
c>            cspcon_coli = 2d0*pi2_6+cspenc_coli(1d0/pz,-ip)
c>     &          -cln_coli(1d0-pz,-ip)*sumln
c>     &        +.5d0*(cln_coli(-pz,-ip))**2
          endif
        else
          cspcon_coli =0d0
          if(abs(sumln).gt.calacc) then 
            write(*,*) 
     &          'cspcon_coli: continuation possibly wrong'
            write(*,*) 'cspcon_coli: z1 = ',z1
            write(*,*) 'cspcon_coli: z2 = ',z2
            write(*,*) 'cspcon_coli: pz = ',pz
            write(*,*) 'cspcon_coli:ompz= ',ompz
            write(*,*) 'cspcon_coli: i1 = ',i1
            write(*,*) 'cspcon_coli: i2 = ',i2
            write(*,*) 'cspcon_coli: sumln = ',sumln
          endif
        endif
      endif

#ifdef CHECK
      if(conscheck) then
        if(pz.ne.1d0.and.pz.ne.0d0)then
          cspcon_check = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &        *(cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2))
        elseif(pz.eq.1d0)then
          cspcon_check = 0d0
        else
          cspcon_check = pi2_6
        endif 

c          write(*,*) 'cspcon_coli out:',cspcon_check

        if(cdabs(cspcon_coli-cspcon_check).gt.
     &      testacc*cdabs(cspcon_coli).and.
     &      cdabs(cspcon_coli).gt.1d3*testacc) then
          write(testout,*) 'cspcon_coli inconsistent'
          write(testout,100) 'cspcon_coli z1  = ',z1,i1
          write(testout,100) 'cspcon_coli z2  = ',z2,i2
          write(testout,100) 'cspcon_coli pz  = ',pz,ip
          write(testout,100) 'cspcon_coli ompz= ',ompz
          write(testout,100) 'cspcon_coli   = ',cspcon_coli
          write(testout,100) 'cspcon_check  = ',cspcon_check
          write(testout,100) 'diff        = ',
     &        cdabs(cspcon_coli/cspcon_check-1d0)
          write(*,*) 'cspcon_coli inconsistent'
          write(*,100) 'cspcon_coli z1  = ',z1,i1
          write(*,100) 'cspcon_coli z2  = ',z2,i2
          write(*,100) 'cspcon_coli pz  = ',pz,ip
          write(*,100) 'cspcon_coli ompz= ',ompz
          write(*,100) 'cspcon_coli input = ',z1,z2,pz,ompz
          write(*,100) 'cspcon_coli inputi= ',i1,i2
          write(*,100) 'cspcon_coli   = ',cspcon_coli
          write(*,100) 'cspcon_check  = ',cspcon_check
          write(*,100) 'diff        = ',
     &        cdabs(cspcon_coli/cspcon_check-1d0)
          stop
        endif
      endif
#endif

c      write(*,*) 'cspcon_coli out:',cspcon_coli,cspcon_check

      end

************************************************************************
      function cspcos_coli(z1,z2,i1,i2)
************************************************************************
*  complex spence function   sp(1-z1*z2)                               *
*                          plus continuation terms                     *
*  i1,i2 determine signs of infinitesimal imaginary parts of           *
*  z1,z2                                                               *
*----------------------------------------------------------------------*
*  30.04.08 Ansgar Denner      last changed 28.10.08                   *
************************************************************************
      implicit   none
      complex*16 z1,z2,pz,ompz
      real*8     i1,i2,ip
      complex*16 cspcos_coli,cspenc_coli,cln_coli,sumln
#ifdef CHECK
      complex*16 cspcos_check
#endif


#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) 'cspcos_coli in: ',z1,z2,i1,i2

 100  format(a21,8g25.17)
      pz=z1*z2
      ompz=1d0-pz
c to ensure correct continuation for infinitesimal imaginary parts
      if(i2.ne.0d0)then
        ip = i2*sign(1d0,dreal(z1))
      else
        ip = i1*sign(1d0,dreal(z2))
      endif
      if(pz.eq.0d0)then
        cspcos_coli = pi2_6
      else
        sumln=cln_coli(z1,i1)+cln_coli(z2,i2)
        if(abs(pz).lt.1d0) then
          if(dreal(pz).gt.0.5d0) then
            cspcos_coli = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &          *(cln_coli(pz,ip)-sumln)
          else
            cspcos_coli = pi2_6-cspenc_coli(pz,ip)
     &          -cln_coli(1d0-pz,-ip)*sumln 
          endif
        elseif(pz.ne.1d0) then
          if(dreal(pz).lt.2d0) then
            cspcos_coli = -cspenc_coli(1d0-1d0/pz,ip)
     &          +cln_coli(1d0-1d0/pz,ip)
     &          *(cln_coli(pz,ip)-sumln)
     &          -.5d0*(sumln)**2
c>            cspcos_coli = -cspenc_coli(1d0-1d0/pz,ip)
c>     &          +cln_coli(1d0-pz,-ip)
c>     &          *(cln_coli(pz,ip)-sumln)
c>     &        -.5d0*(cln_coli(pz,ip))**2
          else
            cspcos_coli = -pi2_6+cspenc_coli(1d0/pz,-ip)
     &          -cln_coli(1d0-1d0/pz,ip)*sumln
     &          -.5d0*(sumln)**2
c>            cspcos_coli = 2d0*pi2_6+cspenc_coli(1d0/pz,-ip)
c>     &          -cln_coli(1d0-pz,-ip)*sumln
c>     &        +.5d0*(cln_coli(-pz,-ip))**2
          endif
        else
          cspcos_coli =0d0
          if(sumln.ne.cd0) then 
            write(*,*) 
     &          'cspcos_coli: continuation possibly wrong'
            write(*,*) 'cspcos_coli: z1 = ',z1
            write(*,*) 'cspcos_coli: z2 = ',z2
            write(*,*) 'cspcos_coli: i1 = ',i1
            write(*,*) 'cspcos_coli: i2 = ',i2
          endif
        endif
      endif

#ifdef CHECK
      if(conscheck) then
        if(ip.eq.0) ip=1  ! for z1,z2>0 needed to avoid error message
        if(pz.ne.1d0.and.pz.ne.0d0)then
          cspcos_check = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &        *(cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2))
        elseif(pz.eq.1d0)then
          cspcos_check = 0d0
        else
          cspcos_check = pi2_6
        endif 

c          write(*,*) 'cspcos_coli out:',cspcos_check

        if(abs(cspcos_coli-cspcos_check).gt.testacc*abs(cspcos_coli)
     &      .and.abs(cspcos_coli).gt.testacc) 
     &      then
          write(testout,*) 'cspcos_coli inconsistent'
          write(testout,100) 'cspcos_coli z1  = ',z1,i1
          write(testout,100) 'cspcos_coli z2  = ',z2,i2
          write(testout,100) 'cspcos_coli pz  = ',pz,ip
          write(testout,100) 'cspcos_coli ompz= ',ompz
          write(testout,100) 'cspcos_coli   = ',cspcos_coli
          write(testout,100) 'cspcos_check  = ',cspcos_check
          write(testout,100) 'diff        = ',
     &        cdabs(cspcos_coli/cspcos_check-1d0)
          write(*,*) 'cspcos_coli inconsistent'
          write(*,100) 'cspcos_coli z1  = ',z1,i1
          write(*,100) 'cspcos_coli z2  = ',z2,i2
          write(*,100) 'cspcos_coli pz  = ',pz,ip
          write(*,100) 'cspcos_coli ompz= ',ompz
          write(*,100) 'cspcos_coli input = ',z1,z2
          write(*,100) 'cspcos_coli inputi= ',i1,i2
          write(*,100) 'cspcos_coli   = ',cspcos_coli
          write(*,100) 'cspcos_check  = ',cspcos_check
          write(*,100) 'diff        = ',
     &        cdabs(cspcos_coli/cspcos_check-1d0)
          stop
        endif
      endif
#endif

c      write(*,*) 'cspcos_coli out:',cspcos_coli,cspcos_check

      end

#ifdef CHECK
************************************************************************
      function cspco0_coli(z1,z2,ompz,i1,i2)
************************************************************************
*  complex spence function   sp(1-z1*z2)                               *
*                          plus continuation terms                     *
*  for ompz = 1-pz -> 0                                                * 
*  i1,i2 determine signs of infinitesimal imaginary parts of z1, z2    *
*----------------------------------------------------------------------*
*  04.11.08 Ansgar Denner      last changed 04.11.08                   *
************************************************************************
      implicit   none
      complex*16 z1,z2,ompz
      real*8     i1,i2,ip
      complex*16 cspco0_coli,cspenc_coli,cln_coli,sumln
      complex*16 pz
#ifdef CHECK
      complex*16 cspco0_check
#endif


#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) 'cspco0_coli in: ',z1,z2,ompz,i1,i2

#ifdef CHECK
      if(argcheck) then  
         if(abs(z1*z2-1d0).gt.1d3*calacc.or.
     &        abs(ompz).gt.1d-3) then
            write(*,*) 'cspco0_coli:  arguments inconsistent  '
            write(*,100) 'cspco0_coli: z1 = ',z1
            write(*,100) 'cspco0_coli: z2 = ',z2
            write(*,100) 'cspco0_coli:omzp= ',ompz,1d0-z1*z2
            write(*,*) 'cspco0_coli: test1  ',
     &         abs(z1*z2-1d0).gt.1d3*calacc
            write(*,*) 'cspco0_coli: test2  ',
     &          abs(ompz).gt.1d-3
 100      format(a21,8g25.17)
c            stop
         endif
      endif
#endif

c to ensure correct continuation for infinitesimal imaginary parts
      ip=-sign(1d0,i1+i2)
      sumln=cln_coli(z1,i1)+cln_coli(z2,i2)
      if (sumln.ne.cd0)then
        cspco0_coli = -sumln*cln_coli(ompz,-ip)
c>        if(dimag(ompz).gt.0d0.and.ip.gt.0d0)then
c>          cspco0_coli = cspco0_coli - sumln*dcmplx(0d0,-2d0*pi)
c>        elseif(dimag(ompz).lt.0d0.and.ip.lt.0d0)then
c>          cspco0_coli = cspco0_coli - sumln*dcmplx(0d0,2d0*pi)
c>        endif
      else
        cspco0_coli=cd0
      endif

#ifdef CHECK
      if(conscheck) then
        pz=z1*z2 
        cspco0_check = cspenc_coli(1d0-pz,-ip)+cln_coli(ompz,-ip)
     &      *(cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2))
c>        if(dimag(ompz).gt.0d0.and.ip.gt.0d0)then
c>          cspco0_check = cspco0_check + dcmplx(0d0,-2d0*pi)
c>     &      *(cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2))
c>        elseif(dimag(ompz).lt.0d0.and.ip.lt.0d0)then
c>          cspco0_check = cspco0_check + dcmplx(0d0,2d0*pi)
c>     &      *(cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2))
c>        endif

c          write(*,*) 'cspco0_coli out:',cspco0_check

        if(cdabs(cspco0_coli-cspco0_check).gt.
     &      testacc*cdabs(cspco0_coli).and.
     &      cdabs(cspco0_coli).gt.1d3*testacc) then
          write(testout,*) 'cspco0_coli inconsistent'
          write(testout,100) 'cspco0_coli z1  = ',z1,i1
          write(testout,100) 'cspco0_coli z2  = ',z2,i2
          write(testout,100) 'cspco0_coli pz  = ',pz,ip
          write(testout,100) 'cspco0_coli ompz= ',ompz
          write(testout,100) 'cspco0_coli   = ',cspco0_coli
          write(testout,100) 'cspco0_check  = ',cspco0_check
          write(testout,100) 'diff        = ',
     &        cdabs(cspco0_coli/cspco0_check-1d0)
          write(*,*) 'cspco0_coli inconsistent'
          write(*,100) 'cspco0_coli z1  = ',z1,i1
          write(*,100) 'cspco0_coli z2  = ',z2,i2
          write(*,100) 'cspco0_coli pz  = ',pz,ip
          write(*,100) 'cspco0_coli ompz= ',ompz
          write(*,100) 'cspco0_coli input = ',z1,z2,pz,ompz
          write(*,100) 'cspco0_coli inputi= ',i1,i2
          write(*,100) 'cspco0_coli   = ',cspco0_coli
          write(*,100) 'cspco0_check  = ',cspco0_check
          write(*,100) 'diff        = ',
     &        cdabs(cspco0_coli/cspco0_check-1d0)
          stop
        endif
      endif
#endif

c      write(*,*) 'cspco0_coli out:',cspco0_coli,cspco0_check

      end
#endif

************************************************************************
      function csp2con_coli(z1,z2,pz,ompz,i1,i2)
************************************************************************
*  complex spence function   sp(1-z1*z2)                               *
*                          plus continuation terms                     *
*  pz = z1*z2, ompz = 1-pz                                             * 
*  i1,i2 determine signs of infinitesimal imaginary parts of           *
*  z1,z2                                                               *
*----------------------------------------------------------------------*
*  30.04.08 Ansgar Denner      last changed 03.06.08                   *
************************************************************************
      implicit   none
      complex*16 z1,z2,pz,ompz
      real*8     i1,i2,ip
      complex*16 csp2con_coli,cspenc_coli,cln_coli,sumln
#ifdef CHECK
      complex*16 csp2con_check
#endif


#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) 'csp2con_coli in: ',z1,z2,pz,ompz,i1,i2,ip,pz.ne.1d0

#ifdef CHECK
      if(argcheck) then  
         if(abs(z1*z2-pz).gt.1d2*calacc*abs(pz).or.
     &        abs(1d0-pz-ompz).gt.
     &      1d2*calacc*max(cdabs(pz),cdabs(1d0-pz))) then
            write(*,*) 'csp2con_coli:  arguments inconsistent  '
            write(*,100) 'csp2con_coli: z1 = ',z1
            write(*,100) 'csp2con_coli: z2 = ',z2
            write(*,100) 'csp2con_coli: zp = ',pz,z1*z2
            write(*,100) 'csp2con_coli: omzp = ',ompz,1d0-pz
 100      format(a21,4g25.17)
            stop
         endif
      endif
#endif

      ip=1d0
      if(pz.eq.0d0)then
        csp2con_coli = pi2_6
      else
        sumln=cln_coli(z1,i1)+cln_coli(z2,i2)
        if(abs(pz).lt.1d0) then
          if(dreal(pz).gt.0.5d0) then
            csp2con_coli = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &          *(cln_coli(pz,ip)-sumln)
          else
            csp2con_coli = pi2_6-cspenc_coli(pz,ip)
     &          -cln_coli(1d0-pz,-ip)*sumln 
          endif
        elseif(pz.ne.1d0) then
          if(dreal(pz).lt.2d0) then
#ifdef WIB
            csp2con_coli = -cspenc_coli(1d0-1d0/pz,ip)
     &          +cln_coli(1d0-1d0/pz,ip)
     &          *(cln_coli(pz,ip)-sumln)
     &          -.5d0*(sumln)**2
#else
            csp2con_coli = -cspenc_coli(1d0-1d0/pz,ip)
     &          +cln_coli(1d0-pz,-ip)
     &          *(cln_coli(pz,ip)-sumln)
     &        -.5d0*(cln_coli(pz,ip))**2
#endif
          else
#ifdef WIB
            csp2con_coli = -pi2_6+cspenc_coli(1d0/pz,-ip)
     &          -cln_coli(1d0-1d0/pz,ip)*sumln
     &          -.5d0*(sumln)**2
#else
            csp2con_coli = 2d0*pi2_6+cspenc_coli(1d0/pz,-ip)
     &          -cln_coli(1d0-pz,-ip)*sumln
     &        +.5d0*(cln_coli(-pz,-ip))**2
#endif
          endif
        else
          csp2con_coli =0d0
        endif
      endif

#ifdef CHECK
      if(conscheck) then
        ip=-1d0
        if(pz.ne.1d0.and.pz.ne.0d0)then
          csp2con_check = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &        *(cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2))
#ifdef WIB
          if(cdabs(pz).gt.1d0) then
            csp2con_check = csp2con_check - 
     &          (cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2))*
     &          (cln_coli(-pz,-ip)-.5d0*(cln_coli(pz,ip) 
     &          +cln_coli(z1,i1)+cln_coli(z2,i2)))
          endif 
#endif
        elseif(pz.eq.1d0)then
          csp2con_check = 0d0
        else
          csp2con_check = pi2_6
        endif 

c          write(*,*) 'csp2con_coli out:',csp2con_check

        if(cdabs(csp2con_coli-csp2con_check).gt.
     &      testacc*cdabs(csp2con_coli)) then
          write(testout,*) 'csp2con_coli inconsistent'
          write(testout,100) 'csp2con_coli z1  = ',z1,i1
          write(testout,100) 'csp2con_coli z2  = ',z2,i2
          write(testout,100) 'csp2con_coli pz  = ',pz,ip
          write(testout,100) 'csp2con_coli ompz= ',ompz
          write(testout,100) 'csp2con_coli   = ',csp2con_coli
          write(testout,100) 'csp2con_check  = ',csp2con_check
          write(testout,100) 'diff        = ',
     &        cdabs(csp2con_coli/csp2con_check-1d0)
          write(*,*) 'csp2con_coli inconsistent'
          write(*,100) 'csp2con_coli z1  = ',z1,i1
          write(*,100) 'csp2con_coli z2  = ',z2,i2
          write(*,100) 'csp2con_coli pz  = ',pz,ip
          write(*,100) 'csp2con_coli ompz= ',ompz
          write(*,100) 'csp2con_coli input = ',z1,z2,pz,ompz
          write(*,100) 'csp2con_coli inputi= ',i1,i2
          write(*,100) 'csp2con_coli   = ',csp2con_coli
          write(*,100) 'csp2con_check  = ',csp2con_check
          write(*,100) 'diff        = ',
     &        cdabs(csp2con_coli/csp2con_check-1d0)
          stop
        endif
      endif
#endif

c      write(*,*) 'csp2con_coli out:',csp2con_coli,csp2con_check

      end
************************************************************************
      function csp2cos_coli(z1,z2,i1,i2)
************************************************************************
*  complex spence function   sp(1-z1*z2)                               *
*                          plus continuation terms                     *
*  i1,i2 determine signs of infinitesimal imaginary parts of           *
*  z1,z2                                                               *
*----------------------------------------------------------------------*
*  30.04.08 Ansgar Denner      last changed 22.09.08                   *
************************************************************************
      implicit   none
      complex*16 z1,z2,pz,ompz
      real*8     i1,i2,ip
      complex*16 csp2cos_coli,cspenc_coli,cln_coli,sumln
#ifdef CHECK
      complex*16 csp2cos_check
#endif


#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) 'csp2cos_coli in: ',z1,z2,i1,i2

 100  format(a21,8g25.17)
      pz=z1*z2
      ompz=1d0-pz
      ip=1d0
      if(pz.eq.0d0)then
        csp2cos_coli = pi2_6
      else
        sumln=cln_coli(z1,i1)+cln_coli(z2,i2)
        if(abs(pz).lt.1d0) then
          if(dreal(pz).gt.0.5d0) then
            csp2cos_coli = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &          *(cln_coli(pz,ip)-sumln)
          else
            csp2cos_coli = pi2_6-cspenc_coli(pz,ip)
     &          -cln_coli(1d0-pz,-ip)*sumln 
          endif
        elseif(pz.ne.1d0) then
          if(dreal(pz).lt.2d0) then
#ifdef WIB
            csp2cos_coli = -cspenc_coli(1d0-1d0/pz,ip)
     &          +cln_coli(1d0-1d0/pz,ip)
     &          *(cln_coli(pz,ip)-sumln)
     &          -.5d0*(sumln)**2
#else
            csp2cos_coli = -cspenc_coli(1d0-1d0/pz,ip)
     &          +cln_coli(1d0-pz,-ip)
     &          *(cln_coli(pz,ip)-sumln)
     &        -.5d0*(cln_coli(pz,ip))**2
#endif
          else
#ifdef WIB
            csp2cos_coli = -pi2_6+cspenc_coli(1d0/pz,-ip)
     &          -cln_coli(1d0-1d0/pz,ip)*sumln
     &          -.5d0*(sumln)**2
#else
            csp2cos_coli = 2d0*pi2_6+cspenc_coli(1d0/pz,-ip)
     &          -cln_coli(1d0-pz,-ip)*sumln
     &        +.5d0*(cln_coli(-pz,-ip))**2
#endif
          endif
        else
          csp2cos_coli =0d0
          if(sumln.ne.cd0) write(*,*) 
     &        'csp2cos_coli: continuation possibly wrong'
        endif
      endif

#ifdef CHECK
      if(conscheck) then
        ip=-1d0
        if(pz.ne.1d0.and.pz.ne.0d0)then
          csp2cos_check = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &        *(cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2))
#ifdef WIB
          if(cdabs(pz).gt.1d0) then
            csp2cos_check = csp2cos_check - 
     &          (cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2))*
     &          (cln_coli(-pz,-ip)-.5d0*(cln_coli(pz,ip) 
     &          +cln_coli(z1,i1)+cln_coli(z2,i2)))
          endif 
#endif
        elseif(pz.eq.1d0)then
          csp2cos_check = 0d0
        else
          csp2cos_check = pi2_6
        endif 

c          write(*,*) 'csp2cos_coli out:',csp2cos_check

        if(abs(csp2cos_coli-csp2cos_check).gt.testacc*abs(csp2cos_coli)
     &      .and.abs(csp2cos_coli).gt.testacc) 
     &      then
          write(testout,*) 'csp2cos_coli inconsistent'
          write(testout,100) 'csp2cos_coli z1  = ',z1,i1
          write(testout,100) 'csp2cos_coli z2  = ',z2,i2
          write(testout,100) 'csp2cos_coli pz  = ',pz,ip
          write(testout,100) 'csp2cos_coli ompz= ',ompz
          write(testout,100) 'csp2cos_coli   = ',csp2cos_coli
          write(testout,100) 'csp2cos_check  = ',csp2cos_check
          write(testout,100) 'diff        = ',
     &        cdabs(csp2cos_coli/csp2cos_check-1d0)
          write(*,*) 'csp2cos_coli inconsistent'
          write(*,100) 'csp2cos_coli z1  = ',z1,i1
          write(*,100) 'csp2cos_coli z2  = ',z2,i2
          write(*,100) 'csp2cos_coli pz  = ',pz,ip
          write(*,100) 'csp2cos_coli ompz= ',ompz
          write(*,100) 'csp2cos_coli input = ',z1,z2
          write(*,100) 'csp2cos_coli inputi= ',i1,i2
          write(*,100) 'csp2cos_coli   = ',csp2cos_coli
          write(*,100) 'csp2cos_check  = ',csp2cos_check
          write(*,100) 'diff        = ',
     &        cdabs(csp2cos_coli/csp2cos_check-1d0)
          stop
        endif
      endif
#endif

c      write(*,*) 'csp2cos_coli out:',csp2cos_coli,csp2cos_check

      end

************************************************************************
      function csp3con_coli(z1,z2,z3,pz,ompz,i1,i2,i3)
************************************************************************
*  complex spence function   sp(1-z1*z2*z3)                            *
*                          plus continuation terms                     *
*  pz = z1*z2*z3, ompz = 1-pz                                          * 
*  i1,i2,i3 determine signs of infinitesimal imaginary parts of        *
*  z1,z2,z3                                                            *
*----------------------------------------------------------------------*
*  30.04.08 Ansgar Denner      last changed 03.06.08                   *
************************************************************************
      implicit   none
      complex*16 z1,z2,z3,pz,ompz
      real*8     i1,i2,i3,ip
      complex*16 csp3con_coli,cspenc_coli,cln_coli,sumln

#ifdef CHECK
      complex*16 csp3con_check
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

#ifdef CHECK
      if(argcheck) then
        if(cdabs(z1*z2*z3-pz).gt.1d3*calacc*cdabs(pz).or.
     &         cdabs(1d0-pz-ompz).gt.
     &      1d3*calacc*max(cdabs(pz),cdabs(1d0-pz))) then
          write(*,*) 'csp3con_coli:  arguments inconsistent  '
          write(*,100) 'csp3con_coli: z1 = ',z1
          write(*,100) 'csp3con_coli: z2 = ',z2
          write(*,100) 'csp3con_coli: z3 = ',z3
          write(*,100) 'csp3con_coli: zp = ',pz,z1*z2*z3
          write(*,100) 'csp3con_coli: omzp = ',ompz,1d0-pz
          write(*,100) 'test1 ',cdabs(z1*z2*z3-pz),1d3*calacc*cdabs(pz),
     &        cdabs(z1*z2*z3-pz).gt.1d3*calacc*cdabs(pz)
          write(*,100) 'test2 ',cdabs(1d0-pz-ompz),
     &        1d3*calacc*max(cdabs(pz),cdabs(1d0-pz)),
     &        cdabs(1d0-pz-ompz).gt.
     &      1d3*calacc*max(cdabs(pz),cdabs(1d0-pz))
 100      format(a21,4g25.17)
          stop
        endif
      endif
#endif

c      write(*,*) 'csp3con_coli in: ',z1,z2,z3,pz,ompz,i1,i2,i3,ip

      ip=1d0
      if(pz.eq.0d0)then
        csp3con_coli = pi2_6
      else
        sumln=cln_coli(z1,i1)+cln_coli(z2,i2)+cln_coli(z3,i3)
        if(abs(pz).lt.1d0) then
          if(dreal(pz).gt.0.5d0) then
            csp3con_coli = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &          *(cln_coli(pz,ip)-sumln)
          else
            csp3con_coli = pi2_6-cspenc_coli(pz,ip)
     &          -cln_coli(1d0-pz,-ip)*sumln 
          endif
        elseif(pz.ne.(1d0,0d0)) then
          if(dreal(pz).lt.2d0) then
#ifdef WIB
            csp3con_coli = -cspenc_coli(1d0-1d0/pz,ip)
     &          +cln_coli(1d0-1d0/pz,ip)
     &          *(cln_coli(pz,ip)-sumln)
     &          -.5d0*(sumln)**2
#else
            csp3con_coli = -cspenc_coli(1d0-1d0/pz,ip)
     &          +cln_coli(1d0-pz,-ip)
     &          *(cln_coli(pz,ip)-sumln)
     &        -.5d0*(cln_coli(pz,ip))**2
#endif
          else
#ifdef WIB
            csp3con_coli = -pi2_6+cspenc_coli(1d0/pz,-ip)
     &          -cln_coli(1d0-1d0/pz,ip)*sumln
     &          -.5d0*(sumln)**2
#else
            csp3con_coli = 2d0*pi2_6+cspenc_coli(1d0/pz,-ip)
     &          -cln_coli(1d0-pz,-ip)*sumln
     &        +.5d0*(cln_coli(-pz,-ip))**2 
#endif
          endif
        else
          csp3con_coli =0d0
        endif
      endif

#ifdef CHECK
      if(conscheck) then
        ip=-1d0
        if(pz.ne.1d0.and.pz.ne.0d0)then
          csp3con_check = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &        *(cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2)
     &        -cln_coli(z3,i3))
#ifdef WIB
          if(cdabs(pz).gt.1d0) then
            csp3con_check = csp3con_check -
     &          (cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2)
     &          -cln_coli(z3,i3))*
     &          (cln_coli(-pz,-ip)-.5d0*(cln_coli(pz,ip) 
     &          +cln_coli(z1,i1)+cln_coli(z2,i2)+cln_coli(z3,i3)))
          endif 
#endif
        elseif(pz.eq.1d0)then
          csp3con_check = 0d0
        else
          csp3con_check = pi2_6
        endif 
c      write(*,*) 'csp3con_coli out:',csp3con_check
        
        if(cdabs(csp3con_coli/csp3con_check-1d0).gt.1d1*testacc) then
          write(testout,*) 'csp3con_coli inconsistent'
          write(testout,100) 'csp3con_coli z1  = ',z1,i1
          write(testout,100) 'csp3con_coli z2  = ',z2,i2
          write(testout,100) 'csp3con_coli z3  = ',z3,i3
          write(testout,100) 'csp3con_coli pz  = ',pz,ip
          write(testout,100) 'csp3con_coli ompz= ',ompz
          write(testout,100) 'csp3con_coli     = ',csp3con_coli
          write(testout,100) 'csp3con_check    = ',csp3con_check
          write(testout,100) 'diff        = ',
     &        cdabs(csp3con_coli/csp3con_check-1d0)
          write(*,*) 'csp3con_coli inconsistent'
          write(*,100) 'csp3con_coli z1  = ',z1,i1
          write(*,100) 'csp3con_coli z2  = ',z2,i2
          write(*,100) 'csp3con_coli z3  = ',z3,i3
          write(*,100) 'csp3con_coli pz  = ',pz,ip
          write(*,100) 'csp3con_coli ompz= ',ompz
          write(*,100) 'csp3con_coli   = ',csp3con_coli
          write(*,100) 'csp3con_check  = ',csp3con_check
          write(*,100) 'diff        = ',
     &        cdabs(csp3con_coli/csp3con_check-1d0)
          stop
        endif
      endif
#endif

c      write(*,*) 'csp3con_coli out:',csp3con_coli,csp3con_check

      end

************************************************************************
      function csp4con_coli(z1,z2,z3,z4,pz,ompz,i1,i2,i3,i4)
************************************************************************
*  complex spence function   sp(1-z1*z2*z3)                            *
*                          plus continuation terms                     *
*  pz = z1*z2*z3*z4, ompz = 1-pz                                       * 
*  i1,i2,i3,i4 determine signs of infinitesimal imaginary parts of     *
*  z1,z2,z3,z4                                                         *
*----------------------------------------------------------------------*
*  02.05.08 Ansgar Denner      last changed 03.06.08                   *
************************************************************************
      implicit   none
      complex*16 z1,z2,z3,z4,pz,ompz
      real*8     i1,i2,i3,i4,ip
      complex*16 csp4con_coli,cspenc_coli,cln_coli,sumln

#ifdef CHECK
      complex*16 csp4con_check
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) 'csp4con_coli in: ',z1,z2,z3,z4,pz,ompz,i1,i2,i3,i4,ip

#ifdef CHECK
      if(argcheck) then
        if(cdabs(z1*z2*z3*z4-pz).gt.1d3*calacc*cdabs(pz).or.
     &         cdabs(1d0-pz-ompz).gt.
     &      1d3*calacc*max(cdabs(pz),cdabs(1d0-pz))) then
          write(*,*) 'csp4con_coli:  arguments inconsistent  '
          write(*,100) 'csp4con_coli: z1 = ',z1
          write(*,100) 'csp4con_coli: z2 = ',z2
          write(*,100) 'csp4con_coli: z3 = ',z3
          write(*,100) 'csp4con_coli: z4 = ',z4
          write(*,100) 'csp4con_coli: zp = ',pz,z1*z2*z3*z4
          write(*,100) 'csp4con_coli: omzp = ',ompz,1d0-pz
          write(*,100) 'test1 ',cdabs(z1*z2*z3*z4-pz),
     &        1d3*calacc*cdabs(pz)
          write(*,100) 'test2 ',cdabs(1d0-pz-ompz),
     &        1d3*calacc*max(cdabs(pz),cdabs(1d0-pz))
 100      format(a21,4g25.17)
          stop
        endif
      endif
#endif

      ip=1d0
      if(pz.eq.0d0)then
         csp4con_coli =pi2_6
      else
        sumln=cln_coli(z1,i1)+cln_coli(z2,i2)
     &      +cln_coli(z3,i3)+cln_coli(z4,i4)
        if(abs(pz).lt.1d0) then
          if(dreal(pz).gt.0.5d0) then
            csp4con_coli = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &          *(cln_coli(pz,ip)-sumln)
          else
            csp4con_coli = pi2_6-cspenc_coli(pz,ip)
     &          -cln_coli(1d0-pz,-ip)*sumln 
          endif
        elseif(pz.ne.(1d0,0d0)) then
          if(dreal(pz).lt.2d0) then
#ifdef WIB
            csp4con_coli = -cspenc_coli(1d0-1d0/pz,ip)
     &          +cln_coli(1d0-1d0/pz,ip)
     &          *(cln_coli(pz,ip)-sumln)
     &        -.5d0*(sumln)**2
#else
            csp4con_coli = -cspenc_coli(1d0-1d0/pz,ip)
     &          +cln_coli(1d0-pz,-ip)
     &          *(cln_coli(pz,ip)-sumln)
     &        -.5d0*(cln_coli(pz,ip))**2
#endif
          else
#ifdef WIB
            csp4con_coli = -pi2_6+cspenc_coli(1d0/pz,-ip)
     &          -cln_coli(1d0-1d0/pz,ip)*sumln
     &          -.5d0*(sumln)**2
#else
            csp4con_coli = 2d0*pi2_6+cspenc_coli(1d0/pz,-ip)
     &          -cln_coli(1d0-pz,-ip)*sumln
     &        +.5d0*(cln_coli(-pz,-ip))**2
#endif
          endif
        else
          csp4con_coli =0d0
        endif
      endif

#ifdef CHECK
      if(conscheck) then
        ip=-1d0
        if(pz.ne.1d0.and.pz.ne.0d0)then
          csp4con_check = cspenc_coli(1d0-pz,-ip)+cln_coli(1d0-pz,-ip)
     &        *(cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2)
     &        -cln_coli(z3,i3) -cln_coli(z4,i4))
#ifdef WIB
          if(cdabs(pz).gt.1d0) then
            csp4con_check = csp4con_check -
     &          (cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2)
     &          -cln_coli(z3,i3)-cln_coli(z4,i4) )*
     &          (cln_coli(-pz,-ip)-.5d0*(cln_coli(pz,ip) 
     &          +cln_coli(z1,i1)+cln_coli(z2,i2)+cln_coli(z3,i3)
     &          +cln_coli(z4,i4)))
          endif
#endif
        elseif(pz.eq.1d0)then
          csp4con_check = 0d0
        else
          csp4con_check = pi2_6
        endif 
        
        if(cdabs(csp4con_coli/csp4con_check-1d0).gt.testacc) then
          write(testout,100) 'csp4con_coli inconsistent'
          write(testout,100) 'csp4con_coli z1  = ',z1,i1
          write(testout,100) 'csp4con_coli z2  = ',z2,i2
          write(testout,100) 'csp4con_coli z3  = ',z3,i3
          write(testout,100) 'csp4con_coli z4  = ',z4,i4
          write(testout,100) 'csp4con_coli pz  = ',pz,ip
          write(testout,100) 'csp4con_coli ompz= ',ompz
          write(testout,100) 'csp4con_coli   = ',csp4con_coli
          write(testout,100) 'csp4con_check  = ',csp4con_check
          write(testout,100) 'diff        = ',
     &        cdabs(csp4con_coli/csp4con_check-1d0)
          write(*,*) 'csp4con_coli inconsistent'
          write(*,100) 'csp4con_coli z1  = ',z1,i1
          write(*,100) 'csp4con_coli z2  = ',z2,i2
          write(*,100) 'csp4con_coli z3  = ',z3,i3
          write(*,100) 'csp4con_coli z4  = ',z4,i4
          write(*,100) 'csp4con_coli pz  = ',pz,ip
          write(*,100) 'csp4con_coli ompz= ',ompz
          write(*,100) 'csp4con_coli   = ',csp4con_coli
          write(*,100) 'csp4con_check  = ',csp4con_check
          write(*,100) 'diff        = ',
     &        cdabs(csp4con_coli/csp4con_check-1d0)
          stop
        endif
      endif
#endif

      end

************************************************************************
      function eta2_coli(z1,z2,pz,i1,i2,ip)
************************************************************************
*     complex eta-function with explicit infinitesimal im parts        *
*     pz = z1*z2                                                      *
*----------------------------------------------------------------------*
*     06.10.08    Ansgar Denner       last changed   06.10.08          *
************************************************************************
      implicit   none
      complex*16 eta2_coli,z1,z2,pz
      real*8     i1,i2,ip
      real*8     im1,im2,imp,re1,re2,rep

#ifdef CHECK
      complex*16 eta2_check,cln_coli
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) 'eta2_coli in',z1,z2,pz,i1,i2,ip,calacc,pi

#ifdef CHECK
      if(argcheck) then
        if(cdabs(z1*z2-pz).gt.1d3*calacc*cdabs(pz)) then
          write(*,*) 'eta2_coli:  arguments inconsistent  '
          write(*,100) 'eta2_coli: z1 = ',z1
          write(*,100) 'eta2_coli: z2 = ',z2
          write(*,100) 'eta2_coli: zp = ',pz,z1*z2
          write(*,100) 'test1 ',cdabs(z1*z2-pz),
     &        1d3*calacc*cdabs(pz)
 100      format(a21,4g25.17)
          stop
        endif
      endif
#endif

      im1   = dimag(z1)
      im2   = dimag(z2)
      imp   = dimag(pz)
      re1   = dreal(z1)
      re2   = dreal(z2)
      rep   = dreal(pz)

c      write(*,*) 'eta2_coli im',im1,im2,im12

      if (abs(im1).lt.impacc*abs(re1)) im1 = 0d0
      if (abs(im2).lt.impacc*abs(re2)) im2 = 0d0
      if (abs(imp).lt.impacc*abs(rep)) imp = 0d0

c      write(*,*) 'eta2_coli im',im1,im2,im12

      if (im1.ne.0d0.or.im2.ne.0d0.or.imp.ne.0d0) then
        if(im1.eq.0d0) im1 = i1
        if(im2.eq.0d0) im2 = i2
        if(imp.eq.0d0) imp = ip
        if(im1.lt.0d0.and.im2.lt.0d0.and.imp.gt.0d0) then
          eta2_coli = dcmplx(0d0,2d0*pi)
        else if (im1.gt.0d0.and.im2.gt.0d0.and.imp.lt.0d0) then
          eta2_coli = dcmplx(0d0,-2d0*pi)
        else
          eta2_coli = dcmplx(0d0)
          if(im1.eq.0.and.re1.lt.0d0 .or.
     &         im2.eq.0.and.re2.lt.0d0 .or.
     &         imp.eq.0.and.rep.lt.0d0) then
            write(*,*) ' eta2_coli not defined '
            write(*,*) ' eta2_coli:  z1 = ',z1,i1
            write(*,*) ' eta2_coli:  z2 = ',z2,i2
            write(*,*) ' eta2_coli:  pz = ',pz,ip
            stop
          end if
        end if
      else 
        eta2_coli = dcmplx(0d0)      
        if (re1.lt.0d0.and.re2.lt.0d0) then
          if(i1.gt.0d0.and.i2.gt.0d0) then
            eta2_coli = dcmplx(0d0,-2d0*pi)
          else   if(i1.lt.0d0.and.i2.lt.0d0) then
            eta2_coli = dcmplx(0d0,2d0*pi)
          end if
        else if (rep.lt.0d0.and.re1.lt.0d0) then
          if(i1.gt.0d0.and.ip.lt.0d0) then
            eta2_coli = dcmplx(0d0,-2d0*pi)
          else   if(i1.lt.0d0.and.ip.gt.0d0) then
            eta2_coli = dcmplx(0d0,2d0*pi)
          end if
        else if (rep.lt.0d0.and.re2.lt.0d0) then
          if(i2.gt.0d0.and.ip.lt.0d0) then
            eta2_coli = dcmplx(0d0,-2d0*pi)
          else   if(i2.lt.0d0.and.ip.gt.0d0) then
            eta2_coli = dcmplx(0d0,2d0*pi)
          end if
        end if
      end if

#ifdef CHECK
      if(conscheck) then
        eta2_check = cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2)
        
        if(cdabs(eta2_coli-eta2_check).gt.testacc) then
          write(testout,100) 'eta2_coli inconsistent'
          write(testout,100) 'eta2_coli z1  = ',z1,i1,im1
          write(testout,100) 'eta2_coli z2  = ',z2,i2,im2
          write(testout,100) 'eta2_coli pz  = ',pz,ip,imp,z1*z2
          write(testout,100) 'eta2_coli   = ',eta2_coli
          write(testout,100) 'eta2_check  = ',eta2_check
          write(testout,100) 'diff        = ',
     &        cdabs(eta2_coli/eta2_check-1d0)
          write(*,*) 'eta2_coli inconsistent'
          write(*,100) 'eta2_coli z1  = ',z1,i1,im1
          write(*,100) 'eta2_coli z2  = ',z2,i2,im2
          write(*,100) 'eta2_coli pz  = ',pz,ip,imp
          write(*,100) 'eta2_coli   = ',eta2_coli
          write(*,100) 'eta2_check  = ',eta2_check
          write(*,100) 'diff        = ',
     &        cdabs(eta2_coli/eta2_check-1d0)
          stop
        endif
      endif
#endif

      end

************************************************************************
      function eta2s_coli(z1,z2,i1,i2,ip)
************************************************************************
*     complex eta-function with explicit infinitesimal im parts        *
*     pz = z1*z2                                                       *
*----------------------------------------------------------------------*
*     06.10.08    Ansgar Denner       last changed   06.10.08          *
************************************************************************
      implicit   none
      complex*16 eta2s_coli,z1,z2,pz
      real*8     i1,i2,ip
      real*8     im1,im2,imp,re1,re2,rep

#ifdef CHECK
      complex*16 eta2s_check,cln_coli
#endif

#include "params_coli.h"
#include "checkparams_coli.h"

c      write(*,*) 'eta2s_coli in',z1,z2,pz,i1,i2,ip,calacc

      pz=z1*z2

      im1 = dimag(z1)
      im2 = dimag(z2)
      imp = dimag(pz)
      re1 = dreal(z1)
      re2 = dreal(z2)
      rep = dreal(pz)

c      write(*,*) 'eta2s_coli im',im1,im2,imp
c      write(*,*) 'eta2s_coli re',re1,re2,rep

      if (abs(im1).lt.impacc*abs(re1)) im1 = 0d0
      if (abs(im2).lt.impacc*abs(re2)) im2 = 0d0
      if (abs(imp).lt.impacc*abs(rep)) imp = 0d0

c      write(*,*) 'eta2s_coli im',im1,im2,imp

      if (im1.ne.0d0.or.im2.ne.0d0.or.imp.ne.0d0) then
        if(im1.eq.0d0) im1 = i1
        if(im2.eq.0d0) im2 = i2
        if(imp.eq.0d0) imp = ip
        if(im1.lt.0d0.and.im2.lt.0d0.and.imp.gt.0d0) then
          eta2s_coli = dcmplx(0d0,2d0*pi)
        else if (im1.gt.0d0.and.im2.gt.0d0.and.imp.lt.0d0) then
          eta2s_coli = dcmplx(0d0,-2d0*pi)
        else
          eta2s_coli = dcmplx(0d0)
          if(im1.eq.0.and.re1.lt.0d0 .or.
     &         im2.eq.0.and.re2.lt.0d0 .or.
     &         imp.eq.0.and.rep.lt.0d0) then
            write(*,*) ' eta2s_coli not defined '
            write(*,*) ' eta2s_coli:  z1 = ',z1,i1
            write(*,*) ' eta2s_coli:  z2 = ',z2,i2
            write(*,*) ' eta2s_coli:  pz = ',pz,ip
            stop
          end if
        end if
      else 
        eta2s_coli = dcmplx(0d0)      
        if (re1.lt.0d0.and.re2.lt.0d0) then
          if(i1.gt.0d0.and.i2.gt.0d0) then
            eta2s_coli = dcmplx(0d0,-2d0*pi)
          else   if(i1.lt.0d0.and.i2.lt.0d0) then
            eta2s_coli = dcmplx(0d0,2d0*pi)
          end if
        else if (rep.lt.0d0.and.re1.lt.0d0) then
          if(i1.gt.0d0.and.ip.lt.0d0) then
            eta2s_coli = dcmplx(0d0,-2d0*pi)
          else   if(i1.lt.0d0.and.ip.gt.0d0) then
            eta2s_coli = dcmplx(0d0,2d0*pi)
          end if
        else if (rep.lt.0d0.and.re2.lt.0d0) then
          if(i2.gt.0d0.and.ip.lt.0d0) then
            eta2s_coli = dcmplx(0d0,-2d0*pi)
          else   if(i2.lt.0d0.and.ip.gt.0d0) then
            eta2s_coli = dcmplx(0d0,2d0*pi)
          end if
        end if
      end if

#ifdef CHECK
      if(conscheck) then
        eta2s_check = cln_coli(pz,ip)-cln_coli(z1,i1)-cln_coli(z2,i2)
        
        if(abs(eta2s_coli-eta2s_check).gt.testacc) then
 100      format(a21,6g25.17)
          write(testout,100) 'eta2s_coli inconsistent'
          write(testout,100) 'eta2s_coli z1  = ',z1,i1,im1
          write(testout,100) 'eta2s_coli z2  = ',z2,i2,im2
          write(testout,100) 'eta2s_coli pz  = ',pz,ip,imp,z1*z2
          write(testout,100) 'eta2s_coli   = ',eta2s_coli
          write(testout,100) 'eta2s_check  = ',eta2s_check
          write(testout,100) 'diff        = ',
     &        abs(eta2s_coli/eta2s_check-1d0)
          write(*,*) 'eta2s_coli inconsistent'
          write(*,100) 'eta2s_coli z1  = ',z1,i1,im1
          write(*,100) 'eta2s_coli z2  = ',z2,i2,im2
          write(*,100) 'eta2s_coli pz  = ',pz,ip,imp,z1*z2
          write(*,100) 'eta2s_coli   = ',eta2s_coli
          write(*,100) 'eta2s_check  = ',eta2s_check
          write(*,100) 'eta2s_check c= ',cln_coli(pz,ip)
          write(*,100) 'eta2s_check 1= ',-cln_coli(z1,i1)
          write(*,100) 'eta2s_check 2= ',-cln_coli(z2,i2)
          write(*,100) 'diff        = ',
     &        abs(eta2s_coli/eta2s_check-1d0)
          stop
        endif
      endif
#endif

      end
************************************************************************
      subroutine setmuuv2_coli(muuv2in)
************************************************************************
*     sets UV-regulator mass squared to muuv2in                        *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      real*8     muuv2in

#include "common_coli.h"

      muuv2 = muuv2in

      if(coliinfo)then
        write(*,*) 'coli: UV regulator mass squared set to muuv2 = ',
     &      muuv2
      endif

      end
************************************************************************
      subroutine getmuuv2_coli(muuv2out)
************************************************************************
*     returns UV-regulator mass squared to muuv2in                     *
*-----------------------------------------------------------------------
*     23.10.08 Ansgar Denner         last changed  23.10.08            *
************************************************************************
      implicit   none
      real*8     muuv2out

#include "common_coli.h"

      muuv2out = muuv2

      end

************************************************************************
      subroutine setmuir2_coli(muir2in)
************************************************************************
*     sets photon mass squared regulator to lambda2in                  *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      real*8     muir2in

#include "common_coli.h"

      muir2 = muir2in

      if(coliinfo)then
        write(*,*) 'coli: IR regularor mass squared set to muir2 = ',
     &      muir2
      endif

      end
************************************************************************
      subroutine getmuir2_coli(muir2out)
************************************************************************
*     returns UV-regulator mass squared to muir2in                     *
*-----------------------------------------------------------------------
*     23.10.08 Ansgar Denner         last changed  23.10.08            *
************************************************************************
      implicit   none
      real*8     muir2out

#include "common_coli.h"

      muir2out = muir2

      end


#ifdef SING
************************************************************************
      subroutine setdeltauv_coli(deltauvin)
************************************************************************
*     sets UV pole 2/(4-D) to deltauv                                  *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      real*8     deltauvin

#include "common_coli.h"

      deltauv = deltauvin
      if(coliinfo)then
        write(*,*) 'coli: UV pole set to                 deltauv = ',
     &    deltauv
      endif

      end
************************************************************************
      subroutine getdeltauv_coli(deltauvout)
************************************************************************
*     returns UV pole 2/(4-D)  deltauv                                 *
*-----------------------------------------------------------------------
*     23.10.08 Ansgar Denner         last changed  23.10.08            *
************************************************************************
      implicit   none
      real*8     deltauvout

#include "common_coli.h"

      deltauvout = deltauv

      end
************************************************************************
      subroutine setdeltair_coli(delta1irin,delta2irin)
************************************************************************
*     sets IR pole 2/(4-D) to delta1ir and squared pole to delta2ir    *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      real*8     delta1irin,delta2irin

#include "common_coli.h"

      delta1ir = delta1irin
      delta2ir = delta2irin
      if(coliinfo)then
        write(*,*) 'coli: IR pole set to                delta1ir = ',
     &      delta1ir
        write(*,*) 'coli: squared IR pole set to        delta2ir = ',
     &      delta2ir
      endif
      end
************************************************************************
      subroutine getdeltair_coli(delta1irout,delta2irout)
************************************************************************
*     returns UV pole 2/(4-D)  deltauv                                 *
*-----------------------------------------------------------------------
*     23.10.08 Ansgar Denner         last changed  23.10.08            *
************************************************************************
      implicit   none
      real*8     delta1irout,delta2irout

#include "common_coli.h"

      delta1irout = delta1ir
      delta2irout = delta2ir

      end
************************************************************************
      subroutine setshiftms2_coli(shiftms2in)
************************************************************************
*     sets  shift for squared mass-singular logs to shiftms2in         *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      real*8     shiftms2in

#include "common_coli.h"

      colishiftms2 = shiftms2in

      if(coliinfo)then
        write(*,*) 'coli: shift for squared mass-singular logs set to' 
        write(*,*) '                               colishiftms2  = ',
     &      colishiftms2
      endif

      end

************************************************************************
      subroutine getshiftms2_coli(shiftms2out)
************************************************************************
*     gets shift for squared mass-singular logs colishiftms2           *
*-----------------------------------------------------------------------
*     15.09.08 Ansgar Denner         last changed  15.09.08            *
************************************************************************
      implicit   none
      real*8     shiftms2out

#include "common_coli.h"

      shiftms2out=colishiftms2 

      end
#endif

************************************************************************
      subroutine setminf2_coli(m2)
************************************************************************
*     adds m2 to list of mass squared to be neglected                  *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      complex*16 m2
      integer    i

#include "common_coli.h"

      data ncoliminf /0/

      do 10 i=1,ncoliminf
        if(m2.eq.coliminffix2(i)) return
 10   continue
      
      ncoliminf = ncoliminf+1
      coliminffix(ncoliminf) = sqrt(m2)
      coliminffix2(ncoliminf) = m2
      coliminf2(ncoliminf) = m2/coliminfscale2
      coliminf(ncoliminf) = sqrt(coliminf2(ncoliminf))

      if(coliinfo)then
        write(*,*) 'coli: added to small masses squared:  m2 =  ',m2
      endif     

      end
************************************************************************
      subroutine clearcoliminf2()
************************************************************************
*     clears list of masses squared to be neglected                    *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      integer    i

#include "common_coli.h"

      do 10 i=1,ncoliminf
        coliminffix2(i) = 0d0
        coliminffix(i) = 0d0
        coliminf2(i) = 0d0
        coliminf(i) = 0d0
 10   continue
      
      ncoliminf = 0
 
      if(coliinfo)then
        write(*,*) 'coli: list of complex small masses cleared'
      endif    

      end

************************************************************************
      subroutine setminfscale2_coli(scale2)
************************************************************************
*     sets factor with which small masses are scaled up in             *
*     mass-singular logarithms                                         *
*     and rescales masses in /coliminf/ /coliminf2/                    *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      real*8     scale2

#include "common_coli.h"

      integer    i
      coliminfscale2 = scale2
      coliminfscale = sqrt(scale2)

      do 10 i=1,ncoliminf
        coliminf(i) = coliminffix(i)/coliminfscale
        coliminf2(i) = coliminffix2(i)/coliminfscale2
 10   continue

      if(coliinfo)then
        write(*,*) 'coli: small masses scaled with coliminfscale = ',
     &      coliminfscale
        write(*,*) '      in mass-singular logarithms'
      endif

      end

************************************************************************
      subroutine modminfscale2_coli(scale2)
************************************************************************
*     sets factor with which small masses are scaled up in             *
*     mass-singular logarithms                                         *
*-----------------------------------------------------------------------
*     15.09.09 Ansgar Denner         last changed  15.09.09            *
************************************************************************
      implicit   none
      real*8     scale2

#include "common_coli.h"

      coliminfscale2 = scale2
      coliminfscale = sqrt(scale2)
      if(coliinfo)then
        write(*,*) 'coli: scalefactor for small masses ',
     &      ' in mass-singular logarithms'
        write(*,*) 
     &      '     changed to coliminfscale = ',
     &      coliminfscale
      endif

      end

************************************************************************
      subroutine getminfscale2_coli(scale2)
************************************************************************
*     gets factor with which small masses are scaled up in             *
*     mass-singular logarithms                                         *
*-----------------------------------------------------------------------
*     15.09.09 Ansgar Denner         last changed  15.09.09            *
************************************************************************
      implicit   none
      real*8     scale2

#include "common_coli.h"

      scale2=coliminfscale2
      end

************************************************************************
      function minfscaledown_coli(m)
************************************************************************
*     scales mass down by coliminfscale if m in /coliminf/             *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      complex*16 minfscaledown_coli,m

#include "common_coli.h"

      integer    i

      minfscaledown_coli = m
      do 10 i=1,ncoliminf
         if(m.eq.coliminffix(i)) then
           minfscaledown_coli = coliminf(i)
           return
         end if 
 10   continue

      end

************************************************************************
      function minfscaledown2_coli(m2) 
************************************************************************
*     scales mass down by coliminfscale^2 if m2 in /coliminf2/         *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      complex*16 minfscaledown2_coli,m2

#include "common_coli.h"

      integer    i
      minfscaledown2_coli = m2
      do 10 i=1,ncoliminf
         if(m2.eq.coliminffix2(i)) then
           minfscaledown2_coli = coliminf2(i)
           return
         end if 
 10   continue

      end

************************************************************************
      function elimminf_coli(m)
************************************************************************
*     if m = coliminfs(i)  then  elimcoliminf = 0                      *
*                          else  elimcoliminf  = m                     *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  08.09.08            *
************************************************************************
      implicit   none
      complex*16 m,elimminf_coli

#include "common_coli.h"

      integer    i

c      write(*,*) 'elimcoliminf',m 
c      write(*,*) 'elimcoliminf',coliminf
      elimminf_coli = m
      do 10 i=1,ncoliminf
         if(m.eq.coliminf(i)) then
           elimminf_coli = 0d0
           return
         end if 
 10   continue
c      write(*,*) 'elimminf_coli',elimminf_coli,m-coliminf(1)
      end

************************************************************************
      function elimminf2_coli(m2)
************************************************************************
*     if m2 = coliminf2(i)  then  elimminf2_coli = 0                   *
*                           else  elimminf2_coli = m2                  *
*-----------------------------------------------------------------------
*     28.08.08 Ansgar Denner         last changed  28.08.08            *
************************************************************************
      implicit   none
      complex*16 m2,elimminf2_coli

#include "common_coli.h"

      integer    i

c      write(*,*) 'elimminf2_coli in =',m2

      elimminf2_coli = m2
      do 10 i=1,ncoliminf           
         if(m2.eq.coliminf2(i)) then
           elimminf2_coli = 0d0
           return
         end if 
 10   continue

c      write(*,*) 'elimminf2_coli out =',m2

      end

************************************************************************
      subroutine setinfo_coli
************************************************************************
*     sets photon mass squared regulator to lambda2in                  *
*-----------------------------------------------------------------------
*     15.09.09 Ansgar Denner         last changed  15.09.09            *
************************************************************************
      implicit   none

#include "common_coli.h"

      coliinfo=.true.
      write(*,*) 'coli: information printed'

      end
************************************************************************
      subroutine unsetinfo_coli
************************************************************************
*     sets photon mass squared regulator to lambda2in                  *
*-----------------------------------------------------------------------
*     15.09.09 Ansgar Denner         last changed  15.09.09            *
************************************************************************
      implicit   none

#include "common_coli.h"

      coliinfo=.false.
      write(*,*) 'coli: no information printed'

      end
