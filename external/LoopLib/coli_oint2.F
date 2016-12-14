************************************************************************
*                                                                      *
*     scalar functions: interface to old routines for testing          *
*                                                                      *
************************************************************************
*                                                                      *
*     last changed  28.10.08  Ansgar Denner                            *
*                                                                      *
************************************************************************
* subroutines:                                                         *
* crootse,celminfiv,cerrB0                                             *
* functions:                                                           *
* cdln, cdlnr, olieta, olietal, olietae, olietae2, olieta2l, ettile    *
* olicspenc, olicspenh, olicspen, cspcon, cspcoe                       *
* cDB00f                                                               * 
* fpve,  yfpve, xyfpve,xlogxe                                          *
************************************************************************

************************************************************************
      function cdlnr(x,eps)
************************************************************************
*     complex logarithm of x + i*eps    real x                         *
*----------------------------------------------------------------------*
*     09.01.90 Ansgar Denner        last changed 21.05.04              *
************************************************************************
      implicit   none
      real*8     x,eps
      complex*16 cdlnr
      real*8     pi,pi2_6
      common /pietc/  pi,pi2_6

      if(x.le.0d0)then
        if (eps.eq.0d0) then
          write(*,*) 'cdlnr:  argument on cut '
          write(*,*) 'cdlnr:  eps = 0'
          write(*,*) 'cdlnr:  x   = ',x 
          cdlnr = log(-x)+dcmplx(0d0,pi)
        else
          cdlnr = log(-x)+dcmplx(0d0,pi)*dsign(1d0,eps)
        end if 
      else
        cdlnr = log(x)
      end if
      end

************************************************************************
      function cdln(cz,eps)
************************************************************************
*     complex logarithm of cz + i*eps                                  *
*----------------------------------------------------------------------*
*     09.01.90 Ansgar Denner         last changed 15.03.04             *
************************************************************************
      implicit   none
      real*8     eps
      complex*16 cdln,cz
      real*8     pi,pi2_6
      common /pietc/  pi,pi2_6
      real*8     impacc
      common /impacc/ impacc
c      real*8     calacc
c      common /calacc/ calacc

c changed 15.03.04
c>      if(dimag(cz).eq.0d0.and.dreal(cz).le.0.d0)then
c>        if (eps.eq.0d0) then
c>          write(80,*) 'cdln:  argument on cut '
c>          write(80,*) 'cdln:  eps = 0'
c>          write(80,*) 'cdln:  cz  = ',cz
c>        end if
c>        cdln=cdlog(-cz)+dcmplx(0d0,pi)*dsign(1d0,eps)
c>      else
c>        cdln=cdlog(cz)
c>      end if

c     if(dreal(cz).lt.0d0.and.
c    &      abs(dimag(cz)/dreal(cz)).lt.impacc) then
c       if (eps.ne.0d0) then
c         cdln=cdlog(-cz)+dcmplx(0d0,pi)*dsign(1d0,eps)
c       else 
c         cdln=cdlog(-cz)+dcmplx(0d0,pi)
c        write(*,*) 'cdln:  argument on cut '
c        write(*,*) 'cdln:  eps = 0'
c        write(*,*) 'cdln:  cz  = ',cz
c       end if
c     else
c       cdln=cdlog(cz)
c     end if

c      write(*,*) 'cdln ',cz,eps,calacc

      if(eps.ne.0d0.and.dreal(cz).lt.0d0.and.
     &      abs(dimag(cz)/dreal(cz)).lt.impacc) then
        cdln=cdlog(-cz)+dcmplx(0d0,pi)*dsign(1d0,eps)
      else if (eps.eq.0d0.and.dimag(cz).eq.0d0.and.dreal(cz).le.0d0)
     &       then
        cdln=cdlog(-cz)+dcmplx(0d0,pi)
        write(*,*) 'cdln:  argument on cut '
        write(*,*) 'cdln:  eps = 0'
        write(*,*) 'cdln:  cz  = ',cz
      else
        cdln=cdlog(cz)
      end if

c      write(*,*) 'cdln = ',cz,eps,cdln
      end

************************************************************************
      function olicspenc(cz,eps)
************************************************************************
*       complex spence function  of cz + i*eps                         *
*       calculated by mapping on the area where there is a quickly     *
*       convergent series                                              *
*----------------------------------------------------------------------*
*     08.01.90 Ansgar Denner        last changed  3.02.97 ad           *
************************************************************************
      implicit   none
      real*8     eps
      complex*16 olicspenc,cz
      real*8     pi,pi2_6
      real*8     az,rz,az1
      complex*16 cz1,olicspenh,cdln

      common /pietc/  pi,pi2_6

c     pi      = 4d0*datan(1d0)
c     pi2_6   = pi*pi/6d0
c     pi2_6   = 1.64493406684822643d0
      cz1     = 1d0-cz
      az1     = cdabs(cz1)
      az      = cdabs(cz)
      rz      = dreal(cz)

      if (az1.lt.1d-15) then
         olicspenc = pi2_6
      else if (rz.lt.0.5d0) then
         if (az.lt.1d0) then
            olicspenc = olicspenh(cz,eps)
         else
            olicspenc = -pi2_6 - .5d0*cdln(-cz,-eps)**2
     &                - olicspenh(1d0/cz,-eps)
         end if
      else
         if (az1.lt.1d0) then
            olicspenc =  pi2_6 - cdln(cz,eps)*cdln(cz1,-eps)
     &                       - olicspenh(cz1,-eps)
         else
            olicspenc = 2d0*pi2_6 + .5d0*cdln(-cz1,-eps)**2
     &              - cdln(cz,eps)*cdln(cz1,-eps)
     &              + olicspenh(1d0/cz1,eps)
         end if
      end if
c      write(*,*) 'olicspenc = ',cz,eps,olicspenc
      end

************************************************************************
      function olicspenh(cz,eps)
************************************************************************
*       complex spence function of cz + i*eps                          *
*       in convergence region                                          *
*       calculation of bernoulli series                                *
*----------------------------------------------------------------------*
*     09.01.90 Ansgar Denner                                           *
************************************************************************
      implicit   none
      complex*16 olicspenh,cdln,cz,x,x2
      real*8     eps
c     real*8     b(11),eps
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
c     beachte:                 b(n)=b2n
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
      x      =  -cdln(1d0-cz,-eps)
c     write(80,*)  'olicspenh'
      x2     =  x*x
      power  =  x
      factor =  1d0
      olicspenh =  x - x2/4d0
      do 10 j=2,22,2
         factor = factor / j / (j+1)
         power  = power * x2
         term   = b(j/2) * factor * power
         csp    = olicspenh + term
         if (csp.eq.olicspenh) return
         olicspenh = csp
10    continue
      if (cdabs(term/csp).gt.1d-15) then
        write(80,*) 'olicspenh converges badly  ',cz,x
        write(80,*) 'olicspenh converges badly  ',csp-term,olicspenh
     &      ,term
      end if 
      end

************************************************************************
      function olicspen(cz)
************************************************************************
*       complex spence function                                        *
*----------------------------------------------------------------------*
*     08.07.94 Ansgar Denner                                           *
************************************************************************
      implicit   none
      complex*16 olicspen,olicspenc,cz
      real*8     impacc
      common /impacc/ impacc

      if((dimag(cz).eq.0d0).and.(dreal(cz-1d0).gt.impacc))then
        write(60,*) 'olicspen:  argument on cut '
        write(60,*) 'olicspen:  cz  = ',cz
        write(61,*) 'olicspen:  argument on cut '
        write(61,*) 'olicspen:  cz  = ',cz
        write(80,*) 'olicspen:  argument on cut '
        write(80,*) 'olicspen:  cz  = ',cz
      end if
      olicspen = olicspenc(cz,0d0)
      end

************************************************************************
      function cspcon(z1,z2)
************************************************************************
*  complex spence function plus continuation terms                     *
*----------------------------------------------------------------------*
*  29.03.92 Ansgar Denner      last changed 9.05.03                    *
************************************************************************
      implicit   none
      complex*16 cspcon,olicspen,olieta,z1,z2
      real*8     pi,pi2_6
      common /pietc/  pi,pi2_6
 
c     write(81,*) 'cspcon0 ',z1,z2,z1*z2
      if(dreal(z1*z2).gt.0d0) then
        cspcon = olicspen(1d0-z1*z2) + olieta(z1,z2)*log(1d0-z1*z2)
c       write(81,*) 'cspcon1 ', olicspen(1d0-z1*z2)
c       write(81,*) 'cspcon1 ', olieta(z1,z2),log(1d0-z1*z2)
c       write(81,*) 'cspcon1 ', olieta(z1,z2)*log(1d0-z1*z2)
c       write(81,*) 'cspcon1 ', cspcon
      else
c---> changed 9.05.03
c       cspcon = pi2_6-olicspen(z1*z2)
c    &           - (log(z1)+log(z2))*log(1d0-z1*z2)
        if(z1*z2.ne.0d0) then
          cspcon = pi2_6-olicspen(z1*z2)
     &             - (log(z1)+log(z2))*log(1d0-z1*z2)          
        else
          cspcon = pi2_6
        end if
c---< changed 9.05.03
      end if
      end
************************************************************************
      function cspcoe(z1,z2,i1,i2)
************************************************************************
*  complex spence function plus continuation terms                     *
*  i2 is assumed to dominate i1                                        *
*----------------------------------------------------------------------*
*  08.07.94 Ansgar Denner      last changed  17.03.95                  *
************************************************************************
      implicit   none
      complex*16 cspcoe,olicspenc,olietae,cdln,z1,z2,z12,olietaa
      real*8     i1,i2,i12
      real*8     pi,pi2_6
      common /pietc/  pi,pi2_6
 
      z12 = z1*z2
c      write(*,*) 'cspcoe ',z1,z2,z12
c --> changed  23.09.96
c     i12 = i2*sign(1d0,real(z1))
      i12 = i2*sign(1d0,dreal(z1))
c      write(*,*) 'cspcoe ',i1,i2,i12
c     write(80,*) 'cspcoe0 ',z1,z2,z12,i1,i2,i12
      if(dreal(z12).gt.0.5d0) then
        cspcoe = olicspenc(1d0-z12,0d0) 
        olietaa   = olietae(z1,z2,i1,i2,i12)
        if(olietaa.ne.0d0) then
          if(z12.ne.1d0)then
            cspcoe = cspcoe + olietaa*cdln(1d0-z12,-i12)
          else
            write(*,*) ' cspcoe: continuation possibly wrong'
          endif
        end if
c       write(80,*) 'cspcoe1 ',olicspenc(1d0-z12,0d0)
c       write(80,*) 'cspcoe1 ', olietaa,cdln(1d0-z12,-i12)
c       write(80,*) 'cspcoe1 ', olietaa*cdln(1d0-z12,-i12)
c added 23.5.02
      else if(abs(z12).eq.0d0) then
        cspcoe = pi2_6
      else if(abs(z12).lt.1d-4) then
        cspcoe = pi2_6-olicspenc(z12,0d0)
     &           + (cdln(z1,i1)+cdln(z2,i2))
     &                  *z12*(1d0+z12/2d0+z12*z12/3d0+z12*z12*z12/4d0)
      else
        cspcoe = pi2_6-olicspenc(z12,0d0)
     &           - (cdln(z1,i1)+cdln(z2,i2))
     &                  *cdln(1d0-z12,-0d0)
      end if
c     write(80,*) 'cspcoe '

      end
************************************************************************
      function ettile(c1,c2,i1,i2)
************************************************************************
*  complex eta-tilde   for 16-spence-d0                                *
*----------------------------------------------------------------------*
*  29.03.92 Ansgar Denner        last changed 12.05.97                 *
************************************************************************
      implicit   none
      complex*16 ettile,olietae,c1,c2
      real*8     i1,i2
      real*8     im1,im2,re2
      real*8     pi,pi2_6
      common /pietc/  pi,pi2_6
 
      im1    = dimag(c1)
      if(im1.eq.0d0) im1 = i1
      im2    = dimag(c2)
      re2    = dreal(c2)
      if(im2.ne.0d0) then
          ettile = olietae(c1,c2,i1,0d0,0d0)
      else if (re2.gt.0d0) then
          ettile = dcmplx(0d0)
      else if (im1.gt.0d0.and.i2.gt.0d0) then
          ettile = dcmplx(0d0,-2d0*pi)
      else if (im1.lt.0d0.and.i2.lt.0d0) then
          ettile = dcmplx(0d0, 2d0*pi)
      else
          ettile = dcmplx(0d0)
          if(im1.eq.0.and.dreal(c1).lt.0d0 .or.
     &       im2.eq.0.and.i2.eq.0d0.and.dreal(c1*c2).lt.0d0) then
             write(80,*) ' ettile not defined '
             write(80,*) ' ettile:  c1  = ',c1
             write(80,*) ' ettile:  c2  = ',c2
             write(80,*) ' ettile:  i1  = ',i1
             write(80,*) ' ettile:  i2  = ',i2
          end if
      end if
      end

************************************************************************
      function olieta(c1,c2)
************************************************************************
*     complex eta-function                                             *
*----------------------------------------------------------------------*
*     8.06.90    Ansgar Denner       last changed   11.07.94           *
************************************************************************
      implicit   none
      complex*16 olieta,c1,c2
      real*8     im1,im2,im12,re1,re2
      real*8     pi,pi2_6
      common /pietc/  pi,pi2_6
      integer    ltest
      common /ltest/ ltest
      complex*16 olietal,olietatest

      re1    = dreal(c1)
      re2    = dreal(c2)
      im1    = dimag(c1)
      im2    = dimag(c2)
      im12   = dimag(c1*c2)
 
      if(im1.lt.0d0.and.im2.lt.0d0.and.im12.gt.0d0) then
          olieta = dcmplx(0d0,2d0*pi)
      else if (im1.gt.0d0.and.im2.gt.0d0.and.im12.lt.0d0) then
          olieta = dcmplx(0d0,-2d0*pi)
      else
          olieta = dcmplx(0d0)
          if(.not.(im2.eq.0d0.and.re2.gt.0d0.or.
     &             im1.eq.0d0.and.re1.gt.0d0).and.
     &       (im1.eq.0.and.re1.lt.0d0 .or.
     &        im2.eq.0.and.re2.lt.0d0 .or.
     &        im12.eq.0.and.dreal(c1*c2).lt.0d0)) then
             write(80,*) ' olieta not defined '
             write(80,*) ' olieta:  c1  = ',c1
             write(80,*) ' olieta:  c2  = ',c2
             write(80,*) ' olieta:  c12 = ',c1*c2
          end if
      end if

      if (ltest.gt.0) then
        olietatest = olietal(c1,c2)
        if (abs(olietatest-olieta).gt.1d-10) then
          write(*,*) ' olieta''s inconsistent'
          write(*,*) ' olieta:  c1  = ',c1
          write(*,*) ' olieta:  c2  = ',c2
          write(*,*) ' olieta:  c12 = ',c1*c2
          write(*,*) ' olieta       = ',olieta
          write(*,*) ' olietal      = ',olietatest
        end if
      end if

      end

************************************************************************
      function olietal(c1,c2)
************************************************************************
*     complex eta-function via logarithms                              *
*----------------------------------------------------------------------*
*     21.05.04    Ansgar Denner       last changed   21.05.04          *
************************************************************************
      implicit   none
      complex*16 olietal,c1,c2

      olietal = cdlog(c1*c2) - cdlog(c1) - cdlog(c2)

      end


************************************************************************
      function olietae(c1,c2,i1,i2,i12)
************************************************************************
*     complex eta-function with explicit infinitesimal im parts        *
*----------------------------------------------------------------------*
*     8.06.90    Ansgar Denner       last changed   22.06.04           *
************************************************************************
      implicit   none
      complex*16 olietae,c1,c2
      real*8     i1,i2,i12
      real*8     im1,im2,im12,re1,re2,re12
      real*8     pi,pi2_6
      common /pietc/  pi,pi2_6
      integer    ltest
      common /ltest/ ltest
      complex*16 olietae2l,olietaetestl
      real*8     impacc
      common /impacc/ impacc
c      real*8     calacc
c      common /calacc/ calacc

c      write(*,*) 'olietae in',c1,c2,c1*c2,i1,i2,i12

      im1    = dimag(c1)
      im2    = dimag(c2)
      im12   = dimag(c1*c2)
      re1    = dreal(c1)
      re2    = dreal(c2)
      re12   = dreal(c1*c2)
c      write(*,*) 'olietae im',im1,im2,im12
c --> added 22.06.04
      if (abs(im1).lt.abs(re1)*impacc) im1 = 0d0
      if (abs(im2).lt.abs(re2)*impacc) im2 = 0d0
      if (abs(im12).lt.abs(re12)*impacc) im12 = 0d0
c <-- added 22.06.04
c      write(*,*) 'olietae im',im1,im2,im12

      if (im1.ne.0d0.or.im2.ne.0d0.or.im12.ne.0d0) then
        if(im1 .eq.0d0) im1  = i1
        if(im2 .eq.0d0) im2  = i2
        if(im12.eq.0d0) im12 = i12
        if(im1.lt.0d0.and.im2.lt.0d0.and.im12.gt.0d0) then
          olietae = dcmplx(0d0,2d0*pi)
        else if (im1.gt.0d0.and.im2.gt.0d0.and.im12.lt.0d0) then
          olietae = dcmplx(0d0,-2d0*pi)
        else
          olietae = dcmplx(0d0)
          if(im1.eq.0.and.re1.lt.0d0 .or.
     &         im2.eq.0.and.re2.lt.0d0 .or.
     &         im12.eq.0.and.re12.lt.0d0) then
            write(80,*) ' olietae not defined '
            write(80,*) ' olietae:  c1  = ',c1,i1
            write(80,*) ' olietae:  c2  = ',c2,i2
            write(80,*) ' olietae:  c12 = ',c1*c2,i12
          end if
        end if
      else 
        olietae = dcmplx(0d0)      
        if (re1.lt.0d0.and.re2.lt.0d0) then
          if(i1.gt.0d0.and.i2.gt.0d0) then
            olietae = dcmplx(0d0,-2d0*pi)
          else   if(i1.lt.0d0.and.i2.lt.0d0) then
            olietae = dcmplx(0d0,2d0*pi)
          end if
        else if (re12.lt.0d0.and.re1.lt.0d0) then
          if(i1.gt.0d0.and.i12.lt.0d0) then
            olietae = dcmplx(0d0,-2d0*pi)
          else   if(i1.lt.0d0.and.i12.gt.0d0) then
            olietae = dcmplx(0d0,2d0*pi)
          end if
        else if (re12.lt.0d0.and.re2.lt.0d0) then
          if(i2.gt.0d0.and.i12.lt.0d0) then
            olietae = dcmplx(0d0,-2d0*pi)
          else   if(i2.lt.0d0.and.i12.gt.0d0) then
            olietae = dcmplx(0d0,2d0*pi)
          end if
        end if
      end if

c      write(*,*) 'olietae = ',olietae

      if (ltest.gt.0) then
        olietaetestl = olietae2l(c1,c2,c1*c2,i1,i2,i12)
        if (abs(olietaetestl-olietae).gt.1d-10) then
          write(*,*) ' olietae''s inconsistent'
          write(*,*) ' olietae:  c1  = ',c1,i1
          write(*,*) ' olietae:  c2  = ',c2,i2
          write(*,*) ' olietae:  c12 = ',c1*c2,i12
          write(*,*) ' olietae       = ',olietae
          write(*,*) ' olietae2l     = ',olietaetestl
        end if
      end if

      end
************************************************************************
      function olietae2(c1,c2,c12,i1,i2,i12)
************************************************************************
*     complex eta-function with explicit infinitesimal im parts        *
*     c12 = c1*c2                                                      *
*----------------------------------------------------------------------*
*     25.06.04    Ansgar Denner       last changed   30.09.08          *
************************************************************************
      implicit   none
      complex*16 olietae2,c1,c2,c12
      real*8     i1,i2,i12
      real*8     im1,im2,im12,re1,re2,re12
      real*8     pi,pi2_6
      common /pietc/  pi,pi2_6
      integer    ltest
      common /ltest/ ltest
      complex*16 olietae2l,olietaetestl
      real*8     impacc
      common /impacc/ impacc
      real*8     calacc
      common /calacc/ calacc

c      write(*,*) 'olietae in',c1,c2,c1*c2,i1,i2,i12,calacc

      im1    = dimag(c1)
      im2    = dimag(c2)
      im12   = dimag(c12)
      re1    = dreal(c1)
      re2    = dreal(c2)
      re12   = dreal(c12)
c      write(*,*) 'olietae im',im1,im2,im12
c --> added 22.06.04
      if (abs(im1).lt.impacc*abs(re1)) im1 = 0d0
      if (abs(im2).lt.impacc*abs(re2)) im2 = 0d0
      if (abs(im12).lt.impacc*abs(re12)) im12 = 0d0
c <-- added 22.06.04
c      write(*,*) 'olietae im',im1,im2,im12

      if (im1.ne.0d0.or.im2.ne.0d0.or.im12.ne.0d0) then
        if(im1 .eq.0d0) im1  = i1
        if(im2 .eq.0d0) im2  = i2
        if(im12.eq.0d0) im12 = i12
        if(im1.lt.0d0.and.im2.lt.0d0.and.im12.gt.0d0) then
          olietae2 = dcmplx(0d0,2d0*pi)
        else if (im1.gt.0d0.and.im2.gt.0d0.and.im12.lt.0d0) then
          olietae2 = dcmplx(0d0,-2d0*pi)
        else
          olietae2 = dcmplx(0d0)
          if(im1.eq.0.and.re1.lt.0d0 .or.
     &         im2.eq.0.and.re2.lt.0d0 .or.
     &         im12.eq.0.and.re12.lt.0d0) then
            write(80,*) ' olietae2 not defined '
            write(80,*) ' olietae2:  c1  = ',c1,i1
            write(80,*) ' olietae2:  c2  = ',c2,i2
            write(80,*) ' olietae2:  c12 = ',c1*c2,i12
          end if
        end if
      else 
        olietae2 = dcmplx(0d0)      
        if (re1.lt.0d0.and.re2.lt.0d0) then
          if(i1.gt.0d0.and.i2.gt.0d0) then
            olietae2 = dcmplx(0d0,-2d0*pi)
          else   if(i1.lt.0d0.and.i2.lt.0d0) then
            olietae2 = dcmplx(0d0,2d0*pi)
          end if
        else if (re12.lt.0d0.and.re1.lt.0d0) then
          if(i1.gt.0d0.and.i12.lt.0d0) then
            olietae2 = dcmplx(0d0,-2d0*pi)
          else   if(i1.lt.0d0.and.i12.gt.0d0) then
            olietae2 = dcmplx(0d0,2d0*pi)
          end if
        else if (re12.lt.0d0.and.re2.lt.0d0) then
          if(i2.gt.0d0.and.i12.lt.0d0) then
            olietae2 = dcmplx(0d0,-2d0*pi)
          else   if(i2.lt.0d0.and.i12.gt.0d0) then
            olietae2 = dcmplx(0d0,2d0*pi)
          end if
        end if
      end if

      if (ltest.gt.0) then
        olietaetestl = olietae2l(c1,c2,c12,i1,i2,i12)
        if (abs(olietaetestl-olietae2).gt.1d-10) then
          write(*,*) ' olietae2''s inconsistent'
          write(*,*) ' olietae2: c1  = ',c1,i1,im1
          write(*,*) ' olietae2: c2  = ',c2,i2,im2
          write(*,*) ' olietae2: c12 = ',c1*c2,i12,im12
          write(*,*) ' olietae2: c12 = ',c12,i12,im12
          write(*,*) ' olietae2      = ',olietae2
          write(*,*) ' olietae2l     = ',olietaetestl
        end if
      end if

      end
************************************************************************
      function olietae2l(c1,c2,c12,i1,i2,i12)
************************************************************************
*     complex eta-function with explicit infinitesimal im parts        *
*     from logarithms                                                  *
*----------------------------------------------------------------------*
*     21.05.04    Ansgar Denner       last changed   28.06.04          *
************************************************************************
      implicit   none
      complex*16 olietae2l,c1,c2,c12
      real*8     i1,i2,i12
      complex*16 cdln

      olietae2l = cdln(c12,i12) - cdln(c1,i1) - cdln(c2,i2)

      end

************************************************************************
      subroutine crootse(p2,m12,m22,x1,x2,y1,y2,r)
************************************************************************
*     roots of quadratic equation                                      *
*     p2*x*x+(m22-m12-p2)*x+m12     = p2*(x-x1)*(x-x2)                 *
*                                     p2*(x-1+y1)*(x-1+y2)             *
*----------------------------------------------------------------------*
*     21.07.93 Ansgar Denner      last changed  22.10.04 ad            *
************************************************************************
      implicit   none
      complex*16 p2,m12,m22
      complex*16 x1,x2,y1,y2,r
      complex*16 q12,q22,m1,m2
 
c --->   case y1=y2=0  <-> cB0f(m2,m,0) makes it difficult!
      m1 = sqrt(m12)
      m2 = sqrt(m22)
      r=sqrt((p2-(m1+m2)**2)*(p2-(m1-m2)**2))
c     r=sqrt(p2*p2+m12*m12+m22*m22-2d0*(m12*m22+m12*p2+m22*p2))
      q12=m12-m22+p2
      q22=m22-m12+p2
      x1=(q12+r)/(2d0*p2)
      x2=(q12-r)/(2d0*p2)

c      write(20,*) 'croot: ', r,x1,x2
c      write(20,*) 'croot: ', p2,q12,q22

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
      end

************************************************************************
      function fpve(n,x,y,eps)
************************************************************************
*     passarino and veltmans f(n,x) function                           *
*     y = 1-x,  x  = x + i eps                                         *
*----------------------------------------------------------------------*
*     21.07.93 Ansgar Denner      last changed  16.03.04 ad            *
************************************************************************
      implicit   none
      complex*16 x,y,fpve,cdln
      real*8     calacc,eps
      integer    n,m
      common /calacc/ calacc
 
      if(abs(x).lt.10d0) then
        if (n.eq.0) then
          fpve = -cdln(-y/x,eps)
        elseif (x.ne.dcmplx(0d0)) then
          fpve = -x**n*cdln(-y/x,eps)
          do m=1,n
            fpve = fpve - x**(n-m)/m
          enddo
        else
          fpve = -1d0/n
        endif
      else
        fpve = 0d0
        do m=1,10
          fpve = fpve + 1d0/((m+n)*x**m)
          if (abs(1d0/(x**m*fpve)).lt.calacc**2) return
        enddo
      endif
      end

************************************************************************
      function yfpve(n,x,y,eps)
************************************************************************
*     passarino and veltmans f(n,x) function times (1-x) = y           *
*     x  = x + i eps                                                   *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  16.03.04 ad            *
************************************************************************
      implicit   none
      complex*16 x,y,fpve,yfpve
      real*8     eps 
      integer    n
 
      if(abs(y).eq.0d0) then
        yfpve = 0d0
      else
        yfpve = y*fpve(n,x,y,eps)
      endif
      end

************************************************************************
      function xyfpve(n,x,y,eps)
************************************************************************
*     passarino and veltmans f(n,x) function times x*(1-x) = x*y       *
*     x  = x + i eps                                                   *
*----------------------------------------------------------------------*
*     20.09.04 Ansgar Denner      last changed  20.09.04 ad            *
************************************************************************
      implicit   none
      complex*16 x,y,fpve,xyfpve
      real*8     eps 
      integer    n
 
      if(abs(y).eq.0d0.or.abs(x).eq.0d0) then
        xyfpve = 0d0
      else
        xyfpve = x*y*fpve(n,x,y,eps)
      endif
      end

************************************************************************
      function xlogxe(x,eps)
************************************************************************
*     x*log(x),    x  = x + i eps                                      *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  16.03.04 ad            *
************************************************************************
      implicit   none
      complex*16 x,xlogxe,cdln
      real*8     eps
 
      if(abs(x).eq.0d0) then
        xlogxe = 0d0
      else
        xlogxe = x*cdln(x,eps)
      endif
      end

************************************************************************
      subroutine elcminfiv(p2,m12,m22,q2,mm12,mm22,flag)
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
      complex*16 elimcminf2
      integer    flag
 
      p2s  = elimcminf2(p2)
      mm12s = elimcminf2(m12)
      mm22s = elimcminf2(m22)

      if(p2s.ne.0d0.or.mm12s.ne.0d0.or.mm22s.ne.0d0) then
        q2 = p2s
        mm12 = mm12s
        mm22 = mm22s
        flag = 0
      else
        q2 = p2
        mm12 = m12
        mm22 = m22
        flag = 1
      end if
      end

************************************************************************
      subroutine cerrB0(q2,m12,m22,value,name,code)
************************************************************************
*     error messages for B0 functions                                  *
*     q2 outer momentum      m12,m22 inner masses squared              *
*----------------------------------------------------------------------*
*     19.08.93 Ansgar Denner      last changed  12.05.04 ad            *
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




