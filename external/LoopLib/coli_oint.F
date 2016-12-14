************************************************************************
*                                                                      *
*     scalar functions: interface to old routines                      *
*                                                                      *
************************************************************************
*                                                                      *
*     last changed  28.10.08  Ansgar Denner                            *
*                                                                      *
************************************************************************
* subroutines:                                                         *
* defcoli,setmudim2,setmphoton2                                        *
* setcminf2,clearcminf2,etcminfscale2                                  *
* functions:                                                           *
* cA0f                                                                 *
* cB0f,  cB1f,   cB00f,  cB11f, cB001f, cB111f                         *
* cDB0f, cDB1f,  cDB00f                                                *
* cBnf,  cB00nf, cB0000nf                                              *
* cC0f,  cD0f                                                          *
* cscaledown, scaledown2, elimcminf, elimcminf2                        *
************************************************************************
      subroutine defcoli
************************************************************************
* fixing of default values for various variables                       *
*----------------------------------------------------------------------*
* 07.02.90  Ansgar Denner         last changed  07.06.04               *
************************************************************************
      implicit   none
      real*8     eps,calacc,impacc,pi,pi2_6
      integer    ltest,sym
 
      common /ltest/  ltest
      real*8     mudim2
      common /uv/     mudim2
      common /infeps/ eps
      common /impacc/ impacc
      common /calacc/ calacc
      common /pietc/  pi,pi2_6
      common /sym/    sym
      integer    ncminf
      common /ncminf/  ncminf
      complex*16 cminf(10)
      common /cminf/   cminf
      complex*16 cminf2(10)
      common /cminf2/  cminf2
      integer    nfw
      common /nfw/    nfw
      real*8     mfw(10)
      common /mfw/    mfw     
      real*8     gfw(10)
      common /gfw/    gfw     
  
      real*8     lambda,lambda2
      common /lambda/ lambda,lambda2
      real*8     minfscale,minfscale2,shiftms2
      common /cmsing/ minfscale,minfscale2,shiftms2

      integer  errout,testout
      common /out/ errout,testout

      real*8     testaccd0
      common /testaccd0/ testaccd0

      integer    init
      data init /0/
      save init

c      data pi /3.1415926535897932384626433832d0/
c     &     pi2_6 /1.6449340668482264364724151666d0/
c     &     eps /1d-20/
c     &     calacc /5d-16/
c     &     mudim2 /1d0/
c     &     lambda,lambda2 /1d0,1d0/
c     &     minfscale,minfscale2,shiftms2 /1d0,1d0,0d0/
c      data ltest,sym,nminf,nfw,lolifw /5*0/
     &    

c      write(*,*) pi,pi2_6,eps,calacc
c      write(*,*) mudim2,lambda,lambda2
c      write(*,*) minfscale,minfscale2,shiftms2
c      write(*,*)  ltest,sym,nminf,nfw,lolifw

      if (init.ne.0) return
      init = 1
C      print *,'======================================================='
C      print *,'  COLI, a package to evaluate one-loop integrals'
C      print *,'written by A. Denner, PSI, Villigen, Switzerland'
C      print *,'======================================================='
C     print *,'THIS CODE IS PRELIMINARY!!!'
C      print *,'There is no guarantee that the results are correct.'
C      print *,'This code has been used and tested for e+e- -> 4f.' 
C      print *,'Cases that were not required there' 
C      print *,'are either not implemented or not tested.'
C      print *,'======================================================='
C      print *,'If you use the code, you must refer to it'
C      print *,'and cite the papers it is based on.'
C      print *,'A. Denner and S. Dittmaier, Nucl.Phys.B734:62-115,2006' 
C      print *,'A. Denner and S. Dittmaier, Nucl.Phys.B658:175-202,2003'
C      print *,'  for the tensor integral reduction,'
C      print *,'A. Denner, U. Nierste, and R. Scharf,'
C      print *,'Nucl.Phys.B367:637-656,1991'
C      print *,'A. Denner and W. Beenakker, Nucl.Phys.B338:349-370,1990'
C      print *,'  for the scalar integrals.'
C      print *,'========================================================'



c mu^2 of dimensional regularization in uv divergent loop integrals
      call setmudim2(1d0)
c photon mass in ir divergent loop integrals
      call setmphoton2(1d0)
c scale factor for mass-singular masses
      call setcminfscale2(1d0)
c extra term for mass-singular squared logarithms
      shiftms2 = 0d0
c infinitesimal parameter
      eps    = 1d-20
c size of imaginary parts below which explicit ieps take over
      impacc = 1d-16
c  ltest=0 --> no testing
c  ltest>0 --> with testing
      ltest  = 0
c  accuracy for testing of 
      testaccd0 =1d-10

c should be set to the accuracy of the computer: n digits -> 1d-n
      calacc = 5d-16
c pi
      pi=3.1415926535897932384626433832d0
      pi2_6 = 1.6449340668482264364724151666d0

c sym = 1 calculation of symmetric tensor coefficients
c sym = 0 no calculation of symmetric tensor coefficients
      sym = 0

c ncminf = 0   no masses are neglected
      ncminf = 0

c errout:  output channel for error messages
      errout = 6
c testout: output channel for testing purposes
      testout = 9

      call initcoli

      end

************************************************************************
      subroutine setmudim2(mudim2in)
************************************************************************
*     sets UV-regulator mass squared to mudim2in                       *
*-----------------------------------------------------------------------
*     28.10.08 Ansgar Denner                                           *
************************************************************************
      implicit   none
      real*8     mudim2in
      real*8     mudim2

      common /uv/     mudim2

      mudim2 = mudim2in
C      write(*,*) 'UV regulator mass set to    mudim2 = ',mudim2

      call setmuuv2_coli(mudim2)

      end

************************************************************************
      subroutine setmphoton2(lambdain2)
************************************************************************
*     sets photon mass squared regulator to lambda2in                  *
*-----------------------------------------------------------------------
*     28.10.08 Ansgar Denner           last changed  28.10.08          *
************************************************************************
      implicit   none
      real*8     lambdain2
      real*8     lambda,lambda2
      common /lambda/ lambda,lambda2

      lambda2 = lambdain2
      lambda  = sqrt(lambda2)
C      write(*,*) 'photon mass squared set to lambda2 = ',lambda2

      call setmuir2_coli(lambdain2)

      end

************************************************************************
      subroutine setcminf2(m2)
************************************************************************
*     adds m2 to list of mass squared to be neglected                  *
*-----------------------------------------------------------------------
*     11.03.04 Ansgar Denner         last changed 04.05.04             *
************************************************************************
      implicit   none
      complex*16 m2
      integer    ncminf
      common /ncminf/     ncminf
      complex*16 cminf(10)
      common /cminf/     cminf    
      complex*16 cminf2(10)
      common /cminf2/    cminf2    
      complex*16 cminffix(10)
      common /cminffix/  cminffix    
      complex*16 cminffix2(10)
      common /cminffix2/ cminffix2    
      real*8     minfscale,minfscale2,shiftms2
      common /cmsing/ minfscale,minfscale2,shiftms2
      integer    i

c      data ncminf /0/

      do 10 i=1,ncminf
        if(m2.eq.cminffix2(i)) return
 10   continue
      
      ncminf = ncminf+1
      cminffix(ncminf) = sqrt(m2)
      cminffix2(ncminf) = m2
      cminf2(ncminf) = m2/minfscale2
      cminf(ncminf) = sqrt(cminf2(ncminf))

C      write(*,*) 'added to small masses squared:  m2 =  ',m2

      call setminf2_coli(m2)
    
      end



************************************************************************
      function cA0f(m2)
************************************************************************
*     scalar one-loop-integral cA0f                                    *
*-----------------------------------------------------------------------
*     28.10.08 Ansgar Denner      last changed 28.10.08                *
************************************************************************
      implicit   none
      real*8     mudim2
      complex*16 cA0f,m2
      complex*16 A0_coli
 
      cA0f = A0_coli(m2)

      end

************************************************************************
      function cB0f(p2,m12,m22)
************************************************************************
*     scalar one loop integral cB0f                                    *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     26.10.08 Ansgar Denner      last changed 28.10.08 ad             *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 cB0f,B0_coli

      cB0f = B0_coli(p2,m12,m22)

      end

************************************************************************
      function cB1f(p2,m12,m22)
************************************************************************
*     scalar one loop integral cB1f                                    *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     26.10.08 Ansgar Denner      last changed 28.10.08 ad             *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 cB1f,B1_coli

      cB1f = B1_coli(p2,m12,m22)

      end

************************************************************************
      function cB00f(p2,m12,m22)
************************************************************************
*     scalar one loop integral cB00f                                   *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     26.10.08 Ansgar Denner      last changed 28.10.08 ad             *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 cB00f,B00_coli

      cB00f = B00_coli(p2,m12,m22)

      end

************************************************************************
      function cB11f(p2,m12,m22)
************************************************************************
*     scalar one loop integral cB11f                                   *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     26.10.08 Ansgar Denner      last changed 28.10.08 ad             *
************************************************************************
       implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 cB11f,B11_coli

      cB11f = B11_coli(p2,m12,m22)

      end


************************************************************************
      function cB001f(p2,m12,m22)
************************************************************************
*     scalar one loop integral cB001f                                  *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     26.10.08 Ansgar Denner      last changed 28.10.08 ad             *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 cB001f,B001_coli

      cB001f = B001_coli(p2,m12,m22)

      end

************************************************************************
      function cB111f(p2,m12,m22)
************************************************************************
*     scalar one loop integral cB111f                                  *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     26.10.08 Ansgar Denner      last changed 28.10.08 ad             *
************************************************************************
       implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 cB111f,B111_coli

      cB111f = B111_coli(p2,m12,m22)

      end


************************************************************************
      function cDB0f(p2,m12,m22)
************************************************************************
*     scalar one loop integral cDB0f                                   *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     26.10.08 Ansgar Denner      last changed 28.10.08 ad             *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 cDB0f,DB0_coli

      cDB0f = DB0_coli(p2,m12,m22)

      end

************************************************************************
      function cDB1f(p2,m12,m22)
************************************************************************
*     scalar one loop integral cDB1f                                   *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     26.10.08 Ansgar Denner      last changed 28.10.08 ad             *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 cDB1f,DB1_coli

      cDB1f = DB1_coli(p2,m12,m22)

      end

************************************************************************
      function cDB00f(p2,m12,m22)
************************************************************************
*     scalar one loop integral cDB00f                                  *
*     p2 outer momentum, m12,m22 inner masses squared                  *
*----------------------------------------------------------------------*
*     26.10.08 Ansgar Denner      last changed 28.10.08 ad             *
************************************************************************
      implicit   none
      complex*16 p2
      complex*16 m12,m22
      complex*16 cDB00f,DB00_coli

      cDB00f = DB00_coli(p2,m12,m22)

      end

************************************************************************
      function cBnf(n,p2,m12,m22)
************************************************************************
*     tensor one loop integral cB001..1 of rank n                      *
*     p2 outer momentum, m12,m22 inner masses                          *
*----------------------------------------------------------------------*
*     30.10.08 Ansgar Denner      last changed 30.10.08 ad             *
************************************************************************
      implicit   none
      integer    n
      complex*16 p2
      complex*16 m12,m22
      complex*16 cBnf,Bn_coli

      cBnf = Bn_coli(n,p2,m12,m22)

      end

************************************************************************
      function cB00nf(n,p2,m12,m22)
************************************************************************
*     tensor one loop integral cB001..1 of rank n+2                    *
*     p2 outer momentum, m12,m22 inner masses                          *
*----------------------------------------------------------------------*
*     30.10.08 Ansgar Denner      last changed 30.10.08 ad             *
************************************************************************
      implicit   none
      integer    n
      complex*16 p2
      complex*16 m12,m22
      complex*16 cB00nf,B00n_coli

      cB00nf = B00n_coli(n,p2,m12,m22)

      end

************************************************************************
      function cB0000nf(n,p2,m12,m22)
************************************************************************
*     tensor one loop integral cB00001..1 of rank n+4                  *
*     p2 outer momentum, m12,m22 inner masses                          *
*----------------------------------------------------------------------*
*     30.10.08 Ansgar Denner      last changed 30.10.08 ad             *
************************************************************************
      implicit   none
      integer    n
      complex*16 p2
      complex*16 m12,m22
      complex*16 cB0000nf,B0000n_coli

      cB0000nf = B0000n_coli(n,p2,m12,m22)

      end

************************************************************************
      function cC0f(p12,p23,p13,m12,m22,m32)
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
*----------------------------------------------------------------------*
*     28.10.08 Ansgar Denner        last changed  28.10.08             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p13
      complex*16 m12,m22,m32
      complex*16 cC0f,C0_coli

      cC0f = C0_coli(p12,p23,p13,m12,m22,m32)

      end

************************************************************************
      function cD0f(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)
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
*     01.03.04 Ansgar Denner        last changed  16.02.07             *
************************************************************************
      implicit   none
      complex*16 p12,p23,p34,p14,p13,p24
      complex*16 m12,m22,m32,m42
      complex*16 cD0f,D0_coli

      cD0f = D0_coli(p12,p23,p34,p14,p13,p24,m12,m22,m32,m42)

      end



************************************************************************
      subroutine clearcminf2()
************************************************************************
*     clears list of masses squared to be neglected                    *
*-----------------------------------------------------------------------
*     28.06.04 Ansgar Denner         last changed 28.06.04             *
************************************************************************
      implicit   none
      integer    ncminf
      common /ncminf/     ncminf
      complex*16 cminf(10)
      common /cminf/     cminf    
      complex*16 cminf2(10)
      common /cminf2/    cminf2    
      complex*16 cminffix(10)
      common /cminffix/  cminffix    
      complex*16 cminffix2(10)
      common /cminffix2/ cminffix2    
      real*8     minfscale,minfscale2,shiftms2
      common /cmsing/ minfscale,minfscale2,shiftms2
      integer    i

c      data ncminf /0/

      do 10 i=1,ncminf
        cminffix2(i) = 0d0
        cminffix(i) = 0d0
        cminf2(i) = 0d0
        cminf(i) = 0d0
 10   continue
      
      ncminf = 0
 
      write(*,*) 'list of complex small masses cleared'
    
      end

************************************************************************
      subroutine setcminfscale2(scale2)
************************************************************************
*     sets factor with with which small masses are scaled up in        *
*     mass-singular logarithms and rescales masses in /cminf/ /cminf2/ *
*-----------------------------------------------------------------------
*     19.03.04 Ansgar Denner         last changed  04.05.04            *
************************************************************************
      implicit   none
      real*8     scale2
      integer    ncminf
      common /ncminf/    ncminf
      complex*16 cminf(10)
      common /cminf/     cminf    
      complex*16 cminf2(10)
      common /cminf2/    cminf2    
      complex*16 cminffix(10)
      common /cminffix/  cminffix    
      complex*16 cminffix2(10)
      common /cminffix2/ cminffix2    
      real*8     minfscale,minfscale2,shiftms2
      common /cmsing/ minfscale,minfscale2,shiftms2

      integer    i
      minfscale2 = scale2
      minfscale = sqrt(scale2)

      do 10 i=1,ncminf
        cminf(i) = cminffix(i)/minfscale
        cminf2(i) = cminffix2(i)/minfscale2
 10   continue

C      write(*,*) 'small masses scaled with minfscale = ',minfscale
C      write(*,*) 'in mass-singular logarithms'

      end

************************************************************************
      function cscaledown(m)
************************************************************************
*     scales mass down by minfscale if m in /cminf/                     *
*-----------------------------------------------------------------------
*     22.03.04 Ansgar Denner                                           *
************************************************************************
      implicit   none
      complex*16 cscaledown,m
      integer    ncminf
      common /ncminf/     ncminf
      complex*16 cminf(10)
      common /cminf/     cminf    
      complex*16 cminffix(10)
      common /cminffix/  cminffix    

      integer    i

      cscaledown = m
      do 10 i=1,ncminf
         if(m.eq.cminffix(i)) then
           cscaledown = cminf(i)
           return
         end if 
 10   continue

      end

************************************************************************
      function cscaledown2(m2)
************************************************************************
*     scales mass down by minfscale^2 if m2 in /cminf2/                *
*-----------------------------------------------------------------------
*     22.03.04 Ansgar Denner                                           *
************************************************************************
      implicit   none
      complex*16 cscaledown2,m2
      integer    ncminf
      common /ncminf/     ncminf
      complex*16 cminf2(10)
      common /cminf2/    cminf2    
      complex*16 cminffix2(10)
      common /cminffix2/ cminffix2    

      integer    i
      cscaledown2 = m2
      do 10 i=1,ncminf
         if(m2.eq.cminffix2(i)) then
           cscaledown2 = cminf2(i)
           return
         end if 
 10   continue

      end

************************************************************************
      function elimcminf(m)
************************************************************************
*     if m = cminfs(i)  then  elimcminf = 0  else  elimcminf  = m      *
*-----------------------------------------------------------------------
*     11.02.04 Ansgar Denner                                           *
************************************************************************
      implicit   none
      complex*16 m,elimcminf

      integer    ncminf
      common /ncminf/     ncminf
      complex*16 cminf(10)
      common /cminf/     cminf 
      integer    i

c      write(*,*) 'elimcminf',m 
c      write(*,*) 'elimcminf',cminf
      elimcminf = m
      do 10 i=1,ncminf
         if(m.eq.cminf(i)) then
           elimcminf = 0d0
           return
         end if 
 10   continue
c      write(*,*) 'elimcminf',elimcminf,m-cminf(1)
      end

************************************************************************
      function elimcminf2(m2)
************************************************************************
*     if m2 = cminf2(i)  then  elimcminf2 = 0  else  elimcminf2  = m2  *
*-----------------------------------------------------------------------
*     11.02.04 Ansgar Denner          last changed  04.05.04           *
************************************************************************
      implicit   none
      complex*16 m2,elimcminf2

      integer    ncminf
      common /ncminf/     ncminf
      complex*16 cminf2(10)
      common /cminf2/    cminf2    
      integer    i

      elimcminf2 = m2
      do 10 i=1,ncminf           
         if(m2.eq.cminf2(i)) then
           elimcminf2 = 0d0
           return
         end if 
 10   continue
      end


