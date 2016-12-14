************************************************************************
*                                                                      *
*     interface to Stefans and Markus scalar integrals                 *
*                                                                      * 
*     last changed 01.07.04 Ansgar Denner                              *
*                                                                      *
************************************************************************
      subroutine initstdints(m12,m22,m32,m42,m52)
************************************************************************
* initialiation of Stefans loop integrals                              *
*----------------------------------------------------------------------*
*     30.06.04  Ansgar Denner         last changed  30.06.04 ad        *
************************************************************************
      implicit none
      complex*16 m12,m22,m32,m42,m52
      complex*16 cme2,cmf2(5)
      real*8     mf2(5),rmf2(5),me2,rme2
      real*8 lambda,lambda2,lambdasd,lambda2sd
      integer    i
c      integer i,n,gen12,gen34,dir(7),zdir(7)

c      common/fstate/qf(4),mf(4),nc4f,dir
      common/fmass/me2,rme2,cme2,mf2,rmf2,cmf2
      common/lambda/lambda,lambda2
      common/ir/lambdasd,lambda2sd

c ir regulators
      lambdasd = lambda
      lambda2sd = lambda2


c*** settings for full e+e- --> 4f reaction
c physical masses
      mf2(1) = m12
      mf2(2) = m22
      mf2(3) = m32
      mf2(4) = m42
      mf2(5) = m52
      me2    = m52

c fermion identifiers
      do i=1,5
        rmf2(i) = i*1d-20
        cmf2(i) = rmf2(i)
      enddo
      rme2 = rmf2(5)
      cme2 = rme2

      end

************************************************************************
      function getmass2std(m2)
************************************************************************
*     evaluation of masses for Stefans loop integrals                  *
*     if m2 = minf2(i)  then  getmass2std = mid  else  getmass2std = m2*
*----------------------------------------------------------------------*
*     30.06.04  Ansgar Denner         last changed  30.06.04 ad        *
************************************************************************
      implicit   none
      complex*16 m2,getmass2std

      integer    ncminf
      common /ncminf/     ncminf
      complex*16 cminf2(10)
      common /cminf2/    cminf2    
      complex*16 cminffix2(10)
      common /cminffix2/ cminffix2    
      complex*16 cme2,cmf2(5)
      real*8     mf2(5),rmf2(5),me2,rme2
      common/fmass/me2,rme2,cme2,mf2,rmf2,cmf2
      integer    i,j

c      write(*,*) 'getmass2std in  ',m2
c      write(*,*) 'getmass2std inf ',cminf2
c      write(*,*) 'getmass2std std ',mf2
      getmass2std = m2
      do 10 i=1,ncminf           
        if(m2.eq.cminf2(i)) then
c          write(*,*) 'i',i
          do 20 j=1,5
            if(dreal(cminffix2(i)).eq.mf2(j)) then
c          write(*,*) 'j',j
              getmass2std = j*1d-20
              return
            end if  
 20       continue
        end if 
 10   continue
      end

************************************************************************
      function cB0sd(p10,m02,m12)
************************************************************************
*	interface to Stefans                                           *
*	scalar 2-point function with complex masses                    *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cB0sd,xB0_
      complex*16 p10,m02,m12
      complex*16 getmass2std,mm02,mm12
      real*8     p1  
       
c      write(*,*) 'cB0sd in ',p10,m02,m12
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      p1   = dreal(getmass2std(p10))
c      write(*,*) 'cB0sd out',p1,mm02,mm12
      cB0sd = xB0_(p1,mm02,mm12,1)
      
      end

************************************************************************
      function cDB0sd(p10,m02,m12)
************************************************************************
*	interface to Stefans                                           *
*	momentum derivative of scalar 2-point function                 *
*       with complex masses                                            *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cDB0sd,xDB0_
      complex*16 p10,m02,m12
      complex*16 getmass2std,mm02,mm12,elimcminf2
      real*8     p1  
      real*8     minfscale,minfscale2,shiftms2
      common /cmsing/ minfscale,minfscale2,shiftms2
       
c      write(*,*) 'cDB0sd in ',p10,m02,m12
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      p1   = dreal(getmass2std(p10))
c      write(*,*) 'cDB0sd out',p1,mm02,mm12
       
      cDB0sd = xDB0_(p1,mm02,mm12,1)
      if (elimcminf2(p10).eq.0d0.and.elimcminf2(m02).eq.0d0.and.
     &    elimcminf2(m12).eq.0d0) then
        cDB0sd = cDB0sd*minfscale2
      end if      
      end

************************************************************************
      function cB1sd(p10,m02,m12)
************************************************************************
*	interface to Stefans                                           *
*	coefficient for vector 2-point function with complex masses    *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cB1sd,xB1_
      complex*16 p10,m02,m12
      complex*16 getmass2std,mm02,mm12
      real*8     p1  
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      p1   = dreal(getmass2std(p10))
       
      cB1sd = xB1_(p1,mm02,mm12,1)
      
      end

************************************************************************
      function cDB1sd(p10,m02,m12)
************************************************************************
*	interface to Stefans                                           *
*	momentum derivative of coefficient for vector 2-point function *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cDB1sd,xDB1_
      complex*16 p10,m02,m12
      complex*16 getmass2std,mm02,mm12,elimcminf2
      real*8     minfscale,minfscale2,shiftms2
      common /cmsing/ minfscale,minfscale2,shiftms2
      real*8     p1  
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      p1   = dreal(getmass2std(p10))
       
      cDB1sd = xDB1_(p1,mm02,mm12,1)
      if (elimcminf2(p10).eq.0d0.and.elimcminf2(m02).eq.0d0.and.
     &    elimcminf2(m12).eq.0d0) then
        cDB1sd = cDB1sd*minfscale2
      end if            
      end

************************************************************************
      function cB00sd(p10,m02,m12)
************************************************************************
*	interface to Stefans                                           *
*	tensor coefficient B00 for 2nd-rank tensor 2-point function    *
*	with complex masses                                            *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cB00sd,xB20_
      complex*16 p10,m02,m12
      complex*16 getmass2std,mm02,mm12
      real*8     p1  
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      p1   = dreal(getmass2std(p10))
       
      cB00sd = xB20_(p1,mm02,mm12,1)
      
      end

************************************************************************
      function cB11sd(p10,m02,m12)
************************************************************************
*	interface to Stefans                                           *
*	tensor coefficient B11 for 2nd-rank tensor 2-point function    *
*	with complex masses                                            *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cB11sd,xB21_
      complex*16 p10,m02,m12
      complex*16 getmass2std,mm02,mm12
      real*8     p1  
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      p1   = dreal(getmass2std(p10))
       
      cB11sd = xB21_(p1,mm02,mm12,1)
      
      end

************************************************************************
      function cB001sd(p10,m02,m12)
************************************************************************
*	interface to Stefans                                           *
*	tensor coefficient B001 for 3rd-rank tensor 2-point function   *
*	with complex masses                                            *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cB001sd,xB001_
      complex*16 p10,m02,m12
      complex*16 getmass2std,mm02,mm12
      real*8     p1  
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      p1   = dreal(getmass2std(p10))
       
      cB001sd = xB001_(p1,mm02,mm12,1)
      
      end

************************************************************************
      function cB111sd(p10,m02,m12)
************************************************************************
*	interface to Stefans                                           *
*	tensor coefficient B111 for 3rd-rank tensor 2-point function   *
*	with complex masses                                            *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cB111sd,xB111_
      complex*16 p10,m02,m12
      complex*16 getmass2std,mm02,mm12
      real*8     p1  
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      p1   = dreal(getmass2std(p10))
       
      cB111sd = xB111_(p1,mm02,mm12,1)
      
      end

c>************************************************************************
c>      subroutine cB0123sd(p10,m02,m12,B0,B1,B2,B3,rank)
c>************************************************************************
c>*	interface to Stefans                                           *
c>*	2-POINT TENSOR COEFFICIENTS with complex masses                *
c>*----------------------------------------------------------------------*
c>*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
c>************************************************************************
c>      implicit   none
c>      integer    rank
c>      complex*16 p10,m02,m12
c>      complex*16 B0sd,B1sd(1),B2sd(0:2,0:2),B3sd(0:2,0:2,2)
c>      complex*16 B0,B1,B2(0:1,0:1),B3(0:1,0:1,0:1)
c>      complex*16 getmass2std,mm02,mm12
c>      real*8     p1  
c>       
c>      mm02 = getmass2std(m02)
c>      mm12 = getmass2std(m12)
c>      p1   = dreal(getmass2std(p10))
c>       
c>      call  xB0123_(p1,mm02,mm12,B0sd,B1sd,B2sd,B3sd,rank,1)
c>      
c>      B0      =  B0sd
c>      if (rank.eq.0) return
c>      B1      =  B1sd(1)
c>      if (rank.eq.1) return
c>      B2(0,0) = B2sd(0,0)
c>      B2(1,1) = B2sd(1,1)
c>      if (rank.eq.2) return
c>      B3(0,0,1) = B3sd(0,0,1)
c>      B3(1,1,1) = B3sd(1,1,1)
c>
c>      end
c>
c>************************************************************************
c>      subroutine cB01234sd(p10,m02,m12,B0,B1,B2,B3,B4,rank)
c>************************************************************************
c>*	interface to Stefans                                           *
c>*	2-POINT TENSOR COEFFICIENTS with complex masses                *
c>*----------------------------------------------------------------------*
c>*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
c>************************************************************************
c>      implicit   none
c>      integer    rank
c>      complex*16 p10,m02,m12
c>      complex*16 B0sd,B1sd(1),B2sd(0:2,0:2),B3sd(0:2,0:2,2)
c>      complex*16 B4sd(0:2,0:2,0:2,0:2)
c>      complex*16 B0,B1,B2(0:1,0:1),B3(0:1,0:1,0:1)
c>      complex*16 B4(0:1,0:1,0:1,0:1)
c>      complex*16 getmass2std,mm02,mm12
c>      real*8     p1  
c>       
c>      mm02 = getmass2std(m02)
c>      mm12 = getmass2std(m12)
c>      p1   = dreal(getmass2std(p10))
c>       
c>      call  xB01234_(p1,mm02,mm12,B0sd,B1sd,B2sd,B3sd,B4sd,rank,1)
c>      
c>      B0      =  B0sd
c>      if (rank.eq.0) return
c>      B1      =  B1sd(1)
c>      if (rank.eq.1) return
c>      B2(0,0) = B2sd(0,0)
c>      B2(1,1) = B2sd(1,1)
c>      if (rank.eq.2) return
c>      B3(0,0,1) = B3sd(0,0,1)
c>      B3(1,1,1) = B3sd(1,1,1)
c>      if (rank.eq.3) return
c>      B4(0,0,0,0) = B4sd(0,0,0,0)
c>      B4(0,0,1,1) = B4sd(0,0,1,1)
c>      B4(1,1,1,1) = B4sd(1,1,1,1)
c>
c>      end

************************************************************************
      subroutine cB012345sd(p10,m02,m12,B0,B1,B2,B3,B4,B5,rank)
************************************************************************
*	interface to Stefans                                           *
*	2-POINT TENSOR COEFFICIENTS with complex masses                *
*----------------------------------------------------------------------*
*     25.11.04  Ansgar Denner         last changed  25.11.04 ad        *
************************************************************************
      implicit   none
      integer    rank
      complex*16 p10,m02,m12
      complex*16 B0sd,B1sd(1),B2sd(0:2,0:2),B3sd(0:2,0:2,2)
      complex*16 B4sd(0:2,0:2,0:2,0:2),B5sd(0:2,0:2,0:2,0:2,2)
      complex*16 B0,B1,B2(0:1,0:1),B3(0:1,0:1,0:1)
      complex*16 B4(0:1,0:1,0:1,0:1),B5(0:1,0:1,0:1,0:1,0:1)
      complex*16 getmass2std,mm02,mm12
      real*8     p1  
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      p1   = dreal(getmass2std(p10))
       
      call  xB012345_(p1,mm02,mm12,B0sd,B1sd,B2sd,B3sd,B4sd,B5sd,rank,1)
      
      B0      =  B0sd
      if (rank.eq.0) return
      B1      =  B1sd(1)
      if (rank.eq.1) return
      B2(0,0) = B2sd(0,0)
      B2(1,1) = B2sd(1,1)
      if (rank.eq.2) return
      B3(0,0,1) = B3sd(0,0,1)
      B3(1,1,1) = B3sd(1,1,1)
      if (rank.eq.3) return
      B4(0,0,0,0) = B4sd(0,0,0,0)
      B4(0,0,1,1) = B4sd(0,0,1,1)
      B4(1,1,1,1) = B4sd(1,1,1,1)
      if (rank.eq.3) return
      B5(0,0,0,0,1) = B5sd(0,0,0,0,1)
      B5(0,0,1,1,1) = B5sd(0,0,1,1,1)
      B5(1,1,1,1,1) = B5sd(1,1,1,1,1)

      end

************************************************************************
      function cC0sd(p10,p21,p20,m02,m12,m22)
************************************************************************
*	interface to Stefans                                           *
*	scalar 3-point function with complex masses                    *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cC0sd,xC0_
      complex*16 p10,p21,p20,m02,m12,m22
      real*8     p1,p2,p3
      complex*16 getmass2std,mm02,mm12,mm22
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      mm22 = getmass2std(m22)
      p1 = dreal(getmass2std(p10))      
      p2 = dreal(getmass2std(p21))      
      p3 = dreal(getmass2std(p20))     
       
c      write(*,*) 'cCstd m0 ',mm02,m02
c      write(*,*) 'cCstd m1 ',mm12,m22
c      write(*,*) 'cCstd m2 ',mm22,m22
c      write(*,*) 'cCstd p1 ',p1,p10
c      write(*,*) 'cCstd p2 ',p2,p21
c      write(*,*) 'cCstd p3 ',p3,p20

      cC0sd = xC0_(p1,p2,p3,mm02,mm12,mm22,1)
      
      end

************************************************************************
      subroutine cCp1234sd(p10,p21,p20,m02,m12,m22,
     &                C0,C1,C2,C3,C4,rank)
************************************************************************
*	interface to Stefans                                           *
*	2-POINT TENSOR COEFFICIENTS with complex masses                *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      integer    rank
      complex*16 p10,p21,p20,m02,m12,m22
      complex*16 C0,C1(2),C2(0:2,0:2),C3(0:2,0:2,0:2),
     &           C4(0:2,0:2,0:2,0:2)
      complex*16 B0sd,B1sd(1),B2sd(0:2,0:2),B3sd(0:2,0:2,2)
      complex*16 C0sd,C1sd(2),C2sd(0:2,0:2),C3sd(0:2,0:2,2),
     &           C4sd(0:2,0:2,0:2,0:2)
      real*8     p1,p2,p3
      complex*16 getmass2std,mm02,mm12,mm22

      integer    sym       
      common /sym/    sym

c      write(*,*) 'cCp1234sd ',p10,p21,p20,m02,m12,m22,rank

      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      mm22 = getmass2std(m22)
      p1 = dreal(getmass2std(p10))      
      p2 = dreal(getmass2std(p21))      
      p3 = dreal(getmass2std(p20))     

      call  xC01234_(p1,p2,p3,mm02,mm12,mm22,
     &    C0sd,C1sd,C2sd,C3sd,C4sd,B0sd,B1sd,B2sd,B3sd,rank,1)
       
      C0      =  C0sd
      if (rank.eq.0) return
      C1(1)    =  C1sd(1)
      C1(2)    =  C1sd(2)
      if (rank.eq.1) return
      C2(0,0) = C2sd(0,0)
      C2(1,1) = C2sd(1,1)
      C2(1,2) = C2sd(1,2)
      C2(2,2) = C2sd(2,2)

      if (sym.gt.0) then
        C2(2,1) =  C2sd(2,1)
      end if

      if (rank.eq.2) return
      C3(0,0,1) = C3sd(0,0,1)
      C3(0,0,2) = C3sd(0,0,2)
      C3(1,1,1) = C3sd(1,1,1)
      C3(1,1,2) = C3sd(1,1,2)
      C3(1,2,2) = C3sd(1,2,2)
      C3(2,2,2) = C3sd(2,2,2)

      if (sym.gt.0) then
        C3(1,2,1) = C3sd(1,2,1)
        C3(2,1,1) = C3sd(2,1,1)
        C3(2,1,2) = C3sd(2,1,2)
        C3(2,2,1) = C3sd(2,2,1)
      end if

      if (rank.eq.3) return
      C4(0,0,0,0) = C4sd(0,0,0,0)
      C4(0,0,1,1) = C4sd(0,0,1,1)
      C4(0,0,1,2) = C4sd(0,0,1,2)
      C4(0,0,2,2) = C4sd(0,0,2,2)
      C4(1,1,1,1) = C4sd(1,1,1,1)
      C4(1,1,1,2) = C4sd(1,1,1,2)
      C4(1,1,2,2) = C4sd(1,1,2,2)
      C4(1,2,2,2) = C4sd(1,2,2,2)
      C4(2,2,2,2) = C4sd(2,2,2,2)

      if (sym.gt.0) then
        C4(0,0,2,1) = C4sd(0,0,2,1)
        C4(2,1,1,1) = C4sd(2,1,1,1)
        C4(1,2,1,1) = C4sd(1,2,1,1)
        C4(1,1,2,1) = C4sd(1,1,2,1)
        C4(1,2,1,2) = C4sd(1,2,1,2)
        C4(2,1,1,2) = C4sd(2,1,1,2)
        C4(1,2,2,1) = C4sd(1,2,2,1)
        C4(2,1,2,1) = C4sd(2,1,2,1)
        C4(2,2,1,1) = C4sd(2,2,1,1)
        C4(2,1,2,2) = C4sd(2,1,2,2)
        C4(2,2,1,2) = C4sd(2,2,1,2)
        C4(2,2,2,1) = C4sd(2,2,2,1)
      end if
      end

************************************************************************
      function cD0sd(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32)
************************************************************************
*	interface to Stefans                                           *
*	scalar 3-point function with complex masses                    *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cD0sd,xD0_
      complex*16 p10,p21,p32,p30,p20,p31,m02,m12,m22,m32
      real*8     p1,p2,p3,p4,s12,s23
      complex*16 getmass2std,mm02,mm12,mm22,mm32
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      mm22 = getmass2std(m22)
      mm32 = getmass2std(m32)
      p1 = dreal(getmass2std(p10))      
      p2 = dreal(getmass2std(p21))      
      p3 = dreal(getmass2std(p32))     
      p4 = dreal(getmass2std(p30))
      s12 = dreal(getmass2std(p20))
      s23 = dreal(getmass2std(p31))

c      write(90,*) 'cD0sd ',p10,p21,p32,p30,p20,p31,m02,m12,m22,m32
c      write(90,*) 'cD0sd ',p1,p2,p3,p4,s12,s23,mm02,mm12,mm22,mm32

c note same order of arguments s12,s23
      cD0sd = xD0_(p1,p2,p3,p4,s12,s23,mm02,mm12,mm22,mm32,1)
      
      end

c>************************************************************************
c>      subroutine cDp1234sd(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
c>     &                D0,D1,D2,D3,D4,rank)
c>************************************************************************
c>*	interface to Stefans                                           *
c>*	2-POINT TENSOR COEFFICIENTS with complex masses                *
c>*----------------------------------------------------------------------*
c>*     01.07.04  Ansgar Denner         last changed  12.10.04 ad        *
c>************************************************************************
c>      implicit   none
c>      integer    rank
c>      complex*16 p10,p21,p32,p30,p20,p31,m02,m12,m22,m32
c>      complex*16 D0,D1(3),D2(0:3,0:3),D3(0:3,0:3,0:3),
c>     &           D4(0:3,0:3,0:3,0:3)
c>      complex*16 C0sd,C1sd(2),C2sd(0:2,0:2),C3sd(0:2,0:2,2)
c>      complex*16 D0sd,D1sd(3),D2sd(0:3,0:3),D3sd(0:3,0:3,3),
c>     &           D4sd(0:3,0:3,0:3,0:3)
c>      real*8     p1,p2,p3,p4,s12,s23   
c>      complex*16 getmass2std,mm02,mm12,mm22,mm32
c>
c>      integer    i,j,k,l       
c>      integer    sym       
c>      common /sym/    sym
c>
c>      mm02 = getmass2std(m02)
c>      mm12 = getmass2std(m12)
c>      mm22 = getmass2std(m22)
c>      mm32 = getmass2std(m32)
c>      p1 = dreal(getmass2std(p10))      
c>      p2 = dreal(getmass2std(p21))      
c>      p3 = dreal(getmass2std(p32))     
c>      p4 = dreal(getmass2std(p30))
c>      s12 = dreal(getmass2std(p20))
c>      s23 = dreal(getmass2std(p31))
c>
c>c note different order of arguments s23,s12
c>      call  xD01234_(p1,p2,p3,p4,s23,s12,mm02,mm12,mm22,mm32,
c>     &    D0sd,D1sd,D2sd,D3sd,D4sd,C0sd,C1sd,C2sd,C3sd,rank,1)
c>       
c>      D0      =  D0sd
c>      if (rank.eq.0) return
c>      D1(1)    =  D1sd(1)
c>      D1(2)    =  D1sd(2)
c>      D1(3)    =  D1sd(3)
c>      if (rank.eq.1) return
c>      D2(0,0) = D2sd(0,0)
c>      D2(1,1) = D2sd(1,1)
c>      D2(1,3) = D2sd(1,3)
c>      D2(1,2) = D2sd(1,2)
c>      D2(2,2) = D2sd(2,2)
c>      D2(2,3) = D2sd(2,3)
c>      D2(3,3) = D2sd(3,3)
c>
c>      if (sym.gt.0) then
c>        D2(3,1) = D2sd(3,1)
c>        D2(2,1) = D2sd(2,1)
c>        D2(3,2) = D2sd(3,2)
c>      end if
c>
c>      if (rank.eq.2) return
c>      D3(0,0,1) = D3sd(0,0,1)
c>      D3(0,0,2) = D3sd(0,0,2)
c>      D3(0,0,3) = D3sd(0,0,3)
c>      D3(1,1,1) = D3sd(1,1,1)
c>      D3(1,1,2) = D3sd(1,1,2)
c>      D3(1,1,3) = D3sd(1,1,3)
c>      D3(1,2,2) = D3sd(1,2,2)
c>      D3(1,2,3) = D3sd(1,2,3)
c>      D3(1,3,3) = D3sd(1,3,3)
c>      D3(2,2,2) = D3sd(2,2,2)
c>      D3(2,2,3) = D3sd(2,2,3)
c>      D3(2,3,3) = D3sd(2,3,3)
c>      D3(3,3,3) = D3sd(3,3,3)
c>
c>      if (sym.gt.0) then
c>        D3(1,2,1) = D3sd(1,2,1)
c>        D3(2,1,1) = D3sd(2,1,1)
c>        D3(1,3,1) = D3sd(1,3,1)
c>        D3(3,1,1) = D3sd(3,1,1)
c>        D3(2,1,2) = D3sd(2,1,2)
c>        D3(2,2,1) = D3sd(2,2,1)
c>        D3(3,3,1) = D3sd(3,3,1)
c>        D3(3,1,3) = D3sd(3,1,3)
c>        D3(2,3,2) = D3sd(2,3,2)
c>        D3(3,2,2) = D3sd(3,2,2)
c>        D3(3,2,3) = D3sd(3,2,3)
c>        D3(3,3,2) = D3sd(3,3,2)
c>        D3(1,3,2) = D3sd(1,3,2)
c>        D3(2,1,3) = D3sd(2,1,3)
c>        D3(2,3,1) = D3sd(2,3,1)
c>        D3(3,2,1) = D3sd(3,2,1)
c>        D3(3,1,2) = D3sd(3,1,2)
c>      end if
c>
c>      if (rank.eq.3) return
c>      D4(0,0,0,0) = D4sd(0,0,0,0)
c>      D4(0,0,1,1) = D4sd(0,0,1,1)
c>      D4(0,0,1,2) = D4sd(0,0,1,2)
c>      D4(0,0,1,3) = D4sd(0,0,1,3)
c>      D4(0,0,2,2) = D4sd(0,0,2,2)
c>      D4(0,0,2,3) = D4sd(0,0,2,3)
c>      D4(0,0,3,3) = D4sd(0,0,3,3)
c>      D4(1,1,1,1) = D4sd(1,1,1,1)
c>      D4(1,1,1,2) = D4sd(1,1,1,2)
c>      D4(1,1,1,3) = D4sd(1,1,1,3)
c>      D4(1,1,2,2) = D4sd(1,1,2,2)
c>      D4(1,1,2,3) = D4sd(1,1,2,3)
c>      D4(1,1,3,3) = D4sd(1,1,3,3)
c>      D4(1,2,2,2) = D4sd(1,2,2,2)
c>      D4(1,2,2,3) = D4sd(1,2,2,3)
c>      D4(1,2,3,3) = D4sd(1,2,3,3)
c>      D4(1,3,3,3) = D4sd(1,3,3,3)
c>      D4(2,2,2,2) = D4sd(2,2,2,2)
c>      D4(2,2,2,3) = D4sd(2,2,2,3)
c>      D4(2,2,3,3) = D4sd(2,2,3,3)
c>      D4(2,3,3,3) = D4sd(2,3,3,3)
c>      D4(3,3,3,3) = D4sd(3,3,3,3)
c>      if (sym.gt.0) then
c>        do i=1,3
c>        do j=1,3
c>          D4(0,0,i,j) = D4sd(0,0,i,j)
c>        do k=1,3
c>        do l=1,3
c>          D4(i,j,k,l) = D4sd(i,j,k,l)
c>        end do
c>        end do
c>        end do
c>        end do
c>      end if
c>      end

************************************************************************
      subroutine cDp12345sd(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,
     &                D0,D1,D2,D3,D4,D5,rank)
************************************************************************
*	interface to Stefans                                           *
*	2-POINT TENSOR COEFFICIENTS with complex masses                *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  12.10.04 ad        *
************************************************************************
      implicit   none
      integer    rank
      complex*16 p10,p21,p32,p30,p20,p31,m02,m12,m22,m32
      complex*16 D0,D1(3),D2(0:3,0:3),D3(0:3,0:3,0:3),
     &           D4(0:3,0:3,0:3,0:3),D5(0:3,0:3,0:3,0:3,0:3)
      complex*16 C0sd,C1sd(2),C2sd(0:2,0:2),C3sd(0:2,0:2,2)
     &           ,C4sd(0:2,0:2,0:2,0:2)
      complex*16 D0sd,D1sd(3),D2sd(0:3,0:3),D3sd(0:3,0:3,3),
     &           D4sd(0:3,0:3,0:3,0:3),D5sd(0:3,0:3,0:3,0:3,3)
      real*8     p1,p2,p3,p4,s12,s23   
      complex*16 getmass2std,mm02,mm12,mm22,mm32
       
      integer    i,j,k,l,m
      integer    sym       
      common /sym/    sym

c      write(*,*) 'cDp12345sd in ',
c     &  p10,p21,p32,p30,p20,p31,m02,m12,m22,m32

      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      mm22 = getmass2std(m22)
      mm32 = getmass2std(m32)
      p1 = dreal(getmass2std(p10))      
      p2 = dreal(getmass2std(p21))      
      p3 = dreal(getmass2std(p32))     
      p4 = dreal(getmass2std(p30))
      s12 = dreal(getmass2std(p20))
      s23 = dreal(getmass2std(p31))

c note different order of arguments s23,s12
      call  xD012345_(p1,p2,p3,p4,s23,s12,mm02,mm12,mm22,mm32,
     &    D0sd,D1sd,D2sd,D3sd,D4sd,D5sd,C0sd,C1sd,C2sd,C3sd,C4sd,rank,1)
       
      D0      =  D0sd
      if (rank.eq.0) return
      D1(1)    =  D1sd(1)
      D1(2)    =  D1sd(2)
      D1(3)    =  D1sd(3)
      if (rank.eq.1) return
      D2(0,0) = D2sd(0,0)
      D2(1,1) = D2sd(1,1)
      D2(1,3) = D2sd(1,3)
      D2(1,2) = D2sd(1,2)
      D2(2,2) = D2sd(2,2)
      D2(2,3) = D2sd(2,3)
      D2(3,3) = D2sd(3,3)

      if (sym.gt.0) then
        D2(3,1) = D2sd(3,1)
        D2(2,1) = D2sd(2,1)
        D2(3,2) = D2sd(3,2)
      end if

      if (rank.eq.2) return
      D3(0,0,1) = D3sd(0,0,1)
      D3(0,0,2) = D3sd(0,0,2)
      D3(0,0,3) = D3sd(0,0,3)
      D3(1,1,1) = D3sd(1,1,1)
      D3(1,1,2) = D3sd(1,1,2)
      D3(1,1,3) = D3sd(1,1,3)
      D3(1,2,2) = D3sd(1,2,2)
      D3(1,2,3) = D3sd(1,2,3)
      D3(1,3,3) = D3sd(1,3,3)
      D3(2,2,2) = D3sd(2,2,2)
      D3(2,2,3) = D3sd(2,2,3)
      D3(2,3,3) = D3sd(2,3,3)
      D3(3,3,3) = D3sd(3,3,3)

      if (sym.gt.0) then
        D3(1,2,1) = D3sd(1,2,1)
        D3(2,1,1) = D3sd(2,1,1)
        D3(1,3,1) = D3sd(1,3,1)
        D3(3,1,1) = D3sd(3,1,1)
        D3(2,1,2) = D3sd(2,1,2)
        D3(2,2,1) = D3sd(2,2,1)
        D3(3,3,1) = D3sd(3,3,1)
        D3(3,1,3) = D3sd(3,1,3)
        D3(2,3,2) = D3sd(2,3,2)
        D3(3,2,2) = D3sd(3,2,2)
        D3(3,2,3) = D3sd(3,2,3)
        D3(3,3,2) = D3sd(3,3,2)
        D3(1,3,2) = D3sd(1,3,2)
        D3(2,1,3) = D3sd(2,1,3)
        D3(2,3,1) = D3sd(2,3,1)
        D3(3,2,1) = D3sd(3,2,1)
        D3(3,1,2) = D3sd(3,1,2)
      end if

      if (rank.eq.3) return
      D4(0,0,0,0) = D4sd(0,0,0,0)
      D4(0,0,1,1) = D4sd(0,0,1,1)
      D4(0,0,1,2) = D4sd(0,0,1,2)
      D4(0,0,1,3) = D4sd(0,0,1,3)
      D4(0,0,2,2) = D4sd(0,0,2,2)
      D4(0,0,2,3) = D4sd(0,0,2,3)
      D4(0,0,3,3) = D4sd(0,0,3,3)
      D4(1,1,1,1) = D4sd(1,1,1,1)
      D4(1,1,1,2) = D4sd(1,1,1,2)
      D4(1,1,1,3) = D4sd(1,1,1,3)
      D4(1,1,2,2) = D4sd(1,1,2,2)
      D4(1,1,2,3) = D4sd(1,1,2,3)
      D4(1,1,3,3) = D4sd(1,1,3,3)
      D4(1,2,2,2) = D4sd(1,2,2,2)
      D4(1,2,2,3) = D4sd(1,2,2,3)
      D4(1,2,3,3) = D4sd(1,2,3,3)
      D4(1,3,3,3) = D4sd(1,3,3,3)
      D4(2,2,2,2) = D4sd(2,2,2,2)
      D4(2,2,2,3) = D4sd(2,2,2,3)
      D4(2,2,3,3) = D4sd(2,2,3,3)
      D4(2,3,3,3) = D4sd(2,3,3,3)
      D4(3,3,3,3) = D4sd(3,3,3,3)

      if (sym.gt.0) then
        do i=1,3
        do j=1,3
          D4(0,0,i,j) = D4sd(0,0,i,j)
        do k=1,3
        do l=1,3
          D4(i,j,k,l) = D4sd(i,j,k,l)
        end do
        end do
        end do
        end do
      end if

      if (rank.eq.4) return
      D5(0,0,0,0,1) = D5sd(0,0,0,0,1)
      D5(0,0,0,0,2) = D5sd(0,0,0,0,2)
      D5(0,0,0,0,3) = D5sd(0,0,0,0,3)
      D5(0,0,1,1,1) = D5sd(0,0,1,1,1)
      D5(0,0,1,1,2) = D5sd(0,0,1,1,2)
      D5(0,0,1,1,3) = D5sd(0,0,1,1,3)
      D5(0,0,1,2,2) = D5sd(0,0,1,2,2)
      D5(0,0,1,2,3) = D5sd(0,0,1,2,3)
      D5(0,0,1,3,3) = D5sd(0,0,1,3,3)
      D5(0,0,2,2,2) = D5sd(0,0,2,2,2)
      D5(0,0,2,2,3) = D5sd(0,0,2,2,3)
      D5(0,0,2,3,3) = D5sd(0,0,2,3,3)
      D5(0,0,3,3,3) = D5sd(0,0,3,3,3)
      D5(1,1,1,1,1) = D5sd(1,1,1,1,1)
      D5(1,1,1,1,2) = D5sd(1,1,1,1,2)
      D5(1,1,1,1,3) = D5sd(1,1,1,1,3)
      D5(1,1,1,2,2) = D5sd(1,1,1,2,2)
      D5(1,1,1,2,3) = D5sd(1,1,1,2,3)
      D5(1,1,1,3,3) = D5sd(1,1,1,3,3)
      D5(1,1,2,2,2) = D5sd(1,1,2,2,2)
      D5(1,1,2,2,3) = D5sd(1,1,2,2,3)
      D5(1,1,2,3,3) = D5sd(1,1,2,3,3)
      D5(1,1,3,3,3) = D5sd(1,1,3,3,3)
      D5(1,2,2,2,2) = D5sd(1,2,2,2,2)
      D5(1,2,2,2,3) = D5sd(1,2,2,2,3)
      D5(1,2,2,3,3) = D5sd(1,2,2,3,3)
      D5(1,2,3,3,3) = D5sd(1,2,3,3,3)
      D5(1,3,3,3,3) = D5sd(1,3,3,3,3)
      D5(2,2,2,2,2) = D5sd(2,2,2,2,2)
      D5(2,2,2,2,3) = D5sd(2,2,2,2,3)
      D5(2,2,2,3,3) = D5sd(2,2,2,3,3)
      D5(2,2,3,3,3) = D5sd(2,2,3,3,3)
      D5(2,3,3,3,3) = D5sd(2,3,3,3,3)
      D5(3,3,3,3,3) = D5sd(3,3,3,3,3)
      if (sym.gt.0) then
        do i=1,3
          D5(0,0,0,0,i) = D5sd(0,0,0,0,i)
        do j=1,3
        do k=1,3
          D5(0,0,i,j,k) = D5sd(0,0,i,j,k)
        do l=1,3
        do m=1,3
          D5(i,j,k,l,m) = D5sd(i,j,k,l,m)
        end do
        end do
        end do
        end do
        end do
      end if
      end

************************************************************************
      function cE0sd(p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
     &                                    m02,m12,m22,m32,m42)
************************************************************************
*	interface to Stefans                                           *
*	5-POINT TENSOR COEFFICIENTS with complex masses                *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cE0sd
      complex*16 p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
     &           m02,m12,m22,m32,m42
      complex*16 D0sd,D1sd(3),D2sd(0:3,0:3),D3sd(0:3,0:3,3)
      complex*16 E0sd,E1sd(4),E2sd(0:4,0:4),E3sd(0:4,0:4,4),
     &           E4sd(0:4,0:4,0:4,0:4)
      real*8     p1,p2,p3,p4,p5,s12,s23,s34,s45,s15
      complex*16 getmass2std,mm02,mm12,mm22,mm32,mm42
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      mm22 = getmass2std(m22)
      mm32 = getmass2std(m32)
      mm42 = getmass2std(m42)
      p1 = dreal(getmass2std(p10))      
      p2 = dreal(getmass2std(p21))      
      p3 = dreal(getmass2std(p32))     
      p4 = dreal(getmass2std(p43))      
      p5 = dreal(getmass2std(p40))
      s12 = dreal(getmass2std(p20))
      s23 = dreal(getmass2std(p31))
      s34 = dreal(getmass2std(p42))
      s45 = dreal(getmass2std(p30))
      s15 = dreal(getmass2std(p41))

      call  xE01234_new(p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,
     &    mm02,mm12,mm22,mm32,mm42,
     &    E0sd,E1sd,E2sd,E3sd,E4sd,D0sd,D1sd,D2sd,D3sd,0,1)
       
      cE0sd      =  E0sd
      end

************************************************************************
      subroutine cEp1234sd(p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
     &                                    m02,m12,m22,m32,m42,
     &                E0,E1,E2,E3,E4,rank)
************************************************************************
*	interface to Stefans                                           *
*	5-POINT TENSOR COEFFICIENTS with complex masses                *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  13.10.04 ad        *
************************************************************************
      implicit   none
      integer    rank
      complex*16 p10,p21,p32,p43,p40,p20,p31,p42,p30,p41,
     &           m02,m12,m22,m32,m42
      complex*16 E0,E1(4),E2(0:4,0:4),E3(0:4,0:4,0:4),
     &           E4(0:4,0:4,0:4,0:4)
      complex*16 D0sd,D1sd(3),D2sd(0:3,0:3),D3sd(0:3,0:3,3)
      complex*16 E0sd,E1sd(4),E2sd(0:4,0:4),E3sd(0:4,0:4,4),
     &           E4sd(0:4,0:4,0:4,0:4)
      real*8     p1,p2,p3,p4,p5,s12,s23,s34,s45,s15
      integer    i,j,k,l
      complex*16 getmass2std,mm02,mm12,mm22,mm32,mm42
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      mm22 = getmass2std(m22)
      mm32 = getmass2std(m32)
      mm42 = getmass2std(m42)
      p1 = dreal(getmass2std(p10))      
      p2 = dreal(getmass2std(p21))      
      p3 = dreal(getmass2std(p32))     
      p4 = dreal(getmass2std(p43))      
      p5 = dreal(getmass2std(p40))
      s12 = dreal(getmass2std(p20))
      s23 = dreal(getmass2std(p31))
      s34 = dreal(getmass2std(p42))
      s45 = dreal(getmass2std(p30))
      s15 = dreal(getmass2std(p41))

      call  xE01234_new(p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,
     &    mm02,mm12,mm22,mm32,mm42,
     &    E0sd,E1sd,E2sd,E3sd,E4sd,D0sd,D1sd,D2sd,D3sd,rank,1)
       
      E0      =  E0sd
      if (rank.eq.0) return
      do 1 i=1,4
        E1(i)    =  E1sd(i)
 1    continue
      if (rank.eq.1) return
      
      E2(0,0)  =  E2sd(0,0)
      do 2 i=1,4
      do 2 j=1,4
        E2(i,j)  =  E2sd(i,j)
 2    continue

      if (rank.eq.2) return
      do 3 i=1,4
        E3(0,0,i)  =  E3sd(0,0,i)
      do 3 j=1,4
      do 3 k=1,4
        E3(i,j,k)  =  E3sd(i,j,k)
 3    continue

      if (rank.eq.3) return
        E4(0,0,0,0)=  E4sd(0,0,0,0)
      do 4 i=1,4
      do 4 j=1,4
        E4(0,0,i,j)=  E4sd(0,0,i,j)
      do 4 k=1,4
      do 4 l=1,4
        E4(i,j,k,l)=  E4sd(i,j,k,l)
 4    continue
      end

************************************************************************
      function cF0sd( p10,p21,p32,p43,p54,p50
     &          ,p20,p31,p42,p53,p40,p51,p30,p41,p52
     &          ,m02,m12,m22,m32,m42,m52)
************************************************************************
*	interface to Stefans                                           *
*	6-POINT TENSOR COEFFICIENTS with complex masses                *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      complex*16 cF0sd
      complex*16 p10,p21,p32,p43,p54,p50
     &          ,p20,p31,p42,p53,p40,p51,p30,p41,p52
     &          ,m02,m12,m22,m32,m42,m52
      complex*16 E0sd,E1sd(4),E2sd(0:4,0:4)
c     &    ,E3sd(0:4,0:4,4)
      complex*16 F0sd,F1sd(5),F2sd(0:5,0:5),F3sd(0:5,0:5,5)
c     &           F4sd(0:5,0:5,0:5,0:5)
      real*8     p1,p2,p3,p4,p5,p6,s12,s23,s34,s45,s56,s16,
     &	 s123,s234,s345
      complex*16 getmass2std,mm02,mm12,mm22,mm32,mm42,mm52
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      mm22 = getmass2std(m22)
      mm32 = getmass2std(m32)
      mm42 = getmass2std(m42)
      mm52 = getmass2std(m52)
      p1 = dreal(getmass2std(p10))      
      p2 = dreal(getmass2std(p21))      
      p3 = dreal(getmass2std(p32))     
      p4 = dreal(getmass2std(p43))      
      p5 = dreal(getmass2std(p54))
      p6 = dreal(getmass2std(p50))
      s12 = dreal(getmass2std(p20))
      s23 = dreal(getmass2std(p31))
      s34 = dreal(getmass2std(p42))
      s45 = dreal(getmass2std(p53))
      s56 = dreal(getmass2std(p40))
      s16 = dreal(getmass2std(p51))
      s123 = dreal(getmass2std(p30))
      s234 = dreal(getmass2std(p41))
      s345 = dreal(getmass2std(p52))

      call  xF0123_new(p1,p2,p3,p4,p5,p6,s12,s23,s34,s45,s56,s16,
     &	 s123,s234,s345,mm02,mm12,mm22,mm32,mm42,mm52,
     &    F0sd,F1sd,F2sd,F3sd,E0sd,E1sd,E2sd,0,1)
       
      cF0sd  =  F0sd

      end
************************************************************************
      subroutine cFp123sd( p10,p21,p32,p43,p54,p50
     &          ,p20,p31,p42,p53,p40,p51,p30,p41,p52
     &          ,m02,m12,m22,m32,m42,m52,
     &                F0,F1,F2,F3,rank)
************************************************************************
*	interface to Stefans                                           *
*	6-POINT TENSOR COEFFICIENTS with complex masses                *
*----------------------------------------------------------------------*
*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
************************************************************************
      implicit   none
      integer    rank
      complex*16 p10,p21,p32,p43,p54,p50
     &          ,p20,p31,p42,p53,p40,p51,p30,p41,p52
     &          ,m02,m12,m22,m32,m42,m52
      complex*16 F0,F1(5),F2(0:5,0:5),F3(0:5,0:5,0:5)
c     &           F4(0:5,0:5,0:5,0:5)
      complex*16 E0sd,E1sd(4),E2sd(0:4,0:4)
c     &      ,E3sd(0:4,0:4,4)
      complex*16 F0sd,F1sd(5),F2sd(0:5,0:5),F3sd(0:5,0:5,5)
c     &           F4sd(0:5,0:5,0:5,0:5)
      real*8     p1,p2,p3,p4,p5,p6,s12,s23,s34,s45,s56,s16,
     &	 s123,s234,s345
      integer    i,j,k
      complex*16 getmass2std,mm02,mm12,mm22,mm32,mm42,mm52
       
      mm02 = getmass2std(m02)
      mm12 = getmass2std(m12)
      mm22 = getmass2std(m22)
      mm32 = getmass2std(m32)
      mm42 = getmass2std(m42)
      mm52 = getmass2std(m52)
      p1 = dreal(getmass2std(p10))      
      p2 = dreal(getmass2std(p21))      
      p3 = dreal(getmass2std(p32))     
      p4 = dreal(getmass2std(p43))      
      p5 = dreal(getmass2std(p54))
      p6 = dreal(getmass2std(p50))
      s12 = dreal(getmass2std(p20))
      s23 = dreal(getmass2std(p31))
      s34 = dreal(getmass2std(p42))
      s45 = dreal(getmass2std(p53))
      s56 = dreal(getmass2std(p40))
      s16 = dreal(getmass2std(p51))
      s123 = dreal(getmass2std(p30))
      s234 = dreal(getmass2std(p41))
      s345 = dreal(getmass2std(p52))

      call  xF0123_new(p1,p2,p3,p4,p5,p6,s12,s23,s34,s45,s56,s16,
     &	 s123,s234,s345,mm02,mm12,mm22,mm32,mm42,mm52,
     &    F0sd,F1sd,F2sd,F3sd,E0sd,E1sd,E2sd,rank,1)
       
      F0      =  F0sd

      if (rank.eq.0) return
      do 1 i=1,5
        F1(i)    =  F1sd(i)
 1    continue

      if (rank.eq.1) return
      F2(0,0) = F2sd(0,0)
      do 2 i=1,5
      do 2 j=1,5
        F2(i,j)  =  F2sd(i,j)
 2    continue

      if (rank.eq.2) return
      do 3 i=1,5
        F3(0,0,i)  =  F3sd(0,0,i)
      do 3 j=1,5
      do 3 k=1,5
        F3(i,j,k)  =  F3sd(i,j,k)
 3    continue

      end



c>************************************************************************
c>      function cB0mr(p10,m02,m12)
c>************************************************************************
c>*	interface to Stefans                                           *
c>*	scalar 2-point function with complex masses                    *
c>*----------------------------------------------------------------------*
c>*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
c>************************************************************************
c>      implicit   none
c>      complex*16 cB0mr,yB0_
c>      complex*16 p10,m02,m12
c>      complex*16 getmass2std,mm02,mm12
c>      real*8     p1  
c>       
c>c      write(*,*) 'cB0mr in ',p10,m02,m12
c>      mm02 = getmass2std(m02)
c>      mm12 = getmass2std(m12)
c>      p1   = dreal(getmass2std(p10))
c>c      write(*,*) 'cB0mr out',p1,mm02,mm12
c>      cB0mr = yB0_(p1,mm02,mm12,1)
c>      
c>      end
c>
c>************************************************************************
c>      function cDB0mr(p10,m02,m12)
c>************************************************************************
c>*	interface to Stefans                                           *
c>*	momentum derivative of scalar 2-point function                 *
c>*       with complex masses                                            *
c>*----------------------------------------------------------------------*
c>*     01.07.04  Ansgar Denner         last changed  01.07.04 ad        *
c>************************************************************************
c>      implicit   none
c>      complex*16 cDB0mr,yDB0_
c>      complex*16 p10,m02,m12
c>      complex*16 getmass2std,mm02,mm12
c>      real*8     p1  
c>       
c>c      write(*,*) 'cDB0mr in ',p10,m02,m12
c>      mm02 = getmass2std(m02)
c>      mm12 = getmass2std(m12)
c>      p1   = dreal(getmass2std(p10))
c>c      write(*,*) 'cDB0mr out',p1,mm02,mm12
c>       
c>      cDB0mr = yDB0_(p1,mm02,mm12,1)
c>      
c>      end 
