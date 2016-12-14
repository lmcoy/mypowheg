* -*-Fortran-*-

***********************************************************************
*     file common_coli.h                                              *
*     contains global common blocks for COLI                          *
*---------------------------------------------------------------------*
*     04.08.08  Ansgar Denner     last changed  05.08.08              *
***********************************************************************


c information output switch
      logical   coliinfo
      common/info_coli/coliinfo

c output channels

#ifdef SING
c regularization parameters
      real*8    deltauv,delta2ir,delta1ir,colishiftms2
      common/sing_coli/deltauv,delta2ir,delta1ir,colishiftms2
#endif
c regularization parameters
      real*8    muuv2,muir2
      common/dimreg_coli/muuv2,muir2

c infinitesimal masses
      integer    ncoliminf
      common /ncoliminf/     ncoliminf
      complex*16 coliminf(10)
      common /coliminf/     coliminf    
      complex*16 coliminf2(10)
      common /coliminf2/    coliminf2    
      complex*16 coliminffix(10)
      common /coliminffix/  coliminffix    
      complex*16 coliminffix2(10)
      common /coliminffix2/ coliminffix2    

c scaling of infinitesimal masses
      real*8     coliminfscale,coliminfscale2
      common /colimsing/ coliminfscale,coliminfscale2

c infinitesimal parameters for D0
c      complex*16 ps12,ps23,ps34,ps14,ps13,ps24,ms12,ms22,ms32,ms42

c      common /d0infpar/ ps12,ps23,ps34,ps14,ps13,ps24,
c     &    ms12,ms22,ms32,ms42

