* -*-Fortran-*-

***********************************************************************
*     file checkparams_coli.h                                         *
*     global input parameters for check option of  COLI               *
*---------------------------------------------------------------------*
*     02.05.08  Ansgar Denner     last changed  02.05.08              *
***********************************************************************

#ifdef CHECK

      logical argcheck,conscheck,oldcheck
      common /check/argcheck,conscheck,oldcheck

      integer   testout,errout	
      parameter (testout=20, errout=6)

      real*8 testacc
      parameter (testacc=1d-5)   
c      parameter (testacc=1d-11)   

      argcheck=.false.
      conscheck=.true.
      oldcheck=.true.
#endif
