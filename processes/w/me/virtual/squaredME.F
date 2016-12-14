#define COLILib
C#define DDLib 

******************************************************************************************+

      Subroutine init_all(alfa_in,mw_in,gw_in,mz_in,gz_in,
     & me_in,mm_in,ml_in,mu_in,mc_in,mt_in,md_in,ms_in,mb_in,mh_in,gh_in,
     & inveps_in, lambda_in, dr_in)
      Implicit None     
*
* Call before first call to Msquared to provide necessary input 
* for common blocks and initialize loop routines
*
#include "params.h"
#include "amps/ampcommon.h"

      double precision alfa_in
      double precision mw_in,gw_in,mz_in,gz_in
      double precision me_in,mm_in,ml_in
      double precision mu_in,mc_in,mt_in
      double precision md_in,ms_in,mb_in
      double precision mh_in,gh_in
      integer dr_in
      integer dr
      double precision lambda_in, inveps_in

      double precision nflavour_out,mureg_out,Qren2
      common /ct_common/nflavour_out,mureg_out,Qren2
      
      integer i1,gen,procnr
      character*5 name(np)
#include "regcommon.h"

      PI=3.1415926535897932384D0
      sqrt2=sqrt(2D0)

      nexternal(1)=np-1

*     renormalization scale
      mureg=1D0

*     deltareg=2/epsilon-gamma_E+log(4*pi)
      deltareg=0D0

      epsfac=inveps_in
      epsfac2=inveps_in*inveps_in

      dr = dr_in
*     photon mass
      if (dr.eq.1) then
	     lambda=1d0
      else
	     lambda=lambda_in
      endif

*     number of light flavors
      nflavour=5D0
      
      nflavour_out=nflavour
      mureg_out=mureg

      MWr=MW_in
      MZr=MZ_in
      MH=MH_in
      GW=GW_in
      GZ=GZ_in
      GH=GH_in

      ME=ME_in
      MM=MM_in
      ML=ML_in

      MU=MU_in
      MC=MC_in
      MT=MT_in

      MD=MD_in
      MS=MS_in
      MB=MB_in

      alfa=alfa_in
      
      cwr = mwr/mzr
      cw2r = cwr**2
      sw2r = 1d0-cw2r
      swr=dsqrt(sw2r) 

      MW=MWr
      MZ=MZr

      MW2r=MWr**2
      MZ2r=MZr**2
      MH2=MH**2

      ME2=ME**2
      MM2=MM**2
      ML2=ML**2

      MU2=MU**2
      MC2=MC**2
      MT2=MT**2

      MD2=MD**2
      MS2=MS**2
      MB2=MB**2

      MLE(1)=ME
      MLE(2)=MM
      MLE(3)=ML

      MLE2(1)=ME2
      MLE2(2)=MM2
      MLE2(3)=ML2

      MQU(1)=MU
      MQU(2)=MC
      MQU(3)=MT

      MQU2(1)=MU2
      MQU2(2)=MC2
      MQU2(3)=MT2

      MQD(1)=MD
      MQD(2)=MS
      MQD(3)=MB

      MQD2(1)=MD2
      MQD2(2)=MS2
      MQD2(3)=MB2
     

      alfa2=alfa**2

      el=dsqrt(4D0*pi*alfa)

      MW2c=MW2r-(0D0,1D0)*GW*MWr
      MZ2c=MZ2r-(0D0,1D0)*GZ*MZr

      MWc=sqrt(MW2c)
      MZc=sqrt(MZ2c)

      CW2c=MW2c/MZ2c
      SW2c=1D0-CW2c
      CWc=sqrt(CW2c)
      SWc=sqrt(SW2c)

      do gen=1,1
       do procnr=1,nproctree
        do i1=1,np						   
           name(i1)=subprocs(i1,procnr,gen)		   
           procmass(gen,procnr,i1)=0D0 			   
           if(name(i1).eq." up  ".or.name(i1).eq."~up  ")then     
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." do  ".or.name(i1).eq."~do  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." ch  ".or.name(i1).eq."~ch  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." st  ".or.name(i1).eq."~st  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." bo  ".or.name(i1).eq."~bo  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." to  ".or.name(i1).eq."~to  ")then 
              procmass(gen,procnr,i1)=MT
           elseif(name(i1).eq." el  ".or.name(i1).eq."~el  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." ve  ".or.name(i1).eq."~ve  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." mu  ".or.name(i1).eq."~mu  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." ta  ".or.name(i1).eq."~ta  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." W+  ".or.name(i1).eq."~W+  ")then 
              procmass(gen,procnr,i1)=MWr
           elseif(name(i1).eq." Z0  ")then			   
              procmass(gen,procnr,i1)=MZr
           elseif(name(i1).eq." H0  ")then			   
              procmass(gen,procnr,i1)=MH
           elseif(name(i1).eq." gl  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq." ga  ")then 
              procmass(gen,procnr,i1)=0D0
           elseif(name(i1).eq."     ")then			   
              procmass(gen,procnr,i1)=0D0 		   
           else
	      Print *,"Particle identification failed in mass assignment"
	      Stop
	   endif						   
        enddo  						   
       enddo
      enddo

*initialize amplitudes
      call init_amps

      End

******************************************************************************************+
      Subroutine getdZR(CT)
      Implicit None
#include "params.h"
#include "amps/ampcommon.h"
#include "amps/ampvars.h" 
#include "amps/ampcom.h"
      double precision CT

c     call abbr_const()
      CT = dZR1
c      write(*,*) dZR1
      End

******************************************************************************************+

      Subroutine set_alphas(Alfas_in,Qren2_in)
      Implicit None     
*
* Call to provide Alphas and the squared renormalization scale
* (at each phase-space point)
*
      integer dr
      double precision Alfas,Alfas2,GS,Alfas_in,Pi
      double precision nflavour,mureg,Qren2,Qren2_in
#include "regcommon.h"

      common /my_alfas_common/Alfas,Alfas2,GS
      common /ct_common/nflavour,mureg,Qren2 

      PI=3.1415926535897932384D0
      
      Alfas=Alfas_in
      Alfas2=Alfas_in**2
      gs=dsqrt(4D0*pi*Alfas)
      
      Qren2=Qren2_in
    
      if (dr.eq.1) then
          call abbrconst_mu(Qren2)
      endif

      End

      subroutine abbrconst_mu(mu2_in)
#include "params.h" 
      integer dr
      double precision mu2_in
#include "amps/ampvars.h" 
#include "amps/ampcommon.h"
#include "amps/ampcom.h" 
#include "regcommon.h"

      lambda=sqrt(mu2_in)
      Cpl105 = 1/Pi*Log(lambda**2)

      dZfL1(2,1,1) = 2*Cpl9*LInt185 + Cpl9*LInt186 - 
     -   0.25D0*Cpl104 + 0.0625D0*Cpl37 + 0.25D0*Cpl38 + 
     -   LInt184*(-(0.5D0*Cpl104) + 0.125D0*Cpl37 + 0.5D0*Cpl38) + 
     -   0.125D0*Cpl8 + 0.25D0*Cpl9 + 0.25D0*(Cpl8*LInt182) + 
     -   0.5D0*(Cpl9*LInt183)

	If(dr.Eq.1) Then
	 dZfL1(2,1,1) = dZfL1(2,1,1) -
     -  (2*Cpl9*LInt185 + Cpl9*LInt186 + 0.25D0*Cpl9 + 
     -   0.5D0*(Cpl9*LInt183)) +
     -   Alfa*(0.25D0*Cpl105 - 0.25D0*Cpl106 + 
     -        0.25D0*Cpl107 - 0.25D0*Cpl108)
	Endif

        dZfL1(3,1,1) = -(0.16666666666666666D0*Cpl104) + 
     -   0.0625D0*Cpl37 + 0.1111111111111111D0*Cpl38 + 
     -   LInt189*(-(0.3333333333333333D0*Cpl104) + 0.125D0*Cpl37 + 
     -      0.2222222222222222D0*Cpl38) + 0.125D0*Cpl8 + 
     -   0.1111111111111111D0*Cpl9 + 0.25D0*(Cpl8*LInt187) + 
     -   0.2222222222222222D0*(Cpl9*LInt188) + 
     -   0.8888888888888888D0*(Cpl9*LInt190) + 
     -   0.4444444444444444D0*(Cpl9*LInt191)
       If(dr.Eq.1) Then
        dZfL1(3,1,1) = dZfL1(3,1,1) -
     - (0.1111111111111111D0*Cpl9+ 
     -   0.2222222222222222D0*(Cpl9*LInt188) + 
     -   0.8888888888888888D0*(Cpl9*LInt190) + 
     -   0.4444444444444444D0*(Cpl9*LInt191)) +
     -   4D0/9D0*Alfa*(0.25D0*Cpl105 - 0.25D0*Cpl106 + 
     -        0.25D0*Cpl107 - 0.25D0*Cpl108)
       Endif

#ifdef DEBUG
	print *, 'dZfL1(3,1,1) =', dZfL1(3,1,1)
#endif

        dZfL1(4,1,1) = -(0.08333333333333333D0*Cpl104) + 
     -   0.0625D0*Cpl37 + 0.027777777777777776D0*Cpl38 + 
     -   LInt193*(-(0.16666666666666666D0*Cpl104) + 0.125D0*Cpl37 + 
     -      0.05555555555555555D0*Cpl38) + 0.125D0*Cpl8 + 
     -   0.027777777777777776D0*Cpl9 + 
     -   0.05555555555555555D0*(Cpl9*LInt192) + 
     -   0.25D0*(Cpl8*LInt194) + 0.2222222222222222D0*(Cpl9*LInt195) + 
     -   0.1111111111111111D0*(Cpl9*LInt196)
       If(dr.Eq.1) Then
        dZfL1(4,1,1) = dZfL1(4,1,1) -
     - (0.027777777777777776D0*Cpl9 + 
     -   0.05555555555555555D0*(Cpl9*LInt192) + 
     -   0.2222222222222222D0*(Cpl9*LInt195) + 
     -   0.1111111111111111D0*(Cpl9*LInt196)) + 
     -   1D0/9D0*Alfa*(0.25D0*Cpl105 - 0.25D0*Cpl106 + 
     -        0.25D0*Cpl107 - 0.25D0*Cpl108)
       Endif
       
#ifdef DEBUG
	print *, 'dZfL1(4,1,1) =', dZfL1(4,1,1)
#endif

        dZfL1dimreg(3,1,1) = 
     -  0.3333333333333333D0*Cpl105 - 0.3333333333333333D0*Cpl106 + 
     -   0.3333333333333333D0*Cpl107 - 0.3333333333333333D0*Cpl108

#ifdef DEBUG
	print *, 'dZfL1dimreg(3,1,1) =', dZfL1dimreg(3,1,1)
#endif

        dZfL1dimreg(4,1,1) = 
     -  0.3333333333333333D0*Cpl105 - 0.3333333333333333D0*Cpl106 + 
     -   0.3333333333333333D0*Cpl107 - 0.3333333333333333D0*Cpl108

#ifdef COLILib
      call setmphoton2(lambda**2)           
#endif
#ifdef DDLib
*Initialize Stefan's libraries
      Call DDsetparam(deltareg,mureg,epsfac2,epsfac,
     &			      lambda**2,smallm2)
      cacc = 1d-10
      dacc = 1d-10
      mode34 = 0
      mode5  = 0
      mode6  = 0
      outlevel = 0
      call DDsetmode(cacc,dacc,mode34,mode5,mode6,outlevel)
#endif
      end

******************************************************************************************+

      subroutine Msquared(msqu0,msqu1,k, QCDorEW)
      implicit none                   
* Calculates:
* msqu0(nproc): Born matrix element squared
* msqu1(nproc): virtual 1-loop QCD corrected matrix element squared
*       nproc=1: ubar u
*       nproc=2: dbar d
* for momentum configuration k(np,0:3) 
* where np is the number of particles at born + 1
* np=1: anti-quark
* np=2: quark
* np=3: anti-lepton
* np=4: lepton
* np=5: empty
* You have to call "init_all" first!
* You have to provide alphas and the renormalization scale 
* with "set_alphas"!

#include "params.h"

      integer procnr,helnr,nc1
      integer resnr,nrproc      
      double precision k(np,0:3),msqu0(nproc),msqu1(nproc)
      double complex M0(nproc,ncol,nhel),M1(nproc,ncol,nhel)
      double complex M1virtQCD1(nproc,ncol,nhel)
      double complex M1ampsfacctQCD1(nproc,ncol,nhel)
      double complex M1virtEW1(nproc,ncol,nhel)
      double complex M1ampsfacctEW1(nproc,ncol,nhel)
      integer QCDorEW


      resnr=0
      
      do procnr=1,nproctree
         msqu0(procnr)=0D0
         msqu1(procnr)=0D0
      enddo
 
      Do procnr=1,nproctree
	 do helnr=1,nheltree
	    do nc1=1,ncol
               M0(procnr,nc1,helnr)=(0D0,0D0)
               M1(procnr,nc1,helnr)=(0D0,0D0)
            enddo
         enddo
      enddo
      
      
       
      If(QCDorEW.eq.1) Then
	call setzero(M1virtQCD1)
	call setzero(M1ampsfacctQCD1)
      Else if(QCDorEW.eq.2) Then
	call setzero(M1virtEW1)
	call setzero(M1ampsfacctEW1)
      Endif
       
#ifdef COLILib
      call cacheinit(0)
#endif

      call setprods(k)

      call abbrborn1(resnr)
      call born1(M0,resnr)

      If(QCDorEW.eq.1) Then

	call abbrvirtQCD1(resnr)
	call virtQCD1(M1virtQCD1,resnr)
	call sumamplitudes(M1,M1virtQCD1)

	call abbrampsfacctQCD1(resnr)
	call ampsfacctQCD1(M1ampsfacctQCD1,resnr)
	call sumamplitudes(M1,M1ampsfacctQCD1)

      Else if(QCDorEW.eq.2) Then

	call abbrvirtEW1(resnr)
	call virtEW1(M1virtEW1,resnr)
	call sumamplitudes(M1,M1virtEW1)

	call abbrampsfacctEW1(resnr)
	call ampsfacctEW1(M1ampsfacctEW1,resnr)
	call sumamplitudes(M1,M1ampsfacctEW1)

      Endif

      call squareamp(msqu0,msqu1,M0,M1)

      end  

****************************************************************************************************************

      subroutine squareamp(msqu0,msqu1,amp0,amp1)
      implicit none

#include "params.h"

      integer nrproc
      double precision msqu0(nproc),msqu1(nproc)
      double complex amp0(nproc,ncol,nhel),amp1(nproc,ncol,nhel)
     
      integer procnr,helnr,nc1,nc2
      double complex res0,res1
      double precision re

#include "amps/ampcommon.h"

      Do procnr=1,nproctree							    
         do helnr=1,nheltree							    
            do nc1=1,ncol							    
               res0=0D0 							    
               res1=0D0 							    
               do nc2=1,ncol							  
        	  res0=res0+Cfac(nc2,nc1,procnr,1)   
     &  		     *amp0(procnr,nc2,helnr)				    
        	  res1=res1+Cfac(nc2,nc1,procnr,1)   
     &  		     *amp1(procnr,nc2,helnr)				    
               enddo								    
               Msqu0(procnr)=Msqu0(procnr)+res0*conjg(res0)			    
               Msqu1(procnr)=Msqu1(procnr)+2D0*Re(res0,res1)			 
            enddo								    
         enddo  								    
      enddo									    
      
      end

**************************************************************************************************** 
 
      subroutine sumamplitudes(Amp,M)
      implicit none

#include "params.h"

      double complex Amp(nproc,ncol,nhel),M(nproc,ncol,nhel)
      integer procnr,helnr,colnr

       Do procnr=1,nproctree
	  do helnr=1,nheltree
	     do colnr=1,ncol
                amp(procnr,colnr,helnr)
     &           =amp(procnr,colnr,helnr)+M(procnr,colnr,helnr)
             enddo
          enddo 
       enddo
      end

******************************************************************************************+

      subroutine setzero(M)       
      implicit none

#include "params.h"

      integer procnr,helnr,nc1
      double complex M(nproc,ncol,nhel)

       Do procnr=1,nproctree
	  do helnr=1,nheltree
	     do nc1=1,ncol
                M(procnr,nc1,helnr)=(0D0,0D0)
             enddo
          enddo
       enddo
      end  

*************************************************************************************************

      subroutine init_amps
      implicit none
*Intializes the amplitudes

#include "params.h"
#include "amps/ampcommon.h"

      integer gen,ii,mode34,mode5,mode6,outlevel,iimax
      double precision getlambda,smallm2(10),cacc,dacc
      double complex commass 

      common/smallm2_common/smallm2,iimax
      
#ifdef COLILib
*Initialize Ansgar's libraries
      call defcoli(1)
c      call cacheon
      call countinit
#endif

*Initialize UV and infrared regulators

#ifdef COLILib
      call setmudim2(mureg)
      call setmphoton2(lambda**2)
      call setdeltair_coli(epsfac,epsfac2)
      call setdeltauv_coli(deltareg)            
#endif
      
*set small masses for Ansgar's and Stefan's loop integrals
      do ii=1,10
         smallm2(ii)=0D0
      enddo
      ii=0
      commass=MU2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=MU2

      commass=MD2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=MD2

C      commass=MS2
#ifdef COLILib
C      call setcminf2(commass)
#endif
C      ii=ii+1
C      smallm2(ii)=MS2

C      commass=MC2
#ifdef COLILib
C      call setcminf2(commass)
#endif
C      ii=ii+1
C      smallm2(ii)=MC2

      commass=MB2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=MB2

      commass=ME2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=ME2

      commass=ML2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=ML2

      commass=MM2
#ifdef COLILib
      call setcminf2(commass)
#endif
      ii=ii+1
      smallm2(ii)=MM2
            
      iimax=ii

      call compress_smallm2(smallm2,iimax)

#ifdef DDLib
*Initialize Stefan's libraries
      Call DDsetparam(deltareg,mureg,epsfac2,epsfac,
     &			      lambda**2,smallm2)
      cacc = 1d-10
      dacc = 1d-10
      mode34 = 0
      mode5  = 0
      mode6  = 0
      outlevel = 0
      call DDsetmode(cacc,dacc,mode34,mode5,mode6,outlevel)
#endif

      end
 
******************************************************************************************+

      Subroutine compress_smallm2(smallm2,iimax)
      Implicit None
      
      double precision smallm2(10)
      integer iimax,ii,jj,kk,keep_index

      keep_index=0
      do ii=1,10
 99     continue
        If(smallm2(ii).NE.0) Then      
	  If(ii.GT.0) Then
	    do jj=1,ii-1
	       If(smallm2(ii).EQ.smallm2(jj)) Then
	         do kk=ii,9
	           smallm2(kk)=smallm2(kk+1)
	         enddo
	         smallm2(10)=0D0
C	         Print ('(A11,I2,A42)'),"small mass ",ii+keep_index,
C     &                   " is eliminated (equal to other small mass)"
		 keep_index=keep_index+1
		 goto 99
	       Endif
	    enddo
	  Endif
	Else 
	  iimax=ii-1
	  Return
	Endif
      Enddo
      
      End      

      function delta(s1,s2)
      implicit none
      integer delta,s1,s2
c*Kronecker delta
      if (s1.eq.s2) then
        delta=1D0
      else
         delta=0D0
      endif
      end


      function Re(a,b)
      implicit none
c*real part
      double precision Re
      double complex a,b

      Re=(conjg(a)*b+a*conjg(b))/2D0
      end

      subroutine getreg(dr_out)
           integer dr
           integer dr_out
#include "regcommon.h"
           dr_out = dr
      end
